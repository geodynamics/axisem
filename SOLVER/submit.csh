#!/bin/csh -f

set homedir = $PWD

if ( ${#argv} < 1 || "$1" == "-h" ) then
    echo ""
    echo "=================================================="
    echo " Argument 1:  directory name for the simulation"
    echo ""
    echo " Optional arguments after directory name: "
    echo " -q <queue>, where <queue> can be:"
    echo "      'lsf': submit to lsf queue using bsub"
    echo "       Default: submit locally"
    echo "=================================================="
    echo""
    exit
else if ( -d $1) then
    echo " Run or directory" $1 "exists....... its content:"
    ls $1
    exit
endif

set datapath = `grep "DATA_DIR" inparam_advanced  |awk '{print $2}'| sed 's/\"//g'`
set infopath = `grep "INFO_DIR" inparam_advanced |awk '{print $2}'| sed 's/\"//g'`
set meshdir = "MESHES/"`grep "MESHNAME" inparam_basic | awk '{print $2}'`

set svnrevision = `svnversion`
echo $svnrevision "SVN_VERSION      " > runinfo
set username = `whoami`
echo $username "USER_NAME        " >> runinfo
set hostname = `hostname`
echo $hostname "HOST_NAME        " >> runinfo


if ( -d $meshdir) then
    echo "Using mesh " $meshdir
else
    echo "Mesh " $meshdir " not found."
    echo "Available meshes:"
    ls MESHES
    exit
endif

set bgmodel = `grep BACKGROUND_MODEL $meshdir/inparam_mesh | awk '{print $2}'`

if ( ! -f inparam_hetero) then 
  cp inparam_hetero.TEMPLATE inparam_hetero
endif

if ( ! -f inparam_xdmf) then 
  cp inparam_xdmf.TEMPLATE inparam_xdmf
endif


# if the mesh has different mesh_params.h, copy here
if ( ! -f mesh_params.h || `diff mesh_params.h $meshdir/mesh_params.h | wc -l` != "0" ) then
  echo 'copying mesh_params.h from ' $meshdir
  cp $meshdir/mesh_params.h .
endif

# if the mesh has different background_models.f90, copy over
if ( `diff background_models.f90 $meshdir/background_models.f90 | wc -l` != "0" ) then
  echo 'copying background_models.f90 from ' $meshdir
  cp $meshdir/background_models.f90 .
endif

# Check arguments: source types and submission queues
set newqueue = 'false'
if ( "$2" == '-q') then
    set queue = $3
    set newqueue = 'true'
endif

set multisrc = 'false'
set srctype = `grep "SIMULATION_TYPE" inparam_basic |awk '{print $2}'`
if ( $srctype == 'single') then
    set multisrc = 'false'
    set src_file_type = 'sourceparams'
else if ( $srctype == 'force') then
    set multisrc = 'true'
    set src_file_type = 'sourceparams'
else if ( $srctype == 'moment') then
    set multisrc = 'true'
    set src_file_type = 'cmtsolut'
endif

if ( $newqueue == 'true' ) then 
	echo "Submitting to queue type" $queue
endif

# Run make to see whether the code has to be rebuilt and if so, do it.
# If 'make' returns an Error (something >0), then exit.

# MvD: I do not get this: it checks == 0 where 0 is the status when exited without
# problems???
if ( { make -j } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif

if ( $src_file_type == 'cmtsolut' ) then
    set srcfile = 'CMTSOLUTION'
else if ( $src_file_type == 'sourceparams' ) then
    set srcfile = 'sourceparams.dat'
else if ( $src_file_type == 'finfault' ) then
    set srcfile = 'finite_fault.dat'
endif 

if ( ! -f $homedir/$srcfile ) then 
    echo "file $srcfile does not exist"
    exit
endif

# identify receiver input file
set rec_file_type = `grep "RECFILE_TYPE" inparam_basic |awk '{print $2}'`
echo "Receiver file type:" $rec_file_type

if ( $rec_file_type == 'colatlon' ) then
    set recfile = 'receivers.dat'
else if ( $rec_file_type == 'stations' ) then
    set recfile = 'STATIONS'
else if ( $rec_file_type == 'database' ) then
    set recfile = 'database'
    echo "this is a dummy database receiver file" >! $homedir/$recfile
endif

if ( ! -f $homedir/$recfile ) then 
    echo "file $recfile does not exist"
    exit
endif
echo "Source file:" $srcfile, "Receiver file:" $recfile

set num_src = 1
set num_src_arr = ( 1 )
if ( $multisrc == 'true' ) then
    # multiple simulations
    echo "setting up multiple simulations for full" $srctype "source type"
    if ( $srctype == 'moment' ) then 
        set mij_sourceparams = ( 0. 0. 0. 0. 0. 0. )
        set map_mij = ( 1 2 4 6 )
        set numsim = 4
        set srcapp = ( MZZ MXX_P_MYY MXZ_MYZ MXY_MXX_M_MYY )
        set srctype1 = ( "'monopole'" "'monopole'" "'dipole'" "'quadpole'")
        set srctype2 = ( "'mzz'" "'mxx_p_myy'" "'mxz'" "'mxy'" )

    else if ( $srctype == 'force' ) then 
        set numsim = 2
        set srcapp = ( PZ PX )
        set srctype1 = ( "'monopole'" "'dipole'" )
        set srctype2 = ( "'vertforce'" "'xforce'" )

    else if ( $srctype == 'finfault') then 
        set num_src = `cat finite_fault.dat | grep fflt_num |awk '{print $1}'`
        set numsim = 4
        set num_src_arr = `tail -n $num_src finite_fault.dat | awk '{print $1}'`
        set srcapp = ( MZZ MXX_P_MYY MXZ_MYZ MXY_MXX_M_MYY )      
        set srctype1 = ( "'monopole'" "'monopole'" "'dipole'" "'quadpole'")
        set srctype2 = ( "'mzz'" "'mxx_p_myy'" "'mxz'" "'mxy'" )
    else
        echo " Unrecognized source type" $srctype
        echo " Choose either 'moment', 'force', 'finfault' or leave blank for one simulation as in sourceparams.dat"
        exit
    endif

   if ( $src_file_type == 'cmtsolut' ) then 
       set mom_tens1 = ( "Mrr:" "Mtt:" "Mpp:" "Mrt:" "Mrp:" "Mtp:" ) 
       set map_mij = ( 1 2 4 6 )
       set map_mij_name = ( "Mrr:" "Mtt:" "Mrt:" "Mtp:" ) 
       head -n 7 $homedir/CMTSOLUTION > 'CMTSOLUTION.MIJ' 
       foreach el (${mom_tens1}) 
           echo $el "     0.000e+27" >> 'CMTSOLUTION.MIJ'
       end
   endif
else if ( $multisrc == 'false' ) then
    # one simulation
    set numsim = 1; 
    set srctype1 = `grep "excitation type" sourceparams.dat |awk '{print $1}'`
    set srctype2 = `head -n 3 sourceparams.dat  | tail -n 1 |awk '{print $1}'` 
    set srcapp = ( "./"  )
endif 

echo 'source names:' $srcapp
echo 'source radiation types:' $srctype1
echo 'source components:' $srctype2

mkdir $1
cd $1
set mainrundir = $PWD

# make sure moment tensor is copied correctly
cp -p $homedir/$srcfile $mainrundir/



# Prepare and copy relevant files for each simulation
foreach isrc (${num_src_arr})
    set i = 0
    foreach isim  (${srcapp})

        @ i ++

        set num = 6
        echo ""
        echo "Setting up simulation" $isim
        # construct different source file for each simulation
        if  ( $multisrc == 'true' ) then
            echo "constructing separate source files for" $isim 

            if ( $src_file_type == 'sourceparams' ) then
                set mijtmp = `echo $mij_sourceparams`
                set mijtmp[$map_mij[$i]] = '1.E20'
                echo $mijtmp > $srcfile.$isrc.$isim
                echo $srctype1[$i] >> $srcfile.$isrc.$isim
                echo $srctype2[$i] >> $srcfile.$isrc.$isim
                tail -n 12 $homedir/$srcfile >> $srcfile.$isrc.$isim
            
            else if ( $src_file_type == 'cmtsolut' ) then
                set num = `echo $num $map_mij[$i] |awk '{print $1+$2}'`
                head -n $num $homedir/$srcfile.MIJ > $srcfile.$isrc.$isim
                echo $map_mij_name[$i] "    1.0000e+27" >> $srcfile.$isrc.$isim
                if ( $isim == 'MXX_P_MYY' ) then 
                    echo "Mpp:     1.0000e+27" >> $srcfile.$isrc.$isim       
                    @ num++
                endif
                set num = `echo 12 $num |awk '{print $1-$2}'`
                tail -n $num $homedir/$srcfile.MIJ >> $srcfile.$isrc.$isim
            endif 

            else if ( $src_file_type == 'finfault' ) then	    
                head -n 1 $homedir/$srcfile >! $srcfile.$isrc.$isim
                echo "1              fflt_num : number of individual points below" >> $srcfile.$isrc.$isim
                head -n 9 $homedir/$srcfile |tail -n 7 >> $srcfile.$isrc.$isim
                tail -n $num_src $homedir/$srcfile |head -n $i |tail -n 1 >> $srcfile.$isrc.$isim
                echo $srctype1[$i] >> $srcfile.$isrc.$isim
                echo $srctype2[$i] >> $srcfile.$isrc.$isim
            endif

        endif 
        
        if ( $multisrc == 'false' ) then
            set simdir = './'
        else 
            if ( $num_src == 1 ) then
                set simdir = $isim
                mkdir $simdir
                cd $simdir
            else 
                set simdir = $isrc"_"$isim
                mkdir $simdir
                cd $simdir
            endif
        endif 
        
        if ( -d $datapath) then
            echo " Saving data into $datapath"
        else
            echo "creating $datapath" 
            mkdir $datapath
        endif
        
        if ( -d $infopath) then 
            echo " saving info into $infopath"
        else
            echo "creating $infopath"
            mkdir $infopath
        endif
        
        mkdir Code
        cp -p $homedir/*.f90 Code
        cp -p $homedir/*.F90 Code
        cp -p $homedir/Makefile Code
        
        echo "copying crucial files for the simulation..."
        
        if ( $multisrc == 'true' ) then
            mv ../$srcfile.$isrc.$isim $srcfile
        else 
            cp $homedir/$srcfile $srcfile
        endif
        
        cp $homedir/xsem .
        cp $homedir/mesh_params.h .
        cp $homedir/runinfo .
        cp $homedir/$recfile . 
        cp $homedir/inparam_basic .
        cp $homedir/inparam_advanced .
        cp $homedir/inparam_hetero .
        cp $homedir/inparam_xdmf .

        if ( $multisrc == 'false' ) then
            ln -s ../$meshdir/ Mesh
        else 
            ln -s ../../$meshdir/ Mesh
        endif
        
        if ( $bgmodel == 'external' ) then
            cp Mesh/external_model.bm .
        endif
        cd $mainrundir

        cp $homedir/mesh_params.h .
        cp $homedir/inparam_basic .
        cp $homedir/inparam_advanced .
        cp $homedir/inparam_hetero .
        cp $homedir/inparam_xdmf .
    end
end



########################################################
######### submit the jobs ##############################
########################################################

set nodnum = `grep nproc_mesh $homedir/mesh_params.h |awk '{print $6}'`
echo "preparing job on $nodnum nodes..."

foreach isrc (${num_src_arr})
    foreach isim  (${srcapp})
        if ( $num_src == 1) then 
            cd $isim
        else 
            cd $isrc"_"$isim
        endif

        if ( $multisrc == 'true' ) then
            set outputname = "OUTPUT_"`echo $isim |sed 's/\//_/g'`
        else
            set outputname = "OUTPUT_"`echo $1 |sed 's/\//_/g'`
        endif

        if ( $newqueue == 'true' ) then 

            ########## LSF SCHEDULER ######################
            if ( $queue == 'lsf' ) then 
                bsub -R "rusage[mem=2048]" -I -n $nodnum mpirun -n $nodnum ./xsem 2>&1 > $outputname &

            ######## TORQUE/MAUI SCHEDULER #######
            else if ( $queue == 'torque' ) then 
                echo "# Sample PBS for parallel jobs" > run_solver.pbs
                echo "#PBS -l nodes=$nodnum,walltime=7:59:00" >> run_solver.pbs
                echo "ulimit -s unlimited " >> run_solver.pbs
                echo "cd $PWD " >> run_mesh.pbs
                echo "mpirun -n $nodnum ./xsem  > OUTPUT " >> run_solver.pbs
                qsub run_solver.pbs
            endif

        ######## SUBMIT LOCALLY #######
        else 
            mpirun.openmpi -n $nodnum ./xsem >& $outputname &
        endif

        echo "Job running in directory $isim"
        cd $mainrundir
    end
end 

######## post processing ##################################################

set F90 = `grep "FC = " $homedir/Makefile |awk '{print $3}'`

cd $homedir/UTILS
cp $homedir/mesh_params.h .
$F90 post_processing.f90 -o xpost_processing

cd $homedir
cd $1
cp -p $homedir/UTILS/xpost_processing .
cp -p $homedir/UTILS/post_processing.f90 .
cp -p $homedir/UTILS/post_processing.csh .
cp -p $homedir/UTILS/plot_recfile_seis.csh .
cp -p $homedir/UTILS/plot_recs.plot .
cp -p $homedir/UTILS/taup_allrec.csh .
cp -p $homedir/UTILS/plot_record_section.m .

echo "To convolve and sum seismograms, run ./post_processing.csh after the simulations in:" 
echo $mainrundir
echo ".... the post-processing input files param_post_processing and param_snaps are generated in the solver"
echo ".... based on guesses. Edit please."
echo " ~ ~ ~ ~ ~ ~ ~ h a n g   o n   &   l o o s e ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"


