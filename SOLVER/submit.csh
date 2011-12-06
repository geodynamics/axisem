#!/bin/csh -f

set homedir = $PWD

if ( ${#argv} < 1 || "$1" == "-h" ) then
    echo "";echo "=================================================="
    echo " Argument 1:  directory name for the simulation"; echo ""
    echo " Optional arguments after directory name: "
    echo " -s <sourcetype>, where <sourcetype> can be:"
    echo "      'moment': submits 4 simulations for full moment tensor"
    echo "      'force': submits 2 simulations for all single forces"
    echo "       Default: 1 simulation for the source specified in input file"     
    echo " -q <queue>, where <queue> can be:"
    echo "      'lsf': submit to lsf queue using bsub"
    echo "       Default: submit locally"
    echo "==================================================";echo"";    exit
else if ( -d $1) then
    echo " Run or directory" $1 "exists....... its content:"; ls $1; exit
endif

set datapath = `grep "data output path" inparam  |awk '{print $1}'`
set infopath = `grep "info output path" inparam |awk '{print $1}'`

# Check arguments: source types and submission queues
set multisrc = 'false'
set newqueue = 'false'

if ( "$2" == "-s" ) then
    set srctype = `echo $3`
    set multisrc = 'true'
    if ( "$4" == '-q' ) then 
        set queue = `echo $5`
        set newqueue = 'true'
    endif
else if ( "$2" == '-q') then
    set queue = `echo $3`
    set newqueue = 'true'
    if ( "$4" == '-s' ) then 
        set srctype = `echo $5`
        set multisrc = 'true'
    endif
endif

if ( $multisrc == 'true' ) then 
    echo "Submitting multiple jobs for full" $srctype
endif

if ( $newqueue == 'true' ) then 
	echo "Submitting to queue type" $queue
endif

# recompile if mesh_params.h is younger than xsem
if ( ! -f xsem) then 
    echo "==============================================="
    echo "Need to compile since executable xsem not existant"
    echo "==============================================="
    make clean; make
    if ( ! -f xsem ) then 
        echo "Compilation error, please check the code."; exit
    endif
else
    set younger = `ls -lt mesh_params.h xsem |awk '{print $8}'`
    if ( $younger[1] == 'mesh_params.h'  ) then 
        echo "Need to recompile since mesh is younger than executable!"; make clean; make
        if ( ! -f xsem ) then 
          echo "Compilation error, please check the code."; exit
        endif
    endif
endif

# identify source input file 
set src_file_type = `grep "source file type" inparam |awk '{print $1}'`
echo "Source file type:" $src_file_type

if ( $src_file_type == 'cmtsolut' ) then
    set srcfile = 'CMTSOLUTION'
else if ( $src_file_type == 'separate' ) then
    set srcfile = 'sourceparams.dat'
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
set rec_file_type = `grep "receiver file type" inparam |awk '{print $1}'`
echo "Receiver file type:" $rec_file_type

if ( $rec_file_type == 'colatlon' ) then
    set recfile = 'receivers.dat'
else if ( $rec_file_type == 'stations' ) then
    set recfile = 'STATIONS'
endif

if ( ! -f $homedir/$recfile ) then 
    echo "file $recfile does not exist"
    exit
endif
echo "Source file:" $srcfile, "Receiver file:" $recfile

# multiple simulations
set num_src = 1
set num_src_arr = ( 1 )
if ( $multisrc == 'true' ) then
    echo " setting up multiple simulations for full" $srctype "source type"
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
        echo " Choose either 'moment', 'force', 'finfault' or leave blank for one simulation as in sourceparams.dat"; exit
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

# write parameter file for summing seismograms after the simulations
echo $numsim $num_src |awk '{print ( $1 * $2)}'"             number of simulations" > param_sum_seis

# Prepare and copy relevant files for each simulation
foreach isrc (${num_src_arr})
    set i = 0
    ########################
    foreach isim  (${srcapp})
    ########################

        @ i ++

        set num = 6
        echo ""
        echo "Setting up simulation" $isrc,$isim
        # construct different source file for each simulation
        if  ( $multisrc == 'true' ) then
            echo "constructing separate source files for" $isim 

            if ( $src_file_type == 'separate' ) then
        #	    head -n 1 $homedir/$srcfile  > $srcfile.$isrc.$isim
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
                mkdir $simdir; cd $simdir
            else 
                set simdir = $isrc"_"$isim
                mkdir $simdir; cd $simdir
            endif
        endif 
        
        if ( -d $datapath) then
            echo " Saving data into $datapath"
        else
            mkdir $datapath; echo "creating $datapath" 
        endif
        
        if ( -d $infopath) then 
            echo " saving info into $infopath"
        else
            mkdir $infopath; echo "creating $infopath"
        endif
        
        mkdir Code
        cp -p $homedir/*.f90 Code
        cp -p $homedir/Makefile Code
        
        echo "copying crucial files for the simulation..."
        
        if ( $multisrc == 'true' ) then
            cp ../$srcfile.$isrc.$isim $srcfile
        else 
            cp $homedir/$srcfile $srcfile
        endif
        
        cp $homedir/xsem .
        cp $homedir/mesh_params.h .
        cp $homedir/$recfile . 
        cp $homedir/mpif.h .
        cp $homedir/inparam .
        cp $homedir/inparam_hetero .
        cp $homedir/*.bm .

        cd $mainrundir
        
        # write parameter file for summing seismograms after the simulations
        echo '"'$simdir'"' >> param_sum_seis

    ########################
    end
########################
end

########################################################
######### submit the jobs ##############################
########################################################

set nodnum = `grep nproc_mesh $homedir/mesh_params.h |awk '{print $6}'`
echo "preparing job on $nodnum nodes..."

foreach isrc (${num_src_arr})
#-----------------------
    foreach isim  (${srcapp})
    #-----------------------
        if ( $num_src == 1) then 
            cd $isim
        else 
            cd $isrc"_"$isim
        endif

        if ( $numsim > 1 ) then 
            cp $mainrundir/param_sum_seis .
        endif

        if ( $multisrc == 'true' ) then
            set outputname = "OUTPUT_"`echo $isim |sed 's/\//_/g'`
        else
            set outputname = "OUTPUT_"`echo $1 |sed 's/\//_/g'`
        endif

        if ( $newqueue == 'true' ) then 

    ########## LSF SCHEDULER ######################
            if ( $queue == 'lsf' ) then 
                bsub -R "rusage[mem=2048]" -I -n $nodnum mpirun -n $nodnum ./xsem > $outputname &

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
            mpirun -n $nodnum ./xsem > $outputname &
        endif

        echo "Job running in directory $isim"
        cd $mainrundir
    #-----------------------
    end
#-----------------------
end 

######## post processing ##################################################

set F90 = `grep "FC = " $homedir/Makefile |awk '{print $3}'`
cd $homedir/UTILS
cp $homedir/mesh_params.h .
$F90 post_processing.f90 -o xpost_processing
cp -p xpost_processing $homedir/$1
cp -p post_processing.f90 $homedir/$1
cp -p post_processing.csh $homedir/$1
cp -p plot_recfile_seis.csh $homedir/$1
cp -p plot_recs.plot $homedir/$1
cp -p taup_allrec.csh $homedir/$1
cp -p plot_record_section.m $homedir/$1
#echo ""; cd $homedir; make sum_seis; cp sum_seis.x sum_seis.f90 $mainrundir
echo "To convolve and sum seismograms, run ./post_processing.csh after the simulations in:" 
echo $mainrundir
echo ".... the post-processing input files param_post_processing and param_snaps are generated in the solver"
echo ".... based on guesses. Edit please. Input file param_sum_seis is generated with this submit script"
echo "  The moment tensor is given/assumed in the following order:"
echo "Mzz,Mxx,Myy,Mxz,Myz,Mxy / Mrr,Mtt,Mpp,Mtr,Mpr,Mtp"
echo " ~ ~ ~ ~ ~ ~ ~ h a n g   o n   &   l o o s e ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"

