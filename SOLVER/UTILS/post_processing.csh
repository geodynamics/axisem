#!/bin/csh -f

echo ""
echo ""
echo ">>>>>>>>> AXISEM POST PROCESSING <<<<<<<<<"
echo ""
echo ""

if ( ${#argv} < 4 ) then
    echo "FYI -- Optional arguments:"
    echo " EITHER: Argument 1: receiver components: enz,sph,cyl,xyz,src"
    echo " Argument 2: convolution period: 0. if not to be convolved"
    echo " Argument 3: source time function (gauss_0, gauss_1, quheavi)"
    echo " Argument 4: seismogram type (disp, velo)"
    echo " OR: Argument 1: Mij"
    echo " Argument 2-7: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp"
    echo
    if ( $1 == '-h' ) then
        echo " Input for post-processing is retrieved from 4 input files:"
        echo "   1) <simdir>/param_post_processing : contains everything related to summation, convolution, receiver components"
        echo "   2) param_sum_seis : number of simulations (e.g. if moment, then 4)"
        echo "   3) <simdir>/param_snaps : if wavefield snaps were computed, this file defines the 3D geometry"
        echo "   4) <simdir>/simulation.info: parameters from the simulation"
        echo
        echo "  <simdir> = ./  if nsim=1 (i.e. only one simulation)"
        echo "  <simdir> = MZZ, etc if nsim=4 via 'moment' as second argument in submitting simulation"
        echo "  <simdir> = P_Z, etc if nsim=2 via 'force' as second argument in submitting simulation"
        echo " Providing arguments will overwrite the settings from the input files."
        echo
        exit
    endif
endif

if ( $1 == 'Mij' ) then 
    if ( ${#argv} != 7 ) then 
        echo "wrong number of arguments for Mij!"
        exit
    endif
    set Mij_overwrite = ( $2 $3 $4 $5 $6 $7 ) 
    echo "Overwritten Mij:" $Mij_overwrite
endif

set homedir = $PWD
set simdir1 = `/usr/bin/tail -n 1 param_sum_seis |awk '{print $1}' |sed 's/"//g' `

if ( ! -f param_post_processing ) then
#    set outdir = `/usr/bin/tail  -n 1 $simdir1/param_post_processing |awk '{print $1}' |sed "s/'/ /g" `
    set outdir = `grep "Directory for" $simdir1/param_post_processing |awk '{print $1}' |sed "s/'/ /g" `
else
#    set outdir = `/usr/bin/tail  -n 1 param_post_processing |awk '{print $1}' |sed "s/'/ /g" `
    set outdir = `grep "Directory for" param_post_processing |awk '{print $1}' |sed "s/'/ /g" `
endif
echo "All post-processed data in: "$outdir

rm -f $outdir/param_mij

set nsim = `head -n 1 param_sum_seis |awk '{print $1}'`
if ( $nsim == 1 ) then 
    set simlist = "./"
else
    set simlist = `/usr/bin/tail -n $nsim param_sum_seis |awk '{print $1}' |sed 's/"//g' ` 
    if ( $nsim == 4 ) then
        if ( ! -f param_snaps ) then 
            cp $simdir1/param_snaps .
        endif
    endif
endif

echo "number of simulations: $nsim"; echo ""

if ( ! -d $outdir ) then
    mkdir $outdir
    if ( -f param_snaps ) then
        mkdir $outdir/SNAPS
    endif
endif

foreach isimdir (${simlist})
    echo "working in simulation $isimdir "
    if ( ! -d $isimdir/$outdir ) then
        mkdir $isimdir/$outdir
        if ( -f param_snaps && $nsim > 1 ) then
           mkdir $isimdir/$outdir/SNAPS
        endif
    endif
    cp -p param_* $outdir 
    /bin/cp -p -f param_* $outdir
    cp -p post_processing.csh $outdir
    cp -p xpost_processing $outdir
    cp -p post_processing.f90 $outdir 
end
    cp -p plot_record_section.m $outdir

if (! -f mesh_params.h) then 
    cp -p $simlist[1]/mesh_params.h .
endif

echo
echo "%%%%%%%%% PROCESSING seismograms/wavefields %%%%%%%%%"
mkdir $outdir/SEISMOGRAMS
mkdir $outdir/SEISMOGRAMS/UNPROCESSED
if ( ${#argv} > 0 && $1 != 'Mij' )  then
    echo ${#argv} > $outdir/param_post_processing_overwrite
    echo $1 >> $outdir/param_post_processing_overwrite
    if ( ${#argv} > 1 ) then
      echo $2>> $outdir/param_post_processing_overwrite
    endif
    if ( ${#argv} > 2 ) then
      echo $3>> $outdir/param_post_processing_overwrite
    endif
    if ( ${#argv} > 3 ) then
      echo $4>> $outdir/param_post_processing_overwrite
    endif
else if ( $1 == 'Mij' ) then 
    echo $Mij_overwrite > $outdir/param_mij
endif

echo ".... output in "$outdir"/OUTPUT_postprocessing ...."
./xpost_processing > $outdir/OUTPUT_postprocessing
echo "Done with post processing, results in SEISMOGRAMS/ " 

if ( -f param_snaps) then 
    echo " .... and SNAPS/"
endif

set gnu_query = `which gnuplot | wc -l `
if ( $gnu_query == 1 ) then
    echo
    echo "%%%%%%%%% PLOTTING seismograms (gnuplot) %%%%%%%%%%"
    cd $outdir
    set seistype = `tail -n 3 $homedir/$isimdir/param_post_processing |head -n 1 | sed "s/'/ /g" | awk '{print $1}'`
    echo "seismogram type:" $seistype
    set reclist = `cat $homedir/$isimdir/Data/receiver_names.dat |awk '{print $1}'`
    echo "1st receiver:" $reclist[1]
    set colat = `cat $homedir/$isimdir/Data/receiver_names.dat |awk '{print $2}' |sed 's/00000/ /g' |awk '{print $1}'`
    set lon = `cat $homedir/$isimdir/Data/receiver_names.dat |awk '{print $3}' |sed 's/00000/ /g'  |awk '{print $1}'`
    set epidist = `cat $homedir/$isimdir/Data/receiver_pts.dat |awk '{print $1}'`
    echo "1st receiver colatitude/longitude/epidist:" $colat[1] " " $lon[1] " " $epidist[1]

    set reccomp = `ls SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*.dat |sed 's/mij_/ /g ' |awk '{print $2}' | sed 's/_/ /g' |awk '{print $2}'  |sed 's/\.dat/ /g '`

    #set reccomp = `ls SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*.dat |sed 's/mij_/ /g ' |awk '{print $2}' | sed 's/_/ /g' |awk '{print $2}'  |sed 's/\.dat/ /g ' |awk '{print $1}'`
    set conv = `ls SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*.dat |sed 's/_mij_/ /g '  | sed 's/\.dat/ /g ' |awk '{print $2}' `
    set t2 = `tail -n 1 SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*{$reccomp[1]}.dat |awk '{print $1}' `
    echo "convolution:" $conv
    echo "receiver components:" $reccomp
    mkdir GRAPHICS
    set i = 0
    foreach  rec (${reclist})
    @ i++
    set j = 0
    foreach comp (${reccomp}) 
    @ j++
        set recname = `echo $rec"_"$seistype"_post_mij_"{$conv[$j]}`
        echo "Plotting receiver " $recname 
        echo 'set term png linewidth 1  ' >! plot_recs.plot
        echo 'set output "GRAPHICS/'$recname'.png"' >> plot_recs.plot
        echo 'set title "colat,lon: '$colat[$i], $lon[$i]', epidist: '$epidist[$i]'"'>> plot_recs.plot
        echo 'plot "SEISMOGRAMS/'$recname'.dat" with lines' >> plot_recs.plot
        echo "set xrange [ 0: "$t2"];set xlabel 'time [s]';set ylabel 'displacement [m]' " >> plot_recs.plot
        gnuplot plot_recs.plot
        cd GRAPHICS; 
        convert $recname.png $recname.gif
        convert $recname.png $recname.pdf
        rm -f $recname.png
        cd ..
    end
    end
    echo "Done with plotting, results in GRAPHICS/"
endif

cd $homedir
set taup_query = `which taup_time | wc -l `
if ( $taup_query == 1 ) then
    echo
    echo "%%%%%%%%% Computing TRAVELTIMES (taup) %%%%%%%%%"
    set depth = `grep "source depth" $simdir1/simulation.info |awk '{print $1}'`
    set model = `grep Background $simdir1/mesh_params.h |awk '{print $5}'`

    if ( $model == 'prem' || $model == 'iasp91' ) then
	set num_rec = `wc -l $simdir1/Data/receiver_pts.dat |awk '{print $1}'`
#	set epi_list = `tail -n $num_rec $simdir1/Data/receiver_pts.dat | sed 's/999999/9/g' |  sed 's/000000/ /g' | awk '{print $1}'`
	set epi_list = `tail -n $num_rec $simdir1/Data/receiver_pts.dat | awk '{print $1}'`
	set depth_short = `echo $depth |sed 's/\./ /g' |awk '{print $1}'` 
	echo "Earthquake depth:" $depth_short
	set i = 0
	cd $outdir
	foreach rec (${epi_list})
	    @ i++
	    echo "traveltimes for epicentral distance" $rec

	    set tt =  `taup_time -mod $model -h $depth -ph pP -deg $rec | grep " pP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_pP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph PP -deg $rec | grep " PP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_PP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph SS -deg $rec | grep " SS " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_SS_traveltime2.dat ;	    endif

	    if ( $rec < 100. ) then
	    set tt =  `taup_time -mod $model -h $depth -ph P -deg $rec | grep " P " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_P_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph S -deg $rec | grep " S " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_S_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph PcP -deg $rec | grep " PcP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_PcP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph ScS -deg $rec | grep " ScS " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_ScS_traveltime2.dat ;	    endif
	    endif

	    if ( $rec > 95. ) then
	    set tt =  `taup_time -mod $model -h $depth -ph Pdiff -deg $rec | grep " Pdiff " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_Pdiff_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph Sdiff -deg $rec | grep " Sdiff " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_Sdiff_traveltime2.dat ;	    endif

	    endif

	end
	sort -n taup_P_traveltime2.dat |grep -v "==>" > taup_P_traveltime.dat
	sort -n taup_pP_traveltime2.dat |grep -v "==>" > taup_pP_traveltime.dat
	sort -n taup_S_traveltime2.dat |grep -v "==>" > taup_S_traveltime.dat
	sort -n taup_PP_traveltime2.dat |grep -v "==>" > taup_PP_traveltime.dat
	sort -n taup_SS_traveltime2.dat |grep -v "==>" > taup_SS_traveltime.dat
	sort -n taup_PcP_traveltime2.dat |grep -v "==>" > taup_PcP_traveltime.dat
	sort -n taup_ScS_traveltime2.dat |grep -v "==>" > taup_ScS_traveltime.dat
	sort -n taup_Pdiff_traveltime2.dat |grep -v "==>" > taup_Pdiff_traveltime.dat
	sort -n taup_Sdiff_traveltime2.dat |grep -v "==>" > taup_Sdiff_traveltime.dat
	mkdir TAUP
	rm -f taup_*2.dat
	mv taup_*.dat TAUP
	echo "Done with taup, results in TAUP/"
    endif
endif

echo
echo
echo " %%%%%%%%%%% QUICK OVERVIEW OF RESULTS %%%%%%%%%%%"
echo "1) Check earthquake and seismograms: "
echo "        open google earth"
echo "        load file "$outdir"/googleearth_src_rec_seis.kml"
echo "        locate earthquake, click and check parameters"
echo "        click on receiver location, check location and seismograms"
echo
echo " 2) Seismogram record section:"
echo "        open matlab in directory" $outdir
echo "        >> plot_record_section"
echo "        plots of seismograms sorted in epicentral distance with relative amplitudes"
echo "        station name given on the left of each trace"
echo "        traveltimes for some phases (if simulations are upon PREM or IASP91)"
echo
echo " 3) Seismogram time series:"
echo "        load individual seismograms, e.g. with the command:"
echo "        >> xmgrace" $outdir"/SEISMOGRAMS/"$recname".dat"

if ( -f param_snaps ) then
    echo
    echo " 4) Wavefield snapshots:"
    echo "        open paraview"
    echo "        load "$outdir"/SNAPS/snap_mij_cell_*vtk"
    echo "        ... this should be a wavefield movie in 3D. "
    echo "        adjust normalization to e.g. [-0.1,0.1]"
endif

echo
echo "==========================================="
echo "                  DONE with post processing."
echo "==========================================="
