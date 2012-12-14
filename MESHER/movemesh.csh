#!/bin/csh -f

# test if mesher finished without problems:
if (`tail -n 1 OUTPUT |grep 'DONE WITH MESHER'|wc -l` != '1') then
  echo "ERROR: MESHER did not finish yet or encountered a problem."
  echo "       Check 'OUTPUT' file for more details."
  exit
else
  echo 'MESHER finished smoothly...'
endif

set homepath = `echo $PWD`
set meshpath = "../SOLVER/MESHES/$1"

if ( ! -d ../SOLVER/MESHES ) then
  mkdir ../SOLVER/MESHES
endif

if ( ! -d $meshpath ) then
  mkdir $meshpath
else
  echo "ERROR: the mesh folder " $meshpath " already exists!"
  exit
endif

# convert relative to absolute path:
cd $meshpath
set meshpath = `echo $PWD`
cd $homepath

echo "character(len=200) :: meshpath='$meshpath'" >> mesh_params.h

echo "Moving mesh to" $meshpath

mv meshdb.dat0* $meshpath
mv mesh_params.h* $meshpath
mv OUTPUT $meshpath
cp -p inparam_mesh $meshpath
mv unrolled_loops.f90 $meshpath
cp -p background_models.f90 $meshpath
mkdir $meshpath/Code
cp -p *.f90 $meshpath/Code
cp -p Makefile $meshpath/Code
cp -p makemake.pl $meshpath/Code
cp -p inparam_mesh $meshpath/Code
cp -p xmesh $meshpath/Code
cp -p submit.csh $meshpath/Code

#mv Diags/serend_coords_per_proc.dat* $meshpath

mv Diags/* $meshpath

# we should really get rid of this head/tail stuff:
set dump_files = `head -n 12 inparam_mesh |tail -n 1 |awk '{print $1}'`
if ( $dump_files == '.true.') then
  cd UTILS 
  ./plot_meshes.csh
  cd ..

  mv Diags/grid.ps $meshpath
  mv Diags/grid_solid.ps $meshpath
  mv Diags/grid_fluid.ps $meshpath
  mv Diags/grid_central.ps $meshpath
  mv Diags/uppermantle_grid.ps $meshpath
  rm -f Diags/*.dat
endif 

echo "Contents in" $meshpath ":"
ls $meshpath
cd $meshpath
echo "Total size: `du -sh` "
echo "DONE."

#./plot_proc_meshes.csh
#cp -p grid_proc* $meshpath

#./plot_proc_valence.csh
#cp -p valence_proc*ps $meshpath
#cp -p valence_central_proc*ps $meshpath

#./plot_proc_messaging.csh
#cp -p message_send_proc*ps $meshpath
#cp -p message_recv_proc*ps $meshpath
                                                 
