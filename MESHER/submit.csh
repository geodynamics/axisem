#!/bin/csh -f

if ( $1 == '-h' ) then 
  echo "Argument options:"
  echo "   default (no argument): submit xmesh on current machine"
  echo "   lsf: submit to a lsf queue using bsub"
  echo "   torque: submit to a Torque/Mauo queue using qsub"
  exit
endif

# tidy up:
rm -rf OUTPUT
rm -rf meshdb.dat0*
rm -rf mesh_params.h*

if ( ! -d Diags ) then
  mkdir Diags
else
  rm -rf Diags/*
endif

if ( { make all } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif

if ( $1 == 'lsf' ) then 
  ########## LSF SCHEDULER ######################
  bsub -R "rusage[mem=2048]" -I -n 1 ./xmesh > OUTPUT &

else 
  if ( $1 == 'torque' ) then 
    ######## TORQUE/MAUI SCHEDULER #######
    echo "# Sample PBS for serial jobs" > run_mesh.pbs
    echo "#PBS -l nodes=1,walltime=2:00:00" >> run_mesh.pbs
    echo "ulimit -s unlimited " >> run_mesh.pbs
    echo "cd $PWD " >> run_mesh.pbs
    echo "./xmesh > OUTPUT " >> run_mesh.pbs
    qsub run_mesh.pbs

  else
    ######## SUBMIT LOCALLY #######
    nohup ./xmesh > OUTPUT &
  endif
endif

echo 'xmesh submitted, output in "OUTPUT"'
echo "After the run, move the mesh to a new directory <meshdir> via:"
echo "./movemesh <meshdir>"
echo "This will be located in ../SOLVER/MESHES/<meshdir>"

