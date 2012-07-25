#!/bin/csh -f

if ( $1 == '-h' ) then 
  echo "Argument options:"; echo "default (no argument): submit xmesh on current machine"
  echo "lsf: submit to a lsf queue using bsub"
  echo "torque: submit to a Torque/Mauo queue using qsub";  exit
endif

if ( { make all } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif

if ( ! -d Diags ) then
  mkdir Diags
endif

########## LSF SCHEDULER ######################
  if ( $1 == 'lsf' ) then 
    bsub -R "rusage[mem=2048]" -I -n 1 ./xmesh > OUTPUT &

######## TORQUE/MAUI SCHEDULER #######
  else if ( $1 == 'torque' ) then 
    echo "# Sample PBS for serial jobs" > run_mesh.pbs
    echo "#PBS -l nodes=1,walltime=2:00:00" >> run_mesh.pbs
    echo "ulimit -s unlimited " >> run_mesh.pbs
    echo "cd $PWD " >> run_mesh.pbs
    echo "./xmesh > OUTPUT " >> run_mesh.pbs
    qsub run_mesh.pbs

######## SUBMIT LOCALLY #######
  else
    nohup ./xmesh > OUTPUT &
  endif

  echo 'xmesh submitted, output in "OUTPUT"'
  echo "After the run, move the mesh to a new directory <meshdir> via:"
  echo "./copymesh <meshdir>"
  echo "This will be located in ../SOLVER/MESHES/<meshdir>"

endif
