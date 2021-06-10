## Changes to .bashrc

To be able to use parallel NetCDF 4, you need to load the latest LMOD 
modules. To do this, add the following lines to your ~/.bashrc:
```
# User specific aliases and functions
env2lmod
module load gcc/8.2.0
module load openmpi/4.0.2
module load netcdf-fortran
module load python/3.6.5
```

## File storage
Also, you should make sure that the data files are stored on scratch. 
For that, it is easiest to create a run directory in $SCRATCH and 
place a symlink to your axisem parent directory
```
cd axisem
rmdir ./runs
mkdir $SCRATCH/runs
ln -s $SCRATCH/runs .
```

## Run jobs
This can be done using the submit.py script:
```
python ./submit.py Mercury_10s ~/axisem/MESHER/Models/Mercury_Rivoldini_cold_handpicked.bm 10.0 --nrad 2 --ntheta 0 -j Euler 
```
make sure that you have an `inparam_basic` and `inparam_advanced` in your axisem 
parent directory.
