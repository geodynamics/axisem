# AxiSEM 1.4 [![Build Status](https://travis-ci.org/geodynamics/axisem.svg?branch=master)](https://travis-ci.org/geodynamics/axisem)

## Axially symmetric Spectral Element Method

Copyright 2018, Tarje Nissen-Meyer, Martin van Driel, Simon Stähler, Kasra Hosseini, Stefanie Hempel, Lion Krischer, Alexandre Fournier

Webpage and distribution: http://www.axisem.info
Contact and information:  info@axisem.info

April, 13, 2018 

## Citation
If you are publishing results obtained with this code, please cite this paper:

T. Nissen-Meyer, M. van Driel, S. C. Staehler, K. Hosseini, S. Hempel, L. Auer, A. Colombi and A. Fournier:
**"AxiSEM: broadband 3-D seismic wavefields in axisymmetric media"**, *Solid Earth*, 5, 425-445, 2014
doi:10.5194/se-5-425-2014 http://www.solid-earth.net/5/425/2014/

## Content of the repository
`manual_axisem_1.4.pdf` - PDF manual

`MESHER` - The program to generate 2D meshes for the SEM forward solver

`SOLVER` - the SEM forward solver itself

`submit.py` - a Python script to create an Instaseis database

`make_axisem.macros` - macro file to set compiler options

`copytemplates.sh` - reset all input files to default templates 

`COPYING` - The GNU General Public License

`HISTORY` - changelog

`README` - this file

## Basic instructions for running:

More details are found in the manual. For a quick start:

 - The mesher has to be run first to generate a SEM mesh for the solver. 
 - Any changes in the resolution, spherically symmetric background model, or number 
   of processors requires a new mesh. 
 - Changes in the source-receiver settings or 2.5D heterogeneities only need a new solver run.
 - Changes on the moment tensor, receiver components, or filtering only need a new postprocessing run.
 - Using Instaseis, seismograms for any depth or moment tensor can be calculated from the wavefield of one force source at the surface.

General settings and explanations for parameters needed in MESHER and SOLVER 
are found in the `inparam_*` files in the respective directories. 

To create Instaseis databases quickly, run the script submit.py in the main directory.

1) Run `copytemplates.sh` to set up a generic run with pre-set parameters

2) Go into the MESHER directory, run `./submit.csh`. This compiles the code using
gfortran and mpif90 as default compilers, and then submits a job on a single node. 
For high resolution (seismic period below 3s), this requires significant amounts 
of RAM, see manual.

3) Check `OUTPUT`; if finished, then run `./movemesh.csh <MESH_NAME>`

4) go into `SOLVER`, check the vtk files in `/MESHES/<MESH_NAME>` with Paraview, if you want

5) Edit `inparam_basic` to the desired `<MESH_NAME>`

6) Run `./submit.csh <RUN_NAME>` , this compiles and then submits a parallel run.

7) Check `<RUN_NAME>/OUTPUT*`, which will record the progress of the run and any problems.

8) All data-related output is in `<RUN_NAME>/Data/`

9) Convolution with a STF and summation for a moment source is done by running
   `postprocessing.csh` in `<RUN_NAME>`

10) A more modern and efficient way of seismogram retrieval is using Instaseis,
    a Python toolbox to retrieve seismograms for arbitrary depths and moment
    tensors from the stored wavefield of one AxiSEM run (http://www.instaseis.net)

Detailed instructions can be found in the file `manual_axisem_1.4.pdf`

## Installation
Generally, AxiSEM needs
- a C and a Fortran90 compiler (typically gcc and gfortran)
- MPI (typically OpenMPI)
- NetCDF (not absolutely necessary, but needed to create an Instaseis database)

### MacOS X
The easiest way to install the necessary requirements is *homebrew* (https://brew.sh/). If you have not installed homebrew itself yet, do so as described on the homepage. The preferred way (as of May 2022) is 
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
and follow the steps therein.

Afterwards, install gfortran, openmpi and NetCDF with 
```
brew install gfortran gcc openmpi netcdf hdf5-mpi
```

### Linux
The necessary packages can easily be installed using APT
```
sudo apt install gfortran gcc libopenmpi-dev openmpi-bin libnetcdff-dev 
```
