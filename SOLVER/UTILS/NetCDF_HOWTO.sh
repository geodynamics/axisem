#!/bin/bash

export CC=gcc
export FC=gfortran
export LDFLAGS=-L$HOME/local/lib
export CPPFLAGS="-I$HOME/local/include -DgFortran"
export LD_LIBRARY_PATH=$HOME/local/lib 
set -e


#zlib 1.2.7
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.7.tar.gz
tar -xvf zlib-1.2.7.tar.gz
cd zlib-1.2.7
./configure --prefix=$HOME/local
make check
make install
cd ..


#HDF 1.8.9
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.9.tar.gz
tar -xvf hdf5-1.8.9.tar.gz
cd hdf5-1.8.9
# to use older version of gcc (4.6 seems to have a problem, see http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg02983.html):
#export OMPI_CC=gcc-4.5
#export CC=mpicc.openmpi
# there seems to be a problem with PIC, so if make check in netcdf tells you to recompile hdf5 with -fpic, uncomment this line:
#export CC="mpicc.openmpi -fPIC"
#export FC=mpif90.openmpi 

./configure --prefix=$HOME/local --with-zlib=$HOME/local 

# -j parallelizes make -s reduces output
make -sj

#check seems to freeze at test ph5diff. whatever...
make check 
make install
cd ..


# netcdf-4.2.1.1$
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.0.tar.gz 
tar -xvf netcdf-4.3.0.tar.gz
cd netcdf-4.3.0
./configure --enable-netcdf-4 --enable-dap --enable-shared --prefix=$HOME/local
make check # all tests should succeed.
make install
cd ..


#netcdf-fortran-4.2$
export LDFLAGS="-L$HOME/local/lib -lnetcdf "
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.2.tar.gz
tar -xvf netcdf-fortran-4.2.tar.gz
cd netcdf-fortran-4.2
 ./configure --prefix=$HOME/local 
# bug in the netcdf release. By default, the parallel tests are not enabled, whatever you set in configure. Replace the makefile in nf_test with the one attached to run the parallel tests.
make check
make install
cd ..


echo 'You might want to put these to lines in your Makefile to make sure, the'
echo 'libraries just compiled are used:'
echo ''
echo 'LIBS = -lm  -L $(HOME)/local/lib -lnetcdff -Wl,-rpath,$(HOME)/local/lib'
echo 'INCLUDE = -I $(HOME)/local/include'

