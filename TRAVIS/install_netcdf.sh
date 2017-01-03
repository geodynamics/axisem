#!/bin/bash

TRAVIS_ROOT="$1"
GCCVERSION="$2"

export CC=gcc-$GCCVERSION
export FC=gfortran-$GCCVERSION
export TRAVIS_ROOT=$TRAVIS_ROOT
export LDFLAGS=-L$TRAVIS_ROOT/lib
export CPPFLAGS="-I$TRAVIS_ROOT/include -DgFortran"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TRAVIS_ROOT/lib 
set -e


if [ ! -f "$TRAVIS_ROOT/bin/h5stat" ]; then
    ##zlib 1.2.8
    rm -rf zlib-1.2.8 zlib-1.2.8.tar.gz
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.8.tar.gz
    tar -xvf zlib-1.2.8.tar.gz
    cd zlib-1.2.8
    ./configure --prefix=$TRAVIS_ROOT
    make check
    make install
    cd ..
    #
    #
    ##HDF 1.8.17
    rm -rf hdf5-*
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.17/src/hdf5-1.8.17.tar.bz2
    tar -xf hdf5-1.8.17.tar.bz2
    cd hdf5-1.8.17/

    ./configure --prefix=$TRAVIS_ROOT --with-zlib=$TRAVIS_ROOT

    # -j parallelizes make;  -s reduces output
    make -sj4
    #make check 
    make install
    cd ..

else
    echo "HDF5 installed..."
    find $TRAVIS_ROOT -name h5stat
fi


if [ ! -f "$TRAVIS_ROOT/bin/nc-config" ]; then
    # netcdf-4.3.3.1
    rm -rf netcdf-4.*
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
    tar -xvf netcdf-4.3.3.1.tar.gz
    cd netcdf-4.3.3.1
    ./configure --enable-netcdf-4 --enable-dap --enable-shared --prefix=$TRAVIS_ROOT
    make -sj4
    #make check # all tests should succeed.
    make install
    cd ..
else
    echo "HDF5 installed..."
    find $TRAVIS_ROOT -name nc-config
fi


if [ ! -f "$TRAVIS_ROOT/bin/nf-config" ]; then
    #netcdf-fortran-4.4$
    rm -rf netcdf-fortran-4.*
    export LDFLAGS="-L$TRAVIS_ROOT/lib -lnetcdf "
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.2.tar.gz
    tar -xvf netcdf-fortran-4.4.2.tar.gz
    cd netcdf-fortran-4.4.2
    ./configure --prefix=$TRAVIS_ROOT --enable-netcdf-4
    make -sj4
    #make check
    make install
    cd ..
else
    echo "NetCDF-fortran installed..."
    find $TRAVIS_ROOT -name nf-config
fi


# echo 'You might want to put these to lines in your Makefile to make sure, the'
# echo 'libraries just compiled are used:'
# echo ''
# echo 'LIBS = -lm  -L $(TRAVIS_ROOT)/lib -lnetcdff -Wl,-rpath,$(TRAVIS_ROOT)/lib'
# echo 'INCLUDE = -I $(TRAVIS_ROOT)/include'

