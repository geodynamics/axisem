#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Fortran

set -e
set -x

os=`uname`
TRAVIS_ROOT="$1"
MPI_IMPL="$2"
GCCVERSION="$3"

case "$os" in
    Darwin)
        echo "Mac"
        brew update
        case "$MPI_IMPL" in
            mpich)
                brew install mpich
                ;;
            openmpi)
                brew install openmpi
                ;;
            *)
                echo "Unknown MPI implementation: $MPI_IMPL"
                exit 10
                ;;
        esac
        ;;
    Linux)
        echo "Linux"
        export PRK_CC="gcc-$GCCVERSION"
        export PRK_CXX="g++-$GCCVERSION"
        export PRK_FC="gfortran-$GCCVERSION"

        case "$MPI_IMPL" in
            mpich)
                if [ ! -f "$TRAVIS_ROOT/bin/mpichversion" ]; then
                    set +e
                    wget --no-check-certificate -q \
                         http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
                    set -e
                    if [ ! -f "$TRAVIS_ROOT/mpich-3.2.tar.gz" ]; then
                        echo "MPICH download from mpich.org failed - trying Github mirror"
                        wget --no-check-certificate -q \
                             https://github.com/jeffhammond/mpich/archive/v3.2.tar.gz \
                             -O mpich-3.2.tar.gz
                        tar -xzf mpich-3.2.tar.gz
                        cd mpich-3.2
                    else
                        tar -xzf mpich-3.2.tar.gz
                        cd mpich-3.2
                    fi
                    sh $TRAVIS_HOME/travis/install-autotools.sh $TRAVIS_ROOT
                    ./autogen.sh
                    mkdir build ; cd build
                    ../configure --prefix=$TRAVIS_ROOT CC=$PRK_CC CXX=$PRK_CXX FC=$PRK_FC
                    make -sj4 
                    make install
                else
                    echo "MPICH installed..."
                    find $TRAVIS_ROOT -name mpiexec
                    find $TRAVIS_ROOT -name mpif90
                    find $TRAVIS_ROOT -name mpicc
                fi
                ;;
            openmpi)
                if [ ! -f "$TRAVIS_ROOT/bin/ompi_info" ]; then
                    wget --no-check-certificate -q http://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.1.tar.bz2
                    tar -xjf openmpi-1.10.1.tar.bz2
                    cd openmpi-1.10.1
                    mkdir build && cd build
                    ../configure --prefix=$TRAVIS_ROOT CC=$PRK_CC CXX=$PRK_CXX FC=$PRK_FC \
                      --disable-java --disable-cxx --without-verbs --without-fca --without-mxm \
                      --without-ucx --without-portals4 --without-psm --without-psm2 \
                      --without-alps --without-munge --without-sge --without-loadleveler \
                      --without-tm --without-lsf --without-slurm --without-pvfs2 --without-plfs \
                      --without-cuda --disable-oshmem --disable-mpi-io  --disable-io-romio
                    make -sj4
                    make install
                else
                    echo "OpenMPI installed..."
                    find $TRAVIS_ROOT -name mpiexec
                    find $TRAVIS_ROOT -name mpif90
                    find $TRAVIS_ROOT -name mpicc
                fi


                ;;
            *)
                echo "Unknown MPI implementation: $MPI_IMPL"
                exit 20
                ;;
        esac
        ;;
esac

