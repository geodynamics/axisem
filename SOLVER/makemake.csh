#!/bin/csh -f

# check for netcdf library
set netcdf_exists = `which ncdump |wc -l `;
if ( $netcdf_exists>0 ) then
  if ( -f nc_routines_true.f90 ) then
    mv -f nc_routines_true.f90 nc_routines.f90
  endif
  if ( -f nc_routines_ghost.f90 ) then
    mv -f nc_routines_ghost.f90 nc_routines.f90.GHOST
  endif
endif

# check for mpif.h
set mpi_exist = `echo $MPIHOME | wc -l`
if ( $mpi_exist > 0 ) then
/bin/cp -f $MPIHOME/include/mpif.h .
else
echo "Please locate your mpi directory and copy the mpif.h to this solver directory."
endif

./perlmakemake.pl $1 $2 $3 $netcdf_exists
