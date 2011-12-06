#!/bin/csh -f

# check for netcdf library
if ( "$1" == "-netcdf" || "$2" == "-netcdf" || $3 == "-netcdf" ) then
set netcdf_exists = `which ncdump |wc -l `;
  if ( $netcdf_exists>0 ) then
  if ( -f nc_routines_true.f90 ) then
    mv -f nc_routines_true.f90 nc_routines.f90
  endif
  if ( -f nc_routines_ghost.f90 ) then
    mv -f nc_routines_ghost.f90 nc_routines.f90.GHOST
  endif
  endif
endif

# check for mpif.h
set mpi_exist = `echo $MPIHOME | wc -l`

if ( $?MPIHOME ) then
  /bin/cp -f $MPIHOME/include/mpif.h .
else
  echo "Please locate your mpi directory and copy the mpif.h to this solver directory."
endif

./perlmakemake.pl $1 $2 $3 
