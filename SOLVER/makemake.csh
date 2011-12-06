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
/bin/cp -f $MPIHOME/include/mpif.h .

./makemake.pl $1 $2 $3 $netcdf_exists
