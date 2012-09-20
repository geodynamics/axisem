#!/bin/csh -f

# check for netcdf library
if ( "$1" == "-netcdf" || "$2" == "-netcdf" || $3 == "-netcdf" ) then
  set netcdf_exists = `which ncdump |wc -l `;
  if ( $netcdf_exists>0 ) then
    echo "netcdf exists, moving routines"
    /bin/cp -f nc_routines.f90.TRUE nc_routines.f90
  else 
    set netcdf_exists = 0
    /bin/cp -f nc_routines.f90.GHOST nc_routines.f90
  endif
else 
  set netcdf_exists = 0
  /bin/cp -f nc_routines.f90.GHOST nc_routines.f90
endif

echo $1 $2 $3
if (${#argv} < 1) then 
./perlmakemake.pl
else if (${#argv} < 2) then 
./perlmakemake.pl $1 
else if (${#argv} < 3) then 
./perlmakemake.pl $1 $2 
else if (${#argv} < 4) then 
./perlmakemake.pl $1 $2 $3
endif

if ( ! -d MESHES ) then
  mkdir MESHES
endif
