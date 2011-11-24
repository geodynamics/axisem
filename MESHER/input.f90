!===================
module input
!===================

  use data_grid
  use data_diag
  use data_coarse
  use data_bkgrdmodel

  implicit none

  public :: read_params
  private

  contains

!----------------------------------------------------------------------------
subroutine read_params

use global_parameters
use data_mesh
use data_spec

character(len=100) :: junk

open(unit=2,file='inparam_mesh')
read(2,*)junk
read(2,*)bkgrdmodel
!af
lfbkgrdmodel=index(bkgrdmodel,' ')-1 
read(2,*)period
read(2,*)nproc_target
read(2,*)junk
read(2,*)resolve_inner_shear
read(2,*)npol
read(2,*)pts_wavelngth
read(2,*)courant
read(2,*) router
read(2,*) dump_mesh_info_files
read(2,*) dump_mesh_info_screen
read(2,*) diagpath
lfdiag=index(diagpath,' ')-1
read(2,*)nc_init
close(2)


write(6,*) ''
write(6,*)'PREDEFINED MODEL/SIMULATION PARAMETERS'
write(6,*)'Background model                 :',bkgrdmodel(1:lfbkgrdmodel)
if (bkgrdmodel=='prem' .or. bkgrdmodel=='prem_light') &
write(6,*)'Resolve inner core shear wave    :',resolve_inner_shear
write(6,*)'Dominant period [s]              :',period
write(6,*)'Elements per dominant wavelength :',pts_wavelngth
write(6,*)'Courant number                   :',courant
write(6,*)'# coarsening levels              :',nc_init
write(6,*)'# processors used in solver      :',nproc_target
write(6,*)' outer radius [m]:',router
write(6,*)' save mesh info files?',dump_mesh_info_files
write(6,*)' print mesh info to screen?',dump_mesh_info_screen
write(6,*)' path to dump output files:',trim(diagpath)
write(6,*)
call flush(6)

if (realkind==4) then 
!  smallval = smallval_sngl
  smallval = smallval_dble
elseif (realkind==8) then 
  smallval = smallval_dble
endif
smallval=1.E-10
end subroutine read_params
!----------------------------------------------------------------------------

!===================
end module input
!===================
