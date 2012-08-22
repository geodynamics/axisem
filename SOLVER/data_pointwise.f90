!===================
 module data_pointwise
!===================
!
! This module is only known during the time loop if the strain tensor 
! is computed on-the-fly. The fluid section is additionally known if global 
! snapshots are dumped (to compute the displacement in the fluid).

use global_parameters

implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!+++++++++++++++++++++++++++++++++++++++++++++++++++
!  Precomputed matrices for pointwise derivatives
!+++++++++++++++++++++++++++++++++++++++++++++++++++

! Inverse density in fluid - only needed when computing pointwise displacement
  real(kind=realkind), allocatable :: inv_rho_fluid(:,:,:)

! (s rho)^-1 in fluid - only for phi comp. of fluid displacement
  real(kind=realkind), allocatable :: inv_s_rho_fluid(:,:,:)

! (s,z) components of fluid potential
  real(kind=realkind), allocatable :: usz_fluid(:,:,:,:)

  real(kind=realkind), allocatable :: deviator(:,:,:,:)


  real(kind=realkind), allocatable :: inv_s_fluid(:,:,:)

! (s)^-1 in solid - needed for the strain tensor, if computed on-the-fly
  real(kind=realkind), allocatable :: inv_s_solid(:,:,:)

! Note that these matrices may include minus signs where necessary.

! Solid region only
  real(kind=realkind), allocatable :: DsDeta_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DzDeta_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DsDxi_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DzDxi_over_J_sol(:,:,:)

! Fluid region only
  real(kind=realkind), allocatable :: DsDeta_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DzDeta_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DsDxi_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DzDxi_over_J_flu(:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
 end module data_pointwise
!=======================
