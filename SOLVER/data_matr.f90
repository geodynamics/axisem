!===================
module data_matr
!===================
!
! Global arrays (i.e. defined on each GLL point) that are
! needed for the mass, stiffness and boundary terms of the 
! temporal ODE. 

use global_parameters

implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  !	Mass matrix arrays
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: inv_mass_rho
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid) :: inv_mass_fluid
  real(kind=realkind), dimension(:,:,:), allocatable :: unassem_mass_rho_solid
  real(kind=realkind), dimension(:,:,:), allocatable :: unassem_mass_lam_fluid 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  !       Precomputed stiffness matrices
  !++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Static solid matrices for all source types:
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M11s,M21s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M41s,M12s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M22s,M42s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M11z,M21z,M41z
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M32s
  ! for dipole: 
  real(kind=realkind), allocatable :: M13s(:,:,:), M33s(:,:,:), M43s(:,:,:)
  ! for quadrupole: 
  real(kind=realkind), allocatable :: M1phi(:,:,:), M2phi(:,:,:), M4phi(:,:,:)
  
  
  ! for all source types:
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_1, M_2, M_3, M_4
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M_5(:,:,:), M_6(:,:,:), M_7(:,:,:), M_8(:,:,:) 
  
  
  ! for all source types:
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w1
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M_w2(:,:,:), M_w3(:,:,:)
  ! for quadrupole: 
  real(kind=realkind), allocatable :: M_w4(:,:,:), M_w5(:,:,:)
  

  ! for all source types:
  real(kind=realkind), dimension(0:npol,nel_solid) :: M0_w1, M0_w2, M0_w3
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M0_w4(:,:), M0_w5(:,:), M0_w6(:,:)
  ! for dipole
  real(kind=realkind), allocatable :: M0_w7(:,:), M0_w8(:,:), M0_w9(:,:), M0_w10(:,:)

  ! Fluid matrices
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid) :: M1chi_fl,M2chi_fl
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid) :: M4chi_fl
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid) :: M_w_fl
  real(kind=realkind), dimension(0:npol,nel_fluid)        :: M0_w_fl

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  !       Precomputed solid-fluid boundary matrices
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Solid-fluid boundary matrix
  real(kind=realkind), dimension(0:npol,nel_bdry,2)    :: bdry_matr
  real(kind=realkind), dimension(0:npol,nel_bdry,2)    :: bdry_matr_fluid
  real(kind=realkind), dimension(0:npol,nel_bdry,2)    :: bdry_matr_solid
  double precision, dimension(nel_bdry)                :: solflubdry_radius

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=====================
end module data_matr
!=====================
