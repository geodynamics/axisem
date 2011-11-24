!===================
module data_matr
!===================
!
! Global arrays (i.e. defined on each GLL point) that are
! needed for the mass, stiffness and boundary terms of the 
! temporal ODE. 
! Note that these are only the arrays common to all source types,
! additional arrays specific to monopole, dipole, and quadrupole 
! are defined in data_<mono,di,quadru>pole.f90

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
! Static solid matrices for all source types
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M11s,M21s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M41s,M12s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M22s,M42s
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M11z,M21z,M41z
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M32s

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
