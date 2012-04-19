!===================
module data_dipole
!===================

use global_parameters
implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Precomputed stiffness terms solid - dipole
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M13s,M33s,M43s

! Matrices for 2nd,3rd-order Terms
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_2, M_6!Mz_xi_p,Mz_xi_m 
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_1, M_5!Mz_eta_p,Mz_eta_m
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_4!Ms_xi_mu
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_8!Ms_xi_2lam
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_3!Ms_eta_mu
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_7!Ms_eta_2lam
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w1!M_w_4_lam_mu
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w3!M_w_mu
  ! -> moved to data_matr.f90

! Matrices for axial elements (defined globally - only npol*nel sized...)
!  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_w_mu
!  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_s_xi_mu
!  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_z_eta_2_lam_mu
!  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_z_eta_4lm_w_l3m

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_dipole
!=======================
