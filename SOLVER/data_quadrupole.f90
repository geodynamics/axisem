!===================
 module data_quadrupole
!===================

use global_parameters
implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Precomputed stiffness terms solid - quadrupole
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M1phi,M2phi,M4phi
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Ms_eta_lam
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Ms_xi_lam
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_eta_lam
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_xi_lam
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_eta_2mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_xi_2mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Ms_eta_2mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Ms_xi_2mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w_lam_6mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w_min_2lam_6mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w_4lam_9mu
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w_4mu

! z-component: (4mu)M_zeta
  real(kind=realkind), dimension(0:npol,nel_solid) :: M0_w_4mu

! s-component (times us): (2lambda)M_zeta+(lambda+6mu)M_w = 3(lambda+2mu)M_w
  real(kind=realkind), dimension(0:npol,nel_solid) :: M0_z_eta_2l_w_l_6m

! phi-component (times uphi): (-2mu)M_zeta+(4lambda+9mu)M_w = (4lambda+7mu)M_w
  real(kind=realkind), dimension(0:npol,nel_solid) :: M0_z_eta_3m_w_4l_9m 

! s- and phi components: 2(mu-lambda)M_zeta-(2lambda+6mu)M_w = -4(lambda+mu)M_w
  real(kind=realkind), dimension(0:npol,nel_solid) :: M0_z_eta_2lm_min_w_2l_6m

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_quadrupole
!=======================
