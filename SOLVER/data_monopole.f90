!===================
module data_monopole
!===================

use global_parameters

implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Precomputed stiffness terms solid - monopole
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_1, M_2, M_3, M_4
  !real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w1
  !real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_w1, M0_w2
  ! -> moved to data_matr (as it is also needed for quadrupole now)
  !real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_4 ! is now M0_w3

! In the isotropic case (non-axial):
! M_1 =  lambda_ij * sigma_i * sigma_j * z_eta_ij
! M_2 = -lambda_ij * sigma_i * sigma_j * z_xi_ij
! M_3 = -lambda_ij * sigma_i * sigma_j * s_eta_ij
! M_4 =  lambda_ij * sigma_i * sigma_j * s_xi_ij
!
! In the isotropic case (axial, i > 0):
! M_1 =  lambda_ij * bar{sigma}_i / (1 + xi_i) * sigma_j * z_eta_ij
! M_2 = -lambda_ij * bar{sigma}_i / (1 + xi_i) * sigma_j * z_xi _ij
! M_3 = -lambda_ij * bar{sigma}_i / (1 + xi_i) * sigma_j * s_eta_ij
! M_4 =  lambda_ij * bar{sigma}_i / (1 + xi_i) * sigma_j * s_xi _ij
!
! In the isotropic case (axial, i = 0):
! M_1 = 0
! M_2 = 0
! M_3 = 0
! M_4 = 0
!
! M0_w1 = (3 * lambda_ij + 2 * mu) * bar{sigma}_i * sigma_j * J_ij / s_xi_ij
! M0_w2 = 0
! M0_w3 = lambda_ij * bar{sigma}_i * sigma_j * s_xi _ij


! In the anisotropic case (non-axial):
! M_1 = sigma_i * sigma_j * ( C12_ij * z_eta_ij - C25_ij * s_eta_ij)
! M_2 = sigma_i * sigma_j * (-C12_ij * z_xi_ij  + C25_ij * s_xi_ij )
! M_3 = sigma_i * sigma_j * (-C23_ij * s_eta_ij + C25_ij * z_eta_ij)
! M_4 = sigma_i * sigma_j * ( C23_ij * s_xi_ij  - C25_ij * z_xi_ij )
!
! In the anisotropic case (axial, i > 0):
! M_1 = bar{sigma}_i / (1 + xi_i) * sigma_j * ( C12_ij * z_eta_ij - C25_ij * s_eta_ij)
! M_2 = bar{sigma}_i / (1 + xi_i) * sigma_j * (-C12_ij * z_xi_ij  + C25_ij * s_xi_ij )
! M_3 = bar{sigma}_i / (1 + xi_i) * sigma_j * (-C23_ij * s_eta_ij + C25_ij * z_eta_ij)
! M_4 = bar{sigma}_i / (1 + xi_i) * sigma_j * ( C23_ij * s_xi_ij  - C25_ij * z_xi_ij )
!
! In the anisotropic case (axial, i = 0):
! M_1 = 0
! M_2 = 0
! M_3 = 0
! M_4 = 0
!
! M0_w1 = bar{sigma}_i * sigma_j * J_ij / s_xi_ij * (2 * C12 + C22)
! M0_w2 = bar{sigma}_i * sigma_j * J_ij / s_xi_ij * C25
! M0_w3 = bar{sigma}_i * sigma_j * (C23 s_xi _ij - C25 z_xi_ij)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=====================
end module data_monopole
!=======================
