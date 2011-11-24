!===================
module data_monopole
!===================

use global_parameters

implicit none
public 
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Precomputed stiffness terms solid - monopole
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_eta,Ms_eta
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: Mz_xi,Ms_xi
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid) :: M_w
  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_w
  real(kind=realkind), dimension(0:npol,nel_solid)        :: M0_s_xi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=====================
end module data_monopole
!=======================
