!===================
module data_heterogeneous
!===================

use global_parameters
implicit none
include "mesh_params.h"
public 

! Heterogeneous region
   integer  :: num_het
   character(len=10), allocatable :: het_format(:),het_funct_type(:)
   character(len=200), allocatable :: het_file_discr(:)
   logical, allocatable :: rdep(:),grad(:)
   logical :: add_up, ani_hetero
   double precision, allocatable :: gradrdep1(:),gradrdep2(:)
   double precision, allocatable :: p_inv_dist(:), R_inv_dist(:)
   double precision, allocatable :: r_het1(:),r_het2(:),th_het1(:),th_het2(:)
   double precision, allocatable :: delta_rho(:), delta_vp(:), delta_vs(:)
   double precision, allocatable :: delta_rho2(:), delta_vp2(:), delta_vs2(:)
   double precision, allocatable :: a_ica(:), b_ica(:), c_ica(:), fa_theta_ica(:), &
                                    fa_phi_ica(:), theta_slices(:)
   integer, allocatable :: inverseshape(:)
   integer :: num_slices
   double precision :: rhetmin, rhetmax, thhetmin, thhetmax



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_heterogeneous
!=======================
