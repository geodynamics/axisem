!===================
module data_heterogeneous
!===================

use global_parameters
implicit none
include "mesh_params.h"
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Heterogeneous region
! Edited by E. Vanacore on 28/10/2011 to allow for multiple heterogeneities
! Edited by fanie on 20/03/2012 to allow for multiple heterogeneities of different shapes
   integer  :: num_het,nhet_pts
   character(len=10), allocatable :: het_format(:),het_funct_type(:)
   character(len=200), allocatable :: het_file_discr(:)
   logical, allocatable :: rdep(:),grad(:)
   double precision, allocatable :: gradrdep1(:),gradrdep2(:)
   double precision, allocatable :: r_het1(:),r_het2(:),th_het1(:),th_het2(:)
   double precision, allocatable :: ph_het1(:),ph_het2(:)
   double precision, allocatable :: delta_rho(:),delta_vp(:),delta_vs(:)
   double precision, allocatable :: delta_rho2(:,:),delta_vp2(:,:),delta_vs2(:,:)
   integer, allocatable :: inverseshape(:)
!!fanie: only needed in load_het_funct - clean up?
   double precision, allocatable :: rhet(:),thhet(:),phhet(:) 
!!fanie: only needed in load_het_discr and plot_discrete_input - clean up?
   double precision, allocatable :: rhet2(:,:),thhet2(:,:),phhet2(:,:) 
   double precision :: rhetmin,rhetmax,thhetmin,thhetmax
! End of Edits 28/10/2011


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_heterogeneous
!=======================
