!==================
  module data_pdb
!==================

  implicit none
  public

  integer                                   :: nproc
  integer, dimension(:), allocatable        :: nel
  integer, dimension(:), allocatable        :: nel_fluid
  integer, dimension(:), allocatable        :: nel_solid
  integer, dimension(:,:), allocatable      :: procel, inv_procel
  integer, dimension(:,:), allocatable      :: procel_fluid, procel_fluidp
  integer, dimension(:,:), allocatable      :: procel_solid, procel_solidp
  integer, dimension(:,:), allocatable      :: inv_procel_solidp
  integer, dimension(:,:), allocatable      :: inv_procel_fluidp
  integer, dimension(:,:,:,:), allocatable  :: global_index
  integer, dimension(:), allocatable        :: nbelong
  integer, dimension(:), allocatable        :: nprocb
  integer, dimension(:,:), allocatable      :: lprocb
  integer, dimension(:), allocatable        :: el2proc 

! Glocal message passing...redundant eventually!
  integer, dimension(:), allocatable        :: sizerecvp, sizesendp
  integer, dimension(:,:), allocatable      :: listrecvp, listsendp
  integer, dimension(:,:), allocatable      :: sizemsgrecvp, sizemsgsendp
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp

! Slocal - Solid message passing, these all go into the database
  integer, dimension(:), allocatable        :: sizerecvp_solid, sizesendp_solid
  integer, dimension(:,:), allocatable      :: listrecvp_solid, listsendp_solid
  integer, dimension(:,:), allocatable      :: sizemsgrecvp_solid, sizemsgsendp_solid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp_solid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp_solid

! Flocal - Fluid message passing, these all go into the database
  integer, dimension(:), allocatable        :: sizerecvp_fluid, sizesendp_fluid
  integer, dimension(:,:), allocatable      :: listrecvp_fluid, listsendp_fluid
  integer, dimension(:,:), allocatable      :: sizemsgrecvp_fluid, sizemsgsendp_fluid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp_fluid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp_fluid

! Solid-fluid
  integer, dimension(:,:,:,:), allocatable  :: slobal_index, flobal_index
  integer, dimension(:), allocatable        :: nbelong_solid, nbelong_fluid
  integer, dimension(:), allocatable        :: nprocb_solid, nprocb_fluid
  integer, dimension(:,:), allocatable      :: lprocb_solid, lprocb_fluid

  integer, dimension(:,:), allocatable      :: igloc_solid, igloc_fluid
  integer, dimension(:), allocatable        :: nglobp_solid, nglobp_fluid

  integer, dimension(:,:), allocatable      :: slob2sloc, flob2floc

  character(len=6), dimension(:,:), allocatable :: eltypep_solid, eltypep_fluid


! Solid-fluid boundary
  integer, allocatable                      :: nbdry_el(:)
  integer, allocatable                      :: belemp(:,:)
  integer, allocatable                      :: bdry_solid_elp(:,:), bdry_fluid_elp(:,:)
  integer, dimension(:,:), allocatable      :: bdry_jpol_solidp, bdry_jpol_fluidp
  logical, allocatable                      :: have_bdry_elemp(:)

! glocal arrays
  integer, dimension(:,:), allocatable      :: igloc
  integer, dimension(:), allocatable        :: nglobp

! Serendipity arrays
  integer, dimension(:), allocatable        :: nglobmeshp
  double precision, dimension(:,:), allocatable :: scpp, zcpp
  integer, dimension(:,:), allocatable      :: iglobcp
  integer, dimension(:,:,:), allocatable    :: lnodescp
  character(len=6), dimension(:,:), allocatable :: eltypep
  logical, dimension(:,:), allocatable      :: coarsingp

  double precision, dimension(:), allocatable :: theta_min_proc, theta_max_proc

! global to glocal mapping
  integer, dimension(:,:), allocatable      :: glob2gloc

!======================
  end module data_pdb
!======================
