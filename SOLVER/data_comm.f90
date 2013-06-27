!===================
!> Message-passing communication variables
!! Note: For the easy serialization of the code (see commun.f90 and commpi.f90),
!! one still needs this module as these quantities are read in from the mesher.
module data_comm
!===================


  use global_parameters
  
  implicit none
  public 

  integer                              :: comm
  integer                              :: mpi_realkind

  ! Solid domain message passing

  integer                              :: num_send_gll, num_recv_gll
  integer                              :: sizemsgrecvmax_solid
  integer                              :: sizemsgsendmax_solid
  integer                              :: sizemsgmax_solid
  integer, allocatable                 :: glob2el_send(:,:), glob2el_recv(:,:)
  real(kind=realkind), allocatable     :: buffs_solid(:,:), buffr_solid(:,:)
  real(kind=realkind), allocatable     :: buffs_all(:,:,:), buffr_all(:,:,:)

  integer                              :: sizerecv_solid, sizesend_solid
  integer, dimension(:),   allocatable :: listrecv_solid, sizemsgrecv_solid
  integer, dimension(:),   allocatable :: listsend_solid, sizemsgsend_solid
  integer, dimension(:,:), allocatable :: glocal_index_msg_recv_solid
  integer, dimension(:,:), allocatable :: glocal_index_msg_send_solid

  ! Fluid domain message passing
  integer                              :: sizemsgrecvmax_fluid
  integer                              :: sizemsgsendmax_fluid
  integer                              :: sizemsgmax_fluid
  real(kind=realkind), allocatable     :: buffs_fluid(:), buffr_fluid(:)
  integer                              :: sizerecv_fluid, sizesend_fluid
  integer, dimension(:),   allocatable :: listrecv_fluid, sizemsgrecv_fluid
  integer, dimension(:),   allocatable :: listsend_fluid, sizemsgsend_fluid
  integer, dimension(:,:), allocatable :: glocal_index_msg_recv_fluid
  integer, dimension(:,:), allocatable :: glocal_index_msg_send_fluid

!=======================
end module data_comm
!=======================
