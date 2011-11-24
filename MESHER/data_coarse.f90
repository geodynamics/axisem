module data_coarse
  implicit none
  public
! integer, parameter         :: nc=2
! integer, dimension(0:nc+1) :: iclev
  integer  nc
  integer, dimension(:),allocatable :: iclev
  integer :: ns_ib
! central region :
  integer, dimension(:), allocatable :: iclevc
end module data_coarse
