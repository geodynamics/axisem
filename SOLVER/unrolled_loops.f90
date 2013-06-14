!========================
module unrolled_loops
!========================

  implicit none
  public :: mxm, vxm

  contains

!============================================================

subroutine mxm(a,b,c)

  use global_parameters, only: realkind
  include "mesh_params.h" 
  
  real(kind=realkind), intent(in)  :: a(0: npol,0: npol),b(0: npol,0: npol)
  real(kind=realkind), intent(out) :: c(0: npol,0: npol)
  integer i,j

  do j = 0, npol
     do i = 0, npol
        c(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end subroutine mxm
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine vxm(a,b,c)
  
  use global_parameters, only: realkind
  include "mesh_params.h" 

  real(kind=realkind), intent(in)  :: a(0: npol),b(0: npol,0: npol)
  real(kind=realkind), intent(out) :: c(0: npol)
  integer j

  do j = 0, npol
     c(j) = sum(a * b(:,j))
  end do 

end subroutine vxm
!-------------------------------------------------------------

!========================
end module unrolled_loops
!========================
