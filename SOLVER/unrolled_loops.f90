!========================
!> Routines for general matrix-matrix and matrix-vector multiplication. Called a 
!! bazillion times, presumably fast.
module unrolled_loops
!========================

  implicit none
  public :: mxm, vxm

  contains

!============================================================

!-------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! Size is fixed to npol x npol
subroutine mxm(a,b,c)

  use global_parameters, only: realkind
  include "mesh_params.h" 
  
  real(kind=realkind), intent(in)  :: a(0: npol,0: npol),b(0: npol,0: npol) !< Input matrices
  real(kind=realkind), intent(out) :: c(0: npol,0: npol)                    !< Result
  integer                          :: i, j

  do j = 0, npol
     do i = 0, npol
        c(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end subroutine mxm
!-------------------------------------------------------------

!-------------------------------------------------------------
!> Multiplies vector a leftwise to matrix b to have vector c.
!! Size is fixed to npol x npol
subroutine vxm(a,b,c)
  
  use global_parameters, only: realkind
  include "mesh_params.h" 

  real(kind=realkind), intent(in)  :: a(0: npol)         !< Vector a
  real(kind=realkind), intent(in)  :: b(0: npol,0: npol) !< Matrix b
  real(kind=realkind), intent(out) :: c(0: npol)         !< Resulting vector c
  integer                          :: j

  do j = 0, npol
     c(j) = sum(a * b(:,j))
  end do 

end subroutine vxm
!-------------------------------------------------------------

!========================
end module unrolled_loops
!========================
