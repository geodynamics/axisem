 !========================
 module unrolled_loops
 !========================

 ! This module has been created by the mesher
 ! for polynomial order  4

   use global_parameters
   use data_proc, only: mynum

   implicit none

   public :: mxm, vxm
   public :: npol_unrolled_loops
   private

   integer, parameter :: npol_unrolled_loops =            4

   contains

 !============================================================

 subroutine mxm(a,b,c)

   include "mesh_params.h" 


   real(kind=realkind), intent(in)  :: a(0: npol,0: npol),b(0: npol,0: npol)
   real(kind=realkind), intent(out) :: c(0: npol,0: npol)
   integer i,j

   do j = 0, npol
     do i = 0, npol
        c(i,j) = sum(a(i,:) * b(:,j))
   !    c(i,j) = & 
   !       + a(i, 0) * b( 0,j) &
   !       + a(i, 1) * b( 1,j) &
   !       + a(i, 2) * b( 2,j) &
   !       + a(i, 3) * b( 3,j) &
   !       + a(i, 4) * b( 4,j)
     end do
   end do 
   return

 end subroutine mxm
 !-------------------------------------------------------------


 !-------------------------------------------------------------
 subroutine vxm(a,b,c)

   include "mesh_params.h" 

   real(kind=realkind), intent(in)  :: a(0: npol),b(0: npol,0: npol)
   real(kind=realkind), intent(out) :: c(0: npol)
   integer j

   do j = 0, npol
       c(j) = sum(a * b(:,j))
   !    c(j) = & 
   !       + a( 0) * b( 0,j) &
   !       + a( 1) * b( 1,j) &
   !       + a( 2) * b( 2,j) &
   !       + a( 3) * b( 3,j) &
   !       + a( 4) * b( 4,j)
   end do 
   !return

 end subroutine vxm
 !-------------------------------------------------------------

 !========================
 end module unrolled_loops
 !========================
