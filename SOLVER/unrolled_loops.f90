 !========================
 module unrolled_loops
 !========================

 ! This module has been created by the mesher
 ! for polynomial order  4
 ! on 10/13/2012, at 21h 16min

   use global_parameters

   implicit none

   public :: mxm,vxm
   private

   contains

 !============================================================

 subroutine mxm(a,b,c)

   include "mesh_params.h" 

   real(kind=realkind), intent(in)  :: a(0: 4,0: 4),b(0: 4,0: 4)
   real(kind=realkind), intent(out) :: c(0: 4,0: 4)
   integer i,j

   if ( npol /=  4 ) then
      write(6,*)"Problem: unrolled_loops.f90 has different" 
      write(6,*)"         polynomial order than mesh_params.h:" 
      write(6,*)"         4",npol
      stop
   endif

   do j = 0, 4
     do i = 0, 4
       c(i,j) = & 
          + a(i, 0) * b( 0,j) &
          + a(i, 1) * b( 1,j) &
          + a(i, 2) * b( 2,j) &
          + a(i, 3) * b( 3,j) &
          + a(i, 4) * b( 4,j)
     end do
   end do 
   return

 end subroutine mxm
 !-------------------------------------------------------------


 !-------------------------------------------------------------
 subroutine vxm(a,b,c)

   include "mesh_params.h" 

   real(kind=realkind), intent(in)  :: a(0: 4),b(0: 4,0: 4)
   real(kind=realkind), intent(out) :: c(0: 4)
   integer j

   if ( npol /=  4 ) then
      write(6,*)"Problem: unrolled_loops.f90 has different" 
      write(6,*)"         polynomial order than mesh_params.h:" 
      write(6,*)"         4",npol
      stop
   endif

   do j = 0, 4
       c(j) = & 
          + a( 0) * b( 0,j) &
          + a( 1) * b( 1,j) &
          + a( 2) * b( 2,j) &
          + a( 3) * b( 3,j) &
          + a( 4) * b( 4,j)
   end do 
   return

 end subroutine vxm
 !-------------------------------------------------------------

 !========================
 end module unrolled_loops
 !========================
