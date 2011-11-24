!=====================
  module apply_masks 
!=====================

  use global_parameters

  implicit none

  public :: apply_axis_mask_onecomp,apply_axis_mask_twocomp
  public :: apply_axis_mask_threecomp,apply_axis_mask_scal
  private

  contains

!//////////////////////////////////////////////////////


!dk apply_axis_mask_scal--------------------------------------
subroutine apply_axis_mask_scal(u,nel,ax_array,nax_array)
!
!       This routine applies a so-called velocity mask
! by retaining in array v those components  of array u which
! do not belong to the axis of rotation. It sets to zero
! the axial components of the array, for the non-axisymmetric
! components of the variables have to vanish on the axis
! of rotation
! TNM
! TARJE: essentially the same, but FOR ONE COMPONENT 
! and avoiding many loops by using 
! predefined array containing the indices of axial elements
! Also input/output the same array.

  include "mesh_params.h"
  integer, intent(in)               :: nel,nax_array
  real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
    u(0,0:npol,ax_array(ielem))=zero
  end do

end subroutine apply_axis_mask_scal
!---------------------------------------------------------

!dk apply_axis_mask_onecomp--------------------------------------
subroutine apply_axis_mask_onecomp(u,nel,ax_array,nax_array)
!
!       This routine applies a so-called velocity mask
! by retaining in array v those components  of array u which
! do not belong to the axis of rotation. It sets to zero
! the axial components of the array, for the non-axisymmetric
! components of the variables have to vanish on the axis
! of rotation
! TNM
! TARJE: essentially the same, but FOR ONE COMPONENT 
! and avoiding many loops by using 
! predefined array containing the indices of axial elements
! Also input/output the same array.
  
  include "mesh_params.h"
  integer, intent(in)               :: nel,nax_array
  real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel,3)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,0:npol,ax_array(ielem),1)=zero
  end do

end subroutine apply_axis_mask_onecomp
!---------------------------------------------------------

!dk apply_axis_mask_twocomp--------------------------------------
subroutine apply_axis_mask_twocomp(u,nel,ax_array,nax_array)
!
!       This routine applies a so-called velocity mask
! by retaining in array v those components  of array u which
! do not belong to the axis of rotation. It sets to zero
! the axial components of the array, for the non-axisymmetric
! components of the variables have to vanish on the axis
! of rotation
! TNM
! TARJE: essentially the same, but FOR ONE COMPONENT 
! and avoiding many loops by using 
! predefined array containing the indices of axial elements
! Also input/output the same array.
  
  include "mesh_params.h"
  integer, intent(in)               :: nel,nax_array
  real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel,3)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,0:npol,ax_array(ielem),2:3)=zero
  end do

end subroutine apply_axis_mask_twocomp
!---------------------------------------------------------

!dk apply_axis_mask_threecomp--------------------------------------
subroutine apply_axis_mask_threecomp(u,nel,ax_array,nax_array)

!       This routine applies a so-called velocity mask
! by retaining in array v those components  of array u which
! do not belong to the axis of rotation. It sets to zero
! the axial components of the array, for the non-axisymmetric
! components of the variables have to vanish on the axis
! of rotation
! TNM
! TARJE: essentially the same, but FOR ONE COMPONENT 
! and avoiding many loops by using 
! predefined array containing the indices of axial elements
! Also input/output the same array.

  include "mesh_params.h"
  integer, intent(in)               :: nel,nax_array
  real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel,3)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,0:npol,ax_array(ielem),1:3)=zero
  end do

end subroutine apply_axis_mask_threecomp
!---------------------------------------------------------

!//////////////////////////////////////////////////////
!
!========================
  end module apply_masks 
!========================
