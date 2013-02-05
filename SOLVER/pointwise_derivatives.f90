!========================
MODULE pointwise_derivatives
!========================
!
! Various forms of the two basic spatial derivatives d/ds and d/dz.
! Pointwise refers to the notion that these derivatives are not embedded 
! into any integral, but merely the spectral-element based derivative. 
! These are needed to compute the source term, the displacement in the fluid,
! the strain tensor, and the axial expression f/s=df/ds (L'Hospital's rule).

USE global_parameters
USE data_mesh
USE data_spec

IMPLICIT NONE

PUBLIC :: axisym_laplacian_solid,axisym_laplacian_solid_add
PUBLIC :: axisym_laplacian_fluid,axisym_laplacian_fluid_add
PUBLIC :: dsdf_elem_solid,dzdf_elem_solid
PUBLIC :: dsdf_fluid_axis,dsdf_fluid_allaxis,dsdf_solid_allaxis

PRIVATE
CONTAINS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------------
subroutine axisym_laplacian_solid(f,lap)
!
! Computes the axisymmetric laplacian of scalar field f in the solid region:
! lap = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
!
! MvD: I do not understand the nomenclature here, as IMHO this is NOT a
!      laplacian (= divergence of the gradient), but rather the gradient itself
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
use data_pointwise, ONLY: DsDeta_over_J_sol,DsDxi_over_J_sol
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_solid)
REAL(kind=realkind),intent(out)              :: lap(0:npol,0:npol,nel_solid,2)
INTEGER                                      :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta

  do iel = 1, nel_solid

    dzdeta = DzDeta_over_J_sol(:,:,iel)
    dzdxi  = DzDxi_over_J_sol(:,:,iel)
    dsdeta = DsDeta_over_J_sol(:,:,iel)
    dsdxi  = DsDxi_over_J_sol(:,:,iel)

    if (axis_solid(iel)) then 
       call mxm(G1T,f(:,:,iel),mxm1) ! axial elements
    else 
       call mxm(G2T,f(:,:,iel),mxm1) ! non-axial elements
    endif 
    call mxm(f(:,:,iel),G2,mxm2)
    call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf,nsize)
    call collocate2_sum_1d(dsdeta,mxm1,dsdxi,mxm2,dzdf,nsize)
    lap(:,:,iel,1)=dsdf
    lap(:,:,iel,2)=dzdf
 enddo

end subroutine axisym_laplacian_solid
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_laplacian_solid_add(f,lap)
!
! Computes the axisymmetric laplacian of scalar field f in the solid region:
! lap = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
! This routine takes a previously calculated derivative and adds it
! to the result computed here in a permuted fashion.
! This saves the strain dump output two global fields, as the strain 
! trace will hereby be dumped as well as the entire E_31 term instead 
! of its two cross-derivative contributions.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
use data_pointwise, ONLY: DsDeta_over_J_sol,DsDxi_over_J_sol
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_solid)
REAL(kind=realkind),intent(inout)            :: lap(0:npol,0:npol,nel_solid,2)
INTEGER                                      :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta
REAL(kind=realkind),DIMENSION(0:npol,0:npol,2) :: lap_old

  do iel = 1, nel_solid

    dzdeta = DzDeta_over_J_sol(:,:,iel)
    dzdxi  = DzDxi_over_J_sol(:,:,iel)
    dsdeta = DsDeta_over_J_sol(:,:,iel)
    dsdxi  = DsDxi_over_J_sol(:,:,iel)

    if (axis_solid(iel)) then 
       call mxm(G1T,f(:,:,iel),mxm1) ! axial elements
    else 
       call mxm(G2T,f(:,:,iel),mxm1) ! non-axial elements
    endif 
    call mxm(f(:,:,iel),G2,mxm2)
    call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf,nsize)
    call collocate2_sum_1d(dsdeta,mxm1,dsdxi,mxm2,dzdf,nsize)

    lap_old(0:npol,0:npol,1) = lap(0:npol,0:npol,iel,2)
    lap_old(0:npol,0:npol,2) = lap(0:npol,0:npol,iel,1)

    lap(0:npol,0:npol,iel,1) = lap_old(0:npol,0:npol,1) + dsdf(0:npol,0:npol)
    lap(0:npol,0:npol,iel,2) = lap_old(0:npol,0:npol,2) + dzdf(0:npol,0:npol)

 enddo

end subroutine axisym_laplacian_solid_add
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_elem_solid(dsdf,f,iel)
!
! Computes the elemental s-derivative of scalar field f in the solid region.
! This is used to compute the source term within the source element only.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol)
REAL(kind=realkind),intent(out)              :: dsdf(0:npol,0:npol)
INTEGER,intent(in)                           :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dzdxi,dzdeta
  
  dzdeta = DzDeta_over_J_sol(:,:,iel)
  dzdxi  = DzDxi_over_J_sol(:,:,iel)

  if (axis_solid(iel)) then 
     call mxm(G1T,f,mxm1) ! axial elements
  else 
     call mxm(G2T,f,mxm1) ! non-axial elements
  endif
  call mxm(f,G2,mxm2)

  call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf,nsize)

end subroutine dsdf_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dzdf_elem_solid(dzdf,f,iel)
!
! Computes the elemental z-derivative of scalar field f in the solid region.
! This is used to compute the source term within the source element only.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DsDeta_over_J_sol,DsDxi_over_J_sol
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol)
REAL(kind=realkind),intent(out)              :: dzdf(0:npol,0:npol)
INTEGER,intent(in)                           :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dsdxi,dsdeta
  
  dsdeta = DsDeta_over_J_sol(:,:,iel)
  dsdxi  = DsDxi_over_J_sol(:,:,iel)

  if (axis_solid(iel)) then 
     call mxm(G1T,f,mxm1) ! axial elements
  else 
     call mxm(G2T,f,mxm1) ! non-axial elements
  endif
  call mxm(f,G2,mxm2)
  call collocate2_sum_1d(dsdeta,mxm1,dsdxi,mxm2,dzdf,nsize)

end subroutine dzdf_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_solid_allaxis(f,dsdf)
!
! Computes the pointwise derivative of scalar f in the s-direction 
! within the solid region, ONLY AT THE AXIS (needed for solid displacement)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_solid)
REAL(kind=realkind),intent(out)              :: dsdf(0:npol,naxel_solid)
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dzdxi,dzdeta,dsdf_el
integer                                      :: ielem,iel

  do ielem=1, naxel_solid
    iel = ax_el_solid(ielem) 
    dzdeta = DzDeta_over_J_sol(:,:,iel)
    dzdxi  = DzDxi_over_J_sol(:,:,iel)
    call mxm(G1T, f(:,:,iel), mxm1) 
    call mxm(f(:,:,iel), G2, mxm2)
    call collocate2_sum_1d(dzdeta, mxm1, dzdxi, mxm2, dsdf_el, nsize)
    dsdf(:,ielem) = dsdf_el(0,:)
  enddo

end subroutine dsdf_solid_allaxis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine axisym_laplacian_fluid(f,lap)
!
! Computes the axisymmetric laplacian of scalar field f in the fluid region:
! lap = \nabla {f}  = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
!
! MvD: I do not understand the nomenclature here, as IMHO this is NOT a
!      laplacian (= divergence of the gradient), but rather the gradient itself
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
use data_pointwise, ONLY: DsDeta_over_J_flu,DsDxi_over_J_flu
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_fluid)
REAL(kind=realkind),intent(out)              :: lap(0:npol,0:npol,nel_fluid,2)
INTEGER                                      :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta

  do iel = 1, nel_fluid

    dzdeta = DzDeta_over_J_flu(:,:,iel)
    dzdxi  = DzDxi_over_J_flu(:,:,iel)
    dsdeta = DsDeta_over_J_flu(:,:,iel)
    dsdxi  = DsDxi_over_J_flu(:,:,iel)

    if (axis_fluid(iel)) then 
       call mxm(G1T,f(:,:,iel),mxm1) ! axial elements
    else 
       call mxm(G2T,f(:,:,iel),mxm1) ! non-axial elements
    endif 
    call mxm(f(:,:,iel),G2,mxm2)
    call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf,nsize)
    call collocate2_sum_1d(dsdeta,mxm1,dsdxi,mxm2,dzdf,nsize)
    lap(:,:,iel,1)=dsdf
    lap(:,:,iel,2)=dzdf
 enddo

end subroutine axisym_laplacian_fluid
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_laplacian_fluid_add(f,lap)
!
! Computes the axisymmetric laplacian of scalar field f in the fluid region:
! lap = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
! This routine takes a previously calculated derivative and adds it
! to the result computed here in a permuted fashion.
! This saves the strain dump output two global fields, as the strain 
! trace will hereby be dumped as well as the entire E_31 term instead 
! of its two cross-derivative contributions.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
use data_pointwise, ONLY: DsDeta_over_J_flu,DsDxi_over_J_flu
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_fluid)
REAL(kind=realkind),intent(inout)            :: lap(0:npol,0:npol,nel_fluid,2)
INTEGER                                      :: iel
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta
REAL(kind=realkind),DIMENSION(0:npol,0:npol,2) :: lap_old

  do iel = 1, nel_fluid

    dzdeta = DzDeta_over_J_flu(:,:,iel)
    dzdxi  = DzDxi_over_J_flu(:,:,iel)
    dsdeta = DsDeta_over_J_flu(:,:,iel)
    dsdxi  = DsDxi_over_J_flu(:,:,iel)

    if (axis_fluid(iel)) then 
       call mxm(G1T,f(:,:,iel),mxm1) ! axial elements
    else 
       call mxm(G2T,f(:,:,iel),mxm1) ! non-axial elements
    endif 
    call mxm(f(:,:,iel),G2,mxm2)
    call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf,nsize)
    call collocate2_sum_1d(dsdeta,mxm1,dsdxi,mxm2,dzdf,nsize)

    lap_old(0:npol,0:npol,1) = lap(0:npol,0:npol,iel,2)
    lap_old(0:npol,0:npol,2) = lap(0:npol,0:npol,iel,1)

    lap(0:npol,0:npol,iel,1) = lap_old(0:npol,0:npol,1) + dsdf(0:npol,0:npol)
    lap(0:npol,0:npol,iel,2) = lap_old(0:npol,0:npol,2) + dzdf(0:npol,0:npol)

 enddo

end subroutine axisym_laplacian_fluid_add
!=============================================================================



!-----------------------------------------------------------------------------
subroutine dsdf_fluid_axis(f,iel,jpol,dsdf)
!
! Computes the pointwise derivative of scalar f in the s-direction 
! within the fluid region, ONLY AT THE AXIS (needed for fluid displacement)
! and for a specific element iel and etsa coordinate index jpol.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

INTEGER,INTENT(in)                           :: iel,jpol
REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol)
REAL(kind=realkind),intent(out)              :: dsdf
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dzdxi,dzdeta,dsdf_el

  dzdeta = DzDeta_over_J_flu(:,:,iel)
  dzdxi  = DzDxi_over_J_flu(:,:,iel)
  call mxm(G1T,f,mxm1)
  call mxm(f,G2,mxm2)
  call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf_el,nsize)
  dsdf=dsdf_el(0,jpol)

end subroutine dsdf_fluid_axis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_fluid_allaxis(f,dsdf)
!
! Computes the pointwise derivative of scalar f in the s-direction 
! within the fluid region, ONLY AT THE AXIS (needed for fluid displacement)
! for all axial elements.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
use unrolled_loops
use unit_stride_colloc

include 'mesh_params.h'

REAL(kind=realkind),INTENT(in)               :: f(0:npol,0:npol,nel_fluid)
REAL(kind=realkind),intent(out)              :: dsdf(0:npol,naxel_fluid)
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: mxm1,mxm2
REAL(kind=realkind),DIMENSION(0:npol,0:npol) :: dzdxi,dzdeta,dsdf_el
integer                                      :: ielem,iel

  do ielem=1,naxel_fluid
    iel=ax_el_fluid(ielem) 
    dzdeta = DzDeta_over_J_flu(:,:,iel)
    dzdxi  = DzDxi_over_J_flu(:,:,iel)
    call mxm(G1T,f(:,:,iel),mxm1) 
    call mxm(f(:,:,iel),G2,mxm2)
    call collocate2_sum_1d(dzdeta,mxm1,dzdxi,mxm2,dsdf_el,nsize)
    dsdf(:,ielem)=dsdf_el(0,:)
  enddo

end subroutine dsdf_fluid_allaxis
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!========================
end MODULE pointwise_derivatives
!========================
