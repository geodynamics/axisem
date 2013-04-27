!========================
MODULE pointwise_derivatives
!========================
  !
  ! Various forms of the two basic spatial derivatives d/ds and d/dz.
  ! Pointwise refers to the notion that these derivatives are not embedded 
  ! into any integral, but merely the spectral-element based derivative. 
  ! These are needed to compute the source term, the displacement in the fluid,
  ! the strain tensor, and the axial expression f/s=df/ds (L'Hospital's rule).
  
  use global_parameters
  use data_mesh
  use data_spec
  
  implicit none
  
  public :: axisym_gradient_solid, axisym_gradient_solid_add
  public :: axisym_gradient_solid_el
  public :: axisym_gradient_solid_el_cg4
  public :: axisym_gradient_fluid, axisym_gradient_fluid_add
  public :: dsdf_elem_solid, dzdf_elem_solid
  public :: dsdf_fluid_axis, dsdf_fluid_allaxis, dsdf_solid_allaxis
  public :: axisym_dsdf_solid
  public :: f_over_s_solid
  public :: f_over_s_solid_el
  public :: f_over_s_solid_el_cg4
  public :: f_over_s_fluid
  
  private

contains

!-----------------------------------------------------------------------------------------
function f_over_s_solid_el_cg4(f, iel)
  !
  ! computes f/s using L'Hospital's rule lim f/s = lim df/ds at the axis (s = 0)
  !
  use data_pointwise,           ONLY: inv_s_solid
  use data_mesh,                ONLY: naxel_solid, ax_el_solid

  include 'mesh_params.h'
  
  real(kind=realkind),intent(in) :: f(0:npol,0:npol)
  integer,intent(in)             :: iel
  real(kind=realkind)            :: f_over_s_solid_el_cg4(1:4)
  real(kind=realkind)            :: dsdf(0:npol,0:npol)
  
  ! in the bulk:
  f_over_s_solid_el_cg4(1) = inv_s_solid(1,1,iel) * f(1,1)
  f_over_s_solid_el_cg4(2) = inv_s_solid(1,3,iel) * f(1,3)
  f_over_s_solid_el_cg4(3) = inv_s_solid(3,1,iel) * f(3,1)
  f_over_s_solid_el_cg4(4) = inv_s_solid(3,3,iel) * f(3,3)

  ! there is no axis (s=0)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function f_over_s_solid_el(f, iel)
  !
  ! computes f/s using L'Hospital's rule lim f/s = lim df/ds at the axis (s = 0)
  !
  use data_pointwise,           ONLY: inv_s_solid
  use data_mesh,                ONLY: naxel_solid, ax_el_solid

  include 'mesh_params.h'
  
  real(kind=realkind),intent(in) :: f(0:npol,0:npol)
  integer,intent(in)             :: iel
  real(kind=realkind)            :: f_over_s_solid_el(0:npol,0:npol)
  real(kind=realkind)            :: dsdf(0:npol,0:npol)
  
  ! in the bulk:
  f_over_s_solid_el = inv_s_solid(:,:,iel) * f

  ! at the axis:
  if (axis_solid(iel)) then
     ! TODO: Optimize this computing the derivative only for i=0
     call dsdf_elem_solid(dsdf,f,iel)
     f_over_s_solid_el(0,:) = dsdf(0,:)
  endif

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function f_over_s_solid(f)
  !
  ! computes f/s using L'Hospital's rule lim f/s = lim df/ds at the axis (s = 0)
  !
  use data_pointwise,           ONLY: inv_s_solid
  use data_mesh,                ONLY: naxel_solid, ax_el_solid

  include 'mesh_params.h'
  
  real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind)            :: f_over_s_solid(0:npol,0:npol,nel_solid)
  real(kind=realkind)            :: dsdf(0:npol,naxel_solid)
  integer                        :: iel
  
  ! in the bulk:
  f_over_s_solid = inv_s_solid * f

  ! at the axis:
  call dsdf_solid_allaxis(f, dsdf) ! axial f/s
  do iel=1, naxel_solid
     f_over_s_solid(0,:,ax_el_solid(iel)) = dsdf(:,iel)
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function f_over_s_fluid(f)
  !
  ! computes f/s using L'Hospital's rule lim f/s = lim df/ds at the axis (s = 0)
  !
  use data_pointwise,           ONLY: inv_s_fluid
  use data_mesh,                ONLY: naxel_fluid, ax_el_fluid

  include 'mesh_params.h'
  
  real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_fluid)
  real(kind=realkind)            :: f_over_s_fluid(0:npol,0:npol,nel_fluid)
  real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
  integer                        :: iel
  
  ! in the bulk:
  f_over_s_fluid = inv_s_fluid * f

  ! at the axis:
  call dsdf_fluid_allaxis(f, dsdf) ! axial f/s
  do iel=1, naxel_fluid
     f_over_s_fluid(0,:,ax_el_fluid(iel)) = dsdf(:,iel)
  enddo

end function
!-----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------
subroutine axisym_dsdf_solid(f, dsdf)
  !
  ! Computes the partial derivative
  ! dsdf = \partial_s(f)
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol, DzDxi_over_J_sol
  use unrolled_loops
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind),intent(out)              :: dsdf(0:npol,0:npol,nel_solid)
  integer                                      :: iel
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1, mxm2
  real(kind=realkind),dimension(0:npol,0:npol) :: dzdeta, dzdxi

  do iel = 1, nel_solid

     dzdeta = DzDeta_over_J_sol(:,:,iel)
     dzdxi  = DzDxi_over_J_sol(:,:,iel)

     if (axis_solid(iel)) then 
        call mxm(G1T,f(:,:,iel),mxm1) ! axial elements
     else 
        call mxm(G2T,f(:,:,iel),mxm1) ! non-axial elements
     endif 
     call mxm(f(:,:,iel),G2,mxm2)

     dsdf(:,:,iel) = dzdeta * mxm1 + dzdxi * mxm2
  enddo

end subroutine
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_gradient_solid_el_cg4(f,grad,iel)
  !
  ! Computes the axisymmetric gradient of scalar field f in the solid region:
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  !use data_pointwise, ONLY: DzDeta_over_J_sol, DzDxi_over_J_sol
  !use data_pointwise, ONLY: DsDeta_over_J_sol, DsDxi_over_J_sol
  use data_pointwise, ONLY: DzDeta_over_J_sol_cg4, DzDxi_over_J_sol_cg4
  use data_pointwise, ONLY: DsDeta_over_J_sol_cg4, DsDxi_over_J_sol_cg4
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)        :: f(0:npol,0:npol)
  real(kind=realkind),intent(out)       :: grad(1:4,2)
  integer,intent(in)                    :: iel
  real(kind=realkind),dimension(1:4)    :: mxm1, mxm2
  real(kind=realkind),dimension(1:4)    :: dsdxi, dzdxi, dsdeta, dzdeta

  ! less memory
  !dzdeta(1) = DzDeta_over_J_sol(1,1,iel)
  !dzdeta(2) = DzDeta_over_J_sol(1,3,iel)
  !dzdeta(3) = DzDeta_over_J_sol(3,1,iel)
  !dzdeta(4) = DzDeta_over_J_sol(3,3,iel)

  !dzdxi(1) = DzDxi_over_J_sol(1,1,iel)
  !dzdxi(2) = DzDxi_over_J_sol(1,3,iel)
  !dzdxi(3) = DzDxi_over_J_sol(3,1,iel)
  !dzdxi(4) = DzDxi_over_J_sol(3,3,iel)

  !dsdeta(1) = DsDeta_over_J_sol(1,1,iel)
  !dsdeta(2) = DsDeta_over_J_sol(1,3,iel)
  !dsdeta(3) = DsDeta_over_J_sol(3,1,iel)
  !dsdeta(4) = DsDeta_over_J_sol(3,3,iel)

  !dsdxi(1) = DsDxi_over_J_sol(1,1,iel)
  !dsdxi(2) = DsDxi_over_J_sol(1,3,iel)
  !dsdxi(3) = DsDxi_over_J_sol(3,1,iel)
  !dsdxi(4) = DsDxi_over_J_sol(3,3,iel)
  
  ! 10% faster
  dzdeta(:) = DzDeta_over_J_sol_cg4(:,iel)
  dzdxi(:)  =  DzDxi_over_J_sol_cg4(:,iel)
  dsdeta(:) = DsDeta_over_J_sol_cg4(:,iel)
  dsdxi(:)  =  DsDxi_over_J_sol_cg4(:,iel)

  if (axis_solid(iel)) then 
     call mxm_cg4_sparse_c(G1T,f(:,:), mxm1) ! axial elements
  else 
     call mxm_cg4_sparse_c(G2T,f(:,:), mxm1) ! non-axial elements
  endif 

  call mxm_cg4_sparse_c(f(:,:), G2, mxm2)

  grad(:,1) = dzdeta * mxm1 + dzdxi * mxm2 ! dsdf
  grad(:,2) = dsdeta * mxm1 + dsdxi * mxm2 ! dzdf

end subroutine axisym_gradient_solid_el_cg4
!=============================================================================
 
!-------------------------------------------------------------
subroutine mxm_cg4_sparse_c(a,b,c)

  include "mesh_params.h" 

  real(kind=realkind), intent(in)  :: a(0:4,0:4), b(0:4,0:4)
  real(kind=realkind), intent(out) :: c(1:4)
  integer i,j

  c(1) = & 
     + a(1,0) * b(0,1) &
     + a(1,1) * b(1,1) &
     + a(1,2) * b(2,1) &
     + a(1,3) * b(3,1) &
     + a(1,4) * b(4,1)

  c(2) = & 
     + a(1,0) * b(0,3) &
     + a(1,1) * b(1,3) &
     + a(1,2) * b(2,3) &
     + a(1,3) * b(3,3) &
     + a(1,4) * b(4,3)

  c(3) = & 
     + a(3,0) * b(0,1) &
     + a(3,1) * b(1,1) &
     + a(3,2) * b(2,1) &
     + a(3,3) * b(3,1) &
     + a(3,4) * b(4,1)

  c(4) = & 
     + a(3,0) * b(0,3) &
     + a(3,1) * b(1,3) &
     + a(3,2) * b(2,3) &
     + a(3,3) * b(3,3) &
     + a(3,4) * b(4,3)

end subroutine mxm_cg4_sparse_c
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_gradient_solid_el(f,grad,iel)
  !
  ! Computes the axisymmetric gradient of scalar field f in the solid region:
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol, DzDxi_over_J_sol
  use data_pointwise, ONLY: DsDeta_over_J_sol, DsDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol)
  real(kind=realkind),intent(out)              :: grad(0:npol,0:npol,2)
  integer,intent(in)                           :: iel
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1, mxm2, dsdf, dzdf
  real(kind=realkind),dimension(0:npol,0:npol) :: dsdxi, dzdxi, dsdeta, dzdeta


  dzdeta = DzDeta_over_J_sol(:,:,iel)
  dzdxi  = DzDxi_over_J_sol(:,:,iel)
  dsdeta = DsDeta_over_J_sol(:,:,iel)
  dsdxi  = DsDxi_over_J_sol(:,:,iel)

  if (axis_solid(iel)) then 
     call mxm(G1T,f(:,:),mxm1) ! axial elements
  else 
     call mxm(G2T,f(:,:),mxm1) ! non-axial elements
  endif 

  call mxm(f(:,:),G2,mxm2)
  dsdf = dzdeta * mxm1 + dzdxi * mxm2
  dzdf = dsdeta * mxm1 + dsdxi * mxm2

  grad(:,:,1) = dsdf
  grad(:,:,2) = dzdf

end subroutine axisym_gradient_solid_el
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_gradient_solid(f,grad)
  !
  ! Computes the axisymmetric gradient of scalar field f in the solid region:
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol, DzDxi_over_J_sol
  use data_pointwise, ONLY: DsDeta_over_J_sol, DsDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind),intent(out)              :: grad(0:npol,0:npol,nel_solid,2)
  integer                                      :: iel
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
  real(kind=realkind),dimension(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta

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
     dsdf = dzdeta * mxm1 + dzdxi * mxm2
     dzdf = dsdeta * mxm1 + dsdxi * mxm2

     grad(:,:,iel,1) = dsdf
     grad(:,:,iel,2) = dzdf
  enddo

end subroutine axisym_gradient_solid
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_gradient_solid_add(f,grad)
  !
  ! Computes the axisymmetric gradient of scalar field f in the solid region:
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  ! This routine takes a previously calculated derivative and adds it
  ! to the result computed here in a permuted fashion.
  ! This saves the strain dump output two global fields, as the strain 
  ! trace will hereby be dumped as well as the entire E_31 term instead 
  ! of its two cross-derivative contributions.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
  use data_pointwise, ONLY: DsDeta_over_J_sol,DsDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)                    :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind),intent(inout)                 :: grad(0:npol,0:npol,nel_solid,2)
  integer                                           :: iel
  real(kind=realkind),dimension(0:npol,0:npol)      :: mxm1, mxm2, dsdf, dzdf
  real(kind=realkind),dimension(0:npol,0:npol)      :: dsdxi, dzdxi, dsdeta, dzdeta
  real(kind=realkind),dimension(0:npol,0:npol,2)    :: grad_old

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
    dsdf = dzdeta * mxm1 + dzdxi * mxm2
    dzdf = dsdeta * mxm1 + dsdxi * mxm2

    grad_old(0:npol,0:npol,1) = grad(0:npol,0:npol,iel,2)
    grad_old(0:npol,0:npol,2) = grad(0:npol,0:npol,iel,1)

    grad(0:npol,0:npol,iel,1) = grad_old(0:npol,0:npol,1) + dsdf(0:npol,0:npol)
    grad(0:npol,0:npol,iel,2) = grad_old(0:npol,0:npol,2) + dzdf(0:npol,0:npol)

 enddo

end subroutine axisym_gradient_solid_add
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_elem_solid(dsdf,f,iel)
  !
  ! Computes the elemental s-derivative of scalar field f in the solid region.
  ! This is used to compute the source term within the source element only.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol,DzDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)               :: f(0:npol,0:npol)
  real(kind=realkind), intent(out)              :: dsdf(0:npol,0:npol)
  integer,intent(in)                            :: iel
  real(kind=realkind), dimension(0:npol,0:npol) :: mxm1, mxm2
  real(kind=realkind), dimension(0:npol,0:npol) :: dzdxi, dzdeta
  
  dzdeta = DzDeta_over_J_sol(:,:,iel)
  dzdxi  = DzDxi_over_J_sol(:,:,iel)

  if (axis_solid(iel)) then 
     call mxm(G1T, f, mxm1) ! axial elements
  else 
     call mxm(G2T, f, mxm1) ! non-axial elements
  endif
  call mxm(f,G2,mxm2)

  dsdf = dzdeta * mxm1 + dzdxi * mxm2

end subroutine dsdf_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dzdf_elem_solid(dzdf,f,iel)
  !
  ! Computes the elemental z-derivative of scalar field f in the solid region.
  ! This is used to compute the source term within the source element only.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DsDeta_over_J_sol,DsDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)               :: f(0:npol,0:npol)
  real(kind=realkind), intent(out)              :: dzdf(0:npol,0:npol)
  integer,intent(in)                            :: iel
  real(kind=realkind), dimension(0:npol,0:npol) :: mxm1, mxm2
  real(kind=realkind), dimension(0:npol,0:npol) :: dsdxi, dsdeta
  
  dsdeta = DsDeta_over_J_sol(:,:,iel)
  dsdxi  = DsDxi_over_J_sol(:,:,iel)

  if (axis_solid(iel)) then 
     call mxm(G1T, f, mxm1) ! axial elements
  else 
     call mxm(G2T, f, mxm1) ! non-axial elements
  endif
  call mxm(f, G2, mxm2)
  dzdf = dsdeta * mxm1 + dsdxi * mxm2

end subroutine dzdf_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_solid_allaxis(f,dsdf)
  !
  ! Computes the pointwise derivative of scalar f in the s-direction 
  ! within the solid region, ONLY AT THE AXIS (needed for solid displacement)
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_sol, DzDxi_over_J_sol
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind),intent(out)              :: dsdf(0:npol,naxel_solid)
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1,mxm2
  real(kind=realkind),dimension(0:npol,0:npol) :: dzdxi,dzdeta,dsdf_el
  integer                                      :: ielem,iel

  do ielem=1, naxel_solid
    iel = ax_el_solid(ielem) 
    dzdeta = DzDeta_over_J_sol(:,:,iel)
    dzdxi  = DzDxi_over_J_sol(:,:,iel)
    call mxm(G1T, f(:,:,iel), mxm1) 
    call mxm(f(:,:,iel), G2, mxm2)
    dsdf_el = dzdeta * mxm1 + dzdxi * mxm2
    dsdf(:,ielem) = dsdf_el(0,:)
  enddo

end subroutine dsdf_solid_allaxis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine axisym_gradient_fluid(f,grad)
  !
  ! Computes the axisymmetric gradient of scalar field f in the fluid region:
  ! grad = \nabla {f}  = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
  use data_pointwise, ONLY: DsDeta_over_J_flu,DsDxi_over_J_flu
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol,nel_fluid)
  real(kind=realkind),intent(out)              :: grad(0:npol,0:npol,nel_fluid,2)
  integer                                      :: iel
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1,mxm2,dsdf,dzdf
  real(kind=realkind),dimension(0:npol,0:npol) :: dsdxi,dzdxi,dsdeta,dzdeta

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
    dsdf = dzdeta * mxm1 + dzdxi * mxm2
    dzdf = dsdeta * mxm1 + dsdxi * mxm2
    grad(:,:,iel,1) = dsdf
    grad(:,:,iel,2) = dzdf
 enddo

end subroutine axisym_gradient_fluid
!=============================================================================

!----------------------------------------------------------------------------
subroutine axisym_gradient_fluid_add(f,grad)
  !
  ! Computes the axisymmetric gradient of scalar field f in the fluid region:
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  ! This routine takes a previously calculated derivative and adds it
  ! to the result computed here in a permuted fashion.
  ! This saves the strain dump output two global fields, as the strain 
  ! trace will hereby be dumped as well as the entire E_31 term instead 
  ! of its two cross-derivative contributions.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_flu,DzDxi_over_J_flu
  use data_pointwise, ONLY: DsDeta_over_J_flu,DsDxi_over_J_flu
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)                 :: f(0:npol,0:npol,nel_fluid)
  real(kind=realkind), intent(inout)              :: grad(0:npol,0:npol,nel_fluid,2)
  integer                                         :: iel
  real(kind=realkind), dimension(0:npol,0:npol)   :: mxm1, mxm2, dsdf, dzdf
  real(kind=realkind), dimension(0:npol,0:npol)   :: dsdxi, dzdxi, dsdeta, dzdeta
  real(kind=realkind), dimension(0:npol,0:npol,2) :: grad_old

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
    dsdf = dzdeta * mxm1 + dzdxi * mxm2
    dzdf = dsdeta * mxm1 + dsdxi * mxm2

    grad_old(0:npol,0:npol,1) = grad(0:npol,0:npol,iel,2)
    grad_old(0:npol,0:npol,2) = grad(0:npol,0:npol,iel,1)

    grad(0:npol,0:npol,iel,1) = grad_old(0:npol,0:npol,1) + dsdf(0:npol,0:npol)
    grad(0:npol,0:npol,iel,2) = grad_old(0:npol,0:npol,2) + dzdf(0:npol,0:npol)

 enddo

end subroutine axisym_gradient_fluid_add
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_fluid_axis(f, iel, jpol, dsdf)
  !
  ! Computes the pointwise derivative of scalar f in the s-direction 
  ! within the fluid region, ONLY AT THE AXIS (needed for fluid displacement)
  ! and for a specific element iel and etsa coordinate index jpol.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_flu, DzDxi_over_J_flu
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  integer,intent(in)                           :: iel, jpol
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol)
  real(kind=realkind),intent(out)              :: dsdf
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1, mxm2
  real(kind=realkind),dimension(0:npol,0:npol) :: dzdxi, dzdeta, dsdf_el

  dzdeta = DzDeta_over_J_flu(:,:,iel)
  dzdxi  = DzDxi_over_J_flu(:,:,iel)
  call mxm(G1T, f, mxm1)
  call mxm(f, G2, mxm2)
  dsdf_el = dzdeta * mxm1 + dzdxi * mxm2
  dsdf = dsdf_el(0,jpol)

end subroutine dsdf_fluid_axis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dsdf_fluid_allaxis(f,dsdf)
  !
  ! Computes the pointwise derivative of scalar f in the s-direction 
  ! within the fluid region, ONLY AT THE AXIS (needed for fluid displacement)
  ! for all axial elements.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_pointwise, ONLY: DzDeta_over_J_flu, DzDxi_over_J_flu
  use unrolled_loops
  use unit_stride_colloc
  
  include 'mesh_params.h'
  
  real(kind=realkind),intent(in)               :: f(0:npol,0:npol,nel_fluid)
  real(kind=realkind),intent(out)              :: dsdf(0:npol,naxel_fluid)
  real(kind=realkind),dimension(0:npol,0:npol) :: mxm1, mxm2
  real(kind=realkind),dimension(0:npol,0:npol) :: dzdxi, dzdeta, dsdf_el
  integer                                      :: ielem, iel

  do ielem=1, naxel_fluid
    iel = ax_el_fluid(ielem) 
    dzdeta = DzDeta_over_J_flu(:,:,iel)
    dzdxi  = DzDxi_over_J_flu(:,:,iel)
    call mxm(G1T, f(:,:,iel), mxm1) 
    call mxm(f(:,:,iel), G2, mxm2)
    dsdf_el = dzdeta * mxm1 + dzdxi * mxm2
    dsdf(:,ielem) = dsdf_el(0,:)
  enddo

end subroutine dsdf_fluid_allaxis
!=============================================================================

!========================
end module pointwise_derivatives
!========================
