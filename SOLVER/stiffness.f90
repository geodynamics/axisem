!====================
module stiffness
!====================

  use global_parameters
  use data_matr
  use data_mesh, ONLY: axis_solid, axis_fluid, nsize, ani_true
  use data_spec
  use data_source
  
  use unrolled_loops
  use unit_stride_colloc
  
  implicit none
  
  public :: glob_stiffness_mono, glob_stiffness_di, glob_stiffness_quad
  public :: glob_anel_stiffness_mono
  public :: glob_fluid_stiffness
  private

contains


!-----------------------------------------------------------------------------
function outerprod(a,b) 
  ! outer product (dyadic) from numerical recipes
  real(kind=realkind), dimension(:), intent(in)     :: a, b
  real(kind=realkind), dimension(size(a),size(b))   :: outerprod

  outerprod = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
end function outerprod
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine glob_anel_stiffness_mono(glob_stiffness, R)

  use attenuation, ONLY: n_sls_attenuation
  include "mesh_params.h"
  
  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)
  real(kind=realkind), intent(in)    :: R(0:npol,0:npol,6,n_sls_attenuation,nel_solid)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  
  real(kind=realkind), dimension(0:npol,0:npol) :: r1, r2, r3, r5

  real(kind=realkind), dimension(0:npol,0:npol) :: yl
  real(kind=realkind), dimension(0:npol,0:npol) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:npol,0:npol) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s
  real(kind=realkind), dimension(0:npol,0:npol) :: S1z, S2z
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4
  
  real(kind=realkind), dimension(0:npol) :: y0l
  real(kind=realkind), dimension(0:npol) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:npol) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4
  
  integer :: ielem, j

  do ielem = 1, nel_solid

     loc_stiffness_s = 0
     loc_stiffness_z = 0

     yl(:,:) = Y(:,:,ielem)
     v_s_etal(:,:) = V_s_eta(:,:,ielem)
     v_s_xil(:,:)  = V_s_xi(:,:,ielem)
     v_z_etal(:,:) = V_z_eta(:,:,ielem)
     v_z_xil(:,:)  = V_z_xi(:,:,ielem)

     r1(:,:) = 0
     r2(:,:) = 0
     r3(:,:) = 0
     r5(:,:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:,:) = r1(:,:) + R(:,:,1,j,ielem)
        r2(:,:) = r2(:,:) + R(:,:,2,j,ielem)
        r3(:,:) = r3(:,:) + R(:,:,3,j,ielem)
        r5(:,:) = r5(:,:) + R(:,:,5,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5
     
     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if ( .not. axis_solid(ielem) ) then
        call mxm(G2,  S1s, X1)
        call mxm(G2,  S1z, X3)
     else
        call mxm(G1,  S1s, X1)
        call mxm(G1,  S1z, X3)
     endif

     call mxm(S2s, G2T, X2)
     call mxm(S2z, G2T, X4)

     loc_stiffness_s = X1 + X2 + yl * r2
     loc_stiffness_z = X3 + X4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! s - component
        V1 = v0_z_etal * r1(0,:) + v0_s_etal * r5(0,:) + y0l * r2(0,:)
        loc_stiffness_s = loc_stiffness_s + outerprod(G0, V1)

        ! z - component
        V2 = v0_z_etal * r5(0,:) + v0_s_etal * r3(0,:)
        V3 = v0_z_xil  * r5(0,:) + v0_s_xil  * r3(0,:)
        call vxm(V3, G2T, V4)

        loc_stiffness_z = loc_stiffness_z + outerprod(G0, V2)
        loc_stiffness_z(0,:) = loc_stiffness_z(0,:) + V4
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:npol,0:npol,ielem,1) = &
            glob_stiffness(0:npol,0:npol,ielem,1) - loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,3) = &
            glob_stiffness(0:npol,0:npol,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_mono
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_stiffness_mono(glob_stiffness,u)

  use global_parameters
  include "mesh_params.h"
  
  ! I/O global arrays
  real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,1:3)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: us, uz
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl, m12sl, m22sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m32sl, m42sl, m11zl, m21zl, m41zl
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4     ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s, S1z, S2z ! Sum arrays
  
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4
  real(kind=realkind), dimension(0:npol) :: uz0
  
  integer :: ielem

  glob_stiffness = zero
  
  do ielem = 1, nel_solid

     us(0:npol,0:npol) = u(0:npol,0:npol,ielem,1)
     uz(0:npol,0:npol) = u(0:npol,0:npol,ielem,3)

     m_1l(0:npol,0:npol) = M_1(:,:,ielem)
     m_2l(0:npol,0:npol) = M_2(:,:,ielem)
     m_3l(0:npol,0:npol) = M_3(:,:,ielem)
     m_4l(0:npol,0:npol) = M_4(:,:,ielem)
     
     m_w1l(0:npol,0:npol) = M_w1(:,:,ielem)

     m11sl(0:npol,0:npol) = M11s(:,:,ielem)
     m21sl(0:npol,0:npol) = M21s(:,:,ielem)
     m41sl(0:npol,0:npol) = M41s(:,:,ielem)
     m12sl(0:npol,0:npol) = M12s(:,:,ielem)
     m22sl(0:npol,0:npol) = M22s(:,:,ielem)
     m32sl(0:npol,0:npol) = M32s(:,:,ielem)
     m42sl(0:npol,0:npol) = M42s(:,:,ielem)
     m11zl(0:npol,0:npol) = M11z(:,:,ielem)
     m21zl(0:npol,0:npol) = M21z(:,:,ielem)
     m41zl(0:npol,0:npol) = M41z(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then
        call mxm(G2T, us, X1)
        call mxm(G2T, uz, X2)
     else 
        call mxm(G1T, us, X1)
        call mxm(G1T, uz, X2)
     endif

     call mxm(us, G2, X3)
     call mxm(uz, G2, X4)

     ! lower order terms in s
     loc_stiffness_s = m_4l * X4 + m_2l * X3 + m_1l * X1 + m_3l * X2 + us * m_w1l

     ! higher order terms + lower order terms with D_xi mxm ()
     S1s = m11sl * X3 + m21sl * X1 + m12sl * X4 + m22sl * X2 + m_1l * us
     S2s = m11sl * X1 + m41sl * X3 + m32sl * X2 + m42sl * X4 + m_2l * us
     S1z = m11zl * X4 + m21zl * X2 + m32sl * X3 + m22sl * X1 + m_3l * us
     S2z = m11zl * X2 + m41zl * X4 + m12sl * X1 + m42sl * X3 + m_4l * us
     
     call mxm(S2s, G2T, X2)
     call mxm(S2z, G2T, X4)

     if ( .not. axis_solid(ielem) ) then
        call mxm(G2, S1s, X1)
        call mxm(G2, S1z, X3)
     else 
        call mxm(G1, S1s, X1)
        call mxm(G1, S1z, X3)
     endif

     loc_stiffness_s = loc_stiffness_s + X1 + X2
     loc_stiffness_z = X3 + X4 

     ! additional axis terms
     if (axis_solid(ielem) ) then
        m0_w1l(0:npol) = M0_w1(0:npol,ielem)
        m0_w2l(0:npol) = M0_w2(0:npol,ielem)
        m0_w3l(0:npol) = M0_w3(0:npol,ielem)
        
        uz0 = uz(0,:)
        X2 = 0

        call vxm(G0, us, V1)
        call vxm(uz0, G2, V2)

        V4 = m0_w1l * V1 + m0_w3l * V2

        if (ani_true) then
           ! additional anisotropic terms
           call vxm(G0, uz, V3)

           V4 = V4 + m0_w2l * V3
           X2 = outerprod(G0, m0_w2l * V1)
        endif
           
        V2 = m0_w3l * V1

        call vxm(V2, G2T, V1)

        X2(0,:) = X2(0,:) + V1
                                
        loc_stiffness_s = loc_stiffness_s + outerprod(G0, V4)
        loc_stiffness_z = X2 + loc_stiffness_z

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_mono
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_stiffness_di(glob_stiffness,u)

  use global_parameters
  include "mesh_params.h"
  !use data_dipole
  
  ! I/O for global arrays
  real(kind=realkind),intent(in)  :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind),intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_1
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_2
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_3
  real(kind=realkind), dimension(0:npol,0:npol) :: u1,u2,u3
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_6l, m_2l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_5l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_4l, m_8l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_3l, m_7l
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m12sl, m22sl, m42sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m13sl, m23sl, m33sl, m43sl
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11zl, m21zl, m41zl
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l
  real(kind=realkind), dimension(0:npol) :: m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l
  
  integer :: ielem

  glob_stiffness = zero

  do ielem = 1, nel_solid

     loc_stiffness_1 = zero
     loc_stiffness_2 = zero
     loc_stiffness_3 = zero
     u1(0:npol,0:npol) = u(0:npol,0:npol,ielem,1)
     u2(0:npol,0:npol) = u(0:npol,0:npol,ielem,2)
     u3(0:npol,0:npol) = u(0:npol,0:npol,ielem,3)
     
     m_1l(0:npol,0:npol) = M_1(:,:,ielem)
     m_2l(0:npol,0:npol) = M_2(:,:,ielem)
     m_3l(0:npol,0:npol) = M_3(:,:,ielem)
     m_4l(0:npol,0:npol) = M_4(:,:,ielem)
     m_5l(0:npol,0:npol) = M_5(:,:,ielem)
     m_6l(0:npol,0:npol) = M_6(:,:,ielem)
     m_7l(0:npol,0:npol) = M_7(:,:,ielem)
     m_8l(0:npol,0:npol) = M_8(:,:,ielem)
     
     m_w1l(0:npol,0:npol) = M_w1(:,:,ielem)
     m_w2l(0:npol,0:npol) = M_w2(:,:,ielem)
     m_w3l(0:npol,0:npol) = M_w3(:,:,ielem)

     m11sl(0:npol,0:npol) = M11s(:,:,ielem)
     m21sl(0:npol,0:npol) = M21s(:,:,ielem)
     m41sl(0:npol,0:npol) = M41s(:,:,ielem)
     m12sl(0:npol,0:npol) = M12s(:,:,ielem)
     m22sl(0:npol,0:npol) = M22s(:,:,ielem)
     m42sl(0:npol,0:npol) = M42s(:,:,ielem)
     m13sl(0:npol,0:npol) = M13s(:,:,ielem)
     m23sl(0:npol,0:npol) = M32s(:,:,ielem) ! correct!! (static memory reasons,
                                            ! reusing static array from
                                            ! monopole)
     m33sl(0:npol,0:npol) = M33s(:,:,ielem)
     m43sl(0:npol,0:npol) = M43s(:,:,ielem)

     m11zl(0:npol,0:npol) = M11z(:,:,ielem)
     m21zl(0:npol,0:npol) = M21z(:,:,ielem)
     m41zl(0:npol,0:npol) = M41z(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then

        call stiffness_di_non_ax(u1, u2, u3, m_w1l, m_w2l, m_w3l, m_2l, m_6l, m_1l, &
            m_5l, m_4l, m_8l, m_3l, m_7l, m11sl, m21sl, m41sl, m12sl, m22sl, &
            m42sl, m13sl, m23sl, m33sl, m43sl, m11zl, m21zl, m41zl, &
            loc_stiffness_1, loc_stiffness_2, loc_stiffness_3)

     else
        m0_w1l(0:npol) = M0_w1(0:npol,ielem)
        m0_w2l(0:npol) = M0_w2(0:npol,ielem)
        m0_w3l(0:npol) = M0_w3(0:npol,ielem)
        m0_w4l(0:npol) = M0_w4(0:npol,ielem)
        m0_w5l(0:npol) = M0_w5(0:npol,ielem)
        m0_w6l(0:npol) = M0_w6(0:npol,ielem)
        m0_w7l(0:npol) = M0_w7(0:npol,ielem)
        m0_w8l(0:npol) = M0_w8(0:npol,ielem)
        m0_w9l(0:npol) = M0_w9(0:npol,ielem)
        m0_w10l(0:npol) = M0_w10(0:npol,ielem)

        call stiffness_di_ax(u1, u2, u3, m_w1l, m_w2l, m_w3l, m_2l, m_6l, m_1l, m_5l, &
            m_4l, m_8l,  m_3l, m_7l, m11sl, m21sl, m41sl, m12sl, m22sl, m42sl, &
            m13sl, m23sl, m33sl, m43sl, m11zl, m21zl, m41zl, m0_w1l, m0_w2l, &
            m0_w3l, m0_w4l, m0_w5l, m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l, &
            loc_stiffness_1, loc_stiffness_2, loc_stiffness_3)

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_1
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_2
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_3

  enddo

end subroutine glob_stiffness_di
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_stiffness_quad(glob_stiffness,u)

  use global_parameters
  include "mesh_params.h"
  !use data_quadrupole
  
  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,1:3)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_phi
  real(kind=realkind), dimension(0:npol,0:npol) :: us,uz,uphi
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_5l, m_6l, m_7l, m_8l
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l, m_w4l, m_w5l
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl,m21sl,m41sl,m12sl,m22sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m32sl,m42sl,m11zl,m21zl,m41zl
  real(kind=realkind), dimension(0:npol,0:npol) :: m1phil,m2phil,m4phil
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol) :: m0_w4l, m0_w5l, m0_w6l
  
  integer :: ielem

  glob_stiffness = zero
  
  do ielem = 1, nel_solid

     loc_stiffness_s = zero
     loc_stiffness_phi = zero
     loc_stiffness_z = zero
     us(0:npol,0:npol)   = u(0:npol,0:npol,ielem,1)
     uphi(0:npol,0:npol) = u(0:npol,0:npol,ielem,2)
     uz(0:npol,0:npol)   = u(0:npol,0:npol,ielem,3)
    
     m_1l(0:npol,0:npol) = M_1(:,:,ielem)
     m_2l(0:npol,0:npol) = M_2(:,:,ielem)
     m_3l(0:npol,0:npol) = M_3(:,:,ielem)
     m_4l(0:npol,0:npol) = M_4(:,:,ielem)
     m_5l(0:npol,0:npol) = M_5(:,:,ielem)
     m_6l(0:npol,0:npol) = M_6(:,:,ielem)
     m_7l(0:npol,0:npol) = M_7(:,:,ielem)
     m_8l(0:npol,0:npol) = M_8(:,:,ielem)
     
     m_w1l(0:npol,0:npol) = M_w1(:,:,ielem)
     m_w2l(0:npol,0:npol) = M_w2(:,:,ielem)
     m_w3l(0:npol,0:npol) = M_w3(:,:,ielem)
     m_w4l(0:npol,0:npol) = M_w4(:,:,ielem)
     m_w5l(0:npol,0:npol) = M_w5(:,:,ielem)

     m11sl(0:npol,0:npol) = M11s(:,:,ielem)
     m21sl(0:npol,0:npol) = M21s(:,:,ielem)
     m41sl(0:npol,0:npol) = M41s(:,:,ielem)
     m12sl(0:npol,0:npol) = M12s(:,:,ielem)
     m22sl(0:npol,0:npol) = M22s(:,:,ielem)
     m32sl(0:npol,0:npol) = M32s(:,:,ielem)
     m42sl(0:npol,0:npol) = M42s(:,:,ielem)
     m11zl(0:npol,0:npol) = M11z(:,:,ielem)
     m21zl(0:npol,0:npol) = M21z(:,:,ielem)
     m41zl(0:npol,0:npol) = M41z(:,:,ielem)
     m1phil(0:npol,0:npol) = M1phi(:,:,ielem)
     m2phil(0:npol,0:npol) = M2phi(:,:,ielem)
     m4phil(0:npol,0:npol) = M4phi(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then

        call stiffness_quad_non_ax(us,uphi,uz, &
             m_1l, m_2l, m_3l, m_4l, m_5l, m_6l, m_7l, m_8l, &
             m_w1l, m_w2l, m_w3l, m_w4l, m_w5l, &
             m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
             m11zl,m21zl,m41zl,m1phil,m2phil,m4phil, &
             loc_stiffness_s,loc_stiffness_phi,loc_stiffness_z)
     else

        m0_w1l(0:npol) = M0_w1(0:npol,ielem)
        m0_w2l(0:npol) = M0_w2(0:npol,ielem)
        m0_w3l(0:npol) = M0_w3(0:npol,ielem)
        m0_w4l(0:npol) = M0_w4(0:npol,ielem)
        m0_w5l(0:npol) = M0_w5(0:npol,ielem)
        m0_w6l(0:npol) = M0_w6(0:npol,ielem)

        call stiffness_quad_ax(us,uphi,uz, &
             m_1l, m_2l, m_3l, m_4l, m_5l, m_6l, m_7l, m_8l, &
             m_w1l, m_w2l, m_w3l, m_w4l, m_w5l, &
             m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
             m11zl,m21zl,m41zl,m1phil,m2phil,m4phil, &
             m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l, m0_w6l, &
             loc_stiffness_s,loc_stiffness_phi,loc_stiffness_z)

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_phi
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_quad
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_fluid_stiffness(glob_stiffness_fl,chi)

  include "mesh_params.h"
  
  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: chi(0:npol,0:npol,nel_fluid)
  real(kind=realkind), intent(out) :: glob_stiffness_fl(0:npol,0:npol,nel_fluid)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: chi_l,loc_stiffness
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w_fl_l
  real(kind=realkind), dimension(0:npol,0:npol) :: m1chil,m2chil,m4chil
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w_fl_l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol) :: V1
  
  integer :: ielem

  glob_stiffness_fl = zero
  
  do ielem = 1, nel_fluid

     loc_stiffness = zero
     chi_l(0:npol,0:npol)=chi(0:npol,0:npol,ielem)
     m1chil(0:npol,0:npol)=M1chi_fl(:,:,ielem)
     m2chil(0:npol,0:npol)=M2chi_fl(:,:,ielem)
     m4chil(0:npol,0:npol)=M4chi_fl(:,:,ielem)

     call fluid_stiffness_joint_part(chi_l,m1chil,m2chil,m4chil, &
          loc_stiffness)

     ! dipole and quadrupole cases: additional 2nd order term
     if (src_type(1) .ne. 'monopole') then

        m_w_fl_l(0:npol,0:npol)=M_w_fl(:,:,ielem)
        ! collocate and add this to previous term

        call collocate_sum_existent_1d(m_w_fl_l,chi_l,loc_stiffness,nsize)
        if ( axis_fluid(ielem) ) then
           m0_w_fl_l(0:npol)=M0_w_fl(0:npol,ielem)
           call vxm(G0,chi_l,V1)
           call collocate_tensor_1d(m0_w_fl_l,V1,G0,chi_l,npol) !chi_l as dummy 
           call sum2s_1d(chi_l,loc_stiffness,nsize) ! sum, add to existent
        endif

     endif

     glob_stiffness_fl(0:npol,0:npol,ielem) = loc_stiffness

  enddo

end subroutine glob_fluid_stiffness
!=============================================================================

!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!     D I P O L E   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine stiffness_di_non_ax(u1, u2, u3, m_w1l, m_w2l, m_w3l, m_2l, m_6l, m_1l, m_5l, &
        m_4l, m_8l, m_3l, m_7l, m11sl, m21sl, m41sl, m12sl, m22sl, m42sl, m13sl, &
        m23sl, m33sl, m43sl, m11zl, m21zl, m41zl, loc_stiffness_1, &
        loc_stiffness_2, loc_stiffness_3)
  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! displacement (+, - and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1, u2, u3
  
  ! precomputed matrices containing weights, mapping, elastic parameters
  ! precomputed matrices for W_x
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w1l, m_w2l, m_w3l
  
  ! precomputed matrices for W_x^d, D_x^y
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_2l, m_6l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_1l, m_5l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_4l, m_8l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_3l, m_7l
  
  ! precomputed matrices for D_xy^xy: + and - components
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m11sl, m21sl, m41sl
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m12sl, m22sl, m42sl
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m13sl, m23sl, m33sl, m43sl
  
  ! precomputed matrices for D_xy^xy: z component
  real(kind=realkind), dimension(0:npol,0:npol) :: m11zl, m21zl, m41zl
  
  ! result: local (elemental) stiffness term (+, - and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_1
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_2
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_3
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2, loc_stiffness_s3
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6, X7, X8 ! MxM 
  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S1m, S2p, S2m, S1z, S2z ! Sum
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_1 = zero
  loc_stiffness_2 = zero
  loc_stiffness_3 = zero

  S1p = zero
  S2p = zero
  S1m = zero
  S2m = zero
  S1z = zero
  S2z = zero
  
  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero
  X7 = zero
  X8 = zero

  ! First MxM
  call mxm(G2T, u1, X1)
  call mxm(G2T, u2, X2)
  call mxm(G2T, u3, X3)
  call mxm(u1, G2, X4)
  call mxm(u2, G2, X5)
  call mxm(u3, G2, X6)

  ! Sum for the z-component
  call sum2_2_1d(X1, X2, X7, X4, X5, X8, nsize)

  ! Collocations and sums of W_x and W_x^d terms

  ! - component
  call collocate8_sum_1d(m_8l, X6, m_7l, X3, m_1l, X1, &
                         m_5l, X2, m_2l, X4, m_6l, X5, &
                         m_w1l, u2, m_w2l, u3, loc_stiffness_s2, nsize)

  ! z component
  call collocate6z_sum_1d(m_4l, X4, m_4l, X5, m_3l, X1, &
                          m_3l, X2, m_w2l, u2, m_w3l, u3, loc_stiffness_s3, nsize)

  ! Collocations and sums of D terms
  call collocate28s_sum_1d(m11sl, m21sl, m41sl, m12sl, m22sl, m42sl, &
                          m13sl, m23sl, m33sl, m43sl,&
                          X1, X2, X3, X4, X5, X6, u2, u3, &
                          m_1l, m_3l, m_5l, &
                          m_2l, m_4l, m_6l, &
                          S1p, S2p, S1m, S2m, nsize)

  call collocate5_sum_1d(m33sl, X8, m23sl, X7, m11zl, X6, m21zl, X3, &
                         m_7l, u2, S1z, nsize)
  call collocate5_sum_1d(m13sl, X7, m43sl, X8, m11zl, X3, m41zl, X6, &
                         m_8l, u2, S2z, nsize)

  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero

  !Second MxM
  call mxm(G2, S1p, X1)
  call mxm(S2p, G2T, X2)
  call mxm(G2, S1m, X3)
  call mxm(S2m, G2T, X4)
  call mxm(G2, S1z, X5)
  call mxm(S2z, G2T, X6)

  call sum2s_3_3_1d(X1, X2, loc_stiffness_1, X3, X4, loc_stiffness_s2, &
                   loc_stiffness_2, X5, X6, loc_stiffness_s3, &
                   loc_stiffness_3, nsize)

end subroutine stiffness_di_non_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine stiffness_di_ax(u1, u2, u3, m_w1l, m_w2l, m_w3l, m_2l, m_6l, m_1l, m_5l, &
        m_4l, m_8l, m_3l, m_7l, m11sl, m21sl, m41sl, m12sl, m22sl, m42sl, m13sl, &
        m23sl, m33sl, m43sl, m11zl, m21zl, m41zl, m0_w1l, m0_w2l, m0_w3l, m0_w4l, &
        m0_w5l, m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l, loc_stiffness_1, &
        loc_stiffness_2, loc_stiffness_3)

  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! displacement (+, - and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1, u2, u3
  
  ! precomputed matrices containing weights, mapping, elastic parameters
  ! precomputed matrices for W_x
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w1l, m_w2l, m_w3l
  
  ! precomputed matrices for W_x^d, D_x^y
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_2l, m_6l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_1l, m_5l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_4l, m_8l
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m_3l, m_7l
  
  ! precomputed matrices for D_xy^xy: + and - components
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m11sl, m21sl, m41sl
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m12sl, m22sl, m42sl
  real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m13sl, m23sl, m33sl, m43sl
  
  ! precomputed matrices for D_xy^xy: z component
  real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m11zl, m21zl, m41zl
  
  ! result: local (elemental) stiffness term (+, - and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_1
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_2
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_3
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol),intent(in) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol),intent(in) :: m0_w4l, m0_w5l
  real(kind=realkind), dimension(0:npol),intent(in) :: m0_w6l, m0_w7l, m0_w8l
  real(kind=realkind), dimension(0:npol),intent(in) :: m0_w9l, m0_w10l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2, loc_stiffness_s3 
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6, X7, X8 ! MxM 
  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S1m, S2p, S2m, S1z, S2z ! Sum
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_1 = zero
  loc_stiffness_2 = zero
  loc_stiffness_3 = zero

  S1p = zero
  S2p = zero
  S1m = zero
  S2m = zero
  S1z = zero
  S2z = zero

  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero
  X7 = zero
  X8 = zero

  ! First MxM
  call mxm(G1T, u1, X1)
  call mxm(G1T, u2, X2)
  call mxm(G1T, u3, X3)
  call mxm(u1, G2, X4)
  call mxm(u2, G2, X5)
  call mxm(u3, G2, X6)

  ! Sum for the z-component
  call sum2_2_1d(X1, X2, X7, X4, X5, X8, nsize)

  ! Collocations and sums of W_x and W_x^d terms
  ! - component
  call collocate8_sum_1d(m_8l, X6, m_7l, X3, m_1l, X1, &
                         m_5l, X2, m_2l, X4, m_6l, X5, &
                         m_w1l, u2, m_w2l, u3, loc_stiffness_s2, nsize)
  ! z component
  call collocate6z_sum_1d(m_4l, X4, m_4l, X5, m_3l, X1, &
                          m_3l, X2, m_w2l, u2, m_w3l, u3, loc_stiffness_s3, nsize)

  ! Collocations and sums of D terms
  ! + and - component
  call collocate28s_sum_1d(m11sl, m21sl, m41sl, m12sl, m22sl, m42sl, &
                          m13sl, m23sl, m33sl, m43sl,&
                          X1, X2, X3, X4, X5, X6, U2, U3, &
                          m_1l, m_3l, m_5l, &
                          m_2l, m_4l, m_6l,&
                          S1p, S2p, S1m, S2m, nsize)
  ! z component
  call collocate5_sum_1d(m33sl, X8, m23sl, X7, m11zl, X6, m21zl, X3, &
                         m_7l, u2, S1z, nsize)
  call collocate5_sum_1d(m13sl, X7, m43sl, X8, m11zl, X3, m41zl, X6, &
                         m_8l, u2, S2z, nsize)

  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero

  ! Second MxM
  call mxm(G1, S1p, X1)
  call mxm(S2p, G2T, X2)
  call mxm(G1, S1m, X3)
  call mxm(S2m, G2T, X4)
  call mxm(G1, S1z, X5)
  call mxm(S2z, G2T, X6)

  S1p = zero
  S1m = zero
  S1z = zero

  ! Additional terms for the axial elements
  if (ani_true) then
    call additional_di_ax_ani(m0_w1l, m0_w2l, m0_w3l, m0_w3l, m0_w4l, m0_w5l, &
                              m0_w7l, m0_w8l, m0_w9l, m0_w10l,&
                              u1, u2, u3, S1p, S1m, S1z)
  else
    call additional_di_ax(m0_w1l, m0_w2l, m0_w7l, m0_w8l, m0_w9l, &
                          u1, u2, u3, S1p, S1m, S1z)
  endif
  

  call sum3s_4_4_1d(X1, X2, S1p, loc_stiffness_1, X3, X4, S1m, loc_stiffness_s2, &
                   loc_stiffness_2, X5, X6, S1z, loc_stiffness_s3, &
                   loc_stiffness_3, nsize)

end subroutine stiffness_di_ax
!=============================================================================


!-----------------------------------------------------------------------------
subroutine additional_di_ax(m0_w1l, m0_w2l, m0_w7l, m0_w8l, m0_w9l, &
                            u1, u2, u3, S1p, S1m, S1z)
include "mesh_params.h"

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w7l, m0_w8l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w1l, m0_w2l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w9l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1, u2, u3
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1p, S1m, S1z
  
  ! work arrays  
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4, V5
  real(kind=realkind), dimension(0:npol) :: u10, u20
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1 = zero
  V2 = zero
  V3 = zero
  V4 = zero

  S1p = zero
  S1m = zero
  S1z = zero

  u10 = u1(0,:)
  u20 = u2(0,:)

  ! VxM
  call vxm(G0, u1, V1)
  call vxm(G0, u2, V2)
  call vxm(G0, u3, V3)

  call vxm(u10, G2, V4)
  call vxm(u20, G2, V5)

  ! Collocations

  ! + comp
  call collocate_tensor_1d(m0_w1l, V2, G0, S1p, npol)

  ! - comp
  call collocate3_sum_tensor_1d(m0_w1l, V1, m0_w2l, V5, m0_w9l, V2, G0, S1m, npol)

  ! z comp
  call collocate2_sum_tensor_1d(m0_w7l, V3, m0_w8l, V4, G0, S1z, npol)

  ! Final VxM D_z^z in + component
  call collocate02_sum_1d(m0_w2l, V2, m0_w8l, V3, V4, npol) 
  call vxm(V4, G2T, V1)
  call add_to_axis_2d(V1, S1p, npol)

end subroutine additional_di_ax
!=============================================================================


!-----------------------------------------------------------------------------
subroutine additional_di_ax_ani(m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l, &
                                m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l, &
                                u1, u2, u3, S1p, S1m, S1z)
  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w1l, m0_w2l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w3l, m0_w5l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w4l, m0_w6l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w7l, m0_w8l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w9l, m0_w10l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1, u2, u3
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1p, S1m, S1z
  
  ! work arrays  
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4, V5
  real(kind=realkind), dimension(0:npol) :: u10, u20
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1 = zero
  V2 = zero
  V3 = zero
  V4 = zero

  S1p = zero
  S1m = zero
  S1z = zero

  u10 = u1(0,:)
  u20 = u2(0,:)

  ! VxM
  call vxm(G0, u1, V1)
  call vxm(G0, u2, V2)
  call vxm(G0, u3, V3)

  call vxm(u10, G2, V4)
  call vxm(u20, G2, V5)

  ! Collocations

  ! + comp
  call collocate2_sum_tensor_1d(m0_w1l, V2, m0_w3l, V3, G0, S1p, npol)

  ! - comp
  call collocate5_sum_tensor_1d(m0_w1l, V1, m0_w2l, V5, m0_w6l, V4, &
                                m0_w9l, V2, m0_w10l, V3, G0, S1m, npol)

  ! z comp
  call collocate5_sum_tensor_1d(m0_w3l, V1, m0_w4l, V4, m0_w7l, V3, &
                                m0_w8l, V4, m0_w10l, V2, G0, S1z, npol)

  ! Final VxM D_z^z in + component
  call collocate02p_sum_1d(m0_w2l, m0_w6l, V2, m0_w4l, m0_w8l, V3, V4, npol) 

  call vxm(V4, G2T, V1)
  call add_to_axis_2d(V1, S1p, npol)

end subroutine additional_di_ax_ani
!=============================================================================


!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!    Q U A D R U P O L E   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine stiffness_quad_non_ax(us, uphi, uz,  &
           m_1l,  m_2l,  m_3l,  m_4l,  m_5l,  m_6l,  m_7l,  m_8l,  &
           m_w1l, m_w2l, m_w3l, m_w4l, m_w5l, &
           m11sl, m21sl, m41sl, m12sl, m22sl, m32sl, m42sl,  &
           m11zl, m21zl, m41zl, m1phil, m2phil, m4phil, &
           loc_stiffness_s, loc_stiffness_phi, loc_stiffness_z)

include "mesh_params.h"

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! displacement (s, phi and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz,uphi
  
  ! precomputed matrices containing weights, mapping, elastic parameters
  ! precomputed matrix for W_s
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w1l, m_w2l, m_w3l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w4l, m_w5l
  
  ! precomputed matrices for W_s^d, D_x^y
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_5l, m_6l, m_7l, m_8l
  
  ! precomputed matrices for D_xy^xy
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11sl,m21sl,m41sl 
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m12sl,m22sl,m32sl,m42sl
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11zl,m21zl,m41zl
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m1phil,m2phil,m4phil
  
  ! result: local (elemental) stiffness term (s, phi and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_phi
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_z
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4,X5,X6 ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s,S2s,S1phi,S2phi,S1z,S2z ! Sum
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_s = zero
  loc_stiffness_phi = zero
  loc_stiffness_z = zero

  S1s = zero
  S2s = zero
  S1phi = zero
  S2phi = zero
  S1z = zero
  S2z = zero

  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero

  !++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
  ! First MxM
  call mxm(G2T, us, X1)
  call mxm(G2T, uphi, X2)
  call mxm(G2T, uz, X3)
  call mxm(us, G2, X4)
  call mxm(uphi, G2, X5)
  call mxm(uz, G2, X6)

  ! Collocations and sums of W terms: s and phi components
  call collocate12s_sum_1d(m_2l, X4, m_1l, X1, m_6l, X5, &
                           m_5l, X2, m_4l, X6, m_3l, X3, &
                           m_w1l, m_w2l, m_w3l, m_w4l, us, uphi, uz, &
                           loc_stiffness_s, loc_stiffness_phi, nsize)

  ! Collocations and sums of W terms: z component
  call collocate5z_sum_1d(m_8l, X5, m_7l, X2, m_w3l, m_w5l, us, uphi, uz, &
                         loc_stiffness_z, nsize)

  !++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
  ! Collocations and sums of D terms

  ! s component
  call collocate6s_sum_1d(m11sl, X4, m21sl, X1, m12sl, X6, m22sl, X3, &
                          m_1l, us, uphi, S1s, nsize)
  call collocate6s_sum_1d(m11sl, X1, m41sl, X4, m32sl, X3, m42sl, X6, &
                          m_2l, us, uphi, S2s, nsize)
  ! z component
  call collocate6s_sum_1d(m11zl, X6, m21zl, X3, m32sl, X4, m22sl, X1, &
                          m_3l, us, uphi, S1z, nsize)
  call collocate6s_sum_1d(m11zl, X3, m41zl, X6, m12sl, X1, m42sl, X4, &
                          m_4l, us, uphi, S2z, nsize)

  ! phi component
  call collocate5ss_sum_1d(m1phil, X5, m2phil, X2, m_5l, us, &
                           m_5l, uphi, m_7l, uz, S1phi, nsize)
  call collocate5ss_sum_1d(m1phil, X2, m4phil, X5, m_6l, us, &
                           m_6l, uphi, m_8l, uz, S2phi, nsize)

  !Second MxM
  call mxm(G2, S1s, X1)
  call mxm(S2s, G2T, X2)
  call mxm(G2, S1phi, X3)
  call mxm(S2phi, G2T, X4)
  call mxm(G2, S1z, X5)
  call mxm(S2z, G2T, X6)

  ! Final Sum
  call sum3s_3_1d(X1, X2, loc_stiffness_s, X3, X4, loc_stiffness_phi, X5, X6, &
                 loc_stiffness_z, nsize)

end subroutine stiffness_quad_non_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine stiffness_quad_ax(us, uphi, uz,  &
           m_1l,  m_2l,  m_3l,  m_4l,  m_5l,  m_6l,  m_7l,  m_8l,  &
           m_w1l, m_w2l, m_w3l, m_w4l, m_w5l, &
           m11sl, m21sl, m41sl, m12sl, m22sl, m32sl, m42sl,  &
           m11zl, m21zl, m41zl, m1phil, m2phil, m4phil,  &
           m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l, m0_w6l, &
           loc_stiffness_s, loc_stiffness_phi, loc_stiffness_z)

include "mesh_params.h"

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! displacement (s, phi and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz,uphi
  
  ! precomputed matrices containing weights, mapping, elastic parameters
  ! precomputed matrix for W_s
  real(kind=realkind), dimension(0:npol,0:npol), intent(in)  :: m_w1l, m_w2l, m_w3l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in)  :: m_w4l, m_w5l
  
  ! precomputed matrices for W_s^d, D_x^y
  real(kind=realkind), dimension(0:npol,0:npol), intent(in)  :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in)  :: m_5l, m_6l, m_7l, m_8l
  
  ! precomputed matrices for D_xy^xy
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11sl,m21sl,m41sl
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m12sl,m22sl,m32sl,m42sl
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11zl,m21zl,m41zl
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m1phil,m2phil,m4phil
  
  ! axial terms
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w4l, m0_w5l, m0_w6l
  
  ! result: local (elemental) stiffness term (s, phi and z components)
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_phi
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_z
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4,X5,X6 ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s,S2s,S1phi,S2phi,S1z,S2z ! Sum
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_s = zero
  loc_stiffness_phi = zero
  loc_stiffness_z = zero
  
  S1s = zero
  S2s = zero
  S1phi = zero
  S2phi = zero
  S1z = zero
  S2z = zero
  
  X1 = zero
  X2 = zero
  X3 = zero
  X4 = zero
  X5 = zero
  X6 = zero

  !++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
  ! First MxM
  call mxm(G1T, us, X1)
  call mxm(G1T, uphi, X2)
  call mxm(G1T, uz, X3)
  call mxm(us, G2, X4)
  call mxm(uphi, G2, X5)
  call mxm(uz, G2, X6)

  ! Collocations and sums of W terms: s and phi components
  call collocate12s_sum_1d(m_2l, X4, m_1l, X1, m_6l, X5, &
                           m_5l, X2, m_4l, X6, m_3l, X3, &
                           m_w1l, m_w2l, m_w3l, m_w4l, us, uphi, uz, &
                           loc_stiffness_s, loc_stiffness_phi, nsize)

  ! Collocations and sums of W terms: z component
  call collocate5z_sum_1d(m_8l, X5, m_7l, X2, m_w3l, m_w5l, us, uphi, uz, &
                         loc_stiffness_z, nsize)

  !++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
  ! Collocations and sums of D terms

  ! s component
  call collocate6s_sum_1d(m11sl, X4, m21sl, X1, m12sl, X6, m22sl, X3, &
                          m_1l, us, uphi, S1s, nsize)
  call collocate6s_sum_1d(m11sl, X1, m41sl, X4, m32sl, X3, m42sl, X6, &
                          m_2l, us, uphi, S2s, nsize)
  ! z component
  call collocate6s_sum_1d(m11zl, X6, m21zl, X3, m32sl, X4, m22sl, X1, &
                          m_3l, us, uphi, S1z, nsize)
  call collocate6s_sum_1d(m11zl, X3, m41zl, X6, m12sl, X1, m42sl, X4, &
                          m_4l, us, uphi, S2z, nsize)

  ! phi component
  call collocate5ss_sum_1d(m1phil, X5, m2phil, X2, m_5l, us, &
                           m_5l, uphi, m_7l, uz, S1phi, nsize)
  call collocate5ss_sum_1d(m1phil, X2, m4phil, X5, m_6l, us, &
                           m_6l, uphi, m_8l, uz, S2phi, nsize)

  !Second MxM
  call mxm(G1, S1s, X1)
  call mxm(S2s, G2T, X2)
  call mxm(G1, S1phi, X3)
  call mxm(S2phi, G2T, X4)
  call mxm(G1, S1z, X5)
  call mxm(S2z, G2T, X6)

  if (ani_true) then
    call additional_quad_ax_ani(m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l, m0_w6l, &
                                us, uphi, uz, S1s, S1phi, S1z)
  else
    call additional_quad_ax(m0_w1l, m0_w4l, m0_w6l, m0_w2l,  &
                            us, uphi, uz, S1s, S1phi, S1z)
  endif

  ! Final Sum
  call sum4_3_1d(X1, X2, S1s, loc_stiffness_s, X3, X4, S1phi, loc_stiffness_phi, X5,  &
                 X6, S1z, loc_stiffness_z, nsize)

end subroutine stiffness_quad_ax
!=============================================================================


!-----------------------------------------------------------------------------
subroutine additional_quad_ax(m0_z_eta_2l_w_l_6ml, m0_z_eta_3m_w_4l_9ml, &
                              m0_w_4mul, m0_z_eta_2lm_min_w_2l_6ml,  &
                              us, uphi, uz, S1s, S1phi, S1z)

  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_2l_w_l_6ml
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_3m_w_4l_9ml
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w_4mul
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_2lm_min_w_2l_6ml
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uphi,uz
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1s,S1phi,S1z
  
  ! work arrays  
  real(kind=realkind), dimension(0:npol) :: V1,V2,V3
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1 = zero
  V2 = zero
  V3 = zero
  S1s = zero
  S1phi = zero
  S1z = zero

  ! VxM
  call vxm(G0, us, V1)
  call vxm(G0, uphi, V2)
  call vxm(G0, uz, V3)

  ! Collocations, Sums, Tensorization

  ! s comp
  call collocate2_sum_tensor_1d(m0_z_eta_2l_w_l_6ml, V1, &
                                m0_z_eta_2lm_min_w_2l_6ml, V2, G0, S1s, npol)
                                
  ! phi comp
  call collocate2_sum_tensor_1d(m0_z_eta_2lm_min_w_2l_6ml, V1,  &
                                m0_z_eta_3m_w_4l_9ml, V2, G0, S1phi, npol)
                                
  ! z comp
  call collocate_tensor_1d(m0_w_4mul, V3, G0, S1z, npol)

end subroutine additional_quad_ax
!=============================================================================


!-----------------------------------------------------------------------------
subroutine additional_quad_ax_ani(m0_w1l, m0_w2l, m0_w3l, m0_w4l, m0_w5l, m0_w6l, &
                                  us, uphi, uz, S1s, S1phi, S1z)

include "mesh_params.h"

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w4l, m0_w5l, m0_w6l
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us, uphi, uz
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1s, S1phi, S1z
  
  ! work arrays  
  real(kind=realkind), dimension(0:npol) :: V1,V2,V3
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1 = zero
  V2 = zero
  V3 = zero
  S1s = zero
  S1phi = zero
  S1z = zero

  ! VxM
  call vxm(G0, us, V1)
  call vxm(G0, uphi, V2)
  call vxm(G0, uz, V3)

  ! Collocations, Sums, Tensorization

  ! s comp
  call collocate3_sum_tensor_1d(m0_w1l, V1, m0_w2l, V2, m0_w3l, V3, G0, S1s, npol)
                                
  ! phi comp
  call collocate3_sum_tensor_1d(m0_w2l, V1, m0_w4l, V2, m0_w5l, V3, G0, S1phi, npol)
                                
  ! z comp
  call collocate3_sum_tensor_1d(m0_w3l, V1, m0_w5l, V2, m0_w6l, V3, G0, S1z, npol)

end subroutine additional_quad_ax_ani
!=============================================================================


!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!      F L U I D   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine fluid_stiffness_joint_part(chi,m1chil,m2chil,m4chil,loc_stiffness)
           
  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! pressure scalar potential
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: chi
  
  ! precomputed matrices for D_xy^xy
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m1chil,m2chil,m4chil
  
  ! result: local (elemental) scalar stiffness term
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2  ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1,S2  ! Sum
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness=zero
  S1=zero;S2=zero;S1=zero;X1=zero;X2=zero

  !++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
  ! First MxM
  call mxm(G2T,chi,X1)
  call mxm(chi,G2,X2)

  ! Collocations and sums of D terms
  call collocate2_sum_1d(m1chil,X2,m2chil,X1,S1,nsize)
  call collocate2_sum_1d(m1chil,X1,m4chil,X2,S2,nsize)

  !Second MxM
  call mxm(G2,S1,X1)
  call mxm(S2,G2T,X2)
  
  ! Final Sum
  call sum_1d(X1,X2,loc_stiffness,nsize)

end subroutine fluid_stiffness_joint_part
!=============================================================================

!-----------------------------------------------------------------------------
subroutine fluid_stiffness_add_multi(chi,m_w_fl_l,loc_stiffness)

  include "mesh_params.h"
  
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  ! pressure scalar potential
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: chi
  
  ! precomputed matrices containing weights, mapping, elastic parameters
  ! precomputed matrix for W_s
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_fl_l
  
  ! result: local (elemental) scalar stiffness term
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness=zero

  ! dipole and quadrupole cases: collocation of W
  call collocate0_1d(m_w_fl_l,chi,loc_stiffness,nsize)

end subroutine fluid_stiffness_add_multi
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!====================
end module stiffness
!====================
