!====================
module stiffness
!====================

use global_parameters
use data_matr
use data_mesh, ONLY: axis_solid, axis_fluid, nsize
use data_spec
use data_source

use unrolled_loops
use unit_stride_colloc

implicit none

public :: glob_stiffness_mono,glob_stiffness_di,glob_stiffness_quad
public :: glob_fluid_stiffness
private

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine glob_stiffness_mono(glob_stiffness,u)

use data_monopole

! I/O global arrays
real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,1:3)
real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)

! local variables for all elements
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
real(kind=realkind), dimension(0:npol,0:npol) :: us,uz
real(kind=realkind), dimension(0:npol,0:npol) :: m_wl
real(kind=realkind), dimension(0:npol,0:npol) :: ms_xil,ms_etal,mz_xil,mz_etal
real(kind=realkind), dimension(0:npol,0:npol) :: m11sl,m21sl,m41sl,m12sl,m22sl
real(kind=realkind), dimension(0:npol,0:npol) :: m32sl,m42sl,m11zl,m21zl,m41zl

! local variables for axial elements
real(kind=realkind), dimension(0:npol) :: m0_wl
real(kind=realkind), dimension(0:npol) :: m0s_xil

integer :: ielem

  glob_stiffness=zero
  
  do ielem = 1, nel_solid

     loc_stiffness_s=zero; loc_stiffness_z=zero
     us(0:npol,0:npol)=u(0:npol,0:npol,ielem,1)
     uz(0:npol,0:npol)=u(0:npol,0:npol,ielem,3)
     m_wl(0:npol,0:npol)=M_w(:,:,ielem)
     ms_xil(0:npol,0:npol)=Ms_xi(:,:,ielem)
     ms_etal(0:npol,0:npol)=Ms_eta(:,:,ielem)
     mz_xil(0:npol,0:npol)=Mz_xi(:,:,ielem)
     mz_etal(0:npol,0:npol)=Mz_eta(:,:,ielem)
     m11sl(0:npol,0:npol)=M11s(:,:,ielem)
     m21sl(0:npol,0:npol)=M21s(:,:,ielem)
     m41sl(0:npol,0:npol)=M41s(:,:,ielem)
     m12sl(0:npol,0:npol)=M12s(:,:,ielem)
     m22sl(0:npol,0:npol)=M22s(:,:,ielem)
     m32sl(0:npol,0:npol)=M32s(:,:,ielem)
     m42sl(0:npol,0:npol)=M42s(:,:,ielem)
     m11zl(0:npol,0:npol)=M11z(:,:,ielem)
     m21zl(0:npol,0:npol)=M21z(:,:,ielem)
     m41zl(0:npol,0:npol)=M41z(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then

        call stiffness_mono_non_ax(us,uz, &
             m_wl,ms_xil,mz_xil,mz_etal,ms_etal, &
             m11sl,m21sl,m41sl,m12sl,m22sl,m32sl, &
             m42sl,m11zl,m21zl,m41zl, &
             loc_stiffness_s,loc_stiffness_z)

     else

        m0s_xil(0:npol)=M0_s_xi(0:npol,ielem)
        m0_wl(0:npol)=M0_w(0:npol,ielem)

        call stiffness_mono_ax(us,uz, &
             m_wl,ms_xil,mz_xil,mz_etal,ms_etal, &
             m0_wl,m0s_xil, &
             m11sl,m21sl,m41sl,m12sl,m22sl, &
             m32sl,m42sl,m11zl,m21zl,m41zl, &
             loc_stiffness_s,loc_stiffness_z)

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_mono
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_stiffness_di(glob_stiffness,u)

use data_dipole

! I/O for global arrays
real(kind=realkind),intent(in)  :: u(0:npol,0:npol,nel_solid,3)
real(kind=realkind),intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

! local variables for all elements
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_1
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_2
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_3
real(kind=realkind), dimension(0:npol,0:npol) :: u1,u2,u3
real(kind=realkind), dimension(0:npol,0:npol) :: m_w_4_lam_mul,m_w_mul
real(kind=realkind), dimension(0:npol,0:npol) :: mz_xi_pl,mz_xi_ml
real(kind=realkind), dimension(0:npol,0:npol) :: mz_eta_pl,mz_eta_ml
real(kind=realkind), dimension(0:npol,0:npol) :: ms_xi_mul,ms_xi_2laml
real(kind=realkind), dimension(0:npol,0:npol) :: ms_eta_mul,ms_eta_2laml

real(kind=realkind), dimension(0:npol,0:npol) :: m11sl,m21sl,m41sl
real(kind=realkind), dimension(0:npol,0:npol) :: m12sl,m22sl,m42sl
real(kind=realkind), dimension(0:npol,0:npol) :: m13sl,m23sl,m33sl,m43sl

real(kind=realkind), dimension(0:npol,0:npol) :: m11zl,m21zl,m41zl

! local variables for axial elements
real(kind=realkind), dimension(0:npol) :: m0_w_mul
real(kind=realkind), dimension(0:npol) :: m0s_xi_mul
real(kind=realkind), dimension(0:npol) :: m0z_eta_2_lam_mul
real(kind=realkind), dimension(0:npol) :: m0z_eta_4lm_w_l3ml

integer :: ielem

  glob_stiffness=zero

  do ielem = 1, nel_solid

     loc_stiffness_1=zero; loc_stiffness_2=zero; loc_stiffness_3=zero;
     u1(0:npol,0:npol)=u(0:npol,0:npol,ielem,1)
     u2(0:npol,0:npol)=u(0:npol,0:npol,ielem,2)
     u3(0:npol,0:npol)=u(0:npol,0:npol,ielem,3)
     mz_xi_pl(0:npol,0:npol)=Mz_xi_p(:,:,ielem)
     mz_xi_ml(0:npol,0:npol)=Mz_xi_m(:,:,ielem)
     mz_eta_pl(0:npol,0:npol)=Mz_eta_p(:,:,ielem)
     mz_eta_ml(0:npol,0:npol)=Mz_eta_m(:,:,ielem)
     ms_xi_mul(0:npol,0:npol)=Ms_xi_mu(:,:,ielem)
     ms_xi_2laml(0:npol,0:npol)=Ms_xi_2lam(:,:,ielem)
     ms_eta_mul(0:npol,0:npol)=Ms_eta_mu(:,:,ielem)
     ms_eta_2laml(0:npol,0:npol)=Ms_eta_2lam(:,:,ielem)
     m_w_4_lam_mul(0:npol,0:npol)=M_w_4_lam_mu(:,:,ielem)
     m_w_mul(0:npol,0:npol)=M_w_mu(:,:,ielem)

     m11sl(0:npol,0:npol)=M11s(:,:,ielem)
     m21sl(0:npol,0:npol)=M21s(:,:,ielem)
     m41sl(0:npol,0:npol)=M41s(:,:,ielem)
     m12sl(0:npol,0:npol)=M12s(:,:,ielem)
     m22sl(0:npol,0:npol)=M22s(:,:,ielem)
     m42sl(0:npol,0:npol)=M42s(:,:,ielem)
     m13sl(0:npol,0:npol)=M13s(:,:,ielem)
     m23sl(0:npol,0:npol)=M32s(:,:,ielem) ! correct!! (static memory reasons)
     m33sl(0:npol,0:npol)=M33s(:,:,ielem)
     m43sl(0:npol,0:npol)=M43s(:,:,ielem)

     m11zl(0:npol,0:npol)=M11z(:,:,ielem)
     m21zl(0:npol,0:npol)=M21z(:,:,ielem)
     m41zl(0:npol,0:npol)=M41z(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then

        call stiffness_di_non_ax(u1,u2,u3, &
             m_w_4_lam_mul,m_w_mul,mz_xi_pl,mz_xi_ml,&
             mz_eta_pl,mz_eta_ml,ms_xi_mul,ms_xi_2laml,&
             ms_eta_mul,ms_eta_2laml,&
             m11sl,m21sl,m41sl,m12sl,m22sl,m42sl,m13sl,m23sl,m33sl,m43sl,&
             m11zl,m21zl,m41zl, &
             loc_stiffness_1,loc_stiffness_2,loc_stiffness_3)

     else
        m0s_xi_mul(0:npol)=M0_s_xi_mu(0:npol,ielem)
        m0z_eta_2_lam_mul(0:npol)=M0_z_eta_2_lam_mu(0:npol,ielem)
        m0_w_mul(0:npol)=M0_w_mu(0:npol,ielem)
        m0z_eta_4lm_w_l3ml(0:npol)=M0_z_eta_4lm_w_l3m(0:npol,ielem)

        call stiffness_di_ax(u1,u2,u3, &
             m_w_4_lam_mul,m_w_mul,mz_xi_pl,mz_xi_ml,&
             mz_eta_pl,mz_eta_ml,ms_xi_mul,ms_xi_2laml,&
             ms_eta_mul,ms_eta_2laml,&
             m11sl,m21sl,m41sl,m12sl,m22sl,m42sl,m13sl,m23sl,m33sl,m43sl,&
             m11zl,m21zl,m41zl,&
             m0s_xi_mul,m0z_eta_2_lam_mul,m0_w_mul,m0z_eta_4lm_w_l3ml,&
             loc_stiffness_1,loc_stiffness_2,loc_stiffness_3)

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_1
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_2
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_3

  enddo

end subroutine glob_stiffness_di
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_stiffness_quad(glob_stiffness,u)

use data_quadrupole

! I/O for global arrays
real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,1:3)
real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)

! local variables for all elements
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_phi
real(kind=realkind), dimension(0:npol,0:npol) :: us,uz,uphi

real(kind=realkind), dimension(0:npol,0:npol) :: mz_xi_laml,ms_xi_laml
real(kind=realkind), dimension(0:npol,0:npol) :: mz_eta_laml,ms_eta_laml
real(kind=realkind), dimension(0:npol,0:npol) :: mz_xi_2mul,ms_xi_2mul
real(kind=realkind), dimension(0:npol,0:npol) :: mz_eta_2mul,ms_eta_2mul
real(kind=realkind), dimension(0:npol,0:npol) :: m_w_lam_6mul,m_w_min_2lam_6mul
real(kind=realkind), dimension(0:npol,0:npol) :: m_w_4lam_9mul,m_w_4mul
real(kind=realkind), dimension(0:npol,0:npol) :: m11sl,m21sl,m41sl,m12sl,m22sl
real(kind=realkind), dimension(0:npol,0:npol) :: m32sl,m42sl,m11zl,m21zl,m41zl
real(kind=realkind), dimension(0:npol,0:npol) :: m1phil,m2phil,m4phil

! local variables for axial elements
real(kind=realkind), dimension(0:npol) :: m0_w_4mul
real(kind=realkind), dimension(0:npol) :: m0_z_eta_3m_w_4l_9ml
real(kind=realkind), dimension(0:npol) :: m0_z_eta_2l_w_l_6ml
real(kind=realkind), dimension(0:npol) :: m0_z_eta_2lm_min_w_2l_6ml

integer :: ielem

  glob_stiffness=zero
  
  do ielem = 1, nel_solid

     loc_stiffness_s=zero
     loc_stiffness_phi=zero
     loc_stiffness_z=zero
     us(0:npol,0:npol)   = u(0:npol,0:npol,ielem,1)
     uphi(0:npol,0:npol) = u(0:npol,0:npol,ielem,2)
     uz(0:npol,0:npol)   = u(0:npol,0:npol,ielem,3)
     mz_xi_laml(0:npol,0:npol)=Mz_xi_lam(:,:,ielem)
     mz_eta_laml(0:npol,0:npol)=Mz_eta_lam(:,:,ielem)
     ms_xi_laml(0:npol,0:npol)=Ms_xi_lam(:,:,ielem)
     ms_eta_laml(0:npol,0:npol)=Ms_eta_lam(:,:,ielem)
     mz_xi_2mul(0:npol,0:npol)=Mz_xi_2mu(:,:,ielem)
     mz_eta_2mul(0:npol,0:npol)=Mz_eta_2mu(:,:,ielem)
     ms_xi_2mul(0:npol,0:npol)=Ms_xi_2mu(:,:,ielem)
     ms_eta_2mul(0:npol,0:npol)=Ms_eta_2mu(:,:,ielem)
     m_w_lam_6mul(0:npol,0:npol)=M_w_lam_6mu(:,:,ielem)
     m_w_min_2lam_6mul(0:npol,0:npol)=M_w_min_2lam_6mu(:,:,ielem)
     m_w_4lam_9mul(0:npol,0:npol)=M_w_4lam_9mu(:,:,ielem)
     m_w_4mul(0:npol,0:npol)=M_w_4mu(:,:,ielem)
     m11sl(0:npol,0:npol)=M11s(:,:,ielem)
     m21sl(0:npol,0:npol)=M21s(:,:,ielem)
     m41sl(0:npol,0:npol)=M41s(:,:,ielem)
     m12sl(0:npol,0:npol)=M12s(:,:,ielem)
     m22sl(0:npol,0:npol)=M22s(:,:,ielem)
     m32sl(0:npol,0:npol)=M32s(:,:,ielem)
     m42sl(0:npol,0:npol)=M42s(:,:,ielem)
     m11zl(0:npol,0:npol)=M11z(:,:,ielem)
     m21zl(0:npol,0:npol)=M21z(:,:,ielem)
     m41zl(0:npol,0:npol)=M41z(:,:,ielem)
     m1phil(0:npol,0:npol)=M1phi(:,:,ielem)
     m2phil(0:npol,0:npol)=M2phi(:,:,ielem)
     m4phil(0:npol,0:npol)=M4phi(:,:,ielem)

     if ( .not. axis_solid(ielem) ) then

        call stiffness_quad_non_ax(us,uphi,uz, &
             mz_xi_laml,ms_xi_laml,mz_eta_laml,ms_eta_laml, &
             mz_xi_2mul,ms_xi_2mul,mz_eta_2mul,ms_eta_2mul, &
             m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,m_w_4mul, &
             m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
             m11zl,m21zl,m41zl,m1phil,m2phil,m4phil, &
             loc_stiffness_s,loc_stiffness_phi,loc_stiffness_z)
     else

        m0_w_4mul(0:npol)=M0_w_4mu(0:npol,ielem)
        m0_z_eta_2l_w_l_6ml(0:npol)=M0_z_eta_2l_w_l_6m(0:npol,ielem)
        m0_z_eta_3m_w_4l_9ml(0:npol)=M0_z_eta_3m_w_4l_9m(0:npol,ielem)
        m0_z_eta_2lm_min_w_2l_6ml(0:npol)=M0_z_eta_2lm_min_w_2l_6m(0:npol,ielem)

        call stiffness_quad_ax(us,uphi,uz, &
             mz_xi_laml,ms_xi_laml,mz_eta_laml,ms_eta_laml, &
             mz_xi_2mul,ms_xi_2mul,mz_eta_2mul,ms_eta_2mul, &
             m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,m_w_4mul, &
             m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
             m11zl,m21zl,m41zl,m1phil,m2phil,m4phil, &
             m0_w_4mul,m0_z_eta_2l_w_l_6ml, &
             m0_z_eta_3m_w_4l_9ml, m0_z_eta_2lm_min_w_2l_6ml, &
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

  glob_stiffness_fl=zero
  
  do ielem = 1, nel_fluid

     loc_stiffness=zero
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
!    M O N O P O L E   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine stiffness_mono_non_ax(us,uz,m_w,ms_xi,mz_xi, &
                                 mz_eta,ms_eta,m11s,m21s,m41s,m12s,m22s, &
                                 m32s,m42s,m11z,m21z,m41z,&
                                 loc_stiffness_s,loc_stiffness_z)
include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (s and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrix for W_s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w 

! precomputed matrices for W_s^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi,mz_xi   
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_eta,ms_eta 

! precomputed matrices for D_xy^xy
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11s,m21s,m41s 
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m12s,m22s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m32s,m42s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11z,m21z,m41z

! result: local (elemental) stiffness term (s and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_s
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_z

! work arrays
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s1  
real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4      ! MxM arrays
real(kind=realkind), dimension(0:npol,0:npol) :: S1s,S2s,S1z,S2z  ! Sum arrays
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_s=zero
  loc_stiffness_s1=zero
  loc_stiffness_z=zero
  S1s=zero;S2s=zero;S1z=zero;S2z=zero
  X1=zero;X2=zero;X3=zero;X4=zero

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! First MxM
  call mxm(G2T,us,X1)
  call mxm(G2T,uz,X2)
  call mxm(us,G2,X3)
  call mxm(uz,G2,X4)

! Collocations and sums of W_x^d terms
  call collocate4_sum_1d(ms_xi,X4,mz_xi,X3,mz_eta,X1,ms_eta,X2,S1s,nsize)

!++++++++++++++++++MONOPOLE SOURCE+++++++++++++++++++++++++++++++++++++++++
! Collocation of W_s term and sum to s-comp (specific to monopole case)
  call collocate_sum_1d(us,m_w,S1s,loc_stiffness_s1,nsize)

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! Collocations and sums of D terms

  call collocate5_sum_1d(m11s,X3,m21s,X1,m12s,X4,m22s,X2,mz_eta,us,S1s,nsize)
  call collocate5_sum_1d(m11s,X1,m41s,X3,m32s,X2,m42s,X4,mz_xi ,us,S2s,nsize)
  call collocate5_sum_1d(m11z,X4,m21z,X2,m32s,X3,m22s,X1,ms_eta,us,S1z,nsize)
  call collocate5_sum_1d(m11z,X2,m41z,X4,m12s,X1,m42s,X3,ms_xi ,us,S2z,nsize)

!Second MxM
  call mxm(G2,S1s,X1)
  call mxm(S2s,G2T,X2)
  call mxm(G2,S1z,X3)
  call mxm(S2z,G2T,X4)

! Final Sum
  call sum2_3_1d(X3,X4,loc_stiffness_z,loc_stiffness_s1,X1,X2, &
                 loc_stiffness_s,nsize)

end subroutine stiffness_mono_non_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine stiffness_mono_ax(us,uz, &
                             m_w,ms_xi,mz_xi,mz_eta,ms_eta, &
                             m0_w,m0s_xi, &
                             m11s,m21s,m41s,m12s,m22s, &
                             m32s,m42s,m11z,m21z,m41z,&
                             loc_stiffness_s,loc_stiffness_z)

include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (s and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrices for W_s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w 
real(kind=realkind), dimension(0:npol), intent(in) :: m0_w
real(kind=realkind), dimension(0:npol), intent(in) :: m0s_xi

! precomputed matrices for W_s^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi,mz_xi   
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_eta,ms_eta 

! precomputed matrices for D_xy^xy
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11s,m21s,m41s 
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m12s,m22s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m32s,m42s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11z,m21z,m41z

! result: local (elemental) stiffness term (s and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_s
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_z
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s1
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z1

! work arrays
real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4     ! MxM arrays
real(kind=realkind), dimension(0:npol,0:npol) :: S1s,S2s,S1z,S2z ! Sum arrays
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_s=zero;  loc_stiffness_z=zero
  loc_stiffness_s1=zero; loc_stiffness_z1=zero
  X1=zero;X2=zero;X3=zero;X4=zero
  S1s=zero;S2s=zero;S1z=zero;S2z=zero
 
!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! First MxM
  call mxm(G1T,us,X1)
  call mxm(G1T,uz,X2)
  call mxm(us,G2,X3)
  call mxm(uz,G2,X4)

! Collocations and sums of W_x^d terms
  call collocate4_sum_1d(ms_xi,X4,mz_xi,X3,mz_eta,X1,ms_eta,X2,S1s,nsize)

!++++++++++++++++++MONOPOLE SOURCE+++++++++++++++++++++++++++++++++++++++++
! Collocation of W_s term and sum to s-comp (specific to monopole case)
  call collocate_sum_1d(us,m_w,S1s,loc_stiffness_s,nsize)

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! Collocations and sums of D terms
  call collocate5_sum_1d(m11s,X3,m21s,X1,m12s,X4,m22s,X2,mz_eta,us,S1s,nsize)
  call collocate5_sum_1d(m11s,X1,m41s,X3,m32s,X2,m42s,X4,mz_xi ,us,S2s,nsize)
  call collocate5_sum_1d(m11z,X4,m21z,X2,m32s,X3,m22s,X1,ms_eta,us,S1z,nsize)
  call collocate5_sum_1d(m11z,X2,m41z,X4,m12s,X1,m42s,X3,ms_xi ,us,S2z,nsize)

!Second MxM
  call mxm(G1,S1s,X1)
  call mxm(S2s,G2T,X2)
  call mxm(G1,S1z,X3)
  call mxm(S2z,G2T,X4)

! Sum
  call sum2_3_1d(X3,X4,loc_stiffness_z1,loc_stiffness_s, &
                 X1,X2,loc_stiffness_s1,nsize)

!++++++++++++++++++MONOPOLE SOURCE+++++++++++++++++++++++++++++++++++++++++
! additional terms for the axial elements
  call additional_mono_ax(m0_w,m0s_xi,us,uz,X1,X2)
                          
 
! Final Sum::::::::::::::::::::::::::::::::::::::::
  call sum2_2_1d(X1,loc_stiffness_s1,loc_stiffness_s,X2, &
                 loc_stiffness_z1,loc_stiffness_z,nsize)

end subroutine stiffness_mono_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine additional_mono_ax(m0_w,m0s_xi,us,uz,S1s,S1z)
                       
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! I/O
include "mesh_params.h"
  real(kind=realkind), dimension(0:npol), intent(in) :: m0_w
  real(kind=realkind), dimension(0:npol), intent(in) :: m0s_xi
  real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz
  real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1s,S1z

! work arrays
  real(kind=realkind), dimension(0:npol) :: V1,V2
  real(kind=realkind), dimension(0:npol) :: uz0

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1=zero;V2=zero;S1s=zero;S1z=zero
  uz0 = uz(0,:)

! VxM
  call vxm(G0,us,V1)
  call vxm(uz0,G2,V2)

! Collocations for s-component
  call collocate2_sum_tensor_1d(m0_w,V1,m0s_xi,V2,G0,S1s,npol)

! Collocation for z-component 
  call collocate0_1d(m0s_xi,V1,V2,npol)

! Final VxM D_z^z for z-component
  call vxm(V2,G2T,V1)
  S1z(0,:)=V1

end subroutine additional_mono_ax
!=============================================================================


!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!     D I P O L E   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine stiffness_di_non_ax(u1,u2,u3, &
           m_w_4_lam_mul,m_w_mul,mz_xi_pl,mz_xi_ml,&
           mz_eta_pl,mz_eta_ml,ms_xi_mul,ms_xi_2laml,&
           ms_eta_mul,ms_eta_2laml,&
           m11sl,m21sl,m41sl,m12sl,m22sl,m42sl,m13sl,m23sl,m33sl,m43sl,&
           m11zl,m21zl,m41zl, &
           loc_stiffness_1,loc_stiffness_2,loc_stiffness_3)

include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (+, - and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1,u2,u3

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrices for W_x
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4_lam_mul,m_w_mul

! precomputed matrices for W_x^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: mz_xi_pl,mz_xi_ml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: mz_eta_pl,mz_eta_ml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: ms_xi_mul,ms_xi_2laml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: ms_eta_mul,ms_eta_2laml

! precomputed matrices for D_xy^xy: + and - components
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m11sl,m21sl,m41sl
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m12sl,m22sl,m42sl
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m13sl,m23sl,m33sl,m43sl

! precomputed matrices for D_xy^xy: z component
real(kind=realkind), dimension(0:npol,0:npol) :: m11zl,m21zl,m41zl

! result: local (elemental) stiffness term (+, - and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_1
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_2
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_3

! work arrays
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2, loc_stiffness_s3
real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4,X5,X6,X7,X8 ! MxM 
real(kind=realkind), dimension(0:npol,0:npol) :: S1p,S1m,S2p,S2m,S1z,S2z ! Sum
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_1=zero;loc_stiffness_2=zero;loc_stiffness_3=zero
  S1p=zero;S2p=zero;S1m=zero;S2m=zero;S1z=zero;S2z=zero
  X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero;X7=zero;X8=zero

! First MxM
  call mxm(G2T,u1,X1)
  call mxm(G2T,u2,X2)
  call mxm(G2T,u3,X3)
  call mxm(u1,G2,X4)
  call mxm(u2,G2,X5)
  call mxm(u3,G2,X6)

! Sum for the z-component
  call sum2_2_1d(X1,X2,X7,X4,X5,X8,nsize)

! Collocations and sums of W_x and W_x^d terms

! - component
  call collocate7_sum_1d(ms_xi_2laml,X6,ms_eta_2laml,X3,mz_eta_pl,X1, &
                         mz_eta_ml,X2,mz_xi_pl,X4,mz_xi_ml,X5, &
                         m_w_4_lam_mul,u2,loc_stiffness_s2,nsize)


! z component
  call collocate5s_sum_1d(ms_xi_mul,X4,ms_xi_mul,X5,ms_eta_mul,X1, &
                          ms_eta_mul,X2,m_w_mul,u3,loc_stiffness_s3,nsize)

!write(3333,*)'locstiff:',loc_stiffness_s2(3,3)/20.d0,loc_stiffness_s3(3,3), &
!             loc_stiffness_s2(3,3)/20.d0/loc_stiffness_s3(3,3)
!write(3333,*)'M:',mz_eta_pl(3,3),mz_xi_pl(3,3),mz_eta_ml(3,3),mz_xi_ml(3,3), &
!                  m_w_4_lam_mul(3,3)

! Collocations and sums of D terms
  call collocate28s_sum_1d(m11sl,m21sl,m41sl,m12sl,m22sl,m42sl, &
                          m13sl,m23sl,m33sl,m43sl,&
                          X1,X2,X3,X4,X5,X6,u2,u3, &
                          mz_eta_pl,ms_eta_mul,mz_eta_ml, &
                          mz_xi_pl,ms_xi_mul,mz_xi_ml, &
                          S1p,S2p,S1m,S2m,nsize)

!write(3333,*)'S:',S1p(5,2),S1m(5,2),S2p(5,2),S2m(5,2)


  call collocate5_sum_1d(m33sl,X8,m23sl,X7,m11zl,X6,m21zl,X3, &
                         ms_eta_2laml,u2,S1z,nsize)
  call collocate5_sum_1d(m13sl,X7,m43sl,X8,m11zl,X3,m41zl,X6, &
                         ms_xi_2laml,u2,S2z,nsize)

X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero

!Second MxM
  call mxm(G2,S1p,X1)
  call mxm(S2p,G2T,X2)
  call mxm(G2,S1m,X3)
  call mxm(S2m,G2T,X4)
  call mxm(G2,S1z,X5)
  call mxm(S2z,G2T,X6)

! Final Sum
!  call sum2_3_3_1d(X1,X2,loc_stiffness_1,X3,X4,loc_stiffness_2,X5,X6, &
!                   loc_stiffness_3,nsize)

  call sum2s_3_3_1d(X1,X2,loc_stiffness_1,X3,X4,loc_stiffness_s2, &
                   loc_stiffness_2,X5,X6,loc_stiffness_s3, &
                   loc_stiffness_3,nsize)

end subroutine stiffness_di_non_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine stiffness_di_ax(u1,u2,u3, &
           m_w_4_lam_mul,m_w_mul,mz_xi_pl,mz_xi_ml,&
           mz_eta_pl,mz_eta_ml,ms_xi_mul,ms_xi_2laml,&
           ms_eta_mul,ms_eta_2laml,&
           m11sl,m21sl,m41sl,m12sl,m22sl,m42sl,m13sl,m23sl,m33sl,m43sl,&
           m11zl,m21zl,m41zl,&
           m0s_xi_mul,m0z_eta_2_lam_mul,m0_w_mul,m0z_eta_4lm_w_l3ml,&
           loc_stiffness_1,loc_stiffness_2,loc_stiffness_3)

include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (+, - and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1,u2,u3

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrices for W_x
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4_lam_mul,m_w_mul

! precomputed matrices for W_x^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: mz_xi_pl,mz_xi_ml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: mz_eta_pl,mz_eta_ml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: ms_xi_mul,ms_xi_2laml
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: ms_eta_mul,ms_eta_2laml

! precomputed matrices for D_xy^xy: + and - components
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m11sl,m21sl,m41sl
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m12sl,m22sl,m42sl
real(kind=realkind),dimension(0:npol,0:npol),intent(in):: m13sl,m23sl,m33sl,m43sl

! precomputed matrices for D_xy^xy: z component
real(kind=realkind), dimension(0:npol,0:npol),intent(in) :: m11zl,m21zl,m41zl

! result: local (elemental) stiffness term (+, - and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_1
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_2
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_3

! local variables for axial elements
real(kind=realkind), dimension(0:npol),intent(in) :: m0_w_mul
real(kind=realkind), dimension(0:npol),intent(in) :: m0s_xi_mul
real(kind=realkind), dimension(0:npol),intent(in) :: m0z_eta_2_lam_mul
real(kind=realkind), dimension(0:npol),intent(in) :: m0z_eta_4lm_w_l3ml

! work arrays
real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2,loc_stiffness_s3 
real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4,X5,X6,X7,X8 ! MxM 
real(kind=realkind), dimension(0:npol,0:npol) :: S1p,S1m,S2p,S2m,S1z,S2z ! Sum
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_1=zero;loc_stiffness_2=zero;loc_stiffness_3=zero
  S1p=zero;S2p=zero;S1m=zero;S2m=zero;S1z=zero;S2z=zero
  X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero;X7=zero;X8=zero

! First MxM
  call mxm(G1T,u1,X1)
  call mxm(G1T,u2,X2)
  call mxm(G1T,u3,X3)
  call mxm(u1,G2,X4)
  call mxm(u2,G2,X5)
  call mxm(u3,G2,X6)

! Sum for the z-component
  call sum2_2_1d(X1,X2,X7,X4,X5,X8,nsize)

! Collocations and sums of W_x and W_x^d terms
! - component
  call collocate7_sum_1d(ms_xi_2laml,X6,ms_eta_2laml,X3,mz_eta_pl,X1, &
                         mz_eta_ml,X2,mz_xi_pl,X4,mz_xi_ml,X5, &
                         m_w_4_lam_mul,u2,loc_stiffness_s2,nsize)
! z component
  call collocate5s_sum_1d(ms_xi_mul,X4,ms_xi_mul,X5,ms_eta_mul,X1, &
                        ms_eta_mul,X2,m_w_mul,u3,loc_stiffness_s3,nsize)

! Collocations and sums of D terms
! + and - component
  call collocate28s_sum_1d(m11sl,m21sl,m41sl,m12sl,m22sl,m42sl, &
                          m13sl,m23sl,m33sl,m43sl,&
                          X1,X2,X3,X4,X5,X6,U2,U3, &
                          mz_eta_pl,ms_eta_mul,mz_eta_ml, &
                          mz_xi_pl,ms_xi_mul,mz_xi_ml,&
                          S1p,S2p,S1m,S2m,nsize)
! z component
  call collocate5_sum_1d(m33sl,X8,m23sl,X7,m11zl,X6,m21zl,X3, &
                         ms_eta_2laml,u2,S1z,nsize)
  call collocate5_sum_1d(m13sl,X7,m43sl,X8,m11zl,X3,m41zl,X6, &
                         ms_xi_2laml,u2,S2z,nsize)

X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero

! Second MxM
  call mxm(G1,S1p,X1)
  call mxm(S2p,G2T,X2)
  call mxm(G1,S1m,X3)
  call mxm(S2m,G2T,X4)
  call mxm(G1,S1z,X5)
  call mxm(S2z,G2T,X6)

S1p=zero;S1m=zero;S1z=zero

! Additional terms for the axial elements
  call additional_di_ax(m0_w_mul,m0s_xi_mul, &
                        m0z_eta_2_lam_mul,m0z_eta_4lm_w_l3ml, &
                        u1,u2,u3,S1p,S1m,S1z)

! Final Sum
!  call sum3_4_4_1d(X1,X2,S1p,loc_stiffness_1,X3,X4,S1m,loc_stiffness_2, &
!                   X5,X6,S1z,loc_stiffness_3,nsize)

  call sum3s_4_4_1d(X1,X2,S1p,loc_stiffness_1,X3,X4,S1m,loc_stiffness_s2, &
                   loc_stiffness_2,X5,X6,S1z,loc_stiffness_s3, &
                   loc_stiffness_3,nsize)

end subroutine stiffness_di_ax
!=============================================================================


!-----------------------------------------------------------------------------
subroutine additional_di_ax(m0_w_mul,m0s_xi_mul, &
                            m0z_eta_2_lam_mul,m0z_eta_4lm_w_l3ml, &
                            u1,u2,u3,S1p,S1m,S1z)
include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
real(kind=realkind), dimension(0:npol), intent(in) :: m0_w_mul,m0s_xi_mul
real(kind=realkind), dimension(0:npol), intent(in) :: m0z_eta_2_lam_mul
real(kind=realkind), dimension(0:npol), intent(in) :: m0z_eta_4lm_w_l3ml
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: u1,u2,u3
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: S1p,S1m,S1z

! work arrays  
real(kind=realkind), dimension(0:npol) :: V1,V2,V3,V4
real(kind=realkind), dimension(0:npol) :: u10
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  V1=zero;V2=zero;V3=zero;V4=zero
  S1p=zero;S1m=zero;S1z=zero

  u10 = u1(0,:)

! VxM
  call vxm(G0,u1,V1)
  call vxm(G0,u2,V2)
  call vxm(G0,u3,V3)
!
  call vxm(u10,G2,V4)

! Collocations

! + comp
  call collocate_tensor_1d(m0z_eta_2_lam_mul,V2,G0,S1p,npol)
! - comp
  call collocate2_sum_tensor_1d(m0z_eta_4lm_w_l3ml,V2, &
                                m0z_eta_2_lam_mul,V1,G0,S1m,npol)
! z comp
  call collocate2_sum_tensor_1d(m0_w_mul,V3,m0s_xi_mul,V4,G0,S1z,npol)

! Final VxM D_z^z in + component
  call collocate0_1d(m0s_xi_mul,V3,V4,npol) 
  call vxm(V4,G2T,V1)
  call add_to_axis_2d(V1,S1p,npol)

end subroutine additional_di_ax
!=============================================================================


!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!    Q U A D R U P O L E   E L E M E N T A L   S T I F F N E S S   T E R M S 
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

!-----------------------------------------------------------------------------
subroutine stiffness_quad_non_ax(us,uphi,uz, &
           mz_xi_laml,ms_xi_laml,mz_eta_laml,ms_eta_laml, &
           mz_xi_2mul,ms_xi_2mul,mz_eta_2mul,ms_eta_2mul, &
           m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,m_w_4mul, &
           m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
           m11zl,m21zl,m41zl,m1phil,m2phil,m4phil,&
           loc_stiffness_s,loc_stiffness_phi,loc_stiffness_z)

include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (s, phi and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz,uphi

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrix for W_s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_lam_6mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_min_2lam_6mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4lam_9mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4mul

! precomputed matrices for W_s^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_xi_laml,mz_eta_laml
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi_laml,ms_eta_laml
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_xi_2mul,mz_eta_2mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi_2mul,ms_eta_2mul

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

  loc_stiffness_s=zero
  loc_stiffness_phi=zero
  loc_stiffness_z=zero
  S1s=zero;S2s=zero;S1phi=zero;S2phi=zero;S1z=zero;S2z=zero
  X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! First MxM
  call mxm(G2T,us,X1)
  call mxm(G2T,uphi,X2)
  call mxm(G2T,uz,X3)
  call mxm(us,G2,X4)
  call mxm(uphi,G2,X5)
  call mxm(uz,G2,X6)

! Collocations and sums of W terms: s and phi components
  call collocate10s_sum_1d(mz_xi_laml,X4,mz_eta_laml,X1,mz_xi_2mul,X5, &
                           mz_eta_2mul,X2,ms_xi_laml,X6,ms_eta_laml,X3, &
                           m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,us,uphi, &
                           loc_stiffness_s,loc_stiffness_phi,nsize)

! Collocations and sums of W terms: z component
  call collocate3_sum_1d(ms_xi_2mul,X5,ms_eta_2mul,X2,m_w_4mul,uz, &
                         loc_stiffness_z,nsize)

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! Collocations and sums of D terms

! s component
  call collocate6s_sum_1d(m11sl,X4,m21sl,X1,m12sl,X6,m22sl,X3, &
                          mz_eta_laml,us,uphi,S1s,nsize)
  call collocate6s_sum_1d(m11sl,X1,m41sl,X4,m32sl,X3,m42sl,X6, &
                          mz_xi_laml,us,uphi,S2s,nsize)
! z component
  call collocate6s_sum_1d(m11zl,X6,m21zl,X3,m32sl,X4,m22sl,X1, &
                          ms_eta_laml,us,uphi,S1z,nsize)
  call collocate6s_sum_1d(m11zl,X3,m41zl,X6,m12sl,X1,m42sl,X4, &
                          ms_xi_laml,us,uphi,S2z,nsize)

! phi component
  call collocate5ss_sum_1d(m1phil,X5,m2phil,X2,mz_eta_2mul,us, &
                           mz_eta_2mul,uphi,ms_eta_2mul,uz,S1phi,nsize)
  call collocate5ss_sum_1d(m1phil,X2,m4phil,X5,mz_xi_2mul,us, &
                           mz_xi_2mul,uphi,ms_xi_2mul,uz,S2phi,nsize)

!Second MxM
  call mxm(G2,S1s,X1)
  call mxm(S2s,G2T,X2)
  call mxm(G2,S1phi,X3)
  call mxm(S2phi,G2T,X4)
  call mxm(G2,S1z,X5)
  call mxm(S2z,G2T,X6)

! Final Sum
  call sum3s_3_1d(X1,X2,loc_stiffness_s,X3,X4,loc_stiffness_phi,X5,X6, &
                 loc_stiffness_z,nsize)

end subroutine stiffness_quad_non_ax
!=============================================================================

!-----------------------------------------------------------------------------
subroutine stiffness_quad_ax(us,uphi,uz, &
           mz_xi_laml,ms_xi_laml,mz_eta_laml,ms_eta_laml, &
           mz_xi_2mul,ms_xi_2mul,mz_eta_2mul,ms_eta_2mul, &
           m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,m_w_4mul, &
           m11sl,m21sl,m41sl,m12sl,m22sl,m32sl,m42sl, &
           m11zl,m21zl,m41zl,m1phil,m2phil,m4phil, &
           m0_w_4mul,m0_z_eta_2l_w_l_6ml, &
           m0_z_eta_3m_w_4l_9ml,m0_z_eta_2lm_min_w_2l_6ml,&
           loc_stiffness_s,loc_stiffness_phi,loc_stiffness_z)

include "mesh_params.h"

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
! displacement (s, phi and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: us,uz,uphi

! precomputed matrices containing weights, mapping, elastic parameters
! precomputed matrix for W_s
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_lam_6mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_min_2lam_6mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4lam_9mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m_w_4mul

! precomputed matrices for W_s^d, D_x^y
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_xi_laml,mz_eta_laml
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi_laml,ms_eta_laml
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: mz_xi_2mul,mz_eta_2mul
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: ms_xi_2mul,ms_eta_2mul

! precomputed matrices for D_xy^xy
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11sl,m21sl,m41sl
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m12sl,m22sl,m32sl,m42sl
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m11zl,m21zl,m41zl
real(kind=realkind), dimension(0:npol,0:npol), intent(in) :: m1phil,m2phil,m4phil

! axial terms
real(kind=realkind), dimension(0:npol), intent(in) :: m0_w_4mul
real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_2l_w_l_6ml
real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_3m_w_4l_9ml
real(kind=realkind), dimension(0:npol), intent(in) :: m0_z_eta_2lm_min_w_2l_6ml

! result: local (elemental) stiffness term (s, phi and z components)
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_s
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_phi
real(kind=realkind), dimension(0:npol,0:npol), intent(out) :: loc_stiffness_z

! work arrays
real(kind=realkind), dimension(0:npol,0:npol) :: X1,X2,X3,X4,X5,X6 ! MxM arrays
real(kind=realkind), dimension(0:npol,0:npol) :: S1s,S2s,S1phi,S2phi,S1z,S2z ! Sum
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  loc_stiffness_s=zero
  loc_stiffness_phi=zero
  loc_stiffness_z=zero
  S1s=zero;S2s=zero;S1phi=zero;S2phi=zero;S1z=zero;S2z=zero
  X1=zero;X2=zero;X3=zero;X4=zero;X5=zero;X6=zero

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! First MxM
  call mxm(G1T,us,X1)
  call mxm(G1T,uphi,X2)
  call mxm(G1T,uz,X3)
  call mxm(us,G2,X4)
  call mxm(uphi,G2,X5)
  call mxm(uz,G2,X6)

! Collocations and sums of W terms: s and phi components
  call collocate10s_sum_1d(mz_xi_laml,X4,mz_eta_laml,X1,mz_xi_2mul,X5, &
                           mz_eta_2mul,X2,ms_xi_laml,X6,ms_eta_laml,X3, &
                           m_w_lam_6mul,m_w_min_2lam_6mul,m_w_4lam_9mul,us,uphi, &
                           loc_stiffness_s,loc_stiffness_phi,nsize)

! Collocations and sums of W terms: z component
  call collocate3_sum_1d(ms_xi_2mul,X5,ms_eta_2mul,X2,m_w_4mul,uz, &
                         loc_stiffness_z,nsize)

!++++++++++++++++++COMMON TO ALL SOURCES+++++++++++++++++++++++++++++++++++
! Collocations and sums of D terms

! s component
  call collocate6s_sum_1d(m11sl,X4,m21sl,X1,m12sl,X6,m22sl,X3, &
                          mz_eta_laml,us,uphi,S1s,nsize)
  call collocate6s_sum_1d(m11sl,X1,m41sl,X4,m32sl,X3,m42sl,X6, &
                          mz_xi_laml,us,uphi,S2s,nsize)
! z component
  call collocate6s_sum_1d(m11zl,X6,m21zl,X3,m32sl,X4,m22sl,X1, &
                          ms_eta_laml,us,uphi,S1z,nsize)
  call collocate6s_sum_1d(m11zl,X3,m41zl,X6,m12sl,X1,m42sl,X4, &
                          ms_xi_laml,us,uphi,S2z,nsize)

! phi component
  call collocate5ss_sum_1d(m1phil,X5,m2phil,X2,mz_eta_2mul,us, &
                           mz_eta_2mul,uphi,ms_eta_2mul,uz,S1phi,nsize)
  call collocate5ss_sum_1d(m1phil,X2,m4phil,X5,mz_xi_2mul,us, &
                           mz_xi_2mul,uphi,ms_xi_2mul,uz,S2phi,nsize)

!Second MxM
  call mxm(G1,S1s,X1)
  call mxm(S2s,G2T,X2)
  call mxm(G1,S1phi,X3)
  call mxm(S2phi,G2T,X4)
  call mxm(G1,S1z,X5)
  call mxm(S2z,G2T,X6)

  call additional_quad_ax(m0_z_eta_2l_w_l_6ml,m0_z_eta_3m_w_4l_9ml,&
                          m0_w_4mul,m0_z_eta_2lm_min_w_2l_6ml, &
                          us,uphi,uz,S1s,S1phi,S1z)

! Final Sum
  call sum4_3_1d(X1,X2,S1s,loc_stiffness_s,X3,X4,S1phi,loc_stiffness_phi,X5, &
                 X6,S1z,loc_stiffness_z,nsize)

end subroutine stiffness_quad_ax
!=============================================================================



!-----------------------------------------------------------------------------
subroutine additional_quad_ax(m0_z_eta_2l_w_l_6ml,m0_z_eta_3m_w_4l_9ml,&
                              m0_w_4mul,m0_z_eta_2lm_min_w_2l_6ml, &
                              us,uphi,uz,S1s,S1phi,S1z)

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

  V1=zero;V2=zero;V3=zero
  S1s=zero;S1phi=zero;S1z=zero

! VxM
  call vxm(G0,us,V1)
  call vxm(G0,uphi,V2)
  call vxm(G0,uz,V3)

! Collocations, Sums, Tensorization

! s comp
  call collocate2_sum_tensor_1d(m0_z_eta_2l_w_l_6ml,V1,&
                                m0_z_eta_2lm_min_w_2l_6ml,V2,G0,S1s,npol)
                                
! phi comp
  call collocate2_sum_tensor_1d(m0_z_eta_2lm_min_w_2l_6ml,V1, &
                                m0_z_eta_3m_w_4l_9ml,V2,G0,S1phi,npol)
                                
! z comp
  call collocate_tensor_1d(m0_w_4mul,V3,G0,S1z,npol)

end subroutine additional_quad_ax
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
