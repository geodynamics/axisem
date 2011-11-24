!==========================
  module analytic_mapping
!==========================
!  
  use global_parameters
  use analytic_spheroid_mapping
  use analytic_semi_mapping
  use subpar_mapping
  use data_mesh

  implicit none

! public :: jacobian_anal,alpha_anal,beta_anal
! public :: gamma_anal,delta_anal,epsilon_anal,zeta_anal
! public :: alphak_anal,betak_anal,gammak_anal
! public :: deltak_anal,epsilonk_anal,zetak_anal
! public :: jacobian_srf_anal, quadfunc_map_anal,grad_quadfunc_map_anal
! public :: mgrad_pointwise_anal,mgrad_pointwisek_anal
  public :: mapping_anal!,s_over_oneplusxi_axis_anal
 
  private

  contains
!/////////////////////////////////////////////////////////////
!
!dk mapping_anal---------------------------------
  double precision function mapping_anal(xi,eta,nodes_crd,iaxis,ielem0)
!
!
  integer        :: iaxis,ielem0
  double precision :: xi, eta, nodes_crd(8,2)
 !  write(61,*)'IEL MAPPING ANAL:',ielem0,eltypeg(ielem0)
  if (eltypeg(ielem0) == 'curved') mapping_anal = map_spheroid(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'linear') mapping_anal = mapping_subpar(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'semino') mapping_anal = map_semino(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'semiso') mapping_anal = map_semiso(xi,eta,nodes_crd,iaxis)

  end function mapping_anal
!-----------------------------------------------
!
!dk quadfunc_map_anal--------------------------------------------
  double precision function quadfunc_map_anal(p,s,z,nodes_crd,ielem0)
!
!        This routines computes the 
!quadratic functional (s-s(xi,eta))**2 + (z-z(xi,eta))**2
!
  integer :: ielem0
  double precision :: p(2), xi,eta,s,z, nodes_crd(8,2)
!
!
  xi  = p(1)
  eta = p(2)
!
  quadfunc_map_anal = (s-mapping_anal(xi,eta,nodes_crd,1,ielem0))**2 &
               +  (z-mapping_anal(xi,eta,nodes_crd,2,ielem0))**2

  end function quadfunc_map_anal
!-----------------------------------------------------------------
!
!dk grad_quadfunc_map_anal------------------------------------------
  subroutine grad_quadfunc_map_anal(grd,p,s,z,nodes_crd,ielem0)
!
!       This routine returns the gradient of the quadratic
!functional associated with the mapping.
!
  integer :: ielem0
  double precision :: grd(2),p(2),xi, eta, s,z, nodes_crd(8,2)
  double precision :: dsdxi,dzdxi,dsdeta,dzdeta
!
  xi  = p(1)
  eta = p(2)

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)

  grd(1) = -((s-mapping_anal(xi,eta,nodes_crd,1,ielem0))*dsdxi&
            +(z-mapping_anal(xi,eta,nodes_crd,2,ielem0))*dzdxi)
  grd(2) = -((s-mapping_anal(xi,eta,nodes_crd,1,ielem0))*dsdeta&
            +(z-mapping_anal(xi,eta,nodes_crd,2,ielem0))*dzdeta)

  end subroutine grad_quadfunc_map_anal
!--------------------------------------------------------------
!
!dk s_over_oneplusxi_axis_anal--------------------------------------------
  double precision function s_over_oneplusxi_axis_anal(xi,eta,nodes_crd,ielem0)
! 
! This routine returns the value of the quantity
!  
!              s/(1+xi) 
!
! when the associated element lies along the axis of 
! symmetry, in the case of an analytical transformation. 
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdxi,dsdeta,dzdeta

  if ( xi == -one ) then 
!    Apply L'Hopital's rule
     call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
!    call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
     s_over_oneplusxi_axis_anal = dsdxi
  else
     s_over_oneplusxi_axis_anal = mapping_anal(xi,eta,nodes_crd,1,ielem0)/(one+xi)
  end if

  end function s_over_oneplusxi_axis_anal
!--------------------------------------------------------------------
!
!dk jacobian_anal-----------------------------------------------------------
  double precision function jacobian_anal(xi,eta,nodes_crd,ielem0)
!
!	This function returns the value of the jacobian of the
! analytical mapping between the reference square [-1,1]^2 and
! the deformed element in the spheroid. 

  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdxi,dsdeta,dzdeta
!

  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)

  jacobian_anal = dsdxi*dzdeta-dsdeta*dzdxi

  end function jacobian_anal
!
!---------------------------------------------------------------------------
!
!dk jacobian_srf_anal-----------------------------------------------
  double precision function jacobian_srf_anal(xi,crdedge)
!
!       This routine computes the Jacobian of the transformation
!that maps [-1,+1] into a portion of the boundary of domain.  
!
!         xi(or eta)
!        ---->
! 1 - - - 2 - - - 3 .
!
  double precision :: xi, crdedge(3,2)
  double precision :: dsdxi,dzdxi,s1,s2,s3,z1,z2,z3
  double precision :: thetabar,deltatheta,a,b
  double precision :: arg,dist
!
  s1 = crdedge(1,1) ; s2 = crdedge(2,1) ; s3 = crdedge(3,1)
  z1 = crdedge(1,2) ; z2 = crdedge(2,2) ; z3 = crdedge(3,2)
!
! write(6,*) ' parameters ',  s1,s3,z1,z3
  call compute_parameters_srf(s1,s3,z1,z3,a,b,deltatheta,thetabar)
!
  if (dabs(deltatheta) > 1.0d-10 ) then 
     arg = xi*deltatheta + thetabar
!
     dsdxi = -deltatheta*a*dsin(arg)
     dzdxi =  deltatheta*b*dcos(arg) 

     jacobian_srf_anal = dsqrt(dsdxi**2+dzdxi**2)
  
!    write(6,*) ielem0,dsqrt(dsdxi**2+dzdxi**2),dsdxi,dzdxi,s1,s3,z1,z3,a,b,deltatheta,thetabar
!    write(6,*) ielem0,jacobian_srf_anal,deltatheta
  else 
     dist = dsqrt((z3-z1)**2+(s3-s1)**2)
     jacobian_srf_anal = half*dist
!    write(6,*) ielem0,jacobian_srf_anal, dist,deltatheta
  end if

!
  end function jacobian_srf_anal
!---------------------------------------------------------------------------
!
!dk alphak_anal------------------------------------------------------
  double precision function alphak_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    alphak =  ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  alphak_anal  = -inv_jacob*dsdxi*dsdeta
!
  end function alphak_anal
!---------------------------------------------------------------------------
!
!dk betak_anal------------------------------------------------------
  double precision function betak_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    betak =  ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.

  implicit none
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  betak_anal  = inv_jacob*dsdxi**2
!
  end function betak_anal
!---------------------------------------------------------------------------
!
!dk gammak_anal------------------------------------------------------
  double precision function gammak_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    gammak =  ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  gammak_anal  = inv_jacob*dsdeta**2
!
  end function gammak_anal
!---------------------------------------------------------------------------
!
!dk deltak_anal------------------------------------------------------
  double precision function deltak_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    deltak = -( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  deltak_anal  = -inv_jacob*dzdxi*dzdeta
!
  end function deltak_anal
!--------------------------------------------------------------------
!
!dk epsilonk_anal------------------------------------------------------
  double precision function epsilonk_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    epsilonk = ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  epsilonk_anal  = inv_jacob*dzdxi**2
!
  end function epsilonk_anal
!---------------------------------------------------------------------------
!
!dk zetak_anal------------------------------------------------------
  double precision function zetak_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    zetak_anal = ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  zetak_anal  = inv_jacob*dzdeta**2
!
  end function zetak_anal
!---------------------------------------------------------------------------
!
!dk alpha_anal------------------------------------------------------
  double precision function alpha_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    alpha = s(xi,eta) * ( -ds/dxi ) * ( ds/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  alpha_anal  = -inv_jacob*dsdxi*dsdeta*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function alpha_anal
!---------------------------------------------------------------------------
!
!dk beta_anal------------------------------------------------------
  double precision function beta_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    beta =  s(xi,eta) * ( ds/dxi ) * ( ds/dxi) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.

  implicit none
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  beta_anal  = inv_jacob*dsdxi**2*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function beta_anal
!---------------------------------------------------------------------------
!
!dk gamma_anal------------------------------------------------------
  double precision function gamma_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    gamma = s(xi,eta) * ( ds/deta ) * ( ds/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  gamma_anal  = inv_jacob*dsdeta**2*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function gamma_anal
!---------------------------------------------------------------------------
!
!dk delta_anal------------------------------------------------------
  double precision function delta_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    delta = -s(xi,eta) * ( dz/dxi ) * ( dz/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  delta_anal  = -inv_jacob*dzdxi*dzdeta*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function delta_anal
!--------------------------------------------------------------------
!
!dk epsilon_anal------------------------------------------------------
  double precision function epsilon_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    epsilon = s(xi,eta) * ( dz/dxi ) * ( dz/dxi) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  epsilon_anal  = inv_jacob*dzdxi**2*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function epsilon_anal
!---------------------------------------------------------------------------
!
!dk zeta_anal------------------------------------------------------
  double precision function zeta_anal(xi,eta,nodes_crd,ielem0)
!
! This routines returns the value of 
!
!    zeta_anal = s(xi,eta) * ( dz/deta ) * ( dz/deta) / J(xi,eta),
!
!a quantity that is needed in the calculation of the laplacian
!operator. alpha is defined within an element, and s(xi,eta) is 
!defined by the analytic transformation .J is the determinant of 
!the Jacobian matrix of the transformation.
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta,inv_jacob
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  inv_jacob  = one/(dsdxi*dzdeta - dsdeta*dzdxi)
!
  zeta_anal  = inv_jacob*dzdeta**2*mapping_anal(xi,eta,nodes_crd,1,ielem0)
!
  end function zeta_anal
!---------------------------------------------------------------------------
!
!dk mgrad_pointwise_anal------------------------------------------------------
  subroutine mgrad_pointwise_anal(mg,xi,eta,nodes_crd,ielem0)
!
! This routines returns the following matrix:
!                      +                     +
!                      |(ds/dxi)  | (ds/deta)|
!    mg =  s(xi,eta) * | ---------|--------- |(xi,eta)
!                      |(dz/dxi ) | (dz/deta)|
!                      +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.
!
  implicit none
!
  integer :: ielem0
  double precision :: mg(2,2) 
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dsdeta,dzdxi,dzdeta

  mg(:,:) = zero
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  mg(1,1)  = mapping_anal(xi,eta,nodes_crd,1,ielem0)*dsdxi
  mg(1,2)  = mapping_anal(xi,eta,nodes_crd,1,ielem0)*dsdeta
  mg(2,1)  = mapping_anal(xi,eta,nodes_crd,1,ielem0)*dzdxi
  mg(2,2)  = mapping_anal(xi,eta,nodes_crd,1,ielem0)*dzdeta
!
  end subroutine mgrad_pointwise_anal
!---------------------------------------------------------------------------
!
!dk mgrad_pointwisek_anal------------------------------------------------------
  subroutine mgrad_pointwisek_anal(mg,xi,eta,nodes_crd,ielem0)
!
! This routines returns the following matrix:
!           +                     +
!           |(ds/dxi)  | (ds/deta)|
!    mg =   | ---------|--------- |(xi,eta)
!           |(dz/dxi ) | (dz/deta)|
!           +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.
!
  implicit none
!
  integer :: ielem0
  double precision :: mg(2,2) 
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dsdeta,dzdxi,dzdeta

  mg(:,:) = zero
!
  call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
! call compute_partial_derivatives_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  mg(1,1)  = dsdxi
  mg(1,2)  = dsdeta
  mg(2,1)  = dzdxi
  mg(2,2)  = dzdeta
!
  end subroutine mgrad_pointwisek_anal
!---------------------------------------------------------------------------
!
!dk compute_partial_derivatives------------------------------------------------------
  subroutine compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd,ielem0)
!
!	This routine returns the analytical values of the partial derivatives
! of the analytic spheroidal mapping. 
!
  integer :: ielem0
  double precision :: xi, eta, nodes_crd(8,2)
  double precision :: dsdxi,dzdxi,dsdeta,dzdeta

  if (eltypeg(ielem0) == 'curved') & 
  call comp_partial_deriv_spheroid(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
  if (eltypeg(ielem0) == 'linear') & 
  call comp_partial_deriv_subpar(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
  if (eltypeg(ielem0) == 'semino') & 
  call comp_partial_deriv_semino(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
  if (eltypeg(ielem0) == 'semiso') & 
  call comp_partial_deriv_semiso(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)

  end subroutine compute_partial_derivatives
!------------------------------------------------------------------------------------
!
!dk compute_parameters---------------------------------------------------------------
  subroutine compute_parameters(nodes_crd,a1,a2,b1,b2,deltatheta,thetabar)
!
  double precision :: nodes_crd(8,2)
  double precision :: a1,a2,b1,b2,deltatheta,thetabar,theta3,theta1
  double precision ::  s1,z1,s3,z3,s5,z5,s7,z7

  s1 = nodes_crd(1,1) ; z1 = nodes_crd(1,2)
  s3 = nodes_crd(3,1) ; z3 = nodes_crd(3,2)
  s5 = nodes_crd(5,1) ; z5 = nodes_crd(5,2)
  s7 = nodes_crd(7,1) ; z7 = nodes_crd(7,2)
! 
  theta1 = datan(s1/(z1+epsi)) ; theta3 = datan(s3/(z3+epsi))
!
  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0.d0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0.d0) theta3 = pi
!
  a1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2))
  a2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(z5**2-z7**2))
! 
  b1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2))
  b2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(s7**2-s5**2))

! b1 = a1  ! Quick and dirty fix for the time being
! b2 = a2

! write(6,*) a1,a2,b1,b2

  deltatheta = half*(theta3-theta1)
  thetabar   = half*(theta3+theta1)   

! write(6,*) deltatheta,thetabar

  end subroutine compute_parameters
!--------------------------------------------------------------------------
!
!dk compute_parameters_new------------------------------------------------------------
  subroutine compute_parameters_new(nodes_crd,a1,a2,b1,b2,deltatheta1,thetabar1,deltatheta2,thetabar2)
!
  double precision :: nodes_crd(8,2)
  double precision :: a1,a2,b1,b2,deltatheta1,thetabar1,deltatheta2,thetabar2
  double precision :: theta3,theta1,theta5,theta7
  double precision ::  s1,z1,s3,z3,s5,z5,s7,z7

  s1 = nodes_crd(1,1) ; z1 = nodes_crd(1,2)
  s3 = nodes_crd(3,1) ; z3 = nodes_crd(3,2)
  s5 = nodes_crd(5,1) ; z5 = nodes_crd(5,2)
  s7 = nodes_crd(7,1) ; z7 = nodes_crd(7,2)
! 
  theta1 = datan(s1/(z1+epsi)) ; theta3 = datan(s3/(z3+epsi))
  theta7 = datan(s7/(z7+epsi)) ; theta5 = datan(s5/(z5+epsi))
!
  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0.d0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0.d0) theta3 = pi
  if ( zero > theta5 ) theta5 = pi + theta5
  if (theta5 == zero .and. z5 < 0.d0) theta5 = pi
  if ( zero > theta7 ) theta7 = pi + theta7
  if (theta7 == zero .and. z7 < 0.d0) theta7 = pi
!
  a1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2))
  a2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(z5**2-z7**2))
! 
  b1 = dsqrt(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2))
  b2 = dsqrt(((s7*z5)**2 - (s5*z7)**2)/(s7**2-s5**2))

! b1 = a1  ! Quick and dirty fix for the time being
! b2 = a2

! write(6,*) a1,a2,b1,b2

  deltatheta1 = half*(theta3-theta1)
  thetabar1   = half*(theta3+theta1)
  deltatheta2 = half*(theta5-theta7)
  thetabar2   = half*(theta5+theta7)

! write(6,*) deltatheta,thetabar

  end subroutine compute_parameters_new
!--------------------------------------------------------------------------
!
!dk compute_parameters_srf-----------------------------------------------------------
  subroutine compute_parameters_srf(s1,s3,z1,z3,a,b,deltatheta,thetabar)
!
  double precision,intent(out) :: a,b,deltatheta,thetabar
  double precision :: theta3,theta1
  double precision,intent(in) ::  s1,z1,s3,z3
! 
  a= zero ; b = zero ; deltatheta = zero ; thetabar = zero

  if (z1/=z3) a = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(z3**2-z1**2)))
!  
  if (s1/=s3) b = dsqrt(dabs(((s1*z3)**2 - (s3*z1)**2)/(s1**2-s3**2)))
!
  theta1 = datan(s1*b/(z1*a+epsi)) ; theta3 = datan(s3*b/(z3*a+epsi))
!
  if ( zero > theta1 ) theta1 = pi + theta1
  if (theta1 == zero .and. z1 < 0.d0) theta1 = pi
  if ( zero > theta3 ) theta3 = pi + theta3
  if (theta3 == zero .and. z3 < 0.d0) theta3 = pi
!
  deltatheta = half*(theta3-theta1)
  thetabar   = half*(theta3+theta1)   

  end subroutine compute_parameters_srf
!--------------------------------------------------------------------------
!
!/////////////////////////////////////////////////////////////
!
!=============================
  end module analytic_mapping
!=============================
