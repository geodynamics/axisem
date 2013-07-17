!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stephanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!==================================
  module analytic_semi_mapping
!==================================
!
!	10/01/2002: This module contains the 
! machinery necessary to describe analytically
! the transformation of the reference element
! into its deformed image in a semi-curved
! enveloppe. 
! 
  use global_parameters 
  use data_mesh, only : smallval

  implicit none
  public :: map_semino,comp_partial_deriv_semino
  public :: map_semiso,comp_partial_deriv_semiso
  private
  contains
!//////////////////////////////////////////

!dk map_semino------------------------------------------------------
  double precision function map_semino(xi,eta,crd_nodes,idir)
!
!	We are working in polar coordinates here: theta
! is the latitude. 
!
  double precision :: xi, eta
  double precision, dimension(8,2),intent(in) :: crd_nodes
  integer :: idir

  double precision :: atop,btop
  double precision :: thetabartop,dthetatop 
  double precision :: sbot,zbot,stop,ztop
  double precision :: sbar,ds,slope,intersect
!
  map_semino = zero
  call comp_parameters_semino(crd_nodes,atop,btop,thetabartop,dthetatop)

  call comp_sz_xi_line_no(sbot,zbot,xi,crd_nodes)
  call comp_sz_xi(stop,ztop,xi,atop,btop,thetabartop,dthetatop)

  sbar = half*(sbot+stop); ds = stop-sbot

  if (dabs(ds)>smallval) then
     intersect = (zbot*stop-ztop*sbot)/ds   
     slope = (ztop-zbot)/ds
  end if   
  
  if (idir == 1) then
     map_semino = sbar+ds*eta*half
  elseif (idir == 2) then
     if (dabs(ds)>smallval) then
     map_semino = slope*(sbar+half*ds*eta)+intersect 
     else
     map_semino = half*(zbot+ztop)+eta*(ztop-zbot)*half
     end if
  end if

  end function map_semino
!----------------------------------------------------------------------
!
!dk map_semiso------------------------------------------------------
  double precision function map_semiso(xi,eta,crd_nodes,idir)
!
!       We are working in polar coordinates here: theta is the latitude. 
!
  double precision :: xi, eta
  double precision, dimension(8,2),intent(in) :: crd_nodes
  integer :: idir
  
  double precision :: abot,bbot
  double precision :: thetabarbot,dthetabot
  double precision :: sbot,zbot,stop,ztop
  double precision :: sbar,ds,slope,intersect
! 
  map_semiso = zero
  call comp_parameters_semiso(crd_nodes,abot,bbot,thetabarbot,dthetabot)
  
  call comp_sz_xi_line_so(stop,ztop,xi,crd_nodes)
  call comp_sz_xi(sbot,zbot,xi,abot,bbot,thetabarbot,dthetabot)
  
  sbar = half*(sbot+stop); ds = stop-sbot
  
  if (dabs(ds)>smallval) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = (ztop-zbot)/ds
  end if
  
  if (idir == 1) then
     map_semiso = sbar+ds*eta*half
  elseif (idir == 2) then
     if (dabs(ds)>smallval) then
     map_semiso = slope*(sbar+half*ds*eta)+intersect
     else
     map_semiso = half*(zbot+ztop)+eta*(ztop-zbot)*half
     end if
  end if

  end function map_semiso
!----------------------------------------------------------------------
!
!dk comp_partial_deriv_semino-----------------------
  subroutine comp_partial_deriv_semino(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)

  double precision, intent(out) :: dsdxi,dzdxi,dsdeta,dzdeta
  double precision, intent(in) :: xi,eta
  double precision, dimension(8,2),intent(in) :: nodes_crd

  double precision :: atop,btop
  double precision :: thetabartop,dthetatop
  double precision :: sbot,zbot,stop,ztop
  double precision :: sbar,ds,dz,slope,intersect,sxieta
  double precision :: dsbotdxi,dzbotdxi
  double precision :: dstopdxi,dztopdxi
  double precision :: dsbardxi,ddsdxi
  double precision :: dzbardxi,ddzdxi
  double precision :: dslopedxi,dintersectdxi
!
! call comp_parameters_sph(nodes_crd,abot,bbot,atop,btop,&
!                             thetabarbot,dthetabot,thetabartop,dthetatop)

! call comp_sz_xi(sbot,zbot,xi,abot,bbot,thetabarbot,dthetabot)
! call comp_sz_xi(stop,ztop,xi,atop,btop,thetabartop,dthetatop)
  call comp_parameters_semino(nodes_crd,atop,btop,thetabartop,dthetatop)

  call comp_sz_xi_line_no(sbot,zbot,xi,nodes_crd)
  call comp_sz_xi(stop,ztop,xi,atop,btop,thetabartop,dthetatop)

  sbar = half*(sbot+stop); ds = stop-sbot ; dz = ztop - zbot
  if (dabs(ds)>smallval) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = dz/ds
  end if

! call comp_dsdxi_dzdxi(dsbotdxi,dzbotdxi,xi,abot,bbot,thetabarbot,dthetabot)
  call comp_dsdxi_dzdxi_line_no(dsbotdxi,dzbotdxi,nodes_crd)
  call comp_dsdxi_dzdxi(dstopdxi,dztopdxi,xi,atop,btop,thetabartop,dthetatop)
!
  dsbardxi = half*(dsbotdxi+dstopdxi)
  ddsdxi = (dstopdxi-dsbotdxi)
 
  dzbardxi = half*(dzbotdxi+dztopdxi)
  ddzdxi = (dztopdxi-dzbotdxi) 
!
  sxieta = sbar+ds*eta*half
!
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   dsdxi = dsbardxi + half*eta*ddsdxi
  dsdeta = half*ds
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
  if (dabs(ds)>smallval) then
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          dslopedxi = (ddzdxi*ds-ddsdxi*dz)/ds**2
      dintersectdxi = ((dzbotdxi*stop-dztopdxi*sbot     &
                       +zbot*dstopdxi-ztop*dsbotdxi)*ds &
                      -ddsdxi*(zbot*stop-ztop*sbot))/ds**2
     dzdxi = (dz/ds)*dsdxi+ sxieta*dslopedxi + dintersectdxi
    dzdeta = (dz/ds)*dsdeta  
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  else
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dzdxi = zero
    dzdeta = half*dz
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  end if

  end subroutine comp_partial_deriv_semino
!--------------------------------------------------------------
!
!dk comp_partial_derivat_semiso-----------------------
  subroutine comp_partial_deriv_semiso(dsdxi,dzdxi,dsdeta,dzdeta,xi,eta,nodes_crd)
!
  double precision, intent(out) :: dsdxi,dzdxi,dsdeta,dzdeta
  double precision, intent(in) :: xi,eta
  double precision, dimension(8,2),intent(in) :: nodes_crd

  double precision :: abot,bbot
  double precision :: thetabarbot,dthetabot
  double precision :: sbot,zbot,stop,ztop
  double precision :: sbar,ds,dz,slope,intersect,sxieta
  double precision :: dsbotdxi,dzbotdxi
  double precision :: dstopdxi,dztopdxi
  double precision :: dsbardxi,ddsdxi
  double precision :: dzbardxi,ddzdxi
  double precision :: dslopedxi,dintersectdxi
!
! call comp_parameters_sph(nodes_crd,abot,bbot,atop,btop,&
!                             thetabarbot,dthetabot,thetabartop,dthetatop)

! call comp_sz_xi(sbot,zbot,xi,abot,bbot,thetabarbot,dthetabot)
! call comp_sz_xi(stop,ztop,xi,atop,btop,thetabartop,dthetatop)
  call comp_parameters_semiso(nodes_crd,abot,bbot,thetabarbot,dthetabot)
  call comp_sz_xi_line_so(stop,ztop,xi,nodes_crd)
  call comp_sz_xi(sbot,zbot,xi,abot,bbot,thetabarbot,dthetabot)

  sbar = half*(sbot+stop); ds = stop-sbot ; dz = ztop - zbot
  if (dabs(ds)>smallval) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = dz/ds
  end if

  call comp_dsdxi_dzdxi(dsbotdxi,dzbotdxi,xi,abot,bbot,thetabarbot,dthetabot)
  call comp_dsdxi_dzdxi_line_so(dstopdxi,dztopdxi,nodes_crd)
! call comp_dsdxi_dzdxi(dstopdxi,dztopdxi,xi,atop,btop,thetabartop,dthetatop)
!
  dsbardxi = half*(dsbotdxi+dstopdxi)
  ddsdxi = (dstopdxi-dsbotdxi)
 
  dzbardxi = half*(dzbotdxi+dztopdxi)
  ddzdxi = (dztopdxi-dzbotdxi) 
!
  sxieta = sbar+ds*eta*half
!
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   dsdxi = dsbardxi + half*eta*ddsdxi
  dsdeta = half*ds
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
  if (dabs(ds)>smallval) then
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          dslopedxi = (ddzdxi*ds-ddsdxi*dz)/ds**2
      dintersectdxi = ((dzbotdxi*stop-dztopdxi*sbot     &
                       +zbot*dstopdxi-ztop*dsbotdxi)*ds &
                      -ddsdxi*(zbot*stop-ztop*sbot))/ds**2
     dzdxi = (dz/ds)*dsdxi+ sxieta*dslopedxi + dintersectdxi
    dzdeta = (dz/ds)*dsdeta  
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  else
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     dzdxi = zero
    dzdeta = half*dz
!--------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  end if

  end subroutine comp_partial_deriv_semiso
!--------------------------------------------------------------
!
!dk comp_sz_xi------------------------------------------------------
  subroutine comp_sz_xi(s,z,xi,a,b,thetabar,dtheta)

  double precision, intent(out) :: s,z
  double precision, intent(in) :: xi, a, b, thetabar,dtheta
  
  s = a*dcos(thetabar+xi*half*dtheta)
  z = b*dsin(thetabar+xi*half*dtheta)

  end subroutine comp_sz_xi
!----------------------------------------------------------------------
!
!dk comp_sz_xi_line_no-------------------------------------
  subroutine comp_sz_xi_line_no(s,z,xi,crd_nodes)

  double precision, intent(out) :: s,z
  double precision, intent(in) :: xi,crd_nodes(8,2)
  double precision :: s1,z1,s3,z3
  
  s1 = crd_nodes(1,1) ; z1 = crd_nodes(1,2)
  s3 = crd_nodes(3,1) ; z3 = crd_nodes(3,2)

  s = half*((one+xi)*s3+(one-xi)*s1)
  z = half*((one+xi)*z3+(one-xi)*z1)

  end subroutine comp_sz_xi_line_no
!----------------------------------------------------------------------
!
!dk comp_sz_xi_line_so-------------------------------------
  subroutine comp_sz_xi_line_so(s,z,xi,crd_nodes)

  double precision, intent(out) :: s,z
  double precision, intent(in) :: xi,crd_nodes(8,2)
  double precision :: s7,z7,s5,z5

  s7 = crd_nodes(7,1) ; z7 = crd_nodes(7,2)
  s5 = crd_nodes(5,1) ; z5 = crd_nodes(5,2)

  s = half*((one+xi)*s5+(one-xi)*s7)
  z = half*((one+xi)*z5+(one-xi)*z7)

  end subroutine comp_sz_xi_line_so
!----------------------------------------------------------------------

!
!dk comp_dsdxi_dzdxi---------------------------
  subroutine comp_dsdxi_dzdxi(dsdxi,dzdxi,xi,a,b,thetabar,dtheta)
!
!
  double precision, intent(out) :: dsdxi,dzdxi
  double precision, intent(in) :: xi,a,b,thetabar,dtheta 

  dsdxi =-a*half*dtheta*dsin(thetabar+xi*half*dtheta)
  dzdxi = b*half*dtheta*dcos(thetabar+xi*half*dtheta)  
 
  end subroutine comp_dsdxi_dzdxi
!-------------------------------------------------
!
!dk comp_dsdxi_dzdxi_line_no------------------------------------------
  subroutine comp_dsdxi_dzdxi_line_no(dsdxi,dzdxi,crd_nodes)
!
!
  double precision, intent(out) :: dsdxi,dzdxi
  double precision, intent(in) :: crd_nodes(8,2)
  double precision :: s1,z1,s3,z3
  
  s1 = crd_nodes(1,1) ; z1 = crd_nodes(1,2)
  s3 = crd_nodes(3,1) ; z3 = crd_nodes(3,2)
  
  dsdxi = half*(s3-s1) ; dzdxi = half*(z3-z1)

  end subroutine comp_dsdxi_dzdxi_line_no
!--------------------------------------------------------------------------------
!
!dk comp_dsdxi_dzdxi_line_so------------------------------------------
  subroutine comp_dsdxi_dzdxi_line_so(dsdxi,dzdxi,crd_nodes)
!
!
  double precision, intent(out) :: dsdxi,dzdxi
  double precision, intent(in) :: crd_nodes(8,2)
  double precision :: s7,z7,s5,z5
  
  s7 = crd_nodes(7,1) ; z7 = crd_nodes(7,2)
  s5 = crd_nodes(5,1) ; z5 = crd_nodes(5,2)
  
  dsdxi = half*(s5-s7) ; dzdxi = half*(z5-z7)

  end subroutine comp_dsdxi_dzdxi_line_so
!--------------------------------------------------------------------------------
!
!dk comp_parameters_semino---------------------------------
  subroutine comp_parameters_semino(crd_nodes,atop,btop,thetabartop,dthetatop)

  double precision, dimension(8,2),intent(in) :: crd_nodes
  double precision,intent(out) :: atop,btop
  double precision,intent(out) :: thetabartop,dthetatop
  double precision :: s1,z1,s3,z3,s5,z5,s7,z7
  double precision :: theta5,theta7
!  
  s1 = crd_nodes(1,1) ; z1 = crd_nodes(1,2)
  s3 = crd_nodes(3,1) ; z3 = crd_nodes(3,2)
  s5 = crd_nodes(5,1) ; z5 = crd_nodes(5,2)
  s7 = crd_nodes(7,1) ; z7 = crd_nodes(7,2)
!
  call compute_ab(atop,btop,s7,z7,s5,z5)
  call compute_theta(theta5,s5,z5,atop,btop) 
  call compute_theta(theta7,s7,z7,atop,btop) 
!
  thetabartop = half*(theta7+theta5)
  dthetatop = theta5 -theta7
!
  end subroutine comp_parameters_semino
!----------------------------------------------------------
!
!dk comp_parameters_semiso---------------------------------
  subroutine comp_parameters_semiso(crd_nodes,abot,bbot,thetabarbot,dthetabot)

  double precision, dimension(8,2),intent(in) :: crd_nodes
  double precision,intent(out) :: abot,bbot
  double precision,intent(out) :: thetabarbot,dthetabot
  double precision :: s1,z1,s3,z3,s5,z5,s7,z7
  double precision :: theta1,theta3
!  
  s1 = crd_nodes(1,1) ; z1 = crd_nodes(1,2)
  s3 = crd_nodes(3,1) ; z3 = crd_nodes(3,2)
  s5 = crd_nodes(5,1) ; z5 = crd_nodes(5,2)
  s7 = crd_nodes(7,1) ; z7 = crd_nodes(7,2)
!
  call compute_ab(abot,bbot,s1,z1,s3,z3)
  call compute_theta(theta3,s3,z3,abot,bbot)
  call compute_theta(theta1,s1,z1,abot,bbot)
!
  thetabarbot = half*(theta1+theta3)
  dthetabot = theta3 -theta1
!
  end subroutine comp_parameters_semiso
!----------------------------------------------------------
!
!dk compute_ab--------------------------------------
  subroutine compute_ab(a,b,s1,z1,s2,z2)

  double precision, intent(out) :: a,b
  double precision, intent(in) :: s1,z1,s2,z2
!
  a = dsqrt(dabs((s2**2*z1**2-z2**2*s1**2)/(z1**2-z2**2)))
  b = dsqrt(dabs((z1**2*s2**2-z2**2*s1**2)/(s2**2-s1**2))) 
  end subroutine compute_ab
!---------------------------------------------------
!
!dk compute_theta---------------------------
  subroutine compute_theta(theta,s,z,a,b)
!
!	This routine returns the latitude
! theta, given s and z.
!
  double precision, intent(out) :: theta
  double precision, intent(in) :: s,z,a,b
!
  if (s /= zero) then
     theta=datan(z*a/(s*b))  
  else
     if (z>0) theta=half*pi
     if (z<0) theta=-half*pi
  end if

  end subroutine compute_theta
!-----------------------------------
!//////////////////////////////////////////
!
!=======================================
  end module analytic_semi_mapping
!=======================================
