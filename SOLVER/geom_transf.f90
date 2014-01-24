!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
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

!=====================
module geom_transf
!=====================

    use data_mesh
    use analytic_mapping

    implicit none
    public :: jacobian

    public :: alpha,  beta,  gamma1, delta,  epsilon1, zeta
    public :: alphak, betak, gammak, deltak, epsilonk, zetak

    public :: jacobian_srf,  quadfunc_map, grad_quadfunc_map
    public :: mgrad_pointwise, mgrad_pointwisek
    public :: mapping, s_over_oneplusxi_axis

    private 
contains


!-----------------------------------------------------------------------------------------
pure real(kind=dp) function mapping(xil, etal, nodes_crd, iaxis, ielem0)

    integer, intent(in)       :: iaxis, ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    mapping = mapping_anal(xil,etal,nodes_crd,iaxis,ielem0)

    if ( iaxis == 1 .and. dabs(mapping/router) < 1.d-12 ) mapping = 0.d0 

end function mapping
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function quadfunc_map(p, s, z, nodes_crd, ielem0)
!<  This routines computes the 
!!  quadratic functional (s-s(xi,eta))**2 + (z-z(xi,eta))**2

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: p(2), s, z, nodes_crd(8,2)

    quadfunc_map = quadfunc_map_anal(p,s,z,nodes_crd,ielem0)

end function quadfunc_map
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine grad_quadfunc_map(grd, p, s, z, nodes_crd, ielem0)
!< This routine returns the gradient of the quadratic
!! functional associated with the mapping.
!
    integer, intent(in)        :: ielem0
    real(kind=dp), intent(in)  :: p(2), s, z, nodes_crd(8,2)
    real(kind=dp), intent(out) :: grd(2)

    call grad_quadfunc_map_anal(grd,p,s,z,nodes_crd,ielem0)

end subroutine grad_quadfunc_map
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function s_over_oneplusxi_axis(xil, etal, nodes_crd, ielem0)
!< This routine returns the value of the quantity
!!  
!!              s/(1+xi) 
!!
!! when the associated element lies along the axis of 
!! symmetry, in the case of an analytical transformation. 

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    s_over_oneplusxi_axis=s_over_oneplusxi_axis_anal(xil,etal,nodes_crd,ielem0)

end function s_over_oneplusxi_axis
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
pure real(kind=dp) function jacobian(xil, etal, nodes_crd, ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    jacobian = jacobian_anal(xil,etal,nodes_crd,ielem0)   

end function jacobian
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function jacobian_srf(xil,crdedge,ielem0)
!<  This routine computes the Jacobian of the transformation
!!  that maps [-1,+1] into a portion of the boundary of domain.  
    use subpar_mapping, only: jacobian_srf_subpar

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, crdedge(3,2)

    if (eltype(ielem0) /= 'linear') &
       jacobian_srf = jacobian_srf_anal(xil,crdedge)  
    if (eltype(ielem0) == 'linear') &
       jacobian_srf = jacobian_srf_subpar(xil,crdedge)  


end function jacobian_srf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function alphak(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    alphak = alphak_anal(xil,etal,nodes_crd,ielem0)

end function alphak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function betak(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    betak = betak_anal(xil,etal,nodes_crd,ielem0)

end function betak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function gammak(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    gammak = gammak_anal(xil,etal,nodes_crd,ielem0)

end function gammak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function deltak(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    deltak = deltak_anal(xil,etal,nodes_crd,ielem0)

end function deltak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function epsilonk(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    epsilonk = epsilonk_anal(xil,etal,nodes_crd,ielem0)

end function epsilonk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function zetak(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    zetak = zetak_anal(xil,etal,nodes_crd,ielem0)

end function zetak
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function alpha(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    alpha = alpha_anal(xil,etal,nodes_crd,ielem0)

end function alpha
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function beta(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    beta = beta_anal(xil,etal,nodes_crd,ielem0)

end function beta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function gamma1(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    gamma1 = gamma_anal(xil,etal,nodes_crd,ielem0)

end function gamma1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function delta(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    delta = delta_anal(xil,etal,nodes_crd,ielem0)

end function delta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function epsilon1(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    epsilon1 = epsilon_anal(xil,etal,nodes_crd,ielem0)

end function epsilon1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function zeta(xil,etal,nodes_crd,ielem0)

    integer, intent(in)       :: ielem0
    real(kind=dp), intent(in) :: xil, etal, nodes_crd(8,2)

    zeta = zeta_anal(xil,etal,nodes_crd,ielem0)

end function zeta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine mgrad_pointwise(mg,xil,etal,nodes_crd,ielem0)
!
! This routines returns the following matrix:
!                      +                     +
!                      |(ds/dxi)  | (ds/deta)|
!    mg =  s(xi,eta) * | ---------|--------- |(xi,eta)
!                      |(dz/dxi ) | (dz/deta)|
!                      +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.

    integer, intent(in)        :: ielem0
    real(kind=dp), intent(in)  :: xil, etal, nodes_crd(8,2)
    real(kind=dp), intent(out) :: mg(2,2)

    call mgrad_pointwise_anal(mg,xil,etal,nodes_crd,ielem0)

end subroutine mgrad_pointwise
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure  subroutine mgrad_pointwisek(mg,xil,etal,nodes_crd,ielem0)

! This routines returns the following matrix:
!          +                     +
!          |(ds/dxi)  | (ds/deta)|
!    mg =  | ---------|--------- |(xi,eta)
!          |(dz/dxi ) | (dz/deta)|
!          +                     +
!       This 2*2 matrix is needed when defining and storing
!gradient/divergence related arrays.

    integer, intent(in)        :: ielem0
    real(kind=dp), intent(in)  :: xil, etal, nodes_crd(8,2)
    real(kind=dp), intent(out) :: mg(2,2)

    call mgrad_pointwisek_anal(mg,xil,etal,nodes_crd,ielem0)

end subroutine mgrad_pointwisek
!----------------------------------------------------------------------------------------


!=========================
end module geom_transf
!=========================
