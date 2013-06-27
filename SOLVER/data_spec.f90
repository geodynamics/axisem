!> Variables concerned with elemental & spectral aspects only
!! (e.g. GLL points, Lagrange interpolant derivatives, quadrature weights)
!===================
 module data_spec
!===================

use global_parameters
implicit none
include 'mesh_params.h'
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  double precision, dimension(0:npol) :: xi_k, eta
  double precision, dimension(0:npol) :: dxi
  double precision, dimension(0:npol) :: wt          !Quadrature weights
  double precision, dimension(0:npol) :: wt_axial_k  !Quad. wgts for the   
                                                     !gaus jacobi(0,1) integration

! Lagrange interpolant derivatives
! 3rd index: 1 - non-axial, 2 - axial
! 4th index: 1 - \partial_\xi, 2 - \partial_\eta
  real(kind=realkind), dimension(0:npol,0:npol,2,2) :: shp_deri_k
  real(kind=realkind), dimension(0:npol,0:npol)     :: G1, G1T, G2, G2T
  real(kind=realkind), dimension(0:npol)            :: G0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_spec
!=======================
