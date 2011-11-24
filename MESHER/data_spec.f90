!===================
  module data_spec
!===================
!
  implicit none
!
  public 
!
!
  integer :: npol
  double precision, dimension(:),allocatable   :: xi,xi_k,eta
  double precision, dimension(:),allocatable   :: dxi
  double precision, dimension (:), allocatable :: wt, wt_axial!Quadrature weights
  double precision, dimension (:), allocatable :: wt_axial_k  !Quad. wgts for the 
                                                              !nonaxisymmetric components
!
!=======================
  end module data_spec
!=======================
