!> This determines the precision for the memory-/CPU-intensive time loop. 
!! Set the parameter realkind to either 
!!  sp: single precision (half memory compared to 8, faster on many systems)
!!  dp: double precision (more expensive (double memory), but more precise.
!! The mesher is intrinsically double precision, as are all precomputed, mesh 
!! related variables. This distinction is only relevant for the global 
!! arrays used in the time evolution.
!=========================
 module global_parameters
!=========================
!

implicit  none
public

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  integer, parameter         :: sp = kind(0.0)
  integer, parameter         :: dp = kind(0.0d0)
  integer, parameter         :: realkind = sp  !< Choose solver precision here

! Do not change these unless problems with any of the accuracy tests arise.
! As floating point rounding is system-dependent, there might be different 
! numbers for different systems, but the below values seem generally reasonable. 
  real(kind=sp), parameter    :: smallval_sngl = 1e-6
  real(kind=dp), parameter    :: smallval_dble = 1e-11

! Do not change these.
  double precision, parameter :: zero = 0d0, half = 5d-1, third = 1d0 / 3d0
  double precision, parameter :: quart = 25d-2, one = 1d0, sixth = 1d0 / 6d0
  double precision, parameter :: two = 2d0, three = 3d0, four = 4d0, five = 5.d0
  double precision, parameter :: fifth = 2.d-1
  double precision, parameter :: pi = 3.1415926535898D0
  double precision, parameter :: epsi = 1d-30

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!============================= 
 end module global_parameters
!=============================
