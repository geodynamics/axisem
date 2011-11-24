!===================
module data_mesh_preloop
!===================
!
! This collective has been extracted from data_mesh to separate those 
! large mesh arrays that only need to be known *prior* to the time loop.
! Remaining inevitable mesh info arrays such as axis_solid are in 
! data_mesh and in active memory throughout the simulation.

use global_parameters
implicit none
include "mesh_params.h"
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Global mesh informations
  integer, dimension(:,:), allocatable          :: lnods ! (nelem,1:8)
  character(len=6), dimension(:), allocatable   :: eltype ! (nelem)
  integer                                       :: npoin
  double precision, dimension(:,:), allocatable :: crd_nodes ! (npoin,2)
  logical, dimension(:), allocatable            :: coarsing,north ! (nelem)
  integer                                       :: num_spher_radii
  double precision, dimension(:), allocatable   :: spher_radii
! Axial elements
  logical, dimension(:), allocatable            :: axis ! (nelem)

! Mapping between solid/fluid elements:
! integer array of size nel_fluid containing the glocal (global per-proc)
! element number for 1:nel_solid/fluid
  integer, dimension(:), allocatable            :: ielsolid ! (nel_solid)
  integer, dimension(:), allocatable            :: ielfluid ! (nel_fluid)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_mesh_preloop
!=======================
