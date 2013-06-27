!> Global numbering for the stiffness assembly
!! (size: all unique grid points in respective domains)
!===================
module data_numbering
!===================


implicit none
public
include "mesh_params.h"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! global number in solid varies across procs due to central cube domain decomposition
 integer                          :: nglob,nglob_solid 
 integer, dimension(npoint_solid) :: igloc_solid
 integer, dimension(npoint_fluid) :: igloc_fluid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!===================
end module data_numbering
!===================
