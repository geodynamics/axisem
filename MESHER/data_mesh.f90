module data_mesh

use global_parameters

implicit none
public
!SKELETON PARAMETERS
integer neltot
integer npointot
double precision, dimension(:), allocatable :: sg,zg
integer, dimension(:,:), allocatable :: lnodesg
character(len=6), dimension(:), allocatable :: eltypeg
logical, dimension(:), allocatable :: coarsing
!Region to which the element pertains, in the case of a stratified bkgrd model
integer, dimension(:), allocatable :: region,nel_region
double precision, dimension(:), allocatable :: scom,zcom,thetacom,rcom ! com = center of mass

! axial elements
integer, allocatable :: naxelp(:),naxel_solidp(:),naxel_fluidp(:)
integer, allocatable :: ax_elp(:,:),ax_el_solidp(:,:),ax_el_fluidp(:,:)
integer,allocatable :: axis(:,:),axis_solid(:,:),axis_fluid(:,:)
logical, allocatable :: have_axis(:)

! Solid fluid distinction
logical, dimension(:), allocatable :: fluid
integer :: neltot_fluid
integer, dimension(:), allocatable :: ielem_fluid,inv_ielem_fluid
!
logical, dimension(:), allocatable :: solid
integer :: neltot_solid
integer, dimension(:), allocatable :: ielem_solid,inv_ielem_solid

!Boundary informations (ICB, CMB, from solid and fluid perspectives)
integer                            :: nbcnd    ! Number of boundaries
integer, dimension(:),allocatable  :: nbelem 
integer, dimension(:,:), allocatable :: belem  ! List of boundary elements
integer, dimension(:,:),allocatable :: my_neighbour

! Solid-fluid boundary arrays (needed by solver)
integer, dimension(:,:),allocatable  :: bdry_above_el, bdry_below_el 
integer, dimension(:,:),allocatable  :: bdry_solid_el, bdry_fluid_el 

double precision, dimension(:,:,:), allocatable :: bdry_s,bdry_z
integer, dimension(:,:,:),allocatable  :: bdry_globnum_above
integer, dimension(:,:,:),allocatable  :: bdry_globnum_below
integer, dimension(:,:,:),allocatable  :: bdry_locnum_above
integer, dimension(:,:,:),allocatable  :: bdry_locnum_below

integer, dimension(:,:,:),allocatable  :: bdry_globnum_solid
integer, dimension(:,:,:),allocatable  :: bdry_globnum_fluid
integer, dimension(:,:,:),allocatable  :: bdry_locnum_solid
integer, dimension(:,:,:),allocatable  :: bdry_locnum_fluid

integer, dimension(:,:), allocatable :: bdry_jpol_solid,bdry_jpol_fluid

!real(kind=realkind) :: smallval
double precision :: smallval

! central region
  double precision, allocatable :: s_arr(:,:),z_arr(:,:)
  integer, allocatable :: central_is_iz_to_globiel(:,:)
contains
subroutine empty_data_mesh
  deallocate(rcom,scom,zcom,thetacom)
end subroutine empty_data_mesh

end module data_mesh
