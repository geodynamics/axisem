module data_grid

implicit none
public

integer :: ns, nz
double precision  :: ri,ro
double precision, dimension(:,:,:),allocatable :: crd_grdc,crd_grds
double precision, dimension(:), allocatable :: s_unif,z_unif
double precision, dimension(:), allocatable :: radius
double precision, dimension(:), allocatable :: ndeta

!number of subdivisions for central square
integer :: ndivs 
double precision :: lsq

!TNM
double precision :: router,lsq_fac
integer :: ngllcube
integer, parameter :: ndim=2

double precision, dimension(:), allocatable :: rmax_el, rmin_el

!ALEX
logical :: southern
end module data_grid
