!===================
module data_bkgrdmodel
!===================

  implicit none
  
  integer                       :: ndisc, nfluidregions
  integer, allocatable          :: idom_fluid(:)
  double precision, allocatable :: discont(:)
  double precision, allocatable :: vp(:,:), vs(:,:), rho(:,:)
  logical, allocatable          :: solid_domain(:)
  integer                       :: lfbkgrdmodel
  character(len=100)            :: bkgrdmodel
  logical                       :: resolve_inner_shear, have_fluid, have_solid
  double precision              :: pts_wavelngth
  double precision              :: period, courant
  double precision              :: dt
  integer                       :: nc_init, nproc_target
  
  ! the sole quantities to be created in create_subregions
  ! that are needed by the rest of the mesher
  integer                       :: nz_glob, ns_glob, nc_glob
  integer, allocatable          :: iclev_glob(:)
  double precision, allocatable :: dz_glob(:)
  double precision              :: rmin, minh_ic, maxh_ic, maxh_icb
  double precision              :: minhvp, maxhvs, maxhnsicb

!===================
end module data_bkgrdmodel
!===================
