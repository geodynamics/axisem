module data_gllmesh
  implicit none
  public 
  double precision, dimension(:,:,:), allocatable :: sgll,zgll
  double precision, dimension(:,:,:), allocatable :: sgll_fluid,zgll_fluid
  double precision, dimension(:,:,:), allocatable :: sgll_solid,zgll_solid

  double precision :: hmin_glob, hmax_glob ! global min/max gll spacing
  double precision :: min_distance_dim ! 0.1*hmin_glob [in meters]
  double precision :: min_distance_nondim ! 0.1*hmin_glob [referenced to 1]
  
  double precision :: char_time_max
  integer :: char_time_max_globel
  double precision :: char_time_max_rad,char_time_max_theta
  double precision :: char_time_min
  integer :: char_time_min_globel
  double precision :: char_time_min_rad,char_time_min_theta

end module data_gllmesh
