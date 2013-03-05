!================
module utlity
!================

  use global_parameters
  use data_mesh, only : smallval
  implicit none
  
  public :: compute_coordinates, scoord, zcoord, rcoord, thetacoord
  public :: dblreldiff_small, reldiff_small
  public :: dblereldiff, reldiff
  public :: dbleabsreldiff, absreldiff
  public :: dbleps2zero, eps2zero
  private

contains

!-----------------------------------------------------------------------------
logical function dblreldiff_small(x1,x2)

  double precision, intent(in) :: x1,x2

  dblreldiff_small = .false.

  if (x1 /= zero) then 
     if (abs((x1-x2)/x1) <= smallval) dblreldiff_small = .true.
  elseif (x2 /=zero) then
     if (abs((x1-x2)/x2) <= smallval) dblreldiff_small = .true.
  else
     dblreldiff_small = .true.
  endif

end function dblreldiff_small
!=============================================================================

!-----------------------------------------------------------------------------
logical function reldiff_small(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2
  real(kind=realkind)             ::  smallval1

  if (realkind==4) smallval1 = smallval
  if (realkind==8) smallval1 = smallval_dble

  reldiff_small = .false.

  if (x1 /= zero) then 
     if (abs((x1-x2)/x1) <= smallval1) reldiff_small = .true.
  elseif (x2 /=zero) then
     if (abs((x1-x2)/x2) <= smallval1) reldiff_small = .true.
  else
     reldiff_small = .true.
  endif

end function reldiff_small
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function reldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1/=zero) then
     reldiff=(x1-x2)/x1
  elseif (x2/=zero) then
     reldiff=(x1-x2)/x2
  else
     reldiff=zero
  endif

end function reldiff
!=============================================================================

!-----------------------------------------------------------------------------
double precision function dblereldiff(x1,x2)

  double precision, intent(in) :: x1,x2

  if (x1/=zero) then
     dblereldiff=(x1-x2)/x1
  elseif (x2/=zero) then 
     dblereldiff=(x1-x2)/x2
  else
     dblereldiff=zero
  endif

end function dblereldiff
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function absreldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1/=zero) then
     absreldiff=abs((x1-x2)/x1)
  elseif (x2/=zero) then 
     absreldiff=abs((x1-x2)/x2)
  else
     absreldiff=zero
  endif

end function absreldiff
!=============================================================================

!-----------------------------------------------------------------------------
double precision function dbleabsreldiff(x1,x2)

  double precision, intent(in) :: x1,x2

  if (x1/=zero) then
     dbleabsreldiff=abs((x1-x2)/x1)
  elseif (x2/=zero) then 
     dbleabsreldiff=abs((x1-x2)/x2)
  else
     dbleabsreldiff=zero
  endif

end function dbleabsreldiff
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_coordinates(s,z,r,theta,ielem,ipol,jpol)
  !
  ! Given the elemental grid point index, outputs s,z,r,theta coordinate [m,rad].
  ! These coordinates are by default ALWAYS global (no solid or fluid domains).
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_mesh,            ONLY: min_distance_dim
  use data_mesh_preloop,    ONLY: lnods, crd_nodes, axis
  use data_spec,            ONLY: xi, xi_k, eta
  use geom_transf,          ONLY: mapping
  
  double precision :: s,z,r,theta
  integer          :: ielem,ipol,jpol,ipt,inode
  double precision :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 

     s= mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)

  else 
     s= mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(s) < min_distance_dim) s=zero
  if (abs(z) < min_distance_dim) z=zero

  r = dsqrt(s**2+z**2)
  theta = datan(s/(z+epsi))
  if ( zero > theta ) theta = pi + theta
  if (theta == zero .and. z < 0) theta = pi

end subroutine compute_coordinates
!=============================================================================

!-----------------------------------------------------------------------------
double precision function scoord(ipol,jpol,ielem)
  !
  ! Given the elemental grid point index, outputs the s coordinate [m].
  ! These coordinates are by default ALWAYS global (no solid or fluid domains).
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_mesh,            ONLY: min_distance_dim
  use data_mesh_preloop,    ONLY: lnods, crd_nodes, axis
  use data_spec,            ONLY : xi, xi_k, eta
  use geom_transf,          ONLY: mapping
  
  integer          :: ielem,ipol,jpol,ipt,inode
  double precision :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     scoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
  else 
     scoord = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(scoord) < min_distance_dim) scoord=zero

end function scoord
!=============================================================================

!-----------------------------------------------------------------------------
double precision function zcoord(ipol,jpol,ielem)
  !
  ! Given the elemental grid point index, outputs the z coordinate [m].
  ! These coordinates are by default ALWAYS global (no solid or fluid domains).
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_mesh,            ONLY: min_distance_dim
  use data_mesh_preloop,    ONLY: lnods, crd_nodes, axis
  use data_spec,            ONLY : xi, xi_k, eta
  use geom_transf,          ONLY: mapping
  
  integer          :: ielem,ipol,jpol,ipt,inode
  double precision :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     zcoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else 
     zcoord = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(zcoord) < min_distance_dim) zcoord=zero

end function zcoord
!=============================================================================

!-----------------------------------------------------------------------------
double precision function rcoord(ipol,jpol,ielem)
  !
  ! Given the elemental grid point index, outputs the radius coordinate [m].
  ! These coordinates are by default ALWAYS global (no solid or fluid domains).
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_mesh,            ONLY: min_distance_dim
  use data_mesh_preloop,    ONLY: lnods, crd_nodes, axis
  use data_spec,            ONLY : xi, xi_k, eta
  use geom_transf,          ONLY: mapping
  
  integer          :: ielem,ipol,jpol,ipt,inode
  double precision :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     s = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else 
     s = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  rcoord = sqrt(s**2 + z**2)
  ! Eliminate roundoff errors
  if (abs(rcoord) < min_distance_dim) rcoord=zero

end function rcoord
!=============================================================================

!-----------------------------------------------------------------------------
double precision function thetacoord(ipol,jpol,ielem)
  !
  ! Given the elemental grid point index, outputs the theta coordinate [rad].
  ! These coordinates are by default ALWAYS global (no solid or fluid domains).
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_mesh, ONLY: min_distance_dim
  use data_mesh_preloop, ONLY: lnods,crd_nodes,axis
  use data_spec, ONLY : xi,xi_k,eta
  use geom_transf, ONLY: mapping
  
  integer          :: ielem,ipol,jpol,ipt,inode
  double precision :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     s = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else 
     s = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  thetacoord = datan(s/(z+epsi))
  if ( zero > thetacoord ) thetacoord = pi + thetacoord
  if (thetacoord == zero .and. z < 0) thetacoord = pi

end function thetacoord
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function eps2zero(val)

  real(kind=realkind) :: val

  eps2zero = val
  if (abs(val) < epsi) eps2zero = zero

end function eps2zero
!=============================================================================

!-----------------------------------------------------------------------------
double precision function dbleps2zero(val)

  double precision :: val

  dbleps2zero=val
  if (abs(val) < epsi) dbleps2zero = zero

end function dbleps2zero
!=============================================================================

!====================
end module utlity
!====================
