module gllmeshgen

  use data_mesh
  use analytic_mapping 
  use data_spec
  use data_gllmesh 
  use splib

  implicit none

  public :: test_mapping,create_gllmesh
  public :: extract_fluid_solid_submeshes
  private

  contains

!----------------------------------------------------------------------------
subroutine create_gllmesh

  use data_diag

  double precision :: crd_nodes(8,2)
  integer npoin
!w
  integer iel,jpol,ipol
  double precision stest
!
! Generate collocation points in the two directions of space
  npoin = neltot*(npol+1)**2
  allocate(sgll(0:npol,0:npol,neltot),zgll(0:npol,0:npol,neltot))

!-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! QUADRATURE POINTS and weights
  allocate (eta(0:npol))
  allocate (xi(0:npol),dxi(0:npol))
  allocate (wt(0:npol),wt_axial(0:npol))
  allocate(xi_k(0:npol), wt_axial_k(0:npol))

  call zemngr(npol,xi)                      ! Pseudo Gauss-Radau quadrature
  call get_welegl_axial(npol,xi,wt_axial,3) !
!
  call zemngl2(npol,xi_k)                       ! Gauss-Jacobi(0,1) quadrature
  call get_welegl_axial(npol,xi_k,wt_axial_k,2) !
!
! In the z-direction and in the s-direction for any other element
!
  call ZELEGL(npol,eta,dxi)                 ! Gauss-Lobatto Points
  call get_welegl(npol,eta,wt)              !


!!TNMOUCHE
!do ipol=0,npol
!write(6,*)'TNMOUCHE:',eta(ipol),xi_k(ipol)
!enddo
!!/TNMOUCHE


!-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do iel = 1, neltot
!
! define dummy coordinate arrays
  crd_nodes(:,:) = 0.d0
  crd_nodes(1,1) = sg(lnodesg(1,iel)) ; crd_nodes(1,2) = zg(lnodesg(1,iel))
  crd_nodes(3,1) = sg(lnodesg(2,iel)) ; crd_nodes(3,2) = zg(lnodesg(2,iel))
  crd_nodes(5,1) = sg(lnodesg(3,iel)) ; crd_nodes(5,2) = zg(lnodesg(3,iel))
  crd_nodes(7,1) = sg(lnodesg(4,iel)) ; crd_nodes(7,2) = zg(lnodesg(4,iel))
  crd_nodes(2,:) = .5d0 * ( crd_nodes(1,:) + crd_nodes(3,:) )  ! midpoints are necessary
  crd_nodes(4,:) = .5d0 * ( crd_nodes(3,:) + crd_nodes(5,:) )  ! for subparametric mapping
  crd_nodes(6,:) = .5d0 * ( crd_nodes(5,:) + crd_nodes(7,:) )  ! (Serendipity elements).
  crd_nodes(8,:) = .5d0 * ( crd_nodes(7,:) + crd_nodes(1,:) )  !  
!
!
  stest = minval(sg(lnodesg(1:4,iel)))
  if ( stest < smallval ) then 
!  if (dump_mesh_info_screen) write(6,*) 'Element ',iel,' is an axial element '
   do jpol = 0, npol
    do ipol = 0, npol
     sgll(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,1,iel)
     zgll(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,2,iel)
    end do
   end do  
  else
   do jpol = 0, npol
    do ipol = 0, npol
     sgll(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,1,iel)
     zgll(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,2,iel)
    end do
   end do
  end if

!    if ( stest < smallval ) then 
!    do jpol = 0, npol-1
!         do ipol = 0, npol-1
!            write(111,*)sgll(ipol,jpol,iel),zgll(ipol,jpol,iel)
!         end do
!      end do
!   endif

 end do

end subroutine create_gllmesh
!--------------------------------------------------------------------------
!dk extract_fluid_solid_submeshes------------------------------------------
  subroutine extract_fluid_solid_submeshes
  integer :: iel_fluid,iel
  integer :: iel_solid
!
!
  allocate(sgll_fluid(0:npol,0:npol,neltot_fluid)) ; sgll_fluid(:,:,:) = 0.d0 
  allocate(zgll_fluid(0:npol,0:npol,neltot_fluid)) ; zgll_fluid(:,:,:) = 0.d0 
  do iel_fluid = 1, neltot_fluid
   iel = ielem_fluid(iel_fluid)
   sgll_fluid(:,:,iel_fluid) = sgll(:,:,iel)
   zgll_fluid(:,:,iel_fluid) = zgll(:,:,iel)
  end do
! 
  allocate(sgll_solid(0:npol,0:npol,neltot_solid)) ; sgll_solid(:,:,:) = 0.d0 
  allocate(zgll_solid(0:npol,0:npol,neltot_solid)) ; zgll_solid(:,:,:) = 0.d0 
  do iel_solid = 1, neltot_solid
   iel = ielem_solid(iel_solid)
   sgll_solid(:,:,iel_solid) = sgll(:,:,iel)
   zgll_solid(:,:,iel_solid) = zgll(:,:,iel)
  end do
! 
  end subroutine extract_fluid_solid_submeshes
!--------------------------------------------------------------------------

!dk test_mapping
subroutine test_mapping
  double precision :: crd_nodes(8,2)
  integer npoin
!w
  integer iel,jpol,ipol
  double precision stest

  double precision, dimension(:,:,:), allocatable :: sglltmp,zglltmp
!
write(6,*)'Test mapping...'; call flush(6)

! Generate collocation points in the  two directions of space
  npoin = neltot*(npol+1)**2
  allocate(sglltmp(0:npol,0:npol,neltot),zglltmp(0:npol,0:npol,neltot))

!-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! QUADRATURE POINTS 
!  allocate (eta(0:npol))
!  allocate (xi(0:npol),dxi(0:npol))
!  allocate (wt(0:npol),wt_axial(0:npol))
!  allocate(xi_k(0:npol), wt_axial_k(0:npol))
!
  call zemngr(npol,xi)                      ! Pseudo Gauss-Radau quadrature
  call get_welegl_axial(npol,xi,wt_axial,3) !
!
  call zemngl2(npol,xi_k)                       ! Gauss-Jacobi(0,1) quadrature
  call get_welegl_axial(npol,xi_k,wt_axial_k,2) !
!
! In the z-direction and in the s-direction for any other element
!
  call ZELEGL(npol,eta,dxi)                 ! Gauss-Lobatto Points
  call get_welegl(npol,eta,wt)              !
!
!-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do iel = 1, neltot
!
! define dummy coordinate arrays
  crd_nodes(:,:) = 0. 
  crd_nodes(1,1) = sg(lnodesg(1,iel)) ; crd_nodes(1,2) = zg(lnodesg(1,iel))
  crd_nodes(3,1) = sg(lnodesg(2,iel)) ; crd_nodes(3,2) = zg(lnodesg(2,iel))
  crd_nodes(5,1) = sg(lnodesg(3,iel)) ; crd_nodes(5,2) = zg(lnodesg(3,iel))
  crd_nodes(7,1) = sg(lnodesg(4,iel)) ; crd_nodes(7,2) = zg(lnodesg(4,iel))
  crd_nodes(2,:) = .5d0 * ( crd_nodes(1,:) + crd_nodes(3,:) )  ! midpoints are necessary
  crd_nodes(4,:) = .5d0 * ( crd_nodes(3,:) + crd_nodes(5,:) )  ! for subparametric mapping
  crd_nodes(6,:) = .5d0 * ( crd_nodes(5,:) + crd_nodes(7,:) )  ! (Serendipity elements).
  crd_nodes(8,:) = .5d0 * ( crd_nodes(7,:) + crd_nodes(1,:) )  !  
!
!
  stest = minval(sg(lnodesg(1:4,iel)))
  if ( stest < smallval ) then 
!  write(6,*) 'Element ',iel,' is an axial element '
   do jpol = 0, npol
    do ipol = 0, npol
     sglltmp(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,1,iel)
     zglltmp(ipol,jpol,iel) = mapping_anal(xi_k(ipol),eta(jpol),crd_nodes,2,iel)
    end do
   end do  
  else
   do jpol = 0, npol
    do ipol = 0, npol
     sglltmp(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,1,iel)
     zglltmp(ipol,jpol,iel) = mapping_anal(eta(ipol),eta(jpol),crd_nodes,2,iel)
    end do
   end do
  end if
 end do
!
   open(21,file='colloc_grid.dat')
   open(23,file='mesh.dat')
   do iel = 1,neltot
! do iel = neltot/2 - 16 , neltot/2
    do jpol = 0, npol
     do ipol = 0, npol
        write(21,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
     end do
     write(21,*)
    end do
                                                  
    do ipol = 0, npol
     do jpol = 0, npol
      write(21,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
     end do
     write(21,*)
    end do

    do jpol = 0,npol,npol
     do ipol = 0, npol
      write(23,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
     end do
     write(23,*)
    end do

    do ipol = 0,npol,npol
     do jpol = 0, npol
      write(23,*) sglltmp(ipol,jpol,iel),zglltmp(ipol,jpol,iel)
     end do
     write(23,*)
    end do
   end do
  close(23) ; close(21)
!
  deallocate(sglltmp,zglltmp,eta,xi,dxi,wt,wt_axial,xi_k,wt_axial_k)

end subroutine test_mapping

!----------------------------------------------------------------------
!dk midpoint
subroutine midpoint(s,z,s1,z1,s2,z2)
  double precision, intent(out) :: s,z
  double precision, intent(in ) :: s1,z1,s2,z2
  s = .5d0*(s1+s2)
  z = .5d0*(z1+z2)
end subroutine midpoint

end module gllmeshgen
