!==============================================================================
module meshgen
!==============================================================================
!   
implicit none

public :: generate_skeleton
private

double precision, dimension(:), allocatable :: aspect_ratio
!REGION BY REGION PARAMETERS
integer, dimension(:,:), allocatable :: lnodeso   ! OUTER SHELL
character(len=6), dimension(:), allocatable :: eltypeo
logical, dimension(:), allocatable :: coarsingo

integer :: nelo
double precision, dimension(:), allocatable :: so,zo
!
integer, dimension(:,:), allocatable :: lnodesi   ! INNER SHELL
character(len=6), dimension(:), allocatable :: eltypei
logical, dimension(:), allocatable :: coarsingi
integer :: neli
double precision, dimension(:), allocatable :: si,zi
!
integer, dimension(:,:), allocatable :: lnodessq  ! INNER SQUARE
character(len=6), dimension(:), allocatable :: eltypesq
integer :: nelsq
double precision, dimension(:), allocatable :: ssq,zsq
!
integer, dimension(:,:), allocatable :: lnodesbuf ! BUFFER SHELL
character(len=6), dimension(:), allocatable :: eltypebuf
integer :: nelbuf
double precision, dimension(:), allocatable :: sbuf,zbuf
!
contains

!----------------------------------------------------------------------------
!dk generate_skeleton
subroutine generate_skeleton

  use data_grid
  use data_diag
  use data_mesh
  use data_numbering
  use numbering

! OUTER SHELL GENERATION
!  define reference grid
   if (dump_mesh_info_screen) then
    write(6,*)'generating reference spherical grid....';call flush(6)
   end if
   call def_reference_spherical_grid_discont

! Aspect ratio
   write(6,*)'estimate aspect ratio...';call flush(6)
   call estimate_aspect_ratio

!  final outer shell w/ potential coarsening w/ depth
   write(6,*)'define spherical shell...';call flush(6)
   call define_spherical_shell
   deallocate(crd_grds,crd_grdc)
   deallocate(z_unif,s_unif)

! CENTRAL REGION GENERATION (inner shell + central square + buffer layer)
   write(6,*)'define central region...';call flush(6)
   call define_central_region

! GATHER INFORMATION FROM DIFFERENT REGIONS IN GLOBAL ARRAYS
   write(6,*)'gather skeleton...';call flush(6)
   call gather_skeleton

! must add a flag there
   if (southern) then
    write(6,*)'generate southern hemisphere...';call flush(6)
    call generate_southern_hemisphere
   else
    call donot_generate_southern_hemisphere 
   end if

   write(6,*);write(6,"(10x,' THE TOTAL NUMBER OF ELEMENTS IS ',i10 )"), neltot; write(6,*)

   call generate_serendipity(npointot,neltot,sg,zg)

end subroutine generate_skeleton
!-----------------------------------------------------------------------------
!
!dk generate_serendipity------------------------------------------------------
subroutine generate_serendipity(npoin,nel,sg,zg)
! 08/09/2004: generates mesh database compatible w/ former version of mesher
!            (element topology defined by 8 control points instead of 4)

 use numbering
 use data_grid
 use data_diag

  integer, intent(in) :: npoin, nel
  double precision, intent(in) :: sg(4,nel),zg(4,nel)
  double precision, dimension(:,:), allocatable :: sg2,zg2
! double precision :: sg2w(8,nel),zg2w(8,nel)
! integer  lnods(8,nel)
  integer, dimension(:,:), allocatable :: lnods
  integer  iel, inode, ipt,npoin2
! Numbering arrays
  integer, dimension(:), allocatable :: ipt_iglob,el_iglob,inode_iglob
  integer :: nglob_serend
  integer, dimension(:), allocatable :: iglob_serend,loc_serend
  logical, dimension(:), allocatable :: ifseg_serend

  npoin2 = 8*npoin/4

  if (dump_mesh_info_screen) write(6,*) ' NPOIN 2 is ', npoin2,nel*8
  allocate(sg2(8,nel),zg2(8,nel))

! The even coordinate indices below are linearly interpolated, i.e. DO NOT 
! represent the correct location for any spheroidally shaped element, 
! only for the linear shapes at the center (the only case for which 
! these serendipity nodes are actually needed)!
!
  do iel = 1, nel
   sg2(1,iel) = sg(1,iel) ; zg2(1,iel) = zg(1,iel)
   sg2(3,iel) = sg(2,iel) ; zg2(3,iel) = zg(2,iel)
   sg2(5,iel) = sg(3,iel) ; zg2(5,iel) = zg(3,iel)
   sg2(7,iel) = sg(4,iel) ; zg2(7,iel) = zg(4,iel)

   sg2(2,iel) = .5d0 * ( sg2(1,iel) + sg2(3,iel) ) 
   zg2(2,iel) = .5d0 * ( zg2(1,iel) + zg2(3,iel) ) 
   sg2(4,iel) = .5d0 * ( sg2(3,iel) + sg2(5,iel) ) 
   zg2(4,iel) = .5d0 * ( zg2(3,iel) + zg2(5,iel) ) 
   sg2(6,iel) = .5d0 * ( sg2(5,iel) + sg2(7,iel) ) 
   zg2(6,iel) = .5d0 * ( zg2(5,iel) + zg2(7,iel) ) 
   sg2(8,iel) = .5d0 * ( sg2(7,iel) + sg2(1,iel) ) 
   zg2(8,iel) = .5d0 * ( zg2(7,iel) + zg2(1,iel) ) 

!  sg2(9,iel) = .5 * ( sg2(1,iel) + sg2(5,iel) ) 
!  zg2(9,iel) = .5 * ( zg2(1,iel) + zg2(5,iel) ) 

! TNM: make sure axial points are equal to zero
  if (sg2(1,iel) <= 0.1d0*dabs(sg2(2,iel)-sg2(1,iel))) then 
!    write(6,*)'AXIAL ELEMENT FIX COORDS 1:',iel,sg2(1,iel),dabs(sg2(2,iel)-sg2(1,iel))
     sg2(1,iel)=0.d0
  endif

  if (sg2(7,iel) <= 0.1d0*dabs(sg2(6,iel)-sg2(7,iel))) then 
     sg2(7,iel)=0.d0
!    write(6,*)'AXIAL ELEMENT FIX COORDS 7:',iel,sg2(7,iel),dabs(sg2(6,iel)-sg2(7,iel))
  endif

  if (sg2(8,iel) <= 0.1d0*dabs(sg2(6,iel)-sg2(7,iel))) then 
     sg2(8,iel)=0.d0
!    write(6,*)'AXIAL ELEMENT FIX COORDS 8:',iel,sg2(8,iel),dabs(sg2(6,iel)-sg2(7,iel))
  endif

  end do

! write out meshes: entire domain, central region, crust, coarsening level
 if (dump_mesh_info_files) then 
   call write_serendipity_meshes(nel,sg2,zg2)
 end if
!
  allocate(iglob_serend(npoin2)) ; iglob_serend(:) = 0 
  allocate(loc_serend(npoin2)) ;   loc_serend(:) = 0
  allocate(ifseg_serend(npoin2)) 
!
! Do we really need these ones? 
! sg2w = sg2 ; zg2w = zg2
!
  if (dump_mesh_info_screen) WRITE(6,*) 'CALLING GLOBAL NUMBERING' 
!  write(6,12)nel,maxval(sg2w),maxval(zg2w),minval(sg2w),minval(zg2w), &
!             npoin2,NDIM
!  12 format(i6,4(1pe15.5),i8,i3)

! call get_global(nel,sg2w,zg2w,iglob_serend,loc_serend,ifseg_serend,nglob_serend,npoin2,8,NDIM)
  call get_global(nel,sg2,zg2,iglob_serend,loc_serend,ifseg_serend,nglob_serend,npoin2,8,NDIM)
  if (dump_mesh_info_screen) write(6,*) 'NGLOB SERENDIPITY IS ' , nglob_serend


!21 format ( 'IEL= ',i3,x,'INODE= ',i1,x, 'IPT= ',i4,x, 'LOC= ',i8,x, ' IGLOB=', i4,2x, &
!            'IGLOB LOC =', i4,x, 2(1pe12.5,2x) )
!
!  write(6,*) MAXVAL (iglob_serend(loc_serend(:)))

  allocate (ipt_iglob(nglob_serend),el_iglob(nglob_serend),inode_iglob(nglob_serend))
 
  do iel = 1, nel
   do inode = 1, 8
    ipt = (iel-1)*8 + inode 
    ipt_iglob(iglob_serend(ipt))   = ipt 
    el_iglob(iglob_serend(ipt))    = iel
    inode_iglob(iglob_serend(ipt)) = inode 
   end do
  end do
!
  allocate(lnods(8,nel))
  do iel = 1, nel
   do inode = 1, 8
    ipt = (iel-1) * 8  + inode 
    lnods(inode,iel) = iglob_serend(ipt)
   end do
  end do

  deallocate(inode_iglob,el_iglob)
  deallocate(ifseg_serend,loc_serend,iglob_serend)
  deallocate(sg2,zg2)
!
end subroutine generate_serendipity
!-----------------------------------------------------------------------------

!dk write_serendipity_meshes--------------------------------------------------
subroutine write_serendipity_meshes(nel,sg2,zg2)
use data_bkgrdmodel
use data_diag
use data_grid, only: ri,router

integer, intent(in) :: nel
double precision, intent(in) :: sg2(8,nel),zg2(8,nel)
integer :: iel
!af test
!integer :: ipt
character(len=80) :: fname

!do iel=1,2
! do ipt = 1,8
!  write(6,*) ipt,iel,sg2(ipt,iel),zg2(ipt,iel)
! end do
!end do
!af test on rock
!  return

!  fname=diagpath(1:lfdiag)//'/global_skel.dat'
!  call out_skel(8,nel,sg2(1,1),zg2(1,1),fname)
! call out_skel(8,nel/2,sg2(1,1),zg2(1,1),fname)
! stop
!  return

  write(6,*)'writing all elements....'
  open(unit=157,file=diagpath(1:lfdiag)//'/global_skel.dat')
  do iel = 1, nel
   write(157,*) sg2(1,iel), zg2(1,iel) 
   write(157,*) sg2(2,iel), zg2(2,iel) 
   write(157,*) sg2(3,iel), zg2(3,iel) 
   write(157,*) sg2(4,iel), zg2(4,iel) 
   write(157,*) sg2(5,iel), zg2(5,iel) 
   write(157,*) sg2(6,iel), zg2(6,iel) 
   write(157,*) sg2(7,iel), zg2(7,iel) 
   write(157,*) sg2(8,iel), zg2(8,iel) 
   write(157,*) sg2(1,iel), zg2(1,iel) 
   write(157,*) 
  end do
  close(157)

!  if (bkgrdmodel=='prem') then 
   write(6,*)'writing regions of elements...'; call flush(6)
   open(unit=2157,file=diagpath(1:lfdiag)//'/foc_skel.dat')
   open(unit=3157,file=diagpath(1:lfdiag)//'/smcic_skel.dat')
   open(unit=1157,file=diagpath(1:lfdiag)//'/center_skel.dat')
   open(unit=1257,file=diagpath(1:lfdiag)//'/um_antipode.dat')
   open(unit=1357,file=diagpath(1:lfdiag)//'/midmantle.dat')
   open(unit=1457,file=diagpath(1:lfdiag)//'/uppermantle.dat')
   open(unit=1557,file=diagpath(1:lfdiag)//'/uppermantle_north.dat')

   do iel = 1, nel

! fluid core 
   if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <=ri) then !.and. &
!      sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 1200./6371.  ) then
    write(2157,*) sg2(1,iel), zg2(1,iel) 
    write(2157,*) sg2(2,iel), zg2(2,iel) 
    write(2157,*) sg2(3,iel), zg2(3,iel) 
    write(2157,*) sg2(4,iel), zg2(4,iel) 
    write(2157,*) sg2(5,iel), zg2(5,iel) 
    write(2157,*) sg2(6,iel), zg2(6,iel) 
    write(2157,*) sg2(7,iel), zg2(7,iel) 
    write(2157,*) sg2(8,iel), zg2(8,iel) 
    write(2157,*) sg2(1,iel), zg2(1,iel) 
    write(2157,*) 
   endif


! solid regions: mantle & inner core
  if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >=ri) then ! .or. &
!      sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <= 1240./6371.  ) then
    write(3157,*) sg2(1,iel), zg2(1,iel) 
    write(3157,*) sg2(2,iel), zg2(2,iel) 
    write(3157,*) sg2(3,iel), zg2(3,iel) 
    write(3157,*) sg2(4,iel), zg2(4,iel) 
    write(3157,*) sg2(5,iel), zg2(5,iel) 
    write(3157,*) sg2(6,iel), zg2(6,iel) 
    write(3157,*) sg2(7,iel), zg2(7,iel) 
    write(3157,*) sg2(8,iel), zg2(8,iel) 
    write(3157,*) sg2(1,iel), zg2(1,iel) 
    write(3157,*) 
  endif

! write only IC section
!  if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <=1400./6371. .and. sg2(5,iel)<=1200./6371. &
!       .and. abs(zg2(5,iel)-1000./6371.) <= 500./6371. .and. zg2(5,iel)>=0. ) then
   if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <=discont(ndisc)/router) then 
    write(1157,*) sg2(1,iel), zg2(1,iel) 
    write(1157,*) sg2(2,iel), zg2(2,iel) 
    write(1157,*) sg2(3,iel), zg2(3,iel) 
    write(1157,*) sg2(4,iel), zg2(4,iel) 
    write(1157,*) sg2(5,iel), zg2(5,iel) 
    write(1157,*) sg2(6,iel), zg2(6,iel) 
    write(1157,*) sg2(7,iel), zg2(7,iel) 
    write(1157,*) sg2(8,iel), zg2(8,iel) 
    write(1157,*) sg2(1,iel), zg2(1,iel) 
    write(1157,*) 
   endif

! plot only upper mantle near antipode
   if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.9d0  .and. &
       sg2(5,iel)<=0.12 .and. zg2(5,iel)<0.d0 ) then
    write(1257,*) sg2(1,iel), zg2(1,iel) 
    write(1257,*) sg2(2,iel), zg2(2,iel) 
    write(1257,*) sg2(3,iel), zg2(3,iel) 
    write(1257,*) sg2(4,iel), zg2(4,iel) 
    write(1257,*) sg2(5,iel), zg2(5,iel) 
    write(1257,*) sg2(6,iel), zg2(6,iel) 
    write(1257,*) sg2(7,iel), zg2(7,iel) 
    write(1257,*) sg2(8,iel), zg2(8,iel) 
    write(1257,*) sg2(1,iel), zg2(1,iel) 
    write(1257,*) 
    endif

!!$! plot only 200km around mid mantle coarsening near equator 
!!$   if ( abs(sqrt(sg2(5,iel)**2+zg2(5,iel)**2) - 0.73) <= 250./6371. .and. & 
!!$        abs(zg2(5,iel)) <= 250./6371.) then
!!$    write(1357,*) sg2(1,iel), zg2(1,iel) 
!!$    write(1357,*) sg2(2,iel), zg2(2,iel) 
!!$    write(1357,*) sg2(3,iel), zg2(3,iel) 
!!$    write(1357,*) sg2(4,iel), zg2(4,iel) 
!!$    write(1357,*) sg2(5,iel), zg2(5,iel) 
!!$    write(1357,*) sg2(6,iel), zg2(6,iel) 
!!$    write(1357,*) sg2(7,iel), zg2(7,iel) 
!!$    write(1357,*) sg2(8,iel), zg2(8,iel) 
!!$    write(1357,*) sg2(1,iel), zg2(1,iel) 
!!$    write(1357,*) 
!!$   endif

! plot only upper mantle
   if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.9d0) then
    write(1457,*) sg2(1,iel), zg2(1,iel) 
    write(1457,*) sg2(2,iel), zg2(2,iel) 
    write(1457,*) sg2(3,iel), zg2(3,iel) 
    write(1457,*) sg2(4,iel), zg2(4,iel) 
    write(1457,*) sg2(5,iel), zg2(5,iel) 
    write(1457,*) sg2(6,iel), zg2(6,iel) 
    write(1457,*) sg2(7,iel), zg2(7,iel) 
    write(1457,*) sg2(8,iel), zg2(8,iel) 
    write(1457,*) sg2(1,iel), zg2(1,iel) 
    write(1457,*) 
   endif

! plot only upper mantle in north
    if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.95d0 .and. &
        sg2(5,iel)<= 1000./6371. .and. zg2(5,iel)> 0.d0  ) then
     write(1557,*) sg2(1,iel), zg2(1,iel) 
     write(1557,*) sg2(2,iel), zg2(2,iel) 
     write(1557,*) sg2(3,iel), zg2(3,iel) 
     write(1557,*) sg2(4,iel), zg2(4,iel) 
     write(1557,*) sg2(5,iel), zg2(5,iel) 
     write(1557,*) sg2(6,iel), zg2(6,iel) 
     write(1557,*) sg2(7,iel), zg2(7,iel) 
     write(1557,*) sg2(8,iel), zg2(8,iel) 
     write(1557,*) sg2(1,iel), zg2(1,iel) 
     write(1557,*) 
    endif

   end do
   close(1557)
   close(1457) 
   close(1357)
   close(1257)
   close(1157)
   close(3157)
   close(2157)

!  end if

end subroutine write_serendipity_meshes
!----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal(u2,mesh,rows,filename)

 implicit none
 integer*4 :: i,rows
 real, dimension(1:rows), intent(in) :: u2
 real, dimension(1:rows) :: u1
real, dimension(1:rows,2), intent(in) :: mesh
 integer*4, dimension(1:rows*2) :: cell
 integer*4, dimension(1:rows) :: cell_type
  character (len=55) :: filename;
 character (len=50) :: ss; !stream
 
!points structure
do i=2,rows*2,2
 cell(i-1)=1;
 cell(i)=(i/2)-1;
enddo
do i=1,rows
 cell_type(i)=1
enddo

u1=real(u2)
do i=1,rows
    if (abs(u1(i))<1.e-25) u1(i)=0.0
enddo
write(6789,*)size(u1),maxval(u1),minval(u1)

 write(6,*)'computing vtk file ',trim(filename),' ...'
open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')

write(100) '# vtk DataFile Version 4.0'//char(10)
write(100) 'mittico'//char(10)
write(100) 'BINARY'//char(10)
write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
write(100) ss//char(10)
!points
do i=1,rows
write(100) mesh(i,1),0.0,mesh(i,2)
enddo
write(100) char(10)
!cell topology
write(ss,fmt='(A5,2I10)') 'CELLS',rows,rows*2
write(100) char(10)//ss//char(10)
write(100) cell
write(100) char(10)
!cell type
write(ss,fmt='(A10,2I10)') 'CELL_TYPES',rows
write(100) char(10)//ss//char(10)
write(100) cell_type
write(100) char(10)
!data
write(ss,fmt='(A10,I10)') 'CELL_DATA',rows
write(100) char(10)//ss//char(10)
write(100) 'SCALARS Displ_u1 float 1'//char(10)
write(100) 'LOOKUP_TABLE default'//char(10) !color table?
write(100) real(u1)
 close(100)
write(6,*)'...saved ',trim(filename)
end subroutine write_VTK_bin_scal
!-----------------------------------------------------------------------------





!dk def_reference_spherical_grid_discont-----------------------------------------------
subroutine def_reference_spherical_grid_discont
!
! ALEX 08/04/2004
! We make the assumption of a uniform rectangular spacing 
! in the associated cylindrical grid, which gets
! mapped into a spherical shell grid of inner 
! radius ri and outer radius ro. 
!
  use data_bkgrdmodel
  use data_grid
  use data_diag
!  use mapping_spheroid
  use analytic_spheroid_mapping

! use prem 
  double precision, dimension(8,2) :: crd_control_nodes(8,2)
  double precision, dimension(:), allocatable :: dz
  integer :: npts
  integer :: iz
!
! FIRST DEFINE PARAMETERS FOR ASSOCIATED CYLINDRICAL/CARTESIAN GRID
!  open(20,file='test_dbprem.dat')
!  read(20,*) ns
!  read(20,*) nz
!  read(20,*) ri
!  read(20,*) ro
!  allocate(dz(1:nz))
!  do iz = 1, nz
!   read(20,*) dz(iz)
!  end do
!  close(20)

  ns=ns_glob ! inherited from discont_meshing routine
  nz=nz_glob ! inherited from discont_meshing routine 
  ri=rmin/router
  ro=1.
  allocate(dz(1:nz))
  do iz = 1, nz
    dz(iz)=2.d0*dz_glob(nz-iz+1)/(router-rmin)
  end do
!
  if (dump_mesh_info_screen) then
   write(6,*)'ns,nz,ri,ro:', ns,nz,ri,ro
   write(6,*)'dz', dz(:)
   write(6,*)'SUM(dz):', SUM(dz)
  end if
!
! Total number of points (soon to be corners)  in the non-coarsened shell
  npts=(ns+1)*(nz+1)
!
  if (dump_mesh_info_screen) then
   write(6,*)'CHECK PARAMS 4 ri,ro,router,nz,ns:',ri,ro,router,nz,ns;call flush(6)
  end if
  allocate(s_unif(npts),z_unif(npts)) ; s_unif(:) = 0.d0 ; z_unif(:) = 0.d0
  allocate(crd_grdc(1:ns+1,1:nz+1,2),crd_grds(1:ns+1,1:nz+1,2))
  crd_grdc(:,:,:) = 0.d0 ;  crd_grds(:,:,:) = 0.d0
!
! Define coordinates in the parent square
  call def_ref_cart_coordinates_discont(ns,nz,crd_grdc,dz)
! In order to use analytic mapping to get coordinates 
! 1) Define control nodes to define shell geometry (Northern H)
  call def_control_nodes(crd_control_nodes,ri,ro)
! 2) Use the map_spheroid function
  call def_mapped_coordinates(ns,nz,crd_grds,crd_grdc,crd_control_nodes)
! At this stage we know the coordinates in the physical domain of the
! 8 control points that define a spectral element which belongs to the
! spherical-shelly part of the mesh
! 3) Use lexicographical ordering to create 1d arrays of coordinates
  call def_global_coordinates(npts,s_unif,z_unif,ns,nz,crd_grds)
! should probably add a possible geometrical diag output at this stage
! ...
end subroutine def_reference_spherical_grid_discont
!--------------------------------------------------------------------------
!
!dk def_global_coordinates-----------------------------------
subroutine def_global_coordinates(npts,s,z,ns1,nz1,crd)
  integer, intent(in) :: npts, ns1, nz1
  double precision, dimension(1:npts), intent(out) :: s, z
  double precision, dimension(1:ns1+1,1:nz1+1,2), intent(in) :: crd
!
  integer :: ipt, is, iz
!
  do iz = 1,nz1+1
   do is = 1,ns1+1
    ipt = uniform_nodenumber(is,iz,ns1)
    s(ipt) = crd(is,iz,1) 
    z(ipt) = crd(is,iz,2)
   end do
  end do
!
end subroutine def_global_coordinates
!------------------------------------------------------------
!
!dk def_mapped_coordinates-----------------------------------
subroutine def_mapped_coordinates(ns1,nz1,crds,crdc,crd_cont)
! use mapping_spheroid
  use analytic_spheroid_mapping

 integer, intent(in) :: ns1, nz1

 double precision, dimension(1:ns1+1,1:nz1+1,2), intent(out) :: crds
 double precision, dimension(1:ns1+1,1:nz1+1,2), intent(in) :: crdc
 double precision, dimension(8,2), intent(in) :: crd_cont
!
 integer :: is, iz ; double precision :: xi, eta

!write(6,*)'CHECK PARAMS 1 nz,ns:',nz1,ns1;call flush(6)

 do iz = 1,nz1+1
  do is = 1,ns1+1
   xi = crdc(is,iz,1) ; eta = crdc(is,iz,2)
   crds(is,iz,1) = map_spheroid(xi,eta,crd_cont,1)
   crds(is,iz,2) = map_spheroid(xi,eta,crd_cont,2)
  end do
 end do
!
end subroutine def_mapped_coordinates
!---------------------------------------------------------
!
!dk def_control_nodes----------------------------------
subroutine def_control_nodes(crd,ri1,ro1)
 double precision, intent(in) :: ri1, ro1
 double precision, dimension(8,2), intent(out) :: crd
!hemispherical case
!write(6,*)'CHECK PARAMS 2 ri1,ro1:',ri1,ro1;call flush(6)
 crd(1,1) = 0.d0             ; crd(1,2) = ri1    
 crd(2,1) = ri1*.5*dsqrt(2.d0); crd(2,2) = ri1*.5d0*dsqrt(2.d0)
 crd(3,1) = ri1            ; crd(3,2) = 0.d0
 crd(4,1) = .5d0*(ri1+ro1)   ; crd(4,2) = 0.d0
 crd(5,1) = ro1            ; crd(5,2) = 0.d0
 crd(6,1) = ro1*.5*dsqrt(2.d0); crd(6,2) = ro1*.5d0*dsqrt(2.d0)
 crd(7,1) = 0.d0             ; crd(7,2) = ro1   
 crd(8,1) = 0.d0             ; crd(8,2) = .5d0*(ro1+ri1) 
end subroutine def_control_nodes
!-------------------------------------------------------
!dk estimate_aspect_ratio--------------------------------------------------
subroutine estimate_aspect_ratio
  use data_grid
  use data_diag
  double precision :: s1,z1,s2,z2,hr
  integer :: iz

  allocate(aspect_ratio(nz)) ; aspect_ratio(:) = 0.
!
  if (dump_mesh_info_files) then
   open(unit=16,file=diagpath(1:lfdiag)//'/fort.16')
  end if
  do iz=1,nz
    s1 = .5d0*(crd_grds(1,iz,1)+crd_grds(1,iz+1,1))
    s2 = .5d0*(crd_grds(2,iz,1)+crd_grds(2,iz+1,1))
    z1 = .5d0*(crd_grds(1,iz,2)+crd_grds(1,iz+1,2))
    z2 = .5d0*(crd_grds(2,iz,2)+crd_grds(2,iz+1,2))
    hr = dsqrt( (s2-s1)**2 + (z2-z1)**2 ) ! WIDTH IS COMPUTED
    s1 = .5d0*(crd_grds(1,iz  ,1)+crd_grds(2,iz  ,1))
    s2 = .5d0*(crd_grds(1,iz+1,1)+crd_grds(2,iz+1,1))
    z1 = .5d0*(crd_grds(1,iz  ,2)+crd_grds(2,iz  ,2))
    z2 = .5d0*(crd_grds(1,iz+1,2)+crd_grds(2,iz+1,2))
    hr = hr / dsqrt( (s2-s1)**2 + (z2-z1)**2 ) ! HR = WIDTH / HEIGHT
    if (dump_mesh_info_files) write(16,*) iz, hr
    aspect_ratio(iz) = hr
  end do
  if (dump_mesh_info_files) close(16) 
  
!
end subroutine estimate_aspect_ratio
!-------------------------------------------------------------------------
!
!dk define_spherical_shell------------------------------------------------
subroutine define_spherical_shell

  use data_bkgrdmodel
  use data_grid
  use data_diag
  use data_coarse
  implicit none
  integer :: nel
!
  integer :: inode, iel, ipt,ipto,ic
! This will have to be parameterised later on
!
!  open(30,file='test_icprem.dat')
!  read(30,*) nc
!! Numbering from the top
!  allocate (iclev(0:nc+1))
!! ALEX GRENOBLE
!  do ic = 0, nc+1
!   read(30,*) iclev(ic)
!  end do
!  close(30)

  nc=nc_glob
  allocate (iclev(0:nc+1))
  iclev(0)=nz_glob+1
  do ic = 1, nc
   iclev(ic)=iclev_glob(ic)
  end do
  iclev(nc+1)=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPATIBILITY CONDITIONS HAVE TO BE DEFINED HERE
! IS NS a multiple of 2**nc ? (an even multiple!)
! IS NZ compatible with the coarsening, etc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute new number of elements

  call compute_nelclean(nel,nc,iclev,ns,nz)
! define elements topology
  allocate(lnodeso(4,nel)) ; lnodeso(:,:) = 0
! array to characterize element geometry
  allocate(eltypeo(nel)) 
! call define_lnodes

  allocate(coarsingo(nel))
  coarsingo = .false.

  call define_lnodesclean(nel,lnodeso,eltypeo,coarsingo,nc,iclev,ns,nz)
!
  nelo  = nel
! Fill so and zo arrays 
  allocate(so(4*nelo),zo(4*nelo)) ; so(:) = 0. ; zo(:) = 0.
  do iel = 1, nel
   do inode = 1, 4
    ipt = lnodeso(inode,iel)
    ipto = (iel-1)*4 + inode
    so(ipto) = s_unif(ipt)
    zo(ipto) = z_unif(ipt)
   end do
  end do
! gnuplot dump
  if (dump_mesh_info_files) then
   open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="UNKNOWN",POSITION="REWIND")
   do iel = 1, nel
    do inode = 1, 4
     ipt = lnodeso(inode,iel)
     ipto = (iel-1)*4 + inode
     write(3,*) so(ipto),zo(ipto)
    end do
    ipt = lnodeso(1,iel)
    ipto = (iel-1)*4 + 1
    write(3,*) so(ipto),zo(ipto)
    write(3,*)
   end do
  end if
!
end subroutine define_spherical_shell
!--------------------------------------------------------------------------
!
!dk uniform_nodenumber-----------------------------------------------------
integer function uniform_nodenumber(is,iz,ns,istart)
!uses lexicographical order

  integer, intent(in) :: is,iz,ns
  integer, intent(in), OPTIONAL :: istart
  integer :: i0
  i0 = 0 ; if (PRESENT(istart)) i0 = istart
  uniform_nodenumber = (iz-1)*(ns+1) + is + i0
!
end function uniform_nodenumber
!--------------------------------------------------------------------------
!
!dk compute_nelclean--------------------------------------------------
subroutine compute_nelclean(nel,nc1,iclev1,ns1,nz1)
use data_diag
!Returns the number of elements defining the new spherical shell grid. 
  integer, intent(out) :: nel
  integer, intent(in) :: nc1
  integer, dimension(0:nc1+1) :: iclev1
  integer, intent(in) :: ns1, nz1
  integer :: iz,ic,icc,icold,nelregion(nc1+1),nelabove,nelregionsum

!!$  write(6,*)'NELCLEANINTRO ns,nz,nc:',ns1,nz1,nc1
!!$  do icc=0,nc1+1
!!$     write(6,*)'NELCLEANINTRO iclev:',iclev1(icc)
!!$  enddo
!!$
!!$! TARJE SEPT 2006: new count of nel
!!$
!!$  icc=0
!!$  neltmp(0) = (iclev1(icc)-iclev1(icc+1)-1)*ns1
!!$  write(6,*)'NELCLEANINTRO neltmp:',iclev1(icc),iclev1(icc+1),neltmp(icc) 
!!$  neltmp(icc) = neltmp(icc) + 3/2*ns1
!!$  write(6,*)'NELCLEANINTRO neltmp coarsening:',neltmp(icc)
!!$  do icc=1,nc1
!!$    ic=2**(icc)
!!$    neltmp(icc) = (iclev1(icc)-iclev1(icc+1)-1)*ns1/ic 
!!$    write(6,*)'NELCLEANINTRO neltmp:',iclev1(icc),iclev1(icc+1),neltmp(icc)
!!$    if (icc <nc1) then 
!!$       neltmp(icc) = neltmp(icc) + 3/2*ns1/ic
!!$    write(6,*)'NELCLEANINTRO neltmp coarsening:',neltmp(icc)
!!$    endif
!!$ enddo
!!$
!!$ write(6,*)'NELCLEANINTRO nel sum:',sum(neltmp)

  nel = 0
  nelregion(:)=0
!
  icold = 0
!
  if ( nc1 == 0 ) then 
    nel = ns1*nz1
    return
  end if
 
!=================================
   icc =1
   do iz = nz1, iclev1(icc),-1
!=================================
    if ( iz >  iclev1(icc) ) nel = nel + ns1  
    if ( iz == iclev1(icc) ) nel = nel + 3*ns1/2
   end do
   icold = 1
   nelregion(1)=nel
   nelabove=nel
!    write(6,*)'NELCLEANDOMAIN       :',icc,ic,nel; call flush(6)
!=================================
  do icc = 2,nc1
   ic = 2**(icc-1)
!  do iz = iclev1(icc-1)-2*icold, iclev(icc),-ic
   do iz = iclev1(icc-1)-2, iclev1(icc),-1
!=================================
    if ( iz >  iclev1(icc) ) nel = nel + ns1/ic ! TEST 
    if ( iz == iclev1(icc) ) nel = nel + 3*ns1/(2*ic)
   end do
   nelregion(icc)=nel-nelabove
   nelabove=nel
   icold = ic
!    write(6,*)'NELCLEANDOMAIN inside:',icc,ic,nel-3*ns1/(2*ic); call flush(6)
!    write(6,*)'NELCLEANDOMAIN coarse:',icc,ic,3*ns1/(2*ic); call flush(6)
!    write(6,*)'NELCLEANDOMAIN       :',icc,ic,nel; call flush(6)
  end do 

!=================================
   icc = nc1+1
   ic = 2**(icc-1)
!  do iz = iclev1(icc-1)-ic, 1,-ic
   do iz = iclev1(icc-1)-2, 1,-1
!=================================
      nel = nel + ns1/ic ! TEST 
   end do

!    write(6,*)'NELCLEANDOMAIN       :',icc,ic,nel; call flush(6)

  nelregion(icc)=nel-nelabove
  if (dump_mesh_info_screen) write(6,*) 'nel =',nel, 'instead of', nz1*ns1
  nelregionsum=sum(nelregion)
  if (dump_mesh_info_screen) write(6,*) 'SUM # EL. ALL SPHERICAL REGIONS: ',nelregionsum
  if (dump_mesh_info_screen) then 
   do icc=1,nc1+1
      write(6,*)'# EL. in REGION ',icc,' :',nelregion(icc),' (',& 
                 real(nelregion(icc))/real(nelregionsum)*100.,' percent)'
   enddo   
  end if
  if (dump_mesh_info_files) then
   open(unit=9999,file=diagpath(1:lfdiag)//'/fort.9999')
    do icc=1,nc1+1
     write(9999,*)icc,nelregion(icc)
    enddo   
   close(9999)
  end if

end subroutine compute_nelclean
!----------------------------------------------------------------

!dk define_lnodesclean-----------------------------------------------------
subroutine define_lnodesclean(nel,lnodes,eltype,coarsingloc,nc,iclev,ns,nz)
  
! The topology of the spherical shell and the array characterizing
! the geometry of each element is defined here.

  integer, intent(in) :: nel
  integer, dimension(4,nel),intent(out) :: lnodes
  character(len=6), dimension(nel), intent(out) :: eltype
  logical, dimension(nel), intent(out) :: coarsingloc
  integer, intent(in) :: nc
  integer, dimension(0:nc+1),intent(in) :: iclev
  integer, intent(in) :: ns,nz
  integer :: iz,ic,icc,iel,is,icold

!write(6,*)'LNODES: nel,ns,nz',nel,ns,nz; call flush(6)
  
  iel = 0 
! no coarsening case 
  if (nc == 0) then
   do iz = nz,1,-1
    do is = 1, ns
     iel = iel + 1
     lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
     lnodes(2,iel) = uniform_nodenumber(is+1,iz,ns)
     lnodes(3,iel) = uniform_nodenumber(is+1,iz+1,ns)
     lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns) 
     eltype(iel) = 'curved'
    end do
   end do  
   return
  end if
!
! Going from top to bottom 
!  do icc = 1,1

    icc = 1
    ic = 1
    do iz = nz, iclev(icc),-ic
    if ( iz >  iclev(icc) ) then 
     do is = 1, ns, ic
      iel = iel + 1
      lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
      lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns)
      lnodes(3,iel) = uniform_nodenumber(is+ic,iz+ic,ns)
      lnodes(4,iel) = uniform_nodenumber(is   ,iz+ic,ns)
      eltype(iel) = 'curved'
     end do
    elseif ( iz == iclev(icc) ) then 
!   This takes care of the "conformal mortars"
!   This case here corresponds to the first remeshing/coarsening
!   when going down.
     do is = 1, ns, ic*4
      
         iel = iel + 1                        
         lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns   )
         lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns   )
         lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1,ns)
         lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns)        
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!      
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is     ,iz-1,ns )
         lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
         lnodes(3,iel) = uniform_nodenumber(is+  ic,iz,ns   )
         lnodes(4,iel) = uniform_nodenumber(is     ,iz,ns   )
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!         
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is  +ic,iz,ns   )
         lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
         lnodes(3,iel) = uniform_nodenumber(is+2*ic,iz+1,ns )
         lnodes(4,iel) = uniform_nodenumber(is  +ic,iz+1,ns )
         eltype(iel) = 'semino'
         coarsingloc(iel)=.true.
!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz- 1,ns)
         lnodes(2,iel) = uniform_nodenumber(is+3*ic,iz,ns   )
         lnodes(3,iel) = uniform_nodenumber(is+3*ic,iz+ 1,ns)
         lnodes(4,iel) = uniform_nodenumber(is+2*ic,iz+ 1,ns)
         eltype(iel) = 'semino'
         coarsingloc(iel)=.true.
!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
         lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz-1,ns )
         lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz ,ns  )
         lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz ,ns  )
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.

     end do 
    end if
   end do
!  end do
!
!    write(6,*)'LNODESDOMAIN       :',icc,ic,iel; call flush(6)

! intermediate domains, bound by coarsening levels above & below
  icold = 1
  do icc = 2, nc
   ic = 2**(icc-1)
   do iz = iclev(icc-1)-2, iclev(icc),-1 
    if ( iz >  iclev(icc) ) then 
     do is = 1, ns, ic
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
         lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns)
         lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is   ,iz+1 ,ns)
         eltype(iel) = 'curved'
     end do
!    write(6,*)'LNODESDOMAIN inside:',icc,ic,iel; call flush(6)
    elseif ( iz == iclev(icc) ) then 
!    Remeshing 
     do is = 1, ns, ic*4
!        
!        write(6,*)'LNODESFINDING is:',is,ns,ic*4
         iel = iel + 1                       
         lnodes(1,iel) = uniform_nodenumber(is   ,iz   ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+ic,iz   ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is   ,iz+1 ,ns)        
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!      
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is     ,iz-1 ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+  ic,iz   ,ns)
         lnodes(4,iel) = uniform_nodenumber(is     ,iz   ,ns)
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!         
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is  +ic,iz   ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+2*ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is  +ic,iz+1 ,ns)
         eltype(iel) = 'semino'
         coarsingloc(iel)=.true.
!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is+2*ic,iz+1 ,ns)
         eltype(iel) = 'semino'         
         coarsingloc(iel)=.true.

!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz-1 ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
         lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
         eltype(iel) = 'curved'
         coarsingloc(iel)=.true.
!
         iel = iel + 1
         lnodes(1,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
         lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
         lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz+1 ,ns)
         lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
         eltype(iel) = 'curved'         
         coarsingloc(iel)=.true.

     end do 
!    write(6,*)'LNODESDOMAIN coarse:',icc,ic,iel; call flush(6)
    end if
   end do
   icold = ic
!    write(6,*)'LNODESDOMAIN       :',icc,ic,iel; call flush(6)
  end do
! Last series of layer after last remeshing
!  do icc = nc+1, nc+1
   icc = nc + 1
    ic = 2**(icc-1)
    do iz = iclev(icc-1)-2, 1,-1
     do is = 1, ns, ic
       iel = iel + 1
       lnodes(1,iel) = uniform_nodenumber(is   ,iz ,ns)
       lnodes(2,iel) = uniform_nodenumber(is+ic,iz ,ns)
       lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1,ns)
       lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns)
       eltype(iel) = 'curved'
     end do
    end do
!  end do
!     write(6,*)'LNODESDOMAIN       :',icc,ic,iel; call flush(6)

!    write(6,*)'defined lnodes cleanly.';call flush(6)

end subroutine define_lnodesclean
!----------------------------------------------------------------

!dk define_central_region--------------------------------------------------
subroutine define_central_region

  use data_bkgrdmodel
  use data_coarse
  use data_diag
  use data_grid
  use global_parameters
  use data_spec
  use data_mesh
!
  integer :: nr, nex, nzs, nzss
  double precision, dimension(:,:,:), allocatable:: crd_cyl
  double precision, dimension(:,:,:), allocatable:: crd_cent
  double precision ric, roc
  double precision, dimension(8,2):: crd_control_nodes
  integer npts,nelc,nelctot,nptscc
!  double precision, dimension(:), allocatable :: s_unif, z_unif
  double precision, dimension(:), allocatable :: scc,zcc
!
  integer  ipt,inode,iel,is,iz,ipti,iptsq,iptbuf,ieltest
  double precision ::  test
!
! TNM GRE
! double precision :: minh_ic,maxh_ic ! DEFINED IN data_bkgrdmodel
! /TNM GRE

! NEW METHOD
  integer :: ix,iy,maxy
  double precision :: rad,angle,p,x,y
!/ NEW METHOD

! number of latitude slices at shell - central region boundary (termed ib)
  ns_ib = ns / 2**nc  
  if (dump_mesh_info_screen) write(6,*) 'ns_ib is ', ns_ib

! define dimensions of central square (just one!) 
  lsq = lsq_fac * ri

! define number of divisions in one direction for central square
  ndivs = ns_ib/2

! compute number of necessary extra coarsenings nex
  nex = ns_ib/(2*ndivs)-1
  if (dump_mesh_info_screen) write(6,*) ' NEX IS ', nex
!
! compute number of radial levels from center of the earth to 
! shell boundary (equidistant spacing)
! PREM CASE
! dr = radius(2) - radius(1) ; write(6,*) ' DR IS ' ,dr
! nzs = int((ri-lsq)/dr) 
! ALEX GRE; af 2009, 4 years later, that seems sensible but I do not remember why I 
! picked this value. Might be useful to keep what follows compatible with nzs=/1
  nzs = 1 
! 
!
! check 
  test = (dsqrt(2.d0)*lsq) / (.9d0 *(lsq + (ri-lsq)/dble(nzs)))
  if ( test >= 1. ) then 
   write(6,*) ' STOP: failed test for generation of central region'
   write(6,*) ' test is ',test
   stop
  end if
!
! nzs = int(real(ndivs)*(ri-lsq)/lsq) 
! nzs = 2
  nr  = ndivs + nzs  ! will come from xtoto (newdiscont_meshing)
  if (dump_mesh_info_screen) write(6,*) 'nr is ', nr   
!
! roc = ri ; ric = roc - real((nzs-1))*(roc-lsq)/dble(nzs)
! Assuming nzs = 1
  roc = ri ; ric = roc*(dble(nr-1)/dble(nr))
  if (dump_mesh_info_screen) write(6,*) 'RIC = ', RIC
  if (dump_mesh_info_screen) write(6,*) 'ROC = ', ROC 
  if ( nzs < 2 * nex + 1) then 
   write(6,*) 'incompatibility in central square between nex and nr' 
   stop
  end if
! The '+1' indicates that we want the  bottommost "spherical layer" to 
! be distinguished from the others, as it will be connected to the central square.
  
! We define the (nzs-1) spherical layers of the central region
! as for the exterior spherical shell.  
! we define, again, the size of the associated uniform grid. 
! nzss is the number of purely spherical levels
  nzss = nzs - 1 ; if (dump_mesh_info_screen) write(6,*) 'nzss = ', nzss

  if (nzss /= 0) then  
  
   npts = (nzss+1) * (ns_ib+1) 

!!
!! TNM JULY 2009: Old method, redundant...
!!
   allocate( crd_cyl(1:ns_ib+1,1:nzss+1,2)) ;  crd_cyl(:,:,:)=0.
   allocate(crd_cent(1:ns_ib+1,1:nzss+1,2)) ; crd_cent(:,:,:)=0.

   allocate(s_unif(1:npts), z_unif(1:npts))
!!$
!  define reference cylindrical coordinates
   call def_ref_cart_coordinates(ns_ib,nzss,crd_cyl,.true.)
   call def_control_nodes(crd_control_nodes,ric,roc)
!!$
!!$! 2) Use the map_spheroid function
   call def_mapped_coordinates(ns_ib,nzss,crd_cent,crd_cyl,crd_control_nodes)
!!$
!!$! 3) Use lexicographical ordering to create 1d arrays of coordinates
   call def_global_coordinates(npts,s_unif,z_unif,ns_ib,nzss,crd_cent)
!!$
!!$!quick check
!!$  if (dump_mesh_info_files) then
!!$   open(unit=4,file=diagpath(1:lfdiag)//'/fort.4')
!!$   do ipt = 1, npts
!!$    write(4,*) s_unif(ipt), z_unif(ipt)
!!$   end do
!!$  end if
!!$! 
! Again, we define the levels at which the 
! potential coarsenings will take place
   allocate(iclevc(0:nex+1)) ; iclevc(:) = 0
   iclevc(0)     = nzss + 1
   iclevc(1)     = nzss 
   iclevc(nex+1) = 1  
!!$!
!!$! compute new number of elements
   call compute_nelclean(nelc,nex,iclevc,ns_ib,nzss)
   write(6,*) 'nelc ', nelc , 'npts' 
!af
! stop
!!$! define elements topology and type
   allocate(lnodesi(4,nelc)) ; lnodesi(:,:) = 0
   allocate(eltypei(nelc),coarsingi(nelc)) 
!!$  call define_lnodesclean(nelc,lnodesi,eltypei,coarsingi,nex,iclevc,ns_ib,nzss)
! fill si and zi arrays
   neli = nelc 
   allocate(si(4*neli),zi(4*neli)) ; si(:) = 0. ; zi(:) = 0. 


   do iel = 1, nelc
    do inode = 1, 4
     ipt = lnodesi(inode,iel)
     ipti = (iel-1)*4 + inode
     si(ipti) = s_unif(ipt)
     zi(ipti) = z_unif(ipt) 
    end do
   end do  

!!$! gnuplot dump
!!$  if (dump_mesh_info_files) then 
!!$   open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',POSITION="APPEND")
!  do iel = 1, nelc
!   do inode = 1, 4
!    ipt = lnodesi(inode,iel)
!    write(3,*) s_unif(ipt),z_unif(ipt)
!    end do
!   ipt = lnodesi(1,iel)
!   write(3,*) s_unif(ipt),z_unif(ipt)
!   write(3,*)
!  end do
!!$  end if
!!$

   write(6,*) 'npts = ', npts 

  else
!af october 2009
   npts=0 
   nelc=0
  end if

! define grid points for central region 
! Points numbering and coordinates
  nptscc = (ndivs+1)**2  

  allocate(scc(nptscc+npts),zcc(nptscc+npts)) ; scc(:) = 0.d0 ; zcc(:) = 0.d0
!
!!$    scc(1:npts) = s_unif(1:npts) ; zcc(1:npts) = z_unif(1:npts)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        NEW METHOD: |x|^p + |y|^p = r^p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(s_arr(ndivs+1,ndivs+1),z_arr(ndivs+1,ndivs+1))
  s_arr=zero; z_arr=zero

  if (dump_mesh_info_screen) write(6,*)'CENTR: ndivs,nr',ndivs,nr
  if (dump_mesh_info_screen) write(6,*)'CENTR: ri',ri

  if (dump_mesh_info_files) then  
   open(unit=46,file=diagpath(1:lfdiag)//'/fort.46')
   open(unit=47,file=diagpath(1:lfdiag)//'/fort.47')
   open(unit=48,file=diagpath(1:lfdiag)//'/fort.48')
  end if
  do ix=1,2*ndivs+1
     maxy=int(min(ix,ceiling(dble(ix)/2.)+1))
     if (dump_mesh_info_screen) write(6,*)'CENTR: ix,maxy',ix,maxy
     do iy=1,maxy
        rad=dble(ix)/dble(2*ndivs+1)*(ri-maxh_icb/router)
     !  linear
        p=dble(ix)/dble(2*ndivs)+1.
     ! quadratic
!        p=(dble(ix)/dble(2*ndivs))**2+1.
     ! cubic
!        p=(dble(ix)/dble(2*ndivs))**3+1.

        if (iy>1) then
          angle=tan( pi/2.d0 * ( 1 - dble(iy-1)/dble(ix) ) )
          y=rad/(angle**p+1.d0)**(1.d0/p)
        else 
          angle=0.d0
          y=0.d0;
        endif
        x=dble(( rad**p-abs(y)**p )**(1.d0/p))

! TNM NOV 2006 checking stuff...
!        write(4141,14)ix,iy,maxy,rad,p,angle,x,y
!        14 format(3(i3),5(1pe14.5))

        ! compute s,z coordinates in rotated frame and indices
        if (mod(ix,2)==0) then

!           write(6,*)'CENTR: rad,p',rad,p
!           write(6,*)'CENTR: angle,tan',(1-dble(iy-1)/dble(ix))*90.,angle
!           write(6,*)'CENTR: y,x',y,x; call flush(6)
           !   below diagonal, s>=z
!           write(6,*)'CENTR: x,y,s,z',ix,iy,int(ix/2)+1,maxy-iy+1
           s_arr(int(ix/2)+1,maxy-iy+1)=x+y 
           z_arr(int(ix/2)+1,maxy-iy+1)=x-y
           !   above diagonal, s<z (switched indices, negative y)
           s_arr(maxy-iy+1,int(ix/2)+1)=x-y
           z_arr(maxy-iy+1,int(ix/2)+1)=x+y

!           write(6,*)'XY1:',ix,iy,int(ix/2)+1,maxy-iy+1
!           write(6,*)'XY2:',ix,iy,maxy-iy+1,int(ix/2)+1
           
           if (dump_mesh_info_files) then  
            write(46,*)y,x
            write(47,*)x+y,x-y
            write(48,*)x-y,x+y
           end if

        endif
     enddo
  enddo
  if (dump_mesh_info_files) then
   close(46)
   close(47)
   close(48)
  end if

! earth center
  s_arr(1,1)=zero
  z_arr(1,1)=zero

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! TNM MAY 2007: Apply stretching for 45-deg pathological elements
!
! This is a trial-and-error fix to the triangularly deformed elements 
! along the 45 deg diagonal inside the central cube which lead to 
! grid spacing of diagonal points (is,iz) vs. (is+1,iz+1) of about 
! a factor of 2-4 lower than the predicted/expected value (for higher 
! resolution grids at least). The solver used to blow up as soon as 
! reaching the inner core if this is left out. One should be careful 
! and examine this region for each new high resolution grid before 
! running the solver to be sure it looks properly (i.e. now crossing 
! element boundaries etc).
  if (dump_mesh_info_files) then 
   open(unit=5559,file=diagpath(1:lfdiag)//'/fort.5559') 
  end if
  do is=2,ndivs+1
     do iz=2,ndivs+1
        if ( is==iz )  then
           p=0.5d0*(s_arr(is,iz)-s_arr(is-1,iz-1))*dble(is)/dble(ndivs)
           if (dump_mesh_info_files) write(5559,*)is,iz,p,s_arr(is,iz)
           s_arr(is,iz) = s_arr(is,iz)+p
           z_arr(is,iz) = z_arr(is,iz)+p
           if (dump_mesh_info_files) write(5559,*)is,iz,p,s_arr(is,iz)
        endif
     enddo
  enddo
  if (dump_mesh_info_files) then 
   close(5559)
  end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  write(6,*)'ri,h_icb,r0:',ri,maxh_icb,router
!  write(6,*)'max s_arr,z_arr:',maxval(abs(s_arr)),maxval(abs(z_arr))

! make sure max val is ri
  s_arr=s_arr/maxval(abs(s_arr))*(ri-maxh_icb/router)
  z_arr=z_arr/maxval(abs(z_arr))*(ri-maxh_icb/router)

  if (dump_mesh_info_files) then  
   open(unit=44,file=diagpath(1:lfdiag)//'/fort.44')
   open(unit=45,file=diagpath(1:lfdiag)//'/fort.45')
  end if
  do is=1,ndivs+1
     do iz=1,ndivs+1
        ipt = uniform_nodenumber(is,iz,ndivs,npts)
        scc(ipt) = s_arr(is,iz)
        zcc(ipt) = z_arr(is,iz)
        if (dump_mesh_info_files) write(44,*) scc(ipt), zcc(ipt)
        if (dump_mesh_info_files) write(45,*) is,iz,s_arr(is,iz),z_arr(is,iz)
     enddo
  enddo
  if (dump_mesh_info_files) then  
   close(44)
   close(45)
  end if

! TNM JULY 2009: Need this for gather_skeleton
  if (nelc /=0) then ! af 10 2009
  do iel = 1, nelc
   write(6,*) 'af iel ', iel
   do inode = 1, 4
    ipt = lnodesi(inode,iel)
    ipti = (iel-1)*4 + inode
    si(ipti) = scc(ipt)
    zi(ipti) = zcc(ipt) 
   end do
  end do  
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        END NEW METHOD: |x|^p + |y|^p = r^p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! total number of elements
  nelctot = nelc + 2*ndivs + ndivs**2  
  if (dump_mesh_info_screen) write(6,*) 'nelctot ', nelctot
  nelsq = ndivs**2
  allocate(lnodessq(4, nelsq))
  allocate(eltypesq(nelsq))
  iel = 0
  do iz = 1, ndivs
   do is = 1, ndivs
    iel = iel + 1
    lnodessq(1,iel) = uniform_nodenumber(is  ,iz  ,ndivs,npts)
    lnodessq(2,iel) = uniform_nodenumber(is+1,iz  ,ndivs,npts)
    lnodessq(3,iel) = uniform_nodenumber(is+1,iz+1,ndivs,npts)
    lnodessq(4,iel) = uniform_nodenumber(is  ,iz+1,ndivs,npts)
    eltypesq(iel) = 'linear'
   end do
  end do

! FILL ARRAY SSQ AND ZSQ
  allocate(ssq(4*nelsq),zsq(4*nelsq)) ; ssq(:) = 0. ; zsq(:) = 0.
  do iel = 1, nelsq
   do inode = 1, 4
    ipt = lnodessq(inode,iel)
    iptsq = 4*(iel-1)+inode
    ssq(iptsq) = scc(ipt)
    zsq(iptsq) = zcc(ipt)
   end do
  end do
!
! gnuplot dump
  if (dump_mesh_info_files) then
   open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="OLD",POSITION="APPEND")
   do iel = 1, ndivs**2
    do inode = 1, 4
     ipt = lnodessq(inode,iel)
     write(3,*) scc(ipt),zcc(ipt)
    end do
    ipt = lnodessq(1,iel)
    write(3,*) scc(ipt),zcc(ipt)
    write(3,*)
   end do
  end if
!
! THE REST IN THIS ROUTINE IS ABOUT THE buffer layer**************************
  nelbuf = 2*ndivs
  allocate(lnodesbuf(4,nelbuf)) ; lnodesbuf(:,:) = 0
  allocate(eltypebuf(nelbuf))
  iel = 0 
  do is = 1, ndivs
   iel = iel + 1
   lnodesbuf(1,iel) = uniform_nodenumber(is,ndivs+1,ndivs,npts)
   lnodesbuf(2,iel) = uniform_nodenumber(is+1,ndivs+1,ndivs,npts)
   if ( neli > 0 ) then 
    lnodesbuf(3,iel) = uniform_nodenumber(1+  (is)*2**nex,1,ns_ib)
    lnodesbuf(4,iel) = uniform_nodenumber(1+(is-1)*2**nex,1,ns_ib)
   else ! CONNECT WITH OUTER SHELL  
    ieltest = nelo + is -2*ndivs
    lnodesbuf(3,iel) = (ieltest-1)*4 + 2 
    lnodesbuf(4,iel) = (ieltest-1)*4 + 1  
!   print*, ' CIAO CIAO 1', ieltest, nelo, zo(lnodeso(2,ieltest)),minval(zo)
   end if
   eltypebuf(iel) = 'semino' 
  end do
  do iz = 1, ndivs
   iel = iel + 1
!  lnodesbuf(1,iel) = uniform_nodenumber(ndivs+1,iz,ndivs,npts)
!  lnodesbuf(2,iel) = uniform_nodenumber(ns_ib+1-(iz-1)*2**nex,1,ns_ib)
!  lnodesbuf(3,iel) = uniform_nodenumber(ns_ib+1-(iz)*2**nex,1,ns_ib)
!  lnodesbuf(4,iel) = uniform_nodenumber(ndivs+1,iz+1,ndivs,npts)
   lnodesbuf(4,iel) = uniform_nodenumber(ndivs+1,iz,ndivs,npts)
   if ( neli > 0 ) then 
    lnodesbuf(1,iel) = uniform_nodenumber(ns_ib+1-(iz-1)*2**nex,1,ns_ib)
    lnodesbuf(2,iel) = uniform_nodenumber(ns_ib+1-(iz)*2**nex,1,ns_ib)
   else ! CONNECT WITH OUTER ELEMENT 
    ieltest = nelo -( iz - 1 )  
!   print*, ' CIAO CIAO 2'
    lnodesbuf(1,iel) = (ieltest-1)*4 + 2 
    lnodesbuf(2,iel) = (ieltest-1)*4 + 1 
   end if 
   lnodesbuf(3,iel) = uniform_nodenumber(ndivs+1,iz+1,ndivs,npts)
   eltypebuf(iel) = 'semiso'
  end do
!
! fill array sbuf and zbuf
  allocate(sbuf(4*nelbuf),zbuf(4*nelbuf)) ; sbuf(:)=0.d0 ; zbuf(:)=0.d0
  do iel = 1, nelbuf

!  print*, 'COORDINATE TEST', iel
   if (neli > 0) then 
    do inode = 1, 4
     ipt = lnodesbuf(inode,iel)
     iptbuf = 4*(iel-1) + inode
     sbuf(iptbuf) = scc(ipt)
     zbuf(iptbuf) = zcc(ipt)
    end do
   else ! CONNECT WITH OUTER SHELL 
    if ( iel < (nelbuf/2 + 1) ) then 
     do inode = 1,2
      ipt = lnodesbuf(inode,iel)
      iptbuf = 4*(iel-1) + inode
      sbuf(iptbuf) = scc(ipt)
      zbuf(iptbuf) = zcc(ipt)
     end do 
     do inode = 3,4
!     print*, ' CIAO CIAO 22'
      ipt = lnodesbuf(inode,iel)
      iptbuf = 4*(iel-1) + inode
      sbuf(iptbuf) = so(ipt) 
      zbuf(iptbuf) = zo(ipt) !; write(6,*) zo(ipt), ipt
     end do
    else 
     do inode = 1,2
!     print*, ' CIAO CIAO 222'
      ipt = lnodesbuf(inode,iel)
!     print*, ' 2222 ', ipt
      iptbuf = 4*(iel-1) + inode
      sbuf(iptbuf) = so(ipt)
      zbuf(iptbuf) = zo(ipt) !; write(6,*) zo(ipt), ipt
!     print*, ' 2223 ', ipt
     end do
     do inode = 3,4
      ipt = lnodesbuf(inode,iel)
      iptbuf = 4*(iel-1) + inode
      sbuf(iptbuf) = scc(ipt)
      zbuf(iptbuf) = zcc(ipt)
     end do 
    end if
   end if
  end do
!
 
! gnuplot dump
  if (dump_mesh_info_files) then
   open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="OLD",POSITION="APPEND")
   do iel = 1, nelbuf
    do inode = 1, 4
     ipt = lnodesbuf(inode,iel)
     iptbuf = 4*(iel-1) + inode
     write(3,*) sbuf(iptbuf),zbuf(iptbuf)
    end do
    ipt = lnodesbuf(1,iel)
    iptbuf = 4*(iel-1) + 1
    write(3,*) sbuf(iptbuf),zbuf(iptbuf)
    write(3,*)
   end do
  end if
!
end subroutine define_central_region
!----------------------------------------------------------------
!dk def_ref_cart_coordinates
subroutine def_ref_cart_coordinates(nst,nzt,crd,inner_shell)
  use data_grid
  integer, intent(in) :: nst, nzt
  double precision, dimension(1:nst+1,1:nzt+1,2), intent(out) :: crd
  logical, OPTIONAL :: inner_shell
! 
  integer :: is, iz ; double precision :: ds, dz
! uniform grid in s and z 
  if (nst/=0 .and. nzt/=0) then
     ds = 2./dble(nst)
     dz = 2./dble(nzt)
  else
     write(6,*)'nst and nzt are ZERO!'
     stop
!af 10_2009
  endif
  do iz = 1, nzt+1
   do is = 1, nst+1
    crd(is,iz,1) = -1. + dble(is-1) * ds
    crd(is,iz,2) = -1. + dble(iz-1) * dz
!   crd(is,iz,2) = ndeta(iz) PREM CASE
    if (PRESENT(inner_shell)) crd(is,iz,2) = -1. + dble(iz-1) * dz
   end do
  end do
end subroutine def_ref_cart_coordinates
!---------------------------------------------------
!
!----------------------------------------------------------------
!dk def_ref_cart_coordinates_discont
subroutine def_ref_cart_coordinates_discont(nst,nzt,crd,dz)
  use data_grid
  integer, intent(in) :: nst, nzt
  double precision, dimension(1:nst+1,1:nzt+1,2), intent(out) :: crd
  double precision, dimension(1:nzt) :: dz
! 
  integer :: is, iz ; double precision :: ds
! uniform grid in s, variable spacing in z  
  ds = 2./dble(nst)
  do is = 1, nst+1
   do iz = 1, nzt+1
    crd(is,iz,1) = -1.d0+ dble(is-1) * ds
   end do
   crd(is,1,2) = -1.d0 
   do iz = 2, nzt+1
    crd(is,iz,2) = crd(is,iz-1,2) +  dz(iz-1)
   end do
  end do
end subroutine def_ref_cart_coordinates_discont
!----------------------------------------------------------------
!
!------------------------------------------------------------------------
!dk gather_skeleton
subroutine gather_skeleton
! This routine defines global arrays to assemble skeleton. 
!
  use data_mesh
  use data_grid, ONLY: ndivs

  integer :: ipt,inode,iel,istart,is,iz
! 
  write(6,*)
  write(6,"(10x,'SKELETON INFORMATIONS (NORTHERN HEMISPHERE ONLY)')")
  write(6,*)
  write(6,"(10x,'Number of elements in the outer shell:    ',i10)")  nelo
  write(6,"(10x,'Number of elements in the inner shell:    ',i10)")  neli
  write(6,"(10x,'Number of elements in the central square: ',i10)")  nelsq
  write(6,"(10x,'Number of elements in the buffer layer:   ',i10)")  nelbuf
!
  neltot = nelo + neli + nelsq + nelbuf
  write(6,*)
  write(6,"(10x,'Total num. of elements in northern skel.: ',i10)")  neltot
  write(6,*)
  write(6,*)
! 
  npointot = 4 * neltot
  allocate(sg(npointot),zg(npointot)) ; sg(:) = 0.d0 ; zg(:) = 0.d0 

! outer shell
  istart = 1
  if(allocated(sg)) sg(istart:4*nelo) = so(1:4*nelo) 
  if(allocated(zg)) zg(istart:4*nelo) = zo(1:4*nelo) 

! inner shell
  istart = 4*nelo + 1
  if(allocated(si)) sg(istart:istart+4*neli-1) = si(1:4*neli)
  if(allocated(zi)) zg(istart:istart+4*neli-1) = zi(1:4*neli)

! central square region
  istart = 4*(nelo+neli) + 1
  if(allocated(ssq)) sg(istart:istart+4*nelsq-1) = ssq(1:4*nelsq)
  if(allocated(zsq)) zg(istart:istart+4*nelsq-1) = zsq(1:4*nelsq)

! define a mapping array for domain decomposition: is,iz ==> global iel
  allocate(central_is_iz_to_globiel(1:ndivs,1:ndivs))
  iel=0
  do iz = 1, ndivs
     do is = 1, ndivs
       iel = iel + 1
       central_is_iz_to_globiel(is,iz) = nelo+neli+iel
    enddo
  enddo
  
! buffer
  istart = 4*(nelo+neli+nelsq) + 1
  if (allocated(sbuf)) sg(istart:istart+4*nelbuf-1) = sbuf(1:4*nelbuf)
  if (allocated(zbuf)) zg(istart:istart+4*nelbuf-1) = zbuf(1:4*nelbuf)

! gather element type
  allocate(lnodesg(4,neltot)) ; lnodesg(:,:) = 0 
  allocate(eltypeg(neltot),coarsing(neltot)) ;
  if (allocated(eltypeo)) eltypeg(1:nelo) = eltypeo(1:nelo)
  if (allocated(eltypei)) eltypeg(nelo+1:nelo+neli) = eltypei(1:neli)
  if (allocated(eltypeg)) eltypeg(nelo+neli+1:nelo+neli+nelsq) = eltypesq(1:nelsq)
  if (allocated(eltypebuf)) eltypeg(nelo+neli+nelsq+1:neltot) = eltypebuf(1:nelbuf)

  coarsing = .false.
  coarsing(1:nelo) = coarsingo(1:nelo)
!  coarsing(nelo+1:nelo+neli) = coarsingi(1:neli)
!  coarsing(nelo+neli+1:nelo+neli+nelsq) = coarsingsq(1:nelsq)
!  eltypeg(nelo+neli+nelsq+1:neltot) = coarsingbuf(1:nelbuf)

!
  do iel = 1, neltot
   do inode = 1, 4
    ipt = (iel-1)*4 + inode
    lnodesg(inode,iel) = ipt
   end do
  end do 
!
end subroutine gather_skeleton

!-----------------------------------------------------------------------
!dk generate_southern_hemisphere--------------------
subroutine generate_southern_hemisphere
!
use data_mesh
use data_diag
!
integer neltot2, npointot2, ipt, iel,inode
double precision, dimension(:), allocatable :: sg2,zg2
integer, dimension(:,:), allocatable        :: lnodesg2
character(len=6), dimension(:), allocatable :: eltypeg2
logical, dimension(:), allocatable          :: coarsing2
!
  npointot2 = 2 * npointot
  neltot2 = 2*neltot
  allocate(sg2(npointot2),zg2(npointot2)) ; sg2(:)=0.d0 ; zg2(:)=0.d0
  allocate(lnodesg2(4,neltot2)) ; lnodesg2(:,:) = 0
  allocate(eltypeg2(neltot2),coarsing2(neltot2))
  do ipt = 1, npointot
   sg2(ipt) = sg(ipt)
   zg2(ipt) = zg(ipt)
   sg2(ipt+npointot) =  sg(ipt)
   zg2(ipt+npointot) = -zg(ipt)
  end do
!
  do iel = 1, neltot
!
   do inode = 1, 4
    ipt = 4*(iel-1) + inode
    lnodesg2(inode,iel) = ipt
    lnodesg2(5-inode,iel+neltot) = ipt + npointot
!   lnodesg2(inode,iel+neltot) = ipt + npointot

!   ALEX: not happy w/ this part, but db is consistent w/ solver...
!         in solver subroutine get_mesh_top, see "QUICK AND DIRTY" comment
!   if(eltypeg(iel) == 'semiso' .and. (inode == 3) ) lnodesg2(1,iel+neltot) = ipt + npointot
!   if(eltypeg(iel) == 'semiso' .and. (inode == 4) ) lnodesg2(2,iel+neltot) = ipt + npointot
   end do
!
   eltypeg2(iel) = eltypeg(iel)
   coarsing2(iel) = coarsing(iel)
   coarsing2(iel+neltot) = coarsing(iel)
!
   if     (eltypeg(iel) == 'semiso') then
    eltypeg2(iel+neltot) = 'semino'
   elseif (eltypeg(iel) == 'semino') then
    eltypeg2(iel+neltot) = 'semiso'
   else
    eltypeg2(iel+neltot) = eltypeg(iel)
   end if

  end do
! copy back into original arrays
  deallocate(sg,zg,lnodesg,eltypeg,coarsing)
  npointot = npointot2
  neltot = neltot2
  allocate(sg(npointot),zg(npointot)) ; sg(:)=0.d0 ; zg(:)=0.d0
  allocate(lnodesg(4,neltot)) ; lnodesg(:,:) = 0
  allocate(eltypeg(neltot),coarsing(neltot))
  sg(:) = sg2(:)
  zg(:) = zg2(:)
  lnodesg(:,:) = lnodesg2(:,:)
  eltypeg(:) = eltypeg2(:)
  coarsing(:)=coarsing2(:)

  deallocate(sg2,zg2,lnodesg2,eltypeg2,coarsing2)
!
end subroutine generate_southern_hemisphere
!
!--------------------------------------------------------------------------
!
!dk donot_generate_southern_hemisphere-------------------------------------
subroutine donot_generate_southern_hemisphere
!
use data_mesh
use data_diag
!
integer neltot2
integer npointot2
double precision, dimension(:), allocatable :: sg2,zg2
integer, dimension(:,:), allocatable :: lnodesg2
character(len=6), dimension(:), allocatable :: eltypeg2
integer :: ipt
integer :: iel,inode
!
  npointot2 =  npointot
  neltot2 = neltot
  allocate(sg2(npointot2),zg2(npointot2)) ; sg2(:)=0.d0 ; zg2(:)=0.d0
  allocate(lnodesg2(4,neltot2)) ; lnodesg2(:,:) = 0
  allocate(eltypeg2(neltot2))
  do ipt = 1, npointot
   sg2(ipt) = sg(ipt)
   zg2(ipt) = zg(ipt)
  end do
!
  do iel = 1, neltot
!
   do inode = 1, 4
    ipt = 4*(iel-1) + inode
    lnodesg2(inode,iel) = ipt
   eltypeg2(iel) = eltypeg(iel)
   end do
!
  end do
! copy back into original arrays
  deallocate(sg,zg,lnodesg,eltypeg)
  npointot = npointot2
  neltot = neltot2
  allocate(sg(npointot),zg(npointot)) ; sg(:)=0.d0 ; zg(:)=0.d0
  allocate(lnodesg(4,neltot)) ; lnodesg(:,:) = 0
  allocate(eltypeg(neltot))
  sg(:) = sg2(:)
  zg(:) = zg2(:)
  lnodesg(:,:) = lnodesg2(:,:)
  eltypeg(:) = eltypeg2(:)
  deallocate(sg2,zg2,lnodesg2,eltypeg2)
!
end subroutine donot_generate_southern_hemisphere
!
!==============================================================================
end module meshgen
!==============================================================================
