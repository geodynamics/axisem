!========================
module lateral_heterogeneities
!========================

use global_parameters
!use data_mesh
use data_heterogeneous
use data_io
use data_proc
use utlity, only :  compute_coordinates
use data_source, only : rot_src

implicit none

public :: compute_heterogeneities
private
contains

!----------------------------------------------------------------------------------------
subroutine compute_heterogeneities(rho,lambda,mu)

implicit none
include 'mesh_params.h'
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu

write(6,*)mynum,'read parameter file for heterogeneities: inparam_hetero'
write(6,*)
write(6,*)' !!!!!!!!! W A R N I N G !!!!!!!! '
write(6,*)'These lateral additions have not been thoroughly tested yet!'
write(6,*)

call read_param_hetero

if (add_hetero) then  ! in case the input files don't exist
   ! load and add heterogeneities
   open(unit=20000+mynum,file='Info/hetero_elements_'//appmynum//'.dat')
   if (het_format=='const') then
      call load_het_const(rho,lambda,mu) ! distinct boxes with constant perturbations
   elseif (het_format=='discr') then
      call load_het_discr(rho,lambda,mu) ! interpolate discrete model of arbitrary locations/perturbations
   elseif (het_format=='funct') then
      call load_het_funct(rho,lambda,mu) ! functional perturbations (sine, gauss, error function)
   elseif (het_format=='rndm') then 
      call load_random(rho,lambda,mu) ! add random fluctuations to radial model
   else
      write(6,*)'Unknown heterogeneity input file type!!'; stop
   endif
   close(20000+mynum)

   call plot_hetero_region_vtk(rho,lambda,mu)
endif 

end subroutine compute_heterogeneities
!-------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
subroutine read_param_hetero
implicit none
integer ::  ij,i
character(len=100)         :: junk

inquire(file="inparam_hetero", EXIST=file_exists)
if (file_exists) then 
   open(unit=91,file='inparam_hetero')
   read(91,*)het_format
   read(91,*)het_file_discr
   read(91,*)het_funct_type

   if (het_format=='discr') then 
      close(91)

   else
      ! ADDED by E.Vanacore 28/10/2011 to allow for multiple heterogeneities
      read(91,*)junk
      read(91,*)num_het
      if (lpr) write(6,*)'adding ',num_het,' regions...'
      allocate(r_het1(num_het))
      allocate(r_het2(num_het))
      allocate(th_het1(num_het))
      allocate(th_het2(num_het))
      allocate(ph_het1(num_het))
      allocate(ph_het2(num_het))
      allocate(delta_rho(num_het))
      allocate(delta_vp(num_het))
      allocate(delta_vs(num_het))
      ! ADDED by E.Vanacore 28/10/2011 to allow for multiple heterogeneities  
      ! Do loop added to put heterogeneity information into arrays     
      do ij = 1, num_het
         read(91,*)r_het1(ij),r_het2(ij)
         read(91,*)th_het1(ij),th_het2(ij)
         read(91,*)ph_het1(ij),ph_het2(ij)
         read(91,*)delta_rho(ij)
         read(91,*)delta_vp(ij)
         read(91,*)delta_vs(ij)
      end do
      ! End of edits by E.Vanacore
      close(91)

      th_het1=th_het1/180.*pi; th_het2=th_het2/180.*pi
      ph_het1=ph_het1/180.*pi; ph_het2=ph_het2/180.*pi

      delta_rho = delta_rho/100.
      delta_vp = delta_vp/100.
      delta_vs = delta_vs/100.
      if (lpr) then 
         do ij=1,num_het
            write(6,*)'Specification of heterogeneous region:',ij
            write(6,*)'Radius (lower/upper bound) [km]:',r_het1(ij)/1000.,r_het2(ij)/1000.
            write(6,*)'Colatitude (lower/upper bound) [deg]:',th_het1(ij)*180./pi,th_het2(ij)*180./pi
            write(6,*)'delta rho [%]:',delta_rho(ij)
            write(6,*)'delta vp  [%]:',delta_vp(ij)
            write(6,*)'delta vs  [%]:',delta_vs(ij)
         enddo
      endif

! need to rotate coordinates if source is not along axis (beneath the north pole)
      if (rot_src ) then 
         do i=1,num_het
            write(6,*)'Before rotation r th ph 1:',r_het1(i),th_het1(i)*180./pi,ph_het1(i)*180./pi
            write(6,*)'Before rotation r th ph 2:',r_het2(i),th_het2(i)*180./pi,ph_het2(i)*180./pi
         enddo
         call rotate_hetero(num_het,1,r_het1,th_het1,ph_het1)
         call rotate_hetero(num_het,1,r_het2,th_het2,ph_het2)
         do i=1,num_het
            write(6,*)'After rotation r th ph 1:',r_het1(i),th_het1(i)*180./pi,ph_het1(i)*180./pi
            write(6,*)'After rotation r th ph 2:',r_het2(i),th_het2(i)*180./pi,ph_het2(i)*180./pi
         enddo
      endif

      if (het_format=='const') then 
         nhet_pts=num_het*2
         allocate(rhet(1:nhet_pts),thhet(1:nhet_pts),phhet(1:nhet_pts))
         rhet(1:num_het)=r_het1(1:num_het)
         rhet(num_het+1:2*num_het)=r_het2(1:num_het)
         thhet(1:num_het)=th_het1(1:num_het)
         thhet(num_het+1:2*num_het)=th_het2(1:num_het)
         phhet(1:num_het)=ph_het1(1:num_het)
         phhet(num_het+1:2*num_het)=ph_het2(1:num_het)
      endif
   endif
else 
   if (lpr) then 
      write(6,*)'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(6,*)'want to add heterogeneity but cannot locate inparam_hetero..... IGNORED.'
      add_hetero=.false.
   endif

endif

end subroutine read_param_hetero
!------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------
subroutine rotate_hetero(n,m,r,th,ph)
implicit none

integer, intent(in) :: n,m
double precision,intent(inout), dimension(1:n,1:m) :: r,th,ph
double precision :: x_vec(3),x_vec_rot(3),r_r
integer :: i,j

write(6,*)'need to rotate the heterogeneous domain with the source....'

open(unit=23,file='Info/hetero_rotations_'//appmynum//'.dat')
do i=1,m
   do j=1,n
   write(23,*)'before rot r th ph:',r(j,i)/1000.,th(j,i)*180./pi
   
   x_vec(1)=r(j,i)*dsin(th(j,i))*dcos(ph(j,i))
   x_vec(2)=r(j,i)*dsin(th(j,i))*dsin(ph(j,i))
   x_vec(3)=r(j,i)*dcos(th(j,i))
   
   x_vec_rot=matmul(trans_rot_mat,x_vec)
   
   r_r= dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
   th(j,i) = dacos(x_vec_rot(3)/( r_r +smallval_dble) )
   ph(j,i) = dasin( x_vec_rot(2)/( r_r*dsin(th(j,i)) +smallval_dble) )
   
   write(23,*)'after rot r th ph:',r_r/1000.,th(j,i)*180./pi
   enddo
enddo 

end subroutine rotate_hetero
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
subroutine load_het_discr(rho,lambda,mu)
implicit none
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
double precision, allocatable, dimension(:) :: disconttmp
double precision :: w(4),wsum
integer :: ndisctmp,ind(4),jj,maxpts,iel,ipol,jpol,i,j
integer, allocatable :: num_het_pts_region(:),het_ind(:,:)
double precision :: s,z,r,th,r1,vptmp,vstmp,r2,r3,r4,th1,th2,th3,th4
double precision, allocatable, dimension(:) :: rmin,rmax,thetamin,thetamax
double precision, allocatable, dimension(:,:) :: shet, zhet

write(6,*)mynum,'reading discrete heterogeneity file...'
open(unit=91,file=trim(het_file_discr))
  read(91,*)num_het
  if (lpr) write(6,*)'Number of distinct-discrete regions:',num_het
  allocate(num_het_pts_region(1:num_het))
  do i=1,num_het
     read(91,*)num_het_pts_region(i)
     if (lpr) write(6,*)'Region',i,'has',num_het_pts_region(i),'points.'
  enddo
  allocate(rmin(num_het),rmax(num_het),thetamin(num_het),thetamax(num_het))
  nhet_pts=sum(num_het_pts_region)
  maxpts=maxval(num_het_pts_region)

  if (lpr) write(6,*)'max number of points in one region:',maxpts
  if (lpr) write(6,*)'total number of points:',nhet_pts
  
  allocate(rhet2(1:maxpts,1:num_het),thhet2(1:maxpts,1:num_het),phhet2(1:maxpts,1:num_het))
  allocate(delta_vs2(1:maxpts,1:num_het),delta_vp2(1:maxpts,1:num_het),delta_rho2(1:maxpts,1:num_het))
  allocate(het_ind(1:maxpts,1:num_het))
  rhet2=-5000.;thhet2=-5000.;delta_vp2=-5000.;delta_vs2=-5000.;delta_rho2=-5000.
  jj=0
  write(6,*)mynum,'read coordinates & medium properties...'
  do i=1,num_het
     do j=1,num_het_pts_region(i)
        jj=jj+1;het_ind(j,i)=jj
        read(91,*)rhet2(j,i),thhet2(j,i),phhet2(j,i),delta_vp2(j,i),delta_vs2(j,i),delta_rho2(j,i)
     enddo
  enddo
close(91)

if (lpr) write(6,*)'percent -> fraction'
delta_vp2=delta_vp2/100.; delta_vs2=delta_vs2/100.; delta_rho2=delta_rho2/100.
thhet2=thhet2*pi/180.; phhet2=phhet2*pi/180.; rhet2=rhet2*1000.

! Rotate coordinates if source is not on axis
do i=1,num_het
   rmin(i)=minval(rhet2(1:num_het_pts_region(i),i),1)
   rmax(i)=maxval(rhet2(1:num_het_pts_region(i),i),1)
   thetamin(i)=minval(thhet2(1:num_het_pts_region(i),i),1)
   thetamax(i)=maxval(thhet2(1:num_het_pts_region(i),i),1)
   write(6,*)mynum,'r min/max:',i,rmin(i)/1000.,rmax(i)/1000.
   write(6,*)mynum,'th min/max:',i,thetamin(i)/pi*180.,thetamax(i)/pi*180.

   if (rot_src) then 
      write(6,*)mynum,'rotate since source is not beneath north pole'
      call rotate_hetero(num_het_pts_region(i),1,rhet2(1:num_het_pts_region(i),i), &
                        thhet2(1:num_het_pts_region(i),i),phhet2(1:num_het_pts_region(i),i))

      rmin(i)=minval(rhet2(1:num_het_pts_region(i),i),1)
      rmax(i)=maxval(rhet2(1:num_het_pts_region(i),i),1)
      thetamin(i)=minval(thhet2(1:num_het_pts_region(i),i),1)
      thetamax(i)=maxval(thhet2(1:num_het_pts_region(i),i),1)
      write(6,*)mynum,'r min/max after rotation:',i,rmin(i)/1000.,rmax(i)/1000.
      write(6,*)mynum,'th min/max after rotation:',i,thetamin(i)/pi*180.,thetamax(i)/pi*180.
   endif
enddo

! plot discrete input file in vtk
call plot_discrete_input(num_het_pts_region)

! for plotting discrete points within heterogeneous region
rhetmin=minval(rmin,1); rhetmax=maxval(rmax,1)
thhetmin=minval(thetamin,1); thhetmax=maxval(thetamax,1)
write(6,*)'r het min/max:',rhetmin/1000.,rhetmax/1000.
write(6,*)'th het min/max:',thhetmin/pi*180.,thhetmax/pi*180.

! revert to cylindrical 
allocate (shet(1:maxpts,1:num_het),zhet(1:maxpts,1:num_het))
shet=rhet2*sin(thhet2);  zhet=rhet2*cos(thhet2)
write(6,*)mynum,'locate GLL points within heterogeneous regions & interpolate over 4 adjacent points'
ind=1;wsum=0;w=0;
do iel=1,nelem
   call compute_coordinates(s,z,r1,th1,iel,npol,npol)
   call compute_coordinates(s,z,r2,th2,iel,0,0)
   do i=1,num_het
      r=max(r1,r2); th=max(th1,th2)
      if ( r>=rmin(i) .and. th>=thetamin(i) ) then
            r=min(r1,r2); th=min(th1,th2)
            if ( r<=rmax(i) .and. th<=thetamax(i) ) then
               do ipol=0,npol
                  do jpol=0,npol
                     ! find closest 4 points of discrete heterogeneous mesh
                     call compute_coordinates(s,z,r,th,iel,ipol,jpol)
                     call find_closest4_points(s,z,num_het_pts_region(i),shet(1:num_het_pts_region(i),i),&
                                               zhet(1:num_het_pts_region(i),i),ind,w,wsum)
                     ! interpolate over 4 points
                     vptmp=sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel))/rho(ipol,jpol,iel))
                     vstmp=sqrt( mu(ipol,jpol,iel)/rho(ipol,jpol,iel))
                     rho(ipol,jpol,iel)=rho(ipol,jpol,iel)* (1. + sum(w*delta_rho2(ind,i))*wsum)
                     vptmp=vptmp* (1. + sum(w*delta_vp2(ind,i))*wsum)
                     vstmp=vstmp* (1. + sum(w*delta_vs2(ind,i))*wsum)
                     mu(ipol,jpol,iel)=rho(ipol,jpol,iel)*vstmp**2
                     lambda(ipol,jpol,iel)=rho(ipol,jpol,iel)*(vptmp**2 -2.*vstmp**2)
                  enddo
               enddo
            endif
      endif
   enddo
enddo

write(6,*)mynum,'DONE loading discrete grid'

end subroutine load_het_discr
!---------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------
subroutine find_closest4_points(s0,z0,n,s,z,ind,w,wsum)

integer, intent(in) :: n
double precision, intent(in) :: s0,z0,s(1:n),z(1:n)
integer, intent(out) :: ind(4)
double precision, intent(out) :: w(4),wsum
integer :: i,p,is1,is2,iz1,iz2
double precision :: dr(4),ds1,ds2,dz1,dz2

! choose 4 closest points
ds1= minval(s0-s(1:n),DIM=1,MASK=(s0-s(1:n))>=0.0)
ds2= minval(s(1:n)-s0,DIM=1,MASK=(s(1:n)-s0)>=0.0)
dz1= minval(z0-z(1:n),DIM=1,MASK=(z0-z(1:n))>=0.0)
dz2= minval(z(1:n)-z0,DIM=1,MASK=(z(1:n)-z0)>=0.0)

is1= minloc(s0-s(1:n),DIM=1,MASK=(s0-s(1:n))>=0.0)
is2= minloc(s(1:n)-s0,DIM=1,MASK=(s(1:n)-s0)>=0.0)
iz1= minloc(z0-z(1:n),DIM=1,MASK=(z0-z(1:n))>=0.0)
iz2= minloc(z(1:n)-z0,DIM=1,MASK=(z(1:n)-z0)>=0.0)

ind(1)=is1;ind(2)=is2;ind(3)=iz1; ind(4)=iz2
dr(1)=ds1; dr(2)=ds2; dr(3)=dz1; dr(4)=dz2
do i=1,4
   if (ind(i)==0) then
      if (i>1) then 
         ind(i)=ind(i-1)
!         write(6,*)'index decreased:',ind(i),i,ind(i-1)
         dr(i)=dr(i-1)
      else
         ind(i)=ind(i+1)
!         write(6,*)'index increased:',ind(i),i,ind(i+1)
         dr(i)=dr(i+1)
      endif
   endif
enddo

write(60+mynum,*)mynum,is1,is2,iz1,iz2
write(60+mynum,*)mynum,ds1,ds2,dz1,dz2
write(60+mynum,*)mynum,minval(s0-s),minval(abs(s0-s)),minval(s0-s,MASK=(s0-s)>=0.0)
write(60+mynum,*)mynum,'gll:',s0/1000.,z0/1000.,sqrt(s0**2+z0**2)/1000.,atan(s0/(z0+epsi))*180./pi
do i=1,4
   write(60+mynum,*)mynum,i,s(ind(i))/1000.,z(ind(i))/1000.,dr(i)/1000.
enddo
write(60+mynum,*)
p=1

! inverse distance weighting
do i=1,4
   w(i)=(dr(i))**(-p)
enddo
wsum=1./sum(w)

end subroutine find_closest4_points
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
subroutine plot_discrete_input(num_het_pts_region)

use background_models, only : velocity
use data_mesh, only : discont,bkgrdmodel

integer, intent(in) :: num_het_pts_region(num_het)
integer :: i,j,idom,icount
real, allocatable, dimension(:,:) :: meshtmp
real, allocatable, dimension(:) :: vptmp,vstmp,rhotmp
character(len=80) :: fname

allocate(vptmp(nhet_pts),vstmp(nhet_pts),rhotmp(nhet_pts),meshtmp(nhet_pts,2))
icount=0
do i=1,num_het
   do j=1,num_het_pts_region(i)
      icount = icount + 1
      idom = minloc(abs(discont-rhet2(j,i)),1)

      vptmp(icount)=velocity(rhet2(j,i),'v_p',idom,bkgrdmodel,lfbkgrdmodel)
      vptmp(icount) = vptmp(icount)* (1.+delta_vp2(j,i))

      vstmp(icount)=velocity(rhet2(j,i),'v_s',idom,bkgrdmodel,lfbkgrdmodel)
      vstmp(icount) = vstmp(icount)* (1.+delta_vs2(j,i))

      rhotmp(icount)=velocity(rhet2(j,i),'rho',idom,bkgrdmodel,lfbkgrdmodel)
      rhotmp(icount) = rhotmp(icount)* (1.+delta_rho2(j,i))

      meshtmp(icount,1) = rhet2(j,i)*sin(thhet2(j,i))
      meshtmp(icount,2) = rhet2(j,i)*cos(thhet2(j,i))
   enddo
enddo

fname=trim('Info/model_rho_discr_het'//appmynum)
call write_VTK_bin_scal_pts(rhotmp(1:nhet_pts),meshtmp(1:nhet_pts,1:2),nhet_pts,fname)

fname=trim('Info/model_vp_discr_het'//appmynum)
call write_VTK_bin_scal_pts(vptmp(1:nhet_pts),meshtmp(1:nhet_pts,1:2),nhet_pts,fname)

fname=trim('Info/model_vs_discr_het'//appmynum)
call write_VTK_bin_scal_pts(vstmp(1:nhet_pts),meshtmp(1:nhet_pts,1:2),nhet_pts,fname)

deallocate(meshtmp,vptmp,vstmp,rhotmp)

end subroutine plot_discrete_input
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
subroutine load_random(rho,lambda,mu)

use commun
use data_mesh, only : naxel, ax_el
use utlity, only :  rcoord,zcoord

implicit none 

double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
real(kind=8) :: t,decay,shift_fact,max_delta_vp,max_delta_vs,max_delta_rho
real(kind=8) :: vptmp,vstmp,rhotmp,s,z,r,th,gauss_val
integer :: iel,ipol,jpol,icount,i
real(kind=8) :: rand
real(kind=8), allocatable :: r_rad(:), rand_rad(:),r_radtmp(:), rand_radtmp(:)

write(6,*)'add random anomalies to structure'
!!$! add randomly to each 2D point : laterally heterogeneous and same random number to vp,vs,rho
!!$do iel=1,nelem
!!$   do jpol=0,npol
!!$      do ipol=0,npol
!!$         call random_number(rand)
!!$         rand = 2.*rand-1.
!!$         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
!!$         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
!!$         rho(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(1)*rand)
!!$         vptmp = vptmp*(1. + delta_vp(1)*rand)
!!$         vstmp = vstmp*(1. + delta_vs(1)*rand)
!!$         lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) *( vptmp*vptmp - two*vstmp*vstmp )
!!$         mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
!!$      enddo
!!$   enddo
!!$enddo

!!$! go along axis to find the radial profile
!!$allocate(r_radtmp(naxel*(npol+1)), rand_radtmp(naxel*(npol+1)))
!!$if (mynum==0) then 
!!$icount=0
!!$do iel=1,naxel
!!$   do jpol=0,npol
!!$      if (zcoord(0,jpol,ax_el(iel)) >=0.) then 
!!$         icount=icount+1
!!$         r_radtmp(icount) = rcoord(0,jpol,ax_el(iel))
!!$         call random_number(rand)
!!$         rand = 2.*rand-1.
!!$         rand_radtmp(icount)=rand
!!$      endif
!!$   enddo
!!$enddo
!!$endif 
!!$
!!$! broadcast the profile to all processors
!!$call broadcast_int(icount,0)
!!$write(6,*)mynum,'number of radii:',icount
!!$allocate(r_rad(icount),rand_rad(icount))
!!$do i=1,icount
!!$   call broadcast_dble(r_radtmp(i),0)
!!$   r_rad(i)=r_radtmp(i)
!!$   call broadcast_dble(rand_radtmp(i),0)
!!$   rand_rad(i)=rand_radtmp(i)
!!$enddo
!!$
!!$! add randomly to each radius, i.e. just altering the 1D background model
!!$do iel=1,nelem
!!$   do jpol=0,npol
!!$      do ipol=0,npol
!!$         i = minloc(abs(rcoord(ipol,jpol,iel)-r_rad(1:icount)),1)
!!$         rand = rand_rad(i)
!!$         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
!!$         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
!!$         rho(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(1)*rand)
!!$         vptmp = vptmp*(1. + delta_vp(1)*rand)
!!$         vstmp = vstmp*(1. + delta_vs(1)*rand)
!!$         lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) *( vptmp*vptmp - two*vstmp*vstmp )
!!$         mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
!!$      enddo
!!$   enddo
!!$enddo


! go along axis to find the radial profile, only per element (not GLL point)
allocate(r_radtmp(naxel*(npol+1)), rand_radtmp(naxel*(npol+1)))
if (mynum==0) then 
   icount=0
   do iel=1,naxel
      if (zcoord(0,npol/2,ax_el(iel)) >=0.) then 
         icount=icount+1
         r_radtmp(icount) = rcoord(0,npol/2,ax_el(iel))
         call random_number(rand)
         rand = 2.*rand-1.
         rand_radtmp(icount)=rand
      endif
   enddo
endif 

! broadcast the profile to all processors
call broadcast_int(icount,0)
write(6,*)mynum,'number of radii:',icount
allocate(r_rad(icount),rand_rad(icount))
do i=1,icount
   call broadcast_dble(r_radtmp(i),0)
   r_rad(i)=r_radtmp(i)
   call broadcast_dble(rand_radtmp(i),0)
   rand_rad(i)=rand_radtmp(i)
enddo

! add randomly to each radius, i.e. just altering the 1D background model
do iel=1,nelem
   i = minloc(abs(rcoord(npol/2,npol/2,iel)-r_rad(1:icount)),1)
   rand = rand_rad(i)
   do jpol=0,npol
      do ipol=0,npol   
         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
         rho(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(1)*rand)
         vptmp = vptmp*(1. + delta_vp(1)*rand)
         vstmp = vstmp*(1. + delta_vs(1)*rand)
         lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) *( vptmp*vptmp - two*vstmp*vstmp )
         mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
      enddo
   enddo
enddo

end subroutine load_random
!---------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
subroutine load_het_funct(rho,lambda,mu)

implicit none 

double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
real(kind=8) :: t,decay,shift_fact,max_delta_vp,max_delta_vs,max_delta_rho
real(kind=8) :: vptmp,vstmp,rhotmp,s,z,r,th,gauss_val
real(kind=8),allocatable,dimension(:) :: r_center_gauss,th_center_gauss
real(kind=8),allocatable,dimension(:) ::s_center_gauss,z_center_gauss,halfwidth_r,halfwidth_th
integer :: iel,ipol,jpol,icount,i

if (het_funct_type=='gauss') then 
! Gaussian
   decay=3.5d0
   
   allocate(r_center_gauss(num_het),th_center_gauss(num_het))
   allocate(s_center_gauss(num_het),z_center_gauss(num_het))
   allocate(halfwidth_r(num_het),halfwidth_th(num_het))
   
   do i=1,num_het
      r_center_gauss(i)=(r_het1(i)+r_het2(i))/2.
      th_center_gauss(i)=(th_het1(i)+th_het2(i))/2.*r_center_gauss(i)
      halfwidth_r(i) = abs(r_het1(i)-r_het2(i))
      halfwidth_th(i) = abs(th_het1(i)-th_het2(i))*r_center_gauss(i)
      if (lpr) write(6,*)i,'center r,th gauss [km]:',r_center_gauss(i)/1000.,&
           th_center_gauss(i)/1000.
      if (lpr) write(6,*)i,'halfwidth r,th gauss [km]:',halfwidth_r(i)/1000.,halfwidth_th(i)/1000.
   enddo
   
20 format(5(1pe12.3))
   allocate(rhet(nelem*(npol+1)**2),thhet(nelem*(npol+1)**2))
   icount=0
   rhet=0.;thhet=0.
   do iel=1,nelem
      do jpol=0,npol
         do ipol=0,npol
            do i=1,num_het
               call compute_coordinates(s,z,r,th,iel,ipol,jpol)
               gauss_val = dexp(-( decay* ( ((r-r_center_gauss(i))/halfwidth_r(i))**2 + &
                    ((th*r_center_gauss(i)-th_center_gauss(i))/halfwidth_th(i))**2 )))
               vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
               vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
               vstmp = vstmp  * (1. + delta_vs(i)*gauss_val)
               vptmp = vptmp  * (1. + delta_vp(i)*gauss_val)
               rho(ipol,jpol,iel) = rho(ipol,jpol,iel) * (1. + delta_rho(i)*gauss_val)
               lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) * (vptmp**2 - 2.*vstmp**2)
               mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * (vstmp**2)
               if (gauss_val> 0.01) then 
                  icount = icount + 1; rhet(icount)=r; thhet(icount)=th
               endif
            enddo
         enddo
      enddo
   enddo
   
!min/max of heterogeneous region
   rhetmin=minval(rhet(1:icount)); rhetmax=maxval(rhet(1:icount))
   thhetmin=minval(thhet(1:icount)); thhetmax=maxval(thhet(1:icount))
   write(6,*)mynum,'r het min/max:',rhetmin/1000.,rhetmax/1000.
   write(6,*)mynum,'th het min/max:',thhetmin*180./pi,thhetmax*180./pi
   
else
   write(6,*)'function type ',het_funct_type,' not implemented yet!'
   stop
! Error fct

! Sinus

endif

end subroutine load_het_funct
!----------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------
subroutine load_het_const(rho,lambda,mu)

double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
double precision :: s,z,rmin,rmax,thetamin,thetamax,r1,vptmp,vstmp,r2,r3,r4,th1,th2,th3,th4
integer :: iel,ipol,jpol,iidom,ieldom(nelem),domcount(ndisc),iel_count,ij,jj
logical :: foundit

if (het_format/='const') then 
   r_het1(1:num_het)=rhet(1:num_het)
   r_het2(1:num_het)=rhet(num_het+1:2*num_het)
   th_het1(1:num_het)=thhet(1:num_het)
   th_het2(1:num_het)=thhet(num_het+1:2*num_het)
   ph_het1(1:num_het)=phhet(1:num_het)
   ph_het2(1:num_het)=phhet(num_het+1:2*num_het)
   deallocate(rhet,thhet,phhet)
endif

jj=0
 !========================
 do iel=1,nelem
!========================
    foundit = .false.
     ! Edited by E. Vanacore to allow for multiple heterogeneities on 28/10/2011
     ! do loop/arrays added
     call compute_coordinates(s,z,r1,th1,iel,0,0);   call compute_coordinates(s,z,r2,th2,iel,0,npol)
     call compute_coordinates(s,z,r3,th3,iel,npol,0);  call compute_coordinates(s,z,r4,th4,iel,npol,npol)
     rmin=1.001*min(r1,r2,r3,r4); thetamin=1.001*min(th1,th2,th3,th4)
     rmax=0.999*max(r1,r2,r3,r4);  thetamax=0.999*max(th1,th2,th3,th4)
     do ij = 1, num_het
        if (rmin>=r_het1(ij) .and. thetamin>=th_het1(ij)) then
           if (rmax<=r_het2(ij) .and. thetamax<=th_het2(ij)) then 
              jj=ij;   foundit = .true.;   iel_count=iel_count + 1
              write(6,*)mynum,'found element inside hetero region:',iel,jj
              write(6,*)mynum,'r,th min',rmin/1000.,thetamin*180./pi
              write(6,*)mynum,'r,th max',rmax/1000.,thetamax*180./pi
           endif
        endif
     end do

     ! change some elements in a given region (see inparam_hetero
     ! Edited by E. Vanacore on 28/10/2011 to allow for multiple heterogeneities
     if (foundit) then
        do ipol=0,npol
           do jpol=0,npol
              vptmp = dsqrt( ( lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel) )
              vstmp = dsqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
!              write(20000+mynum,13)iel_count,iel,rmin/1000.,thetamin*180./pi,rho(ipol,jpol,iel),vptmp,vstmp
              rho(ipol,jpol,iel)=rho(ipol,jpol,iel)+delta_rho(jj)*rho(ipol,jpol,iel)
              vptmp = vptmp + delta_vp(jj)*vptmp
              vstmp = vstmp + delta_vs(jj)*vstmp
 !             write(20000+mynum,13)ipol,jpol,rmin/1000.,thetamin*180./pi,rho(ipol,jpol,iel),vptmp,vstmp
              lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                   ( vptmp*vptmp - two*vstmp*vstmp )
              mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
           enddo
        enddo
     endif

 !========================
 enddo
!========================

rhetmin = 0.9*minval(r_het1,1); rhetmax = 1.1*maxval(r_het2,1)
thhetmin = 0.9*minval(th_het1,1); thhetmax = 1.1*maxval(th_het2,1)

13 format(i4,i8,5(1pe13.4))

end subroutine load_het_const
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine plot_hetero_region_vtk(rho,lambda,mu)

double precision, dimension(0:npol,0:npol,nelem), intent(in) :: rho 
double precision, dimension(0:npol,0:npol,nelem), intent(in) :: lambda,mu
real, dimension(:), allocatable :: vp_all,vs_all,rho_all
real, dimension(:,:), allocatable :: mesh2
character(len=80) :: fname
double precision :: s,z,r,th
integer :: iel,ipol,jpol,icount

write(6,*)'plotting heterogeneous region in pointwise vtk'

allocate(mesh2(nelem*npol**2,2),vp_all(nelem*npol**2),vs_all(nelem*npol**2),rho_all(nelem*npol**2))

if (lpr) then
   write(6,*)'Heterogeneous region rmin,rmax [km]:',rhetmin/1000.,rhetmax/1000.
   write(6,*)'Heterogeneous region thmin,thmax [deg]:',thhetmin*180./pi,thhetmax*180./pi
endif
icount=0

do iel=1,nelem
   do ipol=0,npol
      do jpol=0,npol
         call compute_coordinates(s,z,r,th,iel,ipol,jpol)
         if ( r>=rhetmin .and. r<=rhetmax .and. th>=thhetmin .and. th<=thhetmax) then 
            icount=icount+1
            mesh2(icount,1)=real(s)
            mesh2(icount,2)=real(z)
            vp_all(icount) = sqrt( (lambda(ipol,jpol,iel)+2.*mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel)  )
            vs_all(icount) = sqrt( (mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel)  )
            rho_all(icount) = rho(ipol,jpol,iel) 
            write(666+mynum,20)r/1000.,th*180./pi,0.0,-vp_all(icount)*0.2,-vs_all(icount)*0.4,rho(ipol,jpol,iel)*0.3
         endif
      enddo
   enddo
enddo
20 format(3(1pe12.3),3(1pe10.1))

write(6,*)mynum,'number of points inside heterogeneous region:',icount

fname=trim('Info/model_vp_gll_het'//appmynum)
call write_VTK_bin_scal_pts(real(vp_all(1:icount)),real(mesh2(1:icount,1:2)),icount,fname)

fname=trim('Info/model_vs_gll_het'//appmynum)
call write_VTK_bin_scal_pts(real(vs_all(1:icount)),real(mesh2(1:icount,1:2)),icount,fname)

fname=trim('Info/model_rho_gll_het'//appmynum)
call write_VTK_bin_scal_pts(real(rho_all(1:icount)),real(mesh2(1:icount,1:2)),icount,fname)

deallocate(mesh2,vp_all,vs_all,rho_all)

end subroutine plot_hetero_region_vtk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal_pts(u2,mesh1,rows,filename)

 implicit none
 integer :: i
 integer, intent(in) :: rows
real, dimension(1:rows), intent(in) :: u2
real, dimension(1:rows) :: u1
real, dimension(1:rows,1:2), intent(in) :: mesh1
 integer, dimension(1:rows*2) :: cell
 integer, dimension(1:rows) :: cell_type
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
 !write(6789,*)size(u1),maxval(u1),minval(u1)

 !write(6,*)'computing vtk file ',trim(filename),' ...'
 open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')

 write(100) '# vtk DataFile Version 4.0'//char(10)
 write(100) 'mittico'//char(10)
 write(100) 'BINARY'//char(10)
 write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
 write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
 write(100) ss//char(10)
 !points
 do i=1,rows
    write(100) mesh1(i,1),mesh1(i,2),0.0
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
 write(100) 'SCALARS '//trim(filename)//' float 1'//char(10)
 write(100) 'LOOKUP_TABLE default'//char(10) !color table?
 write(100) real(u1)
 close(100)
 write(6,*)'...saved ',trim(filename)//'.vtk'

end subroutine write_vtk_bin_scal_pts
!-----------------------------------------------------------------------------

end module lateral_heterogeneities
