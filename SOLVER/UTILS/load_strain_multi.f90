!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stephanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

program load_snaps

implicit none

include 'mesh_params.h'

integer iproc,npts,npts_top,npts_bot,npts_meri,nphi,snapskip,snap1,snap2

real, dimension(:,:), allocatable :: coord,disp,disptot,azi_xsect2,azi_xsect1
real, dimension(:,:,:), allocatable :: azi,fluid_prefact,azi_meri
real, dimension(:), allocatable :: vp,x,y,z,vptop,vpbot,vpmeri,xtot,ytot,ztot,azi_phi_top,azi_phi_bot,azi_phi_meri
real, dimension(:), allocatable :: xtop,ytop,ztop,xbot,ybot,zbot,xmeri,ymeri,zmeri

real :: dphi,phi0,prem,r,theta_meri,smallval_meri,dr,theta0,theta_bot,theta_top,phi
integer, dimension(:), allocatable :: ind_proc_top_tmp,ind_pts_top_tmp,ind_proc_bot_tmp,ind_pts_bot_tmp
integer, dimension(:), allocatable :: ind_proc_top,ind_pts_top,ind_proc_bot,ind_pts_bot,azi_ind_top,azi_ind_bot,azi_ind_meri
integer, dimension(:), allocatable :: ind_proc_meri_tmp,ind_pts_meri_tmp,ind_proc_meri,ind_pts_meri
real :: smallval_north_top, smallval_south_top,smallval_north_bot, smallval_south_bot,r0
character(len=50) :: filename1
double precision, parameter :: epsi = 1d-30
logical :: use_meri,use_top,use_bot

! FROM POST PROCESSING

character(len=34), allocatable :: simdir(:)
character(len=4) :: appmynum,appidur,appmynum2
integer :: i,it,nsim,isim,iseis,ibeg,iend,nptstot,k1,k2,k3,j,npts_fluid,skipfact
real,dimension(:), allocatable :: colat,lon,time
real, dimension(:,:), allocatable :: seis,seis_fil
real, dimension(:,:,:), allocatable :: seis_sglcomp

real :: mij(6),conv_period,junk,srccolat,srclon,rec_loc_tol,rbot,rtop,dtheta,phi_incr,smallval
real, parameter :: pi = 3.1415926535898
character(len=32) :: conv_stf,stf_file,rot_rec_file,filename,fname
! simulation.info
character(len=34),allocatable, dimension(:) :: bkgrndmodel,stf_type
character(len=34),allocatable, dimension(:,:) :: src_type
character(len=1),dimension(3) :: reccomp
real,allocatable, dimension(:) :: dt,period,dt_seis,dt_strain,dt_snap,magnitude
logical,allocatable, dimension(:) :: correct_azi,rot_rec,rot_rec_post_array
integer,allocatable, dimension(:) :: nt,nrec,nt_seis,nt_strain,nt_snap
logical :: sum_seis_true,rot_src_rec_true,rot_rec_post

skipfact=1
smallval=10000.

! read snap plot parameters

open(unit=99,file='param_snaps')
read(99,*)phi0
write(6,*)'starting azimuth/phi for cross section on the right [deg]:',phi0
read(99,*)dphi
write(6,*)'ending azimuth/phi for cross section on the left [deg]:',dphi
read(99,*)rtop
write(6,*)'top surface [km]:',rtop
read(99,*)rbot
write(6,*)'bottom surface [km]:',rbot
read(99,*)theta_meri
write(6,*)'colatitude of meridional cross section:',theta_meri
read(99,*)snap1,snap2,snapskip
write(6,*)'1st,last snap, skipfactor:',snap1,snap2,snapskip
read(99,*)use_meri
write(6,*)'consider meridional cross section?',use_meri
read(99,*)use_top
write(6,*)'consider top surface?',use_top
read(99,*)use_bot
write(6,*)'consider bottom surface?',use_bot
read(99,*)nsim
if (nsim==1) write(6,*)'no need to sum, one simulation only!'
allocate(simdir(nsim))
do isim=1,nsim
   read(99,*)simdir(isim)
   write(6,*)isim,'simulation directory: ',trim(simdir(isim))
enddo
close(99)

phi0=phi0/180.*pi; dphi=dphi/180.*pi; rtop=rtop*1000.; rbot=rbot*1000.
theta_meri=theta_meri*pi/180.


allocate(bkgrndmodel(nsim),stf_type(nsim))
allocate(src_type(nsim,2))
allocate(dt(nsim),period(nsim),magnitude(nsim),dt_seis(nsim),dt_strain(nsim),dt_snap(nsim))
allocate(correct_azi(nsim),rot_rec(nsim),rot_rec_post_array(nsim))
allocate(nt(nsim),nrec(nsim),nt_seis(nsim),nt_strain(nsim),nt_snap(nsim))
do isim = 1,nsim
  open(unit=99,file=trim(simdir(isim))//'/simulation.info')
   read(99,*)bkgrndmodel(isim)
   read(99,*)dt(isim)
   read(99,*)nt(isim)
   read(99,*)src_type(isim,1)
   read(99,*)src_type(isim,2)
   read(99,*)stf_type(isim)
   read(99,*)period(isim)
   read(99,*)magnitude(isim)
   read(99,*)nrec(isim)
   read(99,*)nt_seis(isim)
   read(99,*)dt_seis(isim)
   read(99,*)correct_azi(isim)
   read(99,*)nt_strain(isim)
   read(99,*)dt_strain(isim)
   read(99,*)nt_snap(isim)
   read(99,*)dt_snap(isim)
   read(99,*)rot_rec(isim)
   read(99,*)ibeg
   read(99,*)iend
   close(99)
enddo

npts=nelem*(iend-ibeg+1)**2
nptstot=npts*nproc_mesh

write(6,*)'number of points per proc, total points:',npts,nptstot

! if all the same or not!!!!!
if (minval(dt)/=maxval(dt) .or. minval(nt)/=maxval(nt) .or. minval(period)/=maxval(period) .or. &
    minval(nrec)/=maxval(nrec)  .or. minval(nt_seis)/=maxval(nt_seis)  .or. &
    minval(dt_seis)/=maxval(dt_seis) .or. minval(nt_strain)/=maxval(nt_strain) .or. &
    minval(dt_strain)/=maxval(dt_strain) .or. minval(nt_snap)/=maxval(nt_snap) .or. &
    minval(dt_snap)/=maxval(dt_snap) ) then
   write(6,*)'PROBLEM with simulation.info parameters in the respective directories:'
   write(6,*)' one or more of the supposedly equal parameters differ!'
   stop
endif

! load and construct global mesh (one semi-disk)
write(6,*)'reading partitioned mesh...'
allocate(coord(nptstot,2))

smallval_north_top = rtop; smallval_south_top = rtop
smallval_north_bot = rbot; smallval_south_bot = rbot

do iproc=0,nproc_mesh-1
   call define_io_appendix(appmynum,iproc)
   open(unit=99,file=trim(simdir(1))//'/Data/glob_grid_'//appmynum//'.dat')
   do i=1,npts,skipfact
      read(99,*)coord(iproc*npts+i,1),coord(iproc*npts+i,2)

! determine minimal distance from rtop and rbot
      r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 
      if (coord(iproc*npts+i,2)>0.0 .and. abs(r0-rtop)< smallval_north_top) then ! north
         smallval_north_top = abs(r0-rtop)
      elseif (coord(iproc*npts+i,2)<0.0 .and. abs(r0-rtop)< smallval_south_top) then ! south
         smallval_south_top = abs(r0-rtop)
      endif
      if (coord(iproc*npts+i,2)>0.0 .and. abs(r0-rbot)< smallval_north_bot) then ! north            
         smallval_north_bot = abs(r0-rbot)
      elseif (coord(iproc*npts+i,2)<0.0 .and. abs(r0 -rbot)< smallval_south_bot) then ! south
         smallval_south_bot = abs(r0-rbot)
      endif
      
   enddo
   close(99)
enddo

smallval_north_top = smallval_north_top * 1.01
smallval_south_top = smallval_south_top * 1.01
smallval_north_bot = smallval_north_bot * 1.01
smallval_south_bot = smallval_south_bot * 1.01

write(6,*)'Smallest distance to rtop (North,South) [km]:', &
     real(smallval_north_top/1000.),real(smallval_south_top/1000.)
write(6,*)'Smallest distance to rbot (North,South) [km]:', &
     real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)

if (use_top .or. use_meri) allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)),ind_pts_top_tmp(floor(real(nptstot)/10.)))
if (use_bot .or. use_meri) allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)),ind_pts_bot_tmp(floor(real(nptstot)/10.)))

k1=0; k2=0; 

do iproc=0,nproc_mesh-1
   do i=1,npts,skipfact
      ! check for top and bottom radii
      r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 
      if (use_top .or. use_meri) then
      if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
              abs(r0-rtop) <= smallval_north_top) .or. &
           (coord(iproc*npts+i,2)<0.0 .and.  &
              abs(r0-rtop) <= smallval_south_top)) then 
         k1=k1+1         
         ind_proc_top_tmp(k1) = iproc
         ind_pts_top_tmp(k1) = i
       endif
       endif
       if (use_bot .or. use_meri) then 
      if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
              abs(r0-rbot) <= smallval_north_bot) .or. &
           (coord(iproc*npts+i,2)<0.0 .and.  &
              abs(r0-rbot) <= smallval_south_bot)) then 
         k2=k2+1
         ind_proc_bot_tmp(k2) = iproc
         ind_pts_bot_tmp(k2) = i
       endif
       endif
    enddo
enddo

npts_top=k1
npts_bot=k2

write(6,*)'# points on top,bottom:',npts_top,npts_bot

write(6,*)'allocating index arrays for surfaces....'
if (use_top .or. use_meri) then
   allocate(ind_proc_top(npts_top),ind_pts_top(npts_top))
   ind_proc_top=ind_proc_top_tmp(1:npts_top); ind_pts_top=ind_pts_top_tmp(1:npts_top)
   deallocate(ind_proc_top_tmp,ind_pts_top_tmp)
endif
if (use_bot .or. use_meri) then
   allocate(ind_proc_bot(npts_bot),ind_pts_bot(npts_bot))
   ind_proc_bot=ind_proc_bot_tmp(1:npts_bot); ind_pts_bot=ind_pts_bot_tmp(1:npts_bot)
   deallocate(ind_proc_bot_tmp,ind_pts_bot_tmp)
endif


! MERIDIONAL based on rtop and rbottom-----------------------------------------------------------------------------
if (use_meri) then

write(6,*)'computing meridional preparameters....'

! find closest theta at rbot
smallval_meri = 2*pi
do i=1,npts_bot
   r0 = sqrt(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)**2 + coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)**2)
   theta0 = atan(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)/coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)+epsi)
   if ( theta0 <0. ) theta0 = pi + theta0
   if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
!   write(6,*)'r theta bot:',r0/1000.,theta0*180./pi
   if (abs(theta_meri-theta0) < smallval_meri) then 
      smallval_meri = abs(theta_meri-theta0) 
      theta_bot = theta0 
   endif
enddo
write(6,*)'theta meri,theta closest at rbot:',theta_meri*180./pi,theta_bot*180./pi
theta_meri = theta_bot

! find theta at rtop closest to theta from rbot
smallval_meri = 2*pi
do i=1,npts_top
   r0 = sqrt(coord(ind_proc_top(i)*npts+ind_pts_top(i),1)**2 + coord(ind_proc_top(i)*npts+ind_pts_top(i),2)**2)
   theta0 = atan(coord(ind_proc_top(i)*npts+ind_pts_top(i),1)/coord(ind_proc_top(i)*npts+ind_pts_top(i),2)+epsi)
   if ( theta0 <0. ) theta0 = pi + theta0
   if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
!   write(6,*)'r theta top:',r0/1000.,theta0*180./pi
   if (abs(theta_bot-theta0) < smallval_meri) then 
      smallval_meri = abs(theta_bot-theta0) 
      theta_top = theta0 
   endif
enddo

smallval_meri=abs(theta_top-theta_bot)
write(6,*)'theta closest at rtop and smallval:',theta_top*180./pi,smallval_meri*180./pi
if (theta_top > theta_bot) then 
   theta_meri = theta_bot+ smallval_meri/2.
elseif (theta_top < theta_bot) then
   theta_meri = theta_bot-smallval_meri/2.
else
   theta_meri = theta_bot
endif
smallval_meri=smallval_meri*2.5

k3=0
allocate(ind_proc_meri_tmp(floor(real(nptstot)/10.)),ind_pts_meri_tmp(floor(real(nptstot)/10.)))

do iproc=0,nproc_mesh-1
   do i=1,npts,skipfact
      r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 
        theta0=atan(coord(iproc*npts+i,1)/(coord(iproc*npts+i,2)+epsi))
       if ( theta0 <0. ) theta0 = pi + theta0
       if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
       if ( r0>=rbot .and. abs(theta0-theta_meri) <= smallval_meri ) then
          k3=k3+1
          ind_proc_meri_tmp(k3) = iproc
          ind_pts_meri_tmp(k3) = i
          write(62,*)r0,theta0*180./pi
       endif
enddo
enddo
npts_meri=k3

write(6,*)'# points on meridional:',npts_meri
allocate(ind_proc_meri(npts_meri),ind_pts_meri(npts_meri))
ind_proc_meri=ind_proc_meri_tmp(1:npts_meri); ind_pts_meri=ind_pts_meri_tmp(1:npts_meri)
deallocate(ind_proc_meri_tmp,ind_pts_meri_tmp)
endif ! use_meri

! xyz coordinates---------------------------------------------------------------------------------------------
write(6,*)'defining xyz...'
allocate(x(1:2*nptstot),y(1:2*nptstot),z(1:2*nptstot));x=0.;y=0.;z=0.

! left cross section---------------------------------------------------------------------------------------------
call sphi2xy(x(1:nptstot),y(1:nptstot),coord(:,1),phi0,nptstot)
z(1:nptstot)=coord(1:nptstot,2)
write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

! right cross section---------------------------------------------------------------------------------------------
call sphi2xy(x(nptstot+1:2*nptstot),y(nptstot+1:2*nptstot),coord(:,1),dphi,nptstot)
z(nptstot+1:2*nptstot)=coord(1:nptstot,2)

write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

! save xyz
allocate(vp(2*nptstot))
!write(6,*)'saving xyz to file Data/xyz_xsections.dat ...'
!open(unit=99,file='Data/xyz_xsections.dat')
do i=1,2*nptstot  
!   write(99,*)x(i),y(i),z(i)
   r= sqrt(x(i)**2+y(i)**2+z(i)**2)
   vp(i)=prem(r,'v_p')
enddo
!close(99)

! rescale vp
vp=vp/maxval(vp)*0.002-0.001

filename1='mesh_xsect'
call write_VTK_bin_scal(x,y,z,vp,2*nptstot,0,filename1)


! top surface---------------------------------------------------------------------------------------------
k1=0
if (use_top) then
write(6,*)'defining top surface...'
dtheta=rtop*pi/npts_top
write(6,*)'number of top surface points in phi cakepiece:',npts_top,180.*dtheta/pi

allocate(xtop(1:npts_top**2),ytop(1:npts_top**2),ztop(1:npts_top**2));xtop=0.;ytop=0.;ztop=0.
allocate(azi_ind_top(npts_top**2),azi_phi_top(npts_top**2))

   do i=1,npts_top
      nphi = max(floor(coord(ind_proc_top(i)*npts+ind_pts_top(i),1)*(2*pi-(dphi-phi0))/dtheta/2.),1)
      phi_incr = (2.*pi-(dphi-phi0))/nphi
!      write(61,*)'i,nphi,phi_incr',i,npts_top,nphi,phi_incr
      do j=1,nphi
         k1=k1+1
         xtop(k1) = coord(ind_proc_top(i)*npts+ind_pts_top(i),1)*cos(phi0-(j-1)*phi_incr)
         ytop(k1) = coord(ind_proc_top(i)*npts+ind_pts_top(i),1)*sin(phi0-(j-1)*phi_incr)
         ztop(k1) = coord(ind_proc_top(i)*npts+ind_pts_top(i),2)
         azi_ind_top(k1) = ind_proc_top(i)*npts+ind_pts_top(i)
         azi_phi_top(k1) = phi0-(j-1)*phi_incr
!         write(67,*)ztop(k1),coord(ind_proc_top(i)*npts+ind_pts_top(i),2)
   enddo
enddo

open(unit=99,file=trim(simdir(1))//'/Data/xyz_top.dat')
do i=1,k1
   write(99,*)xtop(i),ytop(i),ztop(i)
enddo
close(99)

! extract vp
allocate(vptop(k1))
vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))

! save into AVS
filename1='mesh_top'
call write_VTK_bin_scal(xtop,ytop,ztop,vptop,k1,0,filename1)
endif 

! bottom surface ---------------------------------------------------------------------------------------------
k2=0
if (use_bot) then 
write(6,*)'defining bottom surface...'
dtheta=rbot*pi/npts_bot
write(6,*)'# pts on rbot phi cakepiece and average spacing [m]:',npts_bot,dtheta*180./pi
allocate(xbot(1:npts_bot**2),ybot(1:npts_bot**2),zbot(1:npts_bot**2));xbot=0.;ybot=0.;zbot=0.
allocate(azi_ind_bot(npts_bot**2),azi_phi_bot(npts_bot**2))

   do i=1,npts_bot
      nphi=max(floor(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)*(dphi-phi0)/dtheta/1.),1)
      phi_incr = (dphi-phi0)/nphi
     write(62,*)'i,nphi,phi_incr',i,npts_bot,nphi,phi_incr
      do j=1,nphi
         k2=k2+1
         xbot(k2) = coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)*cos(phi0+(j-1)*phi_incr)
         ybot(k2) = coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)*sin(phi0+(j-1)*phi_incr)
         zbot(k2) = coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)
         azi_ind_bot(k2) = ind_proc_bot(i)*npts+ind_pts_bot(i)
         azi_phi_bot(k2) = phi0+(j-1)*phi_incr
   enddo
enddo

open(unit=99,file=trim(simdir(1))//'/Data/xyz_bot.dat')
do i=1,k2
   write(99,*)xbot(i),ybot(i),zbot(i)
enddo
close(99)

! extract vp
allocate(vpbot(k2))
!vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))
vpbot(1:k2) = minval(vp)*0.9

! save into AVS
write(6,*)'writing vtk bin file...'
filename1='mesh_bot'
!call write_avs_file_scal(xbot,ybot,zbot,vpbot,k2,0,filename1)
call write_VTK_bin_scal(xbot,ybot,zbot,vpbot,k2,0,filename1)

endif


! meridional cross section ---------------------------------------------------------------------------------------------
k3=0
if (use_meri) then 
write(6,*)'defining meridional cross section...'
dr=(6371000.-rbot)/npts_meri
write(6,*)'# pts on rmeri and average spacing [km]:',npts_meri,dr/1000.
allocate(xmeri(1:7*npts_meri**2),ymeri(1:7*npts_meri**2),zmeri(1:7*npts_meri**2));xmeri=0.;ymeri=0.;zmeri=0.
allocate(azi_ind_meri(7*npts_meri**2),azi_phi_meri(7*npts_meri**2))
allocate(vpmeri(1:7*npts_meri**2))

   do i=1,npts_meri
      r0 = sqrt( coord(ind_proc_meri(i)*npts+ind_pts_meri(i),1)**2 + coord(ind_proc_meri(i)*npts+ind_pts_meri(i),2)**2)
      nphi=max(floor(r0*pi/dr/1.),1)
      phi_incr = pi/nphi
      write(6,*)i,'r0,nphi,phi_incr [km]:',r0/1000.,nphi,phi_incr*r0/1000.
      do j=1,nphi
         k3=k3+1
         phi=(j-1)*phi_incr
         call rthetaphi2xyz(xmeri(k3),ymeri(k3),zmeri(k3),r0,theta_meri,phi)
         azi_ind_meri(k3) = ind_proc_meri(i)*npts+ind_pts_meri(i)
         azi_phi_meri(k3) = phi
         vpmeri(k3) = vp(ind_proc_meri(i)*npts+ind_pts_meri(i))
   enddo
enddo

! save into AVS
write(6,*)'writing vtk bin file...'
filename1='mesh_meri'
call write_VTK_bin_scal(xmeri,ymeri,zmeri,vpmeri,k3,0,filename1)

endif !use_meri

! assembling everything to one coordinate array-----------------------------------------------------------------
deallocate(coord)
deallocate(vp)

allocate(xtot(2*nptstot+k1+k2+k3),ytot(2*nptstot+k1+k2+k3),ztot(2*nptstot+k1+k2+k3))

xtot(1:2*nptstot)=x; ytot(1:2*nptstot)=y; ztot(1:2*nptstot)=z
if (use_top) then 
xtot(2*nptstot+1:2*nptstot+k1)=xtop
ytot(2*nptstot+1:2*nptstot+k1)=ytop
ztot(2*nptstot+1:2*nptstot+k1)=ztop
endif
if (use_bot) then 
xtot(2*nptstot+k1+1:2*nptstot+k1+k2)=xbot
ytot(2*nptstot+k1+1:2*nptstot+k1+k2)=ybot
ztot(2*nptstot+k1+1:2*nptstot+k1+k2)=zbot
endif
if (use_meri) then 
   xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=xmeri
   ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=ymeri
   ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=zmeri
endif

npts_fluid=nel_fluid*(iend-ibeg+1)**2
write(6,*)'points in fluid:',npts_fluid

! save mesh for kerner
open(unit=99,file='mesh_tot.xyz')
do i=1,2*nptstot+k1+k2
   write(99,*)xtot(i),ytot(i),ztot(i)
enddo
close(99)

if (use_meri) then 
open(unit=99,file='mesh_meri.xyz')
do i=2*nptstot+k1+k2+1,2*nptstot+k1+k2+k3
   write(99,*)xtot(i),ytot(i),ztot(i)
enddo
close(99)
endif

allocate(disp(2*nptstot+k1+k2+k3,3))
allocate(disptot(2*nptstot+k1+k2+k3,3))
allocate(fluid_prefact(2*npts_fluid,2,nsim))
allocate(azi(k1+k2,2,nsim),azi_xsect2(2,nsim),azi_xsect1(2,nsim),azi_meri(k3,2,nsim))

! load azimuthal prefactors--------------------------------------------------------------------------------------------
do isim = 1,nsim

  open(unit=99,file=trim(simdir(isim))//'/simulation.info')
   read(99,*)bkgrndmodel(isim)
   read(99,*)dt(isim)
   read(99,*)nt(isim)
   read(99,*)src_type(isim,1)
   read(99,*)src_type(isim,2)
   read(99,*)stf_type(isim)
   read(99,*)period(isim)
   read(99,*)magnitude(isim)
   read(99,*)nrec(isim)
   read(99,*)nt_seis(isim)
   read(99,*)dt_seis(isim)
   read(99,*)correct_azi(isim)
   read(99,*)nt_strain(isim)
   read(99,*)dt_strain(isim)
   read(99,*)nt_snap(isim)
   read(99,*)dt_snap(isim)
   read(99,*)rot_rec(isim)
   read(99,*)ibeg
   read(99,*)iend
   close(99)
!

! load fluid prefactors---------------------------------------------------------------------------------------------

write(6,*)'loading fluid prefactors...'
do iproc=0,nproc_mesh-1
   open(unit=190,file=trim(simdir(isim))//'/Data/inv_rho_s_fluid_globsnaps_'//appmynum//'.dat')
   do i=1,npts_fluid,skipfact
      read(190,*)fluid_prefact(iproc*npts_fluid+i,1,isim),fluid_prefact(iproc*npts_fluid+i,2,isim)
   enddo
   close(190)
enddo

   if (src_type(isim,1)=='monopole') then 
      azi_xsect1(1:2,isim) = 1.
      azi_xsect2(1:2,isim) = 1.
      azi(:,1:2,isim) = 1.
      azi_meri(:,1:2,isim) = 1.
   else
      if (src_type(isim,2)=='mxz' .or. src_type(isim,2)=='xforce') then 
         azi_xsect1(1,isim) = cos(phi0);  azi_xsect1(2,isim) = -sin(phi0)
         azi_xsect2(1:2,isim) = cos(phi0+dphi);  azi_xsect1(2,isim) = -sin(phi0+dphi)
         do i=1,k1
            azi(i,1,isim) = cos(azi_phi_top(i))
            azi(i,2,isim) = -sin(azi_phi_top(i))
         enddo
         do i=1,k2
            azi(k1+i,1,isim) = cos(azi_phi_bot(i))
            azi(k1+i,2,isim) = -sin(azi_phi_bot(i))
         enddo
          do i=1,k3
            azi_meri(i,1,isim) = cos(azi_phi_meri(i))
            azi_meri(i,2,isim) = -sin(azi_phi_meri(i))
         enddo
           
      elseif (src_type(isim,2)=='myz' .or. src_type(isim,2)=='yforce') then 
         azi_xsect1(1,isim) = sin(phi0);  azi_xsect1(2,isim) = cos(phi0)
         azi_xsect2(1:2,isim) = sin(phi0+dphi);  azi_xsect1(2,isim) = cos(phi0+dphi)
         do i=1,k1
            azi(i,1,isim) = sin(azi_phi_top(i))
            azi(i,2,isim) = cos(azi_phi_top(i))
         enddo
         do i=1,k2
            azi(k1+i,1,isim) = sin(azi_phi_bot(i))
            azi(k1+i,2,isim) = cos(azi_phi_bot(i))
         enddo
         do i=1,k3
            azi_meri(i,1,isim) = sin(azi_phi_meri(i))
            azi_meri(i,2,isim) = cos(azi_phi_meri(i))
         enddo

      elseif (src_type(isim,2)=='mxy' ) then 
            azi_xsect1(1,isim) = cos(2.*phi0);  azi_xsect1(2,isim) = -sin(2.*phi0)
            azi_xsect2(1:2,isim) = cos(2.*(phi0+dphi));  azi_xsect1(2,isim) = -sin(2.*(phi0+dphi))
            do i=1,k1
               azi(i,1,isim) = cos(2.*azi_phi_top(i))
               azi(i,2,isim) = -sin(2.*azi_phi_top(i))
            enddo
            do i=1,k2
               azi(k1+i,1,isim) = cos(2.*azi_phi_bot(i))
               azi(k1+i,2,isim) = -sin(2.*azi_phi_bot(i))
            enddo
            do i=1,k3
               azi_meri(i,1,isim) = cos(2.*azi_phi_meri(i))
             enddo
         else
            write(6,*)'canot determine source type!',src_type(isim,2)
         endif
      endif

      write(6,*)'azimuthal factors for phi0 xsections:',simdir(isim),azi_xsect1(1,isim),azi_xsect1(2,isim)
      write(6,*)'azimuthal factors for phi0+dphi xsections:',simdir(isim),azi_xsect2(1,isim),azi_xsect2(2,isim)


enddo !isim

! load snaps ===============================================================
write(6,*)'loading snaps...'

do j=snap1,snap2,snapskip

   disp=0.
   disptot=0.

   do isim=1,nsim

   do iproc=0,nproc_mesh-1
      call define_io_appendix(appmynum,iproc)
      call define_io_appendix(appmynum2,j)      
      open(unit=99,file=trim(simdir(isim))//'/Data/snap_'//appmynum//'_'//appmynum2//'.dat')
      do i=1,npts,skipfact
         read(99,*)disp(iproc*npts+i,1),disp(iproc*npts+i,2),disp(iproc*npts+i,3)
      enddo
      close(99)
      disp(iproc*npts+1:iproc*npts+npts_fluid,1)= &
           disp(iproc*npts+1:iproc*npts+npts_fluid,1) * fluid_prefact(:,1,isim)
      disp(iproc*npts+1:iproc*npts+npts_fluid,2)= &
           disp(iproc*npts+1:iproc*npts+npts_fluid,2) * fluid_prefact(:,1,isim) * fluid_prefact(:,2,isim)
      disp(iproc*npts+1:iproc*npts+npts_fluid,3)= &
           disp(iproc*npts+1:iproc*npts+npts_fluid,3) * fluid_prefact(:,1,isim)

   enddo
   disp(nptstot+1:2*nptstot,1)=disp(1:nptstot,1)*azi_xsect2(1,isim)
   disp(nptstot+1:2*nptstot,2)=disp(1:nptstot,2)*azi_xsect2(2,isim)
   disp(nptstot+1:2*nptstot,3)=disp(1:nptstot,3)*azi_xsect2(1,isim)

   do i=1,k1
      disp(2*nptstot+i,1) = disp(azi_ind_top(i),1)*azi(i,1,isim)
      disp(2*nptstot+i,2) = disp(azi_ind_top(i),2)*azi(i,2,isim)
      disp(2*nptstot+i,3) = disp(azi_ind_top(i),3)*azi(i,1,isim)
   enddo
   do i=1,k2
      disp(2*nptstot+k1+i,1) = disp(azi_ind_bot(i),1)*azi(k1+i,1,isim)      
      disp(2*nptstot+k1+i,2) = disp(azi_ind_bot(i),2)*azi(k1+i,2,isim)
      disp(2*nptstot+k1+i,3) = disp(azi_ind_bot(i),3)*azi(k1+i,1,isim)
   enddo
   if (use_meri) then 
      do i=1,k3
         disp(2*nptstot+k1+k2+i,1) = disp(azi_ind_meri(i),1)*azi_meri(i,1,isim)      
         disp(2*nptstot+k1+k2+i,2) = disp(azi_ind_meri(i),2)*azi_meri(i,2,isim)
         disp(2*nptstot+k1+k2+i,3) = disp(azi_ind_meri(i),3)*azi_meri(i,1,isim)
      enddo
   endif

   disp(1:nptstot,1)=disp(1:nptstot,1)*azi_xsect1(1,isim)
   disp(1:nptstot,2)=disp(1:nptstot,2)*azi_xsect1(2,isim)
   disp(1:nptstot,3)=disp(1:nptstot,3)*azi_xsect1(1,isim)

   disptot = disp+ disptot

enddo
filename1='snap_'//appmynum2//'_z'
call write_VTK_bin_scal(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2),ztot(1:2*nptstot+k1+k2), & 
                                           disptot(1:2*nptstot+k1+k2,3),2*nptstot+k1+k2,0,filename1)

if (use_meri) then
   filename1='meri_snap_'//appmynum2//'_z'
   call write_VTK_bin_scal(xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3),ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                                          ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3),disptot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3,3),k3,0,filename1)
endif


enddo

end program load_snaps

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  subroutine define_io_appendix(app,iproc)
!
! Defines the 4 digit character string appended to any 
! data or io file related to process myid. 
!

integer :: iproc
character(len=4) :: app
character(len=1) :: milp,cenp,dizp,unip

  milp = char(48+    iproc/1000)
  cenp = char(48+mod(iproc/100,10))
  dizp = char(48+mod(iproc/10,10))
  unip = char(48+mod(iproc,10))

  app = milp//cenp//dizp//unip

end subroutine define_io_appendix
!=============================================================================

!dk szphi2xyz----------------------------------------------------------
  subroutine sphi2xy(x,y,s,phi,n)
    integer, intent(in) :: n
  real, dimension(1:n), intent(out) :: x,y
  real, dimension(1:n), intent(in) :: s
  real, intent(in) :: phi
 
  x(1:n)=s(1:n)*cos(phi)
  y(1:n)=s(1:n)*sin(phi)

  end subroutine sphi2xy
!--------------------------------------------------------------------------

subroutine write_avs_file(x,y,z,u1,u2,u3,rows,nelem_disk,filename1)

 implicit none
 
 integer :: i,t,rows,ielem,ipol,jpol,nelem_disk,ioerr,dims
real :: Zcoord,phi
real, dimension(1:rows), intent(in) :: x,y,z,u1,u2,u3
 character (len=30) :: celltype;
 character (len=20) :: filename1, file2

 write(6,*)'computing avs file ',trim(filename1),' ...'

 open(unit=2502,file=trim(filename1),status='replace',form='formatted')

 celltype='pt'
 !avs file headers
 write(2502,*) rows,rows, 1, 0, 0;
 write(6,*) rows
 ! write coordinates and numbering
 dims=SIZE(x)
 write(6,*) dims,size(y),size(z)
 do i=1,rows    
    write(2502,'(i10,3f14.2)') i, x(i),y(i),z(i)
 enddo   
   
 ! connectivity structure
 do i=1,rows
    write(2502,'(i12,i9,a4,i12)') i,1,trim(celltype),i
 enddo                   
 !property headers
    write(2502,*) 3,1,1,1;
     write(2502,*) 'u1,', 'm';
    write(2502,*) 'u2,', 'm';
    write(2502,*) 'u3,', 'm';
 !displacement fields

 do i=1,rows
    write(2502,'(i12,3f14.2)') i,u1(i),u2(i),u3(i);
 enddo               

 close(2502)

write(6,*)'...saved ',trim(filename1)

end subroutine write_avs_file

!-----------------------------------------------------------------------------



subroutine write_avs_file_scal(x,y,z,u1,rows,nelem_disk,filename1)

 implicit none
 
 integer :: i,t,rows,ielem,ipol,jpol,nelem_disk,ioerr,dims
real :: Zcoord,phi
real, dimension(1:rows), intent(in) :: x,y,z,u1
 character (len=30) :: celltype
 character (len=20) :: filename1;

 write(6,*)'computing avs file ',trim(filename1),' ...'

 open(unit=2502,file=trim(filename1),status='replace',form='formatted')

 celltype='pt'
 !avs file headers
 write(2502,*) rows,rows, 1, 0, 0;
 write(6,*) rows
 ! write coordinates and numbering
 dims=SIZE(x)
 write(6,*) dims,size(y),size(z)
 do i=1,rows    
    write(2502,'(i10,3f14.2)') i, x(i),y(i),z(i)
 enddo   
  
 ! connectivity structure
 do i=1,rows
    write(2502,'(i12,i9,a4,i12)') i,1,trim(celltype),i
 enddo                   
 !property headers
    write(2502,*) 1,1;
     write(2502,*) 'u1,', 'm';
 !displacement fields

 do i=1,rows
    write(2502,'(i12,e14.3)') i,u1(i)
 enddo               

 close(2502)

write(6,*)'...saved ',trim(filename1)

end subroutine write_avs_file_scal

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal(x,y,z,u1,rows,nelem_disk,filename1)

 implicit none
 integer :: i,t,rows,nelem_disk,ioerr,dims
 real, dimension(1:rows), intent(in) :: x,y,z,u1
 integer, dimension(1:rows,2) :: cell
 integer, dimension(1:rows) :: cell_type
 real, dimension(1:rows,3) :: W
 
 character (len=30) :: celltype;
 character (len=30) :: filename1;
 character (len=50) :: ss; !stream
 
 W(1:rows,1)=x
 W(1:rows,2)=y
 W(1:rows,3)=z
!points structure
do i=1,rows
 cell(i,1)=1
 cell(i,2)=i
 cell_type(i)=1
enddo

 write(6,*)'computing VTK bin file ',trim(filename1)//'.vtk  ...'

open(100,file=trim(filename1)//'.vtk',access='stream',status='replace',convert='big_endian')

write(100) '# vtk DataFile Version 3.0'//char(10)
write(100) 'Cell Fractions'//char(10)
write(100) 'BINARY'//char(10)
write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
write(ss,fmt='(A8,I12,A10)') 'POINTS',rows,' float'
write(100) ss//char(10)
!points
do i=1,rows
write(100) W(i,1:3)
enddo
write(100) char(10)
!cell topology
write(ss,fmt='(A5,2I12)') 'CELLS',rows,rows*2
write(100) char(10)//ss//char(10)
do i=1,rows
write(100) cell(i,1:2)
enddo
write(100) char(10)
!cell type
write(ss,fmt='(A10,2I12)') 'CELL_TYPES',rows
write(100) char(10)//ss//char(10)
do i=1,rows
write(100) cell_type(i)
enddo
write(100) char(10)
!data
write(ss,fmt='(A10,I12)') 'CELL_DATA',rows
write(100) char(10)//ss//char(10)
write(100) 'SCALARS Displ_u1 float 1'//char(10)
write(100) 'LOOKUP_TABLE default'//char(10) !color table?
do i=1,rows
write(100) u1(i)
enddo
 close(100)
write(6,*)'...saved ',trim(filename1)//'.vtk'
end subroutine write_vtk_bin_scal
!-----------------------------------------------------------------------------


!dk rthetaphi2xyz----------------------------------------------------------
  subroutine rthetaphi2xyz(x,y,z,r,theta,phi)
  real, intent(out) :: x,y,z
  real, intent(in) :: r,theta,phi
  x=r*sin(theta)*cos(phi)
  y=r*sin(theta)*sin(phi)
  z=r*cos(theta)
  end subroutine rthetaphi2xyz
!--------------------------------------------------------------------------

real function prem(r0,param)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

real, intent(in) :: r0
real            :: r,x_prem
real             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(r>6356. .and. r .le. 6371.01)THEN        ! upper crustal layer
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  ELSEIF(r>6346.6 .and. r .le. 6356.)THEN
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  ELSEIF(r>6151. .and. r .le. 6346.6)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(r>5971. .and. r .le. 6151. )THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(r>5771. .and. r .le. 5971.)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(r>5701. .and. r .le. 5771. )THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(r>5600. .and. r .le. 5701. )THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(r>3630. .and. r .le. 5600. )THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(r>3480. .and. r .le. 3630.)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(r>1221.5 .and. r .le. 3480. )THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  ELSEIF(r .le. 1221.5)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ELSE 
     write(6,*)'wrong radius!',r; stop
  ENDIF

  if (param=='rho') then
     prem=ro_prem*1000.
  elseif (param=='v_p') then
     prem=vp_prem*1000.
  elseif (param=='v_s') then
     prem=vs_prem*1000.
  else
     write(6,*)'ERROR IN PREM FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem
!=============================================================================
