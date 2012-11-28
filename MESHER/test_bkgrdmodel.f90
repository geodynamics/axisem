!=============================
  module test_bkgrdmodel
!=============================

  use data_gllmesh
  use data_mesh
  use data_spec
  use data_grid
  use data_diag
  use data_bkgrdmodel
  use global_parameters

  implicit none

  public :: bkgrdmodel_testing, write_VTK_bin_scal_old, write_VTK_bin_scal
  private 

  contains

!-----------------------------------------------------------------
subroutine bkgrdmodel_testing

use background_models

double precision, dimension(:,:,:), allocatable :: h,hmin2
double precision, dimension(:,:), allocatable :: crit,crit_max
double precision, dimension(:), allocatable :: hmin,hmax
integer :: iel, ipol, jpol,ntoobig,ntoosmall,j
double precision :: s1,z1,s2,z2,h1,h2,r,velo,velo_max,theta
double precision,allocatable, dimension(:,:,:) :: v_p,v_s,rho

! vtk
real, dimension(:,:), allocatable :: mesh2
real, dimension(:), allocatable :: vp1,vs1,h_real,rho1
character(len=80) :: fname
integer :: npts_vtk,ct
real, allocatable ::  x(:),y(:),z(:)

allocate(crit(0:npol,0:npol)) ; crit(:,:) = 0.d0 
allocate(crit_max(0:npol,0:npol)) ; crit_max(:,:) = 0.d0

allocate(h(0:npol,0:npol,neltot)) ;   h(0:npol,0:npol,neltot) = 0.d0
allocate(hmin2(0:npol,0:npol,neltot)) ;   hmin2(:,:,:) = 0.d0

allocate(hmin(neltot),hmax(neltot)) ;   hmin(:) = 0.d0; hmax(:)=0.d0

ntoobig=0
ntoosmall=0
j=0

if (dump_mesh_info_screen) write(6,*)''

! construct full arrays for velocities... faster! 
if (bkgrdmodel=='solar') then 
   write(6,*)'pre-assembling media arrays for solar case.... '
   write(6,*)'...much faster due to inperpolation routine!'
   allocate(v_p(0:npol,0:npol,neltot),v_s(0:npol,0:npol,neltot),rho(0:npol,0:npol,neltot))
   write(6,*)'allocated global media arrays'
   call arbitr_sub_solar_arr(sgll*router,zgll*router,v_p,v_s,rho,bkgrdmodel)
   write(6,*)'done with media array definition '
   write(6,*)'minmax vp:',minval(v_p),maxval(v_p),minval(rho),maxval(rho)
endif

! find smallest/largest grid spacing
do iel = 1, neltot

!   if (mod(iel,floor(real(neltot)/100.))==0) write(6,*)'iel, %',iel,real(iel)/real(neltot)*100.

   do jpol = 0, npol-1
      do ipol = 0,npol-1
         s1 = sgll(ipol,jpol,iel)
         z1 = zgll(ipol,jpol,iel) 
         s2 = sgll(ipol+1,jpol,iel)
         z2 = zgll(ipol+1,jpol,iel)
         !   write(669,*)'1: ',ipol,jpol,s2-s1,z2-z1
         !   write(667,*)'1c ',s1,s2,z1,z2
         !   write(665,*)sgll(ipol,jpol,iel),zgll(ipol,jpol,iel)

         !   if (iel==10) write(668,*)sgll(ipol,jpol,iel),zgll(ipol,jpol,iel)

         h1 = dsqrt((s2-s1)**2+(z2-z1)**2)

         s2 = sgll(ipol,jpol+1,iel)
         z2 = zgll(ipol,jpol+1,iel)
         !   write(669,*)'2: ',ipol,jpol,s2-s1,z2-z1
         !   write(667,*)'2c ',s1,s2,z1,z2

         h2 = dsqrt((s2-s1)**2+(z2-z1)**2)
         hmin2(ipol,jpol,iel) = min(h1,h2)
         hmin2(ipol,jpol,iel) = router * hmin2(ipol,jpol,iel)  
         h(ipol,jpol,iel) = max(h1,h2)
         !   write(666,*)h1,h2
         h(ipol,jpol,iel) = router * h(ipol,jpol,iel)  

      end do
   end do

   ! check on element edges
   do ipol = npol, npol
      do jpol = 0, npol -1 
         h(ipol,jpol,iel) = h(ipol-1,jpol,iel)
         hmin2(ipol,jpol,iel) = hmin2(ipol-1,jpol,iel)
      end do
   end do

   do jpol = npol, npol 
      do ipol = 0, npol-1
         h(ipol,jpol,iel) = h(ipol,jpol-1,iel)
         hmin2(ipol,jpol,iel) = hmin2(ipol,jpol-1,iel)
      end do
   end do
   h(npol,npol,iel) = h(npol-1,npol,iel)
   hmin2(npol,npol,iel) = hmin2(npol-1,npol,iel)
   hmin(iel)=minval(hmin2(:,:,iel))
   hmax(iel)=maxval(h(:,:,iel))

end do ! elements

write(6,*)'calculate GLL spacing...'

! global min/max spacing
hmin_glob=minval(hmin)
hmax_glob=maxval(hmax)

! min. distance in global domain, e.g. to find identical points on boundaries
min_distance_dim=hmin_glob*0.1d0
min_distance_nondim=hmin_glob*0.1d0/router

if (dump_mesh_info_screen) then  
   write(6,*)'Minimal spacing in global domain [m]:',hmin_glob
   write(6,*)'Maximal spacing in global domain [m]:',hmax_glob
   write(6,*)'Minimal distances in global domain [m]:', &
        min_distance_dim
   write(6,*)'Minimal distances in global domain [non-dim]:', &
        min_distance_nondim; call flush(6)
end if
if (dump_mesh_info_files) then 
   open(unit=61,file=diagpath(1:lfdiag)//'/gridspacing_toolarge_small.dat')
end if

open(unit=62,file=diagpath(1:lfdiag)//'/radial_velocity.dat')  

! vtk preparations
npts_vtk = neltot * 4

if (dump_mesh_vtk) then
    allocate(mesh2(neltot,2), h_real(neltot))
    allocate(vp1(npts_vtk), vs1(npts_vtk), rho1(npts_vtk))
    allocate(x(npts_vtk), y(npts_vtk), z(npts_vtk))
    z = 0.d0
endif

ct = 0

write(6,*)'starting big loop....'
do iel = 1, neltot
!   if (mod(iel,floor(real(neltot)/100.))==0) write(6,*)'iel, %',iel,real(iel)/real(neltot)*100.
   do jpol = 0, npol
      do ipol = 0,npol
         s1 = sgll(ipol,jpol,iel)
         z1 = zgll(ipol,jpol,iel) 

         if (s1<1.d-5 .and. z1>=0.d0) then 
            if (bkgrdmodel=='solar') then 
                  velo = v_p(ipol,jpol,iel)

            else 
               r = dsqrt(s1**2+z1**2)
!               r=dble(int(r*1.d10))*1.d-10 ! TNM: ADDED BACK IN JAN 2011... NOT SURE IF NECESSARY !!!!!

               if ( solid_domain(region(iel))) then 
                  velo = velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)
               else
                  velo = velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
               endif
               write(62,*)r*router,velo; call flush(62)
            endif
         endif

         ! write element boundary coordinates of central region
         !  if (r*router<=rmin ) then 
         !     if (jpol==0 .or. jpol==npol) write(343,*)s1,z1
         !     if (ipol==0 .or. ipol==npol) write(343,*)s1,z1
         !  endif

         ! DANGEROUS STUFF HAPPENING RIGHT HERE!
         ! at ICB, some values are 1221.4999999, some 1221.5 ....
         ! hence we cap the accuracy of the radius

         !   velo = velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)
         !   if (velo < 2.d3 ) velo = velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
         !   if (.not. resolve_inner_shear .and. region(iel)==ndisc) & 
         !        velo = velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)

!         if (dump_mesh_info_files ) then ! .and. ipol==2 .and. jpol==2) then !!! should be for ALL ipol/jpol!!! Just slow for the sun...

            if (bkgrdmodel=='solar') then 
               velo = v_p(ipol,jpol,iel)
               velo_max = v_p(ipol,jpol,iel)

            else
             r = dint(r*1.d10)*1.d-10
            
               if ( solid_domain(region(iel))) then 
                  velo = velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)
               else
                  velo = velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
               endif
               velo_max = velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
            endif

            crit(ipol,jpol) = h(ipol,jpol,iel) / ( velo * period)*dble(npol)
            crit_max(ipol,jpol) = hmin2(ipol,jpol,iel) / ( velo_max )

            theta = datan(dble(s1/(z1+1.d-30)))
            if ( 0.d0 > theta ) theta = pi + theta
            if (theta == 0.d0 .and. z1 < 0.d0) theta = pi

!         endif
         !  write(344,*) r,theta,crit(ipol,jpol,iel)

         ! TNM SEPT 2006: taken out....
         !   write(344,*) s1,z1,crit(ipol,jpol,iel)

         !   if ( abs(r*6371. -1221.5)<= 80.0) then
         !!   write(6,*)'TNMOUCHAUD: ICB',r*6371.,iel
         !   write(6,*)'TNMOUCHAUD: ICB  ',r*6371.,velo
         !   endif

      end do
   end do

   hmax(iel)=maxval(crit(:,:))*period/dble(npol)
   hmin(iel)=minval(crit_max(:,:))

   s1 = sgll(npol/2,npol/2,iel)
   z1 = zgll(npol/2,npol/2,iel)

   ! write(345,10)s1,z1,hmin(iel),hmax(iel)

   if (dump_mesh_info_files) then 

      ! check if grid spacing is within (numerically) allowable limits
      if (hmax(iel)>period/(pts_wavelngth*dble(npol))) then
         write(61,*)'WARNING +: grid spacing TOO LARGE in element',iel; call flush(61)
         write(61,*)'WARNING +: r,theta [km,deg]:',dsqrt(s1**2+z1**2)*router, &
              datan(s1/z1)*180.d0/pi; call flush(61)
         write(61,*)'WARNING +: max. allowed h/vs [s]:',period/(pts_wavelngth*npol)
         call flush(61)
         write(61,*)'WARNING +: actual h/vs [s] here :',hmax(iel); call flush(61)
         write(61,*)
         ntoobig=ntoobig+1
         ! should add a stop here after complete debugging....    
      endif

      if (hmin(iel)<dt/courant) then
         write(61,*)'WARNING -: grid spacing TOO SMALL in element',iel; call flush(61)
         write(61,*)'WARNING -: r,theta [km,deg]:',dsqrt(s1**2+z1**2)*router, &
              datan(s1/z1)*180.d0/pi; call flush(61)
         write(61,*)'WARNING -: max. allowed h/vp [s]:',dt/courant; call flush(61)
         write(61,*)'WARNING -: actual h/vp [s] here :',hmin(iel); call flush(61)
         write(61,*)  
         ntoosmall=ntoosmall+1
         ! should add a stop here after complete debugging....    
      endif
   end if

   ! save into vtk====================================
   if (dump_mesh_vtk) then

      x(ct+1)=sgll(0,0,iel)
      x(ct+2)=sgll(npol,0,iel)
      x(ct+3)=sgll(npol,npol,iel)
      x(ct+4)=sgll(0,npol,iel)
      y(ct+1)=zgll(0,0,iel)
      y(ct+2)=zgll(npol,0,iel)
      y(ct+3)=zgll(npol,npol,iel)
      y(ct+4)=zgll(0,npol,iel)

      if (bkgrdmodel=='solar') then 
         vp1(ct+1)=v_p(0,0,iel)
         vs1(ct+1)=v_s(0,0,iel)
         rho1(ct+1)=rho(0,0,iel)

         vp1(ct+2)=v_p(npol,0,iel)
         vs1(ct+2)=v_s(npol,0,iel)
         rho1(ct+2)=rho(npol,0,iel)

         vp1(ct+3)=v_p(npol,npol,iel)
         vs1(ct+3)=v_s(npol,npol,iel)
         rho1(ct+3)=rho(npol,npol,iel)

         vp1(ct+4)=v_p(0,npol,iel)
         vs1(ct+4)=v_s(0,npol,iel)
         rho1(ct+4)=rho(0,npol,iel)

      else 
         r=sqrt(  (x(ct+1))**2 + (y(ct+1))**2 )
         vp1(ct+1)=velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
         vs1(ct+1)=velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)      
         rho1(ct+1)=velocity(r*router,'rho',region(iel),bkgrdmodel,lfbkgrdmodel)    

         r=sqrt(  (x(ct+2))**2 + (y(ct+2))**2 )
         vp1(ct+2)=velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
         vs1(ct+2)=velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)      
         rho1(ct+2)=velocity(r*router,'rho',region(iel),bkgrdmodel,lfbkgrdmodel)    

         r=sqrt(  (x(ct+3))**2 + (y(ct+3))**2 )
         vp1(ct+3)=velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
         vs1(ct+3)=velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)      
         rho1(ct+3)=velocity(r*router,'rho',region(iel),bkgrdmodel,lfbkgrdmodel)    

         r=sqrt(  (x(ct+4))**2 + (y(ct+4))**2 )
         vp1(ct+4)=velocity(r*router,'v_p',region(iel),bkgrdmodel,lfbkgrdmodel)
         vs1(ct+4)=velocity(r*router,'v_s',region(iel),bkgrdmodel,lfbkgrdmodel)      
         rho1(ct+4)=velocity(r*router,'rho',region(iel),bkgrdmodel,lfbkgrdmodel)     

      endif

      mesh2(iel,1) = real(s1)
      mesh2(iel,2) = real(z1)
   endif

   ct = ct+4

   !=======================================

end do ! iel

if (bkgrdmodel=='solar') deallocate(v_p,v_s,rho)

if (dump_mesh_vtk) then
  fname=trim(diagpath)//'/mesh_vp'
  call write_VTK_bin_scal(x,y,z,vp1,npts_vtk/4,fname)
  deallocate(vp1)
  
  fname=trim(diagpath)//'/mesh_vs'
  call write_VTK_bin_scal(x,y,z,vs1,npts_vtk/4,fname)
  deallocate(vs1)
  
  fname=trim(diagpath)//'/mesh_rho'
  call write_VTK_bin_scal(x,y,z,rho1,npts_vtk/4,fname)
  deallocate(rho1)
    
  deallocate(x,y,z)

  h_real=real(hmax/(period/(pts_wavelngth*real(npol))))
  write(6,*)'minmax hmax:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_hmax'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)

  h_real=real(hmin/(dt/courant))
  write(6,*)'minmax hmin:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_hmin'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)

  h_real=real(period/hmax)
  write(6,*)'minmax pts wavelngth:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_pts_wavelength'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)

  h_real=real(dt/hmin)
  write(6,*)'minmax courant:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_courant'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)

  h_real=real(courant*hmin)
  write(6,*)'minmax dt:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_dt'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)

  h_real=real(pts_wavelngth*real(npol)*hmax)
  write(6,*)'minmax period:',minval(h_real),maxval(h_real)
  fname=trim(diagpath)//'/mesh_period'
  call write_VTK_bin_scal_old(h_real,mesh2,neltot,fname)
  deallocate(mesh2,h_real)
endif

char_time_max=maxval(hmax)
char_time_max_globel=maxloc(hmax,1)
r=dsqrt((sgll(npol/2,npol/2,maxloc(hmax,1)))**2 + &
       (zgll(npol/2,npol/2,maxloc(hmax,1)))**2 )
char_time_max_rad=r
char_time_max_theta=dasin(sgll(npol/2,npol/2,maxloc(hmax,1))/r)*180.d0/pi
write(6,*)'char max:',char_time_max_rad,char_time_max_theta,char_time_max

char_time_min=maxval(hmin)
char_time_min_globel=maxloc(hmin,1)
r=dsqrt(sgll(npol/2,npol/2,maxloc(hmin,1))**2 + &
       zgll(npol/2,npol/2,maxloc(hmin,1))**2)
char_time_min_rad=r
char_time_min_theta=dasin(sgll(npol/2,npol/2,maxloc(hmin,1))/r)*180.d0/pi
write(6,*)'char min:',char_time_min_rad,char_time_min_theta,char_time_min

if (dump_mesh_info_screen) then 
write(6,*)
write(6,*)'Characteristic min/max lead times (ratio h/v):'
write(6,*)'Max value[sec]/el number    :  ',char_time_max,char_time_max_globel
write(6,*)'Max location r[km],theta[deg]: ',char_time_max_rad*router/1000., &
                                            char_time_max_theta
write(6,*)'Min value[sec]/el number     : ',char_time_min,char_time_min_globel
write(6,*)'Min location r[km],theta[deg]: ',char_time_min_rad*router/1000., &
                                            char_time_min_theta
call flush(6)
end if

if (dump_mesh_info_screen) then 
 if (ntoobig > 0) then 
  write(6,*)
  write(6,*)'**********************************************************'
  write(6,*)'SERIOUS WARNING:',ntoobig,'elements are too LARGE!'; call flush(6)
  write(6,*)'                 ...up to', &
            (maxval(hmax)/(period/(pts_wavelngth*dble(npol)))-1.d0)*100.d0, 'percent!'
  write(6,*)'                 ...percent of total elements:',100.d0* &
                                 dble(ntoobig)/dble(neltot)
 
  write(6,*)'**********************************************************'
 endif

 if (ntoosmall > 0) then
  write(6,*)
  write(6,*)'**********************************************************'
  write(6,*)'SERIOUS WARNING:',ntoosmall,'elements are too SMALL!'; call flush(6)
  write(6,*)'                  ...up to', &
            ((dt/courant)/minval(hmin)-1.d0)*100.d0,'percent!'; call flush(6)
  write(6,*)'                 ...percent of total elements:',100.d0* &
                                 dble(ntoosmall)/dble(neltot)
  write(6,*)'**********************************************************'
 endif
 call flush(6)
end if

if (dump_mesh_info_files) then 
 if (ntoobig >0) then 
  write(61,*)'**********************************************************'
  write(61,*)'SERIOUS WARNING:',ntoobig,'elements are too LARGE!'; call flush(6)
  write(61,*)'                 ...up to', &
             (maxval(hmax)/(period/(pts_wavelngth*dble(npol)))-1.d0)*100.d0, 'percent!'
  write(61,*)'                 ...percent of total elements:',100.d0* &
                                dble(ntoobig)/dble(neltot)
  write(61,*)'**********************************************************'
 endif

 if (ntoosmall >0) then
  write(61,*)'**********************************************************'
  write(61,*)'SERIOUS WARNING:',ntoosmall,'elements are too SMALL!'; call flush(6)
  write(61,*)'                  ...up to', &
          ((dt/courant)/minval(hmin)-1.d0)*100.d0,'percent!'; call flush(6)
  write(61,*)'                 ...percent of total elements:',100.d0* &
                                  dble(ntoosmall)/dble(neltot)
  write(61,*)'**********************************************************'
 endif

 close(61)
end if
close(62)
deallocate(h,hmin2,crit,crit_max,hmin,hmax)

end subroutine bkgrdmodel_testing 
!--------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal_old(u2,mesh,rows,filename)

 implicit none
 integer*4 :: i,rows
 real, dimension(1:rows), intent(in) :: u2
 real, dimension(1:rows) :: u1
real, dimension(1:rows,2), intent(in) :: mesh
 integer, dimension(:),allocatable :: cell
 integer, dimension(:),allocatable :: cell_type
 character (len=55) :: filename;
 character (len=50) :: ss; !stream
 
!points structure
allocate(cell(rows*2),cell_type(rows))
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
write(100) mesh(i,1),mesh(i,2),0.0
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
end subroutine write_VTK_bin_scal_old
!-----------------------------------------------------------------------------


subroutine write_VTK_bin_scal(x,y,z,u1,elems,filename)
 implicit none
 integer*4 :: i,t,elems
 real*4, dimension(1:elems*4), intent(in) :: x,y,z,u1
 integer, dimension(:),allocatable :: cell
 integer, dimension(:),allocatable :: cell_type
 character (len=55) :: filename
 character (len=50) :: ss; !stream
!points structure
allocate(cell(elems*5),cell_type(elems))
do i=5,elems*5,5
 cell(i-4)=4;
 enddo
t=0
do i=5,elems*5,5
t=t+4
cell(i-3)=t-4;
cell(i-2)=t-3;
cell(i-1)=t-2;
cell(i)=t-1;
enddo

!do i=1,elems
! cell_type(i)=9
!enddo
cell_type=9
! write(6,*)'computing vtk file ',trim(filename),' ...'
open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')
write(100) '# vtk DataFile Version 4.0'//char(10)
write(100) 'mittico'//char(10)
write(100) 'BINARY'//char(10)
write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
write(ss,fmt='(A6,I10,A5)') 'POINTS',elems*4,'float'
write(100) ss//char(10)
!points
write(100) (x(i),y(i),z(i),i=1,elems*4)
write(100) char(10)
!cell topology
write(ss,fmt='(A5,2I10)') 'CELLS',elems,elems*5
write(100) char(10)//ss//char(10)
write(100) cell
write(100) char(10)
!cell type
write(ss,fmt='(A10,2I10)') 'CELL_TYPES',elems
write(100) char(10)//ss//char(10)
write(100) cell_type
write(100) char(10)
!data
write(ss,fmt='(A10,I10)') 'POINT_DATA',elems*4
write(100) char(10)//ss//char(10)
write(100) 'SCALARS data float 1'//char(10)
write(100) 'LOOKUP_TABLE default'//char(10) !color table?
write(100) u1
 close(100)
write(6,*)'...saved ',trim(filename)//'.vtk'
end subroutine write_VTK_bin_scal


!-----------------------------------------------------------------------------
subroutine arbitr_sub_solar_arr(s,z,v_p,v_s,rho,bkgrdmodel2)
!
! file-based, step-wise model in terms of domains separated by disconts.
! format:
! ndisc
! r vp vs rho
! ...
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double precision, intent(in) :: s(0:npol,0:npol,1:neltot),z(0:npol,0:npol,1:neltot)
character(len=100), intent(in) :: bkgrdmodel2
double precision, dimension(:,:,:), intent(out) :: rho(0:npol,0:npol,1:neltot)
double precision, dimension(:,:,:), intent(out) :: v_s(0:npol,0:npol,1:neltot)
double precision, dimension(:,:,:), intent(out) :: v_p(0:npol,0:npol,1:neltot)
double precision, allocatable, dimension(:) :: disconttmp,rhotmp,vstmp,vptmp
integer :: ndisctmp,i,ndisctmp2,ind(2),ipol,jpol,iel
logical :: bkgrdmodelfile_exists
double precision :: w(2),wsum,r0

! Does the file bkgrdmodel".bm" exist?
  inquire(file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm', &
          exist=bkgrdmodelfile_exists)
  if (bkgrdmodelfile_exists) then
     open(unit=77,file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm')
     read(77,*)ndisctmp
     allocate(disconttmp(1:ndisctmp))
     allocate(vptmp(1:ndisctmp),vstmp(1:ndisctmp),rhotmp(1:ndisctmp))
     do i=1, ndisctmp
        read(77,*)disconttmp(i),rhotmp(i),vptmp(i),vstmp(i)
     enddo
     close(77)
     do iel=1,neltot
        do jpol=0,npol
           do ipol=0,npol
              r0 = dsqrt(s(ipol,jpol,iel)**2 +z(ipol,jpol,iel)**2 )
              call interp_vel(r0,disconttmp(1:ndisctmp),ndisctmp,ind,w,wsum)
              rho(ipol,jpol,iel)=sum(w*rhotmp(ind))*wsum
              v_p(ipol,jpol,iel)=(w(1)*vptmp(ind(1))+w(2)*vptmp(ind(2)))*wsum
              v_s(ipol,jpol,iel)=sum(w*vstmp(ind))*wsum
           enddo
        enddo
     enddo
     deallocate(disconttmp,vstmp,vptmp,rhotmp)
  else 
     write(6,*)'Background model file', &
          bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm','does not exist!!!'
     stop
  endif

end subroutine arbitr_sub_solar_arr
!=============================================================================

!=============================================================================
subroutine interp_vel(r0,r,n,ind,w,wsum)

integer, intent(in) :: n
double precision, intent(in) :: r0,r(1:n)
integer, intent(out) :: ind(2)
double precision, intent(out) :: w(2),wsum
integer :: i,p
double precision :: dr1,dr2

p=1

     i=minloc(dabs(r-r0),1)
!     write(6,*)'INTERP;',r0,i
!     write(6,*)'INTERP:',r(i)

     if (r0>0.d0) then
        if ((r(i)-r0)/r0> 1.d-8) then ! closest discont. at larger radius
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))
        elseif ((r0-r(i))/r0> 1.d-8) then  ! closest discont. at smaller radius
           if (r0>maxval(r)) then ! for round-off errors where mesh is above surface
              ind(1)=i
              ind(2)=i
              dr1=1.d0
              dr2=1.d0
           else
              ind(1)=i-1
              ind(2)=i
              dr1=r(ind(1))-r0
              dr2=r0-r(ind(2))
            endif
        elseif (dabs((r(i)-r0)/r0)< 1.d-8) then ! closest discont identical
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           write(6,*)'problem with round-off errors in interpolating......'
           write(6,*)'r0,r(i),i',r0,r(i),abs((r(i)-r0)/r0),i
           stop
        endif
     else !r0=0
        if (r(i)==0.d0) then ! center of the sun
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))        
        endif
     endif

! inverse distance weighting
     w(1)=(dr1)**(-p)
     w(2)=(dr2)**(-p)
     wsum=1.d0/sum(w)

end subroutine interp_vel
!=============================================================================


!=============================
end module test_bkgrdmodel
!=============================
