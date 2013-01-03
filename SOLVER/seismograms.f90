!> Various subroutines for seismogram preparation and dumping
!========================
module seismograms
!========================

use global_parameters
use data_io
use data_mesh
use data_proc
use data_time

implicit none

public

contains
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine prepare_seismograms

use utlity
use data_mesh_preloop
use data_source, ONLY : have_src
use commun, ONLY : psum_int,barrier
character(len=4)             :: appielem
integer          :: j,iel,ielem,ind,epicount,equcount,anticount,iproc,maxind_glob
double precision :: s,z,r,theta

  if (lpr) &
  write(6,*)'  locating surface elements and generic receivers...'

  ielantipode=1; ielequ=1
  have_epi=.false.
  have_equ=.false.
  have_antipode=.false.
  epicount=0; equcount=0;anticount=0
  ind=0

! test round to determine the amount of surface elements
  do iel=1,nel_solid
     ielem=ielsolid(iel)
     if (north(ielem)) then
        call compute_coordinates(s,z,r,theta,ielem,0,npol)
     else
        call compute_coordinates(s,z,r,theta,ielem,npol,0)
     endif
     if ( dabs(r-router)< min_distance_dim ) then
        ind=ind+1
     endif
  enddo

  num_surf_el = ind
  allocate(surfelem(num_surf_el))

  open(1000+mynum,file=datapath(1:lfdata)//'/surfelem_'//appmynum//'.dat') 
  ind=0
  do iel=1,nel_solid
     ielem=ielsolid(iel)
     if (north(ielem)) then
        call compute_coordinates(s,z,r,theta,ielem,0,npol)
     else
        call compute_coordinates(s,z,r,theta,ielem,npol,0)
     endif

     if ( dabs(r-router)< min_distance_dim ) then
        ind=ind+1
        surfelem(ind)=iel
        write(1000+mynum,*)ind,iel,r,theta*180./pi

! find epicenter
        if (north(ielem)) then
           call compute_coordinates(s,z,r,theta,ielem,0,npol)
           if ( dabs(s) < smallval*router .and. &
                dabs(z-router) < smallval*router)  then
              ielepi=iel
              write(69,*)'Epicenter element:',ielepi
              write(69,*)'Epicenter radius [km], colat [deg]:', &
                         r/1000.,theta/pi*180. 
              write(69,*)''
              have_epi=.true.
              epicount=epicount+1
           end if
        endif

! find antipode at 180deg
       if (.not. north(ielem)) then
          call compute_coordinates(s,z,r,theta,ielem,0,0)
          if ( dabs(theta-pi) < min_distance_nondim*pi) then 
             ielantipode=iel
             write(69,*)'Antipodal element:',ielantipode
             write(69,*)'Antipode radius [km], colat [deg]:',&
                         r/1000.,theta/pi*180.
             write(69,*)''
             have_antipode=.true.
             anticount=anticount+1
          end if
        endif

! find equator (take northern element)
        if (north(ielem)) then
           call compute_coordinates(s,z,r,theta,ielem,npol,npol)
!          if ( dabs(theta-pi/2.d0) < min_distance_nondim*pi/2.d0)  then 
           if ( dabs(z) < 1.e-8)  then 
              ielequ=iel
              write(69,*)'Equatorial element:',ielequ
              write(69,*)'Equatorial radius [km], colat [deg]:',&
                          r/1000.,theta/pi*180.
              write(69,*)''
              have_equ=.true.
              equcount=equcount+1
           end if
        endif

     endif
  enddo
  maxind=ind

  close(1000+mynum)

! make sure only one processor has each location
  if (lpr) &
     write(6,*)'  ensuring uniqueness in generic receiver locations...'

  iel=0
  if ( psum_int(epicount) > 1 ) then 
     do while (iel <= nproc-1)
        call barrier 
        if (mynum==iel) then 
           if (have_epi .and. mynum<nproc-1) then 
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum==iproc) have_epi=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
        iel=iel+1
     enddo
  endif

  if ( psum_int(equcount) > 1 ) then 
     do while (iel<=nproc-1)
        call barrier 
        if (mynum==iel) then 
           if (have_equ .and. mynum<nproc-1) then 
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum==iproc) have_equ=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
        iel=iel+1
     enddo
  endif

  if ( psum_int(anticount) > 1 ) then 
     do while (iel<=nproc-1)
        call barrier 
        if (mynum==iel) then 
           if (have_antipode .and. mynum<nproc-1) then 
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum==iproc) have_antipode=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
     enddo
  endif

  do iel=0,nproc-1
     call barrier
     if (mynum==iel) then
        write(69,*)'  number of surface elements:',maxind
        if (have_epi) write(6,12)procstrg,'epicenter at', &
                                thetacoord(0,npol,ielsolid(ielepi))/pi*180.
        if (have_equ) write(6,12)procstrg,'equator at', &
                                thetacoord(npol,npol,ielsolid(ielequ))/pi*180.
        if (have_antipode) write(6,12)procstrg,'antipode at', &
                                  thetacoord(0,0,ielsolid(ielantipode))/pi*180.
     endif
     call barrier
  enddo
12 format('   ',a8,'has the ',a13,f9.3,' degrees.')

  maxind_glob = psum_int(maxind)

  if (lpr) write(6,*)'  global number of surface elements:',maxind_glob

! open files for displacement and velocity traces in each surface element

  allocate(jsurfel(maxind)) ! for surface strain
  allocate(surfcoord(maxind)) ! theta of surface elements dumped for surface strain

     do iproc=0,nproc-1
        call barrier
        if (mynum==iproc) then 
           call barrier
           open(33333,file=datapath(1:lfdata)// &
                            '/surfelem_coords.dat',position='append')
           open(33334,file=datapath(1:lfdata)// &
                            '/surfelem_coords_jpol.dat',position='append')
           if (mynum==0) write(33334,*)maxind_glob
           if (mynum==0) write(33333,*)maxind_glob
           do iel=1,maxind
              if (thetacoord(npol/2,npol/2,ielsolid(surfelem(iel)))<=pi/2) then
                 jsurfel(iel)=npol
              else
                 jsurfel(iel)=0
              endif
              write(33333,*) 180./pi* &
                   thetacoord(npol/2,jsurfel(iel),ielsolid(surfelem(iel)))

              surfcoord(iel) = 180. / pi * &
                               thetacoord(npol/2,jsurfel(iel),ielsolid(surfelem(iel)))

              write(33334,11) 180./pi* &
                   thetacoord(npol/2,jsurfel(iel),ielsolid(surfelem(iel))),&
                   (rcoord(npol/2,j,ielsolid(surfelem(iel))),j=0,npol)
           enddo
           close(33333); close(33334)
           call barrier
        endif
        call barrier
     enddo

11 format(6(1pe11.4))

  if (dump_wavefields) then
    if(use_netcdf)   then

    else 

       do iel=1,maxind
          call define_io_appendix(appielem,iel+mynum*maxind)
          open(unit=40000000+iel,file=datapath(1:lfdata)// &
                    '/surfelem_disp.dat'//appielem)


          open(unit=50000000+iel,file=datapath(1:lfdata)// &
                                      '/surfelem_velo.dat'//appielem)
       enddo

       do iel=1,maxind
          call define_io_appendix(appielem,iel+mynum*maxind)
          open(unit=60000000+iel,file=datapath(1:lfdata)// &
               '/surfelem_strain.dat'//appielem)
          open(unit=70000000+iel,file=datapath(1:lfdata)// &
               '/surfelem_disp_src.dat'//appielem)
       enddo
    endif
  endif

end subroutine prepare_seismograms
!=============================================================================

!-----------------------------------------------------------------------------
!> Read colatitudes [deg] from a file receivers.dat and locate closest grid 
!! point for seismograms, output grid point locations in 
!! receiver_pts.dat<PROCID>
subroutine prepare_from_recfile_seis
!
! Read colatitudes [deg] from a file receivers.dat and locate closest grid 
! point for seismograms, output grid point locations in 
! receiver_pts.dat<PROCID>
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use utlity
use data_mesh_preloop
use data_mesh, ONLY: loc2globrec
use data_source, ONLY : src_type,rot_src,srccolat,srclon
use commun
use rotations, ONLY : rotate_receivers_recfile,save_google_earth_kml
use nc_routines, ONLY: nc_define_outputfile

integer                      :: i,iel,ipol,irec,num_rec_glob
integer                      :: count_diff_loc,count_procs
double precision             :: s,z,r,theta,recdist,myrecdist
double precision             :: tmprecfile_th,tmprecfile_el(3)
double precision,allocatable :: recfile_readth(:),recfile_readph(:)
double precision,allocatable :: recfile_th_glob(:)
double precision,allocatable :: recfile_th_loc(:),recfile_el_loc(:,:)
double precision,allocatable :: recfile_th(:),recfile_ph_loc(:),recfile_ph_loc2(:)
integer,allocatable          :: rec2proc(:),loc2globrec_loc(:)!, loc2globrec(:)
character(len=4)             :: appielem
character(len=100)           :: junk
character(len=40), allocatable :: receiver_name(:)
double precision             :: maxreclocerr

! Additional arrays within the STATIONS file
character(len=20), allocatable, dimension(:) :: rec_name,rec_network
double precision, allocatable, dimension(:) :: reclat,reclon,recelevation,recbury
integer ierror

    irec=0
    count_diff_loc=0
    dtheta_rec=0. ! default for all cases but database (which has regular spacing)

! Read receiver location file:
! line1: number of receivers
! line2 --> line(number of receivers+1): colatitudes [deg]

    if (rec_file_type=='colatlon') then
       if (lpr) write(6,*)'  reading receiver colatitudes and longitudes from receivers.dat...'
       open(unit=34, file='receivers.dat', iostat=ierror, status='old', position='rewind')
       read(34,*) num_rec_glob
       write(6,*)'  # total receivers:',num_rec_glob; call flush(6)
       allocate( recfile_readth (1:num_rec_glob) )
       allocate( recfile_th_glob(1:num_rec_glob) )
       allocate( recfile_th_loc (1:num_rec_glob) )
       allocate( recfile_el_loc (1:num_rec_glob,3))
       allocate( loc2globrec_loc(1:num_rec_glob) )
       allocate( rec2proc       (1:num_rec_glob) )
       allocate( recfile_readph (1:num_rec_glob) )
       allocate( recfile_ph_loc2(1:num_rec_glob) )   
       allocate( receiver_name  (1:num_rec_glob) )

       open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')
!       call BSORT2(recfile_readth,num_rec_glob) ! need to sort the indices as well....

       do i=1,num_rec_glob
          read(34,*)recfile_readth(i),recfile_readph(i)
          call define_io_appendix(appielem,i)
          receiver_name(i) = 'recfile_'//appielem
          write(30,*) trim(receiver_name(i)), recfile_readth(i), recfile_readph(i)
       enddo
       close(34) 
       close(30)

    elseif (rec_file_type=='database') then 
       if (lpr) write(6,*)'  generating receiver colatitudes at every element edge and midpoint...'
       if (lpr) write(6,*)'  ... which is useful for databases and sinc interpolation'
       dtheta_rec = abs( (thetacoord(npol,npol,ielsolid(surfelem(1))) &
                        - thetacoord(0,   npol,ielsolid(surfelem(1))) ) * 180. / pi ) / 2
       num_rec_glob = ceiling(180./dtheta_rec)+1
       if (lpr) then
         write(6,*) 'delta theta (mid-to-edge), (edge-to-edge) [deg]:', dtheta_rec, &
                    abs( (thetacoord(npol,npol,ielsolid(surfelem(1))) &
                        - thetacoord(0,   npol,ielsolid(surfelem(1))) ) * 180. / pi )
         write(6,*)mynum,'number of surface elements:',maxind
         write(6,*)mynum,'number of global recs (ideal,real):',(180./dtheta_rec)+1,num_rec_glob
       end if
       allocate(recfile_readth(num_rec_glob),recfile_readph(num_rec_glob))
       do i=1,num_rec_glob
          recfile_readth(i) = dtheta_rec*real(i-1)
       enddo
       if(lpr) write(6,*)mynum,'min,max receiver theta [deg]:',minval(recfile_readth),maxval(recfile_readth)
       call flush(6)
       recfile_readph = 0.
       allocate(recfile_th_loc(1:num_rec_glob),recfile_el_loc(1:num_rec_glob,3))
       allocate(loc2globrec_loc(1:num_rec_glob),rec2proc(1:num_rec_glob))
       allocate(recfile_ph_loc2(1:num_rec_glob),recfile_th_glob(1:num_rec_glob))   
       allocate(receiver_name(1:num_rec_glob))
       if (mynum==0) open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')
       do i=1,num_rec_glob
          !call define_io_appendix(appielem,i) Does only work for nrec<9999
          !receiver_name(i) = 'recfile_'//appielem
          write(receiver_name(i),112) i
112       format('recfile_',I6.6)          
          if (mynum==0) write(30,*)trim(receiver_name(i)),recfile_readth(i),recfile_readph(i)
       enddo
       if (mynum==0) close(30)
       if (lpr) write(6,*)mynum,'done with database receiver writing.';call flush(6)

     elseif (rec_file_type=='stations') then 
       if (lpr) write(6,*)'  reading receiver colatitudes and longitudes from STATIONS...'
       ! count number of receivers first
       open(unit=34,file='STATIONS',iostat=ierror,status='old',action='read',position='rewind')
       num_rec_glob = 0
       do while(ierror == 0)
          read(34,*,iostat=ierror) junk
          if(ierror == 0) num_rec_glob = num_rec_glob + 1
       enddo
       close(34)
       if (lpr) write(6,*)'  ...counted number of stations:', num_rec_glob

       allocate(recfile_readth(1:num_rec_glob),recfile_th_glob(1:num_rec_glob))
       allocate(recfile_th_loc(1:num_rec_glob),recfile_el_loc(1:num_rec_glob,3))
       allocate(loc2globrec_loc(1:num_rec_glob),rec2proc(1:num_rec_glob))
       allocate(recfile_readph(1:num_rec_glob))
       allocate(recfile_ph_loc2(1:num_rec_glob))     
       allocate(rec_name(num_rec_glob),rec_network(num_rec_glob))
       allocate(reclat(num_rec_glob),reclon(num_rec_glob),recelevation(num_rec_glob),recbury(num_rec_glob))
       allocate(receiver_name(num_rec_glob))
       open(unit=34,file='STATIONS',iostat=ierror,status='old',action='read',position='rewind')
       open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')

       do i=1,num_rec_glob
          read(34,*)rec_name(i),rec_network(i),reclat(i),reclon(i),recelevation(i),recbury(i)
          if (reclon(i)<=zero) then 
             recfile_readph(i)=reclon(i)+360.d0 
          else
             recfile_readph(i)= reclon(i)
          endif
          recfile_readth(i) = 90.d0 - reclat(i)
          receiver_name(i) = trim(rec_name(i))//'_'//trim(rec_network(i))
          write(30,*)trim(receiver_name(i)),recfile_readth(i),recfile_readph(i)
       enddo
       close(34); close(30)

    elseif (rec_file_type=='deliverg') then 
       if (lpr) write(6,*)'  reading receiver colatitudes/longitudes from deliverg...'
       if (lpr) write(6,*)'  .....not yet implemented!'
       stop
    else 
       write(6,*)procstrg, 'Undefined receiver file format!!'
       stop
    endif !Receiver type rec_file_type

    num_rec_tot = num_rec_glob ! to be known later/globally

! check on consistency of receiver coordinates
    if (maxval(recfile_readph)<smallval .and. correct_azi) then 
       if (lpr) write(6,*)'  WARNING: receivers have zero azimuth but you wish to rotate them..'
!!$    elseif (maxval(recfile_readph)>smallval .and. .not. correct_azi) then 
!!$       if (lpr) write(6,*)'  WARNING: receivers have non-zero azimuth but you chose to ignore it...'
!!$       if (rot_src .and. lpr) then
!!$          write(6,*)'  WARNING: even worse, the receivers will be rotated due to a non-axial source location'
!!$          write(6,*)'                        ... and zero azimuth will then even affect the epicentral distance! Bad!'
!!$       endif
    endif

    if (minval(recfile_readph) < 0.d0) then 
       if (lpr) write(6,*)' ERROR: We do not allow negative receiver longitudes....'
       stop
    endif

    if (maxval(recfile_readph) > 360.001) then 
       if (lpr) write(6,*)' ERROR: We do not allow receiver longitudes larger than 360 degrees....'
       stop
    endif

    if (maxval(recfile_readth) < 0.d0) then 
       if (lpr) write(6,*)' ERROR: We do not allow negative receiver colatitudes....'
       stop
    endif

    if (maxval(recfile_readth) > 180.001) then 
       if (lpr) write(6,*)' ERROR: We do not allow receiver colatitudes larger than 180 degrees....'
       stop
    endif

! rotate receiver locations if source is not located at north pole
    if (rot_src ) then 
       call rotate_receivers_recfile(num_rec_glob,recfile_readth,recfile_readph,receiver_name)
    else
      if ((lpr).and.(.not.(rec_file_type=='database'))) then
        call save_google_earth_kml( real(srccolat*180.0/pi), real(srclon*180.d0/pi), &
                                    real(recfile_readth), real(recfile_readph), &
                                    num_rec_glob, 'original', receiver_name  )
      end if
    endif

  recfile_th_glob(:) = zero
  rec2proc(:) = 0

! find closest grid points

!=====================
  do i=1,num_rec_glob
!=====================
     write(69,*)'  working on receiver #',i,recfile_readth(i)*180./pi; call flush(69)
     recdist=10.d0*router
     do iel=1,maxind
        do ipol=0,npol
           if (north(ielsolid(surfelem(iel)))) then ! NORTH
              call compute_coordinates(s,z,r,theta,ielsolid(surfelem(iel)),&
                                       ipol,npol)
              if (z < zero ) then
               write(6,*)'PROBLEM! north but z<0: ', &
                           north(ielsolid(surfelem(iel))),z
                 write(6,*)'r,theta:',r/1000.,theta*180./pi
                 write(6,*)iel,surfelem(iel),ielsolid(surfelem(iel))
                 stop
              endif

           else ! SOUTH
              call compute_coordinates(s,z,r,theta,ielsolid(surfelem(iel)),&
                                       ipol,0)
              if (z > zero ) then
                 write(6,*)'PROBLEM! south but z>0: ',&
                           north(ielsolid(surfelem(iel))),z
                 write(6,*)'r,theta:',r/1000.,theta*180./pi
                 write(6,*)iel,surfelem(iel),ielsolid(surfelem(iel))
                 stop
              endif

           endif

           if (dabs(theta/pi*180.d0-recfile_readth(i)) < recdist) then
              recdist=dabs(theta/pi*180.d0-recfile_readth(i))
              tmprecfile_th=theta/pi*180.d0
              tmprecfile_el(1)=surfelem(iel) ! only in the solid domain
              tmprecfile_el(2)=ipol
              if (north(ielsolid(surfelem(iel)))) tmprecfile_el(3)=npol
              if (.not. north(ielsolid(surfelem(iel)))) tmprecfile_el(3)=0 
           endif

        enddo
     enddo

! Make sure only one processor takes on each location
     myrecdist=recdist
     recdist=pmin(recdist)
     count_procs=0
     if (dblreldiff_small(myrecdist,recdist)) count_procs=mynum
!    take as default the larger processor ID to take on the receiver
     count_procs=pmax_int(count_procs)
     if (mynum==count_procs) then 
        irec=irec+1
        write(69,*)'found local grid point and processor...',irec,i; call flush(69)
        recfile_th_loc(irec)=tmprecfile_th
        recfile_el_loc(irec,1:3)=tmprecfile_el(1:3)
        loc2globrec_loc(irec)=i
        rec2proc(i)=mynum
        recfile_th_glob(i)=tmprecfile_th
     endif

! longitude
     if (irec > 0)  recfile_ph_loc2(irec) = recfile_readph(i)*pi/180.

! Can do that since NOW only one proc has non-zero values
     recfile_th_glob(i) = psum_dble(recfile_th_glob(i))
     rec2proc(i) = psum_int(rec2proc(i))

!=====================
  enddo ! num_rec_glob
!=====================

! Form local arrays depending on how many receivers each processor has
  num_rec=irec
  allocate(recfile_el(1:num_rec,1:3),loc2globrec(1:num_rec))
  allocate(recfile_th(1:num_rec))
  allocate(recfile_ph_loc(1:num_rec))
!af
  allocate(fname_rec_seis(1:num_rec))
  allocate(fname_rec_velo(1:num_rec))
!end af

  recfile_el(1:num_rec,1:3)=recfile_el_loc(1:num_rec,1:3)
  loc2globrec(1:num_rec)=loc2globrec_loc(1:num_rec)
  recfile_th(1:num_rec)=recfile_th_loc(1:num_rec)
  recfile_ph_loc(1:num_rec)=recfile_ph_loc2(1:num_rec)
  deallocate(recfile_ph_loc2)

! How many receivers does each processor have, do they sum to global number?
  if ( psum_int(num_rec) /= num_rec_glob ) then
     write(6,*)'PROBLEM: sum of local receivers is different than global!'
     if (lpr) write(6,*)'Global number of receivers:',num_rec_glob
     write(6,*)procstrg,'Number of receivers:',num_rec
     stop
  endif

  if (lpr) write(6,*)
  do irec=0,nproc-1
     call barrier
     if (mynum==irec) write(6,14)procstrg,num_rec,num_rec_glob
     call barrier
  enddo
  write(69,*)
  write(69,14)procstrg,num_rec,num_rec_glob
14 format('   ',a8,'has',i4,' out of',i6,' receivers')

! Output colatitudes globally (this is the file needed to plot seismograms)
! NOTE: recfile_th is in degrees!!!
  if (lpr) then
     open(99997,file=datapath(1:lfdata)//'/receiver_pts.dat')
     do i=1,num_rec_glob
        write(99997,*)recfile_th_glob(i),recfile_readph(i),rec2proc(i)
     enddo
     close(99997)
  endif

! Output colatitudes locally to infopath (this file is for info purposes!)

  maxreclocerr=zero
  open(9998+mynum,file=infopath(1:lfinfo)//'/receiver_pts_'//appmynum//'.dat')
  write(9998+mynum,*)num_rec
  do i=1,num_rec ! Only over newly found local receiver locations

     write(9998+mynum,13) i, recfile_readth(loc2globrec(i)), recfile_th(i), &
                          recfile_ph_loc(i),                                &
                          recfile_el(i,1), recfile_el(i,2), recfile_el(i,3)

     call define_io_appendix(appielem,loc2globrec(i))

     if ( pi/180*router*abs(recfile_readth(loc2globrec(i))-recfile_th(i) ) >  &
          min_distance_dim) then
        count_diff_loc=count_diff_loc+1
        write(6,22)procstrg, &
                  recfile_readth(loc2globrec(i)), recfile_th(i)
        if (dabs(recfile_readth(loc2globrec(i))-recfile_th(i))> maxreclocerr) &
             maxreclocerr=dabs(recfile_readth(loc2globrec(i))-recfile_th(i))/ &
                          180.*pi*router

     endif

22 format('   WARNING:',a8,' rec. location file/mesh:',2(f9.3))

! Test: receivers on the surface?
     call compute_coordinates(s,z,r,theta,ielsolid(recfile_el(i,1)), &
          recfile_el(i,2),recfile_el(i,3))
     if ( .not. dblreldiff_small(r,router)) then
        write(6,*)''
        write(6,*)'PROBLEM: receiver is not at the surface!'
        write(6,*)'r [km], colat [deg]:',r/1000.,theta/pi*180.
        stop
     endif


  if(.not.(use_netcdf))   then
     write(6,*)'  ',procstrg,'opening receiver file:',i,appielem
     open(100000+i,file=datapath(1:lfdata)//'/'//trim(receiver_name(loc2globrec(i)))//'_disp.dat')
  endif
 
  enddo
  close(9998+mynum)
  if (use_netcdf) then
    call nc_define_outputfile(num_rec_glob, receiver_name, recfile_th_glob, recfile_readth, recfile_readph, rec2proc)
  end if

  write(69,15)count_diff_loc,num_rec
  write(69,*)'  Maximal receiver location error [m]:',maxreclocerr
  write(69,*)
15 format(i4,' out of',i4,' receivers are located at wrong points.')

  maxreclocerr=pmax(maxreclocerr)
  if (lpr) then 
     write(6,*)
     write(6,*)'  maximal receiver location error [m]:',maxreclocerr
     write(6,*)
  endif

! define general prefactor for all cases
allocate(recfac(num_rec,5))

if (rot_rec=='sph') then ! rotating to r,phi,theta
   if (lpr) write(6,*) &
            '  Calculating prefactors for rotating components to spherical...'
   recfac(:,1) =  sin(recfile_th*pi/180.)
   recfac(:,2) =  cos(recfile_th*pi/180.)
   recfac(:,3) =  1.d0
   recfac(:,4) =  cos(recfile_th*pi/180.)
   recfac(:,5) = -sin(recfile_th*pi/180.)

elseif (rot_rec=='enz') then ! rotating to East North Z
   if (lpr) write(6,*) &
            '  Calculating prefactors for rotating components to NEZ...'
   recfac(:,1) =  sin(recfile_th*pi/180.)
   recfac(:,2) =  cos(recfile_th*pi/180.)
   recfac(:,3) =  1.d0
   recfac(:,4) = -cos(recfile_th*pi/180.)
   recfac(:,5) =  sin(recfile_th*pi/180.)

elseif (rot_rec=='xyz') then ! rotating to Greenwich xyz
   if (lpr) write(6,*) &
            '  Calculating prefactors for rotating components to xyz (global)...NOT DONE YET'

   stop

elseif  (rot_rec=='cyl') then ! keeping s,phi,z
   if (lpr) write(6,*) &
            '  Calculating prefactors for cylindrical components...'
   recfac(:,:) = 1.d0
   recfac(:,2) = 0.d0
   recfac(:,4) = 0.d0
else 
   write(6,*)'Unkown coordinate system for receiver components:',rot_rec
   write(6,*)'Please choose from cyl,sph,enz,xyz'
   stop
endif

! no phi component for monopole
if (src_type(1)=='monopole')  recfac(:,3) = 0.d0

! multiply seismogram with the correct azimuthal factor
if (correct_azi) then 

   if (src_type(2)=='mxz' .or. src_type(2)=='xforce') then 
   if (lpr) write(6,*) &
            '  Calculating prefactors for correct azimuth for',src_type(2)
      recfac(:,1) =  cos(recfile_ph_loc)*recfac(:,1)
      recfac(:,2) =  cos(recfile_ph_loc)*recfac(:,2)
      recfac(:,3) = -sin(recfile_ph_loc)*recfac(:,3)
      recfac(:,4) =  cos(recfile_ph_loc)*recfac(:,4)
      recfac(:,5) =  cos(recfile_ph_loc)*recfac(:,5)

   elseif (src_type(2)=='myz' .or. src_type(2)=='yforce') then 
   if (lpr) write(6,*) &
            '  Calculating prefactors for correct azimuth for',src_type(2)
      recfac(:,1) =  sin(recfile_ph_loc)*recfac(:,1)
      recfac(:,2) =  sin(recfile_ph_loc)*recfac(:,2)
      recfac(:,3) =  cos(recfile_ph_loc)*recfac(:,3)
      recfac(:,4) =  sin(recfile_ph_loc)*recfac(:,4)
      recfac(:,5) =  sin(recfile_ph_loc)*recfac(:,5)
      
   elseif (src_type(2)=='mxx_m_myy') then 
   if (lpr) write(6,*) &
            '  Calculating prefactors for correct azimuth for',src_type(2)
      recfac(:,1) =  cos(2.d0*recfile_ph_loc)*recfac(:,1)
      recfac(:,2) =  cos(2.d0*recfile_ph_loc)*recfac(:,2)
      recfac(:,3) = -sin(2.d0*recfile_ph_loc)*recfac(:,3)
      recfac(:,4) =  cos(2.d0*recfile_ph_loc)*recfac(:,4)
      recfac(:,5) =  cos(2.d0*recfile_ph_loc)*recfac(:,5)
      
   elseif (src_type(2)=='mxy') then 
   if (lpr) write(6,*) &
            '  Calculating prefactors for correct azimuth for',src_type(2)
      recfac(:,1) =  sin(2.d0*recfile_ph_loc)*recfac(:,1)
      recfac(:,2) =  sin(2.d0*recfile_ph_loc)*recfac(:,2)
      recfac(:,3) =  cos(2.d0*recfile_ph_loc)*recfac(:,3)
      recfac(:,4) =  sin(2.d0*recfile_ph_loc)*recfac(:,4)
      recfac(:,5) =  sin(2.d0*recfile_ph_loc)*recfac(:,5)
   endif

endif 

open(300+mynum,file=infopath(1:lfinfo)//'/receiver_recfac_'//appmynum//'.dat')
do i=1,num_rec
   write(300+mynum,*)recfile_th(i),recfile_ph_loc(i)*180./pi,recfac(i,3)
   write(300+mynum,*)recfac(i,1),recfac(i,2),recfac(i,4),recfac(i,5)
   write(300+mynum,*)
enddo
close(300+mynum)

13 format(i3,3(1pe12.4),i8,2(i2))

  deallocate(recfile_ph_loc,recfile_readph)
  deallocate(recfile_readth)
  deallocate(recfile_th_glob,recfile_th,recfile_th_loc)
  deallocate(recfile_el_loc,loc2globrec_loc,rec2proc)

if (rec_file_type=='stations') then 
  deallocate(rec_name,rec_network)
  deallocate(reclat,reclon,recelevation,recbury)
endif

end subroutine prepare_from_recfile_seis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine prepare_from_recfile_cmb
!
! Read colatitudes [deg] from a file receivers.dat and locate closest grid 
! point for seismograms, output grid point locations in 
! receiver_pts.dat<PROCID>
!
! This is the exact same routine as the above, but instead of finding 
! locations on the earth's surface rather locating them in the element 
! above the CMB, with the same colatitude.
! As for the surface, take "above" location, i.e. one element above the CMB
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use utlity
use data_mesh_preloop 
use commun

integer                      :: i,iel,ipol,icmb,num_cmb_glob
integer                      :: count_diff_loc,count_procs,ielglob
double precision             :: s,z,r,theta,cmbdist,mycmbdist
double precision             :: tmpcmbfile_th,tmpcmbfile_r,tmpcmbfile_el(3)
double precision,allocatable :: cmbfile_readth(:),cmbfile_th_glob(:)
double precision,allocatable :: cmbfile_r_glob(:),cmbfile_r_loc(:)
double precision,allocatable :: cmbfile_th_loc(:),cmbfile_el_loc(:,:)
double precision,allocatable :: cmbfile_th(:),cmbfile_r(:)
integer,allocatable          :: cmb2proc(:),loc2globcmb_loc(:),loc2globcmb(:)
character(len=4)             :: appielem
double precision             :: maxcmblocerr

  icmb=0
  count_diff_loc=0

  if (lpr) write(6,*)'  reading cmbrec colatitudes from receivers.dat...'

! Read receiver location file:
! line1: number of receivers
! line2 --> line(number of receivers+1): colatitudes [deg]
  open(99999,file='receivers.dat',status='old',POSITION='REWIND')
    read(99999,*)num_cmb_glob
    allocate(cmbfile_readth(1:num_cmb_glob),cmbfile_th_glob(1:num_cmb_glob))
    allocate(cmbfile_th_loc(1:num_cmb_glob),cmbfile_el_loc(1:num_cmb_glob,3))
    allocate(cmbfile_r_glob(1:num_cmb_glob),cmbfile_r_loc(1:num_cmb_glob))
    allocate(loc2globcmb_loc(1:num_cmb_glob),cmb2proc(1:num_cmb_glob))
    do i=1,num_cmb_glob
       read(99999,*)cmbfile_readth(i)
    enddo
  close(99999)

  cmbfile_th_glob(:) = zero
  cmbfile_r_glob(:) = zero
  cmb2proc(:) = 0

! find closest grid points

!=====================
  do i=1,num_cmb_glob
!=====================
     cmbdist=10.d0*router

  do iel=1,nel_bdry

! Map from boundary to global indexing. Choice of the solid side is random...
     ielglob=ielsolid(bdry_solid_el(iel))

        do ipol=0,npol
           if (north(ielglob)) then ! NORTH
              call compute_coordinates(s,z,r,theta,ielglob,ipol,npol)
              if (z < zero ) then
                 write(6,*)'PROBLEM! north but z<0: ', north(ielglob),z
                 write(6,*)'r,theta:',r/1000.,theta*180/pi
                 write(6,*)iel,bdry_solid_el(iel),ielglob
                 stop
              endif

           else ! SOUTH
              call compute_coordinates(s,z,r,theta,ielglob,ipol,0)
              if (z > zero ) then
                 write(6,*)'PROBLEM! south but z>0: ',north(ielglob),z
                 write(6,*)'r,theta:',r/1000.,theta*180/pi
                 write(6,*)iel,bdry_solid_el(iel),ielglob
                 stop
              endif

           endif
           
! Make sure that we're only dealing with the CMB, not ICB
           if (r>=3480000.) then 

              if (dabs(theta/pi*180.d0-cmbfile_readth(i)) < cmbdist) then
                 cmbdist=dabs(theta/pi*180.d0-cmbfile_readth(i))
                 tmpcmbfile_th=theta/pi*180.d0
                 tmpcmbfile_r=r
                 tmpcmbfile_el(1)=bdry_solid_el(iel) ! only in the solid domain
                 tmpcmbfile_el(2)=ipol
                 if (north(ielglob)) tmpcmbfile_el(3)=npol
                 if (.not. north(ielglob)) tmpcmbfile_el(3)=0 
              endif
              
           endif

        enddo

     enddo ! boundary elements

! Make sure only one processor takes on each location
     mycmbdist=cmbdist
     cmbdist=pmin(cmbdist)
     count_procs=0
     if (dblreldiff_small(mycmbdist,cmbdist)) count_procs=mynum
!    take as default the larger processor ID to take on the receiver
     count_procs=pmax_int(count_procs)
     if (mynum==count_procs) then 
        icmb=icmb+1
        cmbfile_th_loc(icmb)=tmpcmbfile_th
        cmbfile_r_loc(icmb)=tmpcmbfile_r
        cmbfile_el_loc(icmb,1:3)=tmpcmbfile_el(1:3)
        loc2globcmb_loc(icmb)=i
        cmb2proc(i)=mynum
        cmbfile_th_glob(i)=tmpcmbfile_th
        cmbfile_r_glob(i)=tmpcmbfile_r
     endif

! Can do that since NOW only one proc has non-zero values
     cmbfile_th_glob(i) = psum_dble(cmbfile_th_glob(i))
     cmbfile_r_glob(i) = psum_dble(cmbfile_r_glob(i))
     cmb2proc(i) = psum_int(cmb2proc(i))

!=====================
  enddo ! num_cmb_glob
!=====================
  
! Form local arrays depending on how many receivers each processor has
  num_cmb=icmb
  allocate(cmbfile_el(1:num_cmb,1:3),loc2globcmb(1:num_cmb))
  allocate(cmbfile_th(1:num_cmb))
  allocate(cmbfile_r(1:num_cmb))
  cmbfile_el(1:num_cmb,1:3)=cmbfile_el_loc(1:num_cmb,1:3)
  loc2globcmb(1:num_cmb)=loc2globcmb_loc(1:num_cmb)
  cmbfile_th(1:num_cmb)=cmbfile_th_loc(1:num_cmb)
  cmbfile_r(1:num_cmb)=cmbfile_r_loc(1:num_cmb)

! How many receivers does each processor have, do they sum to global number?
  if ( psum_int(num_cmb) /= num_cmb_glob ) then
     write(6,*)'PROBLEM: sum of local cmbcmbs is different than global!'
     if (lpr) write(6,*)'Global number of cmbrecs:',num_cmb_glob
     write(6,*)procstrg,'Number of cmbrecs:',num_cmb
     stop
  endif

  if (lpr) write(6,*)
  do icmb=0,nproc-1
     call barrier
     if (mynum==icmb) write(6,14)procstrg,num_cmb,num_cmb_glob
     call barrier
  enddo
  write(69,*)
  write(69,14)procstrg,num_cmb,num_cmb_glob
14 format('   ',a8,'has',i4,' out of',i6,' cmbrecs')

! Output colatitudes globally (this is the file needed to plot seismograms)
  if (lpr) then
     open(99997,file=datapath(1:lfdata)//'/cmbrec_pts.dat')
     do i=1,num_cmb_glob
        write(99997,*)cmbfile_th_glob(i),cmb2proc(i)
     enddo
     close(99997)
  endif

! Output colatitudes locally to infopath (this file is for info purposes!)
  maxcmblocerr=zero
  open(9998+mynum,file=infopath(1:lfinfo)//'/cmbrec_pts_'//appmynum//'.dat')
  write(9998+mynum,*)num_cmb
  do i=1,num_cmb ! Only over newly found local receiver locations

     write(9998+mynum,12)i,cmbfile_readth(loc2globcmb(i)),cmbfile_th(i), &
                           cmbfile_r(i), &
                           cmbfile_el(i,1),cmbfile_el(i,2),cmbfile_el(i,3)

     call define_io_appendix(appielem,loc2globcmb(i))

     if ( pi/180*router*abs(cmbfile_readth(loc2globcmb(i))-cmbfile_th(i) ) >  &
          min_distance_dim) then
        count_diff_loc=count_diff_loc+1
        write(6,22)procstrg, &
                  cmbfile_readth(loc2globcmb(i)),cmbfile_th(i),cmbfile_r(i)/1.d3

        if (dabs(cmbfile_readth(loc2globcmb(i))-cmbfile_th(i))> maxcmblocerr) &
             maxcmblocerr=dabs(cmbfile_readth(loc2globcmb(i))-cmbfile_th(i))/ &
                          180.*pi*router

     endif

22 format('   WARNING:',a8,' cmbrec location file/mesh:',3(f9.3))

     write(6,*)'  ',procstrg,'opening cmbrec file:',i,appielem

     open(200000+i,file=datapath(1:lfdata)//'/cmbfile_seis.dat'//appielem)
     open(250000+i,file=datapath(1:lfdata)//'/cmbfile_straintrace.dat'//appielem)

  enddo
  close(9998+mynum)

  write(69,15)count_diff_loc,num_cmb
  write(69,*)'  Maximal cmbrec location error [m]:',maxcmblocerr
  write(69,*)
15 format(i4,' out of',i4,' cmbrecs are located at wrong points.')

  maxcmblocerr=pmax(maxcmblocerr)
  if (lpr) then 
     write(6,*)
     write(6,*)'  maximal cmbrec location error [m]:',maxcmblocerr
     write(6,*)
  endif

12 format(i3,3(1pe12.4),i8,2(i2))

  deallocate(cmbfile_readth)
  deallocate(cmbfile_th_glob,cmbfile_th,cmbfile_th_loc)
  deallocate(cmbfile_r_glob,cmbfile_r,cmbfile_r_loc)
  deallocate(cmbfile_el_loc,loc2globcmb_loc,loc2globcmb,cmb2proc)

end subroutine prepare_from_recfile_cmb
!=============================================================================


!-----------------------------------------------------------------------------
subroutine open_hyp_epi_equ_anti
!
! Open files for generic checks: hypocenter,epicenter,equator,antipode.
! File output names: seis<LOCATION>{1,2,3}.dat
!                    where 1=s-component, 2=phi-component, 3=z-component.
! At hypocenter,epicenter,antipode, s-comp = longitudinal comp.
!                                   z-comp = radial/vertical comp.
! At equator                        s-comp = radial/vertical comp.
!                                   z-comp = longitudinal comp.
! Not writing the phi/transverse component for monopole sources.
!   
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_source, ONLY : have_src,src_type

  if (maxind>0) then

    if (have_src) then
       open(10001,file=datapath(1:lfdata)//'/seishypocenter1.dat') 
       if (src_type(1)/='monopole') &
            open(10002,file=datapath(1:lfdata)//'/seishypocenter2.dat') 
       open(10003,file=datapath(1:lfdata)//'/seishypocenter3.dat') 
    endif

    if (have_epi) then
       open(900,file=datapath(1:lfdata)//'/seisepicenter1.dat') 
       if (src_type(1)/='monopole') &
            open(903,file=datapath(1:lfdata)//'/seisepicenter2.dat') 
       open(906,file=datapath(1:lfdata)//'/seisepicenter3.dat') 
    endif

    if (have_equ) then
       open(902,file=datapath(1:lfdata)//'/seisequator1.dat')
       if (src_type(1)/='monopole') &
            open(904,file=datapath(1:lfdata)//'/seisequator2.dat')
       open(907,file=datapath(1:lfdata)//'/seisequator3.dat')
    endif

    if (have_antipode) then
       open(901,file=datapath(1:lfdata)//'/seisantipode1.dat') 
       if (src_type(1)/='monopole') &
            open(905,file=datapath(1:lfdata)//'/seisantipode2.dat')
       open(908,file=datapath(1:lfdata)//'/seisantipode3.dat')
    endif

 endif ! if maxind > 0

end subroutine open_hyp_epi_equ_anti
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_hyp_epi_equ_anti(t,disp)
!
! Generic seismograms for quick checks: hypocenter,epicenter,equator,antipode.
! Not writing the transverse component for monopole sources.
! See open_hyp_epi_equ_anti for component explanation.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_source, ONLY: iel_src,ipol_src,jpol_src,have_src,src_type
double precision :: t
real(kind=realkind), intent(in) :: disp(0:npol,0:npol,nel_solid,3)

  if (maxind>0) then
     if (mynum==0) then
        if (ipol_src /= 0 ) then
           write(6,*)'PROBLEM in hypocenter location!'
           write(6,*)'ipol  is not equal to zero, hence off the axis!',ipol_src
           stop
        endif
     endif

! hypocenter
    if (have_src) then
       if (src_type(1)=='dipole') then 
          write(10001,*)t,disp(ipol_src,jpol_src,iel_src,1)+&
                        disp(ipol_src,jpol_src,iel_src,2) ! s
          write(10002,*)t,disp(ipol_src,jpol_src,iel_src,1)-&
                        disp(ipol_src,jpol_src,iel_src,2) ! phi
       else
          write(10001,*)t,disp(ipol_src,jpol_src,iel_src,1) ! s
          if (src_type(1)=='quadpole') &
               write(10002,*)t,disp(ipol_src,jpol_src,iel_src,2) ! phi
        endif
        write(10003,*)t,disp(ipol_src,jpol_src,iel_src,3)  ! z
     endif
   
! epicenter
     if (have_epi) then
        if (src_type(1)=='dipole') then 
           write(900,*)t,disp(0,npol,ielepi,1)+disp(0,npol,ielepi,2) ! s
           write(903,*)t,disp(0,npol,ielepi,1)-disp(0,npol,ielepi,2) ! phi
        else
           write(900,*)t,disp(0,npol,ielepi,1) ! s
           if (src_type(1)=='quadpole') &
                write(903,*)t,disp(0,npol,ielepi,2) ! phi
        endif
        write(906,*)t,disp(0,npol,ielepi,3)  ! z
     endif

! antipode
    if (have_antipode) then
       if (src_type(1)=='dipole') then 
          write(901,*)t,disp(0,0,ielantipode,1)+disp(0,0,ielantipode,2) ! s
          write(905,*)t,disp(0,0,ielantipode,1)-disp(0,0,ielantipode,2) ! phi
       else
          write(901,*)t,disp(0,0,ielantipode,1) ! s
          if (src_type(1)=='quadpole') &
               write(905,*)t,disp(0,0,ielantipode,2) ! phi
        endif
        write(908,*)t,disp(0,0,ielantipode,3)  ! z
     endif

! equator
    if (have_equ) then
       if (src_type(1)=='dipole') then 
          write(902,*)t,disp(npol,npol,ielequ,1)+disp(npol,npol,ielequ,2) ! s
          write(904,*)t,disp(npol,npol,ielequ,1)-disp(npol,npol,ielequ,2) ! phi
       else
          write(902,*)t,disp(npol,npol,ielequ,1) ! s
          if (src_type(1)=='quadpole') &
               write(904,*)t,disp(npol,npol,ielequ,2) ! phi
        endif
        write(907,*)t,disp(npol,npol,ielequ,3)  ! z
     endif

  endif ! if maxind > 0

end subroutine compute_hyp_epi_equ_anti
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_recfile_seis_bare(disp)

use data_source, ONLY : src_type

include "mesh_params.h"

real(kind=realkind), intent(in) :: disp(0:npol,0:npol,nel_solid,3)
!real(kind=realkind), intent(in) :: velo(0:npol,0:npol,nel_solid,3)

integer :: i

     if (src_type(1)=='monopole') then
        do i=1,num_rec
! order: u_s, u_z
           write(100000+i,*) &
               disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1), &
               disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
        enddo

     elseif (src_type(1)=='dipole') then
        do i=1,num_rec
! order: u_s, u_phi,u_z
! 
           write(100000+i,*) &
                disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) + &
                disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2),  &
                disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) - &
                disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2), &
                disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
        enddo

     elseif (src_type(1)=='quadpole') then
        do i=1,num_rec
! order: u_s, u_phi,u_z 
           write(100000+i,*) &
               disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1), &
               disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2), &
               disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
        enddo

endif !src_type(1)

end subroutine compute_recfile_seis_bare
!=============================================================================
!! Calculate displacement at receiver locations and pass to nc_dump_rec

subroutine nc_compute_recfile_seis_bare(disp) 

use data_source, ONLY : src_type
use nc_routines, ONLY : nc_dump_rec
implicit none
include "mesh_params.h"
real(kind=realkind), intent(in)  :: disp(0:npol,0:npol,nel_solid,3)
real(kind=realkind)              :: disp_rec(3,num_rec)
character(len=50)                :: filename
integer                          :: i

!allocate(disp_surf(3,num_rec))

if (src_type(1)=='monopole') then
! order: u_s, u_z
    do i=1,num_rec
         disp_rec(1,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1))
         disp_rec(2,i)= 0.0 
         disp_rec(3,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3))  
    enddo

elseif (src_type(1)=='dipole') then
    
    do i=1,num_rec
         disp_rec(1,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) + &
                            disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2))
         disp_rec(2,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) - &
                            disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2))
         disp_rec(3,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3))
    enddo

elseif (src_type(1)=='quadpole') then

    do i=1,num_rec
         disp_rec(1,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1))
         disp_rec(2,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2))
         disp_rec(3,i)=real(disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3))  
    enddo

end if !src_type(1)

call nc_dump_rec(disp_rec) 

!deallocate(disp_surf)

end subroutine nc_compute_recfile_seis_bare

!-----------------------------------------------------------------------------
!af
!!$subroutine compute_recfile_seis_binary(time,disp,velo)
!!$
!!$use data_source, ONLY : src_type
!!$!use write_binary
!!$
!!$include "mesh_params.h"
!!$
!!$real(kind=realkind), intent(in) :: time
!!$real(kind=realkind), intent(in) :: disp(0:npol,0:npol,nel_solid,3)
!!$real(kind=realkind), intent(in) :: velo(0:npol,0:npol,nel_solid,3)
!!$
!!$integer :: i
!!$
!!$do i=1,num_rec
!!$  call dump_trace_3comp_dble(time,&
!!$                        disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1),&
!!$                        disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2),&
!!$                        disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3),&
!!$                        fname_rec_seis(i))
!!$  call dump_trace_3comp_dble(time,&
!!$                        velo(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1),&
!!$                        velo(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2),&
!!$                        velo(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3),&
!!$                        fname_rec_velo(i))
!!$enddo
!!$
!!$end subroutine compute_recfile_seis_binary
!!$!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_recfile_cmb(velo,lap_sol)

use data_source, ONLY : src_type

include "mesh_params.h"

real(kind=realkind), intent(in) :: velo(0:npol,0:npol,nel_solid,3)
real(kind=realkind)             :: lap_sol(0:npol,0:npol,nel_solid,1)
integer :: i

  if (src_type(1)=='monopole') then
  do i=1,num_cmb
     write(200000+i,*)velo(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),1),&
                      velo(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),3)

     write(250000+i,*)lap_sol(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),1) 
  enddo

  else
     do i=1,num_cmb
     write(200000+i,*)velo(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),1),&
                      velo(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),2),&
                      velo(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),3)  

     write(250000+i,*)lap_sol(cmbfile_el(i,2),cmbfile_el(i,3),cmbfile_el(i,1),1) 

     enddo
  endif

end subroutine compute_recfile_cmb
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_surfelem(disp,velo)
!
! Save one displacement and velocity trace for each element on the surface 
! which are both needed for kernels (du and v0 inside the cross-correlation)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
use data_io,     ONLY : istrain
use data_source, ONLY : src_type
use nc_routines, ONLY : nc_dump_surface
include "mesh_params.h"

real(kind=realkind), intent(in) :: disp(0:npol,0:npol,nel_solid,3)
real(kind=realkind), intent(in) :: velo(0:npol,0:npol,nel_solid,3)
real                            :: dumpvar(maxind,3)
integer                         :: i

  dumpvar = 0.0

  if (use_netcdf) then
      if (src_type(1)=='monopole') then
          do i=1,maxind
              dumpvar(i,1) = real(disp(npol/2,jsurfel(i),surfelem(i),1))
              dumpvar(i,2) = real(disp(npol/2,jsurfel(i),surfelem(i),3))
          enddo
      else
          do i=1,maxind
              dumpvar(i,:) = real(disp(npol/2,jsurfel(i),surfelem(i),:))
          end do
      end if !monopole
      call nc_dump_surface(dumpvar(:,:), 'disp')
      
      if (src_type(1)=='monopole') then
          do i=1,maxind
              dumpvar(i,1) = real(velo(npol/2,jsurfel(i),surfelem(i),1))
              dumpvar(i,2) = real(velo(npol/2,jsurfel(i),surfelem(i),3))
          end do
      else
          do i=1,maxind    
              dumpvar(i,:) = real(velo(npol/2,jsurfel(i),surfelem(i),:))
          end do
      end if !monopole
      call nc_dump_surface(dumpvar(:,:), 'velo')

  else !use_netcdf
      if (src_type(1)=='monopole') then
      do i=1,maxind
         write(40000000+i,*)real(disp(npol/2,jsurfel(i),surfelem(i),1)),&
                            real(disp(npol/2,jsurfel(i),surfelem(i),3))
         write(50000000+i,*)real(velo(npol/2,jsurfel(i),surfelem(i),1)),&
                            real(velo(npol/2,jsurfel(i),surfelem(i),3))
      enddo

      else
         do i=1,maxind
         write(40000000+i,*)real(disp(npol/2,jsurfel(i),surfelem(i),1)),&
                            real(disp(npol/2,jsurfel(i),surfelem(i),2)),&
                            real(disp(npol/2,jsurfel(i),surfelem(i),3))  
         write(50000000+i,*)real(velo(npol/2,jsurfel(i),surfelem(i),1)),&
                            real(velo(npol/2,jsurfel(i),surfelem(i),2)),&
                            real(velo(npol/2,jsurfel(i),surfelem(i),3))
         enddo
      end if !monopole
  end if !netcdf

end subroutine compute_surfelem
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_surfelem_strain(u)
!
! Save one displacement and velocity trace for each element on the surface 
! which are both needed for kernels (du and v0 inside the cross-correlation)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY: inv_rho_fluid,inv_s_rho_fluid,usz_fluid, inv_s_solid!,deviator
use data_source, ONLY: src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,axisym_laplacian_fluid_add
use pointwise_derivatives, ONLY: axisym_laplacian_solid,axisym_laplacian_solid_add
use unit_stride_colloc, ONLY: collocate0_1d_existent
use nc_routines, ONLY   : nc_dump_surface

include "mesh_params.h"
real(kind=realkind), intent(in) :: u(0:npol,0:npol,nel_solid,3)

real(kind=realkind)             :: lap_sol(0:npol,0:npol,nel_solid,2)
real(kind=realkind)             :: lap_flu(0:npol,0:npol,nel_fluid,2)
real(kind=realkind)             :: dumpvar(maxind, 6)
real(kind=realkind)             :: strain(0:npol,nel_solid,6)
character(len=5)                :: appisnap



integer :: i,jj,j

strain=0.

  if (src_type(1)=='dipole') then
    call axisym_laplacian_solid(u(:,:,:,1)+u(:,:,:,2),lap_sol)
  else
    call axisym_laplacian_solid(u(:,:,:,1),lap_sol) ! 1: dsus, 2: dzus
  endif
  strain(:,:,1)  = lap_sol(npol/2,:,:,1) ! ds us

  call axisym_laplacian_solid_add(u(:,:,:,3),lap_sol) ! 1:dsuz+dzus,2:dzuz+dsus

! calculate entire E31 term: (dsuz+dzus)/2
  strain(:,:,4) = lap_sol(npol/2,:,:,1)* real(.5,kind=realkind) ! ds uz

! Components involving phi....................................................
   if (src_type(1)=='monopole') then
      strain(:,:,2) = inv_s_solid(npol/2,:,:)* u(npol/2,:,:,1) ! dp up
      strain(:,:,3) = lap_sol(npol/2,:,:,2) - strain(:,:,1) ! dz uz

! __________________________________________________________________________

   elseif (src_type(1)=='dipole') then 
      strain(:,:,2) = real(2.,kind=realkind)*inv_s_solid(npol/2,:,:)* u(npol/2,:,:,1) ! dp up
      strain(:,:,3) = lap_sol(npol/2,:,:,2) - strain(:,:,1) ! dz uz

      call axisym_laplacian_solid(u(:,:,:,1)-u(:,:,:,2),lap_sol) !1:dsup,2:dzup
      strain(:,:,5) = ( inv_s_solid(npol/2,:,:) *u(npol/2,:,:,2)  + &
                              real(.5,kind=realkind) * lap_sol(npol/2,:,:,1) ) ! ds up

      strain(:,:,6) = real(.5,kind=realkind) *( inv_s_solid(npol/2,:,:) *u(npol/2,:,:,3)  + &
                               lap_sol(npol/2,:,:,2) ) ! dz up

! __________________________________________________________________________
   elseif (src_type(1)=='quadpole') then

      strain(:,:,2) = inv_s_solid(npol/2,:,:)* &
                             ( u(npol/2,:,:,1)-real(2.,kind=realkind)*u(npol/2,:,:,2)) ! dp up
      strain(:,:,3) = lap_sol(npol/2,:,:,2) - strain(:,:,1) ! dz uz

      call axisym_laplacian_solid(u(:,:,:,2),lap_sol) ! 1: dsup, 2: dzup

      strain(:,:,5) = real(.5,kind=realkind) * ( inv_s_solid(npol/2,:,:) * &
                              (real(2.,kind=realkind)* u(npol/2,:,:,1) - u(npol/2,:,:,2))+ &
                              lap_sol(npol/2,:,:,1) ) ! ds up

      strain(:,:,6) = ( inv_s_solid(npol/2,:,:) *u(npol/2,:,:,3)  + &
                               real(.5,kind=realkind) *lap_sol(npol/2,:,:,2) ) ! dz up

   endif
  if (use_netcdf) then
      do i=1,maxind
        dumpvar(i,:) = real(strain(j,surfelem(i),1:6))
      enddo
      call nc_dump_surface(dumpvar(:,1:6), 'stra')
      do i=1,maxind
        dumpvar(i,1:3) = real(u(npol/2,j,surfelem(i),1:3))
      enddo
      call nc_dump_surface(dumpvar(:,1:3), 'srcd')
  end if
  if (.not. use_netcdf) then

      do i=1,maxind
         do j=0,npol
            write(60000000+i,20)(real(strain(j,surfelem(i),jj)),jj=1,6)
            write(70000000+i,30)(real(u(npol/2,j,surfelem(i),jj)),jj=1,3)
         enddo
      enddo

20 format(6(1pe11.3))
30 format(3(1pe11.3))
  end if

end subroutine compute_surfelem_strain
!=============================================================================



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!========================
end module seismograms
!========================
 
