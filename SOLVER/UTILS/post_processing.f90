!-----------------------------------------------------------------------------
module data_all

  implicit none
  public

  logical,allocatable, dimension(:) :: correct_azi, rot_rec_post_array
  logical                           :: rot_src_rec_true
  character(len=100)                :: stf_file, rot_rec_file, filename, fname

  character(len=100),allocatable, dimension(:)      :: bkgrndmodel, rot_rec
  character(len=7),allocatable, dimension(:)        :: stf_type
  character(len=100),allocatable, dimension(:)      :: recname
  character(len=100),allocatable, dimension(:,:)    :: outname, outname2, rec_full_name
  character(len=10),allocatable, dimension(:,:)     :: src_type
  integer, allocatable, dimension(:)                :: nt, nrec, nt_seis, nt_strain, &
                                                       nt_snap, ibeg, iend
  character(len=1),dimension(3)       :: reccomp
  character(len=100), allocatable     :: simdir(:)
  character(len=4)                    :: appmynum, appidur
  integer                             :: i, it, nsim, isim, iseis, ii, iii, ishift
  real                                :: junk, rec_loc_tol, Mij(6), tshift
  real, dimension(3)                  :: rloc_xyz, rloc_xyz_tmp, rloc_rtp
  real, dimension(:), allocatable     :: time, src_depth
  real, dimension(:,:), allocatable   :: seis, seis_fil, seis_sglcomp
  real, dimension(:,:,:), allocatable :: mij_phi
  real, allocatable, dimension(:)     :: dt, period, dt_seis, dt_strain, dt_snap ,magnitude
  real, allocatable, dimension(:)     :: thr_orig,phr_orig
  real, dimension(3,3)                :: rot_mat
  character(len=3)                    :: rec_comp_sys
  logical                             :: rot_rec_post,any_sum_seis_true,load_snaps
  character(len=12)                   :: src_file_type

  ! discrete dirac sources
  real, allocatable, dimension(:)     :: shift_fact
  integer, allocatable, dimension(:)  :: ishift_deltat, ishift_seisdt, ishift_straindt

  ! parameters from param_post_processing
  logical, dimension(:), allocatable  :: rot_rec_posttmp, sum_seis_true, load_snapstmp, &
                                         negative_time
  character(len=3), allocatable       :: rec_comp_systmp(:)
  real, dimension(:,:), allocatable   :: mij_loc
  real, dimension(:), allocatable     :: conv_period, srccolat, srclon
  character(len=7), allocatable       :: conv_stf(:)
  character(len=100), allocatable     :: outdir(:)
  character(len=4), allocatable       :: seistype(:)

  real, allocatable :: period_final(:)
  real, allocatable :: trans_rot_mat(:,:,:)
  real, allocatable :: colat(:,:),lon(:,:)

end module data_all
!=============================================================================

!-----------------------------------------------------------------------------
module global_par

  double precision, parameter :: pi = 3.1415926535898

  double precision, parameter :: zero = 0d0, half = 5d-1, third = 1d0/3d0
  double precision, parameter :: quart = 25d-2, one = 1d0, sixth = 1d0/6d0
  double precision, parameter :: two = 2d0, three=3d0, four=4d0, five=5.d0
  double precision, parameter :: fifth = 2.d-1
  double precision, parameter :: epsi = 1d-30
  real, parameter             :: epsi_real=1.e-10
  real(kind=8), parameter     :: smallval_dble = 1e-11
  real, parameter             :: decay = 3.5d0
  real, parameter             :: shift_fact1 = 1.5d0
end module global_par
!=============================================================================

!-----------------------------------------------------------------------------
program post_processing_seis

  use data_all
  use global_par
  implicit none 
  double precision :: arg1

  call read_input

  if (rot_rec_post ) then
     rot_rec_post_array = .true.
     if (rec_comp_sys == 'sph') then
        reccomp(1) = 'th'
        reccomp(2) = 'ph'
        reccomp(3) = 'r'
     elseif (rec_comp_sys == 'enz') then
        reccomp(1) = 'N'
        reccomp(2) = 'E'
        reccomp(3) = 'Z'
     elseif (rec_comp_sys == 'cyl') then
        reccomp(1) = 's'
        reccomp(2) = 'ph'
        reccomp(3) = 'z'
     elseif (rec_comp_sys == 'xyz') then
        reccomp(1) = 'x'
        reccomp(2) = 'y'
        reccomp(3) = 'z'
     elseif (rec_comp_sys == 'src') then
        reccomp(1) = 'R'
        reccomp(2) = 'T'
        reccomp(3) = 'Z'
     endif

  else
     rot_rec_post_array = .false.
     reccomp(1) = 's0'
     reccomp(2) = 'ph0'
     reccomp(3) = 'z0' 
     write(6,*) 'WARNING: Receiver components are cylindrical and rotated with the source at the north pole!'
  endif
  write(6,*)'receiver components: ',(reccomp(i),i=1,3)
  write(6,*)
  
  ! input seismogram names
  allocate(recname(nrec(1)), thr_orig(nrec(1)), phr_orig(nrec(1)))
  
  ! these are the original receiver coordinates, having the source at the actual location
  open(unit=61,file=trim(simdir(1))//'/Data/receiver_names.dat')
  do i=1,nrec(1)
     read(61,*)recname(i),thr_orig(i),phr_orig(i)
     if (phr_orig(i)==360.) phr_orig(i)=0.
  enddo
  close(61)
  thr_orig=thr_orig/180.*pi 
  phr_orig=phr_orig/180.*pi

  ! these are the rotated receiver coordinates, having the source at the north pole
  allocate(colat(nrec(1),nsim),lon(nrec(1),nsim))
  do isim=1,nsim
     open(unit=20,file=trim(simdir(isim))//'/Data/receiver_pts.dat')
     open(unit=21,file=trim(outdir(1))//'/receiver_gll_locs.dat')
     write(21,'(4a15)') ' ', 'colat', 'lat', 'lon'
     do i=1,nrec(isim)
        read(20,*) colat(i,isim), lon(i,isim), junk
        if (abs(lon(i,isim)-360.)<0.01) lon(i,isim)=0.
       
        ! transform to xyz
        rloc_xyz(1) = sin(colat(i,isim) * pi / 180.) * cos(lon(i,isim) * pi / 180.) 
        rloc_xyz(2) = sin(colat(i,isim) * pi / 180.) * sin(lon(i,isim) * pi / 180.) 
        rloc_xyz(3) = cos(colat(i,isim) * pi / 180.) 
        
        ! Rotate to the original (i.e. real src-rec coordinate-based) u_xyz
        rloc_xyz_tmp = rloc_xyz
        rloc_xyz(1) =   cos(srccolat(1)) * cos(srclon(1)) * rloc_xyz_tmp(1) &
                    & - sin(srclon(1)) * rloc_xyz_tmp(2) &    
                    & + sin(srccolat(1)) * cos(srclon(1)) * rloc_xyz_tmp(3) 
        rloc_xyz(2) =   cos(srccolat(1)) * sin(srclon(1)) * rloc_xyz_tmp(1) &
                    & + cos(srclon(1)) * rloc_xyz_tmp(2) &
                    & + sin(srccolat(1)) * sin(srclon(1)) * rloc_xyz_tmp(3)
        rloc_xyz(3) = - sin(srccolat(1)) * rloc_xyz_tmp(1) &
                    & + cos(srccolat(1)) * rloc_xyz_tmp(3)

        if (rloc_xyz(1) > 1.) rloc_xyz(1) = 1.
        if (rloc_xyz(1) < -1.) rloc_xyz(1) = -1.
        if (rloc_xyz(2) > 1.) rloc_xyz(2) = 1.
        if (rloc_xyz(2) < -1.) rloc_xyz(2) = -1.
        if (rloc_xyz(3) > 1.) rloc_xyz(3) = 1.
        if (rloc_xyz(3) < -1.) rloc_xyz(3) = -1.
        
        rloc_xyz = rloc_xyz / (rloc_xyz(1)**2 + rloc_xyz(2)**2 + rloc_xyz(3)**2)**.5

        ! compute colat and lon
        rloc_rtp(2) = acos(rloc_xyz(3))

        arg1 = (rloc_xyz(1)  + smallval_dble) / &
               ((rloc_xyz(1)**2 + rloc_xyz(2)**2)**.5 + smallval_dble)
        
        if (arg1 > 1.) arg1 = 1.
        if (arg1 < -1.) arg1 = -1.

        if (rloc_xyz(2) >= 0.) then
           rloc_rtp(3) = acos(arg1)
        else
           rloc_rtp(3) = 2.*pi - acos(arg1)
        end if

        thr_orig(i) = rloc_rtp(2)
        phr_orig(i) = rloc_rtp(3)
        
        ! write exact receiver locations (gll points) in the earth fixed coordinate system to file
        write(21,'(a15,3f15.8)') recname(i), rloc_rtp(2) / pi * 180., &
                                 rloc_rtp(2) / pi * 180. - 90., &
                                 rloc_rtp(3) / pi * 180.
     enddo
     close(20)
     close(21)
  enddo
  if (maxval(lon) > 0.1) any_sum_seis_true = .true. 
  colat = colat / 180. * pi
  lon = lon / 180. * pi
  
  if (any_sum_seis_true .or. rot_rec_post) then
     ! calculate moment tensor and azimuth prefactors/radiation patterns for each simulation
     allocate(mij_phi(nrec(1),nsim,3)); mij_phi = -1.E30
     call compute_radiation_prefactor(mij_phi,nrec(1),nsim,lon)
  endif !sum_seis

  ! convolve with source period?
  allocate(seis_fil(nt_seis(1),3)) ! need this for saving...
  allocate(period_final(nsim))
  
  do isim=1,nsim
     if (conv_period(isim)>0. ) then 
        if (stf_type(isim) /= 'dirac_0' .and. stf_type(isim) /= 'quheavi') then 
           write(6,*) ' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !'
           write(6,*)'    WARNING! You want to convolve although seismograms are upon a ', &
                    trim(stf_type(isim))
           write(6,*) ' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !'
        endif

        if (conv_period(isim) < period(isim)) then
           write(6,*)'Cannot convolve with a period shorter than allowed by the mesh/simulation!'
           write(6,*)'convolution period, min. mesh period [s]:',conv_period(isim),period(isim)
           write(6,*)'Using mesh period instead......'
           conv_period(isim)=period(isim)
        endif

        write(6,*)' Convolve with source period [s]:',conv_period(isim)
        period_final(isim) = conv_period(isim)
     else
        period_final(isim) = period(isim)
     endif
     period(isim) = period_final(isim)
  enddo

  ! define time series
  allocate(time(nt_seis(1)), seis(nt_seis(1),3))
  do iseis=1, nt_seis(1)
     time(iseis) = real(iseis) * dt_seis(1)
  enddo
  allocate(seis_sglcomp(nt_seis(1),3))

  ! output seismogram names
  allocate(outname(nrec(1),nsim), outname2(nrec(1),nsim), rec_full_name(nrec(1),nsim))
  do i=1, nrec(1)
     do isim=1,nsim
        rec_full_name(i,isim) = trim(recname(i))//'_'//seistype(isim)//'_post'
        if (sum_seis_true(isim)) rec_full_name(i,isim) = trim(rec_full_name(i,isim))//'_mij'

        call define_io_appendix(appidur, int(conv_period(isim)))
        rec_full_name(i,isim) = trim(rec_full_name(i,isim))//'_conv'//appidur//''
        outname(i,isim) = trim(outdir(isim))//'/SEISMOGRAMS/'//trim(rec_full_name(i,isim))

        ! if no summation, put output traces in each of the simulation directories
        outname2(i,isim) = trim(outdir(isim))//'/SEISMOGRAMS/UNPROCESSED/'&
                            //trim(recname(i))//'_'//seistype(isim)//'_post_'&
                            //trim(src_type(isim,2))
        
        write(6,*) 'outname: ', i, trim(outname(i,isim))
        write(6,*) 'outname2: ', isim, trim(outname2(i,isim))
     enddo
  enddo

  if (conv_period(1) > 0.)  then 
    tshift = shift_fact1 * conv_period(1)
  else
     tshift = 0.
  endif
  ishift = int((shift_fact(1)+tshift)/(time(2)-time(1)))
  write(6,*) 'ishift1:', shift_fact(1), tshift,time(2)-time(1)
  write(6,*) 'ishift:', (shift_fact(1)+tshift)/(time(2)-time(1)),int((shift_fact(1)+tshift)/(time(2)-time(1)))

  ! ----------------- Start actual post processing -------------------------
  !=================
  do i=1,nrec(1)
  !=================
     call define_io_appendix(appmynum,i)
     seis = 0.
     seis_sglcomp = 0.
     seis_fil = 0.

     !::::::::::::::::::::::::::::::::
     do isim=1,nsim
     !::::::::::::::::::::::::::::::::

        ! load seismograms from all directories
        open(unit=60,file = trim(simdir(isim))//'/Data/'//trim(recname(i))//'_'&
                            //seistype(isim)//'.dat')
        write(6,*) 'opened ', trim(simdir(isim))//'/Data/'//trim(recname(i))//'_'&
                    //seistype(isim)//'.dat'

        if (src_type(isim,1) == 'monopole') then 
           do iseis=1, nt_seis(1) 
              read(60,*) seis_sglcomp(iseis,1), seis_sglcomp(iseis,3)
           enddo
        else 
           do iseis=1, nt_seis(1)
              read(60,*) seis_sglcomp(iseis,1), seis_sglcomp(iseis,2), seis_sglcomp(iseis,3)
           enddo
        endif
        close(60)

        if (nsim>1 .or. .not. sum_seis_true(isim)) then 
        ! output seismograms into each simulation's  data directory. 
        ! In original rotated cylindrical frame, without azimuthal radiation factor
        open(unit=150,file=trim(outname2(i,isim))//'_s0.dat',status='new')
        open(unit=151,file=trim(outname2(i,isim))//'_ph0.dat',status='new')
        open(unit=152,file=trim(outname2(i,isim))//'_z0.dat',status='new')

        if (negative_time(isim)) then
           do it=1,nt_seis(1)
              write(150,*)time(it)-shift_fact(1)-tshift,seis_sglcomp(it,1)
              write(151,*)time(it)-shift_fact(1)-tshift,seis_sglcomp(it,2)
              write(152,*)time(it)-shift_fact(1)-tshift,seis_sglcomp(it,3)
           enddo
        else
           do it=ishift+1,nt_seis(1)
              write(150,*)time(it-ishift-1),seis_sglcomp(it,1)
              write(151,*)time(it-ishift-1),seis_sglcomp(it,2)
              write(152,*)time(it-ishift-1),seis_sglcomp(it,3)
           enddo
        endif
           close(150);close(151);close(152)
        endif 

        ! sum seismograms upon separate moment tensors and include azimuthal radiation factor 
        ! This needs to be done prior to component rotation such that the system is still &
        ! cylindrical.
        if (sum_seis_true(isim)) then
           call sum_individual_wavefields(seis, seis_sglcomp, nt_seis(1), mij_phi(i,isim,:))
        else
           seis = seis_sglcomp
        endif

     !::::::::::::::::::::::::::::::::
     enddo !isim
     !::::::::::::::::::::::::::::::::

     ! convolve with a source time function
     if (conv_period(1)>0.)  then 
        write(6,*)'conv period:',conv_stf(1)
        call convolve_with_stf(conv_period(1),dt_seis(1),nt_seis(1),src_type(1,1),conv_stf(1),&
                                               outdir(1),seis,seis_fil)
     else
        seis_fil=seis
     endif

     ! rotate receiver components to coordinate system of choice
     ! doing this after the nsim loop assumes that all simulations have the same coordinates. 
     ! this is fine within a Mij set, but NOT for a finite fault: In that case, we will have to 
     ! amend this by another loop over fault points. 
     if (rot_rec_post_array(1)) &
            call rotate_receiver_comp(1, rec_comp_sys, srccolat(1), srclon(1), &
                                      colat(i,1), lon(i,1), thr_orig(i), phr_orig(i), &
                                      nt_seis(1), seis_fil)

     ! write processed seismograms into joint outdir
     open(unit=50,file=trim(outname(i,1))//'_'//reccomp(1)//'.dat',status='new')
     open(unit=51,file=trim(outname(i,1))//'_'//reccomp(2)//'.dat',status='new')
     open(unit=52,file=trim(outname(i,1))//'_'//reccomp(3)//'.dat',status='new')
     if (negative_time(1)) then
        write(6,*)' writing seismograms into joint directory: negative time'
        do it=1,nt_seis(1)
           write(50,*)time(it)-shift_fact(1)-tshift,seis_fil(it,1)
           write(51,*)time(it)-shift_fact(1)-tshift,seis_fil(it,2)
           write(52,*)time(it)-shift_fact(1)-tshift,seis_fil(it,3)
        enddo
     else
        write(6,*)' writing seismograms into joint directory: zero time'
        do it=ishift+1,nt_seis(1)
           write(50,*)time(it-ishift-1),seis_fil(it,1)
           write(51,*)time(it-ishift-1),seis_fil(it,2)
           write(52,*)time(it-ishift-1),seis_fil(it,3)
        enddo
     endif
     close(50);close(51);close(52);

  !=================
  enddo ! nrec
  !=================

  ! plot original source and receiver locations in google earth kml file with link 
  ! to seismogram images
  call save_google_earth_kml(srccolat(1), srclon(1), src_depth(1), Mij,period_final(1), &
                             thr_orig,phr_orig, reccomp, src_type(1,1), sum_seis_true(1), &
                             nsim,nrec(1), outdir(1), rec_full_name(:,1))
  
  write(6,*) 'writing matlab input files for record sections...'
  
  open(unit=99,file=trim(outdir(1))//'/info_matlab.dat')
  write(99,*)nrec(1)
  if (negative_time(1)) then 
     write(99,*)nt_seis(1)
  else
     write(99,*)nt_seis(1)-ishift
  endif
  
  write(99,*)time(2)-time(1)
  write(99,*)conv_period(1)
  do i=1,nrec(1)
     write(99,*)colat(i,1)
  enddo
  close(99)
  
  open(unit=99,file=trim(outdir(1))//'/info_seisnames.dat')
  open(unit=98,file=trim(outdir(1))//'/info_seisstations.dat')
  do ii=1,3
     do i=1,nrec(1)
        fname=trim(rec_full_name(i,1))//'_'//reccomp(ii)//".dat"
        write(99,*)trim(fname)
        write(98,*)trim(recname(i))
     enddo
  enddo
  close(99)
  close(98)
  
  open(unit=99,file=trim(outdir(1))//'/info_seiscomp.dat')
  do i=1,3
       write(99,*)reccomp(i)
  enddo
  close(99)
    
  ! compute snapshots of the wavefield in a 3D sphere with given top/bottom surface and 
  ! opening cross section (see param_snaps)
  if (load_snaps) call compute_3d_wavefields

  write(6,*)' DONE with seismogram post processing!'
21 format(a100)
end program post_processing_seis
!=============================================================================

!-----------------------------------------------------------------------------
subroutine read_input

  use data_all
  use global_par
  implicit none 
  logical :: file_exist
  integer :: length

  ! read directories for multiple simulations
  open(unit=99,file='param_sum_seis', status='old',position='rewind')
  read(99,*)nsim
  if (nsim==1) write(6,*)'no need to sum, one simulation only!'
  allocate(simdir(nsim))
  do isim=1,nsim
     read(99,*)simdir(isim)
     write(6,*)'simulation dir:',trim(simdir(isim))
  enddo
  close(99)

  ! read in post processing input file
  allocate(rot_rec_posttmp(nsim))
  allocate(rec_comp_systmp(nsim))
  allocate(sum_seis_true(nsim))
  allocate(mij_loc(6,nsim))
  allocate(conv_period(nsim))

  allocate(conv_stf(nsim))
  allocate(srccolat(nsim))
  allocate(srclon(nsim))
  allocate(load_snapstmp(nsim))
  allocate(outdir(nsim))
  allocate(seistype(nsim))
  allocate(negative_time(nsim))

  any_sum_seis_true = .false.
  do isim=1,nsim
     open(unit=99,file=trim(simdir(isim))//'/param_post_processing')
     read(99,*)rot_rec_posttmp(isim)
     read(99,*)rec_comp_systmp(isim)
     read(99,*)sum_seis_true(isim)
     do i=1,6
        read(99,*)mij_loc(i,isim)
     enddo
     read(99,*)conv_period(isim)
     read(99,*)conv_stf(isim)
     read(99,*)srccolat(isim)
     read(99,*)srclon(isim)
     read(99,*)load_snapstmp(isim)
     read(99,*)seistype(isim)
     read(99,*)outdir(isim)
     read(99,*)negative_time(isim)
     close(99)
   
     write(6,*)'seismogram type:',seistype(isim)
     write(6,*)'Output directory:',trim(outdir(isim))
     write(6,*)'shift time series to negative time (event at zero)?',negative_time(isim)
   
     ! read standard in, if any
     inquire(file=trim(outdir(1))//"/param_post_processing_overwrite",exist=file_exist)
     if (file_exist) then
        write(6,*)'OVERWRITING generic post processing parameters...'
        open(unit=99,file=trim(outdir(1))//'/param_post_processing_overwrite')
        read(99,*)length
        read(99,*)rec_comp_systmp(1)
        write(6,*)'overwrite receiver component system:',rec_comp_systmp(1)
        rec_comp_systmp(:)=rec_comp_systmp(1)
        if (length>1) then 
           read(99,*)conv_period(1)
           write(6,*)'overwrite period:',conv_period(1)
           conv_period(:)=conv_period(1)
        endif
        if (length>2) then 
           read(99,*)conv_stf(1)
           write(6,*)'overwrite convolution type:',conv_stf(1)
           conv_stf(:)=conv_stf(1)
        endif
        if (length>3) then 
           read(99,*)seistype(1)
           write(6,*)'overwrite seismogram type:',seistype(1)
           seistype(:)=seistype(1)
        endif
        close(99)
     endif
   
     if (load_snapstmp(isim)) load_snaps=.true.
   
     if (sum_seis_true(isim)) any_sum_seis_true=.true.
   
     if ( sum_seis_true(isim) .and. nsim ==2)  then
        write(6,'(a,/,a)') ' WARNING: You want to sum seismograms but only have two simulations', &
                           ' ==> are you sure?'
     elseif (.not. sum_seis_true(isim) .and. nsim == 4)  then
        write(6,*) '  WARNING: you have 4 simulations but do not want to sum?'
     endif
   
     write(6,*)'receiver system:',isim,rec_comp_systmp(isim)
     if (isim>1) then
        if (rec_comp_systmp(isim)/=rec_comp_systmp(isim-1) ) then 
            write(6,*) 'inconsistency with receiver component system:'
            write(6,*) simdir(isim),rec_comp_systmp(isim)
            write(6,*) simdir(isim-1),rec_comp_systmp(isim-1)
            write(6,*) 'make sure these are all identical in the subdirectories'
            stop
        end if
     endif
     rec_comp_sys = rec_comp_systmp(1)
   
     write(6,*) 'receiver rotation?',isim,rot_rec_posttmp(isim)
     if (isim>1) then
        if (rot_rec_posttmp(isim) .neqv. rot_rec_posttmp(isim-1))  then 
            write(6,*) 'inconsistency with receiver rotation:'
            write(6,*) simdir(isim),rot_rec_posttmp(isim)
            write(6,*) simdir(isim-1),rot_rec_posttmp(isim-1)
            write(6,*) 'make sure these are all identical in the subdirectories'
            stop
        endif
     endif
     rot_rec_post=rot_rec_posttmp(1)
   
     write(6,*)' Input from param_post_processing:',simdir(isim)
     write(6,*)'  Rotate receivers?',rot_rec_post
     if (rot_rec_post) write(6,*)'  component system: ',rec_comp_sys
     write(6,*)'  sum seismograms to full moment tensor?',sum_seis_true(isim)
     if (sum_seis_true(isim)) then 
        write(6,*)'  Moment tensor of individual run:'
        write(6,3)'','Mrr','Mtt','Mpp','Mtr','Mpr','Mtp'
        write(6,2)'',(Mij_loc(i,isim),i=1,6)
2       format(a2,6(1pe12.2))
3       format(a2,6(a12))
     endif
     write(6,*) '  convolution period:',conv_period(isim)
     write(6,*) '  convolution source time function: ' ,conv_stf(isim)
     write(6,*) '  source colatitude [deg]:',srccolat(isim)*180./pi
     write(6,*) '  source longitude:',srclon(isim)*180./pi
     write(6,*) '  load snaps?',load_snapstmp(isim)
     write(6,*) '  output directory: ',trim(outdir(isim))
     write(6,*)
  enddo 

  tshift = 0.

  allocate(bkgrndmodel(nsim), stf_type(nsim))
  allocate(src_type(nsim,2))
  allocate(dt(nsim), period(nsim), magnitude(nsim), dt_seis(nsim))
  allocate(dt_strain(nsim), dt_snap(nsim))
  allocate(correct_azi(nsim), rot_rec(nsim), rot_rec_post_array(nsim))
  allocate(nt(nsim), nrec(nsim), nt_seis(nsim), nt_strain(nsim), nt_snap(nsim))
  allocate(ibeg(nsim), iend(nsim), src_depth(nsim))
  allocate(shift_fact(nsim))
  allocate(ishift_deltat(nsim), ishift_seisdt(nsim), ishift_straindt(nsim))

  do isim = 1,nsim
     open(unit=99,file=trim(simdir(isim))//'/simulation.info')
     read(99,*) bkgrndmodel(isim)
     read(99,*) dt(isim)
     read(99,*) nt(isim)
     read(99,*) src_type(isim,1)
     read(99,*) src_type(isim,2)
     read(99,*) stf_type(isim)
     read(99,*) src_file_type
     read(99,*) period(isim)
     read(99,*) src_depth(isim)
     read(99,*) magnitude(isim)
     read(99,*) nrec(isim)
     read(99,*) nt_seis(isim)
     read(99,*) dt_seis(isim)
     !read(99,*) correct_azi(isim)
     correct_azi(isim) = .false.
     read(99,*) nt_strain(isim)
     read(99,*) dt_strain(isim)
     read(99,*) nt_snap(isim)
     read(99,*) dt_snap(isim)
     read(99,*) rot_rec(isim)
     read(99,*) ibeg(isim)
     read(99,*) iend(isim)
     read(99,*) shift_fact(isim)
     read(99,*) ishift_deltat(isim)
     read(99,*) ishift_seisdt(isim)
     read(99,*) ishift_straindt(isim)
     close(99)

     if (src_type(isim,2)=='vertforce' .or. src_type(isim,2)=='xforce' .or. &
         src_type(isim,2)=='yforce') sum_seis_true(isim)=.false.

     write(6,*) 'Simulations: ',isim,trim(simdir(isim))
     write(6,*) '  source type:',src_type(isim,1),' ',src_type(isim,2)
     write(6,*) '  source time function:',stf_type(isim)
     write(6,*) '  source file type:',src_file_type
     write(6,*) '  magnitude:',magnitude(isim)
     write(6,*) '  receivers:',nrec(isim)
     write(6,*) '  correct azi?',correct_azi(isim)
     write(6,*) '  source period:',period(isim)
     write(6,*) '  time steps:',nt(isim)
     write(6,*) '  rotate recs?',trim(rot_rec(isim))
     write(6,*) '  shift factor of stf:',shift_fact(isim)
     write(6,*) '  shift factor in samples (dt,dtseis,dtstrain):',&
                    ishift_deltat(isim),ishift_seisdt(isim),ishift_straindt(isim)
     write(6,*) ''  

  enddo

  ! input parameter consistency
  if (minval(dt)/=maxval(dt) .or. minval(nt)/=maxval(nt) .or. minval(period)/=maxval(period) .or. &
       minval(nrec)/=maxval(nrec)  .or. minval(nt_seis)/=maxval(nt_seis)  .or. &
       minval(dt_seis)/=maxval(dt_seis) .or. minval(nt_strain)/=maxval(nt_strain) .or. &
       minval(dt_strain)/=maxval(dt_strain) .or. minval(nt_snap)/=maxval(nt_snap) .or. &
       minval(dt_snap)/=maxval(dt_snap) ) then
     write(6,*)'PROBLEM with simulation.info parameters in the respective directories:'
     write(6,*)' one or more of the supposedly equal parameters differ!'
     stop
  endif


end subroutine read_input
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine compute_radiation_prefactor(mij_prefact,npts,nsim,longit)

  use data_all, only : src_type,magnitude,mij,sum_seis_true,correct_azi,rot_mat
  use data_all, only : trans_rot_mat,src_file_type,srccolat,srclon,any_sum_seis_true,outdir
  use global_par
  
  implicit none
  
  real, dimension(1:npts,1:nsim,1:3), intent(out)   :: mij_prefact
  real, dimension(1:npts), intent(in)               :: longit
  integer, intent(in)   :: npts,nsim
  real, dimension(6)    :: Mij_scale
  character(len=100)    :: junk
  integer               :: isim, i, j
  logical               :: file_exist
  real                  :: transrottmp(1:3,1:3), Mij_matr(3,3)
  
  ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007
  ! to rotate xyz coordinates 

  if (.not. allocated(trans_rot_mat)) allocate(trans_rot_mat(3,3,nsim))

  do isim=1,nsim
     rot_mat(1,1)=cos(srccolat(isim))*cos(srclon(isim))
     rot_mat(2,2)=cos(srclon(isim))
     rot_mat(3,3)=cos(srccolat(isim))
     rot_mat(2,1)=cos(srccolat(isim))*sin(srclon(isim))
     rot_mat(3,1)=-sin(srccolat(isim))
     rot_mat(3,2)=0.d0
     rot_mat(1,2)=-sin(srclon(isim))
     rot_mat(1,3)=sin(srccolat(isim))*cos(srclon(isim))
     rot_mat(2,3)=sin(srccolat(isim))*sin(srclon(isim))
     do i=1,3
        do j=1,3
           if (abs(rot_mat(i,j))<epsi_real) rot_mat(i,j)=0.0
        enddo
     enddo        
     trans_rot_mat(:,:,isim)=transpose(rot_mat)
  enddo

  if (src_file_type=='cmtsolut') then
     write(6,*)'  reading CMTSOLUTION file....'
     open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
     read(20000,*)junk
     read(20000,*)junk
     read(20000,*)junk
     read(20000,*)junk
     read(20000,*)junk
     read(20000,*)junk  
     read(20000,*)junk
     read(20000,*)junk,Mij(1) !Mrr
     read(20000,*)junk,Mij(2) !Mtt
     read(20000,*)junk,Mij(3) !Mpp
     read(20000,*)junk,Mij(4) !Mrt
     read(20000,*)junk,Mij(5) !Mrp
     read(20000,*)junk,Mij(6) !Mtp
     close(20000)

     Mij=Mij/1.E7 ! CMTSOLUTION given in dyn-cm

  elseif (src_file_type=='separate' .or. src_file_type=='sourceparams') then 
     open(unit=20000,file='sourceparams.dat',POSITION='REWIND',status='old')
     read(20000,*)(Mij(i),i=1,6)
     close(20000)        
  else
     write(6,*)'unknown source file type!',src_file_type
  endif

  write(6,*)'Original moment tensor: (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)'
  write(6,*)(Mij(i),i=1,6)
  write(6,*)'magnitudes of each run:'
  write(6,*)(magnitude(isim),isim=1,nsim)

  inquire(file=trim(outdir(1))//"/param_mij",exist=file_exist)
  if (file_exist) then
     write(6,*) 'OVERWRITING moment tensor...'
     open(unit=99,file=trim(outdir(1))//'/param_mij')
     read(99,*) (Mij(i),i=1,6)
     write(6,*) 'overwritten moment tensor:'
     write(6,*) (Mij(i),i=1,6)
     close(99)
  endif

  if (any_sum_seis_true) then
     do isim=1,nsim

        Mij_scale=Mij/magnitude(isim)

        write(6,*)'Mij scaled:',Mij_scale

        if ( (src_file_type=='separate'  .or. src_file_type=='sourceparams') &
                .and. .not. file_exist ) then 
           write(6,*)isim, 'rotating moment tensor from sourceparams..'
           transrottmp(1:3,1:3) = trans_rot_mat(:,:,isim)
           Mij_matr(1,1) = Mij_scale(1)
           Mij_matr(2,2) = Mij_scale(2)
           Mij_matr(3,3) = Mij_scale(3)
           Mij_matr(1,2) = Mij_scale(4)
           Mij_matr(1,3) = Mij_scale(5)
           Mij_matr(2,3) = Mij_scale(6)
           Mij_matr(2,1) = Mij_matr(1,2)
           Mij_matr(3,1) = Mij_matr(1,3)
           Mij_matr(3,2) = Mij_matr(2,3) 

           ! rotate Mij to source at NP system
           Mij_matr = matmul(transrottmp,Mij_matr)
           Mij_matr = matmul(Mij_matr,transpose(transrottmp))

           Mij_scale(1) = Mij_matr(1,1)
           Mij_scale(2) = Mij_matr(2,2)
           Mij_scale(3) = Mij_matr(3,3) 
           Mij_scale(4) = Mij_matr(1,2) 
           Mij_scale(5) = Mij_matr(1,3)
           Mij_scale(6) = Mij_matr(2,3)
            
       endif

        if (.not. correct_azi(isim)) then  ! just double checking the above for each simulation....
           do i=1,npts

              if (src_type(isim,2)=='mzz') then
                 mij_prefact(i,isim,:) = Mij_scale(1)
                 mij_prefact(i,isim,2) = 0.
                 if (i==1) write(6,*) isim, 'Simulation is mzz, prefact:', &
                                      mij_prefact(i,isim,1), mij_prefact(i,isim,2), &
                                      mij_prefact(i,isim,3)

              elseif (src_type(isim,2)=='mxx_p_myy') then
                 mij_prefact(i,isim,:) = Mij_scale(2)+Mij_scale(3)
                 mij_prefact(i,isim,2) = 0.
                 if (i==1) write(6,*) isim, 'Simulation is mxx, prefact:', &
                                      mij_prefact(i,isim,1), mij_prefact(i,isim,2), &
                                      mij_prefact(i,isim,3)

              elseif (src_type(isim,2)=='mxz' .or. src_type(isim,2)=='myz') then
                 mij_prefact(i,isim,:) = Mij_scale(4)*cos(longit(i))+Mij_scale(5)*sin(longit(i))
                 mij_prefact(i,isim,2) = -Mij_scale(4)*sin(longit(i))+Mij_scale(5)*cos(longit(i))

                 if (i==1) write(6,*) isim, 'Simulation is mxz, prefact:', &
                                      mij_prefact(i,isim,1), mij_prefact(i,isim,2), &
                                      mij_prefact(i,isim,3)

              elseif (src_type(isim,2)=='mxy' .or. src_type(isim,2)=='mxx_m_myy') then
                 mij_prefact(i,isim,:) = (Mij_scale(2)-Mij_scale(3))*cos(2.*longit(i))  &
                                                       +2.*Mij_scale(6)*sin(2.*longit(i)) 
                 mij_prefact(i,isim,2) = (Mij_scale(3)-Mij_scale(2))*sin(2.*longit(i)) &
                                                        +2.*Mij_scale(6)*cos(2.*longit(i))
                 if (i==1) write(6,*) isim, 'Simulation is mxy, prefact:', &
                                      mij_prefact(i,isim,1), mij_prefact(i,isim,2), &
                                      mij_prefact(i,isim,3)

              elseif(src_type(isim,2)=='explosion') then 
                 mij_prefact(i,isim,:) = (Mij_scale(1)+Mij_scale(2)+Mij_scale(3)) / 3.
                 if (i==1) write(6,*) isim, 'Simulation is explosion, prefact:', &
                                      mij_prefact(i,isim,1), mij_prefact(i,isim,2), &
                                      mij_prefact(i,isim,3)
              endif
              
           enddo

           write(6,*) 'Mij phi prefactor:', maxval(mij_prefact(:,isim,1)), &
                      maxval(mij_prefact(:,isim,2)), maxval(mij_prefact(:,isim,3))
        else ! correct_azi
           write(6,'(a,a,/,a,/,a,a)') &
                ' .... ISSUE: correct azimuth embedded in simulation, ', &
                'thus various components may be unretrievable!', &
                ' ..... hence we are NOT summing seismograms!', &
                ' .... to be sure to obtain correct results, rerun the ', &
                'solver with correc_azi set to .false. '

           sum_seis_true = .false.
        endif
     enddo ! isim
  endif ! any_sum_seis_true
end subroutine compute_radiation_prefactor
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine sum_individual_wavefields(field_sum, field_in, n, mij_prefact)

  implicit none
  
  integer, intent(in) :: n
  real, dimension(n,3), intent(in) :: field_in
  real, dimension(3), intent(in) :: mij_prefact
  real, dimension(n,3), intent(inout) :: field_sum

  field_sum(:,1) = field_sum(:,1) + mij_prefact(1)*field_in(:,1)
  field_sum(:,2) = field_sum(:,2) + mij_prefact(2)*field_in(:,2)
  field_sum(:,3) = field_sum(:,3) + mij_prefact(3)*field_in(:,3)

end subroutine sum_individual_wavefields
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine rotate_receiver_comp(isim, rec_comp_sys, srccolat, srclon, th_rot, ph_rot, &
                                th_orig, ph_orig, nt, seis)

  use data_all,     only : nsim, trans_rot_mat
  use global_par
  implicit none
  include 'mesh_params.h'
  
  character(len=3)      :: rec_comp_sys
  real, intent(in)      :: th_rot, ph_rot ! coordinates in the rotated (src at pole) system
  real, intent(in)      :: th_orig, ph_orig ! coordinates in the unrotated (actual src) system
  real, intent(in)      :: srccolat, srclon ! orginal source coordinates
  integer, intent(in)   :: nt,isim
  real, intent(inout)   :: seis(nt,3)
  real                  :: seis_tmp(nt,3), seisvec(3), rot(3,3)
  integer               :: i
  
  write(6,*) 'ROTATIONS'
  write(6,*) th_orig*180./pi, ph_orig*180./pi
  write(6,*) th_rot*180./pi, ph_rot*180./pi
  
  ! Need to consider *local* spherical geometry in the first place,
  ! THEN rotate the source-receiver frame to the north pole in the solver.
  ! E.g., consider the difference between source-projected and spherical coordinates for a source 
  ! away from the north pole: they are not the same, but in the framework below would 
  ! be identified as the same.
  
  ! Source projected frame: transform to spherical without any further rotations
  if (rec_comp_sys=='src') then  
     seis_tmp(:,1) = cos(th_rot) * seis(:,1) - sin(th_rot) * seis(:,3)
     seis_tmp(:,2) = seis(:,2)
     seis_tmp(:,3) = sin(th_rot) * seis(:,1) + cos(th_rot) * seis(:,3)
  
  ! Rotate from rotated u_sphiz to rotated u_xyz (both in reference, source-at-pole system) 
  else 
     seis_tmp(:,1) = cos(ph_rot) * seis(:,1) - sin(ph_rot) * seis(:,2) 
     seis_tmp(:,2) = sin(ph_rot) * seis(:,1) + cos(ph_rot) * seis(:,2)
     seis_tmp(:,3) = seis(:,3)
  
     ! Rotate to the original (i.e. real src-rec coordinate-based) u_xyz
     if (srccolat>epsi_real .or. srclon>epsi_real) then 
        rot=transpose(trans_rot_mat(:,:,isim))
        do i=1,nt
           seisvec = seis_tmp(i,:)
           seis_tmp(i,:) = matmul(rot,seisvec)
        enddo
     endif
  
  endif 
  
  ! Rotate to coordinate system of choice
  if (rec_comp_sys=='enz') then
     seis(:,1) = - cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & - cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + sin(th_orig) * seis_tmp(:,3)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) &
               & + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) & 
               & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + cos(th_orig) * seis_tmp(:,3)
     
  
  elseif (rec_comp_sys=='sph') then 
     seis(:,1) =   cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & + cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & - sin(th_orig) * seis_tmp(:,3)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) & 
               & + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) & 
               & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + cos(th_orig) * seis_tmp(:,3)
  
  elseif (rec_comp_sys=='cyl') then 
     seis(:,1) =   cos(ph_orig) * seis_tmp(:,1) + sin(ph_orig) * seis_tmp(:,2)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   seis_tmp(:,3)
  
  elseif (rec_comp_sys=='xyz') then
     seis = seis_tmp
  
  elseif (rec_comp_sys=='src') then
     seis = seis_tmp ! taken from above
  
  else
     write(6,*)'unknown component system',rec_comp_sys
     stop
  endif

end subroutine rotate_receiver_comp
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!! convolve seismograms computed for dirac delta with a Gaussian
subroutine convolve_with_stf(t_0,dt,nt,src_type,stf,outdir,seis,seis_fil)          
  
  use data_all, only : stf_type,shift_fact
  use global_par, only: pi,decay,shift_fact1
  implicit none
  
  integer, intent(in)            :: nt
  real, intent(in)               :: t_0,dt
  character(len=100), intent(in) :: outdir
  real                           :: time(nt)
  real                           :: tau_j,source,sqrt_pi_inv
  integer                        :: i,j,N_j,irec,lffile,ishift
  real, intent(in)               :: seis(nt,3)
  real, intent(out)              :: seis_fil(nt,3)
  real                           :: src_array(nt),temp_expo,alpha
  character(len=7), intent(in)   :: src_type
  character(len=7), intent(in)   :: stf
  character(len=4)               :: appidur,appirec
  logical                        :: monopole

  write(6,*)
  write(6,*)'Convolving with period=',t_0
  write(6,*)'convolve:',stf,stf_type(1),maxval(seis)

  monopole = .false. 
  if (src_type == 'monopole') monopole=.true.
  N_j = int(2.*shift_fact1*t_0/dt)
  if (N_j>nt) N_j = nt

  call define_io_appendix(appidur,int(t_0))
  alpha = decay/t_0
  sqrt_pi_inv = 1./dsqrt(pi)
  do i=1,nt
    time(i) = dt * real(i)
    seis_fil(i,:) = 0.
    do j=1, N_j
       tau_j=dble(j)*dt
       ! convolve with a Gaussian
       if (stf == 'gauss_0') then 
          temp_expo = alpha*(tau_j-shift_fact1*t_0)
          if (temp_expo < 50.) then
             source = alpha*exp(-temp_expo**2 )*sqrt_pi_inv / pi
          else
             source = 0.
          endif
       elseif (stf == 'quheavi') then 
          source = 0.5*(1.0+erf((tau_j-shift_fact1*t_0)/t_0))
       elseif (stf == 'gauss_1') then 
          source = -2.*(decay/t_0)**2*(tau_j-shift_fact1*t_0) * &
                           exp(-( (decay/t_0*(tau_j-shift_fact1*t_0))**2) )
          source=source/( decay/t_0*sqrt(2.)*exp(-2.) )
       else
          write(6,*)' other source time function not implemented yet!',stf
          stop
       endif
       if (i > j .and. i-j <= nt) seis_fil(i,:) = seis_fil(i,:)+seis(i-j,:)*source*dt
       if (i==1 ) src_array(j) = source
    enddo
  enddo

  seis_fil=seis_fil*pi
  write(6,*) 'convolve:', stf, stf_type(1), maxval(seis_fil)

  ! Output source time function as used here
  open(unit=55,file=trim(outdir)//'/stf_'//trim(stf)//'_'//appidur//'sec.dat')
  do i=1,N_j
    write(55,*) time(i), src_array(i)
  enddo
  close(55)

end subroutine convolve_with_stf
!=============================================================================

!-----------------------------------------------------------------------------
subroutine save_google_earth_kml(srccolat1, srclon1, srcdepth, Mij, per, rcvcolat, &
                                 rcvlon, reccomp, src_type, sum_seis_true, nsim, &
                                 num_rec_glob, outdir, receiver_name)

  use global_par, only : pi
  implicit none
  
  integer, intent(in)   :: num_rec_glob,nsim
  real, intent(in)      :: srccolat1, srclon1, srcdepth, rcvcolat(1:num_rec_glob), &
                           rcvlon(1:num_rec_glob)
  real, intent(in)      :: Mij(6),per
  logical, intent(in)   :: sum_seis_true

  character(len=100), intent(in)            :: receiver_name(1:num_rec_glob)
  character(len=100), intent(in)            :: outdir
  character(len=10), intent(in)             :: src_type
  character(len=1),dimension(3), intent(in) :: reccomp

  real                  :: slon, slat, rlon(1:num_rec_glob), rlat(1:num_rec_glob)
  integer               :: i
  character(len=4)      :: app
  character(len=2)      :: comp(3)
  character(len=100)    :: fname2
  
  write(6,*) 'writing google earth kml file for plotting earthquake and receiver locations/seismograms...'
  write(6,*) 'Check it out: '//trim(outdir)//'/googleearth_src_rec_seis.kml'
  
  slat=90.-srccolat1*180./pi
  slon=srclon1*180./pi
  if (slon>180.) slon=slon-360.
  
  rlat=90.-rcvcolat*180./pi
  rlon=rcvlon*180./pi
  do i=1,num_rec_glob
     if (rlon(i)>180.) rlon(i)=rlon(i)-360.
  enddo
  open(unit=88,file=trim(outdir)//'/googleearth_src_rec_seis.kml')
  
  write(88,14)'<?xml version="1.0" encoding="UTF-8"?> '
  write(88,15)'<kml xmlns="http://earth.google.com/kml/2.0"> '
  write(88,16)'<Document> '
  
  write(88,*)
  write(88,*)'  <name> earthquake-receiver configuration</name>'
  write(88,*)'    <LookAt>'
  write(88,12)'     <longitude>',slon,'</longitude><latitude>',slat,'</latitude>'
  write(88,*)'     <range>7000000</range><tilt>0</tilt><heading>0</heading>'
  write(88,*)'    </LookAt>'
  write(88,*)
  write(88,*)'......'
  write(88,*)'  <Placemark>'
  write(88,*)'     <Style id="earthquake">'
  write(88,*)'       <IconStyle>'
  write(88,*)'       <scale>5</scale>'
  write(88,*)'         <Icon>'
  write(88,*)' <href>http://maps.google.com/mapfiles/kml/shapes/earthquake.png</href>'
  write(88,*)'             </Icon>'
  write(88,*)'           </IconStyle>'
  write(88,*)'                  <LabelStyle>'
  write(88,*)'                      <scale>5</scale>'
  write(88,*)'                 </LabelStyle>'
  write(88,*)'        </Style>'
  write(88,*)'    <name> earthquake</name>'
  write(88,*) ' <description> Event details:'
  write(88,20) ' colat,lon [deg]:',srccolat1*180./pi,srclon1*180./pi
  write(88,21)' source depth [km]',srcdepth
  write(88,23)'Mrr=',Mij(1)
  write(88,23)'Mtt=',Mij(2)
  write(88,23)'Mpp=',Mij(3)
  write(88,23)'Mtr=',Mij(4)
  write(88,23)'Mpr=',Mij(5)
  write(88,23)'Mtp=',Mij(6)
  write(88,21)'source period [s]:',per
  write(88,*)'</description>'
  write(88,13)'   <Point><coordinates>',slon,',',slat,'</coordinates></Point>'
  write(88,*)'   </Placemark>'
  
  do i=1,num_rec_glob
     write(88,*)
     write(88,*) ' <Placemark>'
     write(88,*) '     <Style id="cam">'
     write(88,*) '       <IconStyle>'
   write(88,*)'       <scale>2</scale>'
     write(88,*) '         <Icon>'
        write(88,*)' <href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>'
     write(88,*) '         </Icon>'
     write(88,*) '       </IconStyle>'
  write(88,*)'                  <LabelStyle>'
  write(88,*)'                      <scale>2</scale>'
  write(88,*)'                 </LabelStyle>'
     write(88,*) '     </Style>'
     write(88,17) ' <name> ',trim(receiver_name(i)),'  # ',i,'</name>'
     call define_io_appendix(app,i)
     write(88,119) ' <description> station ',trim(receiver_name(i))
     write(88,20) ' colat,lon [deg]:',rcvcolat(i)*180./pi,rcvlon(i)*180./pi
  
     fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(1)//'.gif'
     write(88,*) ' <img src="',trim(fname2),'"></img>'
  
     fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(3)//'.gif'
     write(88,*) ' <img src="',trim(fname2),'"></img>'
  
     if (sum_seis_true .or. nsim>1 .or. src_type/='monopole') then 
        fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(2)//'.gif'
        write(88,*) ' <img src="',trim(fname2),'"></img>'
     endif
     write(88,*) '  </description>'
     write(88,13) '   <Point><coordinates>',rlon(i),',',rlat(i),'</coordinates></Point>'
     write(88,*) ' </Placemark>'
  enddo
  
  write(88,*)'......'
  write(88,*)
  write(88,*)'</Document>'
  write(88,*)'</kml>'
  
  close(88)

12 format(a16,f14.2,a23,f14.2,a12)
13 format(a23,f14.2,a1,f14.2,a23)
14 format(a39)
15 format(a46)
16 format(a11)
17 format(a7,a9,a10,i4,a7)
18 format(a36,a4,a14)
19 format(a24,a15)
119 format(a24,a30)
20 format(A18,f14.2,f14.2)
21 format(A18,f14.2)
23 format(a5,1pe14.2)

end subroutine save_google_earth_kml
!=============================================================================

!-----------------------------------------------------------------------------
subroutine define_io_appendix(app,iproc)
! Defines the 4 digit character string appended to any 
! data or io file related to process myid. 

  integer, intent(in)           :: iproc
  character(len=4), intent(out) :: app
  
  write(app,"(I4.4)") iproc

end subroutine define_io_appendix
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_3d_wavefields

  use data_all
  use global_par
  implicit none

  include 'mesh_params.h'

  integer                               :: iproc, npts, npts_top, npts_bot, npts_meri, &
                                           nphi, snapskip, snap1, snap2
  real, dimension(:,:), allocatable     :: coord, disp, disptot, disptot_sum, azi_xsect2, &
                                           azi_xsect1
  real, dimension(:,:,:), allocatable   :: azi, fluid_prefact, azi_meri, mij_snap
  real, dimension(:), allocatable       :: vp, x, y, z, vptop, vpbot, vpmeri, xtot, ytot, &
                                           ztot, azi_phi_meri
  double precision, dimension(:), allocatable :: xtop, ytop, ztop, xbot, ybot, zbot, &
                                                 azi_phi_top, azi_phi_bot
  real, dimension(:), allocatable       :: xmeri, ymeri, zmeri, longit_snap
  real                                  :: dphi, phi0, prem, r, theta_meri, smallval_meri, &
                                           dr, theta0, theta_bot, theta_top, phi
  integer, dimension(:), allocatable    :: ind_proc_top_tmp, ind_pts_top_tmp, &
                                           ind_proc_bot_tmp, ind_pts_bot_tmp
  integer, dimension(:), allocatable    :: ind_proc_top, ind_pts_top, ind_proc_bot
  integer, dimension(:), allocatable    :: ind_pts_bot, azi_ind_top, azi_ind_bot, &
                                           azi_ind_meri
  integer, dimension(:), allocatable    :: ind_proc_meri_tmp, ind_pts_meri_tmp, &
                                           ind_proc_meri, ind_pts_meri
  real                                  :: smallval_north_top,  smallval_south_top, &
                                           smallval_north_bot,  smallval_south_bot, r0, &
                                           coord9(9, 2), disp9(9, 3)
  character(len=200)                    :: filename1
  logical                               :: use_meri, use_top, use_bot

  character(len=4)                      :: appmynum2
  integer                               :: nptstot, k1, k2, k3, npts_fluid, j, &
                                           naang, npts_read, npts_fluid_read
  real                                  :: rbot, rtop, phi_incr
  double precision                      :: r8, theta8, phi8

  ! read snap plot parameters

  open(unit=99,file=trim(simdir(1))//'/param_snaps')
  read(99,*) phi0
  write(6,*) 'starting azimuth/phi for cross section on the right [deg]:',phi0
  read(99,*) dphi
  write(6,*) 'ending azimuth/phi for cross section on the left [deg]:',dphi
  read(99,*) rtop
  write(6,*) 'top surface [km]:',rtop
  read(99,*) rbot
  write(6,*) 'bottom surface [km]:',rbot
  read(99,*) theta_meri
  write(6,*) 'colatitude of meridional cross section:',theta_meri
  read(99,*) snap1,snap2,snapskip
  write(6,*) '1st,last snap, skipfactor:',snap1,snap2,snapskip
  read(99,*) use_meri
  write(6,*) 'consider meridional cross section?',use_meri
  read(99,*) use_top
  write(6,*) 'consider top surface?',use_top
  read(99,*) use_bot
  write(6,*) 'consider bottom surface?',use_bot
  close(99)

  phi0 = phi0 / 180. * pi
  dphi = dphi / 180. * pi
  rtop = rtop * 1000.
  rbot = rbot * 1000.
  theta_meri = theta_meri * pi / 180.

  ! if all the same or not!!!!!
  if (minval(dt)/=maxval(dt) .or. minval(nt)/=maxval(nt) .or. minval(period)/=maxval(period) .or. &
      minval(nrec)/=maxval(nrec)  .or. minval(nt_seis)/=maxval(nt_seis)  .or. &
      minval(dt_seis)/=maxval(dt_seis) .or. minval(nt_strain)/=maxval(nt_strain) .or. &
      minval(dt_strain)/=maxval(dt_strain) .or. minval(nt_snap)/=maxval(nt_snap) .or. &
      minval(dt_snap)/=maxval(dt_snap) ) then
     write(6,*) 'PROBLEM with simulation.info parameters in the respective directories:'
     write(6,*) ' one or more of the supposedly equal parameters differ!'
     call flush(6)
     stop
  endif

  npts = nelem * 16
  npts_read = nelem * 9

  nptstot = npts * nproc_mesh
  write(6,*)'number of points per proc, total points:',npts,nptstot

  ! load and construct global mesh (one semi-disk)
  write(6,*)'reading partitioned mesh...'
  allocate(coord(nptstot,2))
  smallval_north_top = rtop
  smallval_south_top = rtop
  smallval_north_bot = rbot
  smallval_south_bot = rbot

  write(6,*)'Smallest distance to rtop (North,South) BEFORE [km]:', &
       real(smallval_north_top/1000.),real(smallval_south_top/1000.)
  write(6,*)'Smallest distance to rbot (North,South) BEFORE [km]:', &
       real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)

  do iproc=0,nproc_mesh-1
     call define_io_appendix(appmynum,iproc)
     open(unit=99,file=trim(simdir(1))//'/Data/glob_grid_'//appmynum//'.dat')
     i=0
     do ii=1,npts_read/9
        do iii=1,9
            read(99,*) coord9(iii,1), coord9(iii,2)
         enddo
         coord(iproc*npts+i+1:iproc*npts+i+2,:) = coord9(1:2,:) 
                                                  ! (0,0), (0,npol/2)
         coord(iproc*npts+i+3,:) = coord9(5,:)    ! (npol/2,npol/2)
         coord(iproc*npts+i+4,:) = coord9(4,:)    ! (npol/2,0)
         coord(iproc*npts+i+5:iproc*npts+i+6,:) = coord9(2:3,:)    
                                                  ! (0,npol/2),(0,npol)
         coord(iproc*npts+i+7,:) = coord9(6,:)    ! (npol/2,npol)
         coord(iproc*npts+i+8,:) = coord9(5,:)    ! (npol/2,npol/2)
         coord(iproc*npts+i+9:iproc*npts+i+10,:) = coord9(4:5,:)    
                                                  ! (npol/2,0),(npol/2,npol/2)
         coord(iproc*npts+i+11,:) = coord9(8,:)   ! (npol,npol/2)
         coord(iproc*npts+i+12,:) = coord9(7,:)   ! (npol,0)            
         coord(iproc*npts+i+13:iproc*npts+i+14,:) = coord9(5:6,:)
                                                  ! (npol/2,npol/2),(npol/2,npol)
         coord(iproc*npts+i+15,:) = coord9(9,:)   ! (npol,npol)
         coord(iproc*npts+i+16,:) = coord9(8,:)   ! (npol,npol/2)            

         ! determine minimal distance from rtop and rbot
         do iii=1,16
            r0 = sqrt(coord(iproc*npts+i+iii,1)**2 + coord(iproc*npts+i+iii,2)**2) 
            if (coord(iproc*npts+i+iii,2) > 0.0 &
                    .and. abs(r0-rtop) < smallval_north_top) then ! north
               smallval_north_top = abs(r0-rtop)
            elseif (coord(iproc*npts+i+iii,2) < 0.0 &
                    .and. abs(r0-rtop) < smallval_south_top) then ! south
               smallval_south_top = abs(r0-rtop)
            endif
            
            if (coord(iproc*npts+i+iii,2) > 0.0 &
                    .and. abs(r0-rbot)< smallval_north_bot) then ! north            
               smallval_north_bot = abs(r0-rbot)
            elseif (coord(iproc*npts+i+iii,2) < 0.0 &
                    .and. abs(r0 -rbot)< smallval_south_bot) then ! south
               smallval_south_bot = abs(r0-rbot)
            endif
        enddo
        i=i+16
     enddo
     close(99)
  enddo

  smallval_north_top = smallval_north_top + epsi_real 
  smallval_south_top = smallval_south_top + epsi_real 
  smallval_north_bot = smallval_north_bot + epsi_real 
  smallval_south_bot = smallval_south_bot + epsi_real 

  write(6,*) 'Smallest distance to rtop (North,South) AFTER [km]:', &
       real(smallval_north_top/1000.),real(smallval_south_top/1000.)
  write(6,*) 'Smallest distance to rbot (North,South) AFTER [km]:', &
       real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)

  if (use_top .or. use_meri) then
     allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_top_tmp(floor(real(nptstot)/10.)))
  endif

  if (use_bot .or. use_meri) then
     allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_bot_tmp(floor(real(nptstot)/10.)))
  endif

  k1 = 0
  k2 = 0 

  do iproc=0, nproc_mesh-1
     do i=1, npts
        ! check for top and bottom radii
        r0 = sqrt(coord(iproc*npts+i,1)**2 + coord(iproc*npts+i,2)**2) 

        if (use_top .or. use_meri) then
           if ( (coord(iproc*npts+i,2) >= 0.0 .and. abs(r0-rtop) <= smallval_north_top) .or. &
                (coord(iproc*npts+i,2) < 0.0 .and.  abs(r0-rtop) <= smallval_south_top)) then 
              k1 = k1 + 1         
              ind_proc_top_tmp(k1) = iproc
              ind_pts_top_tmp(k1) = i
           endif
        endif

        if (use_bot .or. use_meri) then 
           if ( (coord(iproc*npts+i,2) >= 0.0 .and. abs(r0-rbot) <= smallval_north_bot) .or. &
                (coord(iproc*npts+i,2) < 0.0 .and. abs(r0-rbot) <= smallval_south_bot)) then 
              k2 = k2 + 1
              ind_proc_bot_tmp(k2) = iproc
              ind_pts_bot_tmp(k2) = i
           endif
        endif
     enddo
  enddo
  
  npts_top = k1
  npts_bot = k2

  write(6,*) '# points on top,bottom:', npts_top, npts_bot

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

  !----------------------------------------------------------------------
  ! MERIDIONAL based on rtop and rbottom
  if (use_meri) then

     write(6,*)'computing meridional preparameters....'

     ! find closest theta at rbot
     smallval_meri = 2.d0*pi
     do i=1,npts_bot
        r0 = sqrt(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)**2 &
                    + coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)**2)
        theta0 = atan(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1) &
                        /coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)+epsi)
        if ( theta0 <0. ) theta0 = pi + theta0
        if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
        if (abs(theta_meri-theta0) < smallval_meri) then 
           smallval_meri = abs(theta_meri-theta0) 
           theta_bot = theta0 
        endif
     enddo
     write(6,*)'theta meri,theta closest at rbot:',theta_meri*180./pi,theta_bot*180./pi
     theta_meri = theta_bot

     ! find theta at rtop closest to theta from rbot
     smallval_meri = 2.d0*pi
     do i=1,npts_top
        r0 = sqrt(coord(ind_proc_top(i)*npts+ind_pts_top(i),1)**2 &
                    + coord(ind_proc_top(i)*npts+ind_pts_top(i),2)**2)
        theta0 = atan(coord(ind_proc_top(i)*npts+ind_pts_top(i),1) &
                        / coord(ind_proc_top(i)*npts+ind_pts_top(i),2)+epsi)
        if ( theta0 <0. ) theta0 = pi + theta0
        if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
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
     allocate(ind_proc_meri_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_meri_tmp(floor(real(nptstot)/10.)))

     do iproc=0,nproc_mesh-1
        do i=1,npts
           r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 
           theta0=atan(coord(iproc*npts+i,1)/(coord(iproc*npts+i,2)+epsi))
           if ( theta0 <0. ) theta0 = pi + theta0
           if (theta0==0.0 .and. coord(iproc*npts+i,2)<0.0) theta0=pi
           if ( r0>=rbot .and. abs(theta0-theta_meri) <= smallval_meri ) then
              k3=k3+1
              ind_proc_meri_tmp(k3) = iproc
              ind_pts_meri_tmp(k3) = i
           endif
        enddo
     enddo
     npts_meri=k3

     write(6,*)'# points on meridional:',npts_meri
     allocate(ind_proc_meri(npts_meri),ind_pts_meri(npts_meri))
     ind_proc_meri = ind_proc_meri_tmp(1:npts_meri)
     ind_pts_meri = ind_pts_meri_tmp(1:npts_meri)
     deallocate(ind_proc_meri_tmp,ind_pts_meri_tmp)
  endif ! use_meri

  ! xyz coordinates-----------------------------------------------------------------------
  write(6,*)'defining xyz...'
  allocate(x(1:2*nptstot),y(1:2*nptstot),z(1:2*nptstot));x=0.;y=0.;z=0.

  ! left cross section--------------------------------------------------------------------
  call sphi2xy(x(1:nptstot),y(1:nptstot),coord(:,1),phi0,nptstot)
  z(1:nptstot)=coord(1:nptstot,2)
  write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

  ! right cross section-------------------------------------------------------------------
  call sphi2xy(x(nptstot+1:2*nptstot),y(nptstot+1:2*nptstot),coord(:,1),phi0+dphi,nptstot)
  z(nptstot+1:2*nptstot)=coord(1:nptstot,2)

  write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

  ! save xyz
  allocate(vp(2*nptstot))
  do i=1,2*nptstot  
     r= sqrt(x(i)**2+y(i)**2+z(i)**2)
     vp(i)=prem(r,'v_p')
  enddo

  ! rescale vp
  vp=vp/maxval(vp)*0.002-0.001

  filename1=trim(outdir(1))//'/SNAPS/mesh_xsect'
  call write_VTK_bin_scal(x,y,z,vp,2*nptstot,0,filename1)
  filename1=trim(outdir(1))//'/SNAPS/mesh_xsect_cell'
  call write_VTK_bin_scal_topology(x,y,z,vp,2*nptstot/4,filename1)


  ! top surface--------------------------------------------------------------------------
  if (use_top) then
     k1 = 0
     write(6,*) 'defining top surface...'

     naang = floor(real(npts_top)/real(2.))**2*6*4
     write(6,*)'points on top surface:',naang
     allocate(xtop(1:naang),ytop(1:naang),ztop(1:naang))
     allocate(azi_phi_top(1:naang),azi_ind_top(1:naang))
     call construct_surface_cubed_sphere(npts_top,npts,dble(rtop),ind_proc_top, &
                                         ind_pts_top,nptstot,dble(coord),k1, &
                                         dble(dphi),dble(phi0),'outside',naang,xtop,&
                                         ytop,ztop,azi_phi_top,azi_ind_top)

     ! extract vp
     write(6,*)'allocating vptop...',k1
     allocate(vptop(k1))
     vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))

     write(6,*)'save into cell vtk...',size(xtop),k1
     call flush(6)
     filename1=trim(outdir(1))//'/SNAPS/mesh_top_cell'
     call write_VTK_bin_scal_topology(real(xtop(1:k1)),real(ytop(1:k1)),real(ztop(1:k1)),&
                                      vptop(1:k1),k1/4,filename1)

     call flush(6)
  endif

  ! bottom surface -----------------------------------------------------------------------
  if (use_bot) then
     k2=0
     write(6,*)'defining bot surface...'

     naang = floor(real(npts_bot)/real(2.))**2*6*4
     write(6,*)'points on bottom surface:',naang
     allocate(xbot(1:naang),ybot(1:naang),zbot(1:naang))
     allocate(azi_phi_bot(1:naang),azi_ind_bot(1:naang))

     call construct_surface_cubed_sphere(npts_bot,npts,dble(rbot),ind_proc_bot, &
                                         ind_pts_bot,nptstot,dble(coord),k2, &
                                         dble(dphi),dble(phi0),'innside',naang,xbot,&
                                         ybot,zbot,azi_phi_bot,azi_ind_bot)

     ! extract vp
     allocate(vpbot(k2))
     vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))

     ! save into cell vtk
     filename1=trim(outdir(1))//'/SNAPS/mesh_bot_cell'
     call write_VTK_bin_scal_topology(real(xbot(1:k2)),real(ybot(1:k2)),real(zbot(1:k2)),&
                                      vpbot,k2/4,filename1)

     call flush(6)
  endif

  ! meridional cross section -------------------------------------------------------------
  k3=0
  if (use_meri) then 
     write(6,*)'defining meridional cross section...'
     dr=(6371000.-rbot)/npts_meri
     write(6,*)'# pts on rmeri and average spacing [km]:',npts_meri,dr/1000.
     allocate(xmeri(1:7*npts_meri**2))
     allocate(ymeri(1:7*npts_meri**2))
     allocate(zmeri(1:7*npts_meri**2))
     xmeri=0.
     ymeri=0.
     zmeri=0.
     allocate(azi_ind_meri(7*npts_meri**2),azi_phi_meri(7*npts_meri**2))
     allocate(vpmeri(1:7*npts_meri**2))

     do i=1,npts_meri
        r0 = sqrt( coord(ind_proc_meri(i)*npts+ind_pts_meri(i),1)**2 &
                    + coord(ind_proc_meri(i)*npts+ind_pts_meri(i),2)**2)
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
     filename1=trim(outdir(1))//'/SNAPS/mesh_meri'
     call write_VTK_bin_scal(xmeri,ymeri,zmeri,vpmeri,k3,0,filename1)

  endif !use_meri

  deallocate(coord)
  deallocate(vp)

  ! assembling everything to one coordinate array-----------------------------------------

  allocate(xtot(2*nptstot+k1+k2+k3),ytot(2*nptstot+k1+k2+k3),ztot(2*nptstot+k1+k2+k3))

  xtot(1:2*nptstot)=x(1:2*nptstot)
  ytot(1:2*nptstot)=y(1:2*nptstot)
  ztot(1:2*nptstot)=z(1:2*nptstot)

  if (use_top) then 
     xtot(2*nptstot+1:2*nptstot+k1)=real(xtop(1:k1))
     ytot(2*nptstot+1:2*nptstot+k1)=real(ytop(1:k1))
     ztot(2*nptstot+1:2*nptstot+k1)=real(ztop(1:k1))
  endif
  if (use_bot) then 
     xtot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(xbot(1:k2))
     ytot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(ybot(1:k2))
     ztot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(zbot(1:k2))
  endif
  if (use_meri) then 
     xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=xmeri(1:k3)
     ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=ymeri(1:k3)
     ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=zmeri(1:k3)
  endif

  ! save mesh for kerner
  open(unit=99,file=trim(outdir(1))//'/SNAPS/mesh_tot.xyz')
  do i=1,2*nptstot+k1+k2
     write(99,*)xtot(i),ytot(i),ztot(i)
  enddo
  close(99)

  if (use_meri) then 
     open(unit=99,file=trim(outdir(1))//'/SNAPS/mesh_meri.xyz')
     do i=2*nptstot+k1+k2+1,2*nptstot+k1+k2+k3
        write(99,*)xtot(i),ytot(i),ztot(i)
     enddo
     close(99)
  endif

  ! load azimuthal prefactors-------------------------------------------------------------
  if (any_sum_seis_true) then 
     allocate(longit_snap(2*nptstot+k1+k2))
     longit_snap(1:nptstot) =phi0
     longit_snap(nptstot+1:2*nptstot) = phi0+dphi
     do i=1,k1
        longit_snap(2*nptstot+i) = azi_phi_top(i)
     enddo
     do i=1,k2
        longit_snap(2*nptstot+k1+i) = azi_phi_bot(i)
     enddo
     do i=1,k3
        longit_snap(2*nptstot+k1+k2+i) = azi_phi_meri(i)
     enddo

     allocate(mij_snap(2*nptstot+k1+k2+k3,nsim,3))
     mij_snap= -1.E30
     call compute_radiation_prefactor(mij_snap,2*nptstot+k1+k2+k3,nsim,longit_snap)
     allocate(disptot_sum(2*nptstot+k1+k2+k3,3))
  endif

  ! FLUID REGION
  npts_fluid=nel_fluid*16
  npts_fluid_read=nel_fluid*9
  write(6,*)'points in fluid:',npts_fluid
  allocate(fluid_prefact(npts_fluid*nproc_mesh,2,nsim))

  do isim = 1,nsim
     ! load fluid prefactors-------------------------------------------------------------
     write(6,*) 'loading fluid prefactors...'
     do iproc=0, nproc_mesh-1
        call define_io_appendix(appmynum,iproc)
        open(unit=190,file=trim(simdir(isim))//'/Data/inv_rho_s_fluid_globsnaps_'//appmynum//'.dat')
        i=0
        do ii=1,npts_fluid_read/9
           do iii=1,9
              read(190,*)coord9(iii,1),coord9(iii,2)
           enddo
           fluid_prefact(iproc*npts_fluid+i+1:iproc*npts_fluid+i+2,:,isim) = coord9(1:2,:)
                                                                       ! (0,0), (0,npol/2)
           fluid_prefact(iproc*npts_fluid+i+3,:,isim) = coord9(5,:)    ! (npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+4,:,isim) = coord9(4,:)    ! (npol/2,0)
           fluid_prefact(iproc*npts_fluid+i+5:iproc*npts_fluid+i+6,:,isim) = coord9(2:3,:)    
                                                                       ! (0,npol/2),(0,npol)
           fluid_prefact(iproc*npts_fluid+i+7,:,isim) = coord9(6,:)    ! (npol/2,npol)
           fluid_prefact(iproc*npts_fluid+i+8,:,isim) = coord9(5,:)    ! (npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+9:iproc*npts_fluid+i+10,:,isim) = coord9(4:5,:)    
                                                                  ! (npol/2,0),(npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+11,:,isim) = coord9(8,:)   ! (npol,npol/2)
           fluid_prefact(iproc*npts_fluid+i+12,:,isim) = coord9(7,:)   ! (npol,0)            
           fluid_prefact(iproc*npts_fluid+i+13:iproc*npts_fluid+i+14,:,isim) = coord9(5:6,:)    
                                                                  ! (npol/2,npol/2),(npol/2,npol)
           fluid_prefact(iproc*npts_fluid+i+15,:,isim) = coord9(9,:)   ! (npol,npol)
           fluid_prefact(iproc*npts_fluid+i+16,:,isim) = coord9(8,:)   ! (npol,npol/2)            
           i=i+16
        enddo  ! npts_fluid_read
        close(190)
     enddo ! nproc_mesh
  enddo !isim

  ! compute longitude and radiation pattern for multiple wavefields

  ! load snaps ===============================================================
  write(6,*)'loading snaps...'
  allocate(disp(2*nptstot+k1+k2+k3,3))

  do j=snap1,snap2,snapskip

     disp=0.
     if (any_sum_seis_true) disptot_sum = 0.

     do isim=1,nsim
        do iproc=0,nproc_mesh-1
           call define_io_appendix(appmynum,iproc)
           call define_io_appendix(appmynum2,j)      
           open(unit=99,file=trim(simdir(isim))//'/Data/snap_'//appmynum//'_'//appmynum2//'.dat', &
                     FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
           i=0
           do ii=1,npts_read/9
              do iii=1,9
                 read(99)disp9(iii,1),disp9(iii,2),disp9(iii,3)
              enddo
              disp(iproc*npts+i+1:iproc*npts+i+2,:) = disp9(1:2,:) ! (0,0), (0,npol/2)
              disp(iproc*npts+i+3,:) = disp9(5,:)    ! (npol/2,npol/2)
              disp(iproc*npts+i+4,:) = disp9(4,:)    ! (npol/2,0)
              disp(iproc*npts+i+5:iproc*npts+i+6,:) = disp9(2:3,:)    ! (0,npol/2),(0,npol)
              disp(iproc*npts+i+7,:) = disp9(6,:)    ! (npol/2,npol)
              disp(iproc*npts+i+8,:) = disp9(5,:)    ! (npol/2,npol/2)
              disp(iproc*npts+i+9:iproc*npts+i+10,:) = disp9(4:5,:)    ! (npol/2,0),(npol/2,npol/2)
              disp(iproc*npts+i+11,:) = disp9(8,:)    ! (npol,npol/2)
              disp(iproc*npts+i+12,:) = disp9(7,:)    ! (npol,0)            
              disp(iproc*npts+i+13:iproc*npts+i+14,:) = disp9(5:6,:)    ! (npol/2,npol/2),(npol/2,npol)
              disp(iproc*npts+i+15,:) = disp9(9,:)    ! (npol,npol)
              disp(iproc*npts+i+16,:) = disp9(8,:)    ! (npol,npol/2)            
              i=i+16
           enddo
           close(99)

           disp(iproc*npts+1:iproc*npts+npts_fluid,1)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,1) * fluid_prefact(:,1,isim)
           disp(iproc*npts+1:iproc*npts+npts_fluid,2)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,2) * fluid_prefact(:,1,isim) &
                    * fluid_prefact(:,2,isim)
           disp(iproc*npts+1:iproc*npts+npts_fluid,3)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,3) * fluid_prefact(:,1,isim)
        enddo

        disp(nptstot+1:2*nptstot,1)=disp(1:nptstot,1)
        disp(nptstot+1:2*nptstot,2)=disp(1:nptstot,2)
        disp(nptstot+1:2*nptstot,3)=disp(1:nptstot,3)

        do i=1,k1
           disp(2*nptstot+i,1) = disp(azi_ind_top(i),1)
           disp(2*nptstot+i,2) = disp(azi_ind_top(i),2)
           disp(2*nptstot+i,3) = disp(azi_ind_top(i),3)
        enddo
        do i=1,k2
           disp(2*nptstot+k1+i,1) = disp(azi_ind_bot(i),1)
           disp(2*nptstot+k1+i,2) = disp(azi_ind_bot(i),2)
           disp(2*nptstot+k1+i,3) = disp(azi_ind_bot(i),3)
        enddo
        if (use_meri) then 
           do i=1,k3
              disp(2*nptstot+k1+k2+i,1) = disp(azi_ind_meri(i),1)   
              disp(2*nptstot+k1+k2+i,2) = disp(azi_ind_meri(i),2)
              disp(2*nptstot+k1+k2+i,3) = disp(azi_ind_meri(i),3)
           enddo
        endif

        disp(1:nptstot,1)=disp(1:nptstot,1)
        disp(1:nptstot,2)=disp(1:nptstot,2)
        disp(1:nptstot,3)=disp(1:nptstot,3)

        filename1 = trim(simdir(isim))//'/'//trim(outdir(1))//'/SNAPS/snap_cell_'&
                    //trim(src_type(isim,2))//'_'//appmynum2//'_z'
        write(6,*)'filename out vtk :',filename1
        call write_VTK_bin_scal_topology(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2),&
                                         ztot(1:2*nptstot+k1+k2), & 
                                         disp(1:2*nptstot+k1+k2,3),(2*nptstot+k1+k2)/4,&
                                         filename1)

        
        if (use_meri) then
           filename1 = trim(simdir(isim))//'/'//trim(outdir(1))//'/SNAPS/meri_snap_'&
                        //trim(src_type(isim,2))//'_'//appmynum2//'_z'
           call write_VTK_bin_scal(xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3),&
                ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                disp(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3,3),k3,0,filename1)
        endif

        if (sum_seis_true(isim)) then 
           disptot_sum(:,1) = disptot_sum(:,1) + mij_snap(:,isim,1)*disp(:,2)
           disptot_sum(:,2) = disptot_sum(:,2) + mij_snap(:,isim,2)*disp(:,2)
           disptot_sum(:,3) = disptot_sum(:,3) + mij_snap(:,isim,3)*disp(:,3)
        endif

     enddo
 
     if (any_sum_seis_true .and. nsim > 1) then 

        filename1=trim(outdir(1))//'/SNAPS/snap_mij_cell_'//appmynum2//'_z'
        write(6,*)'filename out vtk :',filename1
        call write_VTK_bin_scal_topology(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2),&
                                         ztot(1:2*nptstot+k1+k2), & 
                                         disptot_sum(1:2*nptstot+k1+k2,3),&
                                         (2*nptstot+k1+k2)/4,filename1)
        
        if (use_meri) then
           filename1=trim(outdir(1))//'/SNAPS/meri_snap_mij_'//appmynum2//'_z'
           call write_VTK_bin_scal(xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3),&
                ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                disptot_sum(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3,3),k3,0,filename1)
        endif
     endif
  enddo

end subroutine compute_3d_wavefields
!=============================================================================

!-----------------------------------------------------------------------------
subroutine construct_surface_cubed_sphere(npts_surf,npts,rsurf,ind_proc_surf,ind_pts_surf,&
                            nptstot,coord1,kpts,dphi,phi0,in_or_out,n,xsurf,ysurf,zsurf,azi_phi_surf,azi_ind_surf)

!!!!! BEG CUBED SPHERE
! af2tnm: along a great circle, there's the equivalent of two chunks of the cubed
! sphere. According to the convention we use in defining the mapping in the cubed sphere
! module, that amounts to 2*nang spectral elements.
! Assuming we will be using npol_cs=1 in the following, nang has to be even.
! We therefore want nang=.5*(npts_surf-1) if npts_surf is odd
! We therefore want nang=.5*(npts_surf) if npts_surf is even
!
!use data_all
    implicit none
    
    integer, intent(in) :: npts_surf,npts,nptstot,n
    double precision, intent(in) :: rsurf,dphi,phi0
    character(len=7), intent(in) :: in_or_out
    integer, dimension(1:npts_surf), intent(in) :: ind_proc_surf,ind_pts_surf
    double precision, dimension(nptstot,2), intent(in) :: coord1
    double precision, dimension(1:n), intent(out) :: xsurf,ysurf,zsurf,azi_phi_surf
    integer, dimension(1:n), intent(out) :: azi_ind_surf
    integer, intent(out) :: kpts
    
    double precision, allocatable,dimension(:) :: xi_el,eta_el,r_el
    integer, allocatable, dimension(:,:,:,:) :: number_el
    double precision, allocatable,dimension(:) :: xi_i,xi_j,xi_k
    double precision, allocatable,dimension(:,:,:,:) :: xcol,ycol,zcol
    double precision, allocatable,dimension(:,:,:,:) :: x_el,y_el,z_el
    double precision :: dang,C,D,re,ri,teta,tksi,tr,Xe,Ye,delta
    integer :: npol_cs,ii,iii,izone,nang,nelt,nr,nrpol,iel,nel_surf,jj,ipol,jpol,kpol,j,i
    double precision ::  dist,r_ref,th_ref,th,dr,r,xc,yc,zc,ph
    double precision, parameter :: pi = 3.1415926535898

    write(6,*)'computing cubed sphere for surface at r=',rsurf
    write(6,*)'pts surf,total:',npts_surf,npts

    nang=floor(sqrt(real(n))/6./2.)
    write(6,*)'naang,nang,npts_surf:',n,nang,npts_surf

     nr=1 ! one radial layer only
     write(6,*)'NANG:',nang
     allocate(xi_el(0:nang),eta_el(0:nang),r_el(0:nr))
     xi_el = 0.d0;  eta_el = 0.d0 ; r_el = 0.d0
     re=rsurf
     ri=rsurf-1.d0 ! (in km)
     !
     dang=pi/(2.d0*dble(nang))
     dr=(re-ri)/dble(nr)
     do i=0,nang
      xi_el(i)=-pi*0.25d0+dang*dble(i)
      eta_el(i)=-pi*0.25d0+dang*dble(i)
     enddo
     do i=0,nr
      r_el(i)=ri+dr*dble(i)
     enddo
     allocate(x_el(0:nang,0:nang,0:nr,1:6))
     allocate(y_el(0:nang,0:nang,0:nr,1:6))
     allocate(z_el(0:nang,0:nang,0:nr,1:6))
     x_el = 0.d0; y_el=0.d0; z_el = 0.d0
     allocate(number_el(0:nang,0:nang,0:nr,1:6))
     number_el = 0
     do izone = 1, 6
        do iii=0,nr-1     ! loop over  r
            do ii=0,nang-1   ! loop over eta
                do i=0,nang-1   ! loop over ksi
                 number_el(i,ii,iii,izone) = (izone-1)*(nang*nang*nr)+&
                                       iii*(nang*nang)+((ii*nang)+i+1)
                end do
            end do
        end do
     end do
     nelt = maxval(number_el)
     nrpol=1
     npol_cs=1
     allocate(xi_i(0:npol_cs),xi_j(0:npol_cs),xi_k(0:nrpol))

     xi_i(0) = -1
     xi_i(1) = 1
     xi_j(0) = -1
     xi_j(1) = 1
     xi_k(0) = -1
     xi_k(1) = 1

     write(6,*)'Now define collocation points for the SEM-cubed sphere grid'
     allocate(xcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     allocate(ycol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     allocate(zcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     kpts=0
     nel_surf=0

     write(6,*)'ZONE NANG:',nr,nang,npol_cs,nelt
     write(6,*)'ZONE NANG 2 nel_surf:',6*nr*nang**2
     
     do izone=1,6         ! loop over the chunks
      do iii=0,nr-1       ! loop over r   (the global r)
       do ii=0,nang-1     ! loop over eta (the global eta)
        do i=0,nang-1     ! loop over xi  (the global xi)
         iel = number_el(i,ii,iii,izone)
         nel_surf=nel_surf+1
         do kpol = 0, 0 ! TNM nrpol  ! loop over the elemental k index (r direction)
          do jpol = 0, npol_cs  ! loop over the elemental j index (eta direction)
           do ipol = 0, npol_cs ! loop over the elemental i index (xi direction)
            tksi= xi_el(  i) +(1.d0+ xi_i(ipol))*.5d0*dang
            teta=eta_el( ii) +(1.d0 + xi_j(jpol))*.5d0*dang
            tr=    r_el(iii) +(1.d0 + xi_k(kpol))*.5d0*dr
            Xe=tan(tksi)
            Ye=tan(teta)
            C=(1.d0+Xe**2)**(0.5d0)
            D=(1.d0+Ye**2)**(0.5d0)
            delta=1.d0+Xe**2+Ye**2
            kpts=kpts+1
            if(izone==1) then
             Xcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if(izone==2) then
             Xcol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if(izone==3) then
             Xcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if(izone==4) then
             Xcol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if(izone==5) then
             Xcol(ipol,jpol,kpol,iel)=-tr*Ye*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
            endif
            if(izone==6) then
             Xcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
            endif
           end do
          end do
         end do
        end do
       end do
      end do
     end do
!     At this stage we know Xcol, Ycol, Zcol for the cubed sphere
!     These arrays are four dimensional (ipol,jpol,kpol,iel)
!     Their knowledge suffice to define the vtk output that will properly take
!     into account the cubed sphere topology (see the paraview.f90 module)
!!!!! END CUBED SPHERE

      write(6,*)'number of surface pts,elems,tot pts:',npts_surf,nel_surf,kpts
      write(6,*)'max ind_proc, ind_pts:',maxval(ind_proc_surf),maxval(ind_pts_surf)
      write(6,*)size(xsurf)
      xsurf=0.d0; ysurf=0.;zsurf=0.d0
      iii=0
      write(6,*)'constructing 1d array for surface coordinates...'
      do iel=1,nel_surf
         if ( mod(iel,floor(nel_surf/10.))==0  ) then
            write(6,*)'percentage done:',ceiling(real(iel)/real(nel_surf)*100.)
            call flush(6)
         endif
         xc=sum(xcol(:,:,0,iel))/4.d0
         yc=sum(ycol(:,:,0,iel))/4.d0
         zc=sum(zcol(:,:,0,iel))/4.d0
         
         call xyz2rthetaphi(r,th,ph,xc,yc,zc)

         if ( (in_or_out=='innside' .and. ph>=phi0 .and. ph<=phi0+dphi) .or. &
               (in_or_out=='outside' .and. (ph<=phi0 .or. ph>=phi0+dphi) ) ) then

            xsurf(iii+1)=xcol(0,0,0,iel)
            ysurf(iii+1)=ycol(0,0,0,iel)
            zsurf(iii+1)=zcol(0,0,0,iel)
            
            xsurf(iii+2)=xcol(npol_cs,0,0,iel)
            ysurf(iii+2)=ycol(npol_cs,0,0,iel)
            zsurf(iii+2)=zcol(npol_cs,0,0,iel)
            
            xsurf(iii+3)=xcol(npol_cs,npol_cs,0,iel)
            ysurf(iii+3)=ycol(npol_cs,npol_cs,0,iel)
            zsurf(iii+3)=zcol(npol_cs,npol_cs,0,iel)
            
            xsurf(iii+4)=xcol(0,npol_cs,0,iel)
            ysurf(iii+4)=ycol(0,npol_cs,0,iel)
            zsurf(iii+4)=zcol(0,npol_cs,0,iel)
            
            ! determine the corresponding point in the D-shape domain
            do jj=1,4
               call xyz2rthetaphi(r,th, azi_phi_surf(iii+jj),xsurf(iii+jj),ysurf(iii+jj),zsurf(iii+jj))
               dist=2.d0*pi
               do i=1,npts_surf
                  call get_r_theta(coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),1), & 
                       coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),2),r_ref,th_ref)
                  if (abs(th-th_ref)< dist) then 
                     dist=abs(th-th_ref)
                     azi_ind_surf(iii+jj)=ind_proc_surf(i)*npts+ind_pts_surf(i)
                  endif
               enddo
            enddo

            iii=iii+4
            
            endif ! in_or_out

         enddo

         kpts=iii
         write(6,*)'total points in surface:',kpts
         
      write(6,*)'done with cubed sphere for r=',rsurf
      deallocate(xcol,ycol,zcol,x_el,y_el,z_el)

end subroutine construct_surface_cubed_sphere
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine sphi2xy(x,y,s,phi,n)
  integer, intent(in) :: n
  real, dimension(1:n), intent(out) :: x,y
  real, dimension(1:n), intent(in) :: s
  real, intent(in) :: phi
 
  x(1:n)=s(1:n)*cos(phi)
  y(1:n)=s(1:n)*sin(phi)

end subroutine sphi2xy
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal(x,y,z,u1,rows,nelem_disk,filename1)
  use data_all
  implicit none
  integer :: t,rows,nelem_disk,ioerr,dims
  real, dimension(1:rows), intent(in) :: x,y,z,u1
  integer, dimension(1:rows,2) :: cell
  integer, dimension(1:rows) :: cell_type
  real, dimension(1:rows,3) :: W
  
  character (len=30) :: celltype;
  character (len=100) :: filename1;
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
  
  ! 1 IS WRONG FOR OUTDIR !!!!! 
  open(100,file=trim(filename1)//'.vtk',access='stream',&
                          status='replace',convert='big_endian')
  
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
  write(6,*)'...saved ',trim(outdir(1))//'/'//trim(filename1)//'.vtk'
end subroutine write_vtk_bin_scal
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal_topology(x,y,z,u1,elems,filename)
  implicit none
  integer*4 :: i,t,elems
  real*4, dimension(1:elems*4), intent(in) :: x,y,z,u1
  integer*4, dimension(1:elems*5) :: cell
  integer*4, dimension(1:elems) :: cell_type
  character (len=200) :: filename
  character (len=50) :: ss; !stream
  !points structure
  
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
end subroutine write_vtk_bin_scal_topology
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine rthetaphi2xyz(x,y,z,r,theta,phi)
  real, intent(out) :: x,y,z
  real, intent(in) :: r,theta,phi

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta) 

end subroutine rthetaphi2xyz
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
real function prem(r0,param)
  ! prem model in terms of domains separated by discontinuities
  ! MvD: - at discontinuities, upper or lower domain is chosen based on numerical
  !        rounding errors
  !      - only used for mesh plots, nor physical meaning
  
  
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
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine xyz2rthetaphi(r,theta,phi,x,y,z)
  use global_par
  ! Alex2TNM: 
  ! This one you might find useful for other purposes
  ! TNM@Alex: Indeed, done.
  double precision, intent(out) :: r,theta,phi
  double precision, intent(in) :: x,y,z

  r = dsqrt(x**2+y**2+z**2)
  theta = .5d0*pi-dasin(z/r)
  if (y>0) then
    if (x>0) then
      phi = datan(y/(x+1.e-20))
    else
      phi = pi+datan(y/(x+1.e-20))
    endif
  else
    if (x>0) then
      phi = 2*pi+datan(y/(x+1.e-20))
    else
      phi = pi+datan(y/(x+1.e-20))
    end if
  end if
   if (abs(x)<1.e-20) then
     if(y>0.d0) then
       phi = .5d0*pi
     else
       phi = 1.5d0*pi
    end if
  endif
end subroutine xyz2rthetaphi
!--------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine get_r_theta(s,z,r,th)
  use global_par
  double precision, intent(in) :: s,z
  double precision, intent(out) :: r,th
 
  th=datan(s/(z+epsi))
 
  if ( 0.d0 > th ) th = pi + th
  if (th == zero .and. z < 0.d0) th = pi
 
  r=dsqrt(s**2 + z**2)

end subroutine get_r_theta
!==========================================================================
