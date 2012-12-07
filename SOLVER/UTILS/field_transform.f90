program field_transformation

    use netcdf
    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    integer                         :: nvar, ivar
    integer                         :: nsnap, ntimes, ngll
    integer                         :: nmode, nthreads

    integer                         :: ncin_id, ncin_snap_grpid, ncin_snap_dimid, ncin_proc_dimid 
    integer                         :: ncout_id, ncout_fields_grpid, ncout_gll_dimid
    integer                         :: ncout_freq_dimid
    character(len=8)                :: sourcetype
    integer                         :: dimids(2)
    character(len=16), allocatable  :: varnamelist(:)
    integer, dimension(9)           :: ncin_field_varid, ncout_field_dimid
    integer, dimension(9,2)         :: ncout_field_varid
    logical, parameter              :: deflate = .false.
    integer, parameter              :: deflate_level = 0

    integer                         :: rank, istride, ostride, nomega, nextpow2
    integer*8                       :: plan_fftf
    integer                         :: iret
    integer                         :: npointsperstep, istep, nstep

    real(kind=8), dimension(:,:), allocatable       :: datat, datat_t
    complex(kind=8), dimension(:,:), allocatable    :: dataf

    double precision                :: time_fft, time_i, time_o, tick, tack

    npointsperstep = 10000
    nthreads = 4

    ! initialize timer
    time_fft = 0
    time_i = 0
    time_o = 0

    ! initialize multithreading
    if (nthreads > 1) then
        call dfftw_init_threads(iret)
        if (iret /= 1) then
            print *, 'iret = ', iret
            print *, 'WARNING: Problem with initialization of multithreading for fft'
            print *, '         Continuing serially'
        else
            print *, 'setting up with ', nthreads, ' threads'
            call dfftw_plan_with_nthreads(nthreads)
        endif
    endif

    ! open input netcdf file 
    call check( nf90_open(path="../bla_4_v2/Data/axisem_output.nc4", & 
                          mode=NF90_NOWRITE, ncid=ncin_id) )

    ! get Snapshots group id
    call check( nf90_inq_grp_ncid(ncin_id, "Snapshots", ncin_snap_grpid) )

    ! get excitation type (monopole or multipole?)
    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "excitation type", sourcetype))

    print *, sourcetype
    
    if (sourcetype=='monopole')  then
        nvar = 2
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol'/)
        !nvar = 6
        !allocate(varnamelist(nvar))
        !varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
        !                'straintrace', 'velo_s     ', 'velo_z     '/)
    else
        nvar = 9
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                        'strain_dsup', 'strain_dzup', 'straintrace', &
                        'velo_s     ', 'velo_p     ', 'velo_z     '/)

        print *, 'who cares about multipoles?'
        stop
    end if

    ! get variable ids of the fields
    do ivar=1, nvar
        print *, varnamelist(ivar)
        call check( nf90_inq_varid(ncin_snap_grpid, varnamelist(ivar), &
                                   varid=ncin_field_varid(ivar)) )
    end do

    ! get dimension ids (same for all fields)
    ivar = 1
    call check( nf90_inquire_variable(ncin_snap_grpid, ncin_field_varid(ivar), &
                                      dimids=dimids(:)) )

    ! get number of snapshots (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(2), len=nsnap) )
    print *, 'nsnap = ', nsnap

    ! get number of gll points (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(1), len=ngll) )
    print *, 'ngll  = ', ngll

    ! compute optimal length for fft
    nextpow2 = 2
    do while (nextpow2 < nsnap) 
        nextpow2 = nextpow2 * 2
    end do

    nomega = nextpow2 + 1
    ntimes = nextpow2 * 2
    print *, 'nomega = ', nomega
    print *, 'ntimes = ', ntimes


    !! Create output file
    print *, 'Creating output file'
    nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
    call check( nf90_create(path="ordered_output.nc4", cmode=nmode, ncid=ncout_id))

    ! create group for freqdomain fields
    call check( nf90_def_grp(ncout_id, "freqdomain_fields", ncout_fields_grpid) )

    ! create dimensions for freqdomain fields
    call check( nf90_def_dim(ncid=ncout_fields_grpid, name="gllpoints", len=ngll, &
                             dimid=ncout_gll_dimid) )

    call check( nf90_def_dim(ncid=ncout_fields_grpid, name="omega", len=nomega, &
                             dimid=ncout_freq_dimid) )

    ! create variables for real and imaginary part of freq domain fields
    do ivar=1, nvar
        call check( nf90_def_var(ncid=ncout_fields_grpid, name=trim(varnamelist(ivar))//'_real', &
                                 xtype=NF90_FLOAT, &
                                 dimids=(/ncout_freq_dimid, ncout_gll_dimid/),&
                                 varid=ncout_field_varid(ivar, 1), &
                                 chunksizes = (/nomega, npointsperstep/)) )
                                 !chunksizes = (/nomega, 1/)) )

        call check( nf90_def_var(ncid=ncout_fields_grpid, name=trim(varnamelist(ivar))//'_imag', &
                                 xtype=NF90_FLOAT, &
                                 dimids=(/ncout_freq_dimid, ncout_gll_dimid/),&
                                 varid=ncout_field_varid(ivar, 2), &
                                 chunksizes = (/nomega, npointsperstep/)) )
                                 !chunksizes = (/nomega, 1/)) )
       
    end do
    call check( nf90_enddef(ncout_id))

    rank = 1
    istride = 1
    ostride = 1
    print *, 'ntimes = ',  ntimes
    print *, 'nomega = ',  nomega

    ! allocate working arrays for fourier transform
    allocate(datat_t(1:npointsperstep, 1:ntimes))
    allocate(datat(1:ntimes, 1:npointsperstep))
    allocate(dataf(1:nomega, 1:npointsperstep))

    ! generate plan for fft
    call dfftw_plan_many_dft_r2c(plan_fftf, rank, ntimes, npointsperstep, datat, npointsperstep, istride, &
                                 ntimes, dataf, npointsperstep, ostride, nomega, FFTW_ESTIMATE)

    ! loop over fields
    do ivar=1, nvar
        ! loop over subsets of the gll points
        nstep = 0
        do while (nstep + 1 + npointsperstep < ngll)
            !initialize to zero for padding
            datat = 0.
            dataf = 0.
            ! read a chunk of data
            call cpu_time(tick)
            call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(1), values=datat_t(:,1:nsnap), &
                                     start=(/nstep+1, 1/), count=(/npointsperstep, nsnap/)) )!, &
                                     !map=(/nsnap, 1/)) ) 
            datat(1:nsnap,:) = transpose(datat_t(:,1:nsnap))
            call cpu_time(tack)
            time_i = time_i + tack - tick
            print "('read  ', F8.2, ' MB in ', F4.1, ' s => ', F6.2, 'MB/s' )", &
                npointsperstep * nsnap * 4 / 1048576., tack-tick, &
                npointsperstep * nsnap * 4 / 1048576. / (tack-tick)

            ! ADD TAPERING HERE

            ! do fft
            call cpu_time(tick)
            call dfftw_execute(plan_fftf)
            call cpu_time(tack)
            time_fft = time_fft + tack - tick

            ! write real and imaginary parts to output file
            call cpu_time(tick)
            call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar, 1), values=realpart(dataf), &
                                     start=(/1, nstep+1/), count=(/nomega, npointsperstep/)) ) 

            call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar, 2), values=imagpart(dataf), &
                                     start=(/1, nstep+1/), count=(/nomega, npointsperstep/)) ) 
            call cpu_time(tack)
            time_o = time_o + tack - tick
            print "('wrote ', F8.2, ' MB in ', F4.1, ' s => ', F6.2, 'MB/s' )", &
                npointsperstep * nomega * 2 * 4 / 1048576., tack-tick, &
                npointsperstep * nomega * 2 * 4 / 1048576. / (tack-tick)


            nstep = nstep + npointsperstep
        end do

        ! special treatment of the last chunk (having less gll points)
        ! this should be put in the loop above to avoid doubble coding of
        ! read/write statements
        nstep = nstep - npointsperstep
        datat = 0.
        call cpu_time(tick)
        call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(1), values=datat(1:ngll-nstep, 1:nsnap), &
                                 start=(/nstep+1, 1/), count=(/ngll-nstep, nsnap/)) )!, &
                                 !map=(/nsnap, 1/)) ) 
        datat(1:nsnap,1:ngll-nstep) = transpose(datat_t(1:ngll-nstep,1:nsnap))
        call cpu_time(tack)
        time_i = time_i + tack - tick
        print "('read  ', F8.2, ' MB in ', F4.1, ' s => ', F6.2, 'MB/s' )", &
            (ngll - nstep) * nsnap * 4 / 1048576., tack-tick, &
            (ngll - nstep) * nsnap * 4 / 1048576. / (tack-tick)

        ! ADD TAPERING HERE

        call cpu_time(tick)
        call dfftw_execute(plan_fftf)
        call cpu_time(tack)
        time_fft = time_fft + tack - tick

        ! write real and imaginary parts to output file
        call cpu_time(tick)
        ! MVD: npointsperstep should be too much, why no error???
        call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar, 1), values=realpart(dataf), &
                                 start=(/1, nstep+1/), count=(/nomega, npointsperstep/)) ) 

        call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar, 2), values=imagpart(dataf), &
                                 start=(/1, nstep+1/), count=(/nomega, npointsperstep/)) ) 
        call cpu_time(tack)
        time_o = time_o + tack - tick
        print "('wrote ', F8.2, ' MB in ', F4.1, ' s => ', F6.2, 'MB/s' )", &
            (ngll - nstep) * nsnap * 4 / 1048576., tack-tick, &
            (ngll - nstep) * nsnap * 4 / 1048576. / (tack-tick)

    enddo

    call dfftw_destroy_plan(plan_fftf)
    call check( nf90_close(ncin_id))
    call check( nf90_close(ncout_id))

    print *, 'Time spent for FFT: ', time_fft
    print *, 'Time spent for I:   ', time_i
    print *, 'Time spent for O:   ', time_o
    
    contains
!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    implicit none
    integer, intent ( in) :: status !< Error code
    if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
    end if
end subroutine check  
!-----------------------------------------------------------------------------------------


end program
