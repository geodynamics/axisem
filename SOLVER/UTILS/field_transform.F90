!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
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

!=========================================================================================
program field_transformation

#ifdef unc
    use netcdf
#endif
    implicit none

#ifdef unc
    include 'netcdf.inc'
#endif

#ifdef unc
    integer                         :: nvar, ivar
    integer                         :: nsnap, ntimes, ngll, ngllread
    integer                         :: nmode

    integer                         :: ncin_id, ncin_snap_grpid
    integer                         :: ncout_id, ncout_fields_grpid, ncout_gll_dimid
    integer                         :: ncout_freq_dimid, ncout_snap_dimid
    character(len=8)                :: sourcetype
    integer                         :: dimids(2)
    character(len=16), allocatable  :: varnamelist(:)
    integer, dimension(9)           :: ncin_field_varid
    integer, dimension(9,2)         :: ncout_field_varid
    integer                         :: ncin_snaptime_varid

    integer                         :: rank, istride, ostride, nomega, nextpow2
    integer(kind=8)                 :: plan_fftf
    integer                         :: iret
    integer                         :: nstep

    real(kind=8), dimension(:,:), allocatable       :: datat, datat_t
    complex(kind=8), dimension(:,:), allocatable    :: dataf

    double precision                :: time_fft, time_i, time_o, tick, tack
    double precision                :: space_i, space_o

    logical                         :: verbose = .true.
    logical                         :: output_exist = .true.
    integer                         :: npointsperstep = 20000
                                    !< maybe replace this with cache size

    ! lossy compression: setting the fields to numerical zero before p-wave
    ! arrival helps bzip a lot for the compression. The fields are truncated to
    ! the significant digits below the maximum of each time trace.
    ! the idea goes back to point 8) in
    ! http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    ! lossy deflation does not make sense for ffted fields (no blocks of small
    ! numbers expected)
    logical, parameter              :: deflate = .true.
    integer, parameter              :: deflate_level = 1
    logical, parameter              :: deflate_lossy = .true.
    integer, parameter              :: sigdigits =  5       ! significant digits
                                                            ! below max of time trace

    ! initialize timer
    time_fft = 0
    time_i = 0
    time_o = 0

    space_i = 0
    space_o = 0

    ! open input netcdf file 
    call check( nf90_open(path="./Data/axisem_output.nc4", & 
                          mode=NF90_NOWRITE, ncid=ncin_id) )

    ! get Snapshots group id
    call check( nf90_inq_grp_ncid(ncin_id, "Snapshots", ncin_snap_grpid) )

    ! get excitation type (monopole or multipole?)
    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "excitation type", sourcetype))

    if (verbose) &
        print *, sourcetype
    
    if (sourcetype=='monopole')  then
        nvar = 6
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                        'straintrace', 'velo_s     ', 'velo_z     '/)
    else
        nvar = 9
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                        'strain_dsup', 'strain_dzup', 'straintrace', &
                        'velo_s     ', 'velo_p     ', 'velo_z     '/)
    end if

    ! get variable ids of the fields
    do ivar=1, nvar
        call check( nf90_inq_varid(ncin_snap_grpid, varnamelist(ivar), &
                                   varid=ncin_field_varid(ivar)) )
    end do

    ! get dimension ids (same for all fields)
    ivar = 1
    call check( nf90_inquire_variable(ncin_snap_grpid, ncin_field_varid(ivar), &
                                      dimids=dimids(:)) )

    ! get number of snapshots (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(2), len=nsnap) )
    if (verbose) &
        print *, 'nsnap = ', nsnap

    ! get number of gll points (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(1), len=ngll) )
    if (verbose) &
        print *, 'ngll  = ', ngll

    !! Check for output file
    inquire(file="./ordered_output.nc4", exist=output_exist)

    if (output_exist) then
        if (verbose) &
            print *, 'Opening output file'
        nmode = NF90_WRITE 
        call check( nf90_open(path="./ordered_output.nc4", mode=nmode, ncid=ncout_id))
        call check( nf90_redef(ncid=ncout_id))
    else
        !! Create output file
        if (verbose) &
            print *, 'Creating output file'
        nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
        call check( nf90_create(path="./ordered_output.nc4", cmode=nmode, ncid=ncout_id))
    end if

    ! create group for timedomain fields
    call check( nf90_def_grp(ncout_id, "Snapshots", ncout_fields_grpid) )
    print *, 'Defined group "Snapshots"'

    ! copy attributes
    call check( nf90_copy_att( ncin_snap_grpid, NF90_GLOBAL, 'nstrain', &
                               ncout_fields_grpid, NF90_GLOBAL) )
    print *, 'Copyied Attributes'

    ! create dimensions
    call check( nf90_def_dim(ncid=ncout_fields_grpid, name="gllpoints_all", &
                             len=ngll, dimid=ncout_gll_dimid) )

    call check( nf90_def_dim(ncid=ncout_fields_grpid, name="snapshots", len=nsnap, &
                             dimid=ncout_snap_dimid) )

    print *, 'Defined snapshots dimension'


    ! create variables
    do ivar=1, nvar
        call check( nf90_def_var(ncid=ncout_fields_grpid, &
                                 name=trim(varnamelist(ivar)), &
                                 xtype=NF90_FLOAT, &
                                 dimids=(/ncout_snap_dimid, ncout_gll_dimid/),&
                                 varid=ncout_field_varid(ivar, 1), &
                                 chunksizes = (/nsnap, 1/)) )
                                 !chunksizes = (/1, ngll/)) )

        call check( nf90_def_var_fill(ncid=ncout_fields_grpid, &
                                      varid=ncout_field_varid(ivar, 1), &
                                      no_fill=1, fill=0) )


        if (deflate) then
            call check( nf90_def_var_deflate(ncid=ncout_fields_grpid, &
                                             varid=ncout_field_varid(ivar, 1), &
                                             shuffle=1, deflate=1, &
                                             deflate_level=deflate_level) )
        end if
    end do

    call check( nf90_enddef(ncout_id))

    allocate(datat(1:nsnap, 1:npointsperstep))
    allocate(datat_t(1:npointsperstep, 1:nsnap))
   
    ! loop over fields
    do ivar=1, nvar
        if (verbose) &
            print *, varnamelist(ivar)
        ! loop over subsets of the gll points
        nstep = 0
        do while (nstep + 1 < ngll)

            ngllread = min(npointsperstep, ngll - nstep)

            ! read a chunk of data
            call cpu_time(tick)
            call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(1), &
                                     values=datat_t(:,1:nsnap), &
                                     start=(/nstep+1, 1/), &
                                     count=(/ngllread, nsnap/)) )

            datat(1:nsnap,:) = transpose(datat_t(:,1:nsnap))

            call cpu_time(tack)
            time_i = time_i + tack - tick
            space_i = space_i + ngllread * nsnap * 4 / 1048576.
            if (verbose) &
                print "('read  ', F8.2, ' MB in ', F5.2, ' s => ', F6.2, 'MB/s' )", &
                    ngllread * nsnap * 4 / 1048576., tack-tick, &
                    ngllread * nsnap * 4 / 1048576. / (tack-tick)

            ! trunkate for better compression
            if (deflate .and. deflate_lossy) then
                call cpu_time(tick)
                call truncate(datat, sigdigits)
                call cpu_time(tack)
                time_fft = time_fft + tack - tick
            endif

            ! write transposed data to output file
            call cpu_time(tick)
            call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar, 1), &
                                     values=datat, &
                                     start=(/1, nstep+1/), &
                                     count=(/nsnap, ngllread/)) ) 
            call cpu_time(tack)
            time_o = time_o + tack - tick
            space_o = space_o + ngllread * nsnap * 4 / 1048576.
            if (verbose) &
                print "('wrote ', F8.2, ' MB in ', F4.1, ' s => ', F6.2, 'MB/s' )", &
                    ngllread * nsnap * 4 / 1048576., tack-tick, &
                    ngllread * nsnap * 4 / 1048576. / (tack-tick)
                !end if !dofft

            nstep = nstep + npointsperstep
        end do
    enddo

    call check( nf90_close(ncin_id))
    call check( nf90_close(ncout_id))

    print *, 'Time spent for compression/fft: ', time_fft
    print *, 'Time spent for I:               ', time_i
    print *, 'Time spent for O:               ', time_o
    print *, ''
    print *, 'MB I: ', space_i, ' MB, av. speed: ', space_i / time_i, 'MB/s'
    print *, 'MB O: ', space_o, ' MB, av. speed: ', space_o / time_o, 'MB/s'

#else

    stop 'This program can only be run with NetCDF enabled'
#endif
    
contains

!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    implicit none
    integer, intent ( in) :: status !< Error code
    integer, allocatable  :: test(:)
#ifdef unc
    if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        test = test - 100
        stop 0

    end if
#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine truncate(dataIO, sigdigits)
    implicit none
    real(kind=8), intent(inout)     :: dataIO(:,:)
    integer, intent(in)             :: sigdigits
    integer                         :: bits
    double precision                :: scaleit, maxt
    integer                         :: n

    ! find next power of 2 corresponding to sigdigits
    bits = 1
    do while (log10(2.) * bits < real(sigdigits))
        bits = bits + 1
        scaleit = real(2**bits)
    end do

    ! trunkate the data time series wise (each gll point separately)
    do n = lbound(dataIO,2), ubound(dataIO,2)
        maxt = maxval(dataIO(:,n))
        dataIO(:,n) = real(nint(scaleit * dataIO(:,n) / maxt) / scaleit) * maxt
    end do

end subroutine
!-----------------------------------------------------------------------------------------

end program
!=========================================================================================
