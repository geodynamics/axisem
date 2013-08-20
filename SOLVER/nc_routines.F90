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

!> Contains all the routines for NetCDF handling.
module nc_routines

#ifdef unc
    use netcdf
#endif
    use data_io,    only : verbose, deflate_level, nseismo
    use data_proc,  only : mynum, nproc, lpr
    use global_parameters
    use commun,     only : barrier

    implicit none
    save
    private 

    !> Buffer variable for recorder 
    real(sp), allocatable   :: recdumpvar(:,:,:)
    !> Buffer variable for displacement at surface
    real(sp), allocatable   :: surfdumpvar_disp(:,:,:)
    !> Buffer variable for velocity at surface 
    real(sp), allocatable   :: surfdumpvar_velo(:,:,:)
    !> Buffer variable for strain at surface
    real(sp), allocatable   :: surfdumpvar_strain(:,:,:)
    !> Buffer variable for source displacement at surface
    real(sp), allocatable   :: surfdumpvar_srcdisp(:,:,:)
!    real, allocatable   :: oneddumpvar_flu(:,:,:)    !< Buffer variable for everything fluid dumped in nc_dump_field_1d
!    real, allocatable   :: oneddumpvar_sol(:,:,:)    !< Buffer variable for everything solid dumped in nc_dump_field_1d

    !> Buffer variable for everything dumped in nc_dump_field_1d
    real(sp), allocatable   :: oneddumpvar(:,:,:)
    real(sp), allocatable   :: scoord1d(:), zcoord1d(:)
    
    !> Number of steps before kernel specific stuff is dumped 
    integer             :: dumpstepsnap
    !> Number of GLL points per element 
    integer             :: gllperelem
    !> When is this processor supposed to dump. 
    integer             :: outputplan
    !> How many steps since last dump?
    integer             :: stepstodump
    !> Global variables, so that we do not have to pass data to the C subroutine
    integer             :: isnap_global
    !> dito
    integer             :: ndumps
    !> Will any processor dump at this value of isnap?
    logical,allocatable :: dumpposition(:)
    !> Number of GLL points to plot for this processor
    integer             :: npoints
    !> Number of GLL points to plot for all processors
    integer             :: npoints_global
    !> Number of GLL points to plot in solid/fluid domain
    integer             :: npts_sol, npts_flu
    !> Number of GLL points to plot in solid domain for all processors
    integer             :: npts_sol_global
    !> Number of GLL points to plot in fluid domain for all processors
    integer             :: npts_flu_global

    ! Stuff moved from data_io
    integer             :: ncid_out, ncid_recout, ncid_snapout, ncid_surfout, ncid_meshout
    integer             :: nc_snap_dimid, nc_proc_dimid, nc_rec_dimid, nc_recproc_dimid
    integer             :: nc_times_dimid, nc_comp_dimid, nc_disp_varid, nc_stf_seis_varid
    integer             :: nc_time_varid, nc_iter_dimid, nc_stf_iter_varid

    integer             :: nc_strcomp_dimid
    integer             :: nc_surfelem_disp_varid, nc_surfelem_velo_varid
    integer             :: nc_surfelem_strain_varid, nc_surfelem_disp_src_varid
    integer             :: nc_mesh_sol_varid, nc_mesh_flu_varid
    integer             :: nc_point_dimid, nc_pt_sol_dimid, nc_pt_flu_dimid
    integer             :: nc_szcoord_dimid
    integer             :: nc_snaptime_varid, nc_elem_dom_varid, nc_surfelem_theta_varid
    integer,allocatable :: nc_field_varid(:)
    character(len=16), allocatable  :: varnamelist(:)
    character(len=12), allocatable  :: nc_varnamelist(:)
    integer             :: nvar = -1

    !! Variables for dumping of wavefields for plotting purposes
    integer             :: nc_snap_disp_varid, nc_coord_dimid
    integer             :: nc_snap_point_varid, nc_snap_grid_varid
    integer             :: ncid_out_snap
    integer             :: ndim_disp !< 2 for monopole, 3 for rest


    !! Buffer variables for the STF.
    real(kind=sp), allocatable :: stf_seis_dumpvar(:)
    real(kind=sp), allocatable :: stf_dumpvar(:)

    !! @todo These parameters should move to a input file soon
    !> How many snaps should be buffered in RAM?
    integer             :: dumpbuffersize = 256
    
    public              :: nc_dump_strain, nc_dump_rec, nc_dump_surface
    public              :: nc_dump_field_solid, nc_dump_field_fluid
    public              :: nc_write_att_char, nc_write_att_real, nc_write_att_int
    public              :: nc_define_outputfile, nc_finish_prepare, nc_end_output
    public              :: nc_dump_strain_to_disk, nc_dump_mesh_sol, nc_dump_mesh_flu
    public              :: nc_write_el_domains
    public              :: nc_dump_snapshot, nc_dump_snap_points, nc_dump_snap_grid
    public              :: nc_make_snapfile, nc_dump_stf
contains


!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    integer, intent ( in) :: status !< Error code
#ifdef unc
    if (status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        call abort()
    end if
#endif
end subroutine check  
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Write out the model domains for each element
subroutine nc_write_el_domains(idom)
    integer, dimension(:), intent(in) :: idom

#ifdef unc
    call check( nf90_put_var(ncid   = ncid_snapout, &
                             varid  = nc_elem_dom_varid, &
                             values = idom) )
#endif
end subroutine

!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
subroutine nc_dump_field_solid(f, varname)

    real(kind=realkind), intent(in)   :: f(:)     !< Data to dump, size should be npts_sol
    character(len=*), intent(in)      :: varname  !< Internal name of data to dump. 
    !! Is used to identify the NetCDF Variable in question
#ifdef unc
    integer                           :: ivar

    do ivar=1, nvar
        !< Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    end do

    if (ivar > nvar/2) then
        write(6,*) 'nc_dump_field_solid: Trying to access variable: ', trim(varname), &
            ' which is a fluid variable. Contact a developer and shout at him!'
        stop 1
    end if
    
    oneddumpvar(1:npts_sol,stepstodump+1,ivar) = f  !processor specific dump variable
#endif
end subroutine nc_dump_field_solid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
subroutine nc_dump_field_fluid(f, varname)

    real(kind=realkind), intent(in)   :: f(:)  !< Data to dump, size should be npts_flu
    character(len=*), intent(in)      :: varname  !< Internal name of data to dump. 
    !! Is used to identify the NetCDF Variable in question
#ifdef unc
    integer                           :: ivar


    do ivar=1, nvar
        !< Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    end do

    if (ivar <= nvar/2) then
        write(6,*) 'nc_dump_field_fluid: Trying to access variable: ', trim(varname), &
            ' which is a solid variable. Contact a developer and shout at him!'
        stop 1
    end if

    oneddumpvar(npts_sol+1:npoints,stepstodump+1,ivar-nvar/2) = f  
#endif
end subroutine nc_dump_field_fluid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain(isnap_loc)

    use data_io, ONLY    : nstrain
    use clocks_mod, only : tick
    use data_time, only  : iclocknbio, idnbio

    ! explicit interfaces to the c functions to avoide the underscore issues
    ! (fortran 2003 standard)
    interface
        subroutine c_spawn_dumpthread(stepstodump) bind(c, name='c_spawn_dumpthread')
            use, intrinsic      :: iso_c_binding, only : c_int
            integer (c_int)     :: stepstodump
        end subroutine 
    end interface
    
    interface
        subroutine c_wait_for_io() bind(c, name='c_wait_for_io')
        end subroutine 
    end interface

    integer, intent(in) :: isnap_loc
#ifdef unc
    integer             :: iproc
    real                :: tickl, tackl
    
    stepstodump = stepstodump + 1
    if (isnap_loc == 0) return

    if (dumpposition(mod(isnap_loc, dumpstepsnap))) then

        ! wait for other processes to finish writing, measure waiting time and
        ! issue warning in case waiting longer then .5 sec
        ! MvD: I am not sure if cpu_time is a correct measure here, as we idle
        !      until IO is finished.
        !      Therefor testing system_clock, from the clocks module.

        if (mynum == 0) call cpu_time(tickl)

        iclocknbio = tick()
        do iproc=0, nproc-1
            if (iproc == mynum) then
                call c_wait_for_io()
                call flush(6)
            end if
            call barrier
        end do
        iclocknbio = tick(id=idnbio, since=iclocknbio)

        if (mynum == 0) then
            call cpu_time(tackl)
            if ((tackl-tickl) > 0.5 .and. verbose > 0) then
                write(6,"('WARNING: Computation was halted for ', F7.2, ' s to wait for ',&
                         & 'dumping processor. Consider adapting netCDF output variables',&
                         & '(disable compression, increase dumpstepsnap)')") tackl-tickl
            end if
        end if

        ! non blocking write
        if (mod(isnap_loc, dumpstepsnap) == outputplan) then 
            call c_wait_for_io()
            isnap_global = isnap_loc 
            ndumps = stepstodump
            call c_spawn_dumpthread(stepstodump)
            stepstodump = 0
        end if
    end if 

    ! Final and last dump of all remaining 
    ! not done in unblocking fashion using a thread, as there is nothing to
    ! compute anymore
    ! Make sure, nobody is accessing the output file anymore
    if (isnap_loc == nstrain) then 
        do iproc=0, nproc-1
            if (iproc == mynum) then
                call c_wait_for_io()
                call barrier
            end if
        end do
        do iproc=0,nproc-1
            if (iproc == mynum) then
                isnap_global = nstrain 
                ndumps = stepstodump
                call nc_dump_strain_to_disk()
                if (verbose > 1) write(6,*) mynum, 'finished dumping strain'
            end if
            call barrier
        end do
    end if

#endif
end subroutine nc_dump_strain
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain_to_disk() bind(c, name="nc_dump_strain_to_disk")
#ifdef unc

    use data_io
    use data_source,       ONLY: src_type
    use global_parameters, ONLY: realkind
    use data_mesh,         ONLY: loc2globrec, maxind

    integer                           :: ivar, flen, isnap_loc
    real                              :: tick, tack
    integer                           :: dumpsize

    dumpsize = 0
    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                          mode=NF90_WRITE, ncid=ncid_out) )
    
    isnap_loc = isnap_global
    if (verbose > 1) then
        if (ndumps == 0) then 
            write(6,"('  Proc ', I4, ' in dump routine, isnap =', I5, &
                    & ', nothing to dump, returning...')") mynum, isnap_loc
            return
        else
            write(6,"('  Proc ', I4, ' in dump routine, isnap =', I5, &
                    & ', stepstodump = ', I4)") mynum, isnap_loc, ndumps
        end if
    end if


    do ivar=1, nvar/2
        call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
                                 start=(/mynum*npoints+1, isnap_loc-ndumps+1/), &
                                 count=(/npoints, ndumps/), &
                                 values=oneddumpvar(1:npoints,1:ndumps,ivar)) )
        dumpsize = dumpsize + npoints * ndumps
    end do
        
    dumpsize = dumpsize + flen * ndumps
            
    !> Surface dumps 
    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_varid, &
                start = (/isnap_loc-ndumps+1, 1, mynum*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_disp(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_velo_varid, &
                start = (/isnap_loc-ndumps+1, 1, mynum*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_velo(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_strain_varid, &
                start = (/isnap_loc-ndumps+1, 1, mynum*maxind+1/), &
                count = (/ndumps, 6, maxind/), &
                values = surfdumpvar_strain(:,:,:)) )
    dumpsize = dumpsize + 6 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_src_varid, &
                start = (/isnap_loc-ndumps+1, 1, mynum*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_srcdisp(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_close(ncid_out) ) 
    call cpu_time(tack)

    if (verbose > 1) then
        write(6,"('  Proc', I5,': Wrote ', F8.2, ' MB in ', F6.2, 's')") &
            mynum, real(dumpsize) * realkind / 1048576., tack-tick 
        call flush(6)
    end if

#endif
end subroutine nc_dump_strain_to_disk
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine nc_dump_stf(stf)
    use data_io,  only                       : nseismo
    use data_time, only                      : seis_it, niter
    real(kind=sp), intent(in), dimension(:) :: stf   
#ifdef unc
    integer                                 :: it, i

    allocate(stf_dumpvar(niter))
    allocate(stf_seis_dumpvar(nseismo))
    stf_seis_dumpvar = 0.0
    it = 1
    do i = 1, niter
        if ( mod(i,seis_it) == 0) stf_seis_dumpvar(it) = stf(i) 
        it = it + 1
    end do
    stf_dumpvar = stf


#endif
end subroutine nc_dump_stf
!----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Dump receiver specific stuff, especially displacement and velocity
!! N.B.: Works with global indices.
subroutine nc_dump_rec(recfield)
    use data_mesh, ONLY: num_rec
    use data_io,   ONLY: iseismo
    real(sp), intent(in), dimension(3,num_rec) :: recfield
#ifdef unc
   
    recdumpvar(iseismo,:,:) = 0.0
    where(abs(recfield)>epsi) recdumpvar(iseismo,:,:) = recfield(:,:)

#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_rec_to_disk
#ifdef unc
    use data_mesh, ONLY: loc2globrec, num_rec
    use data_io,   ONLY: datapath, lfdata, nseismo

    real                              :: tick, tack
    integer                           :: irec, dumpsize, icomp

    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                          mode=NF90_WRITE, ncid=ncid_out) )
    call check( nf90_inq_varid( ncid_recout, "displacement", nc_disp_varid ) )

    dumpsize = 0
    do irec = 1, num_rec
        do icomp = 1, 3
            call check( nf90_put_var(ncid=ncid_recout, varid=nc_disp_varid, &
                                   start=(/1, icomp, loc2globrec(irec)/), &
                                   count = (/nseismo, 1, 1/), values=recdumpvar(:,icomp,irec)) )
        end do
    end do

    dumpsize = nseismo * 3 * num_rec
    call check( nf90_close(ncid=ncid_out))
    call cpu_time(tack)

    if (verbose > 1) then
        write(6,"(I3,': Receiver data, Wrote ', F8.3, ' MB in ', F6.2, 's')") &
            mynum, real(dumpsize) * realkind / 1048576., tack-tick 
    end if
#endif
end subroutine nc_dump_rec_to_disk
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
!> Dump stuff along surface
subroutine nc_dump_surface(surffield, disporvelo)!, nrec, dim2)
    use data_mesh, ONLY: maxind

    !integer, intent(in)                          :: nrec, dim2
    real(kind=realkind), intent(in), dimension(:,:) :: surffield
    character(len=4), intent(in)                 :: disporvelo
#ifdef unc
   
    select case(disporvelo)
        case('disp')
            surfdumpvar_disp(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('velo')
            surfdumpvar_velo(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('stra')
            surfdumpvar_strain(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('srcd')
            surfdumpvar_srcdisp(stepstodump+1,:,:) = transpose(surffield(:,:))
    end select

#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_sol(scoord_sol, zcoord_sol)

    use data_io, ONLY : ndumppts_el
    real(sp), intent(in), dimension(:,:,:) :: scoord_sol, zcoord_sol
#ifdef unc

    npts_sol = size(scoord_sol)
    if (npts_sol.ne.size(scoord_sol)) then
        write(6,*) 'Inconsistency in mesh size in nc_dump_mesh_solid'
        write(6,*) 'npts_sol=', npts_sol, ', size(scoord_sol)=', size(scoord_sol)
        stop 2
    end if

    npoints = ndumppts_el * nelem
    allocate(scoord1d(npoints))
    allocate(zcoord1d(npoints))

    zcoord1d(1:npts_sol) = pack(zcoord_sol,.true.)
    scoord1d(1:npts_sol) = pack(scoord_sol,.true.)


#endif
end subroutine nc_dump_mesh_sol
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_flu(scoord_flu, zcoord_flu)

    real(sp), intent(in), dimension(:,:,:) :: scoord_flu, zcoord_flu 
#ifdef unc

    npts_flu = size(scoord_flu)
    if (npts_flu.ne.size(scoord_flu))  then
        write(6,*) 'Inconsistency in mesh size in nc_dump_mesh_flu'
        stop 2
    end if

    zcoord1d(npts_sol+1:) = pack(zcoord_flu, .true.)
    scoord1d(npts_sol+1:) = pack(scoord_flu, .true.)


#endif
end subroutine nc_dump_mesh_flu
!-----------------------------------------------------------------------------------------
!subroutine nc_dump_mesh_to_disk()
!    use data_io, only : datapath, lfdata
!    integer iproc, nc_mesh_s_varid, nc_mesh_z_varid
!
!    do iproc=0, nproc-1
!        if (mynum == iproc) then
!            call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
!                                  mode=NF90_WRITE, ncid=ncid_out) )
!
!            call check( nf90_inq_varid( ncid_snapout, "mesh S", nc_mesh_s_varid ) )
!            call check( nf90_put_var(ncid_snapout, varid = nc_mesh_s_varid, &
!                                     values = scoord1d, &
!                                     start  = (/1,mynum*npoints+1/), &
!                                     count  = (/1,npoints/) ))
!            call check( nf90_inq_varid( ncid_snapout, "mesh Z", nc_mesh_z_varid ) )
!            call check( nf90_put_var(ncid_snapout, varid = nc_mesh_z_varid, &
!                                     values = zcoord1d, &
!                                     start  = (/1,mynum*npoints+1/), &
!                                     count  = (/1,npoints/) ))
!            
!            call check( nf90_close(ncid_out))
!        end if
!        call barrier
!    end do
!
!
!end subroutine nc_dump_mesh_to_disk

!-----------------------------------------------------------------------------------------
!> Define the output file variables and dimensions
!! and allocate buffer variables.
subroutine nc_define_outputfile(nrec, rec_names, rec_th, rec_th_req, rec_ph, rec_proc)

    use data_io,     ONLY: nseismo, nstrain, nseismo, ibeg, iend, dump_wavefields
    use data_io,     ONLY: datapath, lfdata, strain_samp
    use data_mesh,   ONLY: maxind, num_rec, discont
    use data_source, ONLY: src_type, t_0
    use data_time,   ONLY: deltat, niter


    integer, intent(in)                  :: nrec              !< Number of receivers
    character(len=40),intent(in)         :: rec_names(nrec)   !< Receiver names
    real(dp), dimension(nrec),intent(in) :: rec_th            !< Receiver theta 
    real(dp), dimension(nrec),intent(in) :: rec_th_req        !< Requested receiver theta
    real(dp), dimension(nrec),intent(in) :: rec_ph            !< Receiver phi
    integer, dimension(nrec),intent(in)  :: rec_proc          !< Receiver processor
#ifdef unc
    real(dp), dimension(nstrain)         :: time_strain
    real(dp), dimension(nseismo)         :: time_seis
    character(len=16), allocatable       :: varname(:)
    integer                              :: ivar, i
    integer                              :: irec, iproc, nmode
    integer                              :: nc_latr_varid, nc_lon_varid 
    integer                              :: nc_lat_varid, nc_ph_varid
    integer                              :: nc_thr_varid, nc_th_varid 
    integer                              :: nc_proc_varid, nc_recnam_dimid
    integer                              :: nc_recnam_varid, nc_surf_dimid
    integer                              :: nc_pt_dimid
    integer                              :: nc_mesh_s_varid, nc_mesh_z_varid   
    integer                              :: nc_disc_dimid, nc_disc_varid

    if ((mynum == 0) .and. (verbose > 1)) then
        write(6,*)
        write(6,*) '************************************************************************'
        write(6,*) '**** Producing netcdf output file with its variables and dimensions ****'
        write(6,*) '************************************************************************'
        write(6,*)
    end if

    call barrier
    dumpstepsnap = int(dumpbuffersize / nproc) * nproc ! Will later be reduced to nstrain, if this is smaller
                                                       ! than value given here
    
    if (lpr .and. verbose > 1) write(6,*) '  Dumping NetCDF file to disk every', dumpstepsnap, ' snaps'

    call barrier ! for nicer output only

    outputplan = mynum * (dumpstepsnap / nproc)

    allocate(dumpposition(0:dumpstepsnap-1))
    dumpposition = .false.
    do iproc=0, nproc-1
        dumpposition(iproc*(dumpstepsnap/nproc)) = .true.
        if ((iproc .eq. mynum) .and. (verbose > 1)) then
            write(6,"('Proc ', I4, ' will dump at position ', I4)") mynum, outputplan
            call flush(6)
        end if
        call barrier ! for nicer output only
    end do


    if (mynum == 0) then
        if (verbose > 1) write (6,*) 'Preparing netcdf file for ', nproc, ' processors'
        nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
        call check( nf90_create(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                                cmode=nmode, ncid=ncid_out) )
        if (verbose > 1) write(6,*) 'Netcdf file with ID ', ncid_out, ' produced.'
    end if


    if (dump_wavefields) then
        if (src_type(1) == 'monopole') then
            nvar = 12
        else
            nvar = 18
        end if
        allocate(varname(nvar))
        allocate(varnamelist(nvar))
        allocate(nc_varnamelist(nvar/2))
        allocate(nc_field_varid(nvar/2))

        if (src_type(1)  ==  'monopole') then 
            varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                            'straintrace_sol', 'velo_sol_s     ', 'velo_sol_z     ', &
                            'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                            'straintrace_flu', 'velo_flu_s     ', 'velo_flu_z     '/)
              
            nc_varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                               'straintrace', 'velo_s     ', 'velo_z     '/)
        else
            varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                            'strain_dsup_sol', 'strain_dzup_sol', 'straintrace_sol', &
                            'velo_sol_s     ', 'velo_sol_p     ', 'velo_sol_z     ', &
                            'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                            'strain_dsup_flu', 'strain_dzup_flu', 'straintrace_flu', &
                            'velo_flu_s     ', 'velo_flu_p     ', 'velo_flu_z     '/)
              
            nc_varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                               'strain_dsup', 'strain_dzup', 'straintrace', &
                               'velo_s     ', 'velo_p     ', 'velo_z     '/)
        end if

        gllperelem = (iend - ibeg + 1)**2
        npoints = nelem * gllperelem
        npoints_global = npoints * nproc
        npts_sol = nel_solid * gllperelem
        npts_flu = nel_fluid * gllperelem 
        npts_sol_global = nel_solid * gllperelem * nproc
        npts_flu_global = nel_fluid * gllperelem * nproc
        
        if (nstrain <= dumpstepsnap) dumpstepsnap = nstrain
    end if ! dump_wavefields
    
    if (mynum == 0) then    
        if (verbose > 1) write(6,*) '  Producing groups for Seismograms and Snapshots'

        call check( nf90_def_grp(ncid_out, "Seismograms", ncid_recout) )
        call check( nf90_def_grp(ncid_out, "Snapshots", ncid_snapout) )
        call check( nf90_def_grp(ncid_out, "Surface", ncid_surfout) )
        call check( nf90_def_grp(ncid_out, "Mesh", ncid_meshout) )
        if (verbose > 1) write(6,*) '  Seismograms group has ID', ncid_recout
        if (verbose > 1) write(6,*) '  Snapshots group has ID', ncid_snapout
        if (verbose > 1) write(6,*) '  Surface group has ID', ncid_surfout
        if (verbose > 1) write(6,*) '  Mesh group has ID', ncid_meshout
        
        if (verbose > 1) write(6,*) 'Define dimensions in ''Seismograms'' group of NetCDF output file'

        if (verbose > 1) write(6,*) '  ''Seismograms'' group has ID ', ncid_recout
110     format(' Dimension ', A20, ' with length ', I8, ' and ID', I6) 
        call check( nf90_def_dim(ncid_out, "seis_timesteps", nseismo, nc_times_dimid) )
        if (verbose > 1) write(6,110) "seis_timesteps", nseismo, nc_times_dimid 

        call check( nf90_def_dim(ncid_out, "sim_timesteps", niter, nc_iter_dimid) )
        if (verbose > 1) write(6,110) "sim_timesteps", niter, nc_iter_dimid 

        call check( nf90_def_dim(ncid_recout, "receivers", nrec, nc_rec_dimid) )
        if (verbose > 1) write(6,110) "receivers", nrec, nc_rec_dimid

        call check( nf90_def_dim(ncid_out, "components", 3, nc_comp_dimid) )
        if (verbose > 1) write(6,110) "components", 3, nc_comp_dimid

        call check( nf90_def_dim(ncid_recout, "recnamlength", 40, nc_recnam_dimid) ) 
        if (verbose > 1) write(6,110) "recnamlength", 40, nc_recnam_dimid

        if (verbose > 1) write(6,*) 'NetCDF dimensions defined'

        if (verbose > 1) write(6,*) 'Define variables in ''Seismograms'' group of NetCDF output file'
        call flush(6)

        call check( nf90_def_var(ncid=ncid_recout, name="displacement", xtype=NF90_FLOAT,&
                                 dimids=(/nc_times_dimid, nc_comp_dimid, nc_rec_dimid/), &
                                 !storage = NF90_CHUNKED, chunksizes=(/nseismo,3,1/), &
                                 !contiguous = .false., chunksizes=(/nseismo,3,1/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_disp_varid) )

        call check( nf90_put_att(ncid_recout, nc_disp_varid, 'units', 'meters') )
        call check( nf90_put_att(ncid_recout, nc_disp_varid, '_FillValue', 0.0) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_seis", xtype=NF90_FLOAT,&
                                 dimids=(/nc_times_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_stf_seis_varid) )
        
        call check( nf90_def_var(ncid=ncid_recout, name="stf_iter", xtype=NF90_FLOAT,&
                                 dimids=(/nc_iter_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_stf_iter_varid) )
        
        call check( nf90_def_var(ncid=ncid_recout, name="time", xtype=NF90_DOUBLE,&
                                 dimids=(/nc_times_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_time_varid) )

        !call check( nf90_def_var(ncid_recout, "Lat_req", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_latr_varid) )
        !call check( nf90_def_var(ncid_recout, "Lon", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_lon_varid) )
        !call check( nf90_def_var(ncid_recout, "Lat", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_lat_varid) )
        call check( nf90_def_var(ncid_recout, "azimuths", NF90_FLOAT, (/nc_rec_dimid/),  &
                                 nc_ph_varid) )
        call check( nf90_def_var(ncid_recout, "distances_requested", NF90_FLOAT, &
                                 (/nc_rec_dimid/), nc_thr_varid) )
        call check( nf90_def_var(ncid_recout, "distances", NF90_FLOAT, &
                                 (/nc_rec_dimid/), nc_th_varid) )
        call check( nf90_def_var(ncid_recout, "processor_of_receiver", NF90_INT, &
                                 (/nc_rec_dimid/), nc_proc_varid) )
        call check( nf90_def_var(ncid_recout, "receiver_name", NF90_CHAR, &
                                 (/nc_rec_dimid, nc_recnam_dimid/), nc_recnam_varid) )

        if (dump_wavefields) then
            ! Wavefields group of output file N.B: Snapshots for kernel calculation
            if (verbose > 1) write(6,*) 'Define variables in ''Snapshots'' group of NetCDF output file', &
                                        '  awaiting', nstrain, ' snapshots'

            call check( nf90_def_dim( ncid  = ncid_out, &
                                      name  = 'snapshots', &
                                      len   = nstrain, &
                                      dimid = nc_snap_dimid) )
            call check( nf90_put_att(ncid   = ncid_snapout, &
                                     varid  = NF90_GLOBAL, &
                                     name   = 'nstrain', &
                                     values = nstrain) )
            
            call check( nf90_def_dim( ncid  = ncid_out, &
                                      name  = 'gllpoints_all', &
                                      len   = npoints_global, &
                                      dimid = nc_pt_dimid) )
            call check( nf90_put_att(ncid   = ncid_out, &
                                     varid  = NF90_GLOBAL, &
                                     name   = 'npoints', &
                                     values = npoints_global) )

            call check( nf90_def_var( ncid   = ncid_out, &
                                      name   = 'snapshot_times', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_snap_dimid,&
                                      varid  = nc_snaptime_varid) )

            call check( nf90_def_dim( ncid  = ncid_meshout, &
                                      name  = 'discontinuities', &
                                      len   = ndisc, &
                                      dimid = nc_disc_dimid) )
            call check( nf90_put_att(ncid   = ncid_meshout, &
                                     varid  = NF90_GLOBAL, &
                                     name   = 'ndisc', &
                                     values = ndisc) )

            call check( nf90_def_var( ncid   = ncid_meshout,  &
                                      name   ='mesh_S', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_s_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_Z', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_z_varid) )

            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'model_domain', &
                                      xtype  = NF90_BYTE, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_elem_dom_varid) )
            
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'disc_depths', &
                                      xtype  = NF90_DOUBLE, &
                                      dimids = nc_disc_dimid,&
                                      varid  = nc_disc_varid) )

            do ivar=1, nvar/2 ! The big snapshot variables for the kerner.
       
                call check( nf90_def_var(ncid=ncid_snapout, name=trim(nc_varnamelist(ivar)), &
                                         xtype = NF90_FLOAT, &
                                         dimids = (/nc_pt_dimid, nc_snap_dimid/),&
                                         varid = nc_field_varid(ivar), &
                                         chunksizes = (/npoints_global, 1/) ))
                call check( nf90_def_var_fill(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
                                              no_fill=1, fill=0) )
                if (verbose > 1) write(6,"('Netcdf variable ', A16,' with ID ', I3, ' and length', &
                        & I8, ' produced.')") &
                          trim(nc_varnamelist(ivar)), nc_field_varid(ivar), npoints_global
            end do

            ! Surface group in output file
            if (verbose > 1) write(6,*) 'Define variables in ''Surface'' group of NetCDF output file'
            !call check( nf90_def_dim( ncid_surfout, name="snapshots", len=nstrain, &
            !                          dimid=nc_snap_dimid) )
            call check( nf90_put_att(ncid   = ncid_surfout, &
                                     name   = 'nstrain', &
                                     varid  = NF90_GLOBAL, &
                                     values = nstrain) )
            call check( nf90_def_dim( ncid_surfout, "straincomponents", len=6, &
                                      dimid=nc_strcomp_dimid) )
            if (verbose > 1) write(6,110) "straincomponents", 6, nc_strcomp_dimid
            
            call check( nf90_def_dim( ncid_surfout, "surf_elems", maxind*nproc, nc_surf_dimid) )     
            call check( nf90_put_att(ncid   = ncid_surfout, &
                                     name   = 'nsurfelem', &
                                     varid  = NF90_GLOBAL, &
                                     values = maxind*nproc) )
            if (verbose > 1) write(6,110) "surf_elems", maxind*nproc, nc_surf_dimid

            call check( nf90_def_var( ncid_surfout, "elem_theta", NF90_FLOAT, &
                                      (/nc_surf_dimid /), &
                                      nc_surfelem_theta_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_theta_varid, 'units', 'degrees'))
           
            call check( nf90_def_var( ncid_surfout, "displacement", NF90_FLOAT, &
                                      (/nc_times_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_disp_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_disp_varid, 'units', 'meters'))
            
            call check( nf90_def_var( ncid_surfout, "velocity", NF90_FLOAT, &
                                      (/nc_times_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_velo_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_velo_varid, 'units', &
                                      'meters per second') )
            
            call check( nf90_def_var( ncid_surfout, "strain", NF90_FLOAT, &
                                      (/nc_times_dimid, nc_strcomp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_strain_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_velo_varid, 'units', &
                                      ' ') )

            call check( nf90_def_var( ncid_surfout, "disp_src", NF90_FLOAT, &
                                      (/nc_times_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_disp_src_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_velo_varid, 'units', &
                                      'meters') )
        end if

        
        if (verbose > 1) write(6,'(a/)') 'NetCDF variables defined'
        ! Leave definition mode
        call check( nf90_enddef(ncid_out))

        if (verbose > 1) write(6,*) 'Writing station info into NetCDF file...'
        call check( nf90_put_var( ncid_recout, nc_th_varid, values = rec_th) )
        call check( nf90_put_var( ncid_recout, nc_ph_varid, values = rec_ph) )
        call check( nf90_put_var( ncid_recout, nc_thr_varid, values = rec_th_req) )
        call check( nf90_put_var( ncid_recout, nc_proc_varid, values = rec_proc) ) 

        do irec=1,nrec
            call check( nf90_put_var( ncid_recout, nc_recnam_varid, start = (/irec, 1/), &
                                      count = (/1, 40/), values = (rec_names(irec))) )
        end do

        ! Write out seismogram dump times
        time_seis = dble((/ (i, i = 1, nseismo) /)) * deltat
        call check( nf90_put_var( ncid_recout, nc_time_varid, values = time_seis ) ) 
        if (verbose > 1) write(6,*) '...done'

        ! Write out STFs
        call check( nf90_put_var(ncid = ncid_recout, varid = nc_stf_iter_varid, &
                                 values = stf_dumpvar) )
        call check( nf90_put_var(ncid = ncid_recout, varid = nc_stf_seis_varid, &
                                 values = stf_seis_dumpvar) )

        ! Write out strain dump times
        if (dump_wavefields) then
            time_strain = dble((/ (i, i = 1, nstrain) /))
            time_strain = time_strain * t_0 / strain_samp
            call check( nf90_put_var(ncid   = ncid_out, &
                                     varid  = nc_snaptime_varid, &
                                     values = time_strain ) ) 
            ! Write out discontinuity depths
            call check( nf90_put_var(ncid   = ncid_snapout, &
                                     varid  = nc_disc_varid, &
                                     values = discont) )
                                     
        end if
    
    end if ! (mynum == 0)
   

! Allocation of Dump buffer variables. Done on all procs
90  format(' Allocated ', A20, ', uses ',F10.3,'MB')
    allocate(recdumpvar(nseismo,3,num_rec))
    recdumpvar = 0.0

    if (dump_wavefields) then
        allocate(surfdumpvar_disp(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_velo(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_strain(dumpstepsnap,6,maxind))
        allocate(surfdumpvar_srcdisp(dumpstepsnap,3,maxind))
       
        if (src_type(1) == 'monopole') then
            !allocate(oneddumpvar_flu(nel_fluid*gllperelem, dumpstepsnap, 6))
            !allocate(oneddumpvar_sol(nel_solid*gllperelem, dumpstepsnap, 6))
            allocate(oneddumpvar(npoints, dumpstepsnap, 6))
        else
            !allocate(oneddumpvar_flu(nel_fluid*gllperelem, dumpstepsnap, 9))
            !allocate(oneddumpvar_sol(nel_solid*gllperelem, dumpstepsnap, 9))
            allocate(oneddumpvar(npoints, dumpstepsnap, 9))
        end if

        !oneddumpvar_flu = 0.0
        !oneddumpvar_sol = 0.0 
        surfdumpvar_disp = 0.0
        surfdumpvar_velo = 0.0
        surfdumpvar_strain = 0.0
        surfdumpvar_srcdisp = 0.0

        if (mynum == 0 .and. verbose > 1) then
            write(6,*)  'Allocating NetCDF buffer variables'
            write(6,90) 'recdumpvar', real(size(recdumpvar))/262144.
            write(6,90) 'surfdumpvar_disp', real(size(surfdumpvar_disp))/262144.
            write(6,90) 'surfdumpvar_velo', real(size(surfdumpvar_velo))/262144.
            write(6,90) 'surfdumpvar_strain', real(size(surfdumpvar_strain))/262144.
            write(6,90) 'surfdumpvar_srcdisp', real(size(surfdumpvar_srcdisp))/262144.
            write(6,90) 'oneddumpvar', real(size(oneddumpvar))/262144.
            !write(6,90) 'oneddumpvar_flu', real(size(oneddumpvar_flu))/262144.
            !write(6,90) 'oneddumpvar_sol', real(size(oneddumpvar_sol))/262144.
            write(6,*)
            write(6,*)  '*********************************************************************'
            write(6,*)  '**** NetCDF output file produced, buffer variables all allocated ****'
            write(6,*)  '*********************************************************************'
            write(6,*)
        end if
    end if
    
#endif
end subroutine nc_define_outputfile
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Character
subroutine nc_write_att_char(attribute_value, attribute_name)
    character(len=*), intent(in)    :: attribute_name, attribute_value

#ifdef unc
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_char
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Real
subroutine nc_write_att_real(attribute_value, attribute_name)
  character(len=*),  intent(in)      :: attribute_name
  real(sp), intent(in)               :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_real
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Integer
subroutine nc_write_att_int(attribute_value, attribute_name)
  character(len=*),  intent(in)     :: attribute_name
  integer, intent(in)               :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! MvD: we are not really opening in parrallel any more, are we? Then this is a
!      missleading name + comment

!> Open the NetCDF output file, check for variable IDs and dump meshes.
subroutine nc_finish_prepare
#ifdef unc
    use data_io,   ONLY  : datapath, lfdata, dump_wavefields
    use data_mesh, ONLY  : maxind, surfcoord
    integer             :: status, ivar, nmode, iproc
    integer             :: nc_mesh_s_varid, nc_mesh_z_varid
    
    if (mynum == 0) then
        call check(nf90_close(ncid_out))
        if (verbose > 1) then
           write(6,*) '  Root process closed netCDF file, waiting for all procs to'
           write(6,*) '  arrive here and then open it to retrieve IDs'
           if (dump_wavefields) write(6,*) '  and dump mesh coordinates.'
        endif
    end if
    call barrier

    do iproc = 0, nproc
        call barrier
        if (iproc.eq.mynum) then
            nmode = ior(NF90_WRITE, NF90_NETCDF4)
            call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                                  mode=nmode, ncid=ncid_out) )

            call check( nf90_inq_grp_ncid(ncid_out, "Seismograms", ncid_recout) )
            call check( nf90_inq_grp_ncid(ncid_out, "Surface", ncid_surfout) )
            call check( nf90_inq_grp_ncid(ncid_out, "Mesh", ncid_meshout) )

            call check( nf90_inq_varid( ncid_recout, "displacement", nc_disp_varid ) )
            
            if (dump_wavefields) then
                call check( nf90_inq_grp_ncid(ncid_out, "Snapshots", ncid_snapout) )
                do ivar=1, nvar/2
                    call check( nf90_inq_varid( ncid_snapout, nc_varnamelist(ivar), &
                                nc_field_varid(ivar)) )
                end do
                
                call check( nf90_inq_varid(ncid_surfout, "elem_theta", &
                                           nc_surfelem_theta_varid) )
                call check( nf90_put_var(ncid_surfout, &
                                         varid  = nc_surfelem_theta_varid, &
                                         values = surfcoord , &
                                         start  = (/mynum*maxind+1/), &
                                         count  = (/maxind/) ))
                
                call check( nf90_inq_varid(ncid_surfout, "displacement", &
                                           nc_surfelem_disp_varid) )

                call check( nf90_inq_varid( ncid_surfout, "velocity", &
                                        nc_surfelem_velo_varid ) )

                call check( nf90_inq_varid( ncid_surfout, "strain", &
                                        nc_surfelem_strain_varid ) )

                call check( nf90_inq_varid( ncid_surfout, "disp_src", &
                                        nc_surfelem_disp_src_varid ) )
            
                call check( nf90_inq_varid( ncid_meshout, "mesh_S", nc_mesh_s_varid ) )
                call check( nf90_put_var(ncid_meshout, varid = nc_mesh_s_varid, &
                                         values = scoord1d, &
                                         start  = (/mynum*npoints+1/), &
                                         count  = (/npoints/) ))
                call check( nf90_inq_varid( ncid_meshout, "mesh_Z", nc_mesh_z_varid ) )
                call check( nf90_put_var(ncid_meshout, varid = nc_mesh_z_varid, &
                                         values = zcoord1d, &
                                         start  = (/mynum*npoints+1/), &
                                         count  = (/npoints/) ))

            end if !dump_wavefields
            call check( nf90_close( ncid_out))
            if (verbose > 1) write(6,"('Proc ', I3, ' dumped its mesh and is ready to rupture')") mynum
        end if !mynum.eq.iproc
    end do
#endif
end subroutine nc_finish_prepare
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Final dumps to netCDF file. In the moment contains only dump of 
!! receiver seismograms.
subroutine nc_end_output
#ifdef unc
    use data_io, only: dump_xdmf
    integer         :: iproc

    call barrier
    do iproc=0, nproc-1
        if (iproc == mynum) then
            if (verbose > 1) write(6,"('Proc ', I3, ' will dump receiver seismograms')") mynum
            call nc_dump_rec_to_disk()
            if (dump_xdmf) then
                call check(nf90_close(ncid_out_snap))
            end if
        end if
        
        call barrier
    end do
    !call check( nf90_close(ncid_out) )

#endif
end subroutine nc_end_output
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_make_snapfile

    use data_mesh,    only: npoint_plot, nelem_plot
    use data_proc,    only: appmynum
    use data_io,      only: datapath, lfdata, nsnap
    use data_source,  only: src_type

#ifdef unc
    integer              :: nmode, nc_snappoint_dimid, nc_snapelem_dimid, nc_snapdim_dimid
    integer              :: nc_snaptime_dimid, nc_snapconnec_dimid
    character(len=120)   :: fname

    if (lpr .and. verbose > 1) write(6,*) '   .... preparing xdmf nc file'

    if (src_type(1) == 'monopole') then
        ndim_disp = 2
    else
        ndim_disp = 3
    end if

    fname = datapath(1:lfdata) // '/netcdf_snap_' // appmynum // '.nc'
    nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
    call check(nf90_create(path=fname, cmode=nmode, ncid=ncid_out_snap) )

    call check(nf90_def_dim(ncid_out_snap, 'points', npoint_plot, nc_snappoint_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'elements', nelem_plot, nc_snapelem_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'dimensions', ndim_disp , nc_snapdim_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 's-z-coordinate', 2 , nc_coord_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'connections', 4 , nc_snapconnec_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'timesteps', nsnap , nc_snaptime_dimid) )

    call flush(6)
    call check(nf90_def_var(ncid   = ncid_out_snap, & 
                            name   = 'displacement',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_snapdim_dimid, nc_snappoint_dimid, &
                                      nc_snaptime_dimid], & 
                            chunksizes = [ndim_disp, npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_disp_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, & 
                            name   = 'points',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_coord_dimid, nc_snappoint_dimid], & 
                            varid  = nc_snap_point_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, & 
                            name   = 'grid',  &
                            xtype  = NF90_INT,     &
                            dimids = [nc_snapconnec_dimid, nc_snapelem_dimid], & 
                            varid  = nc_snap_grid_varid) )

    call check(nf90_enddef(ncid    = ncid_out_snap))

    if (lpr .and. verbose > 1) write(6,*) '   .... DONE'

#endif

end subroutine nc_make_snapfile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snapshot(u)

    use data_mesh,      only: npoint_plot
    use data_io,        only: isnap, nsnap
    use data_source,    only: src_type

    real(kind=realkind), dimension(3,npoint_plot), intent(in)       :: u

#ifdef unc
    if (src_type(1) == 'monopole') then
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [1, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(1,:)) )
       
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [2, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(3,:)) )
    else
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [1, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(1,:)) )
       
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [2, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(2,:)) )

       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [3, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(3,:)) )
    end if

#endif

end subroutine nc_dump_snapshot
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snap_points(points)

    use data_mesh, only: npoint_plot
    real(sp), dimension(2,npoint_plot), intent(in)       :: points

#ifdef unc
    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_point_varid, &
                            count  = [2, npoint_plot], &
                            values = points) )
#endif

end subroutine nc_dump_snap_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snap_grid(grid)

    use data_mesh, only: nelem_plot
    integer, dimension(4, nelem_plot), intent(in)       :: grid 

#ifdef unc
    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_grid_varid, &
                            count  = [4, nelem_plot], &
                            values = grid) )
#endif
end subroutine nc_dump_snap_grid
!-----------------------------------------------------------------------------------------


end module nc_routines
