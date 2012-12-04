!> Contains all the routines for NetCDF handling.
module nc_routines

#ifdef unc
  use netcdf
#endif
  use data_proc, ONLY : mynum, nproc, lpr
  use global_parameters
  real, allocatable   :: recdumpvar(:,:,:)       !< Buffer variable for recorder 
  real, allocatable   :: surfdumpvar_disp(:,:,:)  !< Buffer variable for displacement at surface
  real, allocatable   :: surfdumpvar_velo(:,:,:)  !< Buffer variable for velocity at surface 
  real, allocatable   :: surfdumpvar_strain(:,:,:) !< Buffer variable for strain at surface
  real, allocatable   :: surfdumpvar_srcdisp(:,:,:) !< Buffer variable for source displacement at surface
  real, allocatable   :: oneddumpvar_flu(:,:,:)  !< Buffer variable for everything fluid dumped in nc_dump_field_1d
  real, allocatable   :: oneddumpvar_sol(:,:,:)  !< Buffer variable for everything solid dumped in nc_dump_field_1d
  integer, parameter  :: dumpsteprec = 100       !< Number of steps before receiver specific stuff is dumped
  integer             :: dumpstepsnap            !< Number of steps before kernel specific stuff is dumped 
  integer             :: gllperelem              !< Number of GLL points per element 
  integer             :: outputplan              !< When is this processor supposed to dump. 
  integer             :: stepstodump             !< How many steps since last dump?
  integer,allocatable :: varsfilled(:)           !< Has variable been written to buffer in this step?
  integer,allocatable :: varlength(:)            !< Length of strain dump variables
  integer             :: isnap_global            !< Global variables, so that we do not have to pass data to the C subroutine
  integer             :: ndumps                  !< dito
  logical,allocatable :: dumpposition(:)         !< Will any processor dump at this value of isnap?
  
  ! Stuff moved from data_io
  integer             :: ncid_out, ncid_recout, ncid_snapout, ncid_surfout
  integer             :: nc_snap_dimid, nc_proc_dimid, nc_rec_dimid, nc_recproc_dimid
  integer             :: nc_times_dimid, nc_comp_dimid, nc_disp_varid
  integer             :: nc_strcomp_dimid
  integer             :: nc_surfelem_disp_varid, nc_surfelem_velo_varid
  integer             :: nc_surfelem_strain_varid, nc_surfelem_disp_src_varid
  integer,allocatable :: nc_field_varid(:)
  character(len=16), allocatable  :: varnamelist(:)
  integer             :: nvar = -1



!! @todo These parameters should move to a input file soon
  integer             :: dumpbuffersize = 256    !< How often should each processor dump its buffer to disk?
  logical             :: deflate        = .false. !< Should output be compressed?
  integer             :: deflate_level  = 5      !< Compression level (0 lowest, 9 highest)
  logical             :: verbose        = .true. !< @todo Will hopefully be a global variable one day
contains
!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    implicit none
    integer, intent ( in) :: status !< Error code
#ifdef unc
    if(status /= nf90_noerr) then 
        !if(status.eq.-101) then
         !   write(6,*) 'Nonfatal HDF5 problem' 
        !else
            print *, trim(nf90_strerror(status))
            stop 2
        !end if
    end if
#endif
end subroutine check  
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
!! @todo: Change from local (processor-specific) IO to global.
subroutine nc_dump_field_1d(f, flen, varname, appisnap)
    use data_io,     ONLY : nstrain
    use commun,      ONLY : barrier
    use data_source, ONLY : src_type

    implicit none
    include 'mesh_params.h'
    integer, intent(in)               :: flen     !< Length of f
    real(kind=realkind), intent(in)   :: f(flen)     !< Data to dump
    character(len=*), intent(in)      :: varname  !< Internal name of data to dump. 
    !! Is used to identify the NetCDF Variable in question
    character(len=4), intent(in)      :: appisnap !< String which contains snapshot number
#ifdef unc
    integer                           :: isnap_loc, ivar, iproc
    real                              :: tick, tack

    read(appisnap,*) isnap_loc

    do ivar = 1, nvar
        !< Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar))==trim(varname)) exit
    end do
    if (ivar>nvar) then
        write(6,*) 'nc_dump_field_1d: Trying to access variable: ', trim(varname), &
            ' which is not in varnamelist. Contact a developer and shout at him!'
        stop 1
    end if
20  format('Something is so wrong here', A16, ', flen', I8, ' nel_solid ', I8, /&
           'Probably an error during recompilation. Please run "make clean" and', /&
           'try again. If problem persists, please send a mail to the developers')
21  format('Something is so wrong here', A16, ', flen', I8, ' nel_fluid ', I8, /&
           'Probably an error during recompilation. Please run "make clean" and', /&
           'try again. If problem persists, please send a mail to the developers')

    varsfilled(ivar) = 1
    if (src_type(1).eq.'monopole') then
        if (ivar<=4) then !solid variable
            if (flen .ne. nel_solid*gllperelem) then
                write(6,20)  trim(varname), flen, nel_solid*gllperelem
                stop 2
            end if
            oneddumpvar_sol(1:flen,stepstodump+1,ivar) = f  !processor specific dump variable
        elseif (ivar==5) then ! solid displacement
            if (flen .eq. nel_solid*2*gllperelem) then ! Monopole Source
                oneddumpvar_sol(1:nel_solid, stepstodump+1, 5:6) = &
                                reshape(f,(/nel_solid,2/)) 
            else
                write(6,20)  trim(varname), flen, nel_solid*gllperelem
                stop 2
            end if

        elseif ((ivar>5).and.(ivar<10)) then !fluid variable
            if (flen .ne. nel_fluid*gllperelem) then
                write(6,21)  trim(varname), flen, nel_fluid*gllperelem
                stop 2
            end if
            oneddumpvar_flu(1:flen,stepstodump+1,ivar-5) = f  
        elseif (ivar==10) then ! fluid displacement
            if (flen .ne. nel_fluid*3*gllperelem) then
                write(6,21)  trim(varname), flen, nel_fluid*gllperelem
                stop 2
            end if
            oneddumpvar_flu(1:nel_fluid, stepstodump+1,5:7) = &
                            reshape(f,(/nel_fluid,3/))
        end if

    else  ! Dipole or quadpole source
        if (ivar<=6) then !solid variable
            if (flen .ne. nel_solid*gllperelem) then
                write(6,20)  trim(varname), flen, nel_solid*gllperelem
                stop 2
            end if
            oneddumpvar_sol(1:flen,stepstodump+1,ivar) = f  !processor specific dump variable
        elseif (ivar==7) then ! solid displacement
            if (flen .ne. nel_solid*3*gllperelem) then 
                write(6,20)  trim(varname), flen, nel_solid*gllperelem
                stop 2
            end if
            oneddumpvar_sol(1:nel_solid, stepstodump+1, 7:9) = &
                            reshape(f,(/nel_solid,3/)) 

        elseif ((ivar>7).and.(ivar<14)) then !fluid variable
            if (flen .ne. nel_fluid*gllperelem) then
                write(6,21)  trim(varname), flen, nel_fluid*gllperelem
                stop 2
            end if
            oneddumpvar_flu(1:flen,stepstodump+1,ivar-7) = f  
        elseif (ivar==14) then ! fluid displacement
            if (flen .ne. nel_fluid*3*gllperelem) then
                write(6,21)  trim(varname), flen, nel_fluid*gllperelem
                stop 2
            end if
            oneddumpvar_flu(1:nel_fluid, stepstodump+1,7:9) = &
                            reshape(f,(/nel_fluid,3/))
        end if
    end if
    
#endif
end subroutine nc_dump_field_1d
!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain(isnap_loc)

    use data_io, ONLY    : nstrain
    use commun,  ONLY    : barrier
    implicit none
    integer, intent(in) :: isnap_loc
#ifdef unc
    integer             :: iproc
    real                :: tick, tack
    
    stepstodump = stepstodump + 1
    if (isnap_loc.eq.0) return


    if (dumpposition(mod(isnap_loc, dumpstepsnap))) then

        if (mynum.eq.0) call cpu_time(tick)
        do iproc=0,nproc-1
            if (iproc.eq.mynum) then
                call c_wait_for_io()
                call flush(6)
            end if
            call barrier
        end do
        if (mynum.eq.0) then
            call cpu_time(tack)
            if ((tack-tick).gt.0.5) then
30              format('Computiation was halted for ', F7.2, ' s to wait for ',/&
                       'dumping processor. Consider adapting netCDF output variables',/&
                       '(disable compression, increase dumpstepsnap)')
                write(6,30) tack-tick
            end if
        end if

        if (mod(isnap_loc, dumpstepsnap).eq.outputplan) then 
            call c_wait_for_io()
            isnap_global = isnap_loc 
            ndumps = stepstodump
            call c_spawn_dumpthread(stepstodump)

            stepstodump=0
        end if

    end if 

    
                                        

    if (isnap_loc.eq.nstrain) then !< Final and last dump of all remaining 
        do iproc=0,nproc-1 !< Make sure, nobody is accessing the output file anymore
            if (iproc.eq.mynum) then
                call c_wait_for_io()
                call barrier
            end if
        end do
        !call c_wait_for_io()
        do iproc=0,nproc-1
            if (iproc.eq.mynum) then
                isnap_global = nstrain 
                ndumps = stepstodump
                !call c_spawn_dumpthread(stepstodump)
                !call flush(6)
                !call c_wait_for_io()
                call nc_dump_strain_to_disk()
                write(6,*) mynum, 'finished dumping strain'
            end if
            call barrier
        end do
    end if

#endif
end subroutine nc_dump_strain

!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain_to_disk()
#ifdef unc

    use data_io
    use data_source,       ONLY: src_type
    use global_parameters, ONLY: realkind
    use data_mesh,         ONLY: loc2globrec, maxind

    implicit none
    include 'mesh_params.h'
    !integer, intent(in)               :: stepstodump
    integer                           :: ivar, flen, isnap_loc
    real                              :: tick, tack
    integer                           :: dumpsize
    integer                           :: nvar_deriv
        
    dumpsize = 0
    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                          mode=NF90_WRITE, ncid=ncid_out) )
    isnap_loc = isnap_global
    if (ndumps.eq.0) then 
        write(6,41) mynum, isnap_loc
        return
    else
        write(6,40) mynum, isnap_loc, ndumps
    end if
40  format('  Proc ', I4, ' in dump routine, isnap =', I5, ', stepstodump = ', I4)
41  format('  Proc ', I4, ' in dump routine, isnap =', I5, ', nothing to dump, returning...')
   

    !do ivar = 1, nvar
    if (src_type(1).eq.'monopole') then
        nvar_deriv = 4
    else
        nvar_deriv = 6
    end if

    do ivar=1,nvar_deriv
        flen = varlength(ivar)
        call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar),  &
                                 start=(/1, mynum+1, isnap_loc-ndumps+1/),              &
                                 count=(/flen, 1, ndumps/),                             &
                                 values=oneddumpvar_sol(1:flen,1:ndumps,ivar) ) )
        dumpsize = dumpsize + flen*ndumps
    end do
        
    ivar = nvar_deriv + 1
    flen = varlength(ivar)
    call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar),  &
                             start=(/1, mynum+1, isnap_loc-ndumps+1/), &
                             count=(/flen, 1, ndumps/),                &
                             values=reshape(source=oneddumpvar_sol(1:flen/3,1:ndumps,ivar:ivar+2), &
                                            shape=(/flen/3, 3, ndumps/), &
                                            order=(/1,3,2/) ) ) )
    dumpsize = dumpsize + flen*ndumps
    
    do ivar=nvar_deriv + 2,nvar_deriv*2+1 !< fluid variable
        flen = varlength(ivar)
        call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar),  &
                                 start=(/1, mynum+1, isnap_loc-ndumps+1/),              &
                                 count=(/flen, 1, ndumps/),                             &
                                 values=oneddumpvar_flu(1:flen,1:ndumps,ivar-nvar_deriv-1) ) )
        dumpsize = dumpsize + flen*ndumps
    end do
        
    ivar = nvar_deriv*2+2
    flen = varlength(ivar)
    call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar),  &
                             start=(/1, mynum+1, isnap_loc-ndumps+1/), &
                             count=(/flen, 1, ndumps/),                &
                             values=reshape(source=oneddumpvar_flu(1:flen/3,1:ndumps,&
                                                   nvar_deriv+1:nvar_deriv+3),&
                                            shape=(/flen/3, 3, ndumps/), &
                                            order=(/1,3,2/) ) ) )
    dumpsize = dumpsize + flen*ndumps

            
    !> Surface dumps 
    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_varid, &
                start = (/isnap_loc-ndumps+1, 1, (nproc-1)*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_disp(:,:,:)) )
    dumpsize = dumpsize + 3*maxind*ndumps
    !print *, 'Wrote surfdumpvar_disp'
    call check( nf90_put_var(ncid_surfout, nc_surfelem_velo_varid, &
                start = (/isnap_loc-ndumps+1, 1, (nproc-1)*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_velo(:,:,:)) )
    dumpsize = dumpsize + 3*maxind*ndumps
    !print *, 'Wrote surfdumpvar_velo'
    call check( nf90_put_var(ncid_surfout, nc_surfelem_strain_varid, &
                start = (/isnap_loc-ndumps+1, 1, (nproc-1)*maxind+1/), &
                count = (/ndumps, 6, maxind/), &
                values = surfdumpvar_strain(:,:,:)) )
    dumpsize = dumpsize + 6*maxind*ndumps
    !print *, 'Wrote surfdumpvar_strain'
    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_src_varid, &
                start = (/isnap_loc-ndumps+1, 1, (nproc-1)*maxind+1/), &
                count = (/ndumps, 3, maxind/), &
                values = surfdumpvar_srcdisp(:,:,:)) )
    dumpsize = dumpsize + 3*maxind*ndumps
    !print *, 'Wrote surfdumpvar_srcdisp'


    call cpu_time(tack)
    call check( nf90_close(ncid_out) ) 
    write(6,70) mynum, real(dumpsize) * realkind / 1048576., tack-tick 
70  format('  Proc', I5,': Wrote ', F8.2, ' MB in ', F6.2, 's')    
    call flush(6)
    oneddumpvar_flu = 0.0
    oneddumpvar_sol = 0.0 
    surfdumpvar_disp = 0.0
    surfdumpvar_velo = 0.0
    surfdumpvar_strain = 0.0
    surfdumpvar_srcdisp = 0.0

#endif
end subroutine nc_dump_strain_to_disk

!-----------------------------------------------------------------------------------------
!> Dump receiver specific stuff, especially displacement and velocity
!! N.B.: Works with global indices.
subroutine nc_dump_rec(recfield)
    use data_mesh, ONLY: num_rec
    use data_io,   ONLY: iseismo
    implicit none
    !integer, intent(in)                    :: iseismo 
    real, intent(in), dimension(3,num_rec) :: recfield
#ifdef unc
    
    recdumpvar(iseismo,:,:) = recfield(:,:)

#endif
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_rec_to_disk
#ifdef unc
    use data_mesh, ONLY: loc2globrec, num_rec
    use data_io,   ONLY: datapath, lfdata, nseismo
    implicit none
    real                              :: tick, tack
    integer                           :: irec, status, dumpsize, icomp

    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                          mode=NF90_WRITE, ncid=ncid_out) )
    call check( nf90_inq_varid( ncid_recout, "displacement", nc_disp_varid ) )
    dumpsize = 0
    do irec = 1, num_rec
        do icomp = 1, 3
        !status =( nf90_put_var(ncid=ncid_out, varid=nc_disp_varid, &
            call check( nf90_put_var(ncid=ncid_recout, varid=nc_disp_varid, &
                                   start=(/1, icomp, loc2globrec(irec)/), &
                                   count = (/nseismo, 1, 1/), values=recdumpvar(:,icomp,irec)) )
        end do
        !if(status.ne.0) write(6,*) 'Nonfatal HDF5 problem in nc_dump_rec, idump:'
    end do
    dumpsize = nseismo * 3 * num_rec
    call check( nf90_close(ncid=ncid_out))
    call cpu_time(tack)
    write(6,70) mynum, real(dumpsize) * realkind / 1048576., tack-tick 
70  format(I3,': Receiver data, Wrote ', F8.3, ' MB in ', F6.2, 's')    
#endif
end subroutine nc_dump_rec_to_disk
!-----------------------------------------------------------------------------------------
!


!-----------------------------------------------------------------------------------------
!> Dump stuff along surface
subroutine nc_dump_surface(surffield, disporvelo, nrec, dim2, idump)
    use data_mesh, ONLY: maxind
    implicit none
    integer, intent(in)                          :: nrec, dim2, idump
    real(kind=realkind), intent(in), dimension(nrec,dim2) :: surffield
    character(len=4), intent(in)                 :: disporvelo
#ifdef unc
    integer                                      :: irec, nc_varid
   
    !write(6,*) disporvelo, stepstodump, nrec, dim2, (nproc-1)*maxind, maxind
10  format('In nc_dump_surface to dump ', A, ' no.', I4, ' nrec=', I6,  &
           ' dim2=', I1, ' start=', I6, ' maxind=', I6)
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
!> Define the output file variables and dimensions
!! and allocate buffer variables.
subroutine nc_define_outputfile(nrec, rec_names, rec_th, rec_th_req, rec_ph, rec_proc)

    use data_io,     ONLY: nseismo, nstrain, nseismo, ibeg, iend, dump_wavefields
    use data_time,   ONLY: strain_it
    use data_mesh,   ONLY: maxind, num_rec
    use data_source, ONLY: src_type
    implicit none
    include 'mesh_params.h'

    integer, intent(in)                 :: nrec              !< Number of receivers
    character(len=40),intent(in)        :: rec_names(nrec)   !< Receiver names
    real(8), dimension(nrec),intent(in) :: rec_th            !< Receiver theta 
    real(8), dimension(nrec),intent(in) :: rec_th_req        !< Requested receiver theta
    real(8), dimension(nrec),intent(in) :: rec_ph            !< Receiver phi
    integer, dimension(nrec),intent(in) :: rec_proc          !< Receiver processor
#ifdef unc
    character(len=16), allocatable      :: varname(:)
    integer                             :: ivar
    integer, allocatable                :: nc_f_dimid(:)
    integer                             :: irec
    integer                             :: nc_latr_varid, nc_lon_varid 
    integer                             :: nc_lat_varid, nc_ph_varid
    integer                             :: nc_thr_varid, nc_th_varid 
    integer                             :: nc_proc_varid, nc_recnam_dimid
    integer                             :: nc_recnam_varid, nc_surf_dimid
   
    if (mynum.eq.0) then
        write(6,*) '**** Producing netcdf output file with its variables and dimensions ****'
    end if
    if (dump_wavefields) then
        if (src_type(1).eq.'monopole') then
            nvar = 10
        else
            nvar = 14
        end if
        allocate(varname(nvar))
        allocate(varnamelist(nvar))
        allocate(varlength(nvar))
        allocate(varsfilled(nvar))
        allocate(nc_field_varid(nvar))
        allocate(nc_f_dimid(nvar))

        if (src_type(1) .eq. 'monopole') then 
            varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                            'straintrace_sol', 'velo_sol       ', 'strain_dsus_flu', &
                            'strain_dsuz_flu', 'strain_dpup_flu', 'straintrace_flu', &
                            'velo_flu       '/)
              
            varlength = (/nel_solid, nel_solid,   nel_solid, &
                          nel_solid, nel_solid*3, nel_fluid, &
                          nel_fluid, nel_fluid,   nel_fluid, &
                          nel_fluid*3/)
        else
            varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                            'strain_dsup_sol', 'strain_dzup_sol', 'straintrace_sol', &
                            'velo_sol       ', &
                            'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                            'strain_dsup_flu', 'strain_dzup_flu', 'straintrace_flu', &
                            'velo_flu       '/)
              
            varlength = (/nel_solid, nel_solid,   nel_solid, &
                          nel_solid, nel_solid,   nel_solid, &
                          nel_solid*3, &
                          nel_fluid, nel_fluid,   nel_fluid, &
                          nel_fluid, nel_fluid,   nel_fluid, &
                          nel_fluid*3/)
        end if

        gllperelem = (iend-ibeg+1)**2
        varlength = varlength * gllperelem

        if (mynum.eq.0) then
            !write(6,*) gllperelem, nel_fluid, nel_solid
            !do ivar = 1, nvar
            !    write(6,*) varnamelist(ivar), varlength(ivar)
            !end do
        end if
        if (nstrain<=dumpstepsnap) then
            dumpstepsnap = nstrain
        end if
    end if ! dump_wavefields
    
    if (mynum.eq.0) then    
        write(6,*) '  Producing groups for Seismograms and Snapshots'
        call check( nf90_def_grp( ncid_out, "Seismograms", ncid_recout) )
        call check( nf90_def_grp( ncid_out, "Snapshots", ncid_snapout) )
        call check( nf90_def_grp( ncid_out, "Surface", ncid_surfout) )
        write(6,*) '  Seismograms group has ID', ncid_recout, &
                   ', Snapshots group has ID', ncid_snapout, & 
                   ', Surface group has ID', ncid_surfout
        
        write(6,*) 'Define dimensions in ''Seismograms'' group of NetCDF output file'
        write(6,*) '  ''Seismograms'' group has ID ', ncid_recout
        call check( nf90_def_dim( ncid_out, "timesteps", nseismo, nc_times_dimid) )
        write(6,110) "timesteps", nseismo, nc_times_dimid 
        call check( nf90_def_dim( ncid_recout, "receivers", nrec, nc_rec_dimid) )
        write(6,110) "receivers", nrec, nc_rec_dimid
        !call check( nf90_def_dim( ncid_recout, "processors", nproc, nc_recproc_dimid) )
        !write(6,110) "processors", nproc, nc_recproc_dimid
        call check( nf90_def_dim( ncid_out, "components", 3, nc_comp_dimid) )
        write(6,110) "components", 3, nc_comp_dimid
        call check( nf90_def_dim( ncid_recout, "recnamlength", 40, nc_recnam_dimid) ) 
        write(6,110) "recnamlength", 40, nc_recnam_dimid
        write(6,*) 'NetCDF dimensions defined'
110     format(' Dimension ',A20,' with length ' I8,' and ID', I6) 

        write(6,*) 'Define variables in ''Seismograms'' group of NetCDF output file'
        call flush(6)
        call check( nf90_def_var( ncid=ncid_recout, name="displacement", xtype=NF90_FLOAT,&
                                  dimids=(/nc_times_dimid, nc_comp_dimid, nc_rec_dimid/), &
                                  !storage = NF90_CHUNKED, chunksizes=(/nseismo,3,1/), &
                                  !contiguous = .false., chunksizes=(/nseismo,3,1/), &
                                  varid=nc_disp_varid) )
        call check( nf90_put_att(ncid_recout, nc_disp_varid, 'units', 'meters') )
        call check( nf90_put_att(ncid_recout, nc_disp_varid, '_FillValue', 0.0) )

        call check( nf90_def_var( ncid_recout, "Lat_req", NF90_FLOAT, (/nc_rec_dimid/), &
                                  nc_latr_varid) )

        call check( nf90_def_var( ncid_recout, "Lon", NF90_FLOAT, (/nc_rec_dimid/), &
                                  nc_lon_varid) )
        call check( nf90_def_var( ncid_recout, "Lat", NF90_FLOAT, (/nc_rec_dimid/), &
                                  nc_lat_varid) )
        call check( nf90_def_var( ncid_recout, "azimuths", NF90_FLOAT, (/nc_rec_dimid/),  &
                                  nc_ph_varid) )
        call check( nf90_def_var( ncid_recout, "distances_requested", NF90_FLOAT, &
                                  (/nc_rec_dimid/), nc_thr_varid) )
        call check( nf90_def_var( ncid_recout, "distances", NF90_FLOAT, &
                                  (/nc_rec_dimid/), nc_th_varid) )
        call check( nf90_def_var( ncid_recout, "processor_of_receiver", NF90_INT, &
                                  (/nc_rec_dimid/), nc_proc_varid) )
        call check( nf90_def_var( ncid_recout, "receiver_name", NF90_CHAR, &
                                  (/nc_rec_dimid, nc_recnam_dimid/), nc_recnam_varid) )

    
        
        if (dump_wavefields) then
            ! Wavefields group of output file N.B: Snapshots for kernel calculation
            write(6,*) 'Define variables in ''Snapshots'' group of NetCDF output file'
            write(6,*) '  awaiting', nstrain, ' snapshots'
            call check( nf90_def_dim( ncid=ncid_snapout, name="snapshots", len=nstrain, &
                                      dimid=nc_snap_dimid) )
            call check( nf90_def_dim( ncid=ncid_snapout, name="processors", len=nproc, &
                                      dimid=nc_proc_dimid) )
            do ivar=1, nvar ! The big snapshot variables for the kerner.
                call check( nf90_def_dim(ncid=ncid_snapout, &
                                         name="dim_"//trim(varnamelist(ivar)), &
                                         len=varlength(ivar), dimid=nc_f_dimid(ivar)) )
                call check( nf90_def_var(ncid=ncid_snapout, name=trim(varnamelist(ivar)), &
                                         xtype = NF90_FLOAT, &
                                         dimids = (/nc_f_dimid(ivar), nc_proc_dimid, nc_snap_dimid/),&
                                         varid = nc_field_varid(ivar), &
                                         chunksizes = (/varlength(ivar), 1, 1/) ))
       
                if (deflate) then
                    call check( nf90_def_var_deflate(ncid=ncid_snapout, &
                                                     varid=nc_field_varid(ivar), &
                                                     shuffle=1, deflate=1, &
                                                     deflate_level=deflate_level) )
                end if
                call check( nf90_def_var_fill(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
                                              no_fill=1, fill=0) )
                write(6,100) trim(varnamelist(ivar)), nc_field_varid(ivar), varlength(ivar)
100             format('  Netcdf variable ',A16,' with ID ', I3, ' and length', I8, ' produced.')       
            end do


            varsfilled = 0
        
            ! Surface group in output file
            write(6,*) 'Define variables in ''Surface'' group of NetCDF output file'
            call check( nf90_def_dim( ncid_surfout, name="snapshots", len=nstrain, &
                                      dimid=nc_snap_dimid) )
            call check( nf90_def_dim( ncid_surfout, "straincomponents", len=6, &
                                      dimid=nc_strcomp_dimid) )
            write(6,110) "straincomponents", 6, nc_strcomp_dimid
            
            call check( nf90_def_dim( ncid_surfout, "surf_elems", maxind*nproc, nc_surf_dimid) )     
            write(6,110) "surf_elems", maxind*nproc, nc_surf_dimid

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
        
        
        write(6,*) 'NetCDF variables defined'
        call check( nf90_enddef ( ncid_out))

        write(6,*) 'Writing station info into NetCDF file'
        call check( nf90_put_var( ncid_recout, nc_th_varid, values = rec_th) )
        call check( nf90_put_var( ncid_recout, nc_ph_varid, values = rec_ph) )
        call check( nf90_put_var( ncid_recout, nc_thr_varid, values = rec_th_req) )
        call check( nf90_put_var( ncid_recout, nc_proc_varid, values = rec_proc) ) 

        do irec=1,nrec
            !write(6,*) 'Write name ', rec_names(irec), ' of receiver',irec
            call check( nf90_put_var( ncid_recout, nc_recnam_varid, start = (/irec, 1/), &
                                      count = (/1, 40/), values = (rec_names(irec))) )
        end do
        write(6,*) 'done'
        call check( nf90_redef ( ncid_out))
    end if       
    
90  format(' Allocated ', A20, ', uses ',F10.3,'MB')
    allocate(recdumpvar(nseismo,3,num_rec))
    recdumpvar = 0.0

    if (dump_wavefields) then
        allocate(surfdumpvar_disp(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_velo(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_strain(dumpstepsnap,6,maxind))
        allocate(surfdumpvar_srcdisp(dumpstepsnap,3,maxind))
       
        if (src_type(1).eq.'monopole') then
            allocate(oneddumpvar_flu(nel_fluid*gllperelem, dumpstepsnap, 7) )
            allocate(oneddumpvar_sol(nel_solid*gllperelem, dumpstepsnap, 7) )
        else
            allocate(oneddumpvar_flu(nel_fluid*gllperelem, dumpstepsnap, 9) )
            allocate(oneddumpvar_sol(nel_solid*gllperelem, dumpstepsnap, 9) )
        end if

        oneddumpvar_flu = 0.0
        oneddumpvar_sol = 0.0 
        surfdumpvar_disp = 0.0
        surfdumpvar_velo = 0.0
        surfdumpvar_strain = 0.0
        surfdumpvar_srcdisp = 0.0
        if (mynum.eq.0) then
            write(6,*)  'Allocating NetCDF buffer variables'
            write(6,90) 'recdumpvar', real(size(recdumpvar))/262144.
            write(6,90) 'surfdumpvar_disp', real(size(surfdumpvar_disp))/262144.
            write(6,90) 'surfdumpvar_velo', real(size(surfdumpvar_velo))/262144.
            write(6,90) 'surfdumpvar_strain', real(size(surfdumpvar_strain))/262144.
            write(6,90) 'surfdumpvar_srcdisp', real(size(surfdumpvar_srcdisp))/262144.
            write(6,90) 'oneddumpvar_flu', real(size(oneddumpvar_flu))/262144.
            write(6,90) 'oneddumpvar_sol', real(size(oneddumpvar_sol))/262144.
            write(6,*)  'NetCDF output file produced, buffer variables all allocated'
        end if
    end if
    
#endif
end subroutine nc_define_outputfile
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Processor 0 opens output file so that nc_define_receiverfile can define 
!! the dimensions and variables in it.
subroutine define_netcdf_output
#ifdef unc

    use data_io
    use commpi,      ONLY : MPI_COMM_WORLD, MPI_INFO_NULL
    use commun,      ONLY : barrier
    implicit none    
    integer  :: nmode, ncid_test, test_dimid, test_varid, status
    integer  :: iproc, err

    nmode = ior(NF90_CLOBBER,NF90_NETCDF4)
    ! nmode = ior(nmode, nf90_mpiio)
    ! Testing, whether the netcdf4 library supports parallel IO. The only function
    ! which reliably crashes, when applied to the serial library is nf90_put_var, 
    ! so we try that here.
    call barrier
!    if (nproc>1) then
!        call check( nf90_create (path=datapath(1:lfdata)//"/parallel_test_.nc4", &
!                                 cmode=nmode, ncid=ncid_test, &
!                                 comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))
!        call check( nf90_def_dim(ncid=ncid_test, name='testdim', len=nproc, dimid=test_dimid))
!        call check( nf90_def_var(ncid=ncid_test, name='testvar', xtype=NF90_INT, &
!                                 dimids=test_dimid, varid=test_varid) )
!        call check( nf90_enddef(ncid_test) )
!        status = nf90_var_par_access(ncid_test,test_varid,nf90_collective)
!        if (status.ne.0) then
!            if (mynum.eq.0) then
!                write(6,*) '******************************************************************'
!                write(6,*) '* The NetCDF4 library on the system does not support Parallel    *'
!                write(6,*) '* HDF5 IO. This is a known issue in current Debian and probably  *'
!                write(6,*) '* Ubuntu releases. You probably have to rebuild the library with *'
!                write(6,*) '* support for parallel IO. See UTILS/NetCDF4_HOWTO.              *'
!                write(6,*) '******************************************************************'
!            end if
!            stop 1
!        end if
!        !call check( nf90_put_var(ncid_test, test_varid, values=mynum, start=mynum+1, count=1))
!        call check( nf90_close(ncid_test) )
!        call barrier
!        if (mynum.eq.0) write(6,*) 'Test for parallel NetCDF4 passed'
!    else 
!        if (mynum.eq.0) write(6,*) 'Test for parallel NetCDF4 not necessary.'
!    end if

    !dumpstepsnap = max(dumpbuffersize,nproc) 
    dumpstepsnap = int(dumpbuffersize/nproc+1)*nproc ! Will later be reduced to nstrain, if this is smaller
                                                  ! than value given here
    if (lpr) write(6,*) '  Dumping NetCDF file to disk every', &
                        dumpstepsnap, ' snaps'                                                  
    outputplan = mynum*(dumpstepsnap/nproc)

    allocate(dumpposition(0:dumpstepsnap-1))
    dumpposition = .false.
    do iproc = 0, nproc-1
        dumpposition(iproc*(dumpstepsnap/nproc)) = .true.
    end do

    write(6,60) mynum, outputplan
60  format('Proc ', I4, ' will dump at position ', I4)

    nmode = ior(NF90_CLOBBER,NF90_NETCDF4)
    if (mynum.eq.0) then
        write (6,*) 'Preparing netcdf file' ! for ', nproc, ' processors'
        nmode = ior(NF90_CLOBBER,NF90_NETCDF4)
        !call nc_set_log_level(3,err)
        call check( nf90_create ( path=datapath(1:lfdata)//"/axisem_output.nc4", &
                                  cmode=nmode, ncid=ncid_out))
        write(6,*) 'Netcdf file with ID ',ncid_out,' produced.'
    end if
#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Character
subroutine nc_write_att_char(attribute_value, attribute_name)
    use data_io
    character(len=*), intent(in)  :: attribute_name, attribute_value

#ifdef unc
    !write(6,*) 'Writing ', attribute_value, ' to attr. ', attribute_name, ' in netcdf ID', ncid_out
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )

#endif
end subroutine nc_write_att_char
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Real
subroutine nc_write_att_real(attribute_value, attribute_name)
  use data_io
  character(len=*),  intent(in)  :: attribute_name
  real, intent(in)                :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )

#endif
end subroutine nc_write_att_real
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Integer
subroutine nc_write_att_int(attribute_value, attribute_name)
  use data_io
  character(len=*),  intent(in)  :: attribute_name
  integer, intent(in)                :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )

#endif
end subroutine nc_write_att_int
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Open the NetCDF output file in parallel and check for variable IDs.
subroutine nc_open_parallel
#ifdef unc
    use data_io
    use commpi,     ONLY : MPI_COMM_WORLD, MPI_INFO_NULL
    use commun,     ONLY : barrier
    implicit none
    integer             :: status, ivar, nmode
    !integer, parameter :: par_access_mode = nf90_collective
    integer, parameter  :: par_access_mode = nf90_independent
    
!    varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
!                    'straintrace_sol', 'velo_sol       ', 'strain_dsus_flu', &
!                    'strain_dsuz_flu', 'strain_dpup_flu', 'straintrace_flu', &
!                    'velo_flu       '/)


    if (mynum==0) then
        call check(nf90_close(ncid_out))
        write(6,*) '  Root process closed netCDF file, waiting for all procs to'
        !write(6,*) '  arrive here and then open it for parallel IO'
        write(6,*) '  arrive here and then open it to retrieve IDs'
    end if
    call barrier
    !nmode = ior(NF90_WRITE,NF90_NETCDF4)
    nmode = ior(NF90_NOWRITE,NF90_NETCDF4)
    if (nproc>1) then 
        !nmode = IOR(nmode, nf90_mpiio)
        call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                              mode=nmode, ncid=ncid_out &
                              !,cache_nelems=nel_solid,& !size=500000000, &
                              !comm = MPI_COMM_WORLD, info = MPI_INFO_NULL) )
                              ) )
    else
        call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", & 
                              mode=nmode, ncid=ncid_out) )
    end if

    call check( nf90_inq_grp_ncid(ncid_out, "Seismograms", ncid_recout) )
    call check( nf90_inq_grp_ncid(ncid_out, "Surface", ncid_surfout) )
    if (dump_wavefields) then
        !call check( nf90_open(path=datapath(1:lfdata)//"/axisem_wavefield.nc4", & 
        !                      mode=NF90_NOWRITE, ncid=ncid_snapout) )
        call check( nf90_inq_grp_ncid(ncid_out, "Snapshots", ncid_snapout) )
        do ivar = 1,nvar
            call check( nf90_inq_varid( ncid_snapout, varnamelist(ivar), &
                        nc_field_varid(ivar)) )
        !    if (nproc>1) then
        !        call check( nf90_var_par_access(ncid=ncid_snapout,&
        !                                    varid=nc_field_varid(ivar),&
        !                                    access=par_access_mode))
        !    end if
        end do
        !call check( nf90_close( ncid_snapout))
        
        call check( nf90_inq_varid( ncid_surfout, "displacement", &
                                nc_surfelem_disp_varid ) )
!    if (nproc>1) call check( nf90_var_par_access(ncid=ncid_surfout,&
!                                    varid=nc_surfelem_disp_varid,&
!                                    access=par_access_mode))

        call check( nf90_inq_varid( ncid_surfout, "velocity", &
                                nc_surfelem_velo_varid ) )
!    if (nproc>1) call check( nf90_var_par_access(ncid=ncid_surfout,&
!                                    varid=nc_surfelem_velo_varid,&
!                                    access=par_access_mode))

        call check( nf90_inq_varid( ncid_surfout, "strain", &
                                nc_surfelem_strain_varid ) )
!    if (nproc>1) call check( nf90_var_par_access(ncid=ncid_surfout,&
!                                    varid=nc_surfelem_strain_varid,&
!                                    access=par_access_mode))

        call check( nf90_inq_varid( ncid_surfout, "disp_src", &
                                nc_surfelem_disp_src_varid ) )
!    if (nproc>1) call check( nf90_var_par_access(ncid=ncid_surfout,&
!                                    varid=nc_surfelem_disp_src_varid,&
!                                    access=par_access_mode))
    
    end if
    call check( nf90_inq_varid( ncid_recout, "displacement", nc_disp_varid ) )
!    if (nproc>1) call check( nf90_var_par_access(ncid=ncid_recout,&
!                                    varid=nc_disp_varid,&
!                                    access=par_access_mode))

    call check( nf90_close( ncid_out))
    write(6,70) mynum
70  format('Proc ', I3, ' opened file and is ready to rupture')    
#endif
end subroutine nc_open_parallel
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Close the Output file. Contains barrier.
subroutine end_netcdf_output
#ifdef unc
    !use data_io, only: ncid_out, ncid_recout
    use commun, only: barrier
    integer        :: iproc

    call barrier
    do iproc=0,nproc-1
        if (iproc.eq.mynum) then
            write(6,80) mynum
80          format('Proc ', I3, ' will dump receiver seismograms')    
            call nc_dump_rec_to_disk()
        end if
        call barrier
    end do
    !call check( nf90_close(ncid_out) )

#endif
end subroutine
!-----------------------------------------------------------------------------------------

!
!!> Remainder is just for the compile wo NetCDF case
!#else
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_dump_field_1d(f,flen,varname,appisnap)
!    implicit none
!    integer, intent(in)               :: flen
!    real(kind=realkind), intent(in)   :: f(flen)
!    character(len=*), intent(in)      :: varname
!    character(len=4), intent(in)      :: appisnap
!end subroutine nc_dump_field_1d
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_dump_rec(recfield, nc_varid, nrec, dim2, idump)
!    implicit none
!    integer, intent(in)                          :: nrec, dim2, idump, nc_varid
!    real(kind=realkind), intent(inout), dimension(nrec,dim2) :: recfield
!end subroutine
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_dump_surface(surffield, disporvelo, nrec, dim2, idump)
!    implicit none
!    integer, intent(in)                          :: nrec, dim2, idump
!    real(kind=realkind), intent(inout), dimension(nrec,dim2) :: surffield
!    character(len=4)                             :: disporvelo
!end subroutine nc_dump_surface
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_dump_rec_perproc(recfield, nc_varid, nrec, dim2, idump)
!    implicit none
!    integer, intent(in)                          :: nrec, dim2, idump, nc_varid
!    real(kind=realkind), intent(inout), dimension(nrec,dim2) :: recfield
!end subroutine nc_dump_rec_perproc
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_define_receiverfile(nrec, rec_names, rec_th, rec_th_req, rec_ph, rec_proc)
!    implicit none
!    integer, intent(in)                         :: nrec 
!    character(len=40),intent(in)                :: rec_names(nrec)
!    real(8), dimension(nrec),intent(in) :: rec_th, rec_th_req, rec_ph
!    integer, dimension(nrec),intent(in) :: rec_proc
!end subroutine nc_define_receiverfile
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine define_netcdf_output
!    implicit none    
!end subroutine define_netcdf_output
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_write_att_char(attribute_value, attribute_name)
!    character(len=*), intent(in)  :: attribute_name, attribute_value
!end subroutine nc_write_att_char
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_write_att_real(attribute_value, attribute_name)
!  character(len=*),  intent(in)  :: attribute_name
!  real, intent(in)                :: attribute_value
!end subroutine nc_write_att_real
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_write_att_int(attribute_value, attribute_name)
!  character(len=*),  intent(in)  :: attribute_name
!  integer, intent(in)                :: attribute_value
!end subroutine nc_write_att_int
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine nc_open_parallel
!    implicit none
!end subroutine nc_open_parallel
!!-----------------------------------------------------------------------------------------
!
!
!!-----------------------------------------------------------------------------------------
!subroutine end_netcdf_output
!
!end subroutine end_netcdf_output
!!-----------------------------------------------------------------------------------------
!#endif
end module nc_routines
