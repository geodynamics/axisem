program field_transformation

    use netcdf
    implicit none

    integer                         :: nsnap, nproc, nvar, nmode, ivar, nelemperstep
    integer                         :: ncin_id, ncin_snap_grpid, ncin_snap_dimid, ncin_proc_dimid 
    integer                         :: ncout_id, ncout_snap_grpid, ncout_snap_dimid, ncout_proc_dimid 
    integer                         :: nomega, nextpow2
    character(len=8)                :: sourcetype
    integer                         :: dimids(14,3)
    character(len=16), allocatable  :: varnamelist(:)
    integer, dimension(14)          :: varlength, ncin_field_varid, ncout_field_varid, ncout_field_dimid
    logical, parameter              :: deflate = .false.
    integer, parameter              :: deflate_level = 0

    real, dimension(:,:), allocatable :: strain_dsus_sol, strain_dsuz_sol, strain_dpup_sol, &
                                         straintrace_sol, velo_sol       , strain_dsus_flu, &
                                         strain_dsuz_flu, strain_dpup_flu, straintrace_flu, &
                                         velo_flu       

    complex, dimension(:,:), allocatable :: ft_strain_dsus_sol, ft_strain_dsuz_sol, ft_strain_dpup_sol, &
                                            ft_straintrace_sol, ft_velo_sol       , ft_strain_dsus_flu, &
                                            ft_strain_dsuz_flu, ft_strain_dpup_flu, ft_straintrace_flu, &
                                            ft_velo_flu       
    
    
    call check( nf90_open(path="bla/Data/axisem_output.nc4", & 
                          mode=NF90_NOWRITE, ncid=ncin_id) )

    call check( nf90_inq_grp_ncid(ncin_id, "Snapshots", ncin_snap_grpid) )
    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "excitation type", sourcetype))

    print *, sourcetype
    
    if (sourcetype=='monopole')  then
        nvar = 10
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                        'straintrace_sol', 'velo_sol       ', 'strain_dsus_flu', &
                        'strain_dsuz_flu', 'strain_dpup_flu', 'straintrace_flu', &
                        'velo_flu       '/)
    else
        nvar = 14
        allocate(varnamelist(nvar))
        varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                        'strain_dsup_sol', 'strain_dzup_sol', 'straintrace_sol', &
                        'velo_sol       ', &
                        'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                        'strain_dsup_flu', 'strain_dzup_flu', 'straintrace_flu', &
                        'velo_flu       '/)

        print *, 'who cares about multipoles?'
        stop
    end if

    do ivar=1, nvar
        print *, varnamelist(ivar)
        call check( nf90_inq_varid(ncin_snap_grpid, varnamelist(ivar), &
                                   varid=ncin_field_varid(ivar)) )
        call check( nf90_inquire_variable(ncin_snap_grpid, ncin_field_varid(ivar), dimids=dimids(ivar,:) ))
        call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(ivar,1), len=varlength(ivar)))
        print *, ncin_field_varid(ivar)
        print *, dimids(ivar,:)
        print *, varlength(ivar)
    end do

    call check( nf90_inq_dimid(ncin_snap_grpid, 'snapshots', dimid=ncin_snap_dimid) )
    call check( nf90_inquire_dimension(ncin_snap_grpid, ncin_snap_dimid, len=nsnap))
    print *, 'nsnap = ', nsnap
  
    nelemperstep = 10000
    
    !MvD: I guess varlength needs to be replaced by varlength * processors

    print *, 'allocating arrays for time domain data'
    allocate(strain_dsus_sol(varlength(1), nsnap))
    allocate(strain_dsuz_sol(varlength(2), nsnap))
    allocate(strain_dpup_sol(varlength(3), nsnap))
    allocate(straintrace_sol(varlength(4), nsnap))
    allocate(strain_dsus_flu(varlength(6), nsnap))
    allocate(strain_dsuz_flu(varlength(7), nsnap))
    allocate(strain_dpup_flu(varlength(8), nsnap))
    allocate(straintrace_flu(varlength(9), nsnap))
    allocate(velo_sol(varlength(5),   nsnap))
    allocate(velo_flu(varlength(10),  nsnap))
        

    !! Create output file
    print *, 'creating output file'
    nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
    call check( nf90_create(path="ordered_output.nc4", cmode=nmode, ncid=ncout_id))

    print *, '  Producing groups for Seismograms and Snapshots'
    call check( nf90_def_grp( ncout_id, "Snapshots", ncout_snap_grpid) )


    print *, 'Defining variables in ''Snapshots'' group of NetCDF output file'
    print *, '  awaiting', nsnap, ' snapshots'
    call check( nf90_def_dim(ncid=ncout_snap_grpid, name="snapshots", len=nsnap, &
                             dimid=ncout_snap_dimid) )
    print *, 'done'
    !call check( nf90_def_dim( ncid=ncid_snapout, name="processors", len=nproc, &
    !                          dimid=nc_proc_dimid) )

    do ivar=1, nvar ! The big snapshot variables for the kerner.
        call check( nf90_def_dim(ncid=ncout_snap_grpid, &
                                 name="dim_"//trim(varnamelist(ivar)), &
                                 len=varlength(ivar), dimid=ncout_field_dimid(ivar)) )

        call check( nf90_def_var(ncid=ncout_snap_grpid, name=trim(varnamelist(ivar)), &
                                 xtype=NF90_FLOAT, &
                                 dimids=(/ncout_field_dimid(ivar), ncout_snap_dimid/),&
                                 varid=ncout_field_varid(ivar), &
                                 chunksizes = (/1, nsnap/)) )

        !if (deflate) then
        !    call check( nf90_def_var_deflate(ncid=ncid_snapout, &
        !                                     varid=nc_field_varid(ivar), &
        !                                     shuffle=1, deflate=1, &
        !                                     deflate_level=deflate_level) )
        !end if
        !call check( nf90_def_var_fill(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
        !                              no_fill=1, fill=0) )
        write(6,"('  Netcdf variable ',A16,' with ID ', I3, ' and length', I8, ' produced.')") &
            trim(varnamelist(ivar)), ncout_field_varid(ivar), varlength(ivar)
       
    end do
    call check( nf90_enddef(ncout_id))


    !MvD: I guess the failure here is due to the multiprocessors varlength -> varlength * nprocs
    
    call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(1), values=strain_dsus_sol, &
                             start=(/1, 1/), count=(/varlength(1), nsnap/)) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(2), values = strain_dsuz_sol) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(3), values = strain_dpup_sol) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(4), values = straintrace_sol) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(6), values = strain_dsus_flu) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(7), values = strain_dsuz_flu) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(8), values = strain_dpup_flu) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(9), values = straintrace_flu) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(5), values=velo_sol) ) 
    !call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(10), values=velo_flu) )
    !call check( nf90_close(ncin_id))

    !
    !
    !allocate(ft_strain_dsus_sol(varlength(1), nomega))
    !allocate(ft_strain_dsuz_sol(varlength(2), nomega))
    !allocate(ft_strain_dpup_sol(varlength(3), nomega))
    !allocate(ft_straintrace_sol(varlength(4), nomega))
    !allocate(ft_strain_dsus_flu(varlength(6), nomega))
    !allocate(ft_strain_dsuz_flu(varlength(7), nomega))
    !allocate(ft_strain_dpup_flu(varlength(8), nomega))
    !allocate(ft_straintrace_flu(varlength(9), nomega))
    !allocate(ft_velo_sol(varlength(5),  nomega))
    !allocate(ft_velo_sol(varlength(10), nomega))

    !!! PLACE FOURIER TRANSFORM HERE
    !nextpow2 = 2
    !do while (nextpow2 < nsnap) 
    !    nextpow2 = nextpow2 * 2
    !end do
    !print *, 'nextpow2: ', nextpow2

    !nomega = nextpow2 / 2 + 1



    !call check( nf90_put_var(ncid = ncid_snapout, varid = nc_field_varid_out(1),  & 
    !                         values = ft_strain_dsus_sol, start=(/1,1/)) )
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(2),  & 
    !    values = ft_strain_dsuz_sol) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(3),  & 
    !    values = ft_strain_dpup_sol) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(4),  & 
    !    values = ft_straintrace_sol) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(6),  & 
    !    values = ft_strain_dsus_flu) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(7),  & 
    !    values = ft_strain_dsuz_flu) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(8),  & 
    !    values = ft_strain_dpup_flu) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(9),  & 
    !    values = ft_straintrace_flu) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(5),  & 
    !    values = ft_velo_sol) ) 
    !call check( nf90_put_var(ncid_snapout, nc_field_varid_out(10), & 
    !    values = ft_velo_flu) ) 
    !call check( nf90_close(ncid_snapout))

    
    
    
    
    
    contains

!> Translates NetCDF error code into readable message
subroutine check(status)
    implicit none
    integer, intent ( in) :: status !< Error code
    if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
    end if
end subroutine check  

end program
