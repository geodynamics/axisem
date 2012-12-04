program field_transformation
    use netcdf
    implicit none
    integer                         :: ncid_in, nsnap, nproc, nvar, ncid_snapin, nmode, ivar
    integer                         :: ncid_snapout, ncid_out, nelemperstep, nc_snap_dimid
    integer                         :: nc_proc_dimid
    character(len=8)                :: sourcetype
    integer                         :: dimids(14,3)
    character(len=16), allocatable  :: varnamelist(:)
    integer, dimension(14)          :: varlength, nc_field_varid, nc_field_varid_out, nc_f_dimid
    logical, parameter              :: deflate = .false.
    integer, parameter              :: deflate_level = 0
    real, dimension(:,:,:), allocatable :: strain_dsus_sol, strain_dsuz_sol, strain_dpup_sol, &
                                           straintrace_sol, velo_sol       , strain_dsus_flu, &
                                           strain_dsuz_flu, strain_dpup_flu, straintrace_flu, &
                                           velo_flu       

    
    
    call check( nf90_open(path="test/Data/axisem_output.nc4", & 
                          mode=NF90_NOWRITE, ncid=ncid_in) )


    call check( nf90_inq_grp_ncid(ncid_in, "Snapshots", ncid_snapin) )
    call check( nf90_get_att(ncid_in, NF90_GLOBAL, "excitation type", sourcetype))
    
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
        stop
    end if
    do ivar = 1,nvar
        call check( nf90_inq_varid( ncid_snapin, varnamelist(ivar), &
                    varid=nc_field_varid(ivar)) )
        call check( nf90_inquire_variable(ncid_snapin, nc_field_varid(ivar), dimids=dimids(ivar,:) ))
        call check( nf90_inquire_dimension(ncid_snapin, dimids(ivar,1), len=varlength(ivar)))
    end do

    call check( nf90_inquire_dimension(ncid_snapin, dimids(ivar,3), len=nsnap))
    call check( nf90_inquire_dimension(ncid_snapin, dimids(ivar,2), len=nproc))
  
    nelemperstep = 10000
    allocate(strain_dsus_sol(varlength(1),nproc,nsnap))
    allocate(strain_dsuz_sol(varlength(2),nproc,nsnap))
    allocate(strain_dpup_sol(varlength(3),nproc,nsnap))
    allocate(straintrace_sol(varlength(4),nproc,nsnap))
    allocate(strain_dsus_flu(varlength(6),nproc,nsnap))
    allocate(strain_dsuz_flu(varlength(7),nproc,nsnap))
    allocate(strain_dpup_flu(varlength(8),nproc,nsnap))
    allocate(straintrace_flu(varlength(9),nproc,nsnap))
    allocate(velo_sol(varlength(5),  nproc, nsnap))
    allocate(velo_sol(varlength(10), nproc, nsnap))
        

    !! Create output file
    nmode = ior(NF90_CLOBBER,NF90_NETCDF4)
    call check( nf90_create ( path="ordered_output.nc4", &
                              cmode=nmode, ncid=ncid_out))
    write(6,*) '  Producing groups for Seismograms and Snapshots'
    !call check( nf90_def_grp( ncid_out, "Seismograms", ncid_recout) )
    call check( nf90_def_grp( ncid_out, "Snapshots", ncid_snapin) )
    !call check( nf90_def_grp( ncid_out, "Surface", ncid_surfout) )


    write(6,*) 'Define variables in ''Snapshots'' group of NetCDF output file'
    write(6,*) '  awaiting', nsnap, ' snapshots'
    call check( nf90_def_dim( ncid=ncid_snapout, name="snapshots", len=nsnap, &
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
                                 varid = nc_field_varid_out(ivar), &
                                 chunksizes = (/1, 1, nsnap/) ))

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
    call check( nf90_enddef(ncid_out))


    call check( nf90_get_var(ncid_snapin,  nc_field_varid(1), values=strain_dsus_sol) ) 
    call check( nf90_put_var(ncid_snapout, nc_field_varid(1), values=strain_dsus_sol) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(2), values=strain_dsuz_sol) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(3), values=strain_dpup_sol) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(4), values=straintrace_sol) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(6), values=strain_dsus_flu) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(7), values=strain_dsuz_flu) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(8), values=strain_dpup_flu) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(9), values=straintrace_flu) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(5), values=velo_sol) ) 
    call check( nf90_get_var(ncid_snapin, nc_field_varid(10), values=velo_flu) ) 
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
