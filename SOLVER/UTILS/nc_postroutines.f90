module nc_routines
contains
  subroutine nc_open(dirname, ncid_in_grp, nc_disp_varid)
    use netcdf
    implicit none
    integer, intent(out)          :: ncid_in_grp(4), nc_disp_varid(4)
    character(len=*), intent(in)  :: dirname
    integer                       :: ncid_in(4)

    !status = N
    write(6,*) 'Opening forward file ', trim(dirname)//'/MZZ/Data/axisem_output.nc4'
    call check( nf90_open ( path=trim(dirname)//'/MZZ/Data/axisem_output.nc4', & 
                & mode=NF90_NOWRITE, ncid=ncid_in(1)))
    call check( nf90_inq_grp_ncid(ncid_in(1), "Seismograms", ncid_in_grp(1)) )
    call check( nf90_inq_varid(ncid_in_grp(1), 'displacement', nc_disp_varid(1)) )

    write(6,*) 'Opening forward file ', trim(dirname)//'/MZZ/Data/axisem_output.nc4'
    call check( nf90_open ( path=trim(dirname)//'/MXX_P_MYY/Data/axisem_output.nc4', & 
                & mode=NF90_NOWRITE, ncid=ncid_in(2)))
    call check( nf90_inq_grp_ncid(ncid_in(2), "Seismograms", ncid_in_grp(2)) )
    call check( nf90_inq_varid(ncid_in_grp(2), 'displacement', nc_disp_varid(2)) )
    
    write(6,*) 'Opening forward file ', trim(dirname)//'/MZZ/Data/axisem_output.nc4'
    call check( nf90_open ( path=trim(dirname)//'/MXZ_MYZ/Data/axisem_output.nc4', & 
                & mode=NF90_NOWRITE, ncid=ncid_in(3)))
    call check( nf90_inq_grp_ncid(ncid_in(3), "Seismograms", ncid_in_grp(3)) )
    call check( nf90_inq_varid(ncid_in_grp(3), 'displacement', nc_disp_varid(3)) )
    
    write(6,*) 'Opening forward file ', trim(dirname)//'/MZZ/Data/axisem_output.nc4'
    call check( nf90_open ( path=trim(dirname)//'/MXY_MXX_M_MYY/Data/axisem_output.nc4', & 
                & mode=NF90_NOWRITE, ncid=ncid_in(4)))
    call check( nf90_inq_grp_ncid(ncid_in(4), "Seismograms", ncid_in_grp(4)) )
    call check( nf90_inq_varid(ncid_in_grp(4), 'displacement', nc_disp_varid(4)) )

  end subroutine nc_open

  subroutine nc_close(ncid_in_grp)
    use netcdf
    implicit none
    integer, intent(in) :: ncid_in_grp(4)
    integer             :: ncid_in, isim
    do isim = 1,4
      call check( nf90_inq_grp_parent( ncid_in_grp(isim), ncid_in) ) 
      call check( nf90_close(ncid_in))
    end do
  end subroutine


!--------------------------------------------------------------------
  subroutine nc_read_seis(ncid_in,nc_disp_varid,ind_rec,nt,seis_snglcomp)
    use netcdf
    implicit none
    integer, intent(in)            :: nt,ind_rec, ncid_in(4), nc_disp_varid(4)
  !  character(len=200), intent(in) :: dirname
    real, intent(out)              :: seis_snglcomp(nt,3,4)
    integer                        :: isim
  !  character(len=4)               :: ind_recchar

    do isim = 1,2  ! Monopole sources
    call check( nf90_get_var( ncid_in(isim), nc_disp_varid(isim), start=(/1, 1, ind_rec/), &
               & count = (/nt, 1, 1/), values = seis_snglcomp(:,1,isim)) )
    call check( nf90_get_var( ncid_in(isim), nc_disp_varid(isim), start=(/1, 2, ind_rec/), &
               & count = (/nt, 1, 1/), values = seis_snglcomp(:,3,isim)) )
    end do
    do isim = 3,4 ! Dipole sources
      call check( nf90_get_var( ncid_in(isim), nc_disp_varid(isim), start=(/1, 1, ind_rec/), &
                 & count = (/nt, 3, 1/), values = seis_snglcomp(:,:,isim)) )
    end do
  end subroutine nc_read_seis
!--------------------------------------------------------------------
  subroutine check(status)
  ! Translates netcdf error codes into error messages
    use netcdf
    implicit none
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine check  
end module nc_routines
