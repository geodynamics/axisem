module nc_routines
contains
  subroutine nc_open(dirname, ncid_in_grp, nc_disp_varid)
    implicit none
    integer, intent(out)          :: ncid_in_grp(4), nc_disp_varid(4)
    character(len=*), intent(in)  :: dirname

  end subroutine nc_open

  subroutine nc_close(ncid_in_grp)
    implicit none
    integer, intent(in) :: ncid_in_grp(4)
  end subroutine


!--------------------------------------------------------------------
  subroutine nc_read_seis(ncid_in,nc_disp_varid,ind_rec,nt,seis_snglcomp)
    implicit none
    integer, intent(in)            :: nt,ind_rec, ncid_in(4), nc_disp_varid(4)
  !  character(len=200), intent(in) :: dirname
    real, intent(out)              :: seis_snglcomp(nt,3,4)
  end subroutine nc_read_seis
!--------------------------------------------------------------------
end module nc_routines
