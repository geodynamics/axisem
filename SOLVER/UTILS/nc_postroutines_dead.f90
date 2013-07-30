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
