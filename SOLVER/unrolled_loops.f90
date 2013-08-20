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

!========================
!> Routines for general matrix-matrix and matrix-vector multiplication. Called a 
!! bazillion times, presumably fast.
module unrolled_loops
!========================

  implicit none
  public :: mxm, vxm

  contains

!============================================================

!-------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! Size is fixed to npol x npol
pure subroutine mxm(a,b,c)

  use data_mesh, only: npol
  use global_parameters, only: realkind
  !include "mesh_params.h" 
  
  real(kind=realkind), intent(in)  :: a(0: ,0: ),b(0: ,0: ) !< Input matrices
  real(kind=realkind), intent(out) :: c(0: ,0: )                    !< Result
  integer                          :: i, j

  do j = 0, npol
     do i = 0, npol
        c(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end subroutine mxm
!-------------------------------------------------------------

!-------------------------------------------------------------
!> Multiplies vector a leftwise to matrix b to have vector c.
!! Size is fixed to npol x npol
pure subroutine vxm(a,b,c)
 
  use data_mesh, only: npol
  use global_parameters, only: realkind
  !include "mesh_params.h" 

  real(kind=realkind), intent(in)  :: a(0: )         !< Vector a
  real(kind=realkind), intent(in)  :: b(0: ,0: ) !< Matrix b
  real(kind=realkind), intent(out) :: c(0: )         !< Resulting vector c
  integer                          :: j

  do j = 0, npol
     c(j) = sum(a * b(:,j))
  end do 

end subroutine vxm
!-------------------------------------------------------------

!========================
end module unrolled_loops
!========================
