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

module data_gllmesh
  implicit none
  public 
  double precision, dimension(:,:,:), allocatable :: sgll,zgll
  double precision, dimension(:,:,:), allocatable :: sgll_fluid,zgll_fluid
  double precision, dimension(:,:,:), allocatable :: sgll_solid,zgll_solid

  double precision :: hmin_glob, hmax_glob ! global min/max gll spacing
  double precision :: min_distance_dim     ! 0.1*hmin_glob [in meters]
  double precision :: min_distance_nondim  ! 0.1*hmin_glob [referenced to 1]
  
  double precision :: char_time_max
  integer          :: char_time_max_globel
  double precision :: char_time_max_rad, char_time_max_theta
  double precision :: char_time_min
  integer          :: char_time_min_globel
  double precision :: char_time_min_rad, char_time_min_theta

end module data_gllmesh
