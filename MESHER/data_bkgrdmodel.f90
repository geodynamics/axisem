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

!===================
module data_bkgrdmodel
!===================

  implicit none
  
  integer                       :: ndisc, nfluidregions
  integer, allocatable          :: idom_fluid(:)
  double precision, allocatable :: discont(:)
  double precision, allocatable :: vp(:,:), vs(:,:), rho(:,:)
  logical, allocatable          :: solid_domain(:)
  integer                       :: lfbkgrdmodel
  character(len=100)            :: bkgrdmodel, fnam_ext_model
  logical                       :: resolve_inner_shear, have_fluid, have_solid
  double precision              :: pts_wavelngth
  double precision              :: period, courant
  double precision              :: dt
  integer                       :: nc_init, nproc_target
  
  ! the sole quantities to be created in create_subregions
  ! that are needed by the rest of the mesher
  integer                       :: nz_glob, ns_glob, nc_glob
  integer, allocatable          :: iclev_glob(:)
  double precision, allocatable :: dz_glob(:)
  double precision              :: rmin, minh_ic, maxh_ic, maxh_icb
  double precision              :: minhvp, maxhvs, maxhnsicb

!===================
end module data_bkgrdmodel
!===================
