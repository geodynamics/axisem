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
module data_mesh
!===================

  ! Arrays here pertain to some sort of mesh peculiarities and mainly serve 
  ! as information or parameters for many "if"-decisions such as
  ! axis, north, element type,  solid-fluid boundary mapping, coarsening,
  ! and specifically also related to the background model such as
  ! solid-fluid boundary mapping, discontinuities, and element arrays to 
  ! map solid/fluid to global element domains.
  ! These quantities are in active memory throughout the simulation.
  ! 
  ! Any global arrays containing properties inside elements are defined in data_matr.
  
  use global_parameters
  implicit none
  include "mesh_params.h"
  public 

  ! Misc definitions
  integer                          :: nsize, npoint_solid3
  logical                          :: do_mesh_tests

  ! global numbering array for the solid and fluid assembly
  real(kind=realkind), allocatable :: gvec_solid(:) 
  real(kind=realkind)              :: gvec_fluid(nglob_fluid)

  real(kind=realkind)              :: smallval

  ! Deprecated elemental mesh (radius & colatitude of elemental midpoint)
  ! This is used for blow up localization. Might want to remove this when 
  ! everything is running smoothly in all eternities...
  real(kind=realkind)   :: mean_rad_colat_solid(nel_solid,2)
  real(kind=realkind)   :: mean_rad_colat_fluid(nel_fluid,2)
  
  ! Global mesh informations
  real(kind=dp)         :: router ! Outer radius (surface)

  ! critical mesh parameters (spacing/velocity, characteristic lead time etc)
  real(kind=dp)         :: pts_wavelngth
  real(kind=dp)         :: hmin_glob, hmax_glob
  real(kind=dp)         :: min_distance_dim, min_distance_nondim
  real(kind=dp)         :: char_time_max
  integer               :: char_time_max_globel
  real(kind=dp)         :: char_time_max_rad, char_time_max_theta
  real(kind=dp)         :: char_time_min
  integer               :: char_time_min_globel
  real(kind=dp)         :: char_time_min_rad, char_time_min_theta
  real(kind=dp)         :: vpmin, vsmin, vpmax, vsmax
  real(kind=dp)         :: vpminr, vsminr, vpmaxr, vsmaxr
  integer, dimension(3) :: vpminloc, vsminloc, vpmaxloc, vsmaxloc

  !----------------------------------------------------------------------
  ! Axial elements
  integer               :: naxel, naxel_solid, naxel_fluid
  integer,allocatable   :: ax_el(:), ax_el_solid(:), ax_el_fluid(:)
  logical               :: axis_solid(nel_solid)
  logical               :: axis_fluid(nel_fluid)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Background Model related
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solid-Fluid boundary----------------------------------------------------

  ! mapping nel_bdry elements into the solid/fluid element numbers
  integer, dimension(nel_bdry) :: bdry_solid_el, bdry_fluid_el

  ! mapping the z coordinate of the boundary for each element 
  ! (depending on north/south, above/below)
  integer, dimension(nel_bdry) :: bdry_jpol_solid, bdry_jpol_fluid

  ! Boolean to determine whether proc has solid-fluid boundary elements
  logical :: have_bdry_elem 

  ! integer array of size nel_bdry containing the "global" element number 
  ! for 1:nel_bdry 
  integer, dimension(nel_bdry) :: ibdryel

  ! Background model--------------------------------------------------------
  character(len=100)          :: bkgrdmodel
  character(len=100)          :: meshname
  logical                     :: resolve_inner_shear, have_fluid
  real(kind=dp)               :: discont(ndisc)
  logical                     :: solid_domain(ndisc)
  integer                     :: idom_fluid(ndisc)
  real(kind=dp)               :: rmin, minh_ic, maxh_ic, maxh_icb
  logical                     :: make_homo
  real(kind=dp)               :: vphomo, vshomo, rhohomo
  logical                     :: anel_true ! anelastic model?
  !--------------------------------------------------------------------------

  ! Receiver locations
  integer                      :: maxind, num_rec, num_surf_el, num_rec_tot
  integer, allocatable         :: surfelem(:), jsurfel(:), surfcoord(:)
  integer                      :: ielepi, ielantipode, ielequ
  integer, allocatable         :: recfile_el(:,:), loc2globrec(:)
  logical                      :: have_epi, have_equ, have_antipode
  real                         :: dtheta_rec
  
  ! CMB receivers (same as receivers, just above CMB instead)
  integer                      :: num_cmb
  integer, allocatable         :: cmbfile_el(:,:), loc2globcmb(:)
  !--------------------------------------------------------------------------

  integer                      :: nelem_plot, npoint_plot
  logical, allocatable         :: plotting_mask(:,:,:)
  integer, allocatable         :: mapping_ijel_iplot(:,:,:)

!=======================
end module data_mesh
!=======================
