!========================
module time_evol_wave
!========================

  use global_parameters
  use data_proc  
  use data_mesh
  use data_source
  use data_time
  use seismograms
  use rotations 
  
  implicit none
  public :: prepare_waves, time_loop
  private

contains
 
!-----------------------------------------------------------------------------
subroutine prepare_waves
  !
  ! Contains all the preliminaries to propagate waves; such as the 
  ! background model, the stiffness and mass terms, the source and receiver 
  ! parameters, and preparations for I/O (dumping meshes, opening files).
  ! 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use data_io
  use parameters
  use def_precomp_terms
  use source
  use clocks_mod
  use meshes_io
  use attenuation, only: dump_memory_vars
    
  character(len=120) :: fname

  if (lpr) then
     write(6,*)''
     write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     if (have_fluid) then
        write(6,*)'++++++++    SEISMIC WAVE PROPAGATION: SOLID-FLUID CASE  ++++++++'
     else 
        write(6,*)'+++++++++++  SEISMIC WAVE PROPAGATION: SOLID CASE  +++++++++++++'
     endif
     write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif 

  ! read source parameters from sourceparams.dat
  call read_sourceparams

  ! rotation if source is not on the axis
  if (rot_src ) call def_rot_matrix

  ! Define velocity/density model (velocities in m/s, density in kg/m^3 ) AND 
  ! compute all global matrices (Jacobian, mapping, mass matrix, S/F boundary)
  ! for solid and fluid domains respectively
  call read_model_compute_terms

  ! compute/output some more parameters 
  call compute_numerical_parameters

  ! compute source time function
  call compute_stf

  ! compute source location within mesh and moment tensor/single force field
  call compute_src

  ! Prepare output
  if ( dump_energy .and. lpr) then ! only one proc dumps
     write(6,*)'  opening files for kinetic/potential energy...'
     if (have_fluid) then
        open(unit=4444,file=datapath(1:lfdata)//'/energy_sol.dat')
        open(unit=4445,file=datapath(1:lfdata)//'/energy_flu.dat')
     endif
     open(unit=4446,file=datapath(1:lfdata)//'/energy_glob.dat')
  endif     

  if (dump_wavefields) then 
     if (lpr) write(6,*)'  dumping strain mesh and associated fields...'
     call dump_wavefields_mesh_1d
  endif

  if (dump_snaps_glob) then
     if (lpr) write(6,*)'  dumping global grids for snapshots...'
     call dump_glob_grid_midpoint(ibeg,iend,ibeg,iend)
  endif

  if (dump_xdmf) then
     if (lpr) write(6,*)'  dumping mesh for xdmf snapshots...'
     call dump_xdmf_grid()

     fname = datapath(1:lfdata)//'/xdmf_snap_s_' //appmynum//'.dat'
     open(13100, file=trim(fname), access='stream', status='unknown', &
         convert='little_endian', position='append')

     if (.not. src_type(1)=='monopole') then
         fname = datapath(1:lfdata)//'/xdmf_snap_p_' //appmynum//'.dat'
         open(13101, file=trim(fname), access='stream', status='unknown', &
             convert='little_endian', position='append')
     endif

     fname = datapath(1:lfdata)//'/xdmf_snap_z_' //appmynum//'.dat'
     open(13102, file=trim(fname), access='stream', status='unknown', &
         convert='little_endian', position='append')
  endif

  if (anel_true .and. dump_memory_vars) &
     call prepare_mesh_memoryvar_vtk()

  if (dump_snaps_solflu) then
     if (lpr) write(6,*)'  dumping solid & fluid grids for snapshots...'
     call dump_solid_grid(ibeg,iend,ibeg,iend)
     if (have_fluid) call dump_fluid_grid(ibeg,iend,ibeg,iend)
  endif

  ! Various seismogram output preparations...
  call prepare_seismograms
  call open_hyp_epi_equ_anti

  ! allow for different types of receiver files
  call prepare_from_recfile_seis

  !if (rec_file_type=='colatlon' .and.  &
  !    (trim(bkgrdmodel(1:4))=='prem' .or. trim(bkgrdmodel(1:4))=='iasp')) &
  !   call prepare_from_recfile_cmb

  ! Need to reload old seismograms and add results
  if (isim>1 .and. sum_seis ) then  
     if (lpr) write(6,*)' Running multiple simulations and summing seismograms'
     if (lpr) write(6,*)' ...implementation of multiple simulations not finished'
     stop
  endif

  ! Need to reload old seismograms and add results
  if (isim>1 .and. sum_fields) then    
     if (lpr) write(6,*)' Running multiple simulations and summing wavefields'
     if (lpr) write(6,*)' ...implementation of multiple simulations not finished'
     stop
  endif

  ! Specific file format defined with Karin Sigloch: delivergf.dat
  if (src_file_type=='deliverg') then
     ! MvD: what is supposed to happen here? 
     stop
  endif

  ! write out seismic & numerical information on the simulation
  ! and run some tests on consistency of mesh/spacing/element types/messaging
  call write_parameters

  write(6,*) procstrg, 'done preparing waves.'

end subroutine prepare_waves
!=============================================================================

!-----------------------------------------------------------------------------
subroutine time_loop

  use clocks_mod, ONLY: tick

  iclockold = tick()

  if (time_scheme=='newmark2') then
     call sf_time_loop_newmark
  else
     call symplectic_time_loop
  endif

  iclockold = tick(id=idold, since=iclockold)
 
end subroutine time_loop
!=============================================================================


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           T I M E   E X T R A P O L A T I O N   R O U T I N E S
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!-----------------------------------------------------------------------------
subroutine sf_time_loop_newmark
  !
  ! The conventional explicit, acceleration-driven Newmark scheme of 2nd order.
  ! (e.g. Chaljub & Valette, 2004). The mass matrix is diagonal; we only store 
  ! its pre-assembled inverse at the stage of the time loop.
  ! Explicit axial masking follows Nissen-Meyer et al. 2007, GJI, 
  ! "Spherical-earth Frechet sensitivity kernels" eqs. (80)-(82).
  ! Note that the ordering (starting inside the fluid) is crucial such that no 
  ! iterations for the boundary terms are necessary.
  ! Also note that our definition of the fluid potential is different from 
  ! Chaljub & Valette and the code SPECFEM by an inverse density factor.
  ! This is the correct choice for our case of non-gravitating Earth models, 
  ! but shall be altered once gravity is taken into account.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use commun
  use global_parameters
  use apply_masks
  use stiffness
  use unit_stride_colloc
  use clocks_mod
  use data_matr,            ONLY: inv_mass_rho, inv_mass_fluid
  use attenuation,          ONLY: n_sls_attenuation, time_step_memvars
  
  include 'mesh_params.h'
  
  ! Solid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp, velo
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: acc0, acc1
  
  ! solid memory variables + gradient
  real(kind=realkind), allocatable :: memory_var(:,:,:,:,:)

  ! Fluid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: chi, dchi
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: ddchi0, ddchi1
  
  integer :: iter

  if (lpr) then
     write(6,*)
     write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
     write(6,*)'TTTT  2nd-order, acceleration-driven Newmark time scheme TTTTT'
     write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
     write(6,*)
  endif

  if (anel_true) then
     allocate(memory_var(0:npol,0:npol,6,n_sls_attenuation,nel_solid))
     memory_var = 0
  endif

  ! INITIAL CONDITIONS
  ! initializiation with a small number prevents performance loss (~factor 3) due to
  ! denormal floats in the onset of the p-wave (going from zero to some finite value)
  ! alternatively, compilerflags -ffast-math (gfortran) or -ftz (ifort) might
  ! give the same speedup, but seem to be unstable on some systems
  ! another alternative:
  ! http://software.intel.com/en-us/articles/how-to-avoid-performance-penalties-for-gradual-underflow-behavior
  ! some more reading:
  ! http://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x

  !disp = 1.d-30
  !velo = 1.d-30
  !acc0 = 1.d-30
  !acc1 = 1.d-30
  !
  !chi = 1.d-30
  !dchi = 1.d-30
  !ddchi0 = 1.d-30
  !ddchi1 = 1.d-30

  disp = zero
  velo = zero 
  acc0 = zero
  acc1 = zero
  
  chi = zero
  dchi = zero
  ddchi0 = zero
  ddchi1 = zero
  
  t = zero

  if (lpr) write(6,*)'************ S T A R T I N G   T I M E   L O O P *************'
  write(69,*)'************ S T A R T I N G   T I M E   L O O P *************'

  do iter = 1, niter

     t = t + deltat
     call runtime_info(iter,disp,chi)

     ! ::::::::::::::::::::::::: ACTUAL NEWMARK SOLVER :::::::::::::::::::::::::

     ! FLUID: new fluid potential
     chi = chi +  deltat * dchi + half_dt_sq * ddchi0

     ! SOLID: new displacement
     disp = disp + deltat * velo + half_dt_sq * acc0

     ! FLUID: Axial masking of \chi (dipole,quadrupole)
     if (src_type(1) .ne. 'monopole') & 
          call apply_axis_mask_scal(chi, nel_fluid, ax_el_fluid, naxel_fluid) 

     ! FLUID: stiffness term K_f
     iclockstiff = tick()
     call glob_fluid_stiffness(ddchi1, chi) 
     iclockstiff = tick(id=idstiff, since=iclockstiff)

     ! FLUID: solid-fluid boundary term (solid displ.) added to fluid stiffness
     call bdry_copy2fluid(ddchi1, disp)

     ! FLUID: axial masking of w (dipole,quadrupole)
     if (src_type(1) .ne. 'monopole') & 
          call apply_axis_mask_scal(ddchi1, nel_fluid, ax_el_fluid, naxel_fluid)

     ! FLUID: stiffness term assembly
     iclockcomm = tick()
     call comm2d(ddchi1, nel_fluid, 1, 'fluid')
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     ! FLUID: update 2nd derivative of potential
     ddchi1 = - inv_mass_fluid * ddchi1

     ! SOLID: For each source type, apply the following sequence respectively:
     !        1) masking of u 
     !        2) stiffness term K_s u
     !        3) solid-fluid bdry term (fluid potential) added to solid stiffness
     !        4) masking of w
     ! MvD: masking of w?? you mean acc?

     iclockstiff = tick()
     select case (src_type(1))
        case ('monopole')
           call apply_axis_mask_onecomp(disp, nel_solid, ax_el_solid, naxel_solid)
           call glob_stiffness_mono(acc1, disp)
           if (anel_true) then
              iclockanelst = tick()
              call glob_anel_stiffness_mono(acc1, memory_var)
              iclockanelst = tick(id=idanelst, since=iclockanelst)
           endif
           call bdry_copy2solid(acc1, ddchi1)
           call apply_axis_mask_onecomp(acc1, nel_solid, ax_el_solid, naxel_solid)

        case ('dipole') 
           call apply_axis_mask_twocomp(disp, nel_solid, ax_el_solid, naxel_solid)
           call glob_stiffness_di(acc1,disp) 
           call bdry_copy2solid(acc1,ddchi1)
           call apply_axis_mask_twocomp(acc1, nel_solid, ax_el_solid, naxel_solid)

        case ('quadpole') 
           call apply_axis_mask_threecomp(disp, nel_solid, ax_el_solid, naxel_solid)
           call glob_stiffness_quad(acc1,disp) 
           call bdry_copy2solid(acc1,ddchi1)
           call apply_axis_mask_threecomp(acc1, nel_solid, ax_el_solid, naxel_solid)
     end select
     iclockstiff = tick(id=idstiff, since=iclockstiff)

     ! SOLID: 3-component stiffness term assembly ==> w^T K_s u
     iclockcomm = tick()
     call comm2d(acc1, nel_solid, 3, 'solid') 
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     ! SOLID: add source, only in source elements and for stf/=0
     if (have_src) call add_source(acc1, stf(iter))

     ! SOLID: new acceleration (dipole has factor two due to (+,-,z) coord. system)
     acc1(:,:,:,1) = - inv_mass_rho * acc1(:,:,:,1)

     if (src_type(1)/='monopole') &
          acc1(:,:,:,2) = - inv_mass_rho * acc1(:,:,:,2)

     if (src_type(1)=='dipole') then
        ! for the factor 2 compare eq 32 in TNM (2006)
        acc1(:,:,:,3) = - two * inv_mass_rho * acc1(:,:,:,3)
     else
        acc1(:,:,:,3) = - inv_mass_rho * acc1(:,:,:,3)
     endif

     ! FLUID: new 1st derivative of potential
     dchi = dchi + half_dt * (ddchi0 + ddchi1)

     ! SOLID: new velocity
     velo = velo + half_dt * (acc0 + acc1)

     ! update acceleration & 2nd deriv. of potential
     ddchi0 = ddchi1
     acc0 = acc1

     ! ::::::::::::::::::::::::: END FD SOLVER ::::::::::::::::::::::::::
     
     ! memory variable time evolution with strain as source
     if (anel_true) then
        iclockanelts = tick()
        call time_step_memvars(memory_var, disp)
        iclockanelts = tick(id=idanelts, since=iclockanelts)
     endif

     iclockdump = tick()
     if (anel_true) then
        call dump_stuff(iter, disp, velo, chi, dchi, ddchi0, memory_var)
     else
        call dump_stuff(iter, disp, velo, chi, dchi, ddchi0)
     endif
     iclockdump = tick(id=iddump, since=iclockdump)

  end do ! time loop

end subroutine sf_time_loop_newmark
!=============================================================================

!-----------------------------------------------------------------------------
subroutine symplectic_time_loop
  !
  ! SOLVE coupled solid-fluid system of temporal ODEs:
  !   M*\dot{u}    = -K*u - B*\ddot{\chi} + F (solid)
  !   M*\ddot{chi} = -K*\chi - B*u (fluid)
  ! using symplectic time integration schemes of 4th, 6th, 8th, 10th order 
  !
  ! The time step can be chosen 1.5 times larger than in Newmark, resulting 
  ! in CPU times about 2.5 times longer than Newmark, but considerably more 
  ! accurate. Consult Ampuero & Nissen-Meyer (2007) for examples of when 
  ! this choice might be more appropriate. Generally, for long propagation 
  ! distances (say, > 100 wavelengths), it is worthwhile considering this.
  ! 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use global_parameters
  use commun
  use apply_masks
  use stiffness
  use clocks_mod
  use unit_stride_colloc
  use source,       ONLY: compute_stf_t
  use data_matr,    ONLY: inv_mass_rho,inv_mass_fluid
  
  include 'mesh_params.h'
  
  ! solid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp, velo, acc
  
  ! fluid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: chi, dchi, ddchi
  
  integer          :: iter, i
  
  ! symplectic stuff
  double precision, allocatable, dimension(:) :: coefd, coeff, coefv, subdt
  double precision, allocatable, dimension(:) :: stf_symp
  
  ! choose symplectic scheme (4,6,8,10th order) and compute coefficients
  call symplectic_coefficients(coefd,coeff,coefv)
  allocate(subdt(nstages),stf_symp(nstages))

  ! INITIAL CONDITIONS
  disp = zero
  velo = zero
  acc = zero

  chi = zero
  dchi = zero
  ddchi = zero

  t = zero
  if (lpr) write(6,*)'*********** S T A R T I N G   T I M E   L O O P ************'
  write(69,*)'*********** S T A R T I N G   T I M E   L O O P ************'

  do iter=1, niter

     t = t + deltat
     call runtime_info(iter,disp,chi)

     ! ::::::::::::::::::::::::: ACTUAL SYMPLECTIC SOLVER :::::::::::::::::::::::::

     ! compute external force/source time function at 4 time intervals coeff(1:4)
     subdt = t - deltat + coeff
     call compute_stf_t(nstages,subdt,stf_symp)

     do i = 1, nstages  ! substages 

        ! FLUID: update fluid potential
        chi = chi + dchi * coefd(i)

        ! SOLID: update displacement
        disp = disp + velo * coefd(i)

        ! FLUID: Axial masking of \chi (dipole,quadrupole)
        if (src_type(1) .ne. 'monopole') & 
             call apply_axis_mask_scal(chi,nel_fluid,ax_el_fluid,naxel_fluid) 

        ! FLUID: stiffness term K_f
        iclockstiff = tick()
        call glob_fluid_stiffness(ddchi,chi) 
        iclockstiff = tick(id=idstiff, since=iclockstiff)

        ! FLUID: solid-fluid boundary term (solid displ.) added to fluid stiffness
        call bdry_copy2fluid(ddchi,disp)

        ! FLUID: axial masking of w (dipole,quadrupole)
        if (src_type(1) .ne. 'monopole') & 
             call apply_axis_mask_scal(ddchi,nel_fluid,ax_el_fluid,naxel_fluid)

        ! FLUID: stiffness term assembly
        iclockcomm = tick()
        call comm2d(ddchi,nel_fluid,1,'fluid')
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        ! FLUID: update 1st derivative of potential
        dchi = dchi - coefv(i) * inv_mass_fluid * ddchi

        ! SOLID: For each source type, apply the following sequence respectively:
        !        1) masking of u 
        !        2) stiffness term K_s u
        !        3) solid-fluid bdry term (fluid potential) added to solid stiffness
        !        4) masking of w
        ! MvD: masking of w?? you mean acc?

        call collocate0_neg1d_existent(ddchi,inv_mass_fluid,npoint_fluid)

        iclockstiff = tick()
        select case (src_type(1))
           case ('monopole')
              call apply_axis_mask_onecomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_mono(acc,disp)
              call bdry_copy2solid(acc,ddchi)
              call apply_axis_mask_onecomp(acc,nel_solid, ax_el_solid,naxel_solid)

           case ('dipole') 
              call apply_axis_mask_twocomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_di(acc,disp) 
              call bdry_copy2solid(acc,ddchi)
              call apply_axis_mask_twocomp(acc,nel_solid, ax_el_solid,naxel_solid)

           case ('quadpole') 
              call apply_axis_mask_threecomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_quad(acc,disp) 
              call bdry_copy2solid(acc,ddchi)
              call apply_axis_mask_threecomp(acc,nel_solid, ax_el_solid,naxel_solid)
        end select
        iclockstiff = tick(id=idstiff, since=iclockstiff)

        ! SOLID: stiffness term assembly ==> w^T K_s u
        iclockcomm = tick()
        call comm2d(acc,nel_solid,3,'solid')
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        ! SOLID: add source, only in source elements and for stf/=0
        if (have_src) call add_source(acc,real(stf_symp(i),kind=realkind))

        ! SOLID: new acceleration (dipole has factor two due to (+,-,z) coord. system)
        velo(:,:,:,1) = velo(:,:,:,1) - acc(:,:,:,1) * coefv(i) * inv_mass_rho

        if (src_type(1)/='monopole') &
           velo(:,:,:,2) = velo(:,:,:,2) - acc(:,:,:,2) * coefv(i) * inv_mass_rho
        if (src_type(1)=='dipole') then !factor 2 b/c inv_rho has 1/2 embedded
           velo(:,:,:,3) = velo(:,:,:,3) - two * acc(:,:,:,3) * coefv(i) * inv_mass_rho
        else
           velo(:,:,:,3) = velo(:,:,:,3) - acc(:,:,:,3) * coefv(i) * inv_mass_rho
        endif

     enddo ! ... nstages substages

     ! FLUID: final potential
     chi = chi + dchi * coefd(nstages+1)

     ! SOLID: final displacement
     disp = disp + velo * coefd(nstages+1)

     ! ::::::::::::::::::::::::: END SYMPLECTIC SOLVER ::::::::::::::::::::::::::

     iclockdump = tick()
     call dump_stuff(iter,disp,velo,chi,dchi,ddchi)
     iclockdump = tick(id=iddump, since=iclockdump)

  end do ! time loop

end subroutine symplectic_time_loop
!=============================================================================

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     E N D   O F   T I M E   E X T R A P O L A T I O N   R O U T I N E S
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!-----------------------------------------------------------------------------
subroutine symplectic_coefficients(coefd,coeff,coefv)

  use commun,         ONLY : barrier,pend
  
  double precision, allocatable, dimension(:), intent(out) :: coefd,coeff,coefv
  double precision, allocatable, dimension(:) :: g
  double precision :: zeta_symp,iota_symp,kappa_symp
  double precision :: rho,theta,nu,lambda
  integer :: Q,n,i
  real :: B,C

  select case (time_scheme)
     
  !444444444444444444444444444444444444444444444444444444444444444444444444444
  case('symplec4') ! position extended Forest-Ruth like 
                   ! (Omelyan, Mryglod and Folk, 2002)
     nstages = 4
     allocate(coefd(nstages+1),coeff(nstages),coefv(nstages))

     if (lpr) then
        write(6,*)
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)'TTTT  4th-order symplectic PEFRL scheme  TTTTTTTTTTTTTTT'
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)
     endif
     
     ! symplectic parameters
     zeta_symp=+0.1786178958448091
     iota_symp=-0.2123418310626054
     kappa_symp=-0.06626458266981849
     
     ! extrapolation coefficients
     coefd(1:nstages+1) = (/zeta_symp, kappa_symp, &
          dble(1.d0-2.d0*(zeta_symp+kappa_symp)),  &
                 kappa_symp, zeta_symp/)*deltat
     
     coefv(1:nstages) = (/dble(half-iota_symp), iota_symp, iota_symp, &
          dble(half-iota_symp)/)*deltat

     Q = 4; B=1/12500.d0; C=2.97633; ! empirical


  case('ML_SO4m5')  ! Order 4, S, m=5 in Table 2 of McLachlan (1995), preferred

     nstages=5
     allocate(coefd(nstages+1),coeff(nstages),coefv(nstages))
     
     if (lpr) then
        write(6,*)
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)'TTTT  4th-order symplectic McLachlan scheme  TTTTTTTTTTT'
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)
     endif
     
     rho = (14.d0-dsqrt(19.d0))/108.d0;
     theta = (20.d0-7.d0*dsqrt(19.d0))/108.d0;
     nu = 2.d0/5.d0;
     lambda  = -1.d0/10.d0;

     coefd(1) = rho ; 
     coefd(2) = theta ; 
     coefd(3) = 1.d0/2.d0-rho-theta;  
     coefd(4) = 1.d0/2.d0-rho-theta
     coefd(5) = theta 
     coefd(6) = rho
     
     coefv(1) = nu 
     coefv(2) = lambda 
     coefv(3) = 1.d0-2.d0*(nu+lambda)
     coefv(4) = lambda
     coefv(5) = nu

     coefd = coefd*deltat; coefv = coefv*deltat 
     
     Q = 4; B=1.49e-05; C=3.035

  !6666666666666666666666666666666666666666666666666666666666666666666666666
  case('ML_SO6m7') ! best order 6 so far!
                   ! other order 6 are not better than the best order 4

     nstages=7
     allocate(coefd(nstages+1),coeff(nstages),coefv(nstages))
     
     if (lpr) then
        write(6,*)
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)'TTTT  6th-order symplectic McLachlan scheme  TTTTTTTTTTT'
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)
     endif

     coefd(1) = -1.01308797891717472981
     coefd(2) = 1.18742957373254270702 
     coefd(3) = -0.01833585209646059034
     coefd(4) = 0.34399425728109261313 
     do i=5,8
        coefd(i) = coefd(nstages+2-i)
     enddo

     coefv(1) = 0.00016600692650009894 
     coefv(2) = -0.37962421426377360608 
     coefv(3) = 0.68913741185181063674 
     coefv(4) = 0.38064159097092574080
     do i=5,7
        coefv(i) = coefv(nstages+1-i)
     enddo

     coefd = coefd*deltat; coefv = coefv*deltat; 

     Q = 6; B=1.3e-06; C=3.067;

  !888888888888888888888888888888888888888888888888888888888888888888888888
  case('KL_O8m17') ! Kahan and Li (1997), improved on McLachlan (1995)

     n = 8;  nstages = 2*n+1
     allocate(g(n),coefd(nstages+1),coeff(nstages),coefv(nstages))

     if (lpr) then
        write(6,*)
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)'TTTT  8th-order symplectic Kahan/Li scheme  TTTTTTTTTTTT'
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)
     endif

     g(1) = 0.13020248308889008088;
     g(2) = 0.56116298177510838456;
     g(3) = -0.38947496264484728641;
     g(4) = 0.15884190655515560090;
     g(5) = -0.39590389413323757734;
     g(6) = 0.18453964097831570709;
     g(7) = 0.25837438768632204729;
     g(8) = 0.29501172360931029887;
     
     call SS_scheme(n,coefd,coefv,g)

     coefd = coefd*deltat; coefv = coefv*deltat; 
     
     Q = 8; B=-100000.; C=3.; 
     
  !101010101010101010101010101010101010101010101010101010101010101010101010
  case('SS_35o10') ! Sofroniou and Spaletta (2004), best order 10

     n = 17; nstages = 2*n+1
     allocate(g(n),coefd(nstages+1),coeff(nstages),coefv(nstages))

     if (lpr) then
        write(6,*)
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)'TTTT  10th-order symplectic Sofroniou/Spaletta scheme TT'
        write(6,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(6,*)
     endif

    g(1) = 0.078795722521686419263907679337684;
    g(2) = 0.31309610341510852776481247192647;
    g(3) = 0.027918383235078066109520273275299;
    g(4) =-0.22959284159390709415121339679655;
    g(5) = 0.13096206107716486317465685927961;
    g(6) =-0.26973340565451071434460973222411;
    g(7) = 0.074973343155891435666137105641410;
    g(8) = 0.11199342399981020488957508073640;
    g(9) = 0.36613344954622675119314812353150;
    g(10) =-0.39910563013603589787862981058340;
    g(11) = 0.10308739852747107731580277001372;
    g(12) = 0.41143087395589023782070411897608;
    g(13) =-0.0048663605831352617621956593099771;
    g(14) =-0.39203335370863990644808193642610;
    g(15) = 0.051942502962449647037182904015976;
    g(16) = 0.050665090759924496335874344156866;
    g(17) = 0.049674370639729879054568800279461;

    call SS_scheme(n,coefd,coefv,g)

    coefd = coefd*deltat; coefv = coefv*deltat; 

    Q = 10; B=4.58e-10; C=5.973

  case default
     write(6,*)procstrg,'reporting ERROR ::::::::::'
     write(6,*)procstrg,time_scheme,'Time scheme unknown'; call pend; stop

  end select
 
  do i=1,nstages
     coeff(i) = sum(coefd(1:i))
  enddo

  call barrier
  if (mynum==0) then
     write(6,*)'  :::::::::::::::: Symplectic coefficients :::::::::::::::::'
     write(6,*)'   order,stages:',Q,nstages
     write(6,*)'   dispersion error coeff:',B
     write(6,*)'   CFL factor (wrt Newmark deltat):',C/2.
     do i=1,nstages
        write(6,12)i,coefd(i),coefv(i),coeff(i)
     enddo
     write(6,12)nstages+1,coefd(nstages+1)
     write(6,*)'  :::::::::::::: End Symplectic coefficients :::::::::::::::'
     write(6,*)
  endif
  call barrier

12 format('   ',i3,' coeffd,coeffv,sub_dt:',3(1pe10.2))

end subroutine symplectic_coefficients
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine SS_scheme(n,a,b,g)

  ! coefficients for symmetric compositions of symmetric methods
  
  integer, intent(in) :: n
  double precision, intent(in)  :: g(n)
  double precision, intent(out) :: a(nstages+1), b(nstages)
  integer :: i

  a(1) = g(1) / 2.d0
  a(2:n) = (g(1:n-1) + g(2:n)) / 2.d0
  a(n+1)= 1 / 2 - sum(a(1:n))
  do i=n+2, 2*n+2
     a(i) = a(2*n+3-i)
  enddo

  b(1:n) = g(1:n)
  b(n+1)= 1. - 2.d0 * sum(g)
  do i=n+2,2*n+1
     b(i) = b(2*n+2-i)
  enddo

end subroutine SS_scheme
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine runtime_info(iter,disp,chi)
  ! 
  ! Print time step, time, min/max displacement and potential values
  ! and stop the simulation if displacements blow up beyond acceptable...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use commun, only : pend
  
  include 'mesh_params.h'
  
  integer, intent(in)             :: iter
  real(kind=realkind), intent(in) :: disp(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(in) :: chi(0:npol,0:npol,nel_fluid)
  integer                         :: iblow(4), check_disp, check_iter, time_stamp
  character(len=4)                :: appistamp

  check_iter = 100 ! printing time/percentage done every check_iter-th time step.
  check_disp = 200 ! printing min/max displ. every check_disp-th time step.
  time_stamp = floor(real(niter)/100.) ! a file every time_stamp-th timestep.
  
  ! Stdout time step/percentage announcements
  if (lpr .and. mod(iter,check_iter)==0 ) then
     write(6,13) iter, t, real(iter) / real(niter) * 100.
     call flush(6)
13   format('  time step:',i6,'; t=',f8.2,' s (',f5.1,'%)')
  endif

  ! Time stamp to file
  if (lpr .and. mod(iter,time_stamp)==0) then 
     call define_io_appendix(appistamp, floor(real(iter)/real(time_stamp)))
     open(unit=100, file='timestamp'//appistamp//'.txt')
     write(100,13) iter, t, real(iter) / real(niter) * 100.
     close(100)
  endif

  ! Check on min/max. displacement/potential values globally
  if ( mod(iter,check_disp)==0 ) then
     if (iter==check_disp) then
        write(69,14) 'time', 'absmax(us)', 'absmax(up)', 'absmax(uz)', 'absmax(chi)'
     endif
     write(69,15) t, maxval(abs(disp(:,:,:,1))), maxval(abs(disp(:,:,:,2))), &
                  maxval(abs(disp(:,:,:,3))), maxval(abs(chi))
  endif
14 format(a7,4(a13))
15 format(f7.1,4(1pe12.3))

  ! Stop simulation if displacement exceeds source magnitude
  if ( maxval(abs(disp(1,1,:,:))) > abs(magnitude) ) then
     write(6,*) procstrg,'!!!!!!!!!!!!!!! DISPLACEMENTS BLEW UP !!!!!!!!!!!!!!!'
     write(6,*) procstrg,'  Time step & time:', iter, t 
     write(6,*) procstrg,'  Proc. num, displ value',mynum,maxval(abs(disp))
     iblow(1:4) = maxloc(abs(disp))
     write(6,*) procstrg,'iel,comp   :',iblow(3),iblow(4)
     write(6,*) procstrg,'elem r, th :',mean_rad_colat_solid(iblow(3),1), &
                                       mean_rad_colat_solid(iblow(3),2)
     write(6,*) procstrg,'ipol,jpol  :',iblow(1)-1,iblow(2)-1
     write(6,*) procstrg,'axis       :',axis_solid(iblow(3))
     write(6,*) procstrg,''
     call pend
     stop
  endif

end subroutine runtime_info
!=============================================================================

!-----------------------------------------------------------------------------
subroutine add_source(acc1,stf1)
  !
  ! Add source term inside source elements only if source time function non-zero
  ! and I have the source.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)    :: stf1
  real(kind=realkind), intent(inout) :: acc1(0:npol,0:npol,nel_solid,3)
  integer             :: iel,i

  i=0
  if ( have_src .and. stf1 /= zero) then 
     do iel=1, nelsrc
        i=i+1
        acc1(:,:,ielsrc(iel),:) = acc1(:,:,ielsrc(iel),:) - & 
             source_term_el(:,:,i,:)*stf1
     enddo
  endif

end subroutine add_source
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_stuff(iter, disp, velo, chi, dchi, ddchi, memvar)
  !
  ! Includes all output action done during the time loop such as
  ! various receiver definitions, wavefield snapshots, velocity field & strain 
  ! tensor for 3-D kernels, analytical radiation in a homogeneous model.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use data_io
  use data_mesh
  !use analyt_homog_radiation, ONLY : src_vicinity
  use wavefields_io
  use attenuation,          ONLY: n_sls_attenuation, dump_memory_vars
  
  integer, intent(in)            :: iter
  real(kind=realkind),intent(in) :: disp(0:npol,0:npol,nel_solid,3)
  real(kind=realkind),intent(in) :: velo(0:npol,0:npol,nel_solid,3)
  real(kind=realkind),intent(in) :: chi(0:npol,0:npol,nel_fluid)
  real(kind=realkind),intent(in) :: dchi(0:npol,0:npol,nel_fluid)
  real(kind=realkind),intent(in) :: ddchi(0:npol,0:npol,nel_fluid)
  real(kind=realkind),intent(in), optional :: &
        memvar(0:npol,0:npol,6,n_sls_attenuation,nel_solid)
  real(kind=realkind) :: time
  
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^- Time series^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  
  if ( mod(iter,seis_it)==0) then
     ! receiver locations read in from file (only 3-comp. displacements)

     iseismo = iseismo + 1
     if (use_netcdf) then
        call nc_compute_recfile_seis_bare(disp)
     else
        call compute_recfile_seis_bare(disp)
     endif
  
     time = real(iter)*deltat
     !call compute_recfile_seis_binary(time,disp,velo)
     !!andrea
     !if (dump_wavefields) call compute_surfelem_strain(disp)
     
     ! cmbrec locations read in from file (only velocity & tr(E))
     !    call dump_velo_straintrace_cmb(disp,velo)
     
     ! Generic synthetics at hypo-/epicenter, equator, antipode (including time)
     call compute_hyp_epi_equ_anti(t,disp)

  endif

  ! Analaytical radiation in homogeneous models-^-^-^-^-^-^-^-^-^^-^-^-^
  !  if (srcvic) call src_vicinity(iter,disp)


  ! Compute kinetic and potential energy globally every 5th time step
  !  if (dump_energy .and. mod(iter,5)==0) &
  if (dump_energy) &
       call energy(disp,velo,dchi,ddchi)

  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^ Wavefield snapshots-^-^^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  if (dump_snaps_glob) then
    if (mod(iter,snap_it)==0) then
       isnap = isnap + 1
       if (lpr) then
          write(6,*)
          write(6,*) 'Writing global snap to file: ', isnap
          write(6,*)
       endif
       call glob_snapshot_midpoint(disp, chi, ibeg, iend, ibeg, iend)
     endif
  endif
  
  if (dump_xdmf) then
    if (mod(iter,snap_it)==0) then
        isnap=isnap+1
        if (lpr) then
           write(6,*)
           write(6,*)'Writing global xdmf snap to file:',isnap
           write(6,*)
        endif
        call glob_snapshot_xdmf(disp, chi)
     endif
  endif

  if (anel_true .and. dump_memory_vars) then
    if (mod(iter,snap_it)==0 .or. iter==1) then
       call snapshot_memoryvar_vtk(memvar, iter)
    endif
  endif

  if (dump_snaps_solflu) then
    if (mod(iter,snap_it)==0) then
        isnap=isnap+1
        if (lpr) then
           write(6,*)
           write(6,*)'Writing solid/fluid snap to file:',isnap
           write(6,*)
        endif
       if (have_fluid) call fluid_snapshot(chi, ibeg, iend, ibeg, iend)
       call solid_snapshot(disp, ibeg, iend, ibeg, iend)
     endif
  endif

  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  ! Velocity field and strain tensor wavefields for 3-D kernels^-^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  !
  ! At this point, we offer two end-member versions of dumping the fields  
  ! to eventually calculate waveform kernels.
  !
  ! The FIRST one ('displ_only') dumps a minimal amount and requires extensive 
  ! post-processing when calculating the kernels, but optimizes the SEM 
  ! simulation in terms of memory, storage amount and CPU time.
  ! Note that this method IS ONLY POSSIBLE IF ENTIRE SEM MESH IS DUMPED!!!
  !
  ! The SECOND one ('fullfields') computes the entire strain tensor and velocity 
  ! field on-the-fly, resulting in more output (9 rather than 6 fields), 
  ! more memory and CPU time during the SEM, but no post-processing necessary. 
  ! Any kind of spatial distribution can be dumped, meaning in the long run 
  ! this should be the more effective choice.
  !
  ! Currently, 'fullfields' is the hardcoded choice (parameters.f90:110)
  ! 
  ! Possible cases in between these dumpsters are considered but not yet 
  ! implemented (e.g. dumping 1/s, inverse fluid density, but compute derivatives
  ! on-the-fly to warrant grid flexibility).

  if (dump_wavefields) then

    if (mod(iter,strain_it)==0) then

     ! dump displacement and velocity in each surface element
     !! for netcdf people set .true. in inparam to use it instead of the standard
     
     !!the update of the strain has to preceed the call to the function. 
     !!It starts from 
      istrain=istrain+1

      call compute_surfelem(disp,velo)
       
      select case (dump_type)

        case ('displ_only')
          ! Only dump the 3-comp displacement and velocity fields in solid 
          ! and potential & its derivative in fluid.
          ! Minimal permanent storage, minimal run-time memory, minimal CPU time, 
          ! but extensive post-processing (need to compute strain tensor).
             call dump_disp(disp,chi)       ! displacement in solid, chi in fluid
             call dump_velo_dchi(velo,dchi) ! velocity in solid, dchi in fluid

        case ('fullfields') ! Hardcoded choice
          ! Compute strain tensor on-the-fly here and dump the 6 components.
          ! Also compute corresponding fields in the fluid.
          ! Maximal permanent storage, maximal run-time memory, maximal CPU time, 
          ! but no post-processeing necessary as these are the fields that 
          ! constitute density and elastic kernels.
            call compute_strain(disp,chi)    ! strain globally
            call dump_velo_global(velo,dchi) ! velocity globally

        end select
       
        !> Check, whether it is time to dump the buffer variables to disk and if so,
        !! do so.
        if (use_netcdf) call nc_dump_strain(istrain)

    endif ! dumping interval strain_it

endif   ! dump_wavefields?

end subroutine dump_stuff
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_velo_straintrace_cmb(u,velo)
  ! 
  ! This is for quick checks and singular computations of the bulk moduli kernels
  ! hence only the trace of the strain tensor is needed and computed. 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  use data_source,              ONLY: src_type
  use pointwise_derivatives,    ONLY: axisym_gradient_solid
  use pointwise_derivatives,    ONLY: f_over_s_solid
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)   :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(in)   :: velo(0:npol,0:npol,nel_solid,3)
  real(kind=realkind)               :: grad_sol(0:npol,0:npol,nel_solid,3)
  real(kind=realkind)               :: dsdf(0:npol,naxel_solid)
  integer                           :: iel

  if (src_type(1)=='dipole') then
     call axisym_gradient_solid(u(:,:,:,1)+u(:,:,:,2),grad_sol(:,:,:,1:2))
  else
     call axisym_gradient_solid(u(:,:,:,1),grad_sol(:,:,:,1:2)) ! dsus, dzus
  endif
  call axisym_gradient_solid(u(:,:,:,3),grad_sol(:,:,:,2:3)) ! dsuz, dzuz

  if (src_type(1)=='monopole') then
     grad_sol(:,:,:,2)=u(:,:,:,1)

  elseif (src_type(1)=='dipole') then 
     grad_sol(:,:,:,2)=real(2.,kind=realkind)*u(:,:,:,2)

  elseif (src_type(1)=='quadpole') then
     grad_sol(:,:,:,2)=u(:,:,:,1)-real(2.,kind=realkind)*u(:,:,:,2)

  endif

  grad_sol(:,:,:,2) = f_over_s_solid(grad_sol(:,:,:,2))

  ! Slightly dirty, but memory-cheaper: sum of 3 diag strain elements
  grad_sol(:,:,:,1) = grad_sol(:,:,:,1) + grad_sol(:,:,:,2) + grad_sol(:,:,:,3)

  call compute_recfile_cmb(velo,grad_sol(:,:,:,1))

end subroutine dump_velo_straintrace_cmb
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_strain(u, chi)
  !
  ! Compute the full, global strain tensor on-the-fly. Each of 6 (monopole: 4)
  ! components is stored separately for solid and fluid domains respectively.
  ! The dipole case is transfered to the (s,phi,z) system here.
  !
  ! Dumping Ekk, E11, E22, E13, -E23, and -E12, this has the advantage
  ! that if only lambda/bulk sound speed are of interest, then only a 
  ! scalar Ekk needs to be loaded. Be aware that diuj in the variable names does
  ! NOT stand for partial derivatives, but rather the ij component of the
  ! strain.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use data_pointwise,           ONLY: inv_rho_fluid, inv_s_rho_fluid, usz_fluid
  use data_source,              ONLY: src_type
  use pointwise_derivatives,    ONLY: axisym_gradient_fluid_add
  use pointwise_derivatives,    ONLY: axisym_gradient_fluid
  use pointwise_derivatives,    ONLY: axisym_gradient_solid_add
  use pointwise_derivatives,    ONLY: axisym_gradient_solid
  use pointwise_derivatives,    ONLY: f_over_s_solid
  use pointwise_derivatives,    ONLY: f_over_s_fluid
  !use wavefields_io,            ONLY: dump_half_f1_f2_over_s_fluid
  !use wavefields_io,            ONLY: dump_f1_f2_over_s_fluid
  use wavefields_io,            ONLY: dump_field_1d
  !use wavefields_io,            ONLY: dump_field_over_s_fluid_and_add
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in) :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(in) :: chi(0:npol,0:npol,nel_fluid)
  
  real(kind=realkind)             :: grad_sol(0:npol,0:npol,nel_solid,2)
  real(kind=realkind)             :: buff_solid(0:npol,0:npol,nel_solid)
  real(kind=realkind)             :: grad_flu(0:npol,0:npol,nel_fluid,2)
  character(len=5)                :: appisnap
  real(kind=realkind), parameter  :: two_rk = real(2, kind=realkind)
 
  call define_io_appendix(appisnap,istrain)
 
  ! SSSSSSS Solid region SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 
  ! s,z components, identical for all source types..........................
  if (src_type(1)=='dipole') then
     call axisym_gradient_solid(u(:,:,:,1) + u(:,:,:,2), grad_sol)
  else
     call axisym_gradient_solid(u(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus
  endif
 
  call dump_field_1d(grad_sol(:,:,:,1), '/strain_dsus_sol', appisnap, nel_solid) 
 
  call axisym_gradient_solid_add(u(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus
 
  ! calculate entire E31 term: (dsuz+dzus)/2
  grad_sol(:,:,:,1) = grad_sol(:,:,:,1) / two_rk
  call dump_field_1d(grad_sol(:,:,:,1), '/strain_dsuz_sol', appisnap, nel_solid)
 
  ! Components involving phi....................................................
  if (src_type(1) == 'monopole') then
     buff_solid = f_over_s_solid(u(:,:,:,1))
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !Epp
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
                        nel_solid) !Ekk
 
  elseif (src_type(1) == 'dipole') then 
     buff_solid = two_rk * f_over_s_solid(u(:,:,:,2))
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid)
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
                        nel_solid)
 
     call axisym_gradient_solid(u(:,:,:,1) - u(:,:,:,2), grad_sol) !1:dsup,2:dzup

     ! MvD: E12 and E23 have a minus sign not existent in the gradient, no idea
     !      whether that is on purpose!
     call dump_field_1d(f_over_s_solid(u(:,:,:,2)) + grad_sol(:,:,:,1) / two_rk, &
                        '/strain_dsup_sol', appisnap, nel_solid) !E12
 
     call dump_field_1d((f_over_s_solid(u(:,:,:,3)) +  grad_sol(:,:,:,2)) / two_rk, &
                        '/strain_dzup_sol', appisnap, nel_solid) !E23
 
  elseif (src_type(1) == 'quadpole') then
     buff_solid = f_over_s_solid(u(:,:,:,1) - two_rk * u(:,:,:,2))
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !Epp
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), & !Ekk
                        '/straintrace_sol', appisnap, nel_solid) 
  
     call axisym_gradient_solid(u(:,:,:,2), grad_sol) ! 1: dsup, 2: dzup
 
     ! MvD: E12 and E23 have a minus sign not existent in the gradient, no idea
     !      whether that is on purpose!
     call dump_field_1d(f_over_s_solid(u(:,:,:,1) - u(:,:,:,2) / two_rk) &
                            + grad_sol(:,:,:,1) / two_rk, &
                        '/strain_dsup_sol', appisnap, nel_solid) !E12
  
     call dump_field_1d(f_over_s_solid(u(:,:,:,3)) + grad_sol(:,:,:,2) / two_rk, &
                        '/strain_dzup_sol',appisnap, nel_solid) !E23
  endif
 
  ! FFFFFF Fluid region FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  !
  ! Fluid-region strain tensor is computed just like in the solid but for 
  ! displacement components ds(chi), dz(chi).
 
  if (have_fluid) then
 
     ! construct displacements in the fluid
     call axisym_gradient_fluid(chi, usz_fluid)
     usz_fluid(:,:,:,1) = usz_fluid(:,:,:,1) * inv_rho_fluid
     usz_fluid(:,:,:,2) = usz_fluid(:,:,:,2) * inv_rho_fluid
  
     ! gradient of s component
     call axisym_gradient_fluid(usz_fluid(:,:,:,1), grad_flu)   ! dsus, dzus
 
     ! save Ess
     call dump_field_1d(grad_flu(:,:,:,1), '/strain_dsus_flu', appisnap, nel_fluid) ! E11
 
     ! gradient of z component added to s-comp gradient for strain trace and E13
     call axisym_gradient_fluid_add(usz_fluid(:,:,:,2), grad_flu)  !1:dsuz+dzus, 2:dzuz+dsus
 
     ! calculate entire E31 term: (dsuz+dzus)/2
     grad_flu(:,:,:,1) = grad_flu(:,:,:,1) / two_rk
     call dump_field_1d(grad_flu(:,:,:,1), '/strain_dsuz_flu', appisnap, nel_fluid)
 
     ! Components involving phi................................................
 
     if (src_type(1) == 'monopole') then
        ! Calculate us/s and straintrace
        !call dump_field_over_s_fluid_and_add(usz_fluid(:,:,:,1), grad_flu(:,:,:,2), & !Epp,Eii
        !                                     '/strain_dpup_flu', '/straintrace_flu', &
        !                                     appisnap) 
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)), '/strain_dpup_flu', appisnap, &
                           nel_fluid) 
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)) + grad_flu(:,:,:,2), &
                           '/straintrace_flu', appisnap, nel_fluid) 
 
     elseif (src_type(1) == 'dipole') then
        !call dump_field_over_s_fluid_and_add(usz_fluid(:,:,:,1) - inv_s_rho_fluid * chi, &
        !                                     grad_flu(:,:,:,2), '/strain_dpup_flu', &
        !                                     '/straintrace_flu', appisnap)  !Epp,Eii
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1) - inv_s_rho_fluid * chi), &
                           '/strain_dpup_flu', appisnap, nel_fluid)  !Epp
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1) - inv_s_rho_fluid * chi) &
                            + grad_flu(:,:,:,2), &
                            '/straintrace_flu', appisnap, nel_fluid)  !Eii
 
        ! gradient of phi component
        call axisym_gradient_fluid(inv_s_rho_fluid * chi,grad_flu)   ! dsup, dzup
 
        !call dump_half_f1_f2_over_s_fluid(grad_flu(:,:,:,1), &
        !                                  usz_fluid(:,:,:,1) - inv_s_rho_fluid * chi, &
        !                                  '/strain_dsup_flu', appisnap)   ! E12
        call dump_field_1d((grad_flu(:,:,:,1) &
                            + f_over_s_fluid(usz_fluid(:,:,:,1) &
                                - inv_s_rho_fluid * chi)) / two_rk, &
                           '/strain_dsup_flu', appisnap, nel_fluid)   ! E12
 
        !call dump_half_f1_f2_over_s_fluid(grad_flu(:,:,:,2), usz_fluid(:,:,:,2), &
        !                                  '/strain_dzup_flu', appisnap)  ! E23
        call dump_field_1d((grad_flu(:,:,:,2) + f_over_s_fluid(usz_fluid(:,:,:,2))) / two_rk, &
                            '/strain_dzup_flu', appisnap, nel_fluid)  ! E23
 
     elseif (src_type(1) == 'quadpole') then
        !call dump_field_over_s_fluid_and_add(usz_fluid(:,:,:,1) - &
        !                                     two_rk * inv_s_rho_fluid * chi, &  !Epp,Eii
        !                                     grad_flu(:,:,:,2), & 
        !                                     '/strain_dpup_flu', '/straintrace_flu', &
        !                                     appisnap) 
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1) &
                                           - two_rk * inv_s_rho_fluid * chi), &  !Epp
                           '/strain_dpup_flu', appisnap, nel_fluid) 
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)&
                                           - two_rk * inv_s_rho_fluid * chi) &  !Eii
                            + grad_flu(:,:,:,2), & 
                           '/straintrace_flu', appisnap, nel_fluid) 
 
        ! gradient of phi component
        call axisym_gradient_fluid(inv_s_rho_fluid * chi,grad_flu)   ! dsup, dzup
 
        !call dump_f1_f2_over_s_fluid(grad_flu(:,:,:,1), &
        !                             usz_fluid(:,:,:,1) - inv_s_rho_fluid * chi, &
        !                             '/strain_dsup_flu', appisnap)   !E12
        call dump_field_1d((grad_flu(:,:,:,1) &
                             + f_over_s_fluid(usz_fluid(:,:,:,1) &
                                - inv_s_rho_fluid * chi)), &
                           '/strain_dsup_flu', appisnap, nel_fluid)   !E12
                                         
 
        !call dump_half_f1_f2_over_s_fluid(grad_flu(:,:,:,2), two_rk * usz_fluid(:,:,:,2), &
        !                                   '/strain_dzup_flu', appisnap)   !E23
        call dump_field_1d(grad_flu(:,:,:,2) / two_rk &
                            + f_over_s_fluid(usz_fluid(:,:,:,2)), &
                           '/strain_dzup_flu', appisnap, nel_fluid)   !E23
     endif   !src_type
 
  endif   ! have_fluid
 
end subroutine compute_strain
!=============================================================================


!-----------------------------------------------------------------------------
subroutine energy(disp1,vel,dchi1,ddchi)
  !
  ! Computes the kinetic and potential/elastic/stored energy in the solid and 
  ! fluid subdomains separately. This involves one additional evaluation of 
  ! the stiffness system (for the velocity vector rather than displacement) 
  ! in both domains, but additional cost is rather low if only performed every 
  ! 5th time step such as the default for the time being.
  ! Although non-applicable for accuracy measures of the numerical 
  ! approximation, the preserved total energy over time
  ! (only sufficiently after the source time function is active!) is a useful 
  ! check on the consistency of the used time scheme and also reveals whether 
  ! the system leaks, i.e. whether all boundary conditions are correct.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_source, ONLY: src_type
  use data_matr, ONLY : unassem_mass_rho_solid, unassem_mass_lam_fluid
  use stiffness
  use apply_masks
  use commun
  
  include 'mesh_params.h'
  
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3),intent(in) :: disp1
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3),intent(in) :: vel
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid),intent(in)   :: dchi1
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid),intent(in)   :: ddchi
  
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: dchi
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: stiff
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: stiff_flu
  real(kind=realkind) :: ekin_sol, epot_sol, ekin_flu, epot_flu

  ekin_sol = zero
  epot_sol = zero
  ekin_flu = zero
  epot_flu = zero
  stiff = zero
  stiff_flu=zero

  disp = disp1
  dchi = dchi1

  ! ssssss SOLID REGION sssssssssssssssssssssssssssssssssssssssssssssssssssss

  ! Compute potential/stored/elastic energy in solid region.
  select case (src_type(1))
     case ('monopole')
           call apply_axis_mask_onecomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_mono(stiff,disp)
           call apply_axis_mask_onecomp(stiff,nel_solid, ax_el_solid,naxel_solid)
     case ('dipole')
           call apply_axis_mask_twocomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_di(stiff,disp) 
           call apply_axis_mask_twocomp(stiff,nel_solid, ax_el_solid,naxel_solid)
     case ('quadpole')
           call apply_axis_mask_threecomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_quad(stiff,disp) 
           call apply_axis_mask_threecomp(stiff,nel_solid, ax_el_solid,naxel_solid)
  end select

  stiff = stiff * disp
  epot_sol = sum(stiff)
  epot_sol = two*pi*psum(epot_sol)

  ! Compute kinetic energy in solid region
  stiff(:,:,:,1) = vel(:,:,:,1)**2 * unassem_mass_rho_solid
  stiff(:,:,:,2) = vel(:,:,:,2)**2 * unassem_mass_rho_solid
  stiff(:,:,:,3) = vel(:,:,:,3)**2 * unassem_mass_rho_solid

  if (src_type(1)=='dipole') stiff(:,:,:,3)=two*stiff(:,:,:,3)

  ekin_sol = sum(stiff)
  ekin_sol = two * pi * psum(ekin_sol)

  ! fffff FLUID REGION ffffffffffffffffffffffffffffffffffffffffffffffffffffff
  if (have_fluid) then 

     ! Compute potential/stored/elastic energy in fluid region.
     ! Using the wave equation in the fluid here to avoid spatial derivatives.

     stiff_flu=ddchi**2 * unassem_mass_lam_fluid
     epot_flu=sum(stiff_flu)
     epot_flu=two*pi*psum(epot_flu)
     
     ! Fluid stiffness acting on derivative of potential
     stiff_flu = zero
     if (src_type(1) .ne. 'monopole') &
          call apply_axis_mask_scal(dchi,nel_fluid,ax_el_fluid,naxel_fluid)
     call glob_fluid_stiffness(stiff_flu,dchi)
     if (src_type(1) .ne. 'monopole') & 
         call apply_axis_mask_scal(stiff_flu,nel_fluid,ax_el_fluid,naxel_fluid)

     ! Compute kinetic energy in fluid region
     stiff_flu=stiff_flu*dchi
     ekin_flu=sum(stiff_flu)
     ekin_flu=two*pi*psum(ekin_flu)

  endif ! have_fluid

  ! Save into time series, only one processor needs to do this since global.
  if (lpr) then
     if (have_fluid) then 
        write(4444,8) t, ekin_sol, epot_sol
        write(4445,8) t, epot_flu, ekin_flu
     endif
     write(4446,9) t, epot_sol + epot_flu, ekin_sol + ekin_flu, &
                   half * (epot_sol + epot_flu + ekin_sol + ekin_flu)
  endif
8 format(3(1pe16.6))
9 format(4(1pe16.6))

end subroutine energy
!=============================================================================

!-----------------------------------------------------------------------------
subroutine bdry_copy2fluid(uflu,usol)
  !
  ! FLUID: solid-fluid boundary term (solid displ.) added to fluid stiffness
  !        minus sign in accordance with the definition of bdry_matr
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_matr, ONLY : bdry_matr,solflubdry_radius
  use data_source, ONLY : src_type
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(inout) :: uflu(0:npol,0:npol,1:nel_fluid)
  real(kind=realkind), intent(in)    :: usol(0:npol,0:npol,1:nel_solid,3)
  integer                            :: iel,jpols,jpolf,iels,ielf

  if (src_type(1) == 'dipole') then 
     
     do iel=1,nel_bdry
        jpols = bdry_jpol_solid(iel)
        jpolf = bdry_jpol_fluid(iel)
        iels = bdry_solid_el(iel)
        ielf = bdry_fluid_el(iel)

        uflu(0:npol,jpolf,ielf) = uflu(0:npol,jpolf,ielf) - &
             bdry_matr(0:npol,iel,1) * &
             (usol(0:npol,jpols,iels,1) + usol(0:npol,jpols,iels,2) ) - &
             bdry_matr(0:npol,iel,2) * usol(0:npol,jpols,iels,3)
     enddo

  else

     do iel=1,nel_bdry
        jpols = bdry_jpol_solid(iel)
        jpolf = bdry_jpol_fluid(iel)
        iels = bdry_solid_el(iel)
        ielf = bdry_fluid_el(iel)

        uflu(0:npol,jpolf,ielf) = uflu(0:npol,jpolf,ielf) - &
             bdry_matr(0:npol,iel,1) * usol(0:npol,jpols,iels,1) - &
             bdry_matr(0:npol,iel,2) * usol(0:npol,jpols,iels,3)
     enddo

  endif

end subroutine bdry_copy2fluid
!============================================================================

!----------------------------------------------------------------------------
subroutine bdry_copy2solid(usol,uflu)
  !
  ! SOLID: solid-fluid boundary term (fluid potential) added to solid stiffness
  !        plus sign in accordance with definition of bdry_matr
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  use data_matr, ONLY : bdry_matr,solflubdry_radius
  use data_source, ONLY : src_type
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(inout) :: usol(0:npol,0:npol,1:nel_solid,3)
  real(kind=realkind), intent(in)    :: uflu(0:npol,0:npol,1:nel_fluid)
  integer                            :: iel,jpols,jpolf,iels,ielf

  do iel = 1,nel_bdry
     jpols = bdry_jpol_solid(iel)
     jpolf = bdry_jpol_fluid(iel)
     iels = bdry_solid_el(iel)
     ielf = bdry_fluid_el(iel)  

     if (src_type(1) == 'dipole') then 

     usol(0:npol,jpols,iels,1) = usol(0:npol,jpols,iels,1) + &
          bdry_matr(0:npol,iel,1) * uflu(0:npol,jpolf,ielf)

     usol(0:npol,jpols,iels,2) = usol(0:npol,jpols,iels,2) + &
          bdry_matr(0:npol,iel,1) * uflu(0:npol,jpolf,ielf)

     else
     usol(0:npol,jpols,iels,1) = usol(0:npol,jpols,iels,1) + &
          bdry_matr(0:npol,iel,1) * uflu(0:npol,jpolf,ielf)        
     endif

     usol(0:npol,jpols,iels,3) = usol(0:npol,jpols,iels,3) + &
          bdry_matr(0:npol,iel,2) * uflu(0:npol,jpolf,ielf)

  enddo


end subroutine bdry_copy2solid
!=============================================================================


!============================
end module time_evol_wave
!============================
