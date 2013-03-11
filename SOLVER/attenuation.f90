module attenuation
    use global_parameters,    only: realkind
    implicit none
    include 'mesh_params.h'

    private
    public :: prepare_attenuation
    public :: n_sls_attenuation
    public :: dump_memory_vars
    public :: time_step_memvars

    double precision, allocatable   :: y_j_attenuation(:)
    double precision, allocatable   :: w_j_attenuation(:), exp_w_j_deltat(:)
    double precision, allocatable   :: ts_fac_t(:), ts_fac_tm1(:)
    integer                         :: n_sls_attenuation
    logical                         :: do_corr_lowq, dump_memory_vars = .false.
    real(kind=realkind)             :: src_dev_tm1_glob(0:npol,0:npol,6,nel_solid)  
    real(kind=realkind)             :: src_tr_tm1_glob(0:npol,0:npol,nel_solid)  

contains

!-----------------------------------------------------------------------------------------
subroutine time_step_memvars(memvar, disp)
  !
  ! analytical time integration of memory variables (linear interpolation for
  ! the strain)
  ! MvD, attenutation notes, p 13.2
  !
  use data_time,            only: deltat
  use data_matr,            only: Q_mu, Q_kappa, mu_r, kappa_r
  include 'mesh_params.h'

  real(kind=realkind), intent(inout)    :: memvar(0:npol,0:npol,6,n_sls_attenuation,nel_solid)
  real(kind=realkind), intent(in)       :: disp(0:npol,0:npol,nel_solid,3)
  
  integer               :: iel, l, j, ipol, jpol
  double precision      :: yp_j_mu(n_sls_attenuation)
  double precision      :: yp_j_kappa(n_sls_attenuation)
  real(kind=realkind)   :: grad_t(0:npol,0:npol,nel_solid,6)
  real(kind=realkind)   :: trace_grad_t(0:npol,0:npol)
  real(kind=realkind)   :: trace_grad_tm1(0:npol,0:npol)
  real(kind=realkind)   :: src_tr_t(0:npol,0:npol)
  real(kind=realkind)   :: src_tr_tm1(0:npol,0:npol)
  real(kind=realkind)   :: src_dev_t(0:npol,0:npol,6)
  real(kind=realkind)   :: src_dev_tm1(0:npol,0:npol,6)

  ! compute global strain of current time step
  call compute_strain(disp, grad_t)

  do iel=1, nel_solid
     ! compute local coefficients y_j for kappa and mu
     if (do_corr_lowq) then
        call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
     else
        yp_j_mu = y_j_attenuation / Q_mu(iel)
        yp_j_kappa = y_j_attenuation / Q_kappa(iel)
     endif

     trace_grad_t(:,:) = sum(grad_t(:,:,iel,1:3), dim=3)

     ! analytical time stepping, monople/isotropic hardcoded

     ! compute new source terms (excluding the weighting)
     src_tr_t(:,:) = kappa_r(:,:,iel) * trace_grad_t(:,:)
     src_dev_t(:,:,1) = mu_r(:,:,iel) * 2 * (grad_t(:,:,iel,1)  - trace_grad_t(:,:) / 3)
     src_dev_t(:,:,2) = mu_r(:,:,iel) * 2 * (grad_t(:,:,iel,2)  - trace_grad_t(:,:) / 3)
     src_dev_t(:,:,3) = mu_r(:,:,iel) * 2 * (grad_t(:,:,iel,3)  - trace_grad_t(:,:) / 3)
     src_dev_t(:,:,5) = mu_r(:,:,iel) * grad_t(:,:,iel,5)
     
     ! load old source terms
     src_tr_tm1(:,:) = src_tr_tm1_glob(:,:,iel)
     src_dev_tm1(:,:,:) = src_dev_tm1_glob(:,:,:,iel)
     
     do j=1, n_sls_attenuation
        ! do the timestep
        memvar(:,:,:,j,iel) = memvar(:,:,:,j,iel) * exp_w_j_deltat(j)

        memvar(:,:,1,j,iel) = memvar(:,:,1,j,iel) &
            + ts_fac_t(j) * (yp_j_kappa(j) * src_tr_t(:,:) &
                             +  yp_j_mu(j) * src_dev_t(:,:,1)) &
            + ts_fac_tm1(j) * (yp_j_kappa(j) * src_tr_tm1(:,:) &
                               +  yp_j_mu(j) * src_dev_tm1(:,:,1))

        memvar(:,:,2,j,iel) = memvar(:,:,2,j,iel) &
            + ts_fac_t(j) * (yp_j_kappa(j) * src_tr_t(:,:) &
                             +  yp_j_mu(j) * src_dev_t(:,:,2)) &
            + ts_fac_tm1(j) * (yp_j_kappa(j) * src_tr_tm1(:,:) &
                               +  yp_j_mu(j) * src_dev_tm1(:,:,2))

        memvar(:,:,3,j,iel) = memvar(:,:,3,j,iel) &
            + ts_fac_t(j) * (yp_j_kappa(j) * src_tr_t(:,:) &
                             +  yp_j_mu(j) * src_dev_t(:,:,3)) &
            + ts_fac_tm1(j) * (yp_j_kappa(j) * src_tr_tm1(:,:) &
                               +  yp_j_mu(j) * src_dev_tm1(:,:,3))

        memvar(:,:,5,j,iel) = memvar(:,:,5,j,iel) &
            + ts_fac_t(j) * yp_j_mu(j) * src_dev_t(:,:,5) &
            + ts_fac_tm1(j) * yp_j_mu(j) * src_dev_tm1(:,:,5)

     enddo
     ! save srcs for next iteration
     src_dev_tm1_glob(:,:,:,iel) = src_dev_t(:,:,:)
  enddo
  
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_strain(u, grad_u)
  !
  ! compute strain in Voigt notation 
  ! (i.e. E1 = E11, E2 = E22, E3 = E33, E4 = 2E23, E5 = 2E31, E6 = 2E12)
  !  
  use data_source,              ONLY: src_type
  use pointwise_derivatives,    ONLY: axisym_gradient_solid_add
  use pointwise_derivatives,    ONLY: axisym_gradient_solid
  
  include 'mesh_params.h'
  
  real(kind=realkind), intent(in)   :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(out)  :: grad_u(0:npol,0:npol,nel_solid,6)
  
  real(kind=realkind)             :: grad_buff(0:npol,0:npol,nel_solid,2)
  
  grad_u(:,:,:,:) = 0
  
  ! s,z components, identical for all source types..........................
  if (src_type(1)=='dipole') then
     call axisym_gradient_solid(u(:,:,:,1) + u(:,:,:,2), grad_buff)
  else
     call axisym_gradient_solid(u(:,:,:,1), grad_buff) ! 1: dsus, 2: dzus
  endif

  grad_u(:,:,:,1) = grad_buff(:,:,:,1)  ! dsus

  grad_buff(:,:,:,1) = 0
 
  call axisym_gradient_solid_add(u(:,:,:,3), grad_buff) ! 1:dsuz+dzus, 2:dzuz

  grad_u(:,:,:,5) = grad_buff(:,:,:,1)  ! dsuz + dzus (incl faktor of 2 from voigt notation)
  grad_u(:,:,:,3) = grad_buff(:,:,:,2)  ! dzuz
 
 
  ! Components involving phi....................................................
  ! hardcode monopole for a start

  call field_over_s_solid(u(:,:,:,1), grad_u(:,:,:,2)) ! us / s

 
end subroutine compute_strain
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine field_over_s_solid(f, f_over_s)
  !
  ! computes f/s using L'Hospital's rule lim f/s = lim df/ds at the axis (s = 0)
  !
  use data_pointwise,           ONLY: inv_s_solid
  use pointwise_derivatives,    ONLY: dsdf_solid_allaxis
  use data_mesh,                ONLY: naxel_solid, ax_el_solid
  
  real(kind=realkind),intent(in)      :: f(0:npol,0:npol,nel_solid)
  real(kind=realkind),intent(out)     :: f_over_s(0:npol,0:npol,nel_solid)
  real(kind=realkind)                 :: dsdf(0:npol,naxel_solid)
  integer                             :: iel
  
  f_over_s = f

  call dsdf_solid_allaxis(f_over_s, dsdf) ! axial f/s
  do iel=1, naxel_solid
     inv_s_solid(0,:,ax_el_solid(iel)) = dsdf(:,iel)
     f_over_s(0,:,ax_el_solid(iel)) = 1 ! otherwise this  would result in df/ds * f below
  enddo

  ! construct sum of f/s and g (e.g. straintrace)
  f_over_s = inv_s_solid * f_over_s

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_attenuation(lambda, mu)
  !
  ! read 'inparam_attenuation' file and compute precomputable terms
  !
  use data_io,              only: infopath, lfinfo
  use data_proc,            only: lpr
  use data_time,            only: deltat
  use global_parameters,    only: pi
  use data_matr,            only: Q_mu, Q_kappa, mu_r, kappa_r
  use data_mesh_preloop,    only: ielsolid
  

  include 'mesh_params.h'

  double precision, intent(in)      ::lambda(0:npol,0:npol,1:nelem), mu(0:npol,0:npol,1:nelem)

  double precision                  :: f_min, f_max
  integer                           :: nfsamp, max_it, i, iel
  double precision                  :: Tw, Ty, d
  logical                           :: fixfreq
  double precision, allocatable     :: w_samp(:), q_fit(:), chil(:)
  double precision                  :: yp_j_mu(n_sls_attenuation)
  double precision                  :: yp_j_kappa(n_sls_attenuation)

  if (lpr) print *, '  ...reading inparam_attanuation...'
  open(unit=164, file='inparam_attenuation')

  read(164,*) n_sls_attenuation
  read(164,*) f_min
  read(164,*) f_max
  read(164,*) do_corr_lowq
  
  read(164,*) nfsamp
  read(164,*) max_it
  read(164,*) Tw
  read(164,*) Ty
  read(164,*) d
  read(164,*) fixfreq
  read(164,*) dump_memory_vars

  close(unit=164)
  
  allocate(w_samp(nfsamp))
  allocate(q_fit(nfsamp))
  allocate(chil(max_it))
  
  allocate(w_j_attenuation(n_sls_attenuation))
  allocate(exp_w_j_deltat(n_sls_attenuation))
  allocate(y_j_attenuation(n_sls_attenuation))
  
  allocate(ts_fac_t(n_sls_attenuation))
  allocate(ts_fac_tm1(n_sls_attenuation))
  
  if (lpr) print *, '  ...inverting for standard linear solid parameters...'

  call invert_linear_solids(1.d0, f_min, f_max, n_sls_attenuation, nfsamp, max_it, Tw, &
                            Ty, d, fixfreq, .false., .false., 'maxwell', w_j_attenuation, &
                            y_j_attenuation, w_samp, q_fit, chil)
  
  ! prefactors for the exact time stepping (att nodes p 13.3)
  exp_w_j_deltat = dexp(-w_j_attenuation * deltat)
  ts_fac_tm1 = ((1 - exp_w_j_deltat) / (w_j_attenuation * deltat) - exp_w_j_deltat)
  ts_fac_t = ((exp_w_j_deltat - 1) / (w_j_attenuation * deltat) + 1)

  if (lpr) then
      print *, '  ...log-l2 misfit: ', chil(max_it)
      print *, '  ...frequencies  : ', w_j_attenuation / (2. * pi)
      print *, '  ...exp-frequencies  : ', exp_w_j_deltat

      print *, '  ...writing fitted Q to file...'
      open(unit=165, file=infopath(1:lfinfo)//'/attenuation_q_fitted', status='new')
      write(165,*) (w_samp(i), q_fit(i), char(10), i=1,nfsamp)
      close(unit=165)
      
      print *, '  ...writing convergence of chi to file...'
      open(unit=166, file=infopath(1:lfinfo)//'/attenuation_convergence', status='new')
      write(166,*) (chil(i), char(10), i=1,max_it)
      close(unit=166)
  endif


  if (lpr) print *, '  ...calculating relaxed moduli...'

  allocate(mu_r(0:npol,0:npol,nel_solid))
  allocate(kappa_r(0:npol,0:npol,nel_solid))

  do iel=1, nel_solid
     if (do_corr_lowq) then
        call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
     else
        yp_j_mu = y_j_attenuation / Q_mu(iel)
        yp_j_kappa = y_j_attenuation / Q_kappa(iel)
     endif
     mu_r(:,:,iel) =  mu(:,:,ielsolid(iel)) / (1.d0 + sum(yp_j_mu))
     kappa_r(:,:,iel) =  (lambda(:,:,ielsolid(iel)) &
                            + 2.d0 / 3.d0 * mu(:,:,ielsolid(iel))) &
                            / (1.d0 + sum(yp_j_kappa))
  enddo


end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine q_linear_solid(y_j, w_j, w, exact, Qls)
    !
    ! compute Q after (Emmerich & Korn, inverse of eq 21)
    ! linearized version (exact = false) is eq 22 in E&K
    !
    double precision, intent(in)    :: y_j(:), w_j(:), w(:)
    double precision, intent(out)   :: Qls(size(w))
    integer                         :: j
    
    logical, optional, intent(in)           :: exact
    !f2py logical, optional, intent(in)     :: exact = 0 
    logical                                 :: exact_loc = .false.
    
    double precision                :: Qls_denom(size(w))
    
    if (present(exact)) exact_loc = exact
    
    Qls = 1
    if (exact_loc) then
        do j=1, size(y_j)
            Qls = Qls + y_j(j) *  w**2 / (w**2 + w_j(j)**2)
        enddo
    endif
    
    Qls_denom = 0
    do j=1, size(y_j)
        Qls_denom = Qls_denom + y_j(j) * w * w_j(j) / (w**2 + w_j(j)**2)
    enddo

    Qls = Qls / Qls_denom
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fast_correct(y_j, yp_j)
    !
    ! computes a first order correction to the linearized coefficients:
    ! yp_j_corrected = y_j * delta_j 
    !
    ! MvD Attenuation Notes, p. 17.3 bottom
    !
    double precision, intent(in)    :: y_j(:)
    double precision, intent(out)   :: yp_j(size(y_j))
    
    double precision                :: dy_j(size(y_j))
    integer                         :: k

    dy_j(1) = 1 + .5 * y_j(1)

    do k=2, size(y_j)
        dy_j(k) = dy_j(k-1) + (dy_j(k-1) - .5) * y_j(k-1) + .5 * y_j(k)
    enddo

    yp_j = y_j * dy_j
    
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine l2_error(Q, Qls, lognorm, lse)
    !
    ! returns l2 misfit between constant Q and fitted Q using standard linear solids
    !
    double precision, intent(in)    :: Q, Qls(:)
    
    logical, optional, intent(in)   :: lognorm
    ! optional argument with default value (a bit nasty in f2py)
    !f2py logical, optional, intent(in) :: lognorm = 1
    logical :: lognorm_loc = .true.
    
    double precision, intent(out)   :: lse
    integer                         :: nfsamp, i
    
    if (present(lognorm)) lognorm_loc = lognorm

    lse = 0
    nfsamp = size(Qls)

    if (lognorm_loc) then
        !print *, 'log-l2 norm'
        do i=1, nfsamp
            lse = lse + (log(Q / Qls(i)))**2
        end do
    else
        !print *, 'standard l2 norm'
        do i=1, nfsamp
            lse = lse + (1/Q - 1/Qls(i))**2
        end do
        lse = lse * Q**2
    endif
    lse = lse / float(nfsamp)
    lse = dsqrt(lse)
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine invert_linear_solids(Q, f_min, f_max, N, nfsamp, max_it, Tw, Ty, d, &
                                 fixfreq, verbose, exact, mode, w_j, y_j, w, q_fit, chil)
    !
    ! Inverts for constant Q, minimizing the L2 error for 1/Q using a simulated annealing
    ! approach (varying peak frequencies and amplitudes).

    ! Parameters:
    ! Q:              clear
    ! f_min, fmax:    frequency band (in Hz)
    ! N:              number of standard linear solids
    ! nfsamp:         number of sampling frequencies for computation of the misfit (log
    !                   spaced in freqeuncy band)

    ! max_it:         number of iterations
    ! Tw:             starting temperature for the frequencies
    ! Ty:             starting temperature for the amplitudes
    ! d:              temperature decay
    ! fixfreq:        use log spaced peak frequencies (fixed)
    ! verbose:        clear
    ! exact:          use exact relation for Q and the coefficients (Emmerich & Korn, eq
    !                   21). If false, linearized version is used (large Q approximation, eq
    !                   22).
    ! mode:           'maxwell' (default) oder 'zener' depending on the reology
    !
    ! Returns:
    ! w_j:            relaxation frequencies, equals 1/tau_sigma in zener
    !                   formulation
    ! y_j:            coefficients of the linear solids, (Emmerich & Korn, eq 23 and 24)
    !                   if mode is set to 'zener', this array contains the
    !                   tau_epsilon as defined by Blanch et al, eq 12.
    ! w:              sampling frequencies at which Q(w) is minimized
    ! q_fit:          resulting q(w) at these frequencies
    ! chil:           error as a function of iteration to check convergence,
    !                   Note that this version uses log-l2 norm!
    !

    double precision, intent(in)            :: Q, f_min, f_max
    integer, intent(in)                     :: N, nfsamp, max_it

    double precision, optional, intent(in)          :: Tw, Ty, d
    !f2py double precision, optional, intent(in)    :: Tw=.1, Ty=.1, d=.99995
    double precision                                :: Tw_loc = .1, Ty_loc = .1
    double precision                                :: d_loc = .99995

    logical, optional, intent(in)           :: fixfreq, verbose, exact
    !f2py logical, optional, intent(in)     :: fixfreq = 0, verbose = 0
    !f2py logical, optional, intent(in)     :: exact = 0
    logical                                 :: fixfreq_loc = .false., verbose_loc = .false.
    logical                                 :: exact_loc = .false.
    
    character(len=7), optional, intent(in)  :: mode
    !f2py character(len=7), optional, intent(in) :: mode = 'maxwell'
    character(len=7)                        :: mode_loc = 'maxwell'

    double precision, intent(out)   :: w_j(N)
    double precision, intent(out)   :: y_j(N)
    double precision, intent(out)   :: w(nfsamp) 
    double precision, intent(out)   :: q_fit(nfsamp) 
    double precision, intent(out)   :: chil(max_it) 

    double precision                :: w_j_test(N)
    double precision                :: y_j_test(N)
    double precision                :: expo
    double precision                :: chi, chi_test

    integer             :: j, it, last_it_print

    ! set default values
    if (present(Tw)) Tw_loc = Tw
    if (present(Ty)) Ty_loc = Ty
    if (present(d)) d_loc = d

    if (present(fixfreq)) fixfreq_loc = fixfreq
    if (present(verbose)) verbose_loc = verbose
    if (present(exact)) exact_loc = exact
    
    if (present(mode)) mode_loc = mode
    
    if ((mode_loc .ne. 'maxwell') .and. (mode_loc .ne. 'zener')) then
        print *, "ERROR: mode should be either 'maxwell' or 'zener'"
        return
    endif


    ! Set the starting test frequencies equally spaced in log frequency
    expo = (log10(f_max) - log10(f_min)) / (N - 1.d0)
    do j=1, N
        ! pi = 4 * atan(1)
        w_j_test(j) = datan(1.d0) * 8.d0 * 10**(log10(f_min) + (j - 1) * expo)
    end do

    if (verbose_loc) print *, w_j_test
    
    ! Set the sampling frequencies equally spaced in log frequency
    expo = (log10(f_max) - log10(f_min)) / (nfsamp - 1.d0)
    do j=1, nfsamp
        w(j) = datan(1.d0) * 8.d0 * 10**(log10(f_min) + (j - 1) * expo)
    end do

    if (verbose_loc) print *, w
    
    ! initial weights
    y_j_test = 1.d0 / Q * 1.5
    if (verbose_loc) print *, y_j_test

    ! initial Q(omega)
    call q_linear_solid(y_j=y_j_test, w_j=w_j_test, w=w, exact=exact_loc, Qls=q_fit)
    
    if (verbose_loc) print *, q_fit
   
    ! initial chi
    call l2_error(Q=Q, Qls=q_fit, lognorm=.true., lse=chi)
    if (verbose_loc) print *, 'initital chi: ', chi

    y_j(:) = y_j_test(:)
    w_j(:) = w_j_test(:)
    
    last_it_print = -1
    do it=1, max_it
        do j=1, N
            if (.not. fixfreq_loc) &
                w_j_test(j) = w_j(j) * (1.0 + (0.5 - rand()) * Tw_loc)
            y_j_test(j) = y_j(j) * (1.0 + (0.5 - rand()) * Ty_loc)
        enddo
    
        ! compute Q with test parameters
        call q_linear_solid(y_j=y_j_test, w_j=w_j_test, w=w, exact=exact_loc, Qls=q_fit)
        
        ! compute new misfit and new temperature
        call l2_error(Q=Q, Qls=q_fit, lognorm=.true., lse=chi_test)
        Tw_loc = Tw_loc * d_loc
        Ty_loc = Ty_loc * d_loc
                                        
        ! check if the tested parameters are better
        if (chi_test < chi) then
            y_j(:) = y_j_test(:)
            w_j(:) = w_j_test(:)
            chi = chi_test

            if (verbose_loc) then
                print *, '---------------'
                print *, it, chi
                print *, w_j / (8 * tan(1.))
                print *, y_j
            endif
    
        endif
        chil(it) = chi
    enddo

    if (mode_loc .eq. 'zener') then
        ! compare Attenuation Notes, p 18.1
        ! easy to find from Blanch et al, eq. (12) and Emmerick & Korn, eq. (21)
        y_j = (y_j + 1) / w_j
    endif

end subroutine
!-----------------------------------------------------------------------------------------

end module
