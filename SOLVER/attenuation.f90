! Martin van Driel, 02/2013
! Martin@vanDriel.de
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!

module attenuation
    implicit none

    private
    public :: l2_error, q_linear_solid, fast_correct
    public :: invert_linear_solids

contains

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

    !use progressbar_mod

    double precision, intent(in)            :: Q, f_min, f_max
    integer, intent(in)                     :: N, nfsamp, max_it

    double precision, optional, intent(in)          :: Tw, Ty, d
    !f2py double precision, optional, intent(in)    :: Tw=.1, Ty=.1, d=.9998
    double precision                                :: Tw_loc = .1, Ty_loc = .1
    double precision                                :: d_loc = .9998

    logical, optional, intent(in)           :: fixfreq, verbose, exact
    !f2py logical, optional, intent(in)     :: fixfreq = 0, verbose = 0
    !f2py logical, optional, intent(in)     :: exact = 0
    logical                                 :: fixfreq_loc = .false., verbose_loc = .false.
    logical                                 :: exact_loc = .false.
    
    character(len=10), optional, intent(in) :: mode
    !f2py character(len=10), optional, intent(in) :: mode = 'maxwell'
    character(len=10)                       :: mode_loc = 'maxwell'

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

    !type(progressbar)   :: bar

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
        ! progress bar
        !call bar%printbar(100 * it / max_it)
        
        ! compute perturbed parameters (normal distributed? may gaussian be
        ! better?)
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

end module
