!========================
module lateral_heterogeneities
!========================

  use global_parameters
  use data_heterogeneous
  use data_io
  use data_proc
  use utlity, only :  compute_coordinates
  use data_source, only : rot_src
  
  implicit none
  
  public :: compute_heterogeneities
  public :: write_VTK_bin_scal_pts
  private

contains

!----------------------------------------------------------------------------------------
subroutine compute_heterogeneities(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                   fa_ani_theta, fa_ani_phi, ieldom)
    include 'mesh_params.h'

    integer :: ij
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
    double precision, dimension(0:npol,0:npol,nelem), intent(inout), optional :: &
           xi_ani, phi_ani, eta_ani, fa_ani_theta, fa_ani_phi
    integer, dimension(nelem), intent(in), optional :: ieldom

    double precision, dimension(0:npol,0:npol,nelem) :: rhopost,lambdapost,mupost
    double precision, dimension(:,:,:), allocatable :: &
           xi_ani_post, phi_ani_post, eta_ani_post, fa_ani_theta_post, fa_ani_phi_post

    if (lpr) then
       write(6,*)
       write(6,*) '   +++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*) '   ++++++++    Lateral Heterogeneities  ++++++++'
       write(6,*) '   +++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*)
       write(6,*)
       write(6,*) 'read parameter file for heterogeneities: inparam_hetero'
       write(6,*)
       write(6,*) ' !!!!!!!!! W A R N I N G !!!!!!!! '
       write(6,*) 'These lateral additions have not been thoroughly tested yet!'
       write(6,*)
    endif

    if (present(xi_ani) .and. present(phi_ani) .and.  present(eta_ani) &
          .and. present(fa_ani_theta) .and. present(fa_ani_phi)) then
       allocate(xi_ani_post(0:npol,0:npol,nelem))
       allocate(phi_ani_post(0:npol,0:npol,nelem))
       allocate(eta_ani_post(0:npol,0:npol,nelem))
       allocate(fa_ani_theta_post(0:npol,0:npol,nelem))
       allocate(fa_ani_phi_post(0:npol,0:npol,nelem))
       ani_hetero = .true.
    else
       ani_hetero = .false.
    endif
    
    mupost = mu
    lambdapost = lambda
    rhopost = rho
    if (ani_hetero) then
       xi_ani_post = xi_ani
       phi_ani_post = phi_ani
       eta_ani_post = eta_ani
       fa_ani_theta_post = fa_ani_theta
       fa_ani_phi_post = fa_ani_phi
    endif


    call read_param_hetero

    ! load and add heterogeneities
    ! loop over heterogeneities and call functions for each separately!
 
    do ij = 1, num_het
       if (het_format(ij)=='funct') then
          ! functional perturbations (sinus, gauss, trian, gss1d, inclp, inclr, const)
          call load_het_funct(rho,lambda,mu,rhopost,lambdapost,mupost,ij) 
       elseif (het_format(ij) == 'discr') then
          ! interpolate discrete model of arbitrary locations/perturbations
          if (.not. ani_hetero) then
             call load_het_discr(rho,lambda,mu,rhopost,lambdapost,mupost,ij) 
          else
             call load_het_discr(rho, lambda, mu, rhopost, lambdapost, mupost, ij, &
                           xi_ani, phi_ani, eta_ani, xi_ani_post, phi_ani_post, &
                           eta_ani_post, fa_ani_theta_post, fa_ani_phi_post)
          endif
       elseif (het_format(ij) == 'rndm') then 
          ! add random fluctuations to radial model
          call load_random(rho,lambda,mu,rhopost,lambdapost,mupost,ij) 
       elseif (het_format(ij) == 'ica') then 
          ! add inner core anisotropy
          if (.not. ani_hetero) then
             write(6,*) ''
             write(6,*) 'ERROR: inner core anisotropy need an anisotropic model -'
             write(6,*) '   either choose an anisotropic background model or'
             write(6,*) '   activate force anisotropy in inparam!'
             stop
          endif
          call load_ica(rho, lambda, mu, lambdapost, xi_ani_post, phi_ani_post, &
                        eta_ani_post, fa_ani_theta_post, fa_ani_phi_post, ij, ieldom)
                    
          deallocate(fa_theta_ica, fa_phi_ica)
          deallocate(a_ica, b_ica, c_ica)
       endif

       if (add_up) then
          mu = mupost
          lambda = lambdapost
          rho = rhopost
          if (ani_hetero) then
             xi_ani = xi_ani_post
             phi_ani = phi_ani_post
             eta_ani = eta_ani_post
             fa_ani_theta = fa_ani_theta_post
             fa_ani_phi = fa_ani_phi_post
          endif
       endif
    enddo
 
    if (.not. add_up) then
       mu = mupost
       lambda = lambdapost
       rho = rhopost
       if (ani_hetero) then
          xi_ani = xi_ani_post
          phi_ani = phi_ani_post
          eta_ani = eta_ani_post
          fa_ani_theta = fa_ani_theta_post
          fa_ani_phi = fa_ani_phi_post
       endif
    endif
 
    write(6,*)'final model done, now vtk files...'

    if (ani_hetero) then
       call plot_hetero_region_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                   fa_ani_theta, fa_ani_phi)
    else
       call plot_hetero_region_vtk(rho, lambda, mu)
    endif

    deallocate(het_format, het_file_discr, het_funct_type, rdep, grad, gradrdep1, &
               gradrdep2, r_het1, r_het2, th_het1, th_het2, delta_rho, delta_vp, &
               delta_vs, inverseshape, p_inv_dist, R_inv_dist, het_ani_discr, &
               het_rel_discr) 


    if (lpr) then
       write(6,*)
       write(6,*) '   ++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*) '   ++++++++ done with Lateral Heterogeneities  ++++++++'
       write(6,*) '   ++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*)
    endif

end subroutine compute_heterogeneities
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine read_param_hetero

    integer :: ij, i
    character(len=100) :: junk

    inquire(file="inparam_hetero", EXIST=file_exists)

    if (.not. file_exists) then 
       write(6,*) ''
       write(6,*) 'ERROR: lateral heterogeneity set in inparam, but'
       write(6,*) '       inparam_hetro does not exist!'
       stop
    endif

    write(6,*) 'starting read_param_hetero'
    open(unit=91, file='inparam_hetero')

    read(91,*) num_het
    if (lpr) write(6,*) 'adding ', num_het, ' regions...'

    allocate(het_format(num_het), het_file_discr(num_het), &
             het_funct_type(num_het), rdep(num_het), grad(num_het), &
             gradrdep1(num_het), gradrdep2(num_het), r_het1(num_het), &
             r_het2(num_het), th_het1(num_het), th_het2(num_het), &
             delta_rho(num_het), delta_vp(num_het), delta_vs(num_het), &
             inverseshape(num_het), p_inv_dist(num_het), R_inv_dist(num_het), &
             het_ani_discr(num_het), het_rel_discr(num_het))
    
    read(91,*) add_up

    do ij = 1, num_het
       read(91,*) junk
       read(91,*) het_format(ij)

       if (het_format(ij) == 'discr') then
          read(91,*) het_file_discr(ij)
          read(91,*) het_ani_discr(ij)
          read(91,*) het_rel_discr(ij)
          read(91,*) p_inv_dist(ij)
          read(91,*) R_inv_dist(ij)

       elseif (het_format(ij) == 'rndm') then
          read(91,*) het_funct_type(ij)
          read(91,*) r_het1(ij), r_het2(ij)
          read(91,*) th_het1(ij), th_het2(ij)
          read(91,*) delta_rho(ij)
          read(91,*) delta_vp(ij)
          read(91,*) delta_vs(ij)

          inverseshape(ij) = index(het_funct_type(ij),'_i')-1
          
          if (inverseshape(ij) > 0) then
              het_funct_type(ij) = het_funct_type(ij)(1:inverseshape(ij))
          endif

       elseif (het_format(ij) == 'funct') then
          read(91,*) het_funct_type(ij)
          read(91,*) rdep(ij)
          read(91,*) grad(ij)
          read(91,*) gradrdep1(ij),gradrdep2(ij)
          read(91,*) r_het1(ij),r_het2(ij)
          read(91,*) th_het1(ij),th_het2(ij)
          read(91,*) delta_rho(ij)
          read(91,*) delta_vp(ij)
          read(91,*) delta_vs(ij)

       elseif (het_format(ij) == 'ica') then
          read(91,*) num_slices
          
          allocate(fa_theta_ica(num_slices), fa_phi_ica(num_slices))
          allocate(a_ica(num_slices), b_ica(num_slices), c_ica(num_slices))
          allocate(theta_slices(num_slices + 1))
          
          read(91,*) theta_slices
             theta_slices = theta_slices / 180. * pi

          do i = 1, num_slices
             read(91,*) fa_theta_ica(i), fa_phi_ica(i)
             fa_theta_ica(i) = fa_theta_ica(i) / 180. * pi
             fa_phi_ica(i) = fa_phi_ica(i) / 180. * pi
             read(91,*) a_ica(i), b_ica(i), c_ica(i)
          enddo

       else
          write(6,*)'Unknown heterogeneity input type: ', het_format(ij)
          stop
       endif
     
    enddo

    ! degree to radians
    th_het1 = th_het1 / 180. * pi
    th_het2 = th_het2 / 180. * pi

    ! percent to decimal
    delta_rho = delta_rho / 100.
    delta_vp = delta_vp / 100.
    delta_vs = delta_vs / 100.
            
    R_inv_dist = R_inv_dist * 1000.

    if (lpr) then 
       do ij=1, num_het
          write(6,*) ''
          write(6,*) 'Lateral Heterogeneity No. ', ij, 'of type ', het_format(ij)
          write(6,*) ''
          if (het_format(ij)=='funct' .or. het_format(ij)=='rndm') then
             write(6,*) ''
             write(6,*) '   Radius (lower/upper bound) [km]:', &
                         r_het1(ij) / 1000., r_het2(ij) / 1000.
             write(6,*) '   Colatitude (lower/upper bound) [deg]:', &
                         th_het1(ij) * 180. / pi, th_het2(ij) * 180. / pi
             write(6,*) '   delta rho [%]:', delta_rho(ij)
             write(6,*) '   delta vp  [%]:', delta_vp(ij)
             write(6,*) '   delta vs  [%]:', delta_vs(ij)
             write(6,*) ''
          endif
       enddo
    endif

    ! need to rotate coordinates if source is not along axis (beneath the north pole)
    if (rot_src .and. (het_format(ij) == 'rndm' .or. het_format(ij) == 'funct')) then 
       write(6,*) 'need to rotate the heterogeneous domain with the source....'

       do i=1, num_het
          write(6,*)'Before rotation r th ph 1:', &
             r_het1(i), th_het1(i) * 180. / pi
          write(6,*)'Before rotation r th ph 2:', &
             r_het2(i), th_het2(i) * 180. / pi

          call rotate_hetero(r_het1(i), th_het1(i))
          call rotate_hetero(r_het2(i), th_het2(i))

          write(6,*)'After rotation r th ph 1:', &
             r_het1(i), th_het1(i) * 180. / pi
          write(6,*)'After rotation r th ph 2:', &
             r_het2(i), th_het2(i) * 180. / pi
       enddo
    endif

    rhetmin = 1.d10
    rhetmax = 0.d0
    thhetmin = pi
    thhetmax = 0.d0

    write(6,*)'done with read_param_hetero'

end subroutine read_param_hetero
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine rotate_hetero(r,th)

    double precision,intent(inout) :: r,th
    double precision :: x_vec(3), x_vec_rot(3), r_r, arg1

    x_vec(1) = r * dsin(th)
    x_vec(2) = 0.d0
    x_vec(3) = r * dcos(th) 
    
    x_vec_rot = matmul(trans_rot_mat,x_vec)
    
    write(23,*) x_vec
    write(23,*) x_vec_rot
 
    r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
    th = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )

end subroutine rotate_hetero
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine load_ica(rho, lambda, mu, lambdapost, xi_ani_post, phi_ani_post, & 
                    eta_ani_post, fa_ani_theta_post, fa_ani_phi_post, hetind, ieldom)

    use utlity,     only: thetacoord, rcoord
    use data_mesh,  only: discont

    double precision, dimension(0:npol,0:npol,nelem), intent(in) :: rho, lambda, mu
    double precision, dimension(0:npol,0:npol,nelem), intent(out) :: lambdapost, &
           xi_ani_post, phi_ani_post, eta_ani_post, fa_ani_theta_post, fa_ani_phi_post
    integer, dimension(nelem), intent(in) :: ieldom

    double precision, dimension(1:3) :: fast_axis_np
    integer :: hetind
    double precision, allocatable :: fast_axis_src(:,:)
    double precision :: vptmp, vstmp, arg1
    integer :: iel, ipol, jpol, i
    
    allocate(fast_axis_src(num_slices,1:3))

    if (lpr) then
        write(6,*) ''
        write(6,*) 'Adding Inner Core Anisotropy !!!'
    endif
    
    do i=1, num_slices
        fast_axis_np(1) = sin(fa_theta_ica(i)) * cos(fa_phi_ica(i))
        fast_axis_np(2) = sin(fa_theta_ica(i)) * sin(fa_phi_ica(i))
        fast_axis_np(3) = cos(fa_theta_ica(i))
        !if (rot_src) then 
        !    fast_axis_src(i,:) = matmul(transpose(rot_mat), fast_axis_np)
        !else
            fast_axis_src(i,:) = fast_axis_np
        !endif
        if (lpr) then
            write(6,*) i
            write(6,*) '  Fast Axis        :', fast_axis_np
            write(6,*) '  Fast Axis rotated:', fast_axis_src(i,:)
            write(6,*) '  ROTATION of fastaxis disabled for now - does not make'
            write(6,*) '  sense if not rotating in all three euler angles!'
            write(6,*) ''
        endif
    enddo
          
  
    ! compute theta and phi of the fast axis (phi is not well defined
    ! at the northpole)
    do iel=1, nelem
       if (ieldom(iel) == ndisc) then  
          do ipol=0, npol
             do jpol=0, npol
                
                vptmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                             rho(ipol,jpol,iel))
                vstmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel))

                do i=1, num_slices
                    if (thetacoord(ipol, jpol, iel) <= theta_slices(i + 1) &
                            .and. thetacoord(ipol, jpol, iel) >= theta_slices(i)) then

                        fa_ani_theta_post(ipol,jpol,iel) = acos(fast_axis_src(i,3))
  
                        arg1 = (fast_axis_src(i,1) + smallval_dble) / &
                               ((fast_axis_src(i,1)**2 + fast_axis_src(i,2)**2)**.5 &
                                + smallval_dble)
                        
                        if (fast_axis_src(i,2) >= 0.) then
                           fa_ani_phi_post(ipol,jpol,iel) = acos(arg1)
                        else
                           fa_ani_phi_post(ipol,jpol,iel) = 2. * pi - acos(arg1)
                        end if
                        
                        xi_ani_post(ipol,jpol,iel) = one
                    
                        lambdapost(ipol,jpol,iel) = rho(ipol,jpol,iel) * vptmp**2 * &
                                (1. + a_ica(i))**2 - 2. * mu(ipol,jpol,iel)

                        phi_ani_post(ipol,jpol,iel) = &
                            (1. + 2. * (a_ica(i) + b_ica(i) + c_ica(i))) / (1. + a_ica(i))

                        eta_ani_post(ipol,jpol,iel) = &
                            (vptmp**2 * (1 + a_ica(i)) * (1 + a_ica(i) + b_ica(i)) &
                                                                  - 2. * vstmp**2) &
                            / (vptmp**2 * (1 + a_ica(i))**2 - 2. * vstmp**2)
                    endif
                enddo
                rhetmax = max(rcoord(ipol,jpol,iel), rhetmax)
             enddo
          enddo
       endif
    enddo

    ! for plotting discrete points within heterogeneous region
    rhetmin = 0.
    thhetmin = 0.
    thhetmax = pi

end subroutine load_ica
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine load_het_discr(rho, lambda, mu, rhopost, lambdapost, mupost, hetind, &
                           xi_ani, phi_ani, eta_ani, xi_ani_post, phi_ani_post, &
                           eta_ani_post, fa_ani_theta_post, fa_ani_phi_post)
    use kdtree2_module

    double precision, dimension(0:npol,0:npol,nelem), intent(in)    :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(in)    :: lambda,mu
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rhopost, &
                                                            lambdapost, mupost

    double precision, dimension(0:npol,0:npol,nelem), intent(in), optional :: &
                                                            xi_ani, phi_ani, eta_ani
    double precision, dimension(0:npol,0:npol,nelem), intent(inout), optional :: &
                                               xi_ani_post, phi_ani_post, eta_ani_post
    double precision, dimension(0:npol,0:npol,nelem), intent(inout), optional :: &
                                               fa_ani_theta_post, fa_ani_phi_post

    integer                         :: hetind
    double precision, allocatable   :: w(:)
    integer                         :: iel, ipol, jpol, i, j
    integer                         :: num_het_pts
    double precision                :: s, z, r, th, r1, r2, r3, r4, th1, th2, th3, th4
    double precision                :: vptmp, vstmp, vpvtmp, vsvtmp,vphtmp, vshtmp, etatmp
    double precision                :: fa_theta_tmp, fa_phi_tmp
    double precision                :: rmin, rmax, thetamin, thetamax
    double precision, allocatable   :: szhet(:,:)
    double precision, allocatable   :: rhet2(:), thhet2(:)
    double precision, allocatable   :: delta_rho2(:), delta_vp2(:), delta_vs2(:)
    double precision, allocatable   :: delta_vph2(:), delta_vsh2(:), delta_vpv2(:), &
                                       delta_vsv2(:), delta_eta2(:)
    double precision, allocatable   :: vph2(:), vsh2(:), vpv2(:), vsv2(:), eta2(:)
    double precision, allocatable   :: fa_theta2(:), fa_phi2(:)
    double precision, allocatable   :: rho2(:), vp2(:), vs2(:)

    type(kdtree2), pointer          :: tree
    
    write(6,*) mynum, 'reading discrete heterogeneity file...'

    if (.not. het_ani_discr(hetind) == 'iso') then

        if (.not. ani_hetero) then
           write(6,*) ''
           write(6,*) 'ERROR: anisotropic heterogeneity needs an anisotropic model -'
           write(6,*) '   either choose an anisotropic background model or'
           write(6,*) '   activate force anisotropy in inparam!'
           stop
        endif

        if (.not. (present(xi_ani) .and.  present(phi_ani) .and.  present(eta_ani) .and. &
                present(xi_ani_post) .and.  present(phi_ani_post) &
                .and.  present(eta_ani_post))) then
           write(6,*) ''
           write(6,*) 'ERROR: needs anisotropic parameters in function call!!'
           stop
        endif

    endif

    open(unit=91, file=trim(het_file_discr(hetind)))

    read(91,*) num_het_pts

    if (lpr) write(6,*) 'number of points', num_het_pts
    
    allocate(rhet2(1:num_het_pts), thhet2(1:num_het_pts))

    write(6,*) mynum, 'read coordinates & medium properties...'
    
    if (het_rel_discr(hetind) == 'rel' .and. het_ani_discr(hetind) == 'iso') then
        allocate(delta_vs2(1:num_het_pts), delta_vp2(1:num_het_pts), &
                 delta_rho2(1:num_het_pts), w(1:num_het_pts))

        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), delta_vp2(j), delta_vs2(j), delta_rho2(j)
        enddo
        
        if (lpr) write(6,*) 'percent -> fraction'
       
        delta_vp2 = delta_vp2 / 100.
        delta_vs2 = delta_vs2 / 100.
        delta_rho2 = delta_rho2 / 100.

    elseif (het_rel_discr(hetind) == 'rel' .and. het_ani_discr(hetind) == 'radial') then

        allocate(delta_vsh2(1:num_het_pts), delta_vph2(1:num_het_pts), &
                 delta_vsv2(1:num_het_pts), delta_vpv2(1:num_het_pts), &
                 delta_rho2(1:num_het_pts), delta_eta2(1:num_het_pts), w(1:num_het_pts))

        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), delta_vpv2(j), delta_vsv2(j), &
                       delta_vph2(j), delta_vsh2(j), delta_rho2(j), delta_eta2(j)
        enddo
        
        if (lpr) write(6,*) 'percent -> fraction'
       
        delta_vpv2 = delta_vpv2 / 100.
        delta_vsv2 = delta_vsv2 / 100.
        delta_vph2 = delta_vph2 / 100.
        delta_vsh2 = delta_vsh2 / 100.
        delta_rho2 = delta_rho2 / 100.
        delta_eta2 = delta_eta2 / 100.

    elseif (het_rel_discr(hetind) == 'rel' .and. het_ani_discr(hetind) == 'hex') then

        allocate(delta_vsh2(1:num_het_pts), delta_vph2(1:num_het_pts), &
                 delta_vsv2(1:num_het_pts), delta_vpv2(1:num_het_pts), &
                 delta_rho2(1:num_het_pts), delta_eta2(1:num_het_pts), &
                 fa_theta2(1:num_het_pts), fa_phi2(1:num_het_pts), w(1:num_het_pts))

        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), delta_vpv2(j), delta_vsv2(j), &
                       delta_vph2(j), delta_vsh2(j), delta_rho2(j), delta_eta2(j), &
                       fa_theta2(j), fa_phi2(j)
        enddo
        
        if (lpr) write(6,*) 'percent -> fraction'
       
        delta_vpv2 = delta_vpv2 / 100.
        delta_vsv2 = delta_vsv2 / 100.
        delta_vph2 = delta_vph2 / 100.
        delta_vsh2 = delta_vsh2 / 100.
        delta_rho2 = delta_rho2 / 100.
        delta_eta2 = delta_eta2 / 100.

    elseif (het_rel_discr(hetind) == 'abs' .and. het_ani_discr(hetind) == 'iso') then
        allocate(vs2(1:num_het_pts), vp2(1:num_het_pts), &
             rho2(1:num_het_pts), w(1:num_het_pts))
        
        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), vp2(j), vs2(j), rho2(j)
        enddo
    
    elseif (het_rel_discr(hetind) == 'abs' .and. het_ani_discr(hetind) == 'radial') then

        allocate(vsh2(1:num_het_pts), vph2(1:num_het_pts), &
                 vsv2(1:num_het_pts), vpv2(1:num_het_pts), &
                 rho2(1:num_het_pts), eta2(1:num_het_pts), w(1:num_het_pts))

        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), vpv2(j), vsv2(j), vph2(j), vsh2(j), &
                       rho2(j), eta2(j)
        enddo
    
    elseif (het_rel_discr(hetind) == 'abs' .and. het_ani_discr(hetind) == 'hex') then

        allocate(vsh2(1:num_het_pts), vph2(1:num_het_pts), &
                 vsv2(1:num_het_pts), vpv2(1:num_het_pts), &
                 rho2(1:num_het_pts), eta2(1:num_het_pts), &
                 fa_theta2(1:num_het_pts), fa_phi2(1:num_het_pts), w(1:num_het_pts))

        do j=1, num_het_pts
            read(91,*) rhet2(j), thhet2(j), vpv2(j), vsv2(j), vph2(j), vsh2(j), & 
                       rho2(j), eta2(j), fa_theta2(j), fa_phi2(j)
        enddo

    else
        write(6,*) 'ERROR: combination', het_rel_discr(hetind), 'and ', &
                    het_ani_discr(hetind), 'not yet implemented'
        stop
    endif

    close(91)

    thhet2 = thhet2 * pi / 180.
    rhet2 = rhet2 * 1000.

    ! Rotate coordinates if source is not on axis
    rmin = minval(rhet2(1:num_het_pts))
    rmax = maxval(rhet2(1:num_het_pts))
    thetamin = minval(thhet2(1:num_het_pts))
    thetamax = maxval(thhet2(1:num_het_pts))

    write(6,*) mynum, 'r min/max:', rmin / 1000., rmax / 1000.
    write(6,*) mynum, 'th min/max:', thetamin / pi * 180., thetamax / pi * 180.

    if (rot_src) then 
       write(6,*) mynum, 'rotate since source is not beneath north pole'

       do j=1, num_het_pts
          call rotate_hetero(rhet2(j), thhet2(j))
       enddo

       rmin = minval(rhet2)
       rmax = maxval(rhet2)
       thetamin = minval(thhet2)
       thetamax = maxval(thhet2)

       write(6,*) mynum, 'r min/max after rotation:', rmin / 1000., rmax / 1000.
       write(6,*) mynum, 'th min/max after rotation:', thetamin / pi * 180., &
                  thetamax / pi * 180.
    endif

    ! plot discrete input file in vtk (needs to be done only once)
    if (lpr .and. het_rel_discr(hetind) == 'rel' .and. &
            het_ani_discr(hetind) == 'iso') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, &
                                 delta_vp2=delta_vp2, delta_vs2=delta_vs2, &
                                 delta_rho2=delta_rho2)

    elseif (lpr .and. het_rel_discr(hetind) == 'rel' .and. &
            het_ani_discr(hetind) == 'radial') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, &
                                 delta_vpv2=delta_vpv2, delta_vsv2=delta_vsv2, &
                                 delta_vph2=delta_vph2, delta_vsh2=delta_vsh2, &
                                 delta_rho2=delta_rho2, delta_eta2=delta_eta2)

    elseif (lpr .and. het_rel_discr(hetind) == 'rel' .and. &
            het_ani_discr(hetind) == 'hex') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, &
                                 delta_vpv2=delta_vpv2, delta_vsv2=delta_vsv2, &
                                 delta_vph2=delta_vph2, delta_vsh2=delta_vsh2, &
                                 delta_rho2=delta_rho2, delta_eta2=delta_eta2, &
                                 fa_theta2=fa_theta2, fa_phi2=fa_phi2)

    elseif (lpr .and. het_rel_discr(hetind) == 'abs' .and. & 
            het_ani_discr(hetind) == 'iso') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, &
                                 vp2=vp2, vs2=vs2, rho2=rho2)

    elseif (lpr .and. het_rel_discr(hetind) == 'abs' .and. &
            het_ani_discr(hetind) == 'radial') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, vpv2=vpv2, &
                                 vsv2=vsv2, vph2=vph2, vsh2=vsh2, rho2=rho2, eta2=eta2)
    
    elseif (lpr .and. het_rel_discr(hetind) == 'abs' .and. &
            het_ani_discr(hetind) == 'hex') then

        call plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, vpv2=vpv2, &
                                 vsv2=vsv2, vph2=vph2, vsh2=vsh2, rho2=rho2, eta2=eta2, &
                                 fa_theta2=fa_theta2, fa_phi2=fa_phi2)

    endif


    write(6,*) 'r het min/max:', rhetmin / 1000., rhetmax / 1000.
    write(6,*) 'th het min/max:', thhetmin / pi * 180., thhetmax / pi * 180.

    ! revert to cylindrical 
    allocate (szhet(2,1:num_het_pts))
    szhet(1,:) = rhet2 * sin(thhet2) 
    szhet(2,:) = rhet2 * cos(thhet2)

    tree => kdtree2_create(real(szhet), sort=.false., rearrange=.true.) 

    write(6,*) mynum, 'locate GLL points within heterogeneous regions & '

    do iel=1, nelem
       call compute_coordinates(s, z, r1, th1, iel, npol, npol)
       call compute_coordinates(s, z, r2, th2, iel, 0, 0)
       call compute_coordinates(s, z, r3, th3, iel, 0, npol)
       call compute_coordinates(s, z, r4, th4, iel, npol, 0)

       r = max(max(r1,r2), max(r3, r4))
       th = max(max(th1,th2), max(th3, th4))
    
       if (r >= rmin .and. th >= thetamin) then
          r = min(min(r1,r2), min(r3, r4))
          th = min(min(th1,th2), min(th3, th4))
          if (r <= rmax .and. th <= thetamax) then

             !to plot all gll points in elements effected:
             rhetmax = max(max(max(r1,r2), max(r3, r4)), rhetmax)
             thhetmax = max(max(max(th1,th2), max(th3, th4)), thhetmax)
             rhetmin = min(min(min(r1,r2), min(r3, r4)), rhetmin)
             thhetmin = min(min(min(th1,th2), min(th3, th4)), thhetmin)

             do ipol=0, npol
                do jpol=0, npol
                   
                   call compute_coordinates(s, z, r, th, iel, ipol, jpol)
                   call inverse_distance_weighting(s, z, tree, w, hetind)

                   if (het_ani_discr(hetind) == 'iso' .and. &
                            het_rel_discr(hetind) == 'rel') then
                   
                      vptmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                                   rho(ipol,jpol,iel))
                      vstmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel))

                      rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                                            (1. + sum(w * delta_rho2))

                      vptmp = vptmp * (1. + sum(w * delta_vp2))
                      vstmp = vstmp * (1. + sum(w * delta_vs2))
                   
                   elseif ((het_ani_discr(hetind) == 'radial' .or. &
                            het_ani_discr(hetind) == 'hex') .and. &
                            het_rel_discr(hetind) == 'rel') then
                   
                      vphtmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                                   rho(ipol,jpol,iel))

                      vshtmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel))
                      
                      vpvtmp = sqrt(phi_ani(ipol,jpol,iel)) * vphtmp
                      
                      vsvtmp =  vshtmp / sqrt(xi_ani(ipol,jpol,iel))

                      etatmp =  eta_ani(ipol,jpol,iel)

                      rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                                            (1. + sum(w * delta_rho2))

                      vpvtmp = vpvtmp * (1. + sum(w * delta_vpv2))
                      vsvtmp = vsvtmp * (1. + sum(w * delta_vsv2))
                      vphtmp = vphtmp * (1. + sum(w * delta_vph2))
                      vshtmp = vshtmp * (1. + sum(w * delta_vsh2))
                      etatmp = etatmp * (1. + sum(w * delta_eta2))

                      if (het_ani_discr(hetind) == 'hex') then
                         fa_theta_tmp = sum(w * fa_theta2)
                         fa_phi_tmp = sum(w * fa_phi2)
                      endif
                   
                   elseif (het_ani_discr(hetind) == 'iso' .and. &
                            het_rel_discr(hetind) == 'abs') then
                      
                      rhopost(ipol,jpol,iel) = sum(w * rho2)
                      vptmp = sum(w * vp2)
                      vstmp = sum(w * vs2)

                   elseif ((het_ani_discr(hetind) == 'radial' .or. &
                            het_ani_discr(hetind) == 'hex') .and. &
                            het_rel_discr(hetind) == 'abs') then
                   
                      rhopost(ipol,jpol,iel) = sum(w * rho2)
                      vpvtmp = sum(w * vpv2)
                      vsvtmp = sum(w * vsv2)
                      vphtmp = sum(w * vph2)
                      vshtmp = sum(w * vsh2)
                      etatmp = sum(w * eta2)
                      
                      if (het_ani_discr(hetind) == 'hex') then
                         fa_theta_tmp = sum(w * fa_theta2)
                         fa_phi_tmp = sum(w * fa_phi2)
                      endif
                   
                   endif

                   if (het_ani_discr(hetind) == 'iso') then

                      mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                      lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * &
                                                   (vptmp**2 - 2. * vstmp**2)

                   elseif (het_ani_discr(hetind) == 'radial' .or. &
                           het_ani_discr(hetind) == 'hex') then
                      
                      mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vshtmp**2
                      lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * &
                                                   (vphtmp**2 - 2. * vshtmp**2)

                      xi_ani_post(ipol,jpol,iel) = vshtmp**2 / vsvtmp**2
                      phi_ani_post(ipol,jpol,iel) = vpvtmp**2 / vphtmp**2

                      eta_ani_post(ipol,jpol,iel) = etatmp
                      
                      if (het_ani_discr(hetind) == 'hex') then
                         fa_ani_theta_post(ipol,jpol,iel) = fa_theta_tmp
                         fa_ani_phi_post(ipol,jpol,iel) = fa_theta_tmp
                      endif
                   endif 
                enddo
             enddo
          endif
       endif
    enddo
    
    call kdtree2_destroy(tree)  

    write(6,*) mynum, 'DONE loading discrete grid'
    
    deallocate(rhet2, thhet2)
    !deallocate(shet, zhet)
    deallocate(szhet)

    if (het_ani_discr(hetind) == 'iso' .and. &
            het_rel_discr(hetind) == 'rel') then

        deallocate(delta_vs2, delta_vp2, delta_rho2)

    elseif (het_ani_discr(hetind) == 'radial' .and. &
            het_rel_discr(hetind) == 'rel') then

        deallocate(delta_vsv2, delta_vpv2, delta_rho2, delta_vsh2, delta_vph2, delta_eta2)

    elseif (het_ani_discr(hetind) == 'hex' .and. &
            het_rel_discr(hetind) == 'rel') then

        deallocate(delta_vsv2, delta_vpv2, delta_rho2, delta_vsh2, delta_vph2, &
                   delta_eta2, fa_theta2, fa_phi2)

    elseif (het_ani_discr(hetind) == 'iso' .and. &
            het_rel_discr(hetind) == 'abs') then

        deallocate(vs2, vp2, rho2)

    elseif (het_ani_discr(hetind) == 'radial' .and. &
            het_rel_discr(hetind) == 'abs') then

        deallocate(vsv2, vpv2, rho2, vsh2, vph2, eta2)
    
    elseif (het_ani_discr(hetind) == 'hex' .and. &
            het_rel_discr(hetind) == 'rel') then

        deallocate(vsv2, vpv2, rho2, vsh2, vph2, eta2, fa_theta2, fa_phi2)

    endif

end subroutine load_het_discr
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine inverse_distance_weighting(s0, z0, tree, w, hetind)
    use kdtree2_module

    double precision, intent(in)    :: s0, z0
    type(kdtree2), pointer          :: tree
    double precision, intent(out)   :: w(1:tree%n)
    integer, intent(in)             :: hetind

    double precision                :: d2d
    integer                         :: i, nfound
    real(kdkind), dimension(2)      :: qv

    type(kdtree2_result), allocatable :: results(:) 

    w = 0.

    qv(1) = s0
    qv(2) = z0

    if (p_inv_dist(hetind) <= 100.) then
        if (R_inv_dist(hetind) == 0.) then
            ! http://en.wikipedia.org/wiki/Inverse_distance_weighting#Basic_Form
            do i=1, tree%n
               !d2d = sqrt((s(i) - s0)**2 + (z(i) - z0)**2)
               d2d = sqrt((tree%the_data(1,i) - s0)**2 + (tree%the_data(2,i) - z0)**2)
               w(i) = d2d**(-p_inv_dist(hetind))
            enddo
        else
            ! http://en.wikipedia.org/wiki/Inverse_distance_weighting#Modified_Shepard.27s_Method
            allocate(results(tree%n)) 
            call kdtree2_r_nearest(tp=tree, qv=qv, r2=real(R_inv_dist(hetind)**2), &
                    nfound=nfound, nalloc=tree%n, results=results)
            do i=1, nfound
                w(results(i)%idx) = ((R_inv_dist(hetind) - results(i)%dis) / & 
                                     (R_inv_dist(hetind) * results(i)%dis)) &
                                    **p_inv_dist(hetind)
            enddo
        endif
    else
        allocate(results(1)) 
        call kdtree2_n_nearest(tp=tree, qv=qv, nn=1, results=results)
        w(results(1)%idx) = 1.
    endif
   
    if (sum(w) > 0) then
        w = w / sum(w)
    else
        w = 0.
    endif


end subroutine inverse_distance_weighting
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine plot_discrete_input(hetind, num_het_pts, rhet2, thhet2, delta_vp2, delta_vs2, &
                               delta_rho2, delta_vpv2, delta_vsv2, &
                               delta_vph2, delta_vsh2, delta_eta2, vp2, vs2, rho2, &
                               vpv2, vsv2, vph2, vsh2, eta2, fa_theta2, fa_phi2)


    use background_models,  only : velocity
    use data_mesh,          only : discont, bkgrdmodel

    integer, intent(in)                     :: hetind, num_het_pts
    double precision, intent(in)            :: rhet2(:), thhet2(:)
    double precision, intent(in), optional  :: delta_vp2(:), delta_vs2(:), delta_rho2(:)
    double precision, intent(in), optional  :: delta_vpv2(:), delta_vsv2(:), &
                                               delta_vph2(:), delta_vsh2(:), delta_eta2(:)
    double precision, intent(in), optional  :: vpv2(:), vsv2(:), vph2(:), vsh2(:), eta2(:)
    double precision, intent(in), optional  :: vp2(:), vs2(:), rho2(:)
    double precision, intent(in), optional  :: fa_theta2(:), fa_phi2(:)

    integer                             :: j, idom
    real, allocatable, dimension(:,:)   :: meshtmp
    real, allocatable, dimension(:)     :: vptmp, vstmp, rhotmp
    real, allocatable, dimension(:)     :: vpvtmp, vsvtmp, vphtmp, vshtmp, etatmp
    real, allocatable, dimension(:)     :: fa_theta_tmp, fa_phi_tmp
    character(len=200)                  :: fname

    allocate(meshtmp(num_het_pts,2))

    if (het_ani_discr(hetind) == 'iso' .and. het_rel_discr(hetind) == 'rel') then
        allocate(vptmp(num_het_pts), vstmp(num_het_pts), rhotmp(num_het_pts))

        do j=1, num_het_pts
           idom = minloc((discont - rhet2(j)),1, mask=(discont - rhet2(j))>=0.)

           vptmp(j) = velocity(rhet2(j), 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
           vptmp(j) = vptmp(j) * (1. + delta_vp2(j))

           vstmp(j) = velocity(rhet2(j), 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
           vstmp(j) = vstmp(j) * (1. + delta_vs2(j))

           rhotmp(j) = velocity(rhet2(j), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
           rhotmp(j) = rhotmp(j)* (1. + delta_rho2(j))
        enddo

    elseif ((het_ani_discr(hetind) == 'radial' .or. het_ani_discr(hetind) == 'hex') &
            .and. het_rel_discr(hetind) == 'rel') then

        allocate(vpvtmp(num_het_pts), vsvtmp(num_het_pts), rhotmp(num_het_pts), &
                 vphtmp(num_het_pts), vshtmp(num_het_pts), etatmp(num_het_pts))

        do j=1, num_het_pts
           idom = minloc((discont - rhet2(j)),1, mask=(discont - rhet2(j))>=0.)

           vpvtmp(j) = velocity(rhet2(j), 'vpv', idom, bkgrdmodel, lfbkgrdmodel)
           vpvtmp(j) = vpvtmp(j) * (1. + delta_vpv2(j))
           
           vphtmp(j) = velocity(rhet2(j), 'vph', idom, bkgrdmodel, lfbkgrdmodel)
           vphtmp(j) = vphtmp(j) * (1. + delta_vph2(j))

           vsvtmp(j) = velocity(rhet2(j), 'vsv', idom, bkgrdmodel, lfbkgrdmodel)
           vsvtmp(j) = vsvtmp(j) * (1. + delta_vsv2(j))

           vshtmp(j) = velocity(rhet2(j), 'vsh', idom, bkgrdmodel, lfbkgrdmodel)
           vshtmp(j) = vshtmp(j) * (1. + delta_vsh2(j))

           rhotmp(j) = velocity(rhet2(j), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
           rhotmp(j) = rhotmp(j)* (1. + delta_rho2(j))
           
           etatmp(j) = velocity(rhet2(j), 'eta', idom, bkgrdmodel, lfbkgrdmodel)
           etatmp(j) = etatmp(j)* (1. + delta_eta2(j))

        enddo
        
        if (het_ani_discr(hetind) == 'hex') then
           allocate(fa_theta_tmp(num_het_pts), fa_phi_tmp(num_het_pts))

            fa_theta_tmp = fa_theta2
            fa_phi_tmp = fa_phi2
        endif
        
    elseif (het_ani_discr(hetind) == 'iso' .and. het_rel_discr(hetind) == 'abs') then
        allocate(vptmp(num_het_pts), vstmp(num_het_pts), rhotmp(num_het_pts))

        vptmp = vp2
        vstmp = vs2
        rhotmp = rho2
    
    elseif ((het_ani_discr(hetind) == 'radial' .or. het_ani_discr(hetind) == 'hex') &
            .and. het_rel_discr(hetind) == 'abs') then
        allocate(vpvtmp(num_het_pts), vsvtmp(num_het_pts), rhotmp(num_het_pts), &
                 vphtmp(num_het_pts), vshtmp(num_het_pts), etatmp(num_het_pts))

        vpvtmp = vpv2
        vphtmp = vph2
        vsvtmp = vsv2
        vshtmp = vsh2
        rhotmp = rho2
        etatmp = eta2
        
        if (het_ani_discr(hetind) == 'hex') then
           allocate(fa_theta_tmp(num_het_pts), fa_phi_tmp(num_het_pts))

            fa_theta_tmp = fa_theta2
            fa_phi_tmp = fa_phi2
        endif

    endif

    do j=1, num_het_pts
        meshtmp(j,1) = rhet2(j) * sin(thhet2(j))
        meshtmp(j,2) = rhet2(j) * cos(thhet2(j))
    enddo

    fname = trim(infopath(1:lfinfo)//'/model_rho_discr_het_'//appmynum)
    call write_VTK_bin_scal_pts(rhotmp, meshtmp, num_het_pts, fname)

    if (het_ani_discr(hetind) == 'iso') then
        fname = trim(infopath(1:lfinfo)//'/model_vp_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vptmp, meshtmp, num_het_pts, fname)

        fname = trim(infopath(1:lfinfo)//'/model_vs_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vstmp, meshtmp, num_het_pts, fname)

        if (het_rel_discr(hetind) == 'rel') then
            fname = trim(infopath(1:lfinfo)//'/model_rho_discr_het_rel_'//appmynum)
            call write_VTK_bin_scal_pts(real(delta_rho2), meshtmp, num_het_pts, fname)

            fname = trim(infopath(1:lfinfo)//'/model_vp_discr_het_rel_'//appmynum)
            call write_VTK_bin_scal_pts(real(delta_vp2), meshtmp, num_het_pts, fname)

            fname = trim(infopath(1:lfinfo)//'/model_vs_discr_het_rel_'//appmynum)
            call write_VTK_bin_scal_pts(real(delta_vs2), meshtmp, num_het_pts, fname)
        endif

    elseif (het_ani_discr(hetind) == 'radial' .or. het_ani_discr(hetind) == 'hex') then
        fname = trim(infopath(1:lfinfo)//'/model_vpv_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vpvtmp, meshtmp, num_het_pts, fname)

        fname = trim(infopath(1:lfinfo)//'/model_vsv_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vsvtmp, meshtmp, num_het_pts, fname)

        fname = trim(infopath(1:lfinfo)//'/model_vph_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vphtmp, meshtmp, num_het_pts, fname)

        fname = trim(infopath(1:lfinfo)//'/model_vsh_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(vshtmp, meshtmp, num_het_pts, fname)
        
        fname = trim(infopath(1:lfinfo)//'/model_eta_discr_het_'//appmynum)
        call write_VTK_bin_scal_pts(etatmp, meshtmp, num_het_pts, fname)
        
        if (het_ani_discr(hetind) == 'hex') then
           fname = trim(infopath(1:lfinfo)//'/model_fa_theta_discr_het_'//appmynum)
           call write_VTK_bin_scal_pts(fa_theta_tmp, meshtmp, num_het_pts, fname)
           
           fname = trim(infopath(1:lfinfo)//'/model_fa_phi_discr_het_'//appmynum)
           call write_VTK_bin_scal_pts(fa_phi_tmp, meshtmp, num_het_pts, fname)
        endif
    endif

    if (het_ani_discr(hetind) == 'iso') then
        deallocate(vptmp, vstmp)
    elseif (het_ani_discr(hetind) == 'radial') then
        deallocate(vpvtmp, vsvtmp, vphtmp, vshtmp, etatmp)
        if (het_ani_discr(hetind) == 'hex') then
            deallocate(fa_theta_tmp, fa_phi_tmp)
        endif
    endif

    deallocate(meshtmp, rhotmp)

end subroutine plot_discrete_input
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine load_random(rho,lambda,mu,rhopost,lambdapost,mupost,hetind)

    use commun
    use data_mesh, only : naxel, ax_el
    use utlity, only :  thetacoord, rcoord, zcoord

    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda, mu
    double precision, dimension(0:npol,0:npol,nelem), intent(out)   :: rhopost, &
                                                                       lambdapost, mupost
    integer         :: hetind
    real(kind=8)    :: t, decay, shift_fact, max_delta_vp, max_delta_vs, max_delta_rho
    real(kind=8)    :: vptmp, vstmp, rhotmp, th, r
    integer         :: iel, ipol, jpol, icount, i
    real(kind=8)    :: rand
    real(kind=8), allocatable :: r_rad(:), rand_rad(:), r_radtmp(:), rand_radtmp(:)

    
    if (het_funct_type(hetind) == '2Dgll') then
       write(6,*)'add 2D random anomalies to structure'

       ! add randomly to each 2D point : laterally heterogeneous and same random
       !                                 number to vp,vs,rho
       do iel=1, nelem
          th = thetacoord(npol/2,npol/2,iel)
          r = rcoord(npol/2,npol/2,iel)

          if (th >= th_het1(hetind) .and. th <= th_het2(hetind) .and. &
             r >= r_het1(hetind) .and. r <= r_het2(hetind)) then

             do jpol=0, npol
                do ipol=0, npol
                   call random_number(rand)
                   rand = 2. * rand - 1.
                   vstmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                   vptmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) /&
                                rho(ipol,jpol,iel) )

                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(hetind)*rand)

                   vptmp = vptmp * (1. + delta_vp(hetind) * rand)
                   vstmp = vstmp * (1. + delta_vs(hetind) * rand)
                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * (vptmp**2 - two*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                enddo
             enddo
          endif
       enddo

    elseif (het_funct_type(hetind) == '2Delem') then
       write(6,*)'add 2D random anomalies to structure'

       ! add randomly to each 2D point : laterally heterogeneous and same random
       !                                 number to vp,vs,rho
       do iel=1, nelem
          th = thetacoord(npol/2,npol/2,iel)
          r = rcoord(npol/2,npol/2,iel)

          if (th >= th_het1(hetind) .and. th <= th_het2(hetind) .and. &
             r >= r_het1(hetind) .and. r <= r_het2(hetind)) then

             call random_number(rand)
             rand = 2. * rand - 1.

             do jpol=0, npol
                do ipol=0, npol
                   vstmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                   vptmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) /&
                                rho(ipol,jpol,iel) )

                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(hetind)*rand)

                   vptmp = vptmp*(1. + delta_vp(hetind) * rand)
                   vstmp = vstmp*(1. + delta_vs(hetind) * rand)
                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * (vptmp**2 - two*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                enddo
             enddo
          endif
       enddo

    elseif (het_funct_type(hetind) == '1Dgll') then
       write(6,*)'add 1D random anomalies to structure per GLL point'
       ! go along axis to find the radial profile
       allocate(r_radtmp(naxel*(npol+1)), rand_radtmp(naxel*(npol+1)))
       if (mynum == 0) then 
          icount = 0
          do iel=1, naxel
             do jpol=0, npol
                if (zcoord(0,jpol,ax_el(iel)) >= 0.) then 
                   icount = icount + 1
                   r_radtmp(icount) = rcoord(0,jpol,ax_el(iel))
                   call random_number(rand)
                   rand = 2. * rand - 1.
                   rand_radtmp(icount) = rand
                endif
             enddo
          enddo
       endif 
       
       ! broadcast the profile to all processors
       call broadcast_int(icount,0)
       write(6,*) mynum, 'number of radii: ', icount
       allocate(r_rad(icount),rand_rad(icount))
       do i=1,icount
          call broadcast_dble(r_radtmp(i),0)
          r_rad(i) = r_radtmp(i)
          call broadcast_dble(rand_radtmp(i),0)
          rand_rad(i) = rand_radtmp(i)
       enddo
       
       ! add randomly to each radius, i.e. just altering the 1D background model
       do iel=1, nelem
          do jpol=0, npol
             do ipol=0,npol
                th = thetacoord(ipol,jpol,iel)
                
                i = minloc(abs(rcoord(ipol,jpol,iel)-r_rad(1:icount)),1)

                if (th >= th_het1(hetind) .and. th <= th_het2(hetind) .and. &
                   r_rad(i) >= r_het1(hetind) .and. r_rad(i) <= r_het2(hetind)) then

                   rand = rand_rad(i)
                   vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                   vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(hetind)*rand)

                   vptmp = vptmp*(1. + delta_vp(hetind)*rand)
                   vstmp = vstmp*(1. + delta_vs(hetind)*rand)

                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) *( vptmp**2 - two*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                endif
             enddo
          enddo
       enddo

    elseif (het_funct_type(hetind) == '1Delem') then
       write(6,*)'add 1D random anomalies to structure per element'

       ! go along axis to find the radial profile, only per element (not GLL point)
       allocate(r_radtmp(naxel*(npol+1)), rand_radtmp(naxel*(npol+1)))
       if (mynum==0) then 
          icount = 0
          do iel=1, naxel
             if (zcoord(0,npol/2,ax_el(iel)) >=0.) then 
                icount = icount + 1
                r_radtmp(icount) = rcoord(0,npol/2,ax_el(iel))
                call random_number(rand)
                rand = 2. * rand - 1.
                rand_radtmp(icount) = rand
             endif
          enddo
       endif 

       ! broadcast the profile to all processors
       call broadcast_int(icount,0)
       write(6,*) mynum, 'number of radii:', icount

       allocate(r_rad(icount), rand_rad(icount))
       do i=1, icount
          call broadcast_dble(r_radtmp(i),0)
          r_rad(i) = r_radtmp(i)
          call broadcast_dble(rand_radtmp(i),0)
          rand_rad(i) = rand_radtmp(i)
       enddo
       
       ! add randomly to each radius, i.e. just altering the 1D background model
       do iel=1, nelem
          i = minloc(abs(rcoord(npol/2,npol/2,iel) - r_rad(1:icount)),1)
          th = thetacoord(npol/2,npol/2,iel)

          if (th >= th_het1(hetind) .and. th <= th_het2(hetind) .and. &
             r_rad(i) >= r_het1(hetind) .and. r_rad(i) <= r_het2(hetind)) then

             rand = rand_rad(i)
             do jpol=0, npol
                do ipol=0, npol   
                   vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                   vptmp = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                                  rho(ipol,jpol,iel) )

                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * (1. + delta_rho(hetind) * rand)

                   vptmp = vptmp * (1. + delta_vp(hetind) * rand)
                   vstmp = vstmp * (1. + delta_vs(hetind) * rand)  

                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * ( vptmp**2 - two*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                enddo
             enddo
          endif
       enddo
    else
        write(6,*) 'ERROR: bad parameter: ', het_funct_type(hetind)
        stop
    endif

    !min/max of heterogeneous region
    rhetmin = min(r_het1(hetind), rhetmin)
    rhetmax = max(r_het2(hetind), rhetmax)
    thhetmin = min(th_het1(hetind), thhetmin)
    thhetmax = max(th_het2(hetind), thhetmax)

end subroutine load_random
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine load_het_funct(rho, lambda, mu, rhopost, lambdapost, mupost, hetind)

  ! added by fanie for sharp discontinuites
    use background_models, only : velocity
    use data_mesh, only : discont, bkgrdmodel

    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda, mu
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rhopost, &
                                                                       lambdapost, mupost
    integer, intent(in)                                             :: hetind

    double precision :: t, decay, shift_fact, max_delta_vp, max_delta_vs, max_delta_rho
    double precision :: vptmp, vstmp, rhotmp, s, z, r, th, gauss_val
    double precision :: r_center_gauss, th_center_gauss
    double precision :: s_center_gauss, z_center_gauss, halfwidth_r, halfwidth_th
    double precision, allocatable :: rhet(:), thhet(:), phhet(:) 

    ! start elastic property values
    double precision :: vpst, vsst, rhost
    integer :: iel, ipol, jpol, icount, jj, ij, iel_count, idom
    logical :: foundit
    double precision :: r1, r2, r3, r4, th1, th2, th3, th4
    double precision :: rmin, rmax, thetamin, thetamax

    ! for gradient
    double precision :: grad_halfwidth_r, grad_halfwidth_th
    double precision :: grad_r_het, grad_th_het2, grad_th_het1
    double precision :: dr_outer, dr_inner, dth_outer, dth_inner
    double precision :: val, dro, dru, grad_th_val, gradwidth, grad_val

    if (het_funct_type(hetind) == 'gauss') then 
       decay = 3.5d0
       r_center_gauss = (r_het1(hetind) + r_het2(hetind)) / 2.
       th_center_gauss = (th_het1(hetind) + th_het2(hetind)) / 2. * r_center_gauss
       halfwidth_r = abs(r_het1(hetind) - r_het2(hetind))
       halfwidth_th = abs(th_het1(hetind) - th_het2(hetind)) * r_center_gauss

       if (lpr) then 
          write(6,*) hetind, 'center r,th gauss [km]:', r_center_gauss / 1000., &
                     th_center_gauss / 1000.
          write(6,*) hetind, 'halfwidth r,th gauss [km]:', halfwidth_r / 1000., &
                     halfwidth_th / 1000.
       endif

       allocate(rhet(nelem*(npol+1)**2), thhet(nelem*(npol+1)**2))
       icount = 0
       rhet = 0.
       thhet = 0.
       do iel=1, nelem
          do jpol=0, npol
             do ipol=0, npol
                call compute_coordinates(s,z,r,th,iel,ipol,jpol)
                gauss_val = dexp(-( decay* ( ((r-r_center_gauss)/halfwidth_r)**2 + &
                                  ((th*r_center_gauss-th_center_gauss)/halfwidth_th)**2 )))
                vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / &
                              rho(ipol,jpol,iel) )
                vstmp = vstmp  * (1. + delta_vs(hetind)*gauss_val)
                vptmp = vptmp  * (1. + delta_vp(hetind)*gauss_val)

                if (gauss_val> 0.01) then 
                   icount = icount + 1
                   rhet(icount) = r
                   thhet(icount) = th
                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                                            (1. + delta_rho(hetind)*gauss_val)
                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * &
                                               (vptmp**2 - 2.*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * (vstmp**2)
                endif
             enddo
          enddo
       enddo
       
       !min/max of heterogeneous region
       rhetmin = min(minval(rhet(1:icount)), rhetmin)
       rhetmax = max(maxval(rhet(1:icount)), rhetmax)
       thhetmin = min(minval(thhet(1:icount)), thhetmin)
       thhetmax = max(maxval(thhet(1:icount)), thhetmax)

       write(6,*) mynum, 'r het min/max:', rhetmin / 1000., rhetmax / 1000.
       write(6,*) mynum, 'th het min/max:', thhetmin * 180. / pi, thhetmax * 180. / pi

       deallocate(rhet,thhet)

    elseif (het_funct_type(hetind)=='spher') then 
       ! gaussian with cutoff at 0.9
       decay = 0.5d0
       r_center_gauss = (r_het1(hetind) + r_het2(hetind)) / 2.
       th_center_gauss = (th_het1(hetind) + th_het2(hetind)) / 2. * r_center_gauss
       halfwidth_r = abs(r_het1(hetind) - r_het2(hetind))
       halfwidth_th = abs(th_het1(hetind) - th_het2(hetind)) * r_center_gauss
       if (lpr) then
          write(6,*) hetind, 'center r,th gauss [km]:', r_center_gauss / 1000., &
                     th_center_gauss / 1000.
          write(6,*) hetind, 'halfwidth r,th gauss [km]:', halfwidth_r / 1000., &
                     halfwidth_th / 1000.
       endif

       allocate(rhet(nelem*(npol+1)**2), thhet(nelem*(npol+1)**2))
       icount = 0
       rhet = 0.
       thhet = 0.
       do iel=1, nelem
          do jpol=0, npol
             do ipol=0, npol
                call compute_coordinates(s,z,r,th,iel,ipol,jpol)
                gauss_val = dexp(-( decay* ( ((r-r_center_gauss)/halfwidth_r)**2 + &
                                 ((th*r_center_gauss-th_center_gauss)/halfwidth_th)**2 )))

                if (gauss_val>0.9) then
                   vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                   vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / &
                                 rho(ipol,jpol,iel) )
                   vstmp = vstmp  * (1. + delta_vs(hetind))
                   vptmp = vptmp  * (1. + delta_vp(hetind))
                   rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * (1. + delta_rho(hetind))
                   lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * &
                                               (vptmp**2 - 2.*vstmp**2)
                   mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * (vstmp**2)

                   if (gauss_val > 0.9) then 
                      icount = icount + 1
                      rhet(icount) = r
                      thhet(icount) = th
                   endif
                endif
             enddo
          enddo
       enddo
       
       !min/max of heterogeneous region
       rhetmin = min(minval(rhet(1:icount)), rhetmin)
       rhetmax = max(maxval(rhet(1:icount)), rhetmax)
       thhetmin = min(minval(thhet(1:icount)), thhetmin)
       thhetmax = max(maxval(thhet(1:icount)), thhetmax)

       write(6,*) mynum, 'r het min/max:', rhetmin / 1000., rhetmax / 1000.
       write(6,*) mynum, 'th het min/max:', thhetmin * 180. / pi, thhetmax * 180. / pi

       deallocate(rhet,thhet)

    elseif (het_funct_type(hetind)=='sinus' .or. het_funct_type(hetind)=='trian' &
            .or. het_funct_type(hetind)=='gss1d' .or. het_funct_type(hetind)=='inclp' & 
            .or. het_funct_type(hetind)=='inclr' .or. het_funct_type(hetind)=='const' ) then 
       decay = 1.d0
       ! included by fanie, for a sharp triangle/ sine function/ gauss function/ sphere

       ! if rdep==.false. v_het(z)=const and rho_het(z)=const
       ! determine velocity ontop of heterogeneity and calculate future elastic
       ! properties of the heterogeneity

       if (.not. rdep(hetind) ) then
          if ( inverseshape(hetind)>0 ) then
             idom = minloc(abs(discont - r_het2(hetind)),1)
             !do ij = 1, ndisc
             !   write(6,*)'discmin?',ij,abs(discont(ij)-r_het2(hetind))
             !enddo
             !write(6,*)hetind,'DEBUGGING:',idom,r_het1(hetind),r_het2(hetind)
             !write(6,*)hetind,'discs:',idom,discont(idom),bkgrdmodel,lfbkgrdmodel

!#########################################################################################
! MvD: - calling velocity() here causes problems with anisotropy, anway it is
!        called already in get_model, so why not use the lame parameters here?
!      - vsst / vpst only used if not gradient, so should be easy to copy from
!        the gradient version
!#########################################################################################

             vpst = velocity(r_het2(hetind),'v_p',idom,bkgrdmodel,lfbkgrdmodel)
             !vpst = vpst* (1.+ delta_vp(hetind))
             vsst = velocity(r_het2(hetind),'v_s',idom,bkgrdmodel,lfbkgrdmodel)
             !vsst = vsst* (1.+delta_vs(hetind))
             rhost = velocity(r_het2(hetind),'rho',idom,bkgrdmodel,lfbkgrdmodel)
             !rhost = rhost* (1.+delta_rho(hetind))
             if (lpr) then 
                write(6,*) 'depth-independent variation!'
                write(6,*) hetind, 'elastic properties set to (vp,vs,rho):', &
                           r_het2, vpst, vsst, rhost
             endif
          else
             idom = minloc(abs(discont - r_het1(hetind)),1)
             !do ij = 1, ndisc
             !   write(6,*)'discmin?',ij,abs(discont(ij)-r_het1(hetind))
             !enddo
             !write(6,*)hetind,'DEBUGGING:',idom,r_het1(hetind),r_het2(hetind)
             !write(6,*)hetind,'discs:',idom,discont(idom),bkgrdmodel,lfbkgrdmodel
             vpst = velocity(r_het1(hetind),'v_p',idom,bkgrdmodel,lfbkgrdmodel)
             !vpst = vpst* (1.+ delta_vp(hetind))
             vsst = velocity(r_het1(hetind),'v_s',idom,bkgrdmodel,lfbkgrdmodel)
             !vsst = vsst* (1.+delta_vs(hetind))
             rhost = velocity(r_het1(hetind),'rho',idom,bkgrdmodel,lfbkgrdmodel)
             !rhost = rhost* (1.+delta_rho(hetind))
             if (lpr) then
                write(6,*) 'depth-independent variation!'
                write(6,*) hetind, 'elastic properties set to (vp,vs,rho):', &
                           r_het1, vpst, vsst, rhost
             endif
          endif
       endif

       ! define width and height of structure
       ! if gradient
       if ( grad(hetind) ) then
          ! r_het1: lower boundary of complete thingy, r_het2: upper boundary
          ! r_het2 > grad_r_het: upper boundary of gradient
          ! r_het2 = r_het2 - grad_width: upper boundary of constant part of
          !                               heterogeneity, lower boundary of gradient
          grad_r_het = r_het2(hetind)
          r_het2(hetind) = r_het2(hetind) - gradrdep1(hetind) * 1000.
          ! th_het1 > grad_th_het1: left bound of whole heterogeneity
          ! th_het2 > grad_th_het2: right bound of whole heterogeneity
          ! th_het1 = th_het1 + gradrdep2/pi/r_het1 : 
          !                     left boundary of constant part of het,
          !                     right boundary of left gradient
          ! th_het2 = th_het2 - gradrdep2/pi/r_het1 : 
          !                     right boundary of constant part, 
          !                     left boundary of right gradient
          grad_th_het1 = th_het1(hetind)
          grad_th_het2 = th_het2(hetind)
          th_het1(hetind) = th_het1(hetind) + gradrdep2(hetind) * 1000. / r_het1(hetind)
          th_het2(hetind) = th_het2(hetind) - gradrdep2(hetind) * 1000. / r_het1(hetind)
       endif
    
       ! use this as center of heterogeneity
       r_center_gauss = (r_het1(hetind)+r_het2(hetind)) / 2.
       th_center_gauss = (th_het1(hetind)+th_het2(hetind)) / 2. !*r_center_gauss
       ! height and width of anomaly
       halfwidth_r = abs(r_het1(hetind)-r_het2(hetind))
       halfwidth_th = abs(th_het1(hetind)-th_het2(hetind)) !*r_center_gauss
       
       if ( grad(hetind) ) then
          grad_halfwidth_r = abs(grad_r_het - r_het1(hetind))
          grad_halfwidth_th = abs(grad_th_het1-grad_th_het2) !*r_center_gauss
       endif

       if (lpr) then
          write(6,*) hetind, 'center r [km],th [deg]:', het_funct_type(hetind), &
                     r_center_gauss / 1000., th_center_gauss*180./pi
          write(6,*) hetind, 'halfwidth r [km],th [deg]:', het_funct_type(hetind), &
                     halfwidth_r / 1000., th_het1(hetind) * 180. / pi, &
                     th_het2(hetind) * 180. / pi, halfwidth_th * 180. / pi
       endif

       if ( grad(hetind) .and. lpr ) then
          write(6,*) hetind, 'gradhalfwidth r [km],th [deg]:', het_funct_type(hetind), &
                     grad_halfwidth_r / 1000., grad_halfwidth_th * 180. / pi
       endif

       ! find elements within heterogeneity
       allocate(rhet(nelem*(npol+1)**2), thhet(nelem*(npol+1)**2))
       
       icount = 0
       rhet = 0.
       thhet = 0.
       jj = 0.
       do iel=1, nelem
          foundit = .false.
          call compute_coordinates(s,z,r1,th1,iel,0,0)
          call compute_coordinates(s,z,r2,th2,iel,0,npol)
          call compute_coordinates(s,z,r3,th3,iel,npol,0)
          call compute_coordinates(s,z,r4,th4,iel,npol,npol)

          rmin = 1.001 * min(r1,r2,r3,r4)
          thetamin = 1.001 * min(th1,th2,th3,th4)
          rmax = 0.999 * max(r1,r2,r3,r4)
          thetamax = 0.999 * max(th1,th2,th3,th4)

          if ( rmin >= r_het1(hetind) ) then
             if (((.not. grad(hetind)) .and. rmax<=r_het2(hetind) &
                .and. thetamax<=th_het2(hetind) .and. thetamin>=th_het1(hetind)) .or. &
                ( grad(hetind)  .and. rmax<=grad_r_het .and. thetamax<=grad_th_het2 &
                .and. thetamin>=grad_th_het1)) then 

                jj = hetind
                foundit = .true.
                iel_count = iel_count + 1

                !massive spam of OUTPUT file:
                !write(6,*) mynum, 'found element inside hetero box:', iel, jj
                !write(6,*) mynum, 'r,th min', rmin / 1000., thetamin * 180. / pi
                !write(6,*) mynum, 'r,th max', rmax / 1000., thetamax * 180. / pi
             endif
          endif
            
          !XXX###################################################################
          ! MvD: - is it smart to open the file inside the loop over all elements?
          !      - the file is never closed!
          !      - het_funct_type(hetind) can never be 'spher' in this part of
          !        the code...
          !######################################################################

          open(267, file='hetero_function_adaptions.dat')

          if (het_funct_type(hetind)=='gss1d' .or. het_funct_type(hetind)=='spher') &
             open(268,file='hetero_function_sphere.dat')
          if (foundit) then
             do jpol=0, npol
                do ipol=0, npol
                   call compute_coordinates(s,z,r,th,iel,ipol,jpol)
                   ! check if within function
                   val = 0.
                   dr_outer = 0. ! distance to outer boundary (positive if inside heterogeneity)
                   dr_inner = 0. ! distance to inner boundary (positive if inside constant part)
                   dth_outer = 0.
                   dth_inner = 0.
                   ! horizontal width of gradient at depth r
                   gradwidth = 0.

                   if ( grad(hetind) ) then
                       gradwidth = (gradrdep1(hetind) - gradrdep2(hetind)) * 1000. / &
                                    r_het1(hetind) / &
                                    grad_halfwidth_r * ( r - r_het1(hetind) ) + &
                                    gradrdep2(hetind) * 1000. / r_het1(hetind)
                   endif

                   if (het_funct_type(hetind) == 'const') then
                      dr_inner = r_het2(hetind) - r
                      if ( th < th_center_gauss ) &
                         dth_inner = th - (th_het1(hetind) + gradwidth)

                      if ( th >= th_center_gauss ) &
                         dth_inner = (th_het2(hetind) - gradwidth) - th

                      if ( grad(hetind) ) then
                         dr_outer = grad_r_het - r
                         if (th < th_center_gauss) dth_outer = th - grad_th_het1
                         if (th >= th_center_gauss) dth_outer = grad_th_het2 - th
                      endif
                   endif

                   if (het_funct_type(hetind)=='sinus') then 
                      dr_inner = r_het1(hetind) + halfwidth_r * &
                                 sin( (th - th_het1(hetind)) * pi / halfwidth_th ) - r
                      if ( th < th_center_gauss ) &
                         dth_inner = th - th_het1(hetind) - gradwidth + &
                                     halfwidth_th * cos((th - th_het1(hetind)) * pi / &
                                                         halfwidth_r)
                      if ( th >= th_center_gauss ) &
                         dth_inner = th_het2(hetind) - gradwidth - th + &
                                     halfwidth_th * cos((th_het2(hetind) - th) * pi /  &
                                                        halfwidth_r)
                      if ( grad(hetind) ) then
                          dr_outer = r_het1(hetind) + grad_halfwidth_r * &
                                     sin((th - grad_th_het1) * pi / grad_halfwidth_th ) - r
                          if ( th < th_center_gauss ) &
                             dth_outer = th - grad_th_het1 + grad_halfwidth_th * &
                                         cos((th - grad_th_het1) * pi / grad_halfwidth_r)
                          if ( th >= th_center_gauss ) &
                             dth_outer = grad_th_het2 - th + grad_halfwidth_th * &
                                         cos((grad_th_het2 - th) * pi / grad_halfwidth_r)
                      endif
                   endif

                   if (het_funct_type(hetind) == 'trian') then
                      dr_inner = r_het1(hetind) + 2 * halfwidth_r / halfwidth_th * &
                                 ((th_center_gauss - abs(th - th_center_gauss)) - & 
                                 th_het1(hetind)) - r
                      dth_inner = th - gradwidth + ((th_center_gauss - abs(th - &
                                  th_center_gauss)) - th_het1(hetind))
                      if ( grad(hetind) ) then
                         dr_outer = r_het1(hetind) + & 
                                    2 * grad_halfwidth_r / grad_halfwidth_th * &
                                    ((th_center_gauss - abs(th - th_center_gauss)) - &
                                    grad_th_het1) - r
                         dth_outer = th - gradwidth + &
                                     ((th_center_gauss - abs(th - th_center_gauss)) - &
                                     grad_th_het1)
                      endif
                   endif
                   
                   if (het_funct_type(hetind) == 'inclr') then 
                      dr_inner = r_het1(hetind) + ((halfwidth_r / halfwidth_th) * &
                                 abs(th - th_het1(hetind))) - r
                      dth_inner = th - gradwidth + abs(th - th_het1(hetind))
                      if ( grad(hetind) ) then
                         dr_outer = r_het1(hetind) + grad_halfwidth_r / &
                                    grad_halfwidth_th * abs(th-grad_th_het1) - r
                         dth_outer = th - gradwidth + abs( th - grad_th_het1)
                      endif
                   endif
                   
                   if (het_funct_type(hetind) == 'inclp') then
                      dr_inner = r_het1(hetind) + ((halfwidth_r / halfwidth_th) * &
                                 abs(th - th_het2(hetind))) - r
                      dth_inner = th - gradwidth + &
                                  abs(th - th_het2(hetind))
                      if ( grad(hetind) ) then
                         dr_outer = r_het1(hetind) + grad_halfwidth_r / &
                                    grad_halfwidth_th * abs(th-grad_th_het2) - r
                         dth_outer = th - gradwidth + abs(th - grad_th_het2)
                      endif
                   endif
                   
                   if (het_funct_type(hetind) == 'gss1d') then
                      dr_inner = r_het1(hetind) + halfwidth_r * &
                                 dexp(-(20. * (th - th_center_gauss)**2 / &
                                 halfwidth_th**2 / 2)) - r
                      dth_inner = th - th_het1(hetind) - gradwidth + halfwidth_th * &
                                  cos((th - th_het1(hetind)) * pi / halfwidth_r )
                      if ( grad(hetind) ) then
                          dr_outer = r_het1(hetind) + grad_halfwidth_r * &
                                     dexp(-( 20. * (th - th_center_gauss)**2 / &
                                     grad_halfwidth_th**2 / 2)) - r
                          dth_outer = th - th_het1(hetind) + grad_halfwidth_th * &
                                      cos((th - th_het1(hetind)) * pi / grad_halfwidth_r)
                      endif
                      write(268,*) r, th, dr_inner, dr_outer
                   endif
                   
          !######################################################################
          ! MvD: - het_funct_type(hetind) can never be 'spher' in this part of
          !        the code...
          !      - isnt that implemented above already?
          !######################################################################

                   ! sphere edges with 0.9 of gauss function
                   if (het_funct_type(hetind) == 'spher') then
                      dr_inner = dexp(-((0.5*((r - r_center_gauss)**2 / halfwidth_r**2 / 2) &
                                 + ((th - th_center_gauss)**2 / halfwidth_th**2 / 2))))
                      if ( grad(hetind) ) then
                         dr_outer = dexp(-((0.5*((r - r_center_gauss)**2 / &
                                    grad_halfwidth_r**2 / 2) + &
                                    ((th - th_center_gauss)**2 / grad_halfwidth_th**2 / 2))))

                         if ( dr_outer>=0.9 .and. dr_inner<0.9 ) then
                             dr_outer = dr_outer - dr_inner - 0.9
                             !( (dr_inner/ sqrt(halfwidth_r**2 + halfwidth_th**2 )) + &
                             !(dr_outer/ sqrt(grad_halfwidth_r**2 + grad_halfwidth_th**2 )) ) * &
                             !( sqrt(halfwidth_r**2 + halfwidth_th**2 )/2 + &
                             !	sqrt(grad_halfwidth_r**2 + grad_halfwidth_th**2 )/2 )
                             dth_outer = dr_outer
                         elseif ( dr_outer>=0.9 .and. dr_inner>=0.9 ) then
                             dr_outer = 0.
                             dth_outer = 0.
                         else
                             dr_outer = -1.
                             dth_outer = -1.
                         endif
                      endif
                      if ( dr_inner>=0.9 ) dr_inner = +1.
                      if ( dr_inner<0.9 ) dr_inner = -1.
                      dth_inner = dr_inner
                      write(268,*) r, th, dr_inner, dr_outer
                   endif

                   ! quantify gradual change, if r,th within gradient area
                   if ( grad(hetind) ) then !.and. r>=gauss_val .and. r<=grad_val ) then
                      if ( ( dr_outer>0 .and. dr_inner<0. ) .or. &
                         ( dth_outer>=0. .and. dth_inner<0. ) ) then
                          ! 2d-distance to gradient/constant boundary
                          val = sqrt( dr_outer**2 + dth_outer**2 ) / &
                                ( sqrt( dr_outer**2 + dth_outer**2 ) + &
                                sqrt( dr_inner**2 + dth_inner**2 ) )
                          if (het_funct_type(hetind)=='spher') val = dr_outer
                      else
                          val = 1
                      endif
                   else
                      val = 1
                      dth_outer = dth_inner
                      dr_outer = dr_inner
                   endif

                   write(267,*) hetind, icount,&
                   !r,th,gauss_val,dr_outer,dr_inner,dth_outer,dth_inner,val
                   dth_outer, th, dth_inner, val
                   !write(267,*)hetind,icount,&
                   !(r-r_het1(hetind))/(r_het2(hetind)-r_het1(hetind)),&
                   !(th-th_het1(hetind))/(th_het2(hetind)-th_het1(hetind)),&
                   !(gauss_val-r_het1(hetind))/(r_het2(hetind)-r_het1(hetind)),&
                   !rho(ipol,jpol,iel),lambda(ipol,jpol,iel),mu(ipol,jpol,iel)
        
                   if ( dr_outer>=0. .and. dth_outer>=0.) then
                      if ( rdep(hetind) ) then
                         ! gradual elastic property changes with depth
                         !write(267,*)'depth-dependent variation!'
                         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
                         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                                       rho(ipol,jpol,iel) )
                         vstmp = vstmp * (1. + val * delta_vs(hetind))
                         vptmp = vptmp * (1. + val * delta_vp(hetind))
                         rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * (1. + val * delta_rho(hetind))
                         lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * (vptmp**2 - 2. * vstmp**2)
                         mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                      else
                         ! constant elastic property changes
                         !write(267,*)'depth-independent variation!'
                         rhopost(ipol,jpol,iel) = rhost * (1. + val * delta_rho(hetind))
                         lambdapost(ipol,jpol,iel) = rhost * (1. + val * delta_rho(hetind)) * &
                                                     ((vpst * (1. + val * delta_vp(hetind)))**2 &
                                                      - 2. * (vsst * (1. + val * delta_vs(hetind)))**2)
                         mupost(ipol,jpol,iel) = rhost * (1. + val * delta_rho(hetind)) * &
                                                 ((vsst * (1. + val * delta_vs(hetind)))**2)
                      endif !rdep
                   ! for both...
                   endif !inside heterogeneity
                   !write(267,*)hetind,icount,&
                               !(r-r_het1(hetind))/(r_het2(hetind)-r_het1(hetind)),&
                   !(dr_inner-r_het1(hetind))/(r_het2(hetind)-r_het1(hetind)),&
                   !rho(ipol,jpol,iel),lambda(ipol,jpol,iel),mu(ipol,jpol,iel)
                   icount = icount + 1
                   rhet(icount) = r
                   thhet(icount) = th
                enddo !ipol
             enddo !jpol
          endif !foundit
       enddo
       
       !min/max of heterogeneous region
       rhetmin = min(minval(rhet(1:icount)), rhetmin)
       rhetmax = max(maxval(rhet(1:icount)), rhetmax)
       thhetmin = min(minval(thhet(1:icount)), thhetmin)
       thhetmax = max(maxval(thhet(1:icount)), thhetmax)
       
       write(6,*) mynum, 'r het min/max:', rhetmin/1000., rhetmax/1000.
       write(6,*) mynum, 'th het min/max:', thhetmin*180./pi, thhetmax*180./pi
       deallocate(rhet,thhet)

    else
       write(6,*) 'function type ', het_funct_type(hetind), ' not implemented yet!'
       stop
    endif

end subroutine load_het_funct
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine plot_hetero_region_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                  fa_ani_theta, fa_ani_phi)
    
    double precision, dimension(0:npol,0:npol,nelem), intent(in) :: rho, lambda, mu
    double precision, dimension(0:npol,0:npol,nelem), intent(in), optional :: &
                                    xi_ani, phi_ani, eta_ani, fa_ani_theta, fa_ani_phi

    real, dimension(:), allocatable     :: vp_all, vs_all, rho_all
    real, dimension(:), allocatable     :: vph_all, vsh_all, vpv_all, vsv_all, eta_all
    real, dimension(:), allocatable     :: xi_all, phi_all, fa_theta_all, fa_phi_all
    real, dimension(:,:), allocatable   :: mesh2
    character(len=200)                  :: fname
    double precision                    :: s, z, r, th
    integer                             :: iel, ipol, jpol, icount

    write(6,*) 'plotting heterogeneous region in pointwise vtk'

    allocate(mesh2(nelem * (npol + 1)**2,2), rho_all(nelem * (npol + 1)**2))

    if (.not. ani_hetero) then
       allocate(vp_all(nelem * (npol + 1)**2), vs_all(nelem * (npol + 1)**2))
    else
       allocate(vph_all(nelem * (npol + 1)**2), vsh_all(nelem * (npol + 1)**2), &
                vpv_all(nelem * (npol + 1)**2), vsv_all(nelem * (npol + 1)**2), &
                fa_theta_all(nelem * (npol + 1)**2), fa_phi_all(nelem * (npol + 1)**2), &
                eta_all(nelem * (npol + 1)**2), xi_all(nelem * (npol + 1)**2), &
                phi_all(nelem * (npol + 1)**2))
    endif

    if (lpr) then
       write(6,*) 'Heterogeneous region rmin,rmax [km]:', &
                  rhetmin/1000., rhetmax/1000.
       write(6,*) 'Heterogeneous region thmin,thmax [deg]:', &
                  thhetmin*180./pi, thhetmax*180./pi
    endif

    icount = 0

    do iel=1, nelem
       do ipol=0, npol
          do jpol=0, npol
             call compute_coordinates(s,z,r,th,iel,ipol,jpol)
             if (r >= rhetmin .and. r <= rhetmax .and. th >= thhetmin &
                  .and. th <= thhetmax) then 
                icount = icount + 1
                mesh2(icount,1) = real(s)
                mesh2(icount,2) = real(z)
                
                rho_all(icount) = rho(ipol,jpol,iel) 
                
                if (.not. ani_hetero) then

                   vp_all(icount) = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel) ) / &
                                          rho(ipol,jpol,iel)  )
                   vs_all(icount) = sqrt( (mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel)  )

                else
                   
                   vph_all(icount) = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel) ) / &
                                          rho(ipol,jpol,iel)  )
                   vsh_all(icount) = sqrt( (mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel)  )

                   vpv_all(icount) = sqrt(phi_ani(ipol,jpol,iel)) * vph_all(icount)
                   
                   vsv_all(icount) =  vsh_all(icount) / sqrt(xi_ani(ipol,jpol,iel))

                   eta_all(icount) =  eta_ani(ipol,jpol,iel)
                   xi_all(icount) =  xi_ani(ipol,jpol,iel)
                   phi_all(icount) =  phi_ani(ipol,jpol,iel)
                   fa_theta_all(icount) =  fa_ani_theta(ipol,jpol,iel)
                   fa_phi_all(icount) =  fa_ani_phi(ipol,jpol,iel)
                endif
             endif
          enddo
       enddo
    enddo

    write(6,*) mynum, 'number of points inside heterogeneous region:', icount

    if (icount > 0) then

        fname = trim(infopath(1:lfinfo)//'/model_rho_gll_het'//appmynum)
        call write_VTK_bin_scal_pts(real(rho_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                    icount, fname)
        if (.not. ani_hetero) then
           fname = trim(infopath(1:lfinfo)//'/model_vp_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vp_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)

           fname = trim(infopath(1:lfinfo)//'/model_vs_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vs_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)
        
           deallocate(vp_all, vs_all)
        else
           
           fname = trim(infopath(1:lfinfo)//'/model_vph_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vph_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)

           fname = trim(infopath(1:lfinfo)//'/model_vsh_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vsh_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)
        
           fname = trim(infopath(1:lfinfo)//'/model_vpv_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vpv_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)

           fname = trim(infopath(1:lfinfo)//'/model_vsv_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(vsv_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)
        
           fname = trim(infopath(1:lfinfo)//'/model_eta_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(eta_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)
           
           fname = trim(infopath(1:lfinfo)//'/model_xi_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(xi_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)
           
           fname = trim(infopath(1:lfinfo)//'/model_phi_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(phi_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                       icount, fname)

           fname = trim(infopath(1:lfinfo)//'/model_fa_theta_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(fa_theta_all(1:icount)), &
                                       real(mesh2(1:icount,1:2)), icount, fname)

           fname = trim(infopath(1:lfinfo)//'/model_fa_phi_gll_het'//appmynum)
           call write_VTK_bin_scal_pts(real(fa_phi_all(1:icount)), &
                                       real(mesh2(1:icount,1:2)), icount, fname)

           deallocate(vph_all, vsh_all, vpv_all, vsv_all, fa_theta_all, fa_phi_all, &
                      eta_all, xi_all, phi_all)
        endif
    else
        write(6,*) mynum, 'WARNING: not writing empty vtk file (icount = 0)'
    endif

    deallocate(mesh2, rho_all)

end subroutine plot_hetero_region_vtk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal_pts(u2, mesh1, rows, filename)

   real, dimension(1:rows), intent(in)      :: u2
   real, dimension(1:rows,1:2), intent(in)  :: mesh1
   integer, intent(in)                      :: rows
   character (len=200), intent(in)          :: filename

   integer                      :: i
   real, dimension(1:rows)      :: u1
   integer, dimension(1:rows*2) :: cell
   integer, dimension(1:rows)   :: cell_type
   character (len=50)           :: ss
    
   !points structure
   do i=2, rows*2, 2
      cell(i-1) = 1
      cell(i) = (i / 2) - 1
   enddo

   do i=1, rows
      cell_type(i) = 1
   enddo
  
   u1 = real(u2)
   !do i=1, rows
   !   if (abs(u1(i)) < 1.e-25) u1(i) = 0.0
   !enddo
  
   write(6,*)'computing vtk file ',trim(filename),' ...'

   open(100, file=trim(filename)//'.vtk', access='stream', status='replace', &
        convert='big_endian')
  
   write(100) '# vtk DataFile Version 4.0'//char(10)
   write(100) 'mittico'//char(10)
   write(100) 'BINARY'//char(10)
   write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
   write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
   write(100) ss//char(10)

   !points
   do i=1,rows
      write(100) mesh1(i,1),mesh1(i,2),0.0
   enddo
   write(100) char(10)

   !cell topology
   write(ss,fmt='(A5,2I10)') 'CELLS', rows, rows*2
   write(100) char(10)//ss//char(10)
   write(100) cell
   write(100) char(10)

   !cell type
   write(ss,fmt='(A10,2I10)') 'CELL_TYPES',rows
   write(100) char(10)//ss//char(10)
   write(100) cell_type
   write(100) char(10)
   
   !data
   write(ss,fmt='(A10,I10)') 'CELL_DATA',rows
   write(100) char(10)//ss//char(10)
   write(100) 'SCALARS '//trim(filename)//' float 1'//char(10)
   write(100) 'LOOKUP_TABLE default'//char(10) !color table?
   write(100) real(u1)
   close(100)
   write(6,*)'...saved ',trim(filename)//'.vtk'

end subroutine write_VTK_bin_scal_pts
!-----------------------------------------------------------------------------

end module lateral_heterogeneities
