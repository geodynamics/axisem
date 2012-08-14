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
private
contains

!----------------------------------------------------------------------------------------
subroutine compute_heterogeneities(rho,lambda,mu)
    implicit none

    include 'mesh_params.h'

    integer ::  ij
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
    double precision, dimension(0:npol,0:npol,nelem) :: rhopost,lambdapost,mupost

    mupost = mu
    lambdapost = lambda
    rhopost = rho

    if (lpr) then
       write(6,*)
       write(6,*) '   +++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*) '   ++++++++    Lateral Heterogeneities  ++++++++'
       write(6,*) '   +++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*)
       write(6,*) 'read parameter file for heterogeneities: inparam_hetero'
       write(6,*)
       write(6,*) ' !!!!!!!!! W A R N I N G !!!!!!!! '
       write(6,*) 'These lateral additions have not been thoroughly tested yet!'
       write(6,*)
    endif

    call read_param_hetero

    ! load and add heterogeneities
    ! loop over heterogeneities and call functions for each separately!
 
    do ij = 1, num_het
       if (het_format(ij) == 'const') then
          het_format(ij) = 'funct'
          het_funct_type(ij) = 'const'
          call load_het_funct(rho,lambda,mu,rhopost,lambdapost,mupost,ij)
       elseif (het_format(ij)=='funct') then
          ! functional perturbations (sinus, gauss, trian, gss1d, inclp, inclr, const)
          call load_het_funct(rho,lambda,mu,rhopost,lambdapost,mupost,ij) 
       elseif (het_format(ij) == 'discr') then
          ! interpolate discrete model of arbitrary locations/perturbations
          call load_het_discr(rho,lambda,mu,rhopost,lambdapost,mupost,ij) 
       elseif (het_format(ij) == 'rndm') then 
          ! add random fluctuations to radial model
          call load_random(rho,lambda,mu,ij) 
       else
          write(6,*)'Unknown heterogeneity input file type!!'; stop
       endif
    enddo
 
    !!% write final changes to model!
    mu = mupost
    lambda = lambdapost
    rho = rhopost
 
    write(6,*)'final model done, now vtk files...'
    call plot_hetero_region_vtk(rho,lambda,mu)

    deallocate(het_format, het_file_discr, het_funct_type, rdep, grad, &
               gradrdep1, gradrdep2, r_het1, r_het2, th_het1, th_het2, ph_het1, &
               ph_het2, delta_rho, delta_vp, delta_vs, inverseshape) 

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
    implicit none
    integer :: ij, i
    character(len=100) :: junk

    inquire(file="inparam_hetero", EXIST=file_exists)

    if (.not. file_exists) then 
       write(6,*) 'ERROR: lateral heterogeneity set in inparam, but'
       write(6,*) '       inparam_hetro does not exist!'
       stop
    endif

    write(6,*) 'starting read_param_hetero'
    open(unit=91, file='inparam_hetero')

    read(91,*)num_het
    if (lpr) write(6,*) 'adding ', num_het, ' regions...'

    allocate(het_format(num_het), het_file_discr(num_het), &
             het_funct_type(num_het), rdep(num_het), grad(num_het), &
             gradrdep1(num_het), gradrdep2(num_het), r_het1(num_het), &
             r_het2(num_het), th_het1(num_het), th_het2(num_het), &
             ph_het1(num_het), ph_het2(num_het), delta_rho(num_het), &
             delta_vp(num_het), delta_vs(num_het), inverseshape(num_het)) 

    do ij = 1, num_het
       read(91,*) junk
       read(91,*) het_format(ij)
       read(91,*) het_file_discr(ij)
       read(91,*) het_funct_type(ij)
       read(91,*) rdep(ij)
       read(91,*) grad(ij)
       read(91,*) gradrdep1(ij),gradrdep2(ij)
       read(91,*) r_het1(ij),r_het2(ij)
       read(91,*) th_het1(ij),th_het2(ij)
!#########################################################################################
! MvD: - seems like phi is never ever used, so why have it as parameter?
!#########################################################################################
       read(91,*) ph_het1(ij),ph_het2(ij)
       read(91,*) delta_rho(ij)
       read(91,*) delta_vp(ij)
       read(91,*) delta_vs(ij)
     
!#########################################################################################
! MvD: - what does inverseshape stand for?
!#########################################################################################
       inverseshape(ij) = index(het_funct_type(ij),'_i')-1
       
       if (inverseshape(ij) > 0) then
           het_funct_type(ij) = het_funct_type(ij)(1:inverseshape(ij))
       endif
    enddo

    ! degree to radians
    th_het1 = th_het1 / 180. * pi
    th_het2 = th_het2 / 180. * pi
    ph_het1 = ph_het1 / 180. * pi
    ph_het2 = ph_het2 / 180. * pi

    ! percent to decimal
    delta_rho = delta_rho / 100.
    delta_vp = delta_vp / 100.
    delta_vs = delta_vs / 100.

    if (lpr) then 
       do ij=1,num_het
          if (het_format(ij)=='funct' .or. het_format(ij)=='const') then
             write(6,*) 'Specification of heterogeneous region:', ij
             write(6,*) 'Radius (lower/upper bound) [km]:', &
                         r_het1(ij) / 1000., r_het2(ij) / 1000.
             write(6,*) 'Colatitude (lower/upper bound) [deg]:', &
                         th_het1(ij) * 180. / pi, th_het2(ij) * 180. / pi
             write(6,*) 'delta rho [%]:', delta_rho(ij)
             write(6,*) 'delta vp  [%]:', delta_vp(ij)
             write(6,*) 'delta vs  [%]:', delta_vs(ij)
          endif
       enddo
    endif

    ! need to rotate coordinates if source is not along axis (beneath the north pole)
    if (rot_src ) then 
       write(6,*) 'need to rotate the heterogeneous domain with the source....'

       do i=1, num_het
          write(6,*)'Before rotation r th ph 1:', &
             r_het1(i), th_het1(i) * 180. / pi, ph_het1(i) * 180. / pi
          write(6,*)'Before rotation r th ph 2:', &
             r_het2(i), th_het2(i) * 180. / pi, ph_het2(i) * 180. / pi
       enddo

       call rotate_hetero(num_het, 1, r_het1, th_het1, ph_het1)
       call rotate_hetero(num_het, 1, r_het2, th_het2, ph_het2)

       do i=1, num_het
          write(6,*)'After rotation r th ph 1:', &
             r_het1(i), th_het1(i) * 180. / pi, ph_het1(i) * 180. / pi
          write(6,*)'After rotation r th ph 2:', &
             r_het2(i), th_het2(i) * 180. / pi, ph_het2(i) * 180. / pi
       enddo
    endif

    write(6,*)'done with read_param_hetero'

end subroutine read_param_hetero
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine rotate_hetero(n,m,r,th,ph)
    implicit none

    integer, intent(in) :: n,m
    double precision,intent(inout), dimension(1:n,1:m) :: r,th,ph
    double precision :: x_vec(3),x_vec_rot(3),r_r, arg1
    integer :: i,j

!#########################################################################################
! MvD: - why the double index while in the other functions there is only one dimension?
!#########################################################################################

    do i=1, m
       do j=1, n
           
           x_vec(1) = r(j,i) * dsin(th(j,i)) * dcos(ph(j,i))
           x_vec(2) = r(j,i) * dsin(th(j,i)) * dsin(ph(j,i))
           x_vec(3) = r(j,i) * dcos(th(j,i)) 
           
           x_vec_rot = matmul(trans_rot_mat,x_vec)
           
           write(23,*) x_vec
           write(23,*) x_vec_rot

           r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
           th(j,i) = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )
        
           arg1 = (x_vec_rot(1)  + smallval_dble) / &
                  ((x_vec_rot(1)**2 + x_vec_rot(2)**2)**.5 + smallval_dble)

           if (x_vec_rot(2) >= 0.) then
              ph(j,i) = acos(arg1)
           else
              ph(j,i) = 2. * pi - acos(arg1)
           end if
       enddo
    enddo 

end subroutine rotate_hetero
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_het_discr(rho,lambda,mu,rhopost,lambdapost,mupost,hetind)
    implicit none

    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda,mu
    double precision, dimension(0:npol,0:npol,nelem) :: rhopost,lambdapost,mupost
    integer :: hetind
    double precision, allocatable, dimension(:) :: disconttmp
    double precision :: w(4), wsum, dr(4)
    integer :: ndisctmp, ind(4), maxpts, iel, ipol, jpol, i, j
    integer, allocatable :: num_het_pts_region(:), het_ind(:, :)
    integer :: num_discr_het
    double precision :: s, z, r, th, r1, vptmp, vstmp, r2, r3, r4, th1, th2, th3, th4
    double precision, allocatable, dimension(:) :: rmin, rmax, thetamin, thetamax
    double precision, allocatable, dimension(:,:) :: shet, zhet
!#########################################################################################
! MvD: - rhet2 and r_het2 are quite confusing names..
!#########################################################################################

    double precision, allocatable :: rhet2(:,:), thhet2(:,:), phhet2(:,:) 

!#########################################################################################
! MvD: - stopping the code here until this is fixed
!#########################################################################################
    write(6,*) 'ERROR:'
    write(6,*) '   discrete input: interpolation wrong so far - work in progress...'
    stop

    write(6,*) mynum, 'reading discrete heterogeneity file...'

    open(unit=91, file=trim(het_file_discr(hetind)))

    read(91,*) num_discr_het
    
    if (lpr) write(6,*) 'Number of distinct-discrete regions:', num_discr_het

    allocate(num_het_pts_region(1:num_discr_het))

    do i=1,num_discr_het 
       read(91,*) num_het_pts_region(i)
       if (lpr) write(6,*) 'Region', i, 'has', num_het_pts_region(i), 'points.'
    enddo

    allocate(rmin(num_discr_het), rmax(num_discr_het), thetamin(num_discr_het), thetamax(num_discr_het))

    nhet_pts = sum(num_het_pts_region)
    maxpts = maxval(num_het_pts_region)

    if (lpr) write(6,*) 'max number of points in one region:', maxpts
    if (lpr) write(6,*) 'total number of points:', nhet_pts
    
    allocate(rhet2(1:maxpts,1:num_discr_het), thhet2(1:maxpts,1:num_discr_het), &
             phhet2(1:maxpts,1:num_discr_het))
    allocate(delta_vs2(1:maxpts,1:num_discr_het), delta_vp2(1:maxpts,1:num_discr_het), &
             delta_rho2(1:maxpts,1:num_discr_het))
    allocate(het_ind(1:maxpts,1:num_discr_het))

    rhet2 = -5000.
    thhet2 = -5000.
    phhet2 = -5000.
    delta_vp2 = -5000.
    delta_vs2 = -5000.
    delta_rho2 = -5000.
    
    write(6,*) mynum, 'read coordinates & medium properties...'
    do i=1, num_discr_het
       do j=1, num_het_pts_region(i)
          read(91,*) rhet2(j,i), thhet2(j,i), phhet2(j,i), delta_vp2(j,i), &
                     delta_vs2(j,i), delta_rho2(j,i)
          write(6,*) 'reading line', j, 'of', maxpts, ':', rhet2(j,i), &
                     thhet2(j,i), phhet2(j,i), delta_vp2(j,i), delta_vs2(j,i), &
                     delta_rho2(j,i)
       enddo
    enddo

    close(91)

    if (lpr) write(6,*) 'percent -> fraction'

    delta_vp2 = delta_vp2 / 100.
    delta_vs2 = delta_vs2 / 100.
    delta_rho2 = delta_rho2 / 100.
    thhet2 = thhet2 * pi / 180.
    phhet2 = phhet2 * pi / 180.
    rhet2 = rhet2 * 1000.

    ! Rotate coordinates if source is not on axis
    do i=1, num_discr_het
       rmin(i) = minval(rhet2(1:num_het_pts_region(i),i),1)
       rmax(i) = maxval(rhet2(1:num_het_pts_region(i),i),1)
       thetamin(i) = minval(thhet2(1:num_het_pts_region(i),i),1)
       thetamax(i) = maxval(thhet2(1:num_het_pts_region(i),i),1)

       write(6,*) mynum, 'r min/max:', i,rmin(i) / 1000., rmax(i) / 1000.
       write(6,*) mynum, 'th min/max:', i,thetamin(i) / pi * 180., thetamax(i) / pi * 180.

       if (rot_src) then 
          write(6,*) mynum, 'rotate since source is not beneath north pole'
          call rotate_hetero(num_het_pts_region(i),1,rhet2(1:num_het_pts_region(i),i), &
                             thhet2(1:num_het_pts_region(i),i), &
                             phhet2(1:num_het_pts_region(i),i))

          rmin(i) = minval(rhet2(1:num_het_pts_region(i),i),1)
          rmax(i) = maxval(rhet2(1:num_het_pts_region(i),i),1)
          thetamin(i) = minval(thhet2(1:num_het_pts_region(i),i),1)
          thetamax(i) = maxval(thhet2(1:num_het_pts_region(i),i),1)

          write(6,*) mynum, 'r min/max after rotation:', i, rmin(i) / 1000., &
                     rmax(i) / 1000.
          write(6,*) mynum, 'th min/max after rotation:', i, thetamin(i) / pi * 180., &
                     thetamax(i) / pi * 180.
       endif
    enddo

    ! plot discrete input file in vtk
    call plot_discrete_input(num_discr_het, num_het_pts_region, rhet2, thhet2, phhet2)

    ! for plotting discrete points within heterogeneous region
    rhetmin = minval(rmin,1)
    rhetmax = maxval(rmax,1)
    thhetmin = minval(thetamin,1)
    thhetmax = maxval(thetamax,1)

    write(6,*) 'r het min/max:', rhetmin / 1000., rhetmax / 1000.
    write(6,*) 'th het min/max:', thhetmin / pi * 180., thhetmax / pi * 180.

    ! revert to cylindrical 
    allocate (shet(1:maxpts,1:num_discr_het), zhet(1:maxpts,1:num_discr_het))
    shet = rhet2 * sin(thhet2) 
    zhet = rhet2 * cos(thhet2)

    write(6,*) mynum, 'locate GLL points within heterogeneous regions & '
    write(6,*) mynum, 'interpolate over 4 adjacent points'

    ind = 1
    wsum = 0
    w = 0

    do iel=1, nelem
       call compute_coordinates(s, z, r1, th1, iel, npol, npol)
       call compute_coordinates(s, z, r2, th2, iel, 0, 0)
       do i=1, num_discr_het
          r = max(r1,r2)
          th = max(th1,th2)
          if ( r>=rmin(i) .and. th>=thetamin(i) ) then
             r = min(r1,r2)
             th = min(th1,th2)
             if ( r<=rmax(i) .and. th<=thetamax(i) ) then
                do ipol=0,npol
                   do jpol=0,npol
                      ! find closest 4 points of discrete heterogeneous mesh
                      call compute_coordinates(s, z, r, th, iel, ipol, jpol)
                      call bilinear_interpolation(s, z, num_het_pts_region(i), &
                                                  shet(1:num_het_pts_region(i),i), &
                                                  zhet(1:num_het_pts_region(i),i), &
                                                  ind, w, wsum)
                      ! bilinear interpolationx
                      ! f(x,y) = [x2-x x-x1] [f(11) f(12); f(21) f(22)] [y2-y; y-y1] * &
                      !          1/(x2-x1)(y2-y1)
                      vptmp = sqrt((lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                                   rho(ipol,jpol,iel))
                      vstmp = sqrt(mu(ipol,jpol,iel) / rho(ipol,jpol,iel))
                      rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                                               (1. + sum(w * delta_rho2(ind,i)) * wsum)
                      vptmp = vptmp * (1. + sum(w * delta_vp2(ind,i)) * wsum)
                      vstmp = vstmp * (1. + sum(w * delta_vs2(ind,i)) * wsum)
                      mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
                      lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * &
                                                  (vptmp**2 - 2. * vstmp**2)
                   enddo
                enddo
             endif
          endif
       enddo
    enddo

    write(6,*) mynum, 'DONE loading discrete grid'
    
    deallocate(rhet2, thhet2, phhet2)
    deallocate(delta_vs2, delta_vp2, delta_rho2)
    deallocate(het_ind)
    deallocate (shet, zhet)

end subroutine load_het_discr
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine bilinear_interpolation(s0,z0,n,s,z,ind,w,wsum)
    implicit none

    ! implemented by fanie, find_closest4_points didnt really work for me...
    integer, intent(in) :: n
    double precision, intent(in) :: s0, z0, s(1:n), z(1:n)
    integer, intent(out) :: ind(4)
    double precision, intent(out) :: w(4), wsum
    integer :: i, p, is1, is2, is3, is4, iz1, iz2, iz3, iz4
    double precision :: dr(4), d2d(n,2), ds1, ds2, ds3, ds4

    ! choose 4 closest points, numbered after quadrant seen from s0,z0
    ! better: find minimum in array after calculating 2D-distance d2d in loop over n
    do i=1, n
        if ((s(i)-s0)>=0.0 .and. (z(i)-z0)>=0.0) then
            d2d(i,1) = sqrt((s(i)-s0)**2+(z(i)-z0)**2)
            d2d(i,2) = 1
            write(6,*) 'bilinttest', (s(i)-s0), (z(i)-z0), d2d(i,1), 1
        elseif ((s0-s(i))>=0.0 .and. (z(i)-z0)>=0.0) then
            d2d(i,1) = sqrt((s0-s(i))**2 + (z(i)-z0)**2)
            d2d(i,2) = 2
            write(6,*) 'bilinttest', (s(i)-s0), (z(i)-z0), d2d(i,1), 2
        elseif ((s0-s(i))>=0.0 .and. (z0-z(i))>=0.0) then
            d2d(i,1) = sqrt((s0-s(i))**2 + (z0-z(i))**2)
            d2d(i,2) = 3
            write(6,*) 'bilinttest', (s(i)-s0), (z(i)-z0), d2d(i,1), 3
        elseif ((s(i)-s0)>=0.0 .and. (z0-z(i))>=0.0) then
            d2d(i,1) = sqrt((s(i)-s0)**2+(z0-z(i))**2)
            d2d(i,2) = 4
            write(6,*) 'bilinttest', (s(i)-s0), (z(i)-z0), d2d(i,1), 4
        endif
    enddo

    ! find minimum of d2d for each quadrant
    ds1 = minval(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==1)
    ds2 = minval(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==2)
    ds3 = minval(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==3)
    ds4 = minval(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==4)
    is1 = minloc(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==1)
    is2 = minloc(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==2)
    is3 = minloc(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==3)
    is4 = minloc(d2d(1:n,1), DIM=1, MASK=(d2d(1:n,2))==4)

    dr(1) = ds1
    dr(2) = ds2
    dr(3) = ds3
    dr(4) = ds4

    ind(1) = is1
    ind(2) = is2
    ind(3) = is3
    ind(4) = is4

    write(6,*) 'bil int1', dr(1), ind(1)
    write(6,*) 'bil int2', dr(2), ind(2)
    write(6,*) 'bil int3', dr(3), ind(3)
    write(6,*) 'bil int4', dr(4), ind(4)

    ! inverse distance weighting
    do i=1, 4
       w(i) = (dr(i))**(-p)
    enddo

    wsum = 1. / sum(w)

    ! bilinear interpolation (implemented by fanie)
    ! f(x,y) = [x2-x x-x1] [f(11) f(12); f(21) f(22)] [y2-y; y-y1] * 1/(x2-x1)(y2-y1)

end subroutine bilinear_interpolation
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine plot_discrete_input(num_discr_het, num_het_pts_region, rhet2, thhet2, phhet2)

    use background_models, only : velocity
    use data_mesh, only : discont,bkgrdmodel

    implicit none

    integer, intent(in) :: num_het_pts_region(:), num_discr_het
    integer :: i,j,idom,icount
    real, allocatable, dimension(:,:) :: meshtmp
    real, allocatable, dimension(:) :: vptmp,vstmp,rhotmp
    character(len=80) :: fname
    double precision, intent(in) :: rhet2(:,:), thhet2(:,:), phhet2(:,:)

    allocate(vptmp(nhet_pts), vstmp(nhet_pts), rhotmp(nhet_pts), meshtmp(nhet_pts,2))
    icount = 0
    do i=1, num_discr_het
       do j=1, num_het_pts_region(i)
          icount = icount + 1
          idom = minloc(abs(discont-rhet2(j,i)),1)

!#########################################################################################
! MvD: - calling velocityi() here causes problems with anisotropy, anway it is
!        called already in get_model, so why not use the lame parameters here?
!#########################################################################################

          vptmp(icount) = velocity(rhet2(j,i), 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
          vptmp(icount) = vptmp(icount) * (1. + delta_vp2(j,i))

          vstmp(icount) = velocity(rhet2(j,i), 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
          vstmp(icount) = vstmp(icount) * (1. + delta_vs2(j,i))

          rhotmp(icount) = velocity(rhet2(j,i), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
          rhotmp(icount) = rhotmp(icount)* (1.+delta_rho2(j,i))

          meshtmp(icount,1) = rhet2(j,i) * sin(thhet2(j,i))
          meshtmp(icount,2) = rhet2(j,i) * cos(thhet2(j,i))
       enddo
    enddo

    fname = trim('Info/model_rho_discr_het'//appmynum)
    call write_VTK_bin_scal_pts(rhotmp(1:nhet_pts), meshtmp(1:nhet_pts,1:2), &
                                nhet_pts, fname)

    fname = trim('Info/model_vp_discr_het'//appmynum)
    call write_VTK_bin_scal_pts(vptmp(1:nhet_pts), meshtmp(1:nhet_pts,1:2), &
                                nhet_pts, fname)

    fname = trim('Info/model_vs_discr_het'//appmynum)
    call write_VTK_bin_scal_pts(vstmp(1:nhet_pts), meshtmp(1:nhet_pts,1:2), &
                                nhet_pts, fname)

    deallocate(meshtmp, vptmp, vstmp, rhotmp)

end subroutine plot_discrete_input
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
subroutine load_random(rho,lambda,mu,hetind)

    use commun
    use data_mesh, only : naxel, ax_el
    use utlity, only :  rcoord,zcoord

    implicit none 

    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda, mu
    double precision, dimension(0:npol,0:npol,nelem) :: rhopost, lambdapost, mupost
    integer :: hetind
    real(kind=8) :: t, decay, shift_fact, max_delta_vp, max_delta_vs, max_delta_rho
    real(kind=8) :: vptmp, vstmp, rhotmp, s, z, r, th, gauss_val
    integer :: iel, ipol, jpol, icount, i
    real(kind=8) :: rand
    real(kind=8), allocatable :: r_rad(:), rand_rad(:), r_radtmp(:), rand_radtmp(:)

    write(6,*)'add random anomalies to structure'
    
    !###################################################################### 
    ! MvD: Do we still need this?
    !###################################################################### 

    !!$! add randomly to each 2D point : laterally heterogeneous and same random
    !!                                   number to vp,vs,rho
    !!$do iel=1,nelem
    !!$   do jpol=0,npol
    !!$      do ipol=0,npol
    !!$         call random_number(rand)
    !!$         rand = 2.*rand-1.
    !!$         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
    !!$         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
    !!$         rho(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(1)*rand)
    !!$         vptmp = vptmp*(1. + delta_vp(1)*rand)
    !!$         vstmp = vstmp*(1. + delta_vs(1)*rand)
    !!$         lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) *( vptmp*vptmp - two*vstmp*vstmp )
    !!$         mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
    !!$      enddo
    !!$   enddo
    !!$enddo

    !!$! go along axis to find the radial profile
    !!$allocate(r_radtmp(naxel*(npol+1)), rand_radtmp(naxel*(npol+1)))
    !!$if (mynum==0) then 
    !!$icount=0
    !!$do iel=1,naxel
    !!$   do jpol=0,npol
    !!$      if (zcoord(0,jpol,ax_el(iel)) >=0.) then 
    !!$         icount=icount+1
    !!$         r_radtmp(icount) = rcoord(0,jpol,ax_el(iel))
    !!$         call random_number(rand)
    !!$         rand = 2.*rand-1.
    !!$         rand_radtmp(icount)=rand
    !!$      endif
    !!$   enddo
    !!$enddo
    !!$endif 
    !!$
    !!$! broadcast the profile to all processors
    !!$call broadcast_int(icount,0)
    !!$write(6,*)mynum,'number of radii:',icount
    !!$allocate(r_rad(icount),rand_rad(icount))
    !!$do i=1,icount
    !!$   call broadcast_dble(r_radtmp(i),0)
    !!$   r_rad(i)=r_radtmp(i)
    !!$   call broadcast_dble(rand_radtmp(i),0)
    !!$   rand_rad(i)=rand_radtmp(i)
    !!$enddo
    !!$
    !!$! add randomly to each radius, i.e. just altering the 1D background model
    !!$do iel=1,nelem
    !!$   do jpol=0,npol
    !!$      do ipol=0,npol
    !!$         i = minloc(abs(rcoord(ipol,jpol,iel)-r_rad(1:icount)),1)
    !!$         rand = rand_rad(i)
    !!$         vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
    !!$         vptmp = sqrt( (lambda(ipol,jpol,iel) + 2.*mu(ipol,jpol,iel)) / rho(ipol,jpol,iel) )
    !!$         rho(ipol,jpol,iel) = rho(ipol,jpol,iel)* (1. + delta_rho(1)*rand)
    !!$         vptmp = vptmp*(1. + delta_vp(1)*rand)
    !!$         vstmp = vstmp*(1. + delta_vs(1)*rand)
    !!$         lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) *( vptmp*vptmp - two*vstmp*vstmp )
    !!$         mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
    !!$      enddo
    !!$   enddo
    !!$enddo


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
       i = minloc(abs(rcoord(npol/2,npol/2,iel)-r_rad(1:icount)),1)
       rand = rand_rad(i)
       do jpol=0, npol
          do ipol=0, npol   
             vstmp = sqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
             vptmp = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel)) / &
                            rho(ipol,jpol,iel) )
             rhopost(ipol,jpol,iel) = rho(ipol,jpol,iel) * (1. + delta_rho(1) * rand)
             vptmp = vptmp * (1. + delta_vp(1) * rand)
             vstmp = vstmp * (1. + delta_vs(1) * rand)  
             lambdapost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * ( vptmp**2 - two*vstmp**2)
             mupost(ipol,jpol,iel) = rhopost(ipol,jpol,iel) * vstmp**2
          enddo
       enddo
    enddo

end subroutine load_random
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine load_het_funct(rho,lambda,mu,rhopost,lambdapost,mupost,hetind)

! added by fanie for sharp discontinuites
    use background_models, only : velocity
    use data_mesh, only : discont, bkgrdmodel

    implicit none 

    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: rho
    double precision, dimension(0:npol,0:npol,nelem), intent(inout) :: lambda, mu
    double precision, dimension(0:npol,0:npol,nelem) :: rhopost, lambdapost, mupost
    double precision :: t, decay, shift_fact, max_delta_vp, max_delta_vs, max_delta_rho
    double precision :: vptmp, vstmp, rhotmp, s, z, r, th, gauss_val
    double precision :: r_center_gauss, th_center_gauss
    double precision :: s_center_gauss, z_center_gauss, halfwidth_r, halfwidth_th
    double precision, allocatable :: rhet(:), thhet(:), phhet(:) 

    ! start elastic property values
    double precision :: vpst, vsst, rhost
    integer :: iel, ipol, jpol, icount, hetind, jj, ij, iel_count, idom
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
       rhetmin = minval(rhet(1:icount))
       rhetmax = maxval(rhet(1:icount))
       thhetmin = minval(thhet(1:icount))
       thhetmax = maxval(thhet(1:icount))

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
       rhetmin = minval(rhet(1:icount))
       rhetmax = maxval(rhet(1:icount))
       thhetmin = minval(thhet(1:icount))
       thhetmax = maxval(thhet(1:icount))

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
!#########################################################################################
! MvD : what does the invershape test here?
!#########################################################################################
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
       rhetmin = minval(rhet(1:icount))
       rhetmax = maxval(rhet(1:icount))
       thhetmin = minval(thhet(1:icount))
       thhetmax = maxval(thhet(1:icount))
       
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
subroutine plot_hetero_region_vtk(rho,lambda,mu)
    implicit none

    double precision, dimension(0:npol,0:npol,nelem), intent(in) :: rho 
    double precision, dimension(0:npol,0:npol,nelem), intent(in) :: lambda,mu
    real, dimension(:), allocatable :: vp_all,vs_all,rho_all
    real, dimension(:,:), allocatable :: mesh2
    character(len=80) :: fname
    double precision :: s, z, r, th
    integer :: iel, ipol, jpol, icount

    write(6,*) 'plotting heterogeneous region in pointwise vtk'

    allocate(mesh2(nelem * npol**2,2), vp_all(nelem * npol**2), vs_all(nelem * npol**2), &
             rho_all(nelem * npol**2))

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
                vp_all(icount) = sqrt( (lambda(ipol,jpol,iel) + 2. * mu(ipol,jpol,iel) ) / &
                                       rho(ipol,jpol,iel)  )
                vs_all(icount) = sqrt( (mu(ipol,jpol,iel) ) / rho(ipol,jpol,iel)  )
                rho_all(icount) = rho(ipol,jpol,iel) 
!#########################################################################################
! MvD: - file 666 + mynum is never opened, so how does this work?
!      - why the multiplication whith 0.2 and 0.3 ??
!#########################################################################################
                !write(666+mynum,20) r / 1000., th * 180. / pi, 0.0, -vp_all(icount) * 0.2, &
                !                    -vs_all(icount) * 0.4, rho(ipol,jpol,iel) * 0.3
             endif
          enddo
       enddo
    enddo

    20 format(3(1pe12.3), 3(1pe10.1))

    write(6,*) mynum, 'number of points inside heterogeneous region:', icount

    fname = trim('Info/model_vp_gll_het'//appmynum)
    call write_VTK_bin_scal_pts(real(vp_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                icount, fname)

    fname = trim('Info/model_vs_gll_het'//appmynum)
    call write_VTK_bin_scal_pts(real(vs_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                icount, fname)

    fname = trim('Info/model_rho_gll_het'//appmynum)
    call write_VTK_bin_scal_pts(real(rho_all(1:icount)), real(mesh2(1:icount,1:2)), &
                                icount, fname)

    deallocate(mesh2, vp_all, vs_all, rho_all)

end subroutine plot_hetero_region_vtk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal_pts(u2,mesh1,rows,filename)

   implicit none
   integer :: i
   integer, intent(in) :: rows
   real, dimension(1:rows), intent(in) :: u2
   real, dimension(1:rows) :: u1
   real, dimension(1:rows,1:2), intent(in) :: mesh1
   integer, dimension(1:rows*2) :: cell
   integer, dimension(1:rows) :: cell_type
   character (len=55) :: filename
   character (len=50) :: ss !stream
    
   !points structure
   do i=2, rows*2, 2
      cell(i-1) = 1
      cell(i) = (i / 2) - 1
   enddo

   do i=1, rows
      cell_type(i) = 1
   enddo
  
   u1 = real(u2)
   do i=1, rows
      if (abs(u1(i)) < 1.e-25) u1(i) = 0.0
   enddo
  
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
   write(ss,fmt='(A5,2I10)') 'CELLS',rows,rows*2
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

end subroutine write_vtk_bin_scal_pts
!-----------------------------------------------------------------------------

end module lateral_heterogeneities
