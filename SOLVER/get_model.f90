
!========================
module get_model
!========================

use global_parameters
use data_mesh
use data_io
use data_proc

use background_models
use utlity

implicit none

public :: read_model, read_model_ani, model_output
private
contains
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine read_model(rho,lambda,mu)
!
! First define array ieldom that specifically appoints a respective domain 
! between discontinuities for each element (to avoid issues very close to 
! discontinuities).
! Fill each GLL(GLJ) grid point with rho,lambda,mu on domain basis for given 
! background model as defined in background_models.f90 (NOTE: that routine 
! is shared between solver and mesher, do not edit for one without the other!).
! Compute global min/max values.
!
! Runs various checks and outputs related to the background model:
!
! 1) Discontinuities (given by discont(ndisc) from the mesher):
!    Count points to check on conformity with coarsening and total number of 
!    lateral GLL points, check sphericity and lateral elastic homogeneity, 
!    output velocities just above/below
! 2) Elastic background model:
!    check whether all GLL points assumed reasonable values, whether all 
!    lateral neighbor points share the same elastic properties and 
!    compute maximal variation across an element (issue warning if large)
! 3) Resolution test:
!    Define radial trigonometric functions with wavelengths corresponding 
!    to the source period, or half, or double of it.
!    Compute accuracy of same for each element, using max/min of vp/vs 
!    elementally and globally.
!    This test gives a somewhat reasonable handle on how well the desired 
!    source frequencies are resolved at worst, globally and locally.
!    (of course only spatially: time discretization is usually worse such that 
!     these values should not be taken absolutely, but rather as a relative 
!     heuristic feel for the grid-model melange).
!    Also outputs various forms of this test to files (called resolutionsine).
! 4) Numerical resolution: 
!    computes the characteristic lead time (grid spacing)/velocity which 
!    is the determining factor for the numerical performance, stability 
!    and mesh efficiency. Outputs various descendants such as time step and 
!    source period as a function of radius, and separately writes min/max 
!    characteristic lead times for P- and S-wave velocities of all elements 
!    per processor, and for critical subdomains such as the central cube 
!    and coarsening layers. 
!    On-the-fly verification of respective radial averages:
!       xmgrace timestep_rad.dat or xmgrace period_rad.dat
!    
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use commun, ONLY : barrier
use lateral_heterogeneities
use nc_routines, ONLY: nc_write_el_domains
include 'mesh_params.h'


double precision, dimension(0:npol,0:npol,nelem), intent(out) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(out) :: lambda,mu
double precision :: s,z,r,theta,r1,vptmp,vstmp,r2,r3,r4,th1,th2,th3,th4
integer :: iel,ipol,jpol,iidom,ieldom(nelem),domcount(ndisc),iel_count,ij,jj
logical :: foundit
!af pgf
character(len=100) :: modelstring
double precision, dimension(0:npol,0:npol,nelem) :: stmp,ztmp,v_p,v_s

  if (make_homo ) then 
     if (.not. have_fluid) then
        if (lpr) then
           write(6,*)'  '
           write(6,*)'  '
           write(6,*)'******************** SERIOUS MESSAGE ************************'
           write(6,*)'**                                                         **'
           write(6,*)'**   Overriding background model with HOMOGENEOUS SOLID!   **'
           write(6,*)'**                                                         **'
           write(6,*)'*************************************************************'
           write(6,*)'  '
           write(6,*)'Vp  = ', vphomo
           write(6,*)'Vs  = ', vshomo
           write(6,*)'rho = ', rhohomo
           write(6,*)'  '
        endif
     else
        if (lpr) then
           write(6,*)
           write(6,*)'PROBLEM: Want to override with entirely homogeneous model.'
           write(6,*)'      ...but that is only possible if there is no fluid.'
           write(6,*)'===> Re-run the mesher for a model without a fluid layer'
           write(6,*)'     or change make_homo to false in input parameter file.'
           write(6,*)'NEW ADDITION (NOV 2010):'
           write(6,*)'Keeping the fluid, putting vphomo but vs=0'
           write(6,*)'NOTE: This is obviously NOT homogeneous, as core reflections will occur.'
           write(6,*)
        endif
        !stop
     endif
  endif

! Set elastic parameters to crazy values to later check if all have been filled
  rho(0:npol,0:npol,1:nelem) = -1.E30
  lambda(0:npol,0:npol,1:nelem) = -1.E30
  mu(0:npol,0:npol,1:nelem) = -1.E30

! check for each element which domain it belongs to
  open(unit=65,file=infopath(1:lfinfo)//'/elems_bkgrdmodel_domain.dat'&
                                        //appmynum)

  if (do_mesh_tests) open(60000+mynum,file='Data/model_r_th_rho_vp_vs_'//appmynum//'.dat')

  do iel=1,nelem
     foundit=.false.
     r1=rcoord(int(npol/2),int(npol/2),iel) 
     do iidom=1,ndisc-1
        if (r1<discont(iidom) .and. r1> discont(iidom+1) ) then
           ieldom(iel)=iidom
           foundit=.true.
           write(65,10)iel,r1,iidom,discont(iidom),discont(iidom+1)
        endif
     enddo
     if (r1 < discont(ndisc)) then
        ieldom(iel)=ndisc
        foundit=.true.
        write(65,10)iel,r1,ndisc,discont(ndisc)
     endif
     if (.not. foundit) then 
        write(6,*)'havent found domain for element',iel
        write(6,*)'...of radius',r1
        stop
     endif
  enddo
  close(65)
10 format(i9,1pe11.3,i3,2(1pe11.3))
  
  if (use_netcdf) then
     if (mynum.eq.0) call nc_write_el_domains(ieldom)
  end if


  if (do_mesh_tests) then
    if (lpr) write(6,*)'    checking discontinuity discretization...' 
    call check_mesh_discontinuities(ieldom,domcount)
  endif

! read in respective velocities and density on domain basis
  open(unit=5454,file=infopath(1:lfinfo)//'/background_rad_dom_vel.dat'&
                                          //appmynum)

  if (lpr) write(6,*)'    filling mesh with elastic properties...'   
  call flush(6)

! solar case: interpolations take too long, therefore reading in model parameters here
! construct full arrays for velocities... faster! 
if (bkgrdmodel=='solar') then 
   write(6,*)'pre-assembling media arrays for solar case.... '
   write(6,*)'...much faster due to inperpolation routine!'
   do iel=1,nelem
     do ipol=0,npol
        do jpol=0,npol
           call compute_coordinates(stmp(ipol,jpol,iel),ztmp(ipol,jpol,iel),r,theta,iel,ipol,jpol)
        enddo
     enddo
  enddo
  call arbitr_sub_solar_arr(stmp*router,ztmp*router,v_p,v_s,rho,bkgrdmodel)
  write(6,*)'done with media array definition '
  write(6,*)'minmax vp:',minval(v_p),maxval(v_p),minval(rho),maxval(rho)
endif

!af
  modelstring=bkgrdmodel!(1:(index(bkgrdmodel,' ')-1))
  iel_count=0
!========================
  do iel=1,nelem
!========================

     iidom=ieldom(iel)
     do ipol=0,npol
        do jpol=0,npol

           if (bkgrdmodel/='solar') then

               call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
               vptmp=velocity(r,'v_p',iidom,modelstring,lfbkgrdmodel)
               vstmp=velocity(r,'v_s',iidom,modelstring,lfbkgrdmodel)
               rho(ipol,jpol,iel)= &
                    velocity(r,'rho',iidom,modelstring,lfbkgrdmodel)
               
               if (make_homo) then
                  if (vstmp>0.1) then 
                     vptmp = vphomo; vstmp = vshomo
                     rho(ipol,jpol,iel) = rhohomo
                  else ! in fluid, impose same vp but vs as zero
                     vptmp = vphomo; vstmp = 0.0
                     rho(ipol,jpol,iel) = rhohomo              
                  endif
               endif
           
! HOMOGENEOUS SPHERE, & FLUID CORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         if (iidom == ndisc-1) vstmp = 0.d0
!
!! Keep impedance constant, but vary vp and density respectively across ICB/CMB
!!         if (iidom == ndisc-1) vptmp = 8.d3
!!         if (iidom == ndisc-1) rho = 5.d3

           else ! solar
              vstmp=0.
              vptmp=v_p(ipol,jpol,iel)
           endif

           lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) * &
                                   ( vptmp*vptmp - two*vstmp*vstmp )
           mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vstmp*vstmp
           
           if (save_large_tests) &
                write(5454,12)r,iidom,vptmp,vstmp,rho(ipol,jpol,iel)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! test round-off errors for elastic parameters & velocities
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
           if ( .not. dblreldiff_small(vptmp,dsqrt( (lambda(ipol,jpol,iel) + &
                two*mu(ipol,jpol,iel) )/rho(ipol,jpol,iel) ) ) ) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and p-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vptmp:',&
                   r/1000.,theta*180./pi,vptmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
           
           if ( .not. dblreldiff_small(vstmp,dsqrt( mu(ipol,jpol,iel)/ &
                rho(ipol,jpol,iel) ))) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and s-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vstmp:',&
                   r/1000.,theta*180./pi,vstmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
           
        enddo
     enddo

! write out for later snaps
if (do_mesh_tests) then
do ipol=ibeg,iend
        do jpol=ibeg,iend
                   call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
        write(60000+mynum,14)r,theta,rho(ipol,jpol,iel),&
             sqrt( (lambda(ipol,jpol,iel)+2.*mu(ipol,jpol,iel))/rho(ipol,jpol,iel)), &
             sqrt( mu(ipol,jpol,iel)/rho(ipol,jpol,iel))
     enddo
  enddo
  endif

!========================
  enddo ! nelem
!========================
  close(5454)
write(6,*)mynum,'done with big mesh loop to define model'
 if (do_mesh_tests) close(60000+mynum)
12 format(1pe15.7,i4,3(1pe15.7))
14 format(5(1pe13.4))

! TNM Oct 29: addition of heterogeneities in separate routine
if (add_hetero) call compute_heterogeneities(rho,lambda,mu)

! TNM Oct 28: plot final velocity model in vtk
 write(6,*)mynum,'plotting vtks for the model properties....'
 call plot_model_vtk(rho,lambda,mu)

!@@@@@@@@@@@@@@@@@
! Some tests....
!@@@@@@@@@@@@@@@@@
if (do_mesh_tests) then
  call barrier

  if (.not. add_hetero) then
     if (lpr) write(6,*)'    checking elastic properties on discontinuities...' 
     call check_elastic_discontinuities(ieldom,domcount,lambda,mu,rho)
  else
        if (lpr) write(6,*)'    NOT checking elastic properties on discontinuities since we added hetero...'   
  endif

  if (lpr) write(6,*)'    checking the background model discretization...' 
  call check_background_model(lambda,mu,rho)

  if (lpr) write(6,*)'    testing the mesh/background model resolution...' 
  call test_mesh_model_resolution(lambda,mu,rho)
endif 

! Compute time step, period, courant, points per wavelength throughout the mesh
!af
! if (lpr) write(6,*)'    computing grid spacing vs. velocities, time step...' 
! call compute_numerical_resolution(lambda,mu,rho)

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! compute min/max velocities in whole domain
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  vpmin=minval(dsqrt((lambda+two*mu)/rho))
  vpminloc=minloc(dsqrt((lambda+two*mu)/rho))
  vpmax=maxval(dsqrt((lambda+two*mu)/rho))
  vpmaxloc=maxloc(dsqrt((lambda+two*mu)/rho))
  
  vsmin=minval(dsqrt(mu/rho))
  vsminloc=minloc(dsqrt(mu/rho))
  vsmax=maxval(dsqrt(mu/rho))
  vsmaxloc=maxloc(dsqrt(mu/rho))
  
! since minloc/maxloc start counting at 1...
  vpminloc(1)=vpminloc(1)-1; vpminloc(2)=vpminloc(2)-1;
  vsminloc(1)=vsminloc(1)-1; vsminloc(2)=vsminloc(2)-1;
  vpmaxloc(1)=vpmaxloc(1)-1; vpmaxloc(2)=vpmaxloc(2)-1;
  vsmaxloc(1)=vsmaxloc(1)-1; vsmaxloc(2)=vsmaxloc(2)-1;
  
  call compute_coordinates(s,z,vpminr,theta,vpminloc(3),&
                           vpminloc(1),vpminloc(2))
  call compute_coordinates(s,z,vsminr,theta,vsminloc(3),&
                           vsminloc(1),vsminloc(2))
  call compute_coordinates(s,z,vpmaxr,theta,vpmaxloc(3),&
                           vpmaxloc(1),vpmaxloc(2))
  call compute_coordinates(s,z,vsmaxr,theta,vsmaxloc(3),&
                           vsmaxloc(1),vsmaxloc(2))

end subroutine read_model
!=============================================================================


!-----------------------------------------------------------------------------
subroutine read_model_ani(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                          fa_ani_theta, fa_ani_phi)
!
! same as above, but for anisotropic models
!    
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use commun, ONLY : barrier
use lateral_heterogeneities
use data_mesh_preloop, ONLY : ielsolid
use data_io, ONLY : rot_mat
use data_source, ONLY : rot_src

include 'mesh_params.h'

double precision, dimension(0:npol,0:npol,nelem), intent(out) :: rho
double precision, dimension(0:npol,0:npol,nelem), intent(out) :: lambda, mu
double precision, dimension(0:npol,0:npol,nelem), intent(out) :: xi_ani, phi_ani, eta_ani
double precision, dimension(0:npol,0:npol,nelem), intent(out) :: fa_ani_theta, fa_ani_phi
double precision :: s,z,r,theta,r1,r2,r3,r4,th1,th2,th3,th4
double precision :: vphtmp, vpvtmp, vshtmp, vsvtmp
integer :: iel,ipol,jpol,iidom,ieldom(nelem),domcount(ndisc),iel_count,ij,jj
logical :: foundit
character(len=100) :: modelstring

  if (make_homo ) then 
     write(6,*)'  '
     write(6,*)'ERROR: homogeneous AND anisotropic model does not make '
     write(6,*)'       sense, check input file'
     write(6,*)'  '
     stop
  endif

! Set elastic parameters to crazy values to later check if all have been filled
  rho(0:npol,0:npol,1:nelem) = -1.E30
  lambda(0:npol,0:npol,1:nelem) = -1.E30
  mu(0:npol,0:npol,1:nelem) = -1.E30
  xi_ani(0:npol,0:npol,1:nelem) = -1.E30
  phi_ani(0:npol,0:npol,1:nelem) = -1.E30
  eta_ani(0:npol,0:npol,1:nelem) = -1.E30

! check for each element which domain it belongs to
  open(unit=65,file=infopath(1:lfinfo)//'/elems_bkgrdmodel_domain.dat'&
                                       //appmynum)

  if (do_mesh_tests) open(60000+mynum,file='Data/model_r_th_rho_vp_vs_'&
                                          //appmynum//'.dat')

  do iel=1,nelem
     foundit=.false.
     r1 = rcoord(int(npol/2),int(npol/2),iel) 
     do iidom=1,ndisc-1
        if (r1<discont(iidom) .and. r1> discont(iidom+1)) then
           ieldom(iel) = iidom
           foundit = .true.
           write(65,10)iel,r1,iidom,discont(iidom),discont(iidom+1)
        endif
     enddo
     if (r1 < discont(ndisc)) then
        ieldom(iel) = ndisc
        foundit = .true.
        write(65,10)iel,r1,ndisc,discont(ndisc)
     endif
     if (.not. foundit) then 
        write(6,*)'havent found domain for element',iel
        write(6,*)'...of radius',r1
        stop
     endif
  enddo
  close(65)

10 format(i9,1pe11.3,i3,2(1pe11.3))

  if (do_mesh_tests) then
    if (lpr) write(6,*)'    checking discontinuity discretization...' 
    call check_mesh_discontinuities(ieldom,domcount)
  endif

! read in respective velocities and density on domain basis
  open(unit=5454,file=infopath(1:lfinfo)//'/background_rad_dom_vel.dat'&
                                          //appmynum)

  if (lpr) write(6,*)'    filling mesh with elastic properties...'   
  call flush(6)

  modelstring = bkgrdmodel


  iel_count=0
!========================
  do iel=1,nelem
!========================

     iidom=ieldom(iel)
     do ipol=0,npol
        do jpol=0,npol
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)

           vphtmp = velocity(r,'vph',iidom,modelstring,lfbkgrdmodel)
           vpvtmp = velocity(r,'vpv',iidom,modelstring,lfbkgrdmodel)
           vshtmp = velocity(r,'vsh',iidom,modelstring,lfbkgrdmodel)
           vsvtmp = velocity(r,'vsv',iidom,modelstring,lfbkgrdmodel)
           eta_ani(ipol,jpol,iel) = velocity(r,'eta',iidom,modelstring,lfbkgrdmodel)
           rho(ipol,jpol,iel) = velocity(r,'rho',iidom,modelstring,lfbkgrdmodel)
           
           lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) * (vphtmp**2 - two*vshtmp**2)
           mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vshtmp**2

           if (vsvtmp > (smallval_sngl * vphtmp)) then
              xi_ani(ipol,jpol,iel) = vshtmp**2 / vsvtmp**2
           else
              xi_ani(ipol,jpol,iel) = one
           endif
           
           phi_ani(ipol,jpol,iel) = vpvtmp**2 / vphtmp**2
            
           ! radial anisotropy (otherwise lateral heterogeneity!)
           fa_ani_theta(ipol,jpol,iel) = thetacoord(ipol, jpol, iel)
           fa_ani_phi(ipol,jpol,iel) = 0.

           if (save_large_tests) &
                write(5454,12) r, iidom, vphtmp, vpvtmp, vshtmp, vsvtmp, &
                        eta_ani(ipol,jpol,iel), rho(ipol,jpol,iel)

        ! XXX generalize tests for anisotropic parameters:

        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! test round-off errors for elastic parameters & velocities
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
           if ( .not. dblreldiff_small(vphtmp,dsqrt( (lambda(ipol,jpol,iel) + &
                two*mu(ipol,jpol,iel) )/rho(ipol,jpol,iel) ) ) ) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and p-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vphtmp:',&
                   r/1000.,theta*180./pi,vphtmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
           
           if ( .not. dblreldiff_small(vshtmp,dsqrt( mu(ipol,jpol,iel)/ &
                rho(ipol,jpol,iel) ))) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and s-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vshtmp:',&
                   r/1000.,theta*180./pi,vshtmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
           
        enddo
     enddo

     ! XXX generalize tests for anisotropic parameters:
     ! write out for later snaps
     if (do_mesh_tests) then
       do ipol=ibeg,iend
          do jpol=ibeg,iend
             call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
             write(60000+mynum,14)r,theta,rho(ipol,jpol,iel),&
                  sqrt( (lambda(ipol,jpol,iel)+2.*mu(ipol,jpol,iel))/rho(ipol,jpol,iel)), &
                  sqrt( mu(ipol,jpol,iel)/rho(ipol,jpol,iel))
          enddo
       enddo
     endif

!========================
  enddo ! nelem
!========================
  close(5454)
  write(6,*)mynum,'done with big mesh loop to define model'
  if (do_mesh_tests) close(60000+mynum)

12 format(1pe15.7,i4,6(1pe15.7))
14 format(5(1pe13.4))

  if (add_hetero) call compute_heterogeneities(rho, lambda, mu, xi_ani, phi_ani, &
                         eta_ani, fa_ani_theta, fa_ani_phi, ieldom)

! plot final velocity model in vtk
  write(6,*)mynum,'plotting vtks for the model properties....'

  call plot_model_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, fa_ani_theta, fa_ani_phi)

!@@@@@@@@@@@@@@@@@
! Some tests....
!@@@@@@@@@@@@@@@@
if (do_mesh_tests) then
  call barrier

  if (.not. add_hetero) then
     if (lpr) write(6,*)'    checking elastic properties on discontinuities...' 
     call check_elastic_discontinuities(ieldom,domcount,lambda,mu,rho)
  else
        if (lpr) write(6,*)'    NOT checking elastic properties on discontinuities since we added hetero...'   
  endif

  if (lpr) write(6,*)'    checking the background model discretization...' 
  call check_background_model(lambda,mu,rho)

  if (lpr) write(6,*)'    testing the mesh/background model resolution...' 
  call test_mesh_model_resolution(lambda,mu,rho)
endif 

! Compute time step, period, courant, points per wavelength throughout the mesh
!af
! if (lpr) write(6,*)'    computing grid spacing vs. velocities, time step...' 
! call compute_numerical_resolution(lambda,mu,rho)


!MvD: Do we need the following for each anisotropic velocity can we leave it as is? Do we need it at all? 
! - 2.2.12: just outputted to standard out in parameters.f90

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! compute min/max velocities in whole domain
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  vpmin=minval(dsqrt((lambda+two*mu)/rho))
  vpminloc=minloc(dsqrt((lambda+two*mu)/rho))
  vpmax=maxval(dsqrt((lambda+two*mu)/rho))
  vpmaxloc=maxloc(dsqrt((lambda+two*mu)/rho))
  
  vsmin=minval(dsqrt(mu/rho))
  vsminloc=minloc(dsqrt(mu/rho))
  vsmax=maxval(dsqrt(mu/rho))
  vsmaxloc=maxloc(dsqrt(mu/rho))
  
! since minloc/maxloc start counting at 1...
  vpminloc(1)=vpminloc(1)-1; vpminloc(2)=vpminloc(2)-1;
  vsminloc(1)=vsminloc(1)-1; vsminloc(2)=vsminloc(2)-1;
  vpmaxloc(1)=vpmaxloc(1)-1; vpmaxloc(2)=vpmaxloc(2)-1;
  vsmaxloc(1)=vsmaxloc(1)-1; vsmaxloc(2)=vsmaxloc(2)-1;
  
  call compute_coordinates(s,z,vpminr,theta,vpminloc(3),&
                           vpminloc(1),vpminloc(2))
  call compute_coordinates(s,z,vsminr,theta,vsminloc(3),&
                           vsminloc(1),vsminloc(2))
  call compute_coordinates(s,z,vpmaxr,theta,vpmaxloc(3),&
                           vpmaxloc(1),vpmaxloc(2))
  call compute_coordinates(s,z,vsmaxr,theta,vsmaxloc(3),&
                           vsmaxloc(1),vsmaxloc(2))

end subroutine read_model_ani
!=============================================================================


!-----------------------------------------------------------------------------
subroutine arbitr_sub_solar_arr(s,z,v_p,v_s,rho,bkgrdmodel2)
!
! file-based, step-wise model in terms of domains separated by disconts.
! format:
! ndisc
! r vp vs rho
! ...
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double precision, intent(in) :: s(0:npol,0:npol,1:nelem),z(0:npol,0:npol,1:nelem)
character(len=100), intent(in) :: bkgrdmodel2
double precision, dimension(:,:,:), intent(out) :: rho(0:npol,0:npol,1:nelem)
double precision, dimension(:,:,:), intent(out) :: v_s(0:npol,0:npol,1:nelem)
double precision, dimension(:,:,:), intent(out) :: v_p(0:npol,0:npol,1:nelem)
double precision, allocatable, dimension(:) :: disconttmp,rhotmp,vstmp,vptmp
integer :: ndisctmp,i,ndisctmp2,ind(2),ipol,jpol,iel
logical :: bkgrdmodelfile_exists
double precision :: w(2),wsum,r0

! Does the file bkgrdmodel".bm" exist?
  inquire(file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm', &
          exist=bkgrdmodelfile_exists)
  if (bkgrdmodelfile_exists) then
     open(unit=77,file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm')
     read(77,*)ndisctmp
     allocate(disconttmp(1:ndisctmp))
     allocate(vptmp(1:ndisctmp),vstmp(1:ndisctmp),rhotmp(1:ndisctmp))
     do i=1, ndisctmp
        read(77,*)disconttmp(i),rhotmp(i),vptmp(i),vstmp(i)
     enddo
     close(77)
     do iel=1,nelem
        do jpol=0,npol
           do ipol=0,npol
              r0 = dsqrt(s(ipol,jpol,iel)**2 +z(ipol,jpol,iel)**2 )
              call interp_vel(r0,disconttmp(1:ndisctmp),ndisctmp,ind,w,wsum)
              rho(ipol,jpol,iel)=sum(w*rhotmp(ind))*wsum
              v_p(ipol,jpol,iel)=(w(1)*vptmp(ind(1))+w(2)*vptmp(ind(2)))*wsum
              v_s(ipol,jpol,iel)=sum(w*vstmp(ind))*wsum
           enddo
        enddo
     enddo
     deallocate(disconttmp,vstmp,vptmp,rhotmp)
  else 
     write(6,*)'Background model file', &
          bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm','does not exist!!!'
     stop
  endif

end subroutine arbitr_sub_solar_arr
!=============================================================================

!=============================================================================
subroutine interp_vel(r0,r,n,ind,w,wsum)

integer, intent(in) :: n
double precision, intent(in) :: r0,r(1:n)
integer, intent(out) :: ind(2)
double precision, intent(out) :: w(2),wsum
integer :: i,p
double precision :: dr1,dr2

p=1

     i=minloc(dabs(r-r0),1)
!     write(6,*)'INTERP;',r0,i
!     write(6,*)'INTERP:',r(i)

     if (r0>0.d0) then
        if ((r(i)-r0)/r0> 1.d-8) then ! closest discont. at larger radius
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))
        elseif ((r0-r(i))/r0> 1.d-8) then  ! closest discont. at smaller radius
           if (r0>maxval(r)) then ! for round-off errors where mesh is above surface
              ind(1)=i
              ind(2)=i
              dr1=1.d0
              dr2=1.d0
           else
              ind(1)=i-1
              ind(2)=i
              dr1=r(ind(1))-r0
              dr2=r0-r(ind(2))
            endif
        elseif (dabs((r(i)-r0)/r0)< 1.d-8) then ! closest discont identical
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           write(6,*)'problem with round-off errors in interpolating......'
           write(6,*)'r0,r(i),i',r0,r(i),abs((r(i)-r0)/r0),i
           stop
        endif
     else !r0=0
        if (r(i)==0.d0) then ! center of the sun
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))        
        endif
     endif

! inverse distance weighting
     w(1)=(dr1)**(-p)
     w(2)=(dr2)**(-p)
     wsum=1.d0/sum(w)

end subroutine interp_vel
!=============================================================================



!-----------------------------------------------------------------------------
subroutine check_mesh_discontinuities(ieldom,domcount)

include 'mesh_params.h'

integer, intent(in)  :: ieldom(:)
integer, intent(out) :: domcount(ndisc)
integer              :: iel,ipol,jpol,iidom,eldomcount
double precision     :: s,z,r,theta

  domcount(:)=0
! Count points on discontinuities and check whether right domains are assigned
  do iel=1,nelem
     iidom=ieldom(iel)
     do jpol=0,npol
        eldomcount=-1
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
           if (r>zero) then 
              if ( (r-discont(iidom))/r > dble(smallval) ) then 
                 write(6,*)procstrg,'PROBLEM with domains and discontinuities!'
                 write(6,*)procstrg,'radius > associated discont.:',&
                      iidom,r,discont(iidom)
                 stop
              endif
           endif

           if ( dblreldiff_small(r,discont(iidom)) ) then
              domcount(iidom)=domcount(iidom)+1
              eldomcount=eldomcount+1
           endif
        enddo
        
! Check number of GLL points on discontinuities
        if (dblreldiff_small(r,discont(iidom)) .and. (eldomcount /= npol)) then
           write(6,*)procstrg,&
                'PROBLEM: not every GLL along xi is on the discont!'
           write(6,*)procstrg,'iel,idom,r,theta:',iel,ieldom(iel), &
                r/1000.d0,theta*180.d0/pi
           write(6,*)procstrg,'ipol count, npol:',eldomcount,npol
           stop
        endif
        
     enddo
  enddo
  
! make sure on each discontinuity we have a fraction of # surface points
! (in agreement with conformal coarsening/size doubling)

do iidom=2,ndisc
!     write(6,*)procstrg,'Domain count:',iidom,domcount(1),domcount(iidom)
     if (mod(domcount(1),domcount(iidom)) /= 0 ) then
        write(6,*)' ' 
        write(6,*)procstrg,&
             'PR0BLEM: # points on discontinuity not fraction of surface '
        write(6,*)procstrg,'Domain, discontinuity:',iidom, discont(iidom)
        write(6,*)procstrg,'number of points at discont, surface:', &
             domcount(iidom),domcount(1)
        stop
     endif
  enddo
  
! same as above, but make sure each level is fraction of its above neighbor
  do iidom=2,ndisc
     if (domcount(iidom-1) /= domcount(iidom)) then
        if (mod(domcount(iidom-1),2*domcount(iidom)) /= 0 ) then
           write(6,*) ' '
           write(6,*)procstrg,'PR0BLEM: # points discont. not even fraction', &
                'of discont. above'
           write(6,*)procstrg,'Domain, discontinuity:',iidom, discont(iidom)
           write(6,*)procstrg,'#points at discont above,here:', &
                domcount(iidom-1),domcount(iidom)
           stop
        endif
     endif
  enddo

end subroutine check_mesh_discontinuities
!=============================================================================

!-----------------------------------------------------------------------------
subroutine check_elastic_discontinuities(ieldom,domcount,lambda,mu,rho)

include 'mesh_params.h'

integer, intent(in)          :: ieldom(nelem),domcount(ndisc)
double precision, intent(in) :: rho(0:npol,0:npol,nelem)
double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
double precision, intent(in) :: mu(0:npol,0:npol,nelem)
integer                      :: iel,jpol,iidom
double precision             :: minvpdomabove(ndisc),maxvpdomabove(ndisc)
double precision             :: minvsdomabove(ndisc),maxvsdomabove(ndisc)
double precision             :: minrodomabove(ndisc),maxrodomabove(ndisc)
double precision             :: minvpdombelow(ndisc),maxvpdombelow(ndisc)
double precision             :: minvsdombelow(ndisc),maxvsdombelow(ndisc)
double precision             :: minrodombelow(ndisc),maxrodombelow(ndisc)
double precision             :: s,z,r,theta

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! count points on discontinuities and check whether right domains are assigned
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  open(unit=6464,file=infopath(1:lfinfo)//&
       '/velocities_below_discont.dat'//appmynum)
  open(unit=6465,file=infopath(1:lfinfo)//&
       '/velocities_above_discont.dat'//appmynum)
  
  minvpdombelow=1.d6; maxvpdombelow=0.d0; minvpdomabove=1.d6;maxvpdomabove=0.d0
  minvsdombelow=1.d6; maxvsdombelow=0.d0; minvsdomabove=1.d6;maxvsdomabove=0.d0
  minrodombelow=1.d6; maxrodombelow=0.d0; minrodomabove=1.d6;maxrodomabove=0.d0
  
  do iel=1,nelem
     iidom=ieldom(iel)
     do jpol=0,npol
        
        ! write out velocities for each element that shares a discontinuity
        call compute_coordinates(s,z,r,theta,iel,int(npol/2),jpol)
        if ( dblreldiff_small(r,discont(iidom)) ) then 
           
           if (minvpdombelow(iidom)> sqrt((lambda(int(npol/2),jpol,iel)+ &
                two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))) &
                minvpdombelow(iidom)=sqrt((lambda(int(npol/2),jpol,iel)+ &
                two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))
           
           if (maxvpdombelow(iidom)< sqrt((lambda(int(npol/2),jpol,iel)+ &
                two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))) & 
                maxvpdombelow(iidom)=sqrt((lambda(int(npol/2),jpol,iel)+ &
                two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) 
           
           if (minvsdombelow(iidom)>& 
                sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))) & 
                minvsdombelow(iidom)=&
                sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
           
           if (maxvsdombelow(iidom)<&
                sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) & 
                maxvsdombelow(iidom)=&
                sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
           
           if (minrodombelow(iidom)>rho(int(npol/2),jpol,iel)) & 
                minrodombelow(iidom)=rho(int(npol/2),jpol,iel)
           
           if (maxrodombelow(iidom)<rho(int(npol/2),jpol,iel) ) & 
                maxrodombelow(iidom)=rho(int(npol/2),jpol,iel)
           
! r,theta,vp,vs,rho for one point of each element on the discontinuity
           write(6464,14)r/1000.d0,theta*180.d0/pi, &
                sqrt((lambda(int(npol/2),jpol,iel)+ &
                two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)),&
                sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)),&
                rho(int(npol/2),jpol,iel)
        endif
        
        if (iidom<ndisc) then 
           if ( dblreldiff_small(r,discont(iidom+1)) ) then 
! r,theta,vp,vs,rho for one point of each element on the discontinuity
              
              if (minvpdomabove(iidom+1)> sqrt((lambda(int(npol/2),jpol,iel)+ &
                   two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))) &
                   minvpdomabove(iidom+1)=sqrt((lambda(int(npol/2),jpol,iel)+ &
                   two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))
              
              if (maxvpdomabove(iidom+1)< sqrt((lambda(int(npol/2),jpol,iel)+ &
                   two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))) & 
                   maxvpdomabove(iidom+1)=sqrt((lambda(int(npol/2),jpol,iel)+ &
                   two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) 
              
              if (minvsdomabove(iidom+1)>& 
                   sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))) & 
                   minvsdomabove(iidom+1)=&
                   sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
              
              if (maxvsdomabove(iidom+1)<&
                   sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) &
                   maxvsdomabove(iidom+1)=&
                   sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
              
              if (minrodomabove(iidom+1)>rho(int(npol/2),jpol,iel)) & 
                   minrodomabove(iidom+1)=rho(int(npol/2),jpol,iel)
              
              if (maxrodomabove(iidom+1)<rho(int(npol/2),jpol,iel) ) & 
                   maxrodomabove(iidom+1)=rho(int(npol/2),jpol,iel)
              
              write(6465,14)r/1000.d0,theta*180.d0/pi, &
                   sqrt((lambda(int(npol/2),jpol,iel)+ &
                   two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)),&
                   sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)),&
                   rho(int(npol/2),jpol,iel)
           endif
        endif
        
     enddo
  enddo
  close(6464); close(6465)
14 format(1pe14.5,1pe14.5,3(1pe14.6))
  
! Affix "above" values on the surface to corresponding "below" values
  minvpdomabove(1)=minvpdombelow(1); maxvpdomabove(1)=maxvpdombelow(1)
  minvsdomabove(1)=minvsdombelow(1); maxvsdomabove(1)=maxvsdombelow(1)
  minrodomabove(1)=minrodombelow(1); maxrodomabove(1)=maxrodombelow(1)
  
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! check whether the right number of points resides on each discontinuity
! and compute the velocities above/below the discontinuities
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  write(69,*)' '
  write(69,*)'# points and elastic info around discontinuities'
  write(69,*)
  write(69,133)'','r disc. [km]','# pts','vp [km/s]','vs [km/s]','rho [g/cm^3]'
133 format(a7,a13,a7,3(a13))
  
  write(69,13)'Below',discont(1)/1.d3,domcount(1),minvpdombelow(1)/1.d3,&
       minvsdombelow(1)/1.d3,minrodombelow(1)/1.d3
  
  if (.not. dblreldiff_small(minvpdombelow(1),maxvpdombelow(1))) then
     write(6,*)
     write(6,*)procstrg,&
          'PROBLEM: Vp at the surface is not spherically symmetric!'
     write(6,*)procstrg,'Min/max vp along surface:',&
          minvpdombelow(1),maxvpdombelow(1)
     stop
  endif
  
  if (.not. dblreldiff_small(minvsdombelow(1),maxvsdombelow(1))) then 
     write(6,*)
     write(6,*)procstrg,&
          'PROBLEM: Vs at the surface is not spherically symmetric!'
     write(6,*)procstrg,'Min/max vs along surface:',&
          minvsdombelow(1),maxvsdombelow(1)
     stop
  endif
  
  if (.not. dblreldiff_small(minrodombelow(1),maxrodombelow(1))) then 
     write(6,*)
     write(6,*)procstrg,&
          'PROBLEM: density at the surface is not spherically symmetric!'
     write(6,*)procstrg,&
          'Min/max rho along surface:',minrodombelow(1),maxrodombelow(1)
     stop
  endif
  
!---------------
  do iidom=2,ndisc
!---------------
     write(69,13)'Above',discont(iidom)/1.d3,domcount(iidom), &
          minvpdomabove(iidom)/1.d3,minvsdomabove(iidom)/1.d3, &
          minrodomabove(iidom)/1.d3  
     
     if (.not. dblreldiff_small(minvpdomabove(iidom),maxvpdomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vp not spherically symm. above discont',&
             discont(iidom)
        write(6,*)procstrg,'Min/max vp along discont:',minvpdomabove(iidom), &
             maxvpdomabove(iidom)
        stop
     endif
     
     if (.not. dblreldiff_small(minvsdomabove(iidom),maxvsdomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vs not spherically symm. above discont',&
             discont(iidom)
        write(6,*)procstrg,'Min/max vs along discont:',minvsdomabove(iidom), &
             maxvsdomabove(iidom)
        stop 
     endif
     
     if (.not. dblreldiff_small(minrodomabove(iidom),maxrodomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,&
             'PROBLEM: density not spherically symm. above discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max rho along discont:',minrodomabove(iidom), &
             maxrodomabove(iidom)
        stop
     endif
     
     write(69,13)'Below',discont(iidom)/1.d3,domcount(iidom), &
          minvpdombelow(iidom)/1.d3,minvsdombelow(iidom)/1.d3, &
          minrodombelow(iidom)/1.d3
     
     if (.not. dblreldiff_small(minvpdombelow(iidom),maxvpdombelow(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vp not spherically symm. below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max vp along discont:',minvpdombelow(iidom), &
             maxvpdombelow(iidom)
        stop
     endif
     
     if (.not. dblreldiff_small(minvsdombelow(iidom),maxvsdombelow(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vs not spherically symm. below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max vs along discont:',minvsdombelow(iidom), &
             maxvsdombelow(iidom)
        stop
     endif
     
     if (.not.dblreldiff_small(minrodombelow(iidom),maxrodombelow(iidom))) then
        write(6,*)
        write(6,*)procstrg,&
             'PROBLEM: density not spherically symmetric below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max rho along discont:',minrodombelow(iidom), &
             maxrodombelow(iidom)
        stop
     endif
!---------------
  enddo
!---------------
13 format(a7,1pe12.4,i7,3(1pe13.3))

end subroutine check_elastic_discontinuities
!=============================================================================

!-----------------------------------------------------------------------------
subroutine check_background_model(lambda,mu,rho)

use data_mesh_preloop, ONLY : eltype, coarsing

include 'mesh_params.h'

double precision, intent(in)  :: rho(0:npol,0:npol,nelem)
double precision, intent(in)  :: lambda(0:npol,0:npol,nelem)
double precision, intent(in)  :: mu(0:npol,0:npol,nelem)
integer                       :: iel,ipol,jpol,i
integer, dimension (2)        :: rhominloc,lamminloc,muminloc
integer, dimension (2)        :: rhomaxloc,lammaxloc,mumaxloc
double precision, allocatable :: maxdiffrho(:)
double precision, allocatable :: maxdifflam(:),maxdiffmu(:)

  allocate(maxdiffrho(nelem),maxdifflam(nelem),maxdiffmu(nelem))

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! check whether all grid points assumed non-crazy elastic properties
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  if (minval(rho) < zero) then 
     write(6,*)
     write(6,*)procstrg,'PROBLEM: Not every grid point has a valid density!'
     vsminloc=minloc(rho)
     write(6,*)procstrg,'rho=',minval(rho)
     write(6,*)procstrg,'r,theta:',&
          rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
          thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif
  
  if (minval(mu) < zero) then 
     write(6,*)
     write(6,*)procstrg,&
          'PROBLEM: Not every grid point has a valid mu (Lame) parameter!'
     vsminloc=minloc(mu)
     write(6,*)procstrg,'mu=',minval(mu)
     write(6,*)procstrg,'r,theta:',&
          rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
          thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif
  
  if (minval(lambda) < zero) then 
     write(6,*)
     write(6,*)procstrg,&
            'PROBLEM: Not every grid point has a valid lambda (Lame) parameter'
     vsminloc=minloc(lambda)
     write(6,*)procstrg,'lambda=',minval(lambda)
     write(6,*)procstrg,'r,theta:',&
          rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
          thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif


  do iel = 1, nelem
     if (.not. add_hetero) then 
        if ( eltype(iel)=='curved' .and. .not. coarsing(iel) ) then
           do jpol = 0, npol
              do ipol = 1, npol-1
   
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   ! Test lateral homogeneity: Does each point have same elastic properties 
   ! as its neighbor in the xi-direction?
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   
                 if (.not. dblreldiff_small(rho(ipol,jpol,iel), &
                      rho(ipol-1,jpol,iel)) .or. &
                      .not. dblreldiff_small(rho(ipol,jpol,iel), &
                      rho(ipol+1,jpol,iel))) then
                    write(6,*)
                    write(6,*)procstrg,'PROBLEM: lateral density inhomogeneity!'
                    write(6,*)procstrg,' at r[km], theta[deg]:',&
                         rcoord(ipol,jpol,iel)/1000., &
                          thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*)procstrg,'Lateral ipol,theta,density:'
                    do i=0,npol
                       write(6,*)procstrg,i,thetacoord(i,jpol,iel)*180./pi,&
                            rho(i,jpol,iel)
                    enddo
                    stop
                 endif
                 
                 if (.not. dblreldiff_small(lambda(ipol,jpol,iel), &
                      lambda(ipol-1,jpol,iel)) .or. &
                      .not. dblreldiff_small(lambda(ipol,jpol,iel), &
                      lambda(ipol+1,jpol,iel)) ) then 
                    write(6,*)
                    write(6,*)procstrg,'PROBLEM: lateral elastic inhomogeneity!'
                    write(6,*)procstrg,'at r[km], theta[deg]:',&
                         rcoord(ipol,jpol,iel)/1000., &
                         thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*)procstrg,'Lateral ipol,theta,lambda:'
                    do i=0,npol
                       write(6,*)procstrg,i,thetacoord(i,jpol,iel)*180./pi,&
                            lambda(i,jpol,iel)
                    enddo
                    stop
                 endif
                 
                 if (.not. dblreldiff_small(mu(ipol,jpol,iel), &
                      mu(ipol-1,jpol,iel)) .or. &
                      .not. dblreldiff_small(mu(ipol,jpol,iel), &
                      mu(ipol+1,jpol,iel)) ) then 
                    write(6,*)
                    write(6,*)procstrg,'PROBLEM: lateral elastic inhomogeneity!'
                    write(6,*)procstrg,'at r[km], theta[deg]:',&
                         rcoord(ipol,jpol,iel)/1000., &
                         thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*)'Lateral ipol,theta,mu:'
                    do i=0,npol
                       write(6,*)procstrg,i,thetacoord(i,jpol,iel)*180./pi,&
                            mu(i,jpol,iel)
                    enddo
                    stop
                 endif
        
   
              enddo
           enddo
        endif ! only curved elements
     endif ! add_hetero
           


! Compute maximal elemental variation in medium properties
     rhominloc=minloc(rho(0:npol,0:npol,iel))
     lamminloc=minloc(lambda(0:npol,0:npol,iel))
     muminloc=minloc(mu(0:npol,0:npol,iel))
     rhomaxloc=maxloc(rho(0:npol,0:npol,iel))
     lammaxloc=maxloc(lambda(0:npol,0:npol,iel))
     mumaxloc=maxloc(mu(0:npol,0:npol,iel))
     
     maxdiffrho(iel)=dbleabsreldiff(rho(rhomaxloc(1)-1,rhomaxloc(2)-1,iel), &
          rho(rhominloc(1)-1,rhominloc(2)-1,iel))
     maxdifflam(iel)=dbleabsreldiff(lambda(lammaxloc(1)-1,lammaxloc(2)-1,iel),&
          lambda(lamminloc(1)-1,lamminloc(2)-1,iel) )
     maxdiffmu(iel)=dbleabsreldiff(mu(mumaxloc(1)-1,mumaxloc(2)-1,iel), & 
          mu(muminloc(1)-1,muminloc(2)-1,iel))
  enddo
  
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Maximal elemental variation in medium properties
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  open(unit=6464,file=infopath(1:lfinfo)//&
       '/model_variation_across_elem.dat'//appmynum)
  open(unit=4,file=infopath(1:lfinfo)//&
       '/model_var_morethan_12percent.dat'//appmynum)
  do iel = 1,nelem
     if (maxdiffrho(iel) >= 0.12) then
        write(4,*)
        write(4,*)'WARNING: Density varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     if (maxdifflam(iel) >= 0.12) then 
        write(4,*)
        write(4,*)'WARNING: Lambda varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     if (maxdiffmu(iel) >= 0.12) then 
        write(4,*)
        write(4,*)'WARNING: Mu varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     write(6464,64)rcoord(npol/2,npol/2,iel)/1000., &
          thetacoord(npol/2,npol/2,iel)*180/pi, &
          maxdiffrho(iel),maxdifflam(iel),maxdiffmu(iel)
  enddo
  close(4);close(6464)
  
64 format(1pe14.4,1pe13.3,3(1pe15.5))
  
  write(69,*)
  write(69,*)'Maximal variation in density across one element:', &
              maxval(maxdiffrho)
  write(69,*)'at r [km], theta[deg]:', &
       rcoord(npol/2,npol/2,maxloc(maxdiffrho,1))/1000., &
       thetacoord(npol/2,npol/2,maxloc(maxdiffrho,1))*180./pi
  write(69,*)
  write(69,*)'Maximal variation in lambda across one element:',&
              maxval(maxdifflam)
  write(69,*)'at r [km], theta[deg]:', &
       rcoord(npol/2,npol/2,maxloc(maxdifflam,1))/1000., &
       thetacoord(npol/2,npol/2,maxloc(maxdifflam,1))*180./pi
  write(69,*)
  write(69,*)'Maximal variation in mu across one element:',maxval(maxdiffmu)
  write(69,*)'at r [km], theta[deg]:', &
       rcoord(npol/2,npol/2,maxloc(maxdiffmu,1))/1000., &
       thetacoord(npol/2,npol/2,maxloc(maxdiffmu,1))*180./pi

  deallocate(maxdiffrho,maxdifflam,maxdiffmu)

end subroutine check_background_model
!=============================================================================

!-----------------------------------------------------------------------------
subroutine test_mesh_model_resolution(lambda,mu,rho)
!
! Resolution test: Insert radial sine function with wavenumbers according to 
! actual wavelengths, compute numerical and analytical integrals for various 
! setups by changing the wavenumber by source period and seismic velocities.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use def_grid,  ONLY : massmatrix,massmatrix_dble
use data_time, ONLY : period
use data_mesh_preloop, ONLY : eltype, coarsing, north, axis
use commun, ONLY : psum_dble,broadcast_int,broadcast_dble

include 'mesh_params.h'

double precision, intent(in)  :: rho(0:npol,0:npol,nelem)
double precision, intent(in)  :: lambda(0:npol,0:npol,nelem)
double precision, intent(in)  :: mu(0:npol,0:npol,nelem)
double precision, allocatable :: fp(:,:,:),fs(:,:,:)
double precision, allocatable :: mass(:,:,:)
double precision              :: intp_num,ints_num
double precision              :: intp_ana,ints_ana,arg
integer                       :: iel,ipol,jpol,i,irad,iirad,iidom
double precision              :: tmpradii(naxel,2)
double precision              :: vptmp,vstmp
double precision, allocatable :: radii(:,:),vel(:,:)

  allocate(fp(0:npol,0:npol,nelem),fs(0:npol,0:npol,nelem))
  allocate(mass(0:npol,0:npol,nelem))

  do iidom = 1,2
     write(69,*)
     if (iidom == 1) then
        vptmp=minval(dsqrt( (lambda+2.d0*mu )/rho ))
        !vstmp=minval(sqrt( mu/rho))
        vstmp = vptmp/dsqrt(3.d0)
        write(69,*)':::::::::Resolution test for globally MINIMAL vp, ', &
             'vs::::::::::::::::::'
        write(69,*)'     sin(2 pi/(T0*minv) r) integrated over 3-D volume'
     elseif (iidom ==2) then
        vptmp=maxval(dsqrt( (lambda+2.d0*mu )/rho ))
        !vstmp=minval(sqrt( mu/rho))
        vstmp = vptmp/dsqrt(3.d0)
        write(69,*)':::::::::Resolution test for globally MAXIMAL vp, ', &
             'vs::::::::::::::::::'
        write(69,*)'     sin(2 pi/(T0*maxv) r) integrated over 3-D volume'
     endif

     do i = 1,4 ! 0.5*period, period, 1.5*period, 2*period
        do iel = 1, nelem
           do jpol = 0, npol
              do ipol = 0,npol
                 fp(ipol,jpol,iel) = dsin(two*pi*rcoord(ipol,jpol,iel)/ & 
                      (dble(i)/2.d0*period*vptmp))
                 if (vstmp > zero) then 
                    fs(ipol,jpol,iel) = dsin(two*pi*rcoord(ipol,jpol,iel)/ & 
                         (dble(i)/2.d0*period*vstmp))
                 else
                    fs(ipol,jpol,iel) = fp(ipol,jpol,iel)
                 endif
              enddo
           enddo
        enddo

        if (iidom ==3) fp = one

! Numerical integration:
        call massmatrix_dble(mass,nelem,'total')
        intp_num = zero; ints_num = zero
        do iel = 1, nelem
           do ipol = 0, npol
              do jpol = 0, npol
                 intp_num = intp_num + fp(ipol,jpol,iel)*mass(ipol,jpol,iel)
                 ints_num = ints_num + fs(ipol,jpol,iel)*mass(ipol,jpol,iel)
              end do
           end do
        end do
        intp_num=2.d0*pi*intp_num
        ints_num=2.d0*pi*ints_num

        intp_num = psum_dble(intp_num)
        ints_num = psum_dble(ints_num)

! Analytical integration: 
! \int sin(ar) r^2 dr = 2r/a^2 sin(ar) - (r^2/a - 2/a^3) cos(ar)
!    where a = 2 pi/(T0 v_{s,p})
! Here: analytically for whole earth if v_{s,p} = const

! Constant velocities:
        arg = two*pi/(real(i)/2.d0*period * vptmp)
        intp_ana = four*pi*( two*router/arg**2 * dsin(arg*router) - &
             (router**2/arg-two/arg**3)*dcos(arg*router) - two/arg**3 )
        arg = two*pi/(real(i)/2.*period * vstmp)
        ints_ana = four*pi*( two*router/arg**2 * dsin(arg*router) - &
             (router**2/arg-two/arg**3)*dcos(arg*router) - two/arg**3 )

        write(69,*)'T0=period*',real(i)/2.
        write(69,1234)' Num., ana., diff vp:',intp_num,intp_ana, &
             dbleabsreldiff(intp_num,intp_ana)
        write(69,1234)' Num., ana., diff vs:',ints_num,ints_ana, &
             dbleabsreldiff(ints_num,ints_ana)
1234    format(a22,2(1pe18.8),1pe13.3)

     enddo ! 4 periods
  enddo ! min/max velocities

! Same but for constant function, i.e. the pure volume of the sphere
  write(69,*)
  write(69,*)':::::::::Resolution test for constant function::::', &
       '::::::::::::::::::::'
  write(69,*)'             3-D spherical volume itself'
  call massmatrix_dble(mass,nelem,'total')
  intp_num = zero
  do iel = 1, nelem
     do ipol = 0, npol
        do jpol = 0, npol
           intp_num = intp_num + mass(ipol,jpol,iel)
        end do
     end do
  end do
  intp_num=2.d0*pi*intp_num
  intp_num = psum_dble(intp_num)

  intp_ana = four/three*pi*router**3

  write(69,1234)' Num., ana., diff vp:',intp_num,intp_ana, &
       dbleabsreldiff(intp_num,intp_ana)

! Piecewise constant velocities over each element,
! only for spheroidal element shapes and leaving out coarsening levels

  irad=0
  if (mynum==0) then
     do iel = 1, nelem
        if (eltype(iel)=='curved' .and. .not. coarsing(iel) ) then
           ! radii needed for analytical integration
           if (axis(iel) .and. north(iel)) then
              irad = irad + 1
              tmpradii(irad,1) = rcoord(0,npol,iel)
              tmpradii(irad,2) = rcoord(0,0,iel)
           endif
        endif !eltype curved
     enddo
  endif
  call broadcast_int(irad,0) 
  allocate(radii(irad,2),vel(irad,2))
  if (mynum==0) radii(1:irad,1:2) = tmpradii(1:irad,1:2)
  do iirad=1,irad
     call broadcast_dble(radii(iirad,1),0)
     call broadcast_dble(radii(iirad,2),0)
  enddo

  open(unit=779,file=infopath(1:lfinfo)//'/resolutionsine.dat'//appmynum)

  do iidom = 1,2
     write(69,*)
     if (iidom == 1) then
        open(unit=779,file=infopath(1:lfinfo)// &
             '/resolutionsine_elminvp.dat'//appmynum)      
        write(69,*)':::::::::Resolution test for elementally MINIMAL vp, ', &
             'vs:::::::::::::::'
        write(69,*)'     sin(2 pi/(T0*minv) r) integrated over 3-D volume'
     elseif (iidom ==2) then
        open(unit=779,file=infopath(1:lfinfo)//'/resolutionsine_elmaxvp.dat' &
             //appmynum)  
        write(69,*)':::::::::Resolution test for elementally MAXIMAL vp, ', &
             'vs:::::::::::::::'
        write(69,*)'     sin(2 pi/(T0*maxv) r) integrated over 3-D volume'
     endif

     fp=zero; fs=zero; iirad=0; 
     do iel = 1, nelem
        if (eltype(iel)=='curved' .and. .not. coarsing(iel) ) then

           if (iidom == 1) then
              vptmp=minval(dsqrt((lambda(:,:,iel)+&
                                 two*mu(:,:,iel))/rho(:,:,iel)))
              vstmp= minval( dsqrt( mu(:,:,iel) / rho(:,:,iel) ) )
              if (vstmp==zero) vstmp = vptmp  
           elseif (iidom ==2) then
              vptmp=maxval(dsqrt((lambda(:,:,iel)+&
                                  two*mu(:,:,iel))/rho(:,:,iel)))
              vstmp= maxval( dsqrt( mu(:,:,iel) / rho(:,:,iel) )) 
              if (vstmp==zero) vstmp = vptmp  
           endif
           do jpol = 0, npol
              do ipol = 0,npol
                 fp(ipol,jpol,iel)= dsin(two*pi*rcoord(ipol,jpol,iel)/ &
                      (period*vptmp))
                 fs(ipol,jpol,iel)= dsin(two*pi*rcoord(ipol,jpol,iel)/ &
                      (period*vstmp))
              enddo
              if (axis(iel) .and. north(iel)) then
                 write(779,*) rcoord(0,jpol,iel),fp(0,jpol,iel),fs(0,jpol,iel)
              endif
           enddo
          ! radii needed for analytical integration
           if (axis(iel) .and. north(iel) .and. mynum==0) then
              iirad = iirad + 1
              vel(iirad,1) = vptmp
              vel(iirad,2) = vstmp
           endif

        endif !eltype curved
     enddo

     do iirad=1,irad
        call broadcast_dble(vel(iirad,1),0)
        call broadcast_dble(vel(iirad,2),0)
     enddo

     close(779)

! Numerical integration:
     call massmatrix_dble(mass,nelem,'total')
     intp_num = zero; ints_num = zero
     do iel = 1, nelem
        do ipol = 0, npol
           do jpol = 0, npol
              intp_num = intp_num + fp(ipol,jpol,iel)*mass(ipol,jpol,iel)
              ints_num = ints_num + fs(ipol,jpol,iel)*mass(ipol,jpol,iel)
           end do
        end do
     end do
     intp_num=2.d0*pi*intp_num
     ints_num=2.d0*pi*ints_num

     intp_num = psum_dble(intp_num)
     ints_num = psum_dble(ints_num)

! analytical integration
     intp_ana = zero; ints_ana = zero; 
     do i=1,irad
        arg = two*pi/(period * vel(i,1))
        intp_ana = intp_ana + &
             four*pi*(two*radii(i,1)/arg**2*dsin(arg*radii(i,1)) - &
             (radii(i,1)**2/arg-two/arg**3)*dcos(arg*radii(i,1)) - &
             two*radii(i,2)/arg**2*dsin(arg*radii(i,2)) + &
             (radii(i,2)**2/arg-two/arg**3)*dcos(arg*radii(i,2)) )
        arg = two*pi/(period * vel(i,2))
        ints_ana = ints_ana + &
             four*pi*(two*radii(i,1)/arg**2*dsin(arg*radii(i,1)) - &
             (radii(i,1)**2/arg-two/arg**3)*dcos(arg*radii(i,1)) - &
             two*radii(i,2)/arg**2*dsin(arg*radii(i,2)) + &
             (radii(i,2)**2/arg-two/arg**3)*dcos(arg*radii(i,2)) )
     enddo

     write(69,1234)' Num., ana., diff vp:',intp_num,intp_ana, &
          dbleabsreldiff(intp_num,intp_ana)
     write(69,1234)' Num., ana., diff vs:',ints_num,ints_ana, &
          dbleabsreldiff(ints_num,ints_ana)

  enddo ! min/max velocities

  write(69,*)

  deallocate(fp,fs)
  deallocate(mass)

end subroutine test_mesh_model_resolution
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_numerical_resolution(lambda,mu,rho)
!
! MvD: This function is never called inside the solver!!
!
! Compute critical numerical parameters such as time step, source period, 
! courant number, spacing variations, characteristic spacing lead time, 
! or points per wavelength throughout the domain.
!
! relationships/definitions: 
!   characteristic lead time: spacing/velocity (taup for vp, taus for vs)
!   courant number: (often heuristically determined) stability criterion
!   points per wavelength ptslam: number of gridpoints per DOMINANT wavelength
!
!   dt <= courant * min(taup)   ===> usually set by lower mantle/inner core
!   period >= ptslam * max(taus) ==> usually set by crust
!
! The mesher reads courant and ptslam as input parameters and determines 
! dt and period from the above. Here we validate these computed parameters.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_mesh_preloop
use data_time, ONLY : period,deltat,courant
use commun, ONLY : pmin,pmax

double precision, intent(in)                  :: rho(0:npol,0:npol,nelem)
double precision, intent(in)                  :: lambda(0:npol,0:npol,nelem)
double precision, intent(in)                  :: mu(0:npol,0:npol,nelem)
double precision,dimension(:,:,:),allocatable :: dis
double precision,dimension(:,:),allocatable   :: vs,vp
double precision,dimension(:,:,:),allocatable :: chartime_vp,chartime_vs
double precision,dimension(:),allocatable     :: char_lead_time_vp_min
double precision,dimension(:),allocatable     :: char_lead_time_vp_max
double precision,dimension(:),allocatable     :: char_lead_time_vs_min
double precision,dimension(:),allocatable     :: char_lead_time_vs_max
double precision,dimension(:),allocatable     :: char_lead_time_vp_min_rad
double precision,dimension(:),allocatable     :: char_lead_time_vp_min_lat
double precision,dimension(:),allocatable     :: char_lead_time_vs_max_rad
double precision,dimension(:),allocatable     :: char_lead_time_vs_max_lat
double precision,dimension(:),allocatable     :: char_lead_vp_radmin
double precision,dimension(:),allocatable     :: char_lead_vp_radmax
double precision,dimension(:),allocatable     :: char_lead_vs_radmin
double precision,dimension(:),allocatable     :: char_lead_vs_radmax
double precision,dimension(:),allocatable     :: char_lead_vp_coarsmin
double precision,dimension(:),allocatable     :: char_lead_vp_coarsmax
double precision,dimension(:),allocatable     :: char_lead_vs_coarsmin
double precision,dimension(:),allocatable     :: char_lead_vs_coarsmax
double precision,dimension(:),allocatable     :: char_lead_vp_centralmin
double precision,dimension(:),allocatable     :: char_lead_vp_centralmax
double precision,dimension(:),allocatable     :: char_lead_vs_centralmin
double precision,dimension(:),allocatable     :: char_lead_vs_centralmax
double precision,dimension(:),allocatable     :: scentral,zcentral
double precision,dimension(:),allocatable     :: rcentral,thcentral
double precision,dimension(:),allocatable     :: rad_coars,dt_rad,period_rad
double precision                              :: mincharspher,minradspher 
double precision                              :: minthspher,maxcharspher
double precision                              :: maxradspher,maxthspher
integer                                       :: iel,ipol,jpol,i
integer                                       :: icub,ncub,iirad,idisc,ncoars

! elemental grid spacing
  allocate(dis(0:npol-1,0:npol-1,4))
  allocate(chartime_vp(0:npol-1,0:npol-1,4))
  allocate(chartime_vs(0:npol-1,0:npol-1,4))
  allocate(vs(0:npol-1,0:npol-1))
  allocate(vp(0:npol-1,0:npol-1))

! elemental characteristic lead times
  allocate(char_lead_time_vp_min(nelem))
  allocate(char_lead_time_vp_max(nelem))
  allocate(char_lead_time_vs_min(nelem))
  allocate(char_lead_time_vs_max(nelem))

  allocate(char_lead_time_vp_min_rad(nelem))
  allocate(char_lead_time_vp_min_lat(nelem))
  allocate(char_lead_time_vs_max_rad(nelem))
  allocate(char_lead_time_vs_max_lat(nelem))

! time step and period for radial bins
  allocate(dt_rad(num_spher_radii))
  allocate(period_rad(num_spher_radii))
  dt_rad=10000000.d0
  period_rad=zero

! spherical part
  allocate(char_lead_vp_radmin(num_spher_radii))
  allocate(char_lead_vp_radmax(num_spher_radii))
  allocate(char_lead_vs_radmin(num_spher_radii))
  allocate(char_lead_vs_radmax(num_spher_radii))
  char_lead_vp_radmin = 10*router
  char_lead_vp_radmax = zero
  char_lead_vs_radmin = 10*router
  char_lead_vs_radmax = zero

! coarsening levels
  allocate(rad_coars(ndisc))
  allocate(char_lead_vp_coarsmin(ndisc))
  allocate(char_lead_vp_coarsmax(ndisc))
  allocate(char_lead_vs_coarsmin(ndisc))
  allocate(char_lead_vs_coarsmax(ndisc))
  char_lead_vp_coarsmin = 10*router
  char_lead_vp_coarsmax = zero
  char_lead_vs_coarsmin = 10*router
  char_lead_vs_coarsmax = zero
  ncoars=0

! central cube + buffer/transition level
  ncub=0
  do iel=1,nelem
     if ( rcoord(npol/2,npol/2,iel) <= spher_radii(num_spher_radii) .and. &
          eltype(iel)/='curved') ncub=ncub+1    
  enddo 
  write(69,*)
  write(69,*)'  Elements in central cube and transition level:',ncub

  allocate(char_lead_vp_centralmin(ncub))
  allocate(char_lead_vp_centralmax(ncub))
  allocate(char_lead_vs_centralmin(ncub))
  allocate(char_lead_vs_centralmax(ncub))
  allocate(scentral(ncub),zcentral(ncub))
  char_lead_vp_centralmin = 10*router
  char_lead_vp_centralmax = zero
  char_lead_vs_centralmin = 10*router
  char_lead_vs_centralmax = zero
  icub=0

! Find elemental grid spacing minimum and maximum  
  do iel = 1,nelem
     do ipol=0,npol-1
        do jpol=0,npol-1
           
           dis(ipol,jpol,1) = dsqrt(&
                (scoord(ipol,jpol,iel)-scoord(ipol+1,jpol,iel))**2&
                +(zcoord(ipol,jpol,iel)-zcoord(ipol+1,jpol,iel))**2)

           dis(ipol,jpol,2) = dsqrt(&
                (scoord(ipol,jpol,iel)-scoord(ipol,jpol+1,iel))**2&
                +(zcoord(ipol,jpol,iel)-zcoord(ipol,jpol+1,iel))**2)

! include diagonal distances: crucial for deformed inner cube elements!
           dis(ipol,jpol,3) = dsqrt(&
                (scoord(ipol,jpol,iel)-scoord(ipol+1,jpol+1,iel))**2&
                +(zcoord(ipol,jpol,iel)-zcoord(ipol+1,jpol+1,iel))**2)

           dis(ipol,npol-jpol-1,4) = dsqrt(&
                (scoord(ipol+1,npol-jpol-1,iel)-scoord(ipol,npol-jpol,iel))**2&
               +(zcoord(ipol+1,npol-jpol-1,iel)-zcoord(ipol,npol-jpol,iel))**2)

           vs(ipol,jpol)= dsqrt( mu(ipol,jpol,iel) / rho(ipol,jpol,iel) )
           vp(ipol,jpol)= dsqrt((lambda(ipol,jpol,iel)+two*mu(ipol,jpol,iel))/&
                                rho(ipol,jpol,iel) )
        enddo
     enddo


     do i=1,4
        if (minval(vs)>zero) then 
           chartime_vs(0:npol-1,0:npol-1,i) = dis(0:npol-1,0:npol-1,i)/ &
                                               vs(0:npol-1,0:npol-1)
        else ! fluid
           chartime_vs(0:npol-1,0:npol-1,i) = dis(0:npol-1,0:npol-1,i)/ &
                                               vp(0:npol-1,0:npol-1)
        endif
        chartime_vp(0:npol-1,0:npol-1,i) = dis(0:npol-1,0:npol-1,i)/ &
                                            vp(0:npol-1,0:npol-1)  
     enddo
     
     char_lead_time_vp_min(iel)=minval(chartime_vp(:,:,1:2))
     char_lead_time_vs_min(iel)=minval(chartime_vs(:,:,1:2))

     char_lead_time_vp_min_rad(iel)=minval(chartime_vp(:,:,2))
     char_lead_time_vp_min_lat(iel)=minval(chartime_vp(:,:,1))

!    Accomodate the fact that dis(:,;,3:4) are for diagonal distances, 
!    hence divide by sqrt(2) when checking largest distances since the 
!    diagonal distance does not per say play a role, only in relation to very 
!    distorted shapes (i.e. very small, or when even after this division 
!    larger than ipol- or j-pol parallel distances).
     chartime_vp(:,:,3:4)=chartime_vp(:,:,3:4)/dsqrt(two)
     chartime_vs(:,:,3:4)=chartime_vs(:,:,3:4)/dsqrt(two)

     char_lead_time_vp_max(iel)=maxval(chartime_vp)
     char_lead_time_vs_max(iel)=maxval(chartime_vs)

     char_lead_time_vs_max_rad(iel)=maxval(chartime_vs(:,:,2))
     char_lead_time_vs_max_lat(iel)=maxval(chartime_vs(:,:,1))

!     char_lead_time_vp_max(iel)=maxval(chartime_vp(:,:,1:2))
!     char_lead_time_vs_max(iel)=maxval(chartime_vs(:,:,1:2))

! for quick checks: output deltat and period radially
     do iirad = 1,num_spher_radii-1       
        if ( rcoord(npol/2,npol/2,iel) <= spher_radii(iirad) .and. & 
             rcoord(npol/2,npol/2,iel) > spher_radii(iirad+1) ) then 

           dt_rad(iirad)=min(dt_rad(iirad),char_lead_time_vp_min(iel))
           period_rad(iirad)=max(period_rad(iirad),char_lead_time_vs_max(iel))

        endif
     enddo
     if ( rcoord(npol/2,npol/2,iel) <= spher_radii(num_spher_radii) ) then 
        dt_rad(num_spher_radii) = &
                 min(dt_rad(num_spher_radii),char_lead_time_vp_min(iel))
        period_rad(num_spher_radii) = &
                 max(period_rad(num_spher_radii),char_lead_time_vs_max(iel))
     endif

! check parameters between all spherical radii
     if (.not. coarsing(iel)) then 
        do iirad = 1,num_spher_radii-1
           if ( rcoord(npol/2,npol/2,iel) <= spher_radii(iirad) .and. & 
                rcoord(npol/2,npol/2,iel) > spher_radii(iirad+1) ) then 
              
              char_lead_vp_radmin(iirad) = min(char_lead_vp_radmin(iirad),&
                                               char_lead_time_vp_min(iel))
              
              char_lead_vp_radmax(iirad) = max(char_lead_vp_radmax(iirad),&
                                               char_lead_time_vp_max(iel))
              
              char_lead_vs_radmin(iirad) = min(char_lead_vs_radmin(iirad),&
                                               char_lead_time_vs_min(iel))
              
              char_lead_vs_radmax(iirad) = max(char_lead_vs_radmax(iirad),&
                                               char_lead_time_vs_max(iel))
           endif
        enddo

! put all h/v ratios into bins for each coarsening layer (actually each 
! layer in between discontinuities, i.e. some will remain empty)
     else
        do idisc = 1,ndisc-1
           if ( rcoord(npol/2,npol/2,iel)<=discont(idisc) .and. &
                rcoord(npol/2,npol/2,iel) > discont(idisc+1) ) then 

              rad_coars(idisc) = rcoord(npol/2,npol/2,iel)
              char_lead_vp_coarsmin(idisc) = min(char_lead_vp_coarsmin(idisc),&
                                                 char_lead_time_vp_min(iel))
              
              char_lead_vp_coarsmax(idisc) = max(char_lead_vp_coarsmax(idisc),&
                                                 char_lead_time_vp_max(iel))
              
              char_lead_vs_coarsmin(idisc) = min(char_lead_vs_coarsmin(idisc),&
                                                 char_lead_time_vs_min(iel))
              
              char_lead_vs_coarsmax(idisc) = max(char_lead_vs_coarsmax(idisc),&
                                                 char_lead_time_vs_max(iel))
          endif
       enddo
    endif

! check parameters below spherical radii, i.e. central cube part
     if ( rcoord(npol/2,npol/2,iel) <= spher_radii(num_spher_radii) ) then 

!!$        char_lead_vp_radmin(num_spher_radii) = &
!!$           min(char_lead_vp_radmin(num_spher_radii),char_lead_time_vp_min(iel))
!!$                                               
!!$        char_lead_vp_radmax(num_spher_radii) = &
!!$           max(char_lead_vp_radmax(num_spher_radii),char_lead_time_vp_max(iel))
!!$
!!$        char_lead_vs_radmin(num_spher_radii) = &
!!$           min(char_lead_vs_radmin(num_spher_radii),char_lead_time_vs_min(iel))
!!$        
!!$        char_lead_vs_radmax(num_spher_radii) = &
!!$           max(char_lead_vs_radmax(num_spher_radii),char_lead_time_vs_max(iel))


! additionally save central cube & buffer into separate arrays
        if ( eltype(iel) /= 'curved' ) then 
          icub = icub + 1 
          char_lead_vp_centralmin(icub) = char_lead_time_vp_min(iel)
          char_lead_vp_centralmax(icub) = char_lead_time_vp_max(iel)
          char_lead_vs_centralmin(icub) = char_lead_time_vs_min(iel)
          char_lead_vs_centralmax(icub) = char_lead_time_vs_max(iel)          
          scentral(icub) = scoord(npol/2,npol/2,iel)
          zcentral(icub) = zcoord(npol/2,npol/2,iel)
        endif

     endif

  enddo ! iel

  deallocate(dis)
  deallocate(chartime_vp)
  deallocate(chartime_vs)

! Compare given time step and period to the one calculated here
  write(69,*)
  write(69,*)'  Read-in time step:   ',deltat
  write(69,*)'  Calculated time step:',courant*minval(char_lead_time_vp_min)
  write(69,*)     
  write(69,*)'  Read-in period:   ',period
  write(69,*)'  Calculated period:',pts_wavelngth*dble(npol)*&
                                    maxval(char_lead_time_vs_max)

! Locate minimal and maximal vp- and vs-characteristic lead times
  write(69,*)
  write(69,*)'Minimal/maximal global characteristic spacing lead times:'
  write(69,9)'  min value,r,theta vp:',minval(char_lead_time_vp_min), &
              rcoord(npol/2,npol/2,minloc(char_lead_time_vp_min,1))/1.d3, &
              thetacoord(npol/2,npol/2,minloc(char_lead_time_vp_min,1))*180./pi
  write(69,9)'  min value,r,theta vs:',minval(char_lead_time_vs_min), &
              rcoord(npol/2,npol/2,minloc(char_lead_time_vs_min,1))/1.d3, &
              thetacoord(npol/2,npol/2,minloc(char_lead_time_vs_min,1))*180./pi
  write(69,9)'  max value,r,theta vp:',maxval(char_lead_time_vp_max), &
              rcoord(npol/2,npol/2,maxloc(char_lead_time_vp_max,1))/1.d3, &
              thetacoord(npol/2,npol/2,maxloc(char_lead_time_vp_max,1))*180./pi
  write(69,9)'  max value,r,theta vs:',maxval(char_lead_time_vs_max), &
              rcoord(npol/2,npol/2,maxloc(char_lead_time_vs_max,1))/1.d3, &
              thetacoord(npol/2,npol/2,maxloc(char_lead_time_vs_max,1))*180./pi

  write(69,*)'Min/max spherical radial characteristic spacing lead times:'
  write(69,8)'  min value,r vp:',minval(char_lead_vp_radmin), &
              spher_radii(minloc(char_lead_vp_radmin,1))/1.d3
  write(69,8)'  min value,r vs:',minval(char_lead_vs_radmin), &
              spher_radii(minloc(char_lead_vs_radmin,1))/1.d3
  write(69,8)'  max value,r vp:',maxval(char_lead_vp_radmax), &
              spher_radii(minloc(char_lead_vp_radmax,1))/1.d3
  write(69,8)'  max value,r vs:',maxval(char_lead_vs_radmax), &
              spher_radii(minloc(char_lead_vs_radmax,1))/1.d3

! Find min/max elementally ONLY in spherical part of the domain
  mincharspher = router
  maxcharspher = zero
  do iel=1,nelem
     if (eltype(iel)/='linear') then
        if (char_lead_time_vp_min(iel) < mincharspher) then 
           mincharspher=char_lead_time_vp_min(iel)
           minradspher=rcoord(npol/2,npol/2,iel)/1.d3
           minthspher=thetacoord(npol/2,npol/2,iel)*180./pi
        endif
        if (char_lead_time_vs_max(iel) > maxcharspher) then 
           maxcharspher=char_lead_time_vs_max(iel)
           maxradspher=rcoord(npol/2,npol/2,iel)/1.d3
           maxthspher=thetacoord(npol/2,npol/2,iel)*180./pi
        endif
     endif 
  enddo

  write(69,*)'Min/max characteristic spacing lead times, NON-LINEAR elements:'
  write(69,9)'  min value,r,theta vp:',mincharspher,minradspher,minthspher
  write(69,9)'  max value,r,theta vs:',maxcharspher,maxradspher,maxthspher

! Find min/max characteristic spacing lead times in central cube
  allocate(rcentral(ncub),thcentral(ncub))
  rcentral = dsqrt(scentral**2+zcentral**2)/1.d3
  thcentral = datan(scentral/(zcentral+1.e-30))
  do icub = 1,ncub
     if (thcentral(icub)<zero) thcentral(icub) = pi + thcentral(icub)
     if (thcentral(icub) == zero .and. zcentral(icub) < zero) & 
          thcentral(icub)=pi
  enddo
  thcentral = thcentral*180./pi

  write(69,*)'Minimal/maximal central-cube characteristic spacing lead times:'
  write(69,9)'  min value,r,theta vp:',minval(char_lead_vp_centralmin), &
              rcentral(minloc(char_lead_vp_centralmin,1)), &
              thcentral(minloc(char_lead_vp_centralmin,1))
  write(69,9)'  min value,r,theta vs:',minval(char_lead_vs_centralmin), &
              rcentral(minloc(char_lead_vs_centralmin,1)), &
              thcentral(minloc(char_lead_vs_centralmin,1))
  write(69,9)'  max value,r,theta vp:',maxval(char_lead_vp_centralmax), &
              rcentral(minloc(char_lead_vp_centralmax,1)), &
              thcentral(minloc(char_lead_vp_centralmax,1))
  write(69,9)'  max value,r,theta vs:',maxval(char_lead_vs_centralmax), &
              rcentral(minloc(char_lead_vs_centralmax,1)), &
              thcentral(minloc(char_lead_vs_centralmax,1))

  deallocate(rcentral,thcentral)

!9 format(a25,1pe13.3,f12.1,f7.1)
!8 format(a25,1pe13.3,f12.1)

9 format(a25,3(1pe13.3))
8 format(a25,2(1pe13.3))

!\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\
! save characteristic lead times elementally for each processor
  open(unit=60009,file=infopath(1:lfinfo)//'/char_lead_time_el.dat'//appmynum)
  open(unit=60008,file=infopath(1:lfinfo)//'/char_lead_time_el_radlat.dat'//appmynum)
  do iel=1,nelem
     write(60009,13)scoord(npol/2,npol/2,iel),zcoord(npol/2,npol/2,iel),&
                    char_lead_time_vp_min(iel),char_lead_time_vp_max(iel), &
                    char_lead_time_vs_min(iel),char_lead_time_vs_max(iel)

     write(60008,13)scoord(npol/2,npol/2,iel),zcoord(npol/2,npol/2,iel),&
                    char_lead_time_vp_min_lat(iel),char_lead_time_vp_min_rad(iel),&
                    char_lead_time_vs_max_lat(iel),char_lead_time_vs_max_rad(iel)
  enddo
  close(60009); close(60008)
  deallocate(char_lead_time_vp_min)
  deallocate(char_lead_time_vp_max)
  deallocate(char_lead_time_vs_min)
  deallocate(char_lead_time_vs_max)

  deallocate(char_lead_time_vp_min_rad)
  deallocate(char_lead_time_vp_min_lat)
  deallocate(char_lead_time_vs_max_rad)
  deallocate(char_lead_time_vs_max_lat)

! save spherical average for time step and period for given radius GLOBALLY
   dt_rad=courant*dt_rad
   period_rad=dble(npol)*pts_wavelngth*period_rad
  if (lpr) &
       open(unit=60008,file=infopath(1:lfinfo)//'/period_rad.dat')
  if (lpr) &
       open(unit=60009,file=infopath(1:lfinfo)//'/timestep_rad.dat')

  do iirad=1,num_spher_radii-1
     period_rad(iirad) = pmin(period_rad(iirad))
     dt_rad(iirad) = pmax(dt_rad(iirad))
     if (lpr) then
       write(60008,11) spher_radii(iirad),period_rad(iirad)
       write(60009,11) spher_radii(iirad),dt_rad(iirad)
    endif
  enddo
  if (lpr) close(60009)
  if (lpr) close(60008)
  deallocate(dt_rad)
  deallocate(period_rad)

! save spherical average over given radius GLOBALLY
  if (lpr) &
       open(unit=60009,file=infopath(1:lfinfo)//'/char_lead_time_spher.dat')
  do iirad=1,num_spher_radii-1
     char_lead_vp_radmin(iirad) = pmin(char_lead_vp_radmin(iirad))
     char_lead_vp_radmax(iirad) = pmax(char_lead_vp_radmax(iirad))
     char_lead_vs_radmin(iirad) = pmin(char_lead_vs_radmin(iirad))
     char_lead_vs_radmax(iirad) = pmax(char_lead_vs_radmax(iirad))
     if (lpr) &
       write(60009,12) spher_radii(iirad),char_lead_vp_radmin(iirad),&
                       char_lead_vp_radmax(iirad),char_lead_vs_radmin(iirad),&
                       char_lead_vs_radmax(iirad)
  enddo
  if (lpr) close(60009)
  deallocate(char_lead_vp_radmin)
  deallocate(char_lead_vp_radmax)
  deallocate(char_lead_vs_radmin)
  deallocate(char_lead_vs_radmax)

! save coarsening levels 
  if (lpr) &
       open(unit=60009,file=infopath(1:lfinfo)//'/char_lead_time_coars.dat')
  do idisc=1,ndisc
     char_lead_vp_coarsmin(idisc)=pmin(char_lead_vp_coarsmin(idisc))
     char_lead_vp_coarsmax(idisc)=pmax(char_lead_vp_coarsmax(idisc))
     char_lead_vs_coarsmin(idisc)=pmin(char_lead_vs_coarsmin(idisc))
     char_lead_vs_coarsmax(idisc)=pmax(char_lead_vs_coarsmax(idisc))
     if (lpr .and. char_lead_vp_coarsmin(idisc) < router .and. & 
                   char_lead_vs_coarsmin(idisc) < router .and. & 
                   char_lead_vp_coarsmax(idisc) > zero .and. & 
                   char_lead_vs_coarsmax(idisc) > zero ) then
        ncoars= ncoars+1
        write(60009,12) rad_coars(idisc),char_lead_vp_coarsmin(idisc),&
                    char_lead_vp_coarsmax(idisc),char_lead_vs_coarsmax(idisc),&
                    char_lead_vs_coarsmax(idisc)
     endif
  enddo
  if (lpr) close(60009)
  deallocate(char_lead_vp_coarsmin)
  deallocate(char_lead_vp_coarsmax)
  deallocate(char_lead_vs_coarsmin)
  deallocate(char_lead_vs_coarsmax)

! save central cube per processor
  open(unit=60009,file=infopath(1:lfinfo)//&
                            '/char_lead_time_central.dat'//appmynum)
  do icub=1,ncub
      write(60009,13) scentral(icub),zcentral(icub), &
                char_lead_vp_centralmin(icub),char_lead_vp_centralmax(icub),&
                char_lead_vs_centralmax(icub),char_lead_vs_centralmax(icub)
  enddo
  close(60009)
  deallocate(char_lead_vp_centralmin)
  deallocate(char_lead_vp_centralmax)
  deallocate(char_lead_vs_centralmin)
  deallocate(char_lead_vs_centralmax)
  deallocate(scentral,zcentral)

! save basic parameters to extra file for post-processing of char. lead times
  if (lpr) then
     open(unit=60009,file=infopath(1:lfinfo)//'/char_lead_time_params.dat')
     write(60009,14)'period',period
     write(60009,14)'timestep',deltat
     write(60009,14)'aveGLL_lambda',pts_wavelngth*dble(npol)
     write(60009,14)'courant',courant
     write(60009,15)'ncube',ncub
     write(60009,15)'ncoarse',ncoars
     write(60009,15)'nelem',nelem
     write(60009,15)'nproc',nproc
     write(60009,16)'rmin',rmin
     close(60009)

     open(unit=60009,file=infopath(1:lfinfo)//'/char_lead_time_paramsNEW.dat')
     write(60009,*)period
     write(60009,*)deltat
     write(60009,*)pts_wavelngth*dble(npol)
     write(60009,*)courant
     write(60009,*)ncub
     write(60009,*)ncoars
     write(60009,*)nelem
     write(60009,*)nproc
     write(60009,*)rmin
     close(60009)

  endif

11 format(2(1pe15.5))
12 format(5(1pe13.3))
13 format(2(1pe12.3),4(1pe11.2))
14 format(a15,f10.4)
15 format(a15,i10)
16 format(a15,1pe13.3)
!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

end subroutine compute_numerical_resolution
!=============================================================================

!-----------------------------------------------------------------------------
subroutine model_output(rho,lambda,mu)

include "mesh_params.h"

double precision, intent(in) :: rho(0:npol,0:npol,nelem)
double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
double precision, intent(in) :: mu(0:npol,0:npol,nelem)
double precision             :: s,z,r,theta,vp,vs,ro,pts_per_km
integer                      :: iel,ipol,jpol,idom,iidom

  if (save_large_tests) then
     open(unit=5454,file=infopath(1:lfinfo)//&
          '/backgroundmodel_solver.dat'//appmynum)
     do iel=1,nelem
        do ipol=0,npol
           do jpol=0,npol
              call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
              write(5454,10)s,z,sqrt((lambda(ipol,jpol,iel)+ &
                                      two*mu(ipol,jpol,iel))/&
                   rho(ipol,jpol,iel)), & 
                   sqrt(mu(ipol,jpol,iel)/rho(ipol,jpol,iel)),&
                   rho(ipol,jpol,iel)
           enddo
        enddo
     enddo
     close(5454)
  endif

10 format(2(e15.7),3(e15.7))

! output discontinuities to check if they are discretized correctly
! (superimposed on model output)
  open(unit=5454,file=infopath(1:lfinfo)//&
                      '/discontinuities_solver.dat'//appmynum)
  do idom=1,ndisc
     write(5454,*)discont(idom)
  enddo
  close(5454)

  open(unit=5454,file=infopath(1:lfinfo)//&
                      '/vert_elem_edges_solver.dat'//appmynum)
  do iel=1,nelem
     call compute_coordinates(s,z,r,theta,iel,0,npol)
! only need to output vertical profile, e.g. axis
     if (s==zero .and. z>=zero) then
        write(5454,12)r,sqrt((lambda(0,npol,iel)+ &
             two*mu(0,npol,iel))/rho(0,npol,iel)), & 
             sqrt(mu(0,npol,iel)/rho(0,npol,iel)),rho(0,npol,iel)

        call compute_coordinates(s,z,r,theta,iel,0,0)
        write(5454,12)r,sqrt((lambda(0,0,iel)+ &
             two*mu(0,0,iel))/rho(0,0,iel)), & 
             sqrt(mu(0,0,iel)/rho(0,0,iel)),rho(0,0,iel)
     endif
  enddo
  close(5454)

12 format(4(1pe12.4))

! output radial model with km-scale resolution and points on discontinuities
  if (do_mesh_tests .or. (dump_wavefields .and. mynum == 0)) then
     pts_per_km = 2.
   
     if (lpr) write(6,13)pts_per_km
13 format('     fine-scale radial model output; points per km:',f7.2)

     open(unit=5454,file=infopath(1:lfinfo)//&
          '/backgroundmodel_kmscale.dat'//appmynum)

     write(5454,*) ceiling(router/1000.*pts_per_km)+1

     do iel=1,ceiling(router/1000.*pts_per_km)+1
        r = router-(real(iel)-1)*1000./pts_per_km
        idom = 1000
        do iidom=1,ndisc-1
           if (r<=discont(iidom) .and. r> discont(iidom+1) ) then 
              if (idom == 1000) then 
                 idom=iidom
              else
                 write(6,*)procstrg,'PROBLEM in fine-scale model output!'
                 write(6,*)procstrg,'Found multiple domains for radius=',r
                 write(6,*)procstrg,'Domain 1:',idom
                 write(6,*)procstrg,'Domain 2:',iidom
                 stop
              endif
           endif
        enddo
        if (r <= discont(ndisc)) then 
           if (idom == 1000) then 
              idom=ndisc
           else
              write(6,*)procstrg,'PROBLEM in fine-scale model output!'
              write(6,*)procstrg,'Found multiple domains for radius=',r
              write(6,*)procstrg,'Domain 1:',idom
              write(6,*)procstrg,'Domain 2:',ndisc
              stop
           endif
        endif
        if (idom==1000) then 
           write(6,*)procstrg,'PROBLEM in fine-scale model output!'
           write(6,*)procstrg,'Couldn"t find domain for radius=',r
           stop
        endif
        vp = velocity(r,'v_p',idom,bkgrdmodel,lfbkgrdmodel)
        vs = velocity(r,'v_s',idom,bkgrdmodel,lfbkgrdmodel)
        ro = velocity(r,'rho',idom,bkgrdmodel,lfbkgrdmodel)
        write(5454,12)r,vp,vs,ro
     enddo
  endif
end subroutine model_output
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal(x,y,z,u1,elems,filename)
implicit none
real(4), dimension(1:elems*4), intent(in) :: x,y,z,u1
integer  :: i,t,elems
integer, dimension(1:elems*5) :: cell
integer, dimension(1:elems) :: cell_type
character (len=200) :: filename
character (len=50) :: ss; !stream
!points structure

do i=5, elems*5, 5
   cell(i-4) = 4
enddo

t=0

do i=5,elems*5,5
   t = t + 4
   cell(i-3) = t - 4
   cell(i-2) = t - 3
   cell(i-1) = t - 2
   cell(i) = t - 1
enddo

cell_type = 9

open(100, file=trim(filename)//'.vtk', access='stream', status='replace', &
     convert='big_endian')
write(100) '# vtk DataFile Version 4.0'//char(10)
write(100) 'mittico'//char(10)
write(100) 'BINARY'//char(10)
write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
write(ss,fmt='(A6,I10,A5)') 'POINTS',elems*4,'float'
write(100) ss//char(10)

!points
write(100) (x(i),y(i),z(i),i=1,elems*4)
write(100) char(10)

!cell topology
write(ss,fmt='(A5,2I10)') 'CELLS',elems,elems*5
write(100) char(10)//ss//char(10)
write(100) cell
write(100) char(10)

!cell type
write(ss,fmt='(A10,2I10)') 'CELL_TYPES',elems
write(100) char(10)//ss//char(10)
write(100) cell_type
write(100) char(10)

!data
write(ss,fmt='(A10,I10)') 'POINT_DATA',elems*4
write(100) char(10)//ss//char(10)
write(100) 'SCALARS data float 1'//char(10)
write(100) 'LOOKUP_TABLE default'//char(10) !color table?
write(100) u1

!do i=1,elems*4
!   write(100) u1(i), u1(i) * 2.
!enddo

close(100)
write(6,*)'...saved ',trim(filename)//'.vtk'

end subroutine write_VTK_bin_scal
!-----------------------------------------------------------------------------

!=============================================================================
subroutine plot_model_vtk(rho,lambda,mu, xi_ani, phi_ani, eta_ani, &
                          fa_ani_theta, fa_ani_phi)

double precision, dimension(0:npol,0:npol,nelem), intent(in) :: rho 
double precision, dimension(0:npol,0:npol,nelem), intent(in) :: lambda,mu
double precision, dimension(0:npol,0:npol,nelem), intent(in), optional :: &
                                        xi_ani, phi_ani, eta_ani
double precision, dimension(0:npol,0:npol,nelem), intent(in), optional :: &
                                        fa_ani_theta, fa_ani_phi

real, dimension(:), allocatable :: vp1,vs1,rho1
real, dimension(:), allocatable :: vpv1,vsv1,eta1
real, dimension(:), allocatable :: xi1,phi1
real, dimension(:), allocatable :: fa_ani_theta1, fa_ani_phi1
character(len=200) :: fname
integer :: npts_vtk,ct,iel,i
real, allocatable ::  x(:),y(:),z0(:)
logical :: plot_ani

npts_vtk = nelem * 4

allocate(vp1(npts_vtk),vs1(npts_vtk),rho1(npts_vtk))
if (present(xi_ani) .and. present(phi_ani) .and.  present(eta_ani) &
    .and. present(fa_ani_theta) .and. present(fa_ani_phi))then
    allocate(vpv1(npts_vtk),vsv1(npts_vtk),eta1(npts_vtk))
    allocate(xi1(npts_vtk),phi1(npts_vtk))
    allocate(fa_ani_theta1(npts_vtk), fa_ani_phi1(npts_vtk))
    plot_ani = .true.
else
    plot_ani = .false.
endif

allocate(x(npts_vtk),y(npts_vtk),z0(npts_vtk))

z0 = 0.d0
ct = 0
do iel=1, nelem

   x(ct+1) = scoord(0,0,iel)
   x(ct+2) = scoord(npol,0,iel)
   x(ct+3) = scoord(npol,npol,iel)
   x(ct+4) = scoord(0,npol,iel)
   y(ct+1) = zcoord(0,0,iel)
   y(ct+2) = zcoord(npol,0,iel)
   y(ct+3) = zcoord(npol,npol,iel)
   y(ct+4) = zcoord(0,npol,iel)


   rho1(ct+1) = rho(0,0,iel)
   vp1(ct+1) = sqrt( (lambda(0,0,iel)+2.*mu(0,0,iel) ) / rho1(ct+1)  )
   vs1(ct+1) = sqrt( mu(0,0,iel)  / rho1(ct+1)  )

   if (plot_ani) then
     vpv1(ct+1) = sqrt(phi_ani(0,0,iel)) * vp1(ct+1)
     vsv1(ct+1) = vs1(ct+1) / sqrt(xi_ani(0,0,iel))
     xi1(ct+1) = xi_ani(0,0,iel)
     phi1(ct+1) = phi_ani(0,0,iel)
     eta1(ct+1) = eta_ani(0,0,iel)
     fa_ani_theta1(ct+1) = fa_ani_theta(0,0,iel)
     fa_ani_phi1(ct+1) = fa_ani_phi(0,0,iel)
   endif

   rho1(ct+2) = rho(npol,0,iel)
   vp1(ct+2) = sqrt( (lambda(npol,0,iel)+2.*mu(npol,0,iel) ) / rho1(ct+2)  )
   vs1(ct+2) = sqrt( mu(npol,0,iel)  / rho1(ct+2)  )

   if (plot_ani) then
     vpv1(ct+2) = sqrt(phi_ani(npol,0,iel)) * vp1(ct+2)
     vsv1(ct+2) = vs1(ct+2) / sqrt(xi_ani(npol,0,iel))
     xi1(ct+2) = xi_ani(npol,0,iel)
     phi1(ct+2) = phi_ani(npol,0,iel)
     eta1(ct+2) = eta_ani(npol,0,iel)
     fa_ani_theta1(ct+2) = fa_ani_theta(npol,0,iel)
     fa_ani_phi1(ct+2) = fa_ani_phi(npol,0,iel)
   endif

   rho1(ct+3) = rho(npol,npol,iel)
   vp1(ct+3) = sqrt( (lambda(npol,npol,iel)+2.*mu(npol,npol,iel) ) / rho1(ct+3)  )
   vs1(ct+3) = sqrt( mu(npol,npol,iel)  / rho1(ct+3)  )

   if (plot_ani) then
     vpv1(ct+3) = sqrt(phi_ani(npol,npol,iel)) * vp1(ct+3)
     vsv1(ct+3) = vs1(ct+3) / sqrt(xi_ani(npol,npol,iel))
     xi1(ct+3) = xi_ani(npol,npol,iel)
     phi1(ct+3) = phi_ani(npol,npol,iel)
     eta1(ct+3) = eta_ani(npol,npol,iel)
     fa_ani_theta1(ct+3) = fa_ani_theta(npol,npol,iel)
     fa_ani_phi1(ct+3) = fa_ani_phi(npol,npol,iel)
   endif

   rho1(ct+4) = rho(0,npol,iel)
   vp1(ct+4) = sqrt( (lambda(0,npol,iel)+2.*mu(0,npol,iel) ) / rho1(ct+4)  )
   vs1(ct+4) = sqrt( mu(0,npol,iel)  / rho1(ct+4)  )

   if (plot_ani) then
     vpv1(ct+4) = sqrt(phi_ani(0,npol,iel)) * vp1(ct+4)
     vsv1(ct+4) = vs1(ct+4) / sqrt(xi_ani(0,npol,iel))
     xi1(ct+4) = xi_ani(0,npol,iel)
     phi1(ct+4) = phi_ani(0,npol,iel)
     eta1(ct+4) = eta_ani(0,npol,iel)
     fa_ani_theta1(ct+4) = fa_ani_theta(0,npol,iel)
     fa_ani_phi1(ct+4) = fa_ani_phi(0,npol,iel)
   endif

   ct = ct + 4
enddo

if (plot_ani) then
   fname=trim(infopath(1:lfinfo)//'/model_vph_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vp1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_vsh_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vs1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_rho_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,rho1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_vpv_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vpv1,npts_vtk/4,fname)
   
   fname=trim(infopath(1:lfinfo)//'/model_vsv_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vsv1,npts_vtk/4,fname)
   
   fname=trim(infopath(1:lfinfo)//'/model_eta_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,eta1,npts_vtk/4,fname)
   
   fname=trim(infopath(1:lfinfo)//'/model_xi_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,xi1,npts_vtk/4,fname)
   
   fname=trim(infopath(1:lfinfo)//'/model_phi_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,phi1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_fa_theta_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,fa_ani_theta1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_fa_phi_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,fa_ani_phi1,npts_vtk/4,fname)

   deallocate(vp1,vs1,rho1, vsv1, vpv1, eta1, xi1, phi1, fa_ani_theta1, fa_ani_phi1)
else
   fname=trim(infopath(1:lfinfo)//'/model_vp_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vp1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_vs_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,vs1,npts_vtk/4,fname)

   fname=trim(infopath(1:lfinfo)//'/model_rho_'//appmynum)
   call write_VTK_bin_scal(x,y,z0,rho1,npts_vtk/4,fname)
   deallocate(vp1,vs1,rho1)
endif

end subroutine plot_model_vtk
!=============================================================================

!========================
end module get_model
!========================
