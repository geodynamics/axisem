!======================
module def_precomp_terms
!======================
!
! Read elastic information of the background model, define precomputable 
! matrices for mass, stiffness, boundary terms, pointwise derivatives.
! This is the quintessential module of the code...

use global_parameters
use data_mesh
use data_mesh_preloop
use data_spec
use data_matr
use data_source, ONLY : src_type
use data_proc

use get_mesh,    ONLY : compute_coordinates_mesh
use geom_transf
use utlity
use analytic_mapping

implicit none

public :: read_model_compute_terms
private
contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine read_model_compute_terms
!
! Wrapper routine to contain globally defined large matrices that are not 
! used in the time loop to this module (e.g. rho, lambda, mu).
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use get_model

include 'mesh_params.h'

double precision, dimension(:,:,:),allocatable :: rho,lambda,mu,massmat_kwts2
double precision, dimension(:,:,:),allocatable :: xi_ani, phi_ani, eta_ani

  if(lpr)write(6,*)'  ::::::::: BACKGROUND MODEL & PRECOMPUTED MATRICES:::::::'
  if(lpr)write(6,*)'  allocate elastic fields....';call flush(6)
  allocate(rho(0:npol,0:npol,1:nelem),massmat_kwts2(0:npol,0:npol,1:nelem))
  allocate(lambda(0:npol,0:npol,1:nelem),mu(0:npol,0:npol,1:nelem))
  if (ani_true) then
    allocate(xi_ani(0:npol,0:npol,1:nelem))
    allocate(phi_ani(0:npol,0:npol,1:nelem))
    allocate(eta_ani(0:npol,0:npol,1:nelem))
  endif

! load velocity/density model  (velocities in m/s, density in kg/m^3 )
  if(lpr)write(6,*)'  define background model....';call flush(6)
  if (ani_true) then
    if(lpr)write(6,*)'  background model is anisotropic....';call flush(6)
    call read_model_ani(rho, lambda, mu, xi_ani, phi_ani, eta_ani)
  else 
    call read_model(rho,lambda,mu)
  endif
  call model_output(rho,lambda,mu)

  if(lpr)write(6,*)'  compute Lagrange interpolant derivatives...'
  call lagrange_derivs

  if(lpr)write(6,*)'  define mass matrix....';call flush(6)
  call def_mass_matrix_k(rho,lambda,mu,massmat_kwts2)
   
  if (do_mesh_tests) then
    if(lpr)write(6,*)'  compute mass of the earth model....';call flush(6)
    call compute_mass_earth(rho)
  endif

  if(lpr)write(6,*)'  define precomputed matrices for pointwise derivatives...'
  call compute_pointwisederiv_matrices

  if (do_mesh_tests) then
    if(lpr)write(6,*)'  test pointwise derivatives & Laplacian in solid....'
    call test_pntwsdrvtvs_solid
    if(lpr)write(6,*)'  test pointwise derivatives & Laplacian in fluid....'
    if (have_fluid) call test_pntwsdrvtvs_fluid
  endif


  if(lpr)write(6,*)'  define solid stiffness terms....';call flush(6)
  if (ani_true) then
    call def_solid_stiffness_terms(lambda, mu, massmat_kwts2, xi_ani, phi_ani, eta_ani)
    deallocate(lambda,mu,xi_ani,phi_ani,eta_ani)
  else
    call def_solid_stiffness_terms(lambda,mu,massmat_kwts2)
    deallocate(lambda,mu)
  endif


  if (have_fluid) then
     if(lpr)write(6,*)'  define fluid stiffness terms....';call flush(6)
     call def_fluid_stiffness_terms(rho,massmat_kwts2)

     if(lpr)write(6,*)'  define solid-fluid boundary terms....';call flush(6)
     call def_solid_fluid_boundary_terms
  else
     M_w_fl=zero; M0_w_fl=zero
     M1chi_fl=zero; M2chi_fl=zero; M4chi_fl =zero    
     bdry_matr=zero
  endif

  if (lpr) write(6,*)'  ...defined all precomputed arrays';call flush(6)
  deallocate(rho,massmat_kwts2)

  if (lpr) write(6,*)'  ...deallocated unnecessary elastic arrays'

  if(lpr)write(6,*)'  ::::::: END BACKGROUND MODEL & PRECOMPUTED MATRICES:::::'
  call flush(6)

end subroutine read_model_compute_terms
!=============================================================================

!-----------------------------------------------------------------------------
subroutine lagrange_derivs
!
! Defines elemental arrays for the derivatives of Lagrange interpolating 
! functions either upon 
! Gauss-Lobatto-Legendre (all eta, and xi direction for non-axial elements) or 
! Gauss-Lobatto-Jacobi (0,1) points (axial xi direction):
! G1(i,j) = \partial_\xi ( \bar{l}_i(\xi_j) )  i.e. axial xi direction
! G2(i,j) = \partial_\eta ( l_i(\eta_j) )  i.e. all eta/non-ax xi directions 
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use splib, only : hn_jprime,lag_interp_deriv_wgl

include 'mesh_params.h'

double precision  :: df(0:npol),dg(0:npol)
integer           :: ishp,jpol
character(len=16) :: fmt1
logical           :: tensorwrong

! shp_deri_k only needed for the source and pointwise derivatives,
! otherwise (stiffness terms) we apply G0, G1, G2 and their transposes.

  shp_deri_k(:,:,:,:) = zero

! non-axial elements
        do ishp = 0, npol
           call hn_jprime(eta,ishp,npol,df)
           call hn_jprime(eta,ishp,npol,dg)
           do jpol = 0, npol
              shp_deri_k(ishp,jpol,1,1)=df(jpol)
              shp_deri_k(ishp,jpol,1,2)=dg(jpol)
           end do
        end do

! axial elements
        do ishp = 0, npol
           call lag_interp_deriv_wgl(df,xi_k,ishp,npol)
           call hn_jprime(eta,ishp,npol,dg)
           do jpol = 0, npol 
              shp_deri_k(ishp,jpol,2,1)=df(jpol)
              shp_deri_k(ishp,jpol,2,2)=dg(jpol)
           end do
        end do

! Define elemental Lagrange interpolant derivatives as needed for stiffness

! Derivative in z direction: \partial_\eta (l_j(\eta_q))
  G2 = shp_deri_k(0:npol,0:npol,1,2)
  G2T = transpose(G2)

! Derivative in s-direction: \partial_\xi (\bar{l}_i(\xi_p))
  G1 = shp_deri_k(0:npol,0:npol,2,1)
  G1T = transpose(G1)

! Axial vector
  G0 = shp_deri_k(0:npol,0,2,1)

! Simple test on Lagrange derivative tensor's antisymmetries...
  tensorwrong=.false.
  if ( mod(npol,2) == 0 ) then
     do jpol = 0,npol-1
        do ishp = 0,npol-1
           if (.not.reldiff_small(G2(ishp,jpol),-G2(npol-ishp,npol-jpol)))&
              then 
              write(6,*)procstrg, &
                        'PROBLEM: Lagrange deriv tensor in eta not symmetric!'
              write(6,*)procstrg,'polynomial order:',npol
              write(6,*)procstrg,'ishp,jpol,G2',ishp,jpol,G2(ishp,jpol)
              write(6,*)
              tensorwrong=.true.
           endif
        enddo
     enddo
     if (tensorwrong) then
        fmt1 = "(K(f10.3))"
        write(fmt1(2:2),'(i1.1)') npol+1
        write(6,*)'....Lagrange derivative tensor in eta:'
        do jpol=0,npol
           write(6,fmt1)(G2(ishp,jpol),ishp=0,npol)
        enddo
        stop
     endif
  endif

end subroutine lagrange_derivs
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_pointwisederiv_matrices
!
! The 4 necessary global matrices due to pointwise derivatives d/ds and d/dz:
! dzdeta/J, dzdxi/J, dsdeta/J, dsdxi/J (J: Jacobian).
! These are known during the time loop if the strain is computed on-the-fly.
! This is convenient to avoid recomputing these mapping derivatives at each 
! dumping stage and to avoid knowing the grid itself during the time loop
! (hence the additional 4 global variables in exchange for at least 2 for the
! mesh, i.e. only slightly more memory intensive).
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise
use data_io, ONLY : dump_type,need_fluid_displ

INTEGER          :: iel,inode,ipol,jpol
DOUBLE PRECISION :: dsdxi,dzdxi,dsdeta,dzdeta
double precision :: local_crd_nodes(8,2)

! fluid pointwise derivatives     
  allocate(DsDeta_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DzDeta_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DsDxi_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DzDxi_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(usz_fluid(0:npol,0:npol,1:nel_fluid,2))
  allocate(inv_s_fluid(0:npol,0:npol,1:nel_fluid))
  allocate(inv_s_rho_fluid(0:npol,0:npol,1:nel_fluid))

! solid pointwise derivatives
  allocate(DsDeta_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DzDeta_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DsDxi_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DzDxi_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(inv_s_solid(0:npol,0:npol,1:nel_solid))

! Solid region
  do iel = 1, nel_solid
    do inode = 1, 8
      call compute_coordinates_mesh(local_crd_nodes(inode,1),&
                                  local_crd_nodes(inode,2),ielsolid(iel),inode)
    end do
    if (.not. axis_solid(iel)) then ! non-axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDeta_over_J_sol(ipol,jpol,iel) = -dsdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDeta_over_J_sol(ipol,jpol,iel) = dzdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDxi_over_J_sol(ipol,jpol,iel) = dsdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDxi_over_J_sol(ipol,jpol,iel) = -dzdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         inv_s_solid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielsolid(iel))

       enddo
    enddo
   else ! axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                  xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         DsDeta_over_J_sol(ipol,jpol,iel) = -dsdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDeta_over_J_sol(ipol,jpol,iel) = dzdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDxi_over_J_sol(ipol,jpol,iel) = dsdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDxi_over_J_sol(ipol,jpol,iel) = -dzdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         if (ipol>0) then
            inv_s_solid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielsolid(iel)) 
         else
            inv_s_solid(ipol,jpol,iel) = one
         endif

        enddo
    enddo
    endif !axis
  enddo
 
  write(69,*)'Pointwise derivative precomputed terms in solid:'
  write(69,8)'  min/max DsDeta/J [1/m]:',minval(DsDeta_over_J_sol), &
       maxval(DsDeta_over_J_sol)
  write(69,8)'  min/max DzDeta/J [1/m]:',minval(DzDeta_over_J_sol), &
       maxval(DzDeta_over_J_sol)
  write(69,8)'  min/max DsDxi/J  [1/m]:',minval(DsDxi_over_J_sol), &
       maxval(DsDxi_over_J_sol)
  write(69,8)'  min/max DzDxi/J  [1/m]:',minval(DzDxi_over_J_sol), &
       maxval(DzDxi_over_J_sol)
  write(69,*)

! Fluid region
  do iel = 1, nel_fluid
    do inode = 1, 8
      call compute_coordinates_mesh(local_crd_nodes(inode,1),&
                                  local_crd_nodes(inode,2),ielfluid(iel),inode)
    end do
    if (.not. axis_fluid(iel)) then ! non-axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDeta_over_J_flu(ipol,jpol,iel) = -dsdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDeta_over_J_flu(ipol,jpol,iel) = dzdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDxi_over_J_flu(ipol,jpol,iel) = dsdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDxi_over_J_flu(ipol,jpol,iel) = -dzdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         
         inv_s_rho_fluid(ipol,jpol,iel) = inv_rho_fluid(ipol,jpol,iel)/&
                                          scoord(ipol,jpol,ielfluid(iel))

         inv_s_fluid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielfluid(iel))  
      enddo
    enddo
   else ! axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                  xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))

         DsDeta_over_J_flu(ipol,jpol,iel) = -dsdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDeta_over_J_flu(ipol,jpol,iel) = dzdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDxi_over_J_flu(ipol,jpol,iel) = dsdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDxi_over_J_flu(ipol,jpol,iel) = -dzdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))

         if (ipol>0) then
            inv_s_rho_fluid(ipol,jpol,iel) = inv_rho_fluid(ipol,jpol,iel)/&
                                             scoord(ipol,jpol,ielfluid(iel))
            inv_s_fluid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielfluid(iel))
         else
            inv_s_rho_fluid(ipol,jpol,iel) = inv_rho_fluid(ipol,jpol,iel)
            inv_s_fluid(ipol,jpol,iel) = one
         endif

       enddo
    enddo
    endif !axis
  enddo
  
  write(69,*)'Pointwise derivative precomputed terms in fluid:'
  write(69,8)'  min/max DsDeta/J [1/m]:',minval(DsDeta_over_J_flu),&
       maxval(DsDeta_over_J_flu)
  write(69,8)'  min/max DzDeta/J [1/m]:',minval(DzDeta_over_J_flu),&
       maxval(DzDeta_over_J_flu)
  write(69,8)'  min/max DsDxi/J  [1/m]:',minval(DsDxi_over_J_flu),&
       maxval(DsDxi_over_J_flu)
  write(69,8)'  min/max DzDxi/J  [1/m]:',minval(DzDxi_over_J_flu),&
       maxval(DzDxi_over_J_flu)
  write(69,*)

8 format(a25,2(1pe14.4))

! Prefactor for quadrupole phi-comp of fluid displacement
  if (src_type(1)=='monopole') inv_s_rho_fluid = zero
  if (src_type(1)=='quadpole') inv_s_rho_fluid = two*inv_s_rho_fluid

end subroutine compute_pointwisederiv_matrices
!--------------------------------------------------------------------------
!=============================================================================

!-----------------------------------------------------------------------------
subroutine test_pntwsdrvtvs_solid
!
! Test pointwise derivatives & axisymmetric Laplacian in solid region
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_io
use pointwise_derivatives

include 'mesh_params.h'

real(kind=realkind),allocatable :: tmpsolfield(:,:,:)
real(kind=realkind),allocatable :: tmpsolfieldcomp(:,:,:,:)
real(kind=realkind),allocatable :: tmpsolfielddiff(:,:,:,:)
real(kind=realkind)             :: meandiff(2)
double precision,allocatable    :: elderiv(:,:)
double precision                :: s,z,r,theta
integer                         :: iel,ipol,jpol,iarr(3)
character(len=16)               :: fmt1

  allocate(tmpsolfield(0:npol,0:npol,1:nel_solid))
  allocate(tmpsolfieldcomp(0:npol,0:npol,1:nel_solid,1:3))
  allocate(tmpsolfielddiff(0:npol,0:npol,1:nel_solid,1:3))
  allocate(elderiv(0:npol,0:npol))
  
! Test derivatives: define scalar field inside solid
  do iel=1,nel_solid
     do jpol=0,npol
        do ipol=0,npol
! make tmpsolfield equal to z s^3 + s z^3
           call compute_coordinates(s,z,r,theta,ielsolid(iel),ipol,jpol)
           tmpsolfield(ipol,jpol,iel)= z*s**3+s*z**3 + (s+z)*router
        enddo
     enddo
  enddo

  call axisym_laplacian_solid(tmpsolfield,tmpsolfieldcomp(:,:,:,1:2))

  meandiff=zero

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  open(unit=34,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_sol_num.dat'//appmynum)
  open(unit=36,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_reldiff.dat'//appmynum)
  do iel=1,nel_solid
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielsolid(iel),ipol,jpol)

           tmpsolfielddiff(ipol,jpol,iel,1)=absreldiff(&
                                   tmpsolfieldcomp(ipol,jpol,iel,1), &
                                   real(3.d0*z*s**2+z**3+router,kind=realkind))
           tmpsolfielddiff(ipol,jpol,iel,2)=absreldiff(&
                                   tmpsolfieldcomp(ipol,jpol,iel,2), &
                                   real(3.d0*s*z**2+s**3+router,kind=realkind))

           meandiff(1)=meandiff(1)+tmpsolfielddiff(ipol,jpol,iel,1)
           meandiff(2)=meandiff(2)+tmpsolfielddiff(ipol,jpol,iel,2)

           if (ipol==npol/2 .and. jpol==npol/2) then
              write(34,12)s,z,tmpsolfieldcomp(npol/2,npol/2,iel,1), &
                   tmpsolfieldcomp(npol/2,npol/2,iel,2)
              write(36,13)r/1000.,theta*180./pi, &
                          absreldiff(tmpsolfieldcomp(npol/2,npol/2,iel,1), &
                          real(3.d0*z*s**2+z**3+router,kind=realkind)),&
                          absreldiff(tmpsolfieldcomp(npol/2,npol/2,iel,2), &
                          real(3.d0*s*z**2+s**3+router,kind=realkind))
           endif
        enddo
     enddo
  enddo
  close(34); close(36)
!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  write(69,*)
  write(69,*)'/_/_/_/_/_/_/SOLID pointwise deriv with f=zs^3+sz^3/_/_/_/_/_/_/'
  write(69,9)'  mean error df/ds          :',meandiff(1)/ &
                                            real((npol)**2*nel_solid)
  write(69,8)'  min/max error df/ds       :',minval(tmpsolfielddiff(:,:,:,1)),&
                                            maxval(tmpsolfielddiff(:,:,:,1))
  iarr = maxloc(tmpsolfielddiff(:,:,:,1))
  write(69,8)'  r[m],theta max error df/ds:', &
                     rcoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3))), &
                     thetacoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3)))*180./pi
  write(69,7)'  elem type,coarsening ,axis :',eltype(ielsolid(iarr(3))), &
                                             coarsing(ielsolid(iarr(3))), &
                                             axis_solid(iarr(3))

  fmt1 = "(K(1pe13.4))"
  write(fmt1(2:2),'(i1.1)') npol+1
  write(69,*)' Analytical derivatives df/ds across that element (-->xi):'
  do jpol=0,npol
     do ipol=0,npol
        call compute_coordinates(s,z,r,theta,ielsolid(iarr(3)),ipol,jpol)
        elderiv(ipol,jpol)=3.d0*z*s**2+z**3+router
     enddo
     write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' Numerical derivatives df/ds across that element (-->xi):'
  do jpol=0,npol
     write(69,fmt1)(tmpsolfieldcomp(ipol,jpol,iarr(3),1),ipol=0,npol)
  enddo
  write(69,*)

  write(69,9)'  mean error df/dz          :',meandiff(2)/ &
                                             real((npol)**2*nel_solid)
  write(69,8)'  min/max error df/dz       :',minval(tmpsolfielddiff(:,:,:,2)),&
                                             maxval(tmpsolfielddiff(:,:,:,2))
  iarr = maxloc(tmpsolfielddiff(:,:,:,2))
  write(69,8)'  r[m],theta max error df/dz:', &
                      rcoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3))), &
                      thetacoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3)))*180./pi
  write(69,7)'  elem type,coarsening,axis :',eltype(ielsolid(iarr(3))), &
                                            coarsing(ielsolid(iarr(3))), &
                                            axis_solid(iarr(3))
  write(69,*)' Analytical derivatives df/dz across that element (-->xi):'
  do jpol=0,npol
     do ipol=0,npol
        call compute_coordinates(s,z,r,theta,ielsolid(iarr(3)),ipol,jpol)
                                 elderiv(ipol,jpol)=3.d0*s*z**2+s**3+router
     enddo
     write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' Numerical derivatives df/dz across that element (-->xi):'
  do jpol=0,npol
     write(69,fmt1)(tmpsolfieldcomp(ipol,jpol,iarr(3),2),ipol=0,npol)
  enddo
  write(69,*)'/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'

  write(69,*)
7  format(a30,a10,2(l4))
8  format(a30,2(1pe14.4))
9  format(a30,1pe14.4)
12 format(4(1pe16.7))
13 format(2(1pe15.5),2(1pe16.7))

  deallocate(tmpsolfield)
  deallocate(tmpsolfieldcomp)
  deallocate(tmpsolfielddiff)
  deallocate(elderiv)

end subroutine test_pntwsdrvtvs_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine test_pntwsdrvtvs_fluid
!
! Test pointwise derivatives & axisymmetric Laplacian in fluid
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_io
use pointwise_derivatives

include 'mesh_params.h'

real(kind=realkind),allocatable :: tmpflufield(:,:,:)
real(kind=realkind),allocatable :: tmpflufieldcomp(:,:,:,:)
real(kind=realkind),allocatable :: tmpflufielddiff(:,:,:,:)
real(kind=realkind)             :: meandiff(2)
double precision,allocatable    :: elderiv(:,:)
double precision                :: s,z,r,theta
integer                         :: iel,ipol,jpol,iarr(3)
character(len=16)               :: fmt1

  allocate(tmpflufield(0:npol,0:npol,1:nel_fluid))
  allocate(tmpflufieldcomp(0:npol,0:npol,1:nel_fluid,1:3))
  allocate(tmpflufielddiff(0:npol,0:npol,1:nel_fluid,1:3))
  allocate(elderiv(0:npol,0:npol))
  
! Test derivatives: define scalar field inside fluid
  do iel=1,nel_fluid
     do jpol=0,npol
        do ipol=0,npol
! make tmpflufield equal to z s^3 + s z^3
           call compute_coordinates(s,z,r,theta,ielfluid(iel),ipol,jpol)
           tmpflufield(ipol,jpol,iel)= z*s**3+s*z**3 + (s+z)*router
        enddo
     enddo
  enddo

  call axisym_laplacian_fluid(tmpflufield,tmpflufieldcomp(:,:,:,1:2))

  meandiff=zero

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  open(unit=34,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_flu_num.dat'//appmynum)
  open(unit=36,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_reldiff.dat'//appmynum)
  do iel=1,nel_fluid
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielfluid(iel),ipol,jpol)

           tmpflufielddiff(ipol,jpol,iel,1)=absreldiff(&
                                   tmpflufieldcomp(ipol,jpol,iel,1), &
                                   real(3.d0*z*s**2+z**3+router,kind=realkind))
           tmpflufielddiff(ipol,jpol,iel,2)=absreldiff(&
                                   tmpflufieldcomp(ipol,jpol,iel,2), &
                                   real(3.d0*s*z**2+s**3+router,kind=realkind))

           meandiff(1)=meandiff(1)+tmpflufielddiff(ipol,jpol,iel,1)
           meandiff(2)=meandiff(2)+tmpflufielddiff(ipol,jpol,iel,2)

           if (ipol==npol/2 .and. jpol==npol/2) then
              write(34,12)s,z,tmpflufieldcomp(npol/2,npol/2,iel,1), &
                   tmpflufieldcomp(npol/2,npol/2,iel,2)
              write(36,13)r/1000.,theta*180./pi, &
                          absreldiff(tmpflufieldcomp(npol/2,npol/2,iel,1), &
                          real(3.d0*z*s**2+z**3+router,kind=realkind)),&
                          absreldiff(tmpflufieldcomp(npol/2,npol/2,iel,2), &
                          real(3.d0*s*z**2+s**3+router,kind=realkind))
           endif
        enddo
     enddo
  enddo
  close(34); close(36)
!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  write(69,*)
  write(69,*)'/_/_/_/_/_/_/FLUID pointwise deriv with f=zs^3+sz^3/_/_/_/_/_/_/'
  write(69,9)'  mean error df/ds          :',meandiff(1)/ &
                                            real((npol)**2*nel_fluid)
  write(69,8)'  min/max error df/ds       :',minval(tmpflufielddiff(:,:,:,1)),&
                                            maxval(tmpflufielddiff(:,:,:,1))
  iarr = maxloc(tmpflufielddiff(:,:,:,1))
  write(69,8)'  r[m],theta max error df/ds:', &
                     rcoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3))), &
                     thetacoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3)))*180./pi
  write(69,7)'  elem type,coarsening ,axis :',eltype(ielfluid(iarr(3))), &
                                             coarsing(ielfluid(iarr(3))), &
                                             axis_fluid(iarr(3))

  fmt1 = "(K(1pe13.4))"
  write(fmt1(2:2),'(i1.1)') npol+1
  write(69,*)' Analytical derivatives df/ds across that element (-->xi):'
  do jpol=0,npol
     do ipol=0,npol
        call compute_coordinates(s,z,r,theta,ielfluid(iarr(3)),ipol,jpol)
        elderiv(ipol,jpol)=3.d0*z*s**2+z**3+router
     enddo
     write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' Numerical derivatives df/ds across that element (-->xi):'
  do jpol=0,npol
     write(69,fmt1)(tmpflufieldcomp(ipol,jpol,iarr(3),1),ipol=0,npol)
  enddo
  write(69,*)

  write(69,9)'  mean error df/dz          :',meandiff(2)/ &
                                             real((npol)**2*nel_fluid)
  write(69,8)'  min/max error df/dz       :',minval(tmpflufielddiff(:,:,:,2)),&
                                             maxval(tmpflufielddiff(:,:,:,2))
  iarr = maxloc(tmpflufielddiff(:,:,:,2))
  write(69,8)'  r[m],theta max error df/dz:', &
                      rcoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3))), &
                      thetacoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3)))*180./pi
  write(69,7)'  elem type,coarsening,axis :',eltype(ielfluid(iarr(3))), &
                                            coarsing(ielfluid(iarr(3))), &
                                            axis_fluid(iarr(3))
  write(69,*)' Analytical derivatives df/dz across that element (-->xi):'
  do jpol=0,npol
     do ipol=0,npol
        call compute_coordinates(s,z,r,theta,ielfluid(iarr(3)),ipol,jpol)
                                 elderiv(ipol,jpol)=3.d0*s*z**2+s**3+router
     enddo
     write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' Numerical derivatives df/dz across that element (-->xi):'
  do jpol=0,npol
     write(69,fmt1)(tmpflufieldcomp(ipol,jpol,iarr(3),2),ipol=0,npol)
  enddo
  write(69,*)'/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'

  write(69,*)
7  format(a30,a10,2(l4))
8  format(a30,2(1pe14.4))
9  format(a30,1pe14.4)
12 format(4(1pe16.7))
13 format(2(1pe15.5),2(1pe16.7))

  deallocate(tmpflufield)
  deallocate(tmpflufieldcomp)
  deallocate(tmpflufielddiff)
  deallocate(elderiv)

end subroutine test_pntwsdrvtvs_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine def_mass_matrix_k(rho,lambda,mu,massmat_kwts2)
!
! This routine computes and stores the coefficients of the diagonal 
! mass matrix, when a weighted Gauss-Lobatto quadrature for axial elements.
! It is built here with a factor of volume equal to s * ds * dz, as required
! by our approach.  we also define in this routine the mass matrix weighted
! by 1/s^2, as required by some components of the Laplacian of the
! fields. Note the special contribution arising in the case of an
! axial element. 
! massmat_k    : The actual mass term:    \sigma_I \sigma_J J^{IJ} s^{IJ} 
! massmat_kwts2: Paper 2, table 1, term A=\sigma_I \sigma_J J^{IJ} / s^{IJ}
! jacob: just defined locally for check on extrema further below....
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_io, ONLY : need_fluid_displ,dump_energy
use commun, only : comm2d
use data_pointwise, ONLY: inv_rho_fluid

include "mesh_params.h"

double precision, dimension(0:npol,0:npol,nelem),intent(in)  :: rho,lambda,mu
double precision, dimension(0:npol,0:npol,nelem),intent(out) :: massmat_kwts2

double precision, allocatable    :: massmat_k(:,:,:)   ! Mass matrix
double precision, allocatable    :: jacob (:,:,:)      ! jacobian array
real(kind=realkind), allocatable :: drdxi(:,:,:,:)     ! min/max derivs

double precision  :: local_crd_nodes(8,2),s,z,r,theta
integer           :: iel,inode,iarr(3),ipol,jpol
character(len=16) :: fmt1
double precision  :: dsdxi,dzdeta,dzdxi,dsdeta

  allocate(massmat_k(0:npol,0:npol,1:nelem),jacob(0:npol,0:npol,1:nelem))
  allocate(drdxi(0:npol,0:npol,1:nelem,1:4))
  allocate(inv_rho_fluid(0:npol,0:npol,1:nel_fluid))

  massmat_k(:,:,:) = zero; massmat_kwts2(:,:,:) = zero

!----------------------
  do iel = 1, nelem
!----------------------
     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1),&
             local_crd_nodes(inode,2),iel,inode)
     end do

! computing global arrays for the respective coordinate mapping derivatives
! ONLY needed for min/max write statements below!
     if ( axis(iel) ) then
        do ipol  = 0, npol
           do jpol = 0, npol
              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   xi_k(ipol),eta(jpol),local_crd_nodes,iel)
              drdxi(ipol,jpol,iel,1)=dsdxi; drdxi(ipol,jpol,iel,2)=dzdxi
              drdxi(ipol,jpol,iel,3)=dsdeta; drdxi(ipol,jpol,iel,4)=dzdeta
           enddo
        enddo
     else
        do ipol  = 0, npol
           do jpol = 0, npol
              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,iel)
              drdxi(ipol,jpol,iel,1)=dsdxi; drdxi(ipol,jpol,iel,2)=dzdxi
              drdxi(ipol,jpol,iel,3)=dsdeta; drdxi(ipol,jpol,iel,4)=dzdeta
           enddo
        enddo
     endif

! ::::::::::::::::non-axial elements::::::::::::::::
     if ( .not. axis(iel)) then
        do ipol  = 0, npol
           do jpol = 0, npol
              massmat_k(ipol,jpol,iel) = &
                   jacobian(eta(ipol),eta(jpol),local_crd_nodes,iel) &
                   *scoord(ipol,jpol,iel)*wt(ipol)*wt(jpol)
              massmat_kwts2(ipol,jpol,iel) =  &
                   jacobian(eta(ipol),eta(jpol),local_crd_nodes,iel) &
                   *scoord(ipol,jpol,iel)**(-1)*wt(ipol)*wt(jpol)
              jacob(ipol,jpol,iel) = jacobian(eta(ipol),eta(jpol),&
                                              local_crd_nodes,iel)
           end do
        end do

! ::::::::::::::::axial elements::::::::::::::::
     elseif (axis(iel)) then
        do ipol  = 0, npol ! Be careful here !!!!
           do jpol = 0, npol
              massmat_k(ipol,jpol,iel) = &
                    jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,iel) &
                    *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol),&
                    local_crd_nodes,iel)*wt_axial_k(ipol)*wt(jpol)
              jacob(ipol,jpol,iel) = jacobian(xi_k(ipol),eta(jpol),&
                                              local_crd_nodes,iel)
           end do
        end do
        do ipol = 1, npol
           do jpol = 0, npol
              massmat_kwts2(ipol,jpol,iel) = &
                   jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,iel) / &
                   ( scoord(ipol,jpol,iel) * ( one+xi_k(ipol) ) ) * &
                   wt_axial_k(ipol)*wt(jpol)
           end do
        end do
!TNM ADDITION NOV 2006: for ipol=0
           do jpol = 0, npol
              massmat_kwts2(0,jpol,iel) = &
                   jacobian(xi_k(0),eta(jpol),local_crd_nodes,iel)* &
                   (s_over_oneplusxi_axis(xi_k(0),eta(jpol),&
                   local_crd_nodes,iel))**(-1) * &
                   wt_axial_k(0)*wt(jpol)
           enddo
! END ADDITION

! check consistency of s/(1+xi) definitions...
! axial element but off-axis
        do ipol = 1, npol
           do jpol = 0, npol
              if ( .not. dblreldiff_small( scoord(ipol,jpol,iel)/ &
                                        ( one+xi_k(ipol) ), &
                                       s_over_oneplusxi_axis(xi_k(ipol), &
                                        eta(jpol),local_crd_nodes,iel)) ) then
                 write(6,*)procstrg,'PROBLEM: 2 definitions of s/(1+xi) differ'
                 write(6,*)procstrg,'scoord/(1+xi)=',scoord(ipol,jpol,iel)/ &
                                                     ( one+xi_k(ipol) )
                 write(6,*)procstrg,'s_over_onexi=', &
                                     s_over_oneplusxi_axis(xi_k(ipol), &
                                     eta(jpol),local_crd_nodes,iel)
                 write(6,*)procstrg,'iel,ipol,jpol:',iel,ipol,jpol
                 write(6,*)procstrg,'s,r,theta:',scoord(ipol,jpol,iel),& 
                           rcoord(ipol,jpol,iel),thetacoord(ipol,jpol,iel)
                 stop
              endif
           end do
        end do

     end if ! axial?

!----------------------
  end do
!----------------------

  if(lpr)write(6,*)'  solid mass matrix...'; call flush(6)
! Solid inverse mass term_____________________________________________
  do iel=1,nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
              inv_mass_rho(ipol,jpol,iel)= rho(ipol,jpol,ielsolid(iel))*&
                                           massmat_k(ipol,jpol,ielsolid(iel))
        end do
     end do
  enddo

! unassembled mass matrix in solid for energy.
   if (dump_energy) then
     allocate(unassem_mass_rho_solid(0:npol,0:npol,nel_solid))
     unassem_mass_rho_solid = inv_mass_rho
     if (src_type(1)=='dipole') &
                    unassem_mass_rho_solid = two * unassem_mass_rho_solid
   endif



! Exchange boundary information_______________________________________
  if(lpr)write(6,*)'  assemble solid mass matrix...'; call flush(6)
  call comm2d(inv_mass_rho,nel_solid,1,'solid')

  if(lpr)write(6,*)'  compute inverse solid mass matrix...'; call flush(6)
  do iel=1,nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           if ( inv_mass_rho(ipol,jpol,iel) /= zero) then 
              inv_mass_rho(ipol,jpol,iel) = one / inv_mass_rho(ipol,jpol,iel)
           else
              write(6,*)procstrg,'WARNING: solid mass term zero!', &
                         ipol,jpol,iel,ielsolid(iel)
              inv_mass_rho(ipol,jpol,iel) = one / rho(ipol,jpol,ielsolid(iel))
           endif
        end do
     end do
  enddo

  if (src_type(1)=='dipole') inv_mass_rho = half * inv_mass_rho

! Fluid inverse mass term_____________________________________________
  if(lpr)write(6,*)'  fluid mass matrix...'; call flush(6)
  do iel=1,nel_fluid
! check if fluid element is really fluid throughout
     if (maxval(mu(:,:,ielfluid(iel)))> zero) then
        call compute_coordinates(s,z,r,theta,ielfluid(iel), &
                                 int(npol/2),int(npol/2))
        write(6,*)
        write(6,*)procstrg,'!!!!!!!!!!!!!!!!  PROBLEM  !!!!!!!!!!!!!!!!!!!!!!!'
        write(6,*)procstrg,'Have a non-zero mu in the fluid region!'
        write(6,*)procstrg,'Value & location r[km],theta[deg]:', &
                  maxval(mu(:,:,ielfluid(iel))),r/1000.,theta*180./pi      
        call flush(6)
        stop
     endif

     do ipol = 0, npol
        do jpol = 0, npol
! since mu=0 only consider lambda here for the
! bulk modulus/incompressibility kappa = lambda + 2/3 mu
! in the ideal fluid that harbors the outer core...
           inv_mass_fluid(ipol,jpol,iel)=massmat_k(ipol,jpol,ielfluid(iel))/&
                                         lambda(ipol,jpol,ielfluid(iel))

! inverse density inside the fluid, needed to calculate displacement in fluid
           if (need_fluid_displ) &
              inv_rho_fluid(ipol,jpol,iel) = one / rho(ipol,jpol,ielfluid(iel))
        end do
     end do
  enddo


! unassembled mass matrix in fluid for the energy.
   if (dump_energy) then 
     allocate(unassem_mass_lam_fluid(0:npol,0:npol,nel_fluid))
     unassem_mass_lam_fluid = inv_mass_fluid
   endif


! Exchange boundary information_______________________________________
  if(lpr)write(6,*)'  assemble fluid mass matrix...'; call flush(6)
  call comm2d(inv_mass_fluid,nel_fluid,1,'fluid')

  if(lpr)write(6,*)'  compute inverse fluid mass matrix...'; call flush(6)
  do iel=1,nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           if (inv_mass_fluid(ipol,jpol,iel) /= zero) then
              inv_mass_fluid(ipol,jpol,iel) = one / inv_mass_fluid(ipol,jpol,iel)
           else
              write(6,*)procstrg,'WARNING: Fluid mass term zero!', &
                        ipol,jpol,iel,ielfluid(iel); call flush(6)
              inv_mass_fluid(ipol,jpol,iel)=lambda(ipol,jpol,ielfluid(iel))
           end if
        end do
     end do
  enddo




! In the remainder: document min/max values and locations for velocities, 
!    density, Jacobian, mass terms, GLL points. Lagrange derivatives etc.
  fmt1 = "(a18,K(f9.4))"
  write(fmt1(6:6),'(i1.1)') npol+1
  write(69,*)
  write(69,*)'-+-+-+-+-+-+-+-+-+-+ Integration weights +-+-+-+-+-+-+-+-+-+-+-+'
  write(69,fmt1)'GLJ sigma_I ax   :',(wt_axial_k(ipol),ipol=0,npol)  
  write(69,fmt1)'GLL sigma_J nonax:',(wt(ipol),ipol=0,npol)  

  fmt1 = "(a11,K(f9.4))"
  write(fmt1(6:6),'(i1.1)') npol+1
  write(69,*)
  write(69,*)'-+-+-+- Gauss-Lobatto-Legendre pts, Lagrange derivatives +-+-+-+'
  write(69,fmt1)'GLL eta  :',(eta(ipol),ipol=0,npol)
  write(69,*)
  write(69,*)' Lagrange derivatives eta: \partial_\eta (l_j(\eta_q)), --> li_i'
  do jpol=0,npol
     write(69,fmt1)'',(G2(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' transpose derivatives eta: --> eta_p'
  do jpol=0,npol
     write(69,fmt1)'',(G2T(ipol,jpol),ipol=0,npol)
  enddo

  write(69,*)
  write(69,*)'-+-+-+ Gauss-Lobatto-Jacobi (0,1) pts, Lagrange derivatives +-+-'
  write(69,fmt1)'GLJ xi_ax:',(xi_k(ipol),ipol=0,npol)
  write(69,*)
  write(69,*)' Lagrange derivatives G1 xi : \partial_\xi (l_i(\xi_p)), --> li_i'
  do jpol=0,npol
     write(69,fmt1)'',(G1(ipol,jpol),ipol=0,npol)
  enddo
  write(69,*)' transpose G1T derivatives xi --> xi_p'
  do jpol=0,npol
     write(69,fmt1)'',(G1T(ipol,jpol),ipol=0,npol)
  enddo  

  write(69,*)' Lagrange derivatives axial vector: \partial_\xi (l_i(\xi_0)) '
  write(69,fmt1)'',(G0(ipol),ipol=0,npol)
  write(69,*)' '

  write(69,*)'+-+-+-+-+-+- Coordinate mapping/derivative extrema +-+-+-+-+-+--'
  iarr = minloc(drdxi(:,:,:,1))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th dsdxi [m,m,deg]        :',minval(drdxi(:,:,:,1)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  iarr = maxloc(drdxi(:,:,:,1))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th dsdxi [m,m,deg]        :',maxval(drdxi(:,:,:,1)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

  iarr = minloc(drdxi(:,:,:,2))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th dzdxi [m,m,deg]        :',minval(drdxi(:,:,:,2)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  iarr = maxloc(drdxi(:,:,:,2))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th dzdxi [m,m,deg]        :',maxval(drdxi(:,:,:,2)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

  iarr = minloc(drdxi(:,:,:,3))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th dsdeta [m,m,deg]       :',minval(drdxi(:,:,:,3)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  iarr = maxloc(drdxi(:,:,:,3))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th dsdeta [m,m,deg]       :',maxval(drdxi(:,:,:,3)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

  iarr = minloc(drdxi(:,:,:,4))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th dzdeta [m,m,deg]       :',minval(drdxi(:,:,:,4)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  iarr = maxloc(drdxi(:,:,:,4))
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th dzdeta [m,m,deg]       :',maxval(drdxi(:,:,:,4)), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

  write(69,*)' '
  write(69,*)'+-+-+-+-+-+-+-+- Jacobian & mass term extrema -+-+-+-+-+-+-+-+-+'
  write(69,*)'   Jacobian  = dsdxi dzdeta - dsdeta dzdxi'
  write(69,*)'   (mass term)_{IJ} = Jacobian_{IJ} \sigma_I \sigma_J s_{IJ}))'
  iarr = minloc(jacob)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th Jacobian [m,m,deg]        :',minval(jacob), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  iarr = maxloc(jacob)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th Jacobian [m,m,deg]        :',maxval(jacob), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
  
  iarr = minloc(massmat_k)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th mass term [m^3,m,deg]     :',minval(massmat_k), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 

  iarr = maxloc(massmat_k)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th mass term [m^3,m,deg]     :',maxval(massmat_k), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 
  iarr = minloc(inv_mass_rho)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Min,r,th sol. invmass [1/kg,m,deg] :',minval(inv_mass_rho), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 
  iarr = maxloc(inv_mass_rho)
  theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
  write(69,9)'Max,r,th sol. invmass [1/kg,m,deg] :',maxval(inv_mass_rho), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 
  if (have_fluid) then
     iarr = minloc(inv_mass_fluid)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th flu. invmass [N/m^5,m,deg]:',minval(inv_mass_fluid), &
          rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 
     iarr = maxloc(inv_mass_fluid)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th flu. invmass [N/m^5,m,deg]:',maxval(inv_mass_fluid), &
                                rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta 
  endif


  write(69,*)' '
  write(69,*)'+-+-+-+-+-+-+-+-+-+-+- Elastic extrema -+-+-+-+-+-+-+-+-+-+-+-+-'
  iarr = minloc(sqrt((lambda+2*mu)/rho))
  write(69,8)'Min,r P-vel [m/s,m]   :',minval(dsqrt((lambda+2*mu)/rho)), &
                                     rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
  iarr = minloc(sqrt(mu/rho)) 
  write(69,8)'Min,r S-vel [m/s,m]   :',minval(dsqrt(mu/rho)), &
                                      rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
  iarr = maxloc(sqrt((lambda+2*mu)/rho)) 
  write(69,8)'Max,r P-vel [m/s,m]   :',maxval(dsqrt((lambda+2*mu)/rho)), &
                                      rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
  iarr = maxloc(sqrt(mu/rho))
  write(69,8)'Max,r S-vel [m/s,m]   :',maxval(dsqrt(mu/rho)), &
                                      rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
  iarr = minloc(rho)
  write(69,8)'Min,r rho [kg/m^3,m]  :',minval(rho), &
                                      rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
  iarr = maxloc(rho)
  write(69,8)'Max,r rho [kg/m^3,m]  :',maxval(rho), &
                                       rcoord(iarr(1)-1,iarr(2)-1,iarr(3))/1.d3
  write(69,8)'Min/max mu [N/m^2]    :',minval(mu),maxval(mu)
  write(69,8)'Min/max lambda [N/m^2]:',minval(lambda),maxval(lambda)

  write(69,*)'' 

8 format(a25,2(1pe12.4))
9 format(a37,3(1pe12.4))

  deallocate(massmat_k,jacob)
  deallocate(drdxi)

end subroutine def_mass_matrix_k
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_mass_earth(rho)
!
! A straight computation of the mass of the sphere and that of its 
! solid and fluid sub-volumes. This is the same as computing the volume 
! (see def_grid.f90), but with a multiplicative density factor, i.e. the 
! actual mass matrix. The comparison to the exact value is merely 
! done as an indicator and does not cause any error. One should keep an eye 
! on these values when generating any new kind of background model for which 
! one then needs to dig up the total mass...
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use def_grid, ONLY : massmatrix
use background_models, ONLY: velocity
use commun, ONLY : psum
use data_io, ONLY : infopath,lfinfo

include 'mesh_params.h'

double precision, intent(in)     :: rho(0:npol,0:npol,nelem)
integer                          :: iel,ipol,jpol,idom,iidom,idisc
integer                          :: count_solid,count_fluid,count_sic
double precision                 :: mass_glob,mass_solid,mass_fluid,mass_sic
double precision                 :: mass_glob_num,mass_solid_num
double precision                 :: mass_fluid_num,mass_sic_num
double precision                 :: mass_layer(ndisc)
double precision                 :: r,density,vs,dr
real(kind=realkind), allocatable :: massmat(:,:,:),massmat_solid(:,:,:)
real(kind=realkind), allocatable :: massmat_fluid(:,:,:)
!af pgf
character(len=100) :: modelstring

  allocate(massmat(0:npol,0:npol,1:nelem))
  allocate(massmat_solid(0:npol,0:npol,1:nel_solid))
  allocate(massmat_fluid(0:npol,0:npol,1:nel_fluid))
  
  mass_fluid=zero; mass_glob=zero; mass_solid=zero; mass_sic=zero
  count_fluid=0; count_solid = 0; count_sic = 0

!af
  modelstring=bkgrdmodel!(1:(lfbkgrdmodel))
! actual masses, calculated by summation over layers spaced 100 meters
  do iel=1,int(router/100.)+1
     r = router-(real(iel)-.5)*100.
     dr = (router-(real(iel)-1.)*100.)**3 - (router-(real(iel)*100.))**3
     idom=10*ndisc
     do iidom=1,ndisc-1
        if (r<=discont(iidom) .and. r> discont(iidom+1) ) then 
           idom=iidom
           exit
        endif
        if (idom==iidom) exit
     enddo
     if (r <= discont(ndisc)) idom=ndisc
     density=velocity(r,'rho',idom,modelstring,lfbkgrdmodel)
     mass_glob = mass_glob +  density*dr
     vs=velocity(r,'v_s',idom,modelstring,lfbkgrdmodel)
     if (vs < 0.1d0) then
        mass_fluid = mass_fluid + density*dr
     else
        mass_solid = mass_solid + density*dr
     endif
! Solid inner core only
     if (idom==ndisc) mass_sic = mass_sic + density*dr
  enddo
  mass_glob  = 4.d0/3.d0*pi*mass_glob
  mass_fluid = 4.d0/3.d0*pi*mass_fluid
  mass_solid = 4.d0/3.d0*pi*mass_solid
  mass_sic = 4.d0/3.d0*pi*mass_sic
  
  mass_glob_num=zero; mass_solid_num=zero; mass_fluid_num=zero
  mass_sic_num=zero

! numerically computed masses
  call massmatrix(massmat,nelem,'total')
  call massmatrix(massmat_solid,nel_solid,'solid')
  call massmatrix(massmat_fluid,nel_fluid,'fluid')

  do iel = 1, nelem
     do ipol = 0, npol
        do jpol = 0, npol
           mass_glob_num = mass_glob_num + &
                           rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)
! solid inner core only
           if (rcoord(npol/2,npol/2,iel)< discont(ndisc)) &
                mass_sic_num = mass_sic_num + & 
                               rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)

! compute mass for each layer between discontinuities
           do idisc = 1,ndisc-1
              if (rcoord(npol/2,npol/2,iel) <= discont(idisc) .and. &
                  rcoord(npol/2,npol/2,iel) >  discont(idisc+1) ) then
                 mass_layer(idisc) = mass_layer(idisc) + &
                                     rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)
              endif
           enddo

        end do
     end do
  end do
  mass_layer(ndisc) = mass_sic_num
  mass_glob_num=2.d0*pi*mass_glob_num
  mass_glob_num=psum(real(mass_glob_num,kind=realkind))
  mass_sic_num=2.d0*pi*mass_sic_num
  mass_sic_num=psum(real(mass_sic_num,kind=realkind)) 
  mass_layer=2.d0*pi*mass_layer
  do idisc = 1,ndisc-1
     mass_layer(idisc) = psum(real(mass_layer(idisc),kind=realkind))
  enddo

  do iel = 1, nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           mass_solid_num = mass_solid_num + rho(ipol,jpol,ielsolid(iel))* &
                                             massmat_solid(ipol,jpol,iel)
        end do
     end do
  end do
  mass_solid_num=2.d0*pi*mass_solid_num
  mass_solid_num=psum(real(mass_solid_num,kind=realkind))

  do iel = 1, nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           mass_fluid_num = mass_fluid_num + rho(ipol,jpol,ielfluid(iel))* &
                                             massmat_fluid(ipol,jpol,iel)
        end do
     end do
  end do
  mass_fluid_num=2.d0*pi*mass_fluid_num
  mass_fluid_num=psum(real(mass_fluid_num,kind=realkind))

  if (lpr) then
     write(6,*)'  Calculated masses for earth model: ',&
                                          bkgrdmodel(1:lfbkgrdmodel)
     write(6,10)'  Total mass (real,num,diff)    :',mass_glob,mass_glob_num,&
                                         abs(mass_glob-mass_glob_num)/mass_glob
     write(6,11)'  Sum of layers (real,num,diff) :',sum(mass_layer)
     write(6,10)'    Solid mass (real,num,diff)    :',mass_solid,& 
                       mass_solid_num,abs(mass_solid-mass_solid_num)/mass_solid
     if (have_fluid) then
        write(6,10)'    Fluid mass (real,num,diff)    :',mass_fluid,&
                    mass_fluid_num,abs(mass_fluid-mass_fluid_num)/mass_fluid
     endif
     write(6,10)'    Innercore mass (real,num,diff):',mass_sic, &
                       mass_sic_num,abs(mass_sic-mass_sic_num)/mass_sic
     write(6,*) '   Mass of the Earth             :   5.97 x 10^24 kg'
10   format(a35,2(1pe14.5),1pe12.2)
11   format(a35,1pe14.5)
     write(6,*)
  endif

! write out total masses for each layer
  if (lpr) then 
     open(unit=9876,file=infopath(1:lfinfo)//'/mass_kg_per_discont_layer.dat')
     do idisc=1,ndisc
        write(9876,*)discont(idisc),mass_layer(idisc)
     enddo
     close(9876)
  endif

  deallocate(massmat)
  deallocate(massmat_solid)
  deallocate(massmat_fluid)

end subroutine compute_mass_earth
!=============================================================================

!-----------------------------------------------------------------------------
subroutine def_solid_stiffness_terms(lambda, mu, massmat_kwts2, xi_ani, phi_ani, eta_ani)
!
! This routine is a merged version to minimize global work 
! array definitions. The terms alpha_wt_k etc. are now 
! merely elemental arrays, and defined on the fly when 
! computing the global, final precomputable matrices 
! for the solid stiffness term. The loop is over solid elements only.
! Tarje, Sept 2006.
!
! Adding optional arguments for anisotropy, MvD
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include "mesh_params.h"

double precision, dimension(0:npol,0:npol,nelem), intent(in) :: lambda,mu
double precision, dimension(0:npol,0:npol,nelem), intent(in) :: massmat_kwts2
double precision, dimension(0:npol,0:npol,nelem), intent(in), optional :: xi_ani, phi_ani, eta_ani

double precision :: local_crd_nodes(8,2)
integer          :: ielem,ipol,jpol,inode
double precision :: dsdxi,dzdeta,dzdxi,dsdeta

double precision :: alpha_wt_k(0:npol,0:npol)
double precision :: beta_wt_k(0:npol,0:npol)
double precision :: gamma_wt_k(0:npol,0:npol)
double precision :: delta_wt_k(0:npol,0:npol)
double precision :: epsil_wt_k(0:npol,0:npol)
double precision :: zeta_wt_k(0:npol,0:npol)

double precision :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
double precision :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
double precision :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
double precision :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)
double precision :: M_s_xi_wt_k(0:npol,0:npol)
double precision :: M_z_xi_wt_k(0:npol,0:npol)
double precision :: M_z_eta_wt_k(0:npol,0:npol)
double precision :: M_s_eta_wt_k(0:npol,0:npol)

! non-diagfact
double precision, allocatable :: non_diag_fact(:,:)

!-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->
! NON DIAG FAC------------- --<->-<->-<->-<->-<->-<->-<->-<->-<->-<->
!-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->
  allocate(non_diag_fact(0:npol,1:nel_solid))
  non_diag_fact(:,:) = zero
  do ielem = 1, nel_solid
     if ( axis_solid(ielem) ) then
        do inode = 1, 8
           call compute_coordinates_mesh(local_crd_nodes(inode,1),&
                local_crd_nodes(inode,2),ielsolid(ielem),inode)
        end do
        do jpol = 0,npol 
           non_diag_fact(jpol,ielem) = wt_axial_k(0)*wt(jpol)* &
              jacobian(xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))/&
              s_over_oneplusxi_axis(xi_k(0),eta(jpol),&
                                    local_crd_nodes,ielsolid(ielem))
        end do
     end if
  end do

!-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->
! SOLID STIFFNESS TERMS-------<->-<->-<->-<->-<->-<->-<->-<->-<->-<->
!-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->-<->

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do ielem=1,nel_solid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     alpha_wt_k(:,:) = zero; beta_wt_k(:,:)  = zero; gamma_wt_k(:,:) = zero
     delta_wt_k(:,:) = zero; epsil_wt_k(:,:) = zero; zeta_wt_k(:,:)  = zero
     Ms_z_eta_s_xi_wt_k(:,:) = zero; Ms_z_eta_s_eta_wt_k(:,:) = zero
     Ms_z_xi_s_eta_wt_k(:,:) = zero; Ms_z_xi_s_xi_wt_k(:,:) = zero
     M_s_xi_wt_k(:,:) = zero; M_z_xi_wt_k(:,:) = zero
     M_z_eta_wt_k(:,:) = zero; M_s_eta_wt_k(:,:) = zero

     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1),&
             local_crd_nodes(inode,2),ielsolid(ielem),inode)
     end do

! ::::::::::::::::non-axial elements::::::::::::::::
     if ( .not. axis_solid(ielem) ) then
        do ipol  = 0, npol
           do jpol = 0, npol

              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielsolid(ielem))

              alpha_wt_k(ipol,jpol) = alpha(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              beta_wt_k(ipol,jpol) = beta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              gamma_wt_k(ipol,jpol) = gamma1(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              delta_wt_k(ipol,jpol) = delta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              epsil_wt_k(ipol,jpol) = epsilon1(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              zeta_wt_k(ipol,jpol) = zeta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_eta_s_xi_wt_k(ipol,jpol) = Ms_z_eta_s_xi(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_eta_s_eta_wt_k(ipol,jpol) = Ms_z_eta_s_eta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_xi_s_eta_wt_k(ipol,jpol) = Ms_z_xi_s_eta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_xi_s_xi_wt_k(ipol,jpol) = Ms_z_xi_s_xi(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              M_s_xi_wt_k(ipol,jpol)  =  dsdxi*wt(ipol)*wt(jpol)
              M_z_xi_wt_k(ipol,jpol)  = -dzdxi*wt(ipol)*wt(jpol)
              M_z_eta_wt_k(ipol,jpol) =  dzdeta*wt(ipol)*wt(jpol)
              M_s_eta_wt_k(ipol,jpol) = -dsdeta*wt(ipol)*wt(jpol)
           enddo
        enddo

! ::::::::::::::::axial elements::::::::::::::::
     elseif ( axis_solid(ielem) ) then
        do jpol  = 0, npol
           do ipol = 0, npol
              alpha_wt_k(ipol,jpol) = alphak(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              beta_wt_k(ipol,jpol) = betak(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              gamma_wt_k(ipol,jpol) = gammak(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              delta_wt_k(ipol,jpol) = deltak(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              epsil_wt_k(ipol,jpol) = epsilonk(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              zeta_wt_k(ipol,jpol) = zetak(xi_k(ipol),&
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)

              Ms_z_eta_s_xi_wt_k(ipol,jpol) = Ms_z_eta_s_xi_k(xi_k(ipol), & 
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol) 
              Ms_z_eta_s_eta_wt_k(ipol,jpol) = Ms_z_eta_s_eta_k(xi_k(ipol),& 
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol) 
              Ms_z_xi_s_eta_wt_k(ipol,jpol) = Ms_z_xi_s_eta_k(xi_k(ipol), & 
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol) 

              Ms_z_xi_s_xi_wt_k(ipol,jpol) = Ms_z_xi_s_xi_k(xi_k(ipol), &  
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol) 

              if (ipol>0) then
                 call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                      xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(ielem))
                 M_s_xi_wt_k(ipol,jpol)  = dsdxi*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_z_xi_wt_k(ipol,jpol)  = -dzdxi*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_z_eta_wt_k(ipol,jpol) = dzdeta*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_s_eta_wt_k(ipol,jpol) = -dsdeta*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
              else
                 M_s_xi_wt_k(ipol,jpol) = zero
                 M_z_xi_wt_k(ipol,jpol) = zero
                 M_z_eta_wt_k(ipol,jpol) = zero
                 M_s_eta_wt_k(ipol,jpol) = zero
              endif
           enddo
        enddo
     endif ! axis_solid or not

!    +++++++++++++++++++++
     do jpol=0,npol
!    +++++++++++++++++++++
       select case (src_type(1))
       case ('monopole')
         if (ani_true) then
            call compute_monopole_stiff_terms_ani(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,xi_ani,phi_ani,eta_ani, &
                                       massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)
         else
            call compute_monopole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)
         endif

       case('dipole')
         if (ani_true) then
            write(6,*) 'ERROR: Anisotropy not yet implemented for Dipole sources!'
            stop
         else
            call compute_dipole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)
         endif

       case('quadpole') 
         if (ani_true) then
           write(6,*) 'ERROR: Anisotropy not yet implemented for Qudrupole sources!'
           stop
         else
            call compute_quadrupole_stiff_terms(ielem,jpol,local_crd_nodes,&
                                       lambda,mu,massmat_kwts2,&
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)
         endif

      end select

!    +++++++++++++++++++++
     enddo  ! jpol
!    +++++++++++++++++++++

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  enddo ! solid elements
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  write(69,*)' '
  write(69,*)'Min/max M11s [kg/s^2]:', minval(M11s),maxval(M11s)
  write(69,*)' '

  deallocate(non_diag_fact)

end subroutine def_solid_stiffness_terms
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_monopole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

use data_monopole

integer, intent(in) :: ielem,jpol

double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
double precision, intent(in) :: mu(0:npol,0:npol,nelem)
double precision, intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

double precision, intent(in) :: non_diag_fact(0:npol,nel_solid)
double precision, intent(in) :: local_crd_nodes(8,2)

double precision, intent(in) :: alpha_wt_k(0:npol,0:npol)
double precision, intent(in) :: beta_wt_k(0:npol,0:npol)
double precision, intent(in) :: gamma_wt_k(0:npol,0:npol)
double precision, intent(in) :: delta_wt_k(0:npol,0:npol)
double precision, intent(in) :: epsil_wt_k(0:npol,0:npol)
double precision, intent(in) :: zeta_wt_k(0:npol,0:npol)

double precision, intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

double precision, intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

integer          :: ipol
double precision :: dsdxi,dzdeta,dzdxi,dsdeta

! Clumsy to initialize inside a loop... but hey, the easiest way in this setup.
  if ( ielem==1 .and. jpol==0 ) then
     M0_4(0:npol,1:nel_solid)=zero
     M0_w1(0:npol,1:nel_solid)=zero
     M0_w2(0:npol,1:nel_solid)=zero
  endif
  
! ----------------
  do ipol=0,npol
! ----------------
     M11s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))*delta_wt_k(ipol,jpol)+&
          mu(ipol,jpol,ielsolid(ielem))*alpha_wt_k(ipol,jpol)
     M21s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))*zeta_wt_k(ipol,jpol)+&
          mu(ipol,jpol,ielsolid(ielem))*gamma_wt_k(ipol,jpol)
     M41s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))*epsil_wt_k(ipol,jpol)+&
          mu(ipol,jpol,ielsolid(ielem))*beta_wt_k(ipol,jpol)
     M12s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          Ms_z_eta_s_xi_wt_k(ipol,jpol)+&
          mu(ipol,jpol,ielsolid(ielem))*Ms_z_xi_s_eta_wt_k(ipol,jpol)
     M22s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          mu(ipol,jpol,ielsolid(ielem)))*Ms_z_eta_s_eta_wt_k(ipol,jpol)
     M32s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          Ms_z_xi_s_eta_wt_k(ipol,jpol)+&
          mu(ipol,jpol,ielsolid(ielem))*Ms_z_eta_s_xi_wt_k(ipol,jpol)
     M42s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*Ms_z_xi_s_xi_wt_k(ipol,jpol)
     M11z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
          delta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
          two*mu(ipol,jpol,ielsolid(ielem)))*alpha_wt_k(ipol,jpol)
     M21z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
          zeta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))*gamma_wt_k(ipol,jpol)
     M41z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          epsil_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
          two*mu(ipol,jpol,ielsolid(ielem)))*beta_wt_k(ipol,jpol)

     M_2(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          M_z_xi_wt_k(ipol,jpol)
     M_1(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          M_z_eta_wt_k(ipol,jpol)
     M_4(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          M_s_xi_wt_k(ipol,jpol)
     M_3(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          M_s_eta_wt_k(ipol,jpol)

     M_w(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))* &
          massmat_kwts2(ipol,jpol,ielsolid(ielem))
     if (axis_solid(ielem)) M_w(0,jpol,ielem)=zero
! ----------------
  enddo
! ----------------
  if ( axis_solid(ielem) ) then
     call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
          xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))
! AXIS-------------------
     ipol=0

     M0_4(jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))* &
          dsdxi*wt_axial_k(0)*wt(jpol)
     M0_w1(jpol,ielem)=(three*lambda(ipol,jpol,ielsolid(ielem))+ &
          two*mu(ipol,jpol,ielsolid(ielem)))* &
          non_diag_fact(jpol,ielem)
  endif

end subroutine compute_monopole_stiff_terms
!=============================================================================


!-----------------------------------------------------------------------------
subroutine compute_monopole_stiff_terms_ani(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,xi_ani,phi_ani,eta_ani, &
                                       massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k,&
                                       gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

use data_monopole

implicit none
integer, intent(in) :: ielem,jpol

double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
double precision, intent(in) :: mu(0:npol,0:npol,nelem)
double precision, intent(in) :: xi_ani(0:npol,0:npol,nelem)
double precision, intent(in) :: phi_ani(0:npol,0:npol,nelem)
double precision, intent(in) :: eta_ani(0:npol,0:npol,nelem)
double precision, intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

double precision, intent(in) :: non_diag_fact(0:npol,nel_solid)
double precision, intent(in) :: local_crd_nodes(8,2)

double precision, intent(in) :: alpha_wt_k(0:npol,0:npol)
double precision, intent(in) :: beta_wt_k(0:npol,0:npol)
double precision, intent(in) :: gamma_wt_k(0:npol,0:npol)
double precision, intent(in) :: delta_wt_k(0:npol,0:npol)
double precision, intent(in) :: epsil_wt_k(0:npol,0:npol)
double precision, intent(in) :: zeta_wt_k(0:npol,0:npol)

double precision, intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

double precision, intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

integer          :: ipol
double precision :: dsdxi,dzdeta,dzdxi,dsdeta
double precision :: theta
double precision :: C11, C22, C33, C12, C13, C23, C15, C25, C35, C44, C46, C55, C66, Ctmp
double precision :: lambdal, mul, xil, phil, etal

! Clumsy to initialize inside a loop... but hey, the easiest way in this setup.
  if ( ielem==1 .and. jpol==0 ) then
     M0_4(0:npol,1:nel_solid)=zero
     M0_w1(0:npol,1:nel_solid)=zero
     M0_w2(0:npol,1:nel_solid)=zero
  endif
  
! ----------------
  do ipol=0,npol
! ----------------
     theta = thetacoord(ipol, jpol, ielsolid(ielem))
   
     lambdal = lambda(ipol,jpol,ielsolid(ielem))
     mul = mu(ipol,jpol,ielsolid(ielem))
     xil = xi_ani(ipol,jpol,ielsolid(ielem))
     phil = phi_ani(ipol,jpol,ielsolid(ielem))
     etal = eta_ani(ipol,jpol,ielsolid(ielem))
     
     C11 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 1, 1)
     C12 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 2, 2)
     C13 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 3, 3)
     C15 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 3, 1)
     C22 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 2, 2, 2)
     C23 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 2, 3, 3)
     C25 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 2, 3, 1)
     C33 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 3, 3, 3)
     C35 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 3, 3, 1)
     C44 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 3, 2, 3)
     C46 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 3, 1, 2)
     C55 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 1, 3, 1)
     C66 = c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 2, 1, 2)

    ! Test for the components that should be zero:
     if (do_mesh_tests) then
        if ( ielem==1 .and. jpol==0 .and. ipol==0 ) then
           if (lpr) write(6,*) ' Test for the components of c_ijkl that should be zero in anisotropic case'
        endif
        Ctmp = zero
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 1, 1, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 2, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 2, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 3, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 3, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 2, 3, 3, 1))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, theta, 3, 1, 1, 2))
        
        if (Ctmp > smallval_sngl) then
           write(6,*)procstrg,' ERROR: some stiffness term that should be zero '
           write(6,*)procstrg,'        is not: in compute_monopole_stiff_terms_ani()'
           stop
        endif
     endif

     !Uncomment to test for comparison to isotropic values
     !
     !if (dabs(C11 - lambdal - 2 * mul) / C11 > smallval_sngl) then
     !   write(6,*) 'C11 is wrong'
     !   stop
     !endif
     !if (dabs(C22 - lambdal - 2 * mul) / C22 > smallval_sngl) then
     !   write(6,*) 'C22 is wrong'
     !   stop
     !endif
     !if (dabs(C33 - lambdal - 2 * mul) / C33  > smallval_sngl) then
     !   write(6,*) 'C33 is wrong'
     !   stop
     !endif
     !if (dabs(C44 - mul) / C44 > smallval_sngl) then
     !   write(6,*) 'C44 is wrong'
     !   stop
     !endif
     !if (dabs(C55 - mul) / C55 > smallval_sngl) then
     !   write(6,*) 'C55 is wrong'
     !   stop
     !endif
     !if (dabs(C66 - mul) / C66 > smallval_sngl) then
     !   write(6,*) 'C66 is wrong'
     !   stop
     !endif
     !if (dabs(C12 - lambdal) / C12 > smallval_sngl) then
     !   write(6,*) 'C12 is wrong'
     !   stop
     !endif
     !if (dabs(C13 - lambdal) / C13 > smallval_sngl) then
     !   write(6,*) 'C13 is wrong'
     !   stop
     !endif
     !if (dabs(C23 - lambdal) / C23 > smallval_sngl) then
     !   write(6,*) 'C23 is wrong'
     !   stop
     !endif
     !if (dabs(C15) > smallval_sngl) then
     !   write(6,*) 'C15 is wrong'
     !   stop
     !endif
     !if (dabs(C25) > smallval_sngl) then
     !   write(6,*) 'C25 is wrong'
     !   stop
     !endif
     !if (dabs(C35) > smallval_sngl) then
     !   write(6,*) 'C35 is wrong'
     !   stop
     !endif
     !if (dabs(C46) > smallval_sngl) then
     !   write(6,*) 'C46 is wrong'
     !   stop
     !endif

     !hard coded C for pure vp anisotropy, can be removed after some more
     !testing (MvD)
     !C11 = C11 + (phil - one) * (lambdal + two * mul) * (dsin(theta)**4)
     !C33 = C33 + (phil - one) * (lambdal + two * mul) * (dcos(theta)**4)
     !C13 = C13 + (phil - one) * (lambdal + two * mul) * (dsin(theta)**2) * (dcos(theta)**2)
     !C15 = C15 + (phil - one) * (lambdal + two * mul) * (dsin(theta)**3) * (dcos(theta)**1)
     !C35 = C35 + (phil - one) * (lambdal + two * mul) * (dsin(theta)**1) * (dcos(theta)**3)
     !C55 = C55 + (phil - one) * (lambdal + two * mul) * (dsin(theta)**2) * (dcos(theta)**2)
    
     M11s(ipol,jpol,ielem) = C11 * delta_wt_k(ipol,jpol) &
                           + C15 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C15 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * alpha_wt_k(ipol,jpol)

     M21s(ipol,jpol,ielem) = C11 * zeta_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C55 * gamma_wt_k(ipol,jpol)

     M41s(ipol,jpol,ielem) = C11 * epsil_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C55 * beta_wt_k(ipol,jpol)

     M12s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C55 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M22s(ipol,jpol,ielem) = C15 * zeta_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C35 * gamma_wt_k(ipol,jpol)

     M32s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M42s(ipol,jpol,ielem) = C15 * epsil_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C35 * beta_wt_k(ipol,jpol)

     M11z(ipol,jpol,ielem) = C55 * delta_wt_k(ipol,jpol) &
                           + C35 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C33 * alpha_wt_k(ipol,jpol)

     M21z(ipol,jpol,ielem) = C55 * zeta_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C33 * gamma_wt_k(ipol,jpol)

     M41z(ipol,jpol,ielem) = C55 * epsil_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C33 * beta_wt_k(ipol,jpol)


     M_1(ipol,jpol,ielem) = C12 * M_z_eta_wt_k(ipol,jpol) + C25  * M_s_eta_wt_k(ipol,jpol)
     M_2(ipol,jpol,ielem) = C12 * M_z_xi_wt_k(ipol,jpol)  + C25  * M_s_xi_wt_k(ipol,jpol)
     M_3(ipol,jpol,ielem) = C23 * M_s_eta_wt_k(ipol,jpol) + C25  * M_z_eta_wt_k(ipol,jpol)
     M_4(ipol,jpol,ielem) = C23 * M_s_xi_wt_k(ipol,jpol)  + C25  * M_z_xi_wt_k(ipol,jpol)

     M_w(ipol,jpol,ielem) = C22 * massmat_kwts2(ipol,jpol,ielsolid(ielem))

     if (axis_solid(ielem)) M_w(0,jpol,ielem)=zero
     
     if (axis_solid(ielem) .and. ipol==0) then
        call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
             xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))
   
        M0_4(jpol,ielem) = C23 * dsdxi * wt_axial_k(0) * wt(jpol) &
                         + C25 * dzdxi * wt_axial_k(0) * wt(jpol)
        M0_w1(jpol,ielem) = (2 * C12 + C22) * non_diag_fact(jpol,ielem)
        M0_w2(jpol,ielem) = C25 * non_diag_fact(jpol,ielem)
     endif
  enddo
! ----------------

end subroutine compute_monopole_stiff_terms_ani
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_dipole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                      lambda,mu,massmat_kwts2, &
                                      non_diag_fact,alpha_wt_k,beta_wt_k,&
                                      gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                      zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                      M_z_eta_wt_k,M_s_eta_wt_k, &
                                      Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                      Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

use data_dipole

integer, intent(in) :: ielem,jpol

double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
double precision, intent(in) :: mu(0:npol,0:npol,nelem)
double precision, intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

double precision, intent(in) :: non_diag_fact(0:npol,nel_solid)
double precision, intent(in) :: local_crd_nodes(8,2)

double precision, intent(in) :: alpha_wt_k(0:npol,0:npol)
double precision, intent(in) :: beta_wt_k(0:npol,0:npol)
double precision, intent(in) :: gamma_wt_k(0:npol,0:npol)
double precision, intent(in) :: delta_wt_k(0:npol,0:npol)
double precision, intent(in) :: epsil_wt_k(0:npol,0:npol)
double precision, intent(in) :: zeta_wt_k(0:npol,0:npol)

double precision, intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

double precision, intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
double precision, intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

integer          :: ipol
double precision :: dsdxi,dzdeta,dzdxi,dsdeta

! Clumsy to initialize inside a loop... but hey, the easiest way in this setup.
  if ( ielem==1 .and. jpol==0 ) then
     M0_w_mu(0:npol,1:nel_solid)=zero
     M0_s_xi_mu(0:npol,1:nel_solid)=zero
     M0_z_eta_2_lam_mu(0:npol,1:nel_solid)=zero
     M0_z_eta_4lm_w_l3m(0:npol,1:nel_solid)=zero
  endif
  
! ----------------
  do ipol=0,npol
! ----------------
! + and - components
     M11s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          three*mu(ipol,jpol,ielsolid(ielem)))*delta_wt_k(ipol,jpol)+&
          two*mu(ipol,jpol,ielsolid(ielem))*alpha_wt_k(ipol,jpol)
     M21s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
          three*mu(ipol,jpol,ielsolid(ielem)))*zeta_wt_k(ipol,jpol)+&
          two*mu(ipol,jpol,ielsolid(ielem))*gamma_wt_k(ipol,jpol)
     M41s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          three*mu(ipol,jpol,ielsolid(ielem)))*epsil_wt_k(ipol,jpol)+&
          two*mu(ipol,jpol,ielsolid(ielem))*beta_wt_k(ipol,jpol)

     M12s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*delta_wt_k(ipol,jpol)
     M22s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*zeta_wt_k(ipol,jpol)
     M42s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*epsil_wt_k(ipol,jpol)

     M13s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
          Ms_z_eta_s_xi_wt_k(ipol,jpol)+mu(ipol,jpol,ielsolid(ielem))*&
          Ms_z_xi_s_eta_wt_k(ipol,jpol)
     M32s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*Ms_z_eta_s_eta_wt_k(ipol,jpol)
     M33s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
          Ms_z_xi_s_eta_wt_k(ipol,jpol)+mu(ipol,jpol,ielsolid(ielem))*&
          Ms_z_eta_s_xi_wt_k(ipol,jpol)
     M43s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*Ms_z_xi_s_xi_wt_k(ipol,jpol)

! z component
     M11z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          delta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
          two*mu(ipol,jpol,ielsolid(ielem)))*alpha_wt_k(ipol,jpol)
     M21z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          zeta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
          two*mu(ipol,jpol,ielsolid(ielem)))*gamma_wt_k(ipol,jpol)
     M41z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          epsil_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
          two*mu(ipol,jpol,ielsolid(ielem)))*beta_wt_k(ipol,jpol)

! D^x_y terms
     Mz_xi_p(ipol,jpol,ielem)=two*(lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*M_z_xi_wt_k(ipol,jpol)
     Mz_xi_m(ipol,jpol,ielem)=two*(lambda(ipol,jpol,ielsolid(ielem))-&
          mu(ipol,jpol,ielsolid(ielem)))*M_z_xi_wt_k(ipol,jpol)  
     Mz_eta_p(ipol,jpol,ielem)=two* &
          (lambda(ipol,jpol,ielsolid(ielem))+&
          mu(ipol,jpol,ielsolid(ielem)))*M_z_eta_wt_k(ipol,jpol)
     Mz_eta_m(ipol,jpol,ielem)=two* &
          (lambda(ipol,jpol,ielsolid(ielem))-&
          mu(ipol,jpol,ielsolid(ielem)))*M_z_eta_wt_k(ipol,jpol)
     Ms_xi_mu(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
          M_s_xi_wt_k(ipol,jpol)
     Ms_xi_2lam(ipol,jpol,ielem) = &
          two*lambda(ipol,jpol,ielsolid(ielem))*M_s_xi_wt_k(ipol,jpol)
     Ms_eta_mu(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          M_s_eta_wt_k(ipol,jpol)
     Ms_eta_2lam(ipol,jpol,ielem) = &
          two*lambda(ipol,jpol,ielsolid(ielem))*M_s_eta_wt_k(ipol,jpol)

! AXIS-------------------
! 2nd order terms
     M_w_4_lam_mu(ipol,jpol,ielem)=four* &
          (lambda(ipol,jpol,ielsolid(ielem)) + &
          three*mu(ipol,jpol,ielsolid(ielem)))* &
          massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w_mu(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          massmat_kwts2(ipol,jpol,ielsolid(ielem))
! ----------------
  enddo
! ----------------

! AXIS-------------------
  if ( axis_solid(ielem) ) then

     M11s(0,jpol,ielem)=zero; M12s(0,jpol,ielem)=zero
     M42s(0,jpol,ielem)=zero; M32s(0,jpol,ielem)=zero
     M43s(0,jpol,ielem)=zero; M11z(0,jpol,ielem)=zero

     Mz_xi_p(0,jpol,ielem)=zero; Mz_xi_m(0,jpol,ielem)=zero
     Mz_eta_p(0,jpol,ielem)=zero; Mz_eta_m(0,jpol,ielem)=zero
     Ms_xi_mu(0,jpol,ielem)=zero; Ms_xi_2lam(0,jpol,ielem)=zero
     Ms_eta_mu(0,jpol,ielem)=zero; Ms_eta_2lam(0,jpol,ielem)=zero
     M_w_4_lam_mu(0,jpol,ielem)=zero; M_w_mu(0,jpol,ielem)=zero

     call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
          xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))

     dzdxi=zero; dsdeta=zero; ipol=0
     M0_w_mu(jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          non_diag_fact(jpol,ielem)
     M0_s_xi_mu(jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
          dsdxi*wt_axial_k(0)*wt(jpol)
     M0_z_eta_2_lam_mu(jpol,ielem)=two*&
          (lambda(ipol,jpol,ielsolid(ielem))+ &
          mu(ipol,jpol,ielsolid(ielem)))*dzdeta*wt_axial_k(0)*wt(jpol)
     M0_z_eta_4lm_w_l3m(jpol,ielem)=  &
          (four*(lambda(ipol,jpol,ielsolid(ielem))- &
          mu(ipol,jpol,ielsolid(ielem)))*dzdeta*&
          wt_axial_k(0)*wt(jpol))+&
          (four*(lambda(ipol,jpol,ielsolid(ielem))+ &
          four*mu(ipol,jpol,ielsolid(ielem)))*&
          non_diag_fact(jpol,ielem))
  endif

end subroutine compute_dipole_stiff_terms
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_quadrupole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                      lambda,mu,massmat_kwts2, &
                                      non_diag_fact,alpha_wt_k,beta_wt_k,&
                                      gamma_wt_k,delta_wt_k,epsil_wt_k,&
                                      zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k,&
                                      M_z_eta_wt_k,M_s_eta_wt_k, &
                                      Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k,&
                                      Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

use data_quadrupole

  integer, intent(in) :: ielem,jpol

  double precision, intent(in) :: lambda(0:npol,0:npol,nelem)
  double precision, intent(in) :: mu(0:npol,0:npol,nelem)
  double precision, intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

  double precision, intent(in) :: non_diag_fact(0:npol,nel_solid)
  double precision, intent(in) :: local_crd_nodes(8,2)

  double precision, intent(in) :: alpha_wt_k(0:npol,0:npol)
  double precision, intent(in) :: beta_wt_k(0:npol,0:npol)
  double precision, intent(in) :: gamma_wt_k(0:npol,0:npol)
  double precision, intent(in) :: delta_wt_k(0:npol,0:npol)
  double precision, intent(in) :: epsil_wt_k(0:npol,0:npol)
  double precision, intent(in) :: zeta_wt_k(0:npol,0:npol)

  double precision, intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
  double precision, intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
  double precision, intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
  double precision, intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

  double precision, intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
  double precision, intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
  double precision, intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
  double precision, intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

  integer          :: ipol
  double precision :: dsdxi,dzdeta,dzdxi,dsdeta

!       ----------------
           do ipol=0,npol
!       ----------------
              M11s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
                  two*mu(ipol,jpol,ielsolid(ielem)))*delta_wt_k(ipol,jpol)+&
                  mu(ipol,jpol,ielsolid(ielem))*alpha_wt_k(ipol,jpol)
              M21s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
                  two*mu(ipol,jpol,ielsolid(ielem)))*zeta_wt_k(ipol,jpol)+&
                  mu(ipol,jpol,ielsolid(ielem))*gamma_wt_k(ipol,jpol)
              M41s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
                  two*mu(ipol,jpol,ielsolid(ielem)))*epsil_wt_k(ipol,jpol)+&
                  mu(ipol,jpol,ielsolid(ielem))*beta_wt_k(ipol,jpol)
              M12s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                  Ms_z_eta_s_xi_wt_k(ipol,jpol)+mu(ipol,jpol,ielsolid(ielem))*&
                  Ms_z_xi_s_eta_wt_k(ipol,jpol)
              M22s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
                  mu(ipol,jpol,ielsolid(ielem)))*Ms_z_eta_s_eta_wt_k(ipol,jpol)
              M32s(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                  Ms_z_xi_s_eta_wt_k(ipol,jpol)+mu(ipol,jpol,ielsolid(ielem))*&
                  Ms_z_eta_s_xi_wt_k(ipol,jpol)
              M42s(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+ &
                  mu(ipol,jpol,ielsolid(ielem)))*Ms_z_xi_s_xi_wt_k(ipol,jpol)
              M11z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
                  delta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
                  two*mu(ipol,jpol,ielsolid(ielem)))*alpha_wt_k(ipol,jpol)
              M21z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
                  zeta_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))+&
                  two*mu(ipol,jpol,ielsolid(ielem)))*gamma_wt_k(ipol,jpol)
              M41z(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))* &
                  epsil_wt_k(ipol,jpol)+(lambda(ipol,jpol,ielsolid(ielem))&
                  +two*mu(ipol,jpol,ielsolid(ielem)))*beta_wt_k(ipol,jpol)

              M1phi(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
                  (delta_wt_k(ipol,jpol)+alpha_wt_k(ipol,jpol))
              M2phi(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
                   (zeta_wt_k(ipol,jpol)+ gamma_wt_k(ipol,jpol))
              M4phi(ipol,jpol,ielem)=mu(ipol,jpol,ielsolid(ielem))*&
                   (epsil_wt_k(ipol,jpol)+beta_wt_k(ipol,jpol))
              Mz_xi_lam(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                   M_z_xi_wt_k(ipol,jpol)
              Mz_eta_lam(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                   M_z_eta_wt_k(ipol,jpol)
              Ms_xi_lam(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                   M_s_xi_wt_k(ipol,jpol)
              Ms_eta_lam(ipol,jpol,ielem)=lambda(ipol,jpol,ielsolid(ielem))*&
                   M_s_eta_wt_k(ipol,jpol)

              Mz_xi_2mu(ipol,jpol,ielem)=two*mu(ipol,jpol,ielsolid(ielem))*&
                   M_z_xi_wt_k(ipol,jpol)
              Mz_eta_2mu(ipol,jpol,ielem)=two*mu(ipol,jpol,ielsolid(ielem))*&
                   M_z_eta_wt_k(ipol,jpol)
              Ms_xi_2mu(ipol,jpol,ielem)=two*mu(ipol,jpol,ielsolid(ielem))*&
                   M_s_xi_wt_k(ipol,jpol)
              Ms_eta_2mu(ipol,jpol,ielem)=two*mu(ipol,jpol,ielsolid(ielem))*&
                   M_s_eta_wt_k(ipol,jpol)

              M_w_lam_6mu(ipol,jpol,ielem)=(lambda(ipol,jpol,ielsolid(ielem))+&
                   6.d0*mu(ipol,jpol,ielsolid(ielem)))* &
                   massmat_kwts2(ipol,jpol,ielsolid(ielem))
              M_w_min_2lam_6mu(ipol,jpol,ielem)=&
                   -(two*lambda(ipol,jpol,ielsolid(ielem))+&
                   6.d0*mu(ipol,jpol,ielsolid(ielem)))* &
                   massmat_kwts2(ipol,jpol,ielsolid(ielem))
              M_w_4lam_9mu(ipol,jpol,ielem)=&
                   (four*lambda(ipol,jpol,ielsolid(ielem))+&
                   9.d0*mu(ipol,jpol,ielsolid(ielem)))* &
                   massmat_kwts2(ipol,jpol,ielsolid(ielem))
              M_w_4mu(ipol,jpol,ielem)=four*mu(ipol,jpol,ielsolid(ielem))*&
                   massmat_kwts2(ipol,jpol,ielsolid(ielem))
!       ----------------
           enddo ! ipol
!       ----------------

   ! AXIS-------------------
           if ( axis_solid(ielem) ) then
              ipol=0
              M11s(0,jpol,ielem)=zero; M22s(0,jpol,ielem)=zero
              M42s(0,jpol,ielem)=zero; M11z(0,jpol,ielem)=zero
              M1phi(0,jpol,ielem)=zero
              Mz_xi_lam(0,jpol,ielem)=zero; Mz_eta_lam(0,jpol,ielem)=zero
              Ms_xi_lam(0,jpol,ielem)=zero; Ms_eta_lam(0,jpol,ielem)=zero
              Mz_xi_2mu(0,jpol,ielem)=zero; Mz_eta_2mu(0,jpol,ielem)=zero
              Ms_xi_2mu(0,jpol,ielem)=zero; Ms_eta_2mu(0,jpol,ielem)=zero
              M_w_lam_6mu(0,jpol,ielem)=zero
              M_w_min_2lam_6mu(0,jpol,ielem)=zero
              M_w_4lam_9mu(0,jpol,ielem)=zero; M_w_4mu(0,jpol,ielem)=zero

              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta,&
                   xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))
              dzdxi=zero; dsdeta=zero

              M0_w_4mu(jpol,ielem)=four*mu(ipol,jpol,ielsolid(ielem))*&
                    non_diag_fact(jpol,ielem)

              M0_z_eta_3m_w_4l_9m(jpol,ielem)=&
                   (-two*mu(ipol,jpol,ielsolid(ielem))*dzdeta* &
                   wt_axial_k(0)*wt(jpol))+ &
                   ((four*lambda(ipol,jpol,ielsolid(ielem))+&
                   9.d0*mu(ipol,jpol,ielsolid(ielem)))*&
                   non_diag_fact(jpol,ielem))

              M0_z_eta_2l_w_l_6m(jpol,ielem)=&
                   (two*lambda(ipol,jpol,ielsolid(ielem))*dzdeta* &
                   wt_axial_k(0)*wt(jpol))+ &
                   ((lambda(ipol,jpol,ielsolid(ielem))+&
                   6.d0*mu(ipol,jpol,ielsolid(ielem)))*&
                   non_diag_fact(jpol,ielem))

              M0_z_eta_2lm_min_w_2l_6m(jpol,ielem)=&
                   (two*(mu(ipol,jpol,ielsolid(ielem))-&
                   lambda(ipol,jpol,ielsolid(ielem)))*dzdeta* &
                   wt_axial_k(0)*wt(jpol)) - &
                   ((two*lambda(ipol,jpol,ielsolid(ielem))+&
                   6.d0*mu(ipol,jpol,ielsolid(ielem)))*&
                   non_diag_fact(jpol,ielem))
           endif ! axial element

end subroutine compute_quadrupole_stiff_terms
!=============================================================================

!-----------------------------------------------------------------------------
!double precision function c_ijkl_ani(lambda, mu, epsilon_ani, gamma_ani, delta_ani, theta, i, j, k, l)
double precision function c_ijkl_ani(lambda, mu, xi_ani, phi_ani, eta_ani, theta, i, j, k, l)
!
! returns the stiffness tensor as defined in Nolet(2008), Eq. (16.2)
! i, j, k and l should be in [1,3]
! MvD
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

implicit none

double precision, intent(in) :: lambda, mu, xi_ani, phi_ani, eta_ani
double precision, intent(in) :: theta 
integer, intent(in) :: i, j, k, l
double precision, dimension(1:3, 1:3) :: deltaf
double precision, dimension(1:3) :: s

deltaf = zero
deltaf(1,1) = one
deltaf(2,2) = one
deltaf(3,3) = one

s(1) = dsin(theta)
s(2) = zero
s(3) = dcos(theta)

c_ijkl_ani = zero

! isotropic part:
c_ijkl_ani = c_ijkl_ani + lambda * deltaf(i,j) * deltaf(k,l)

c_ijkl_ani = c_ijkl_ani + mu * (deltaf(i,k) * deltaf(j,l) + deltaf(i,l) * deltaf(j,k))


! anisotropic part:
! in xi, phi, eta

c_ijkl_ani = c_ijkl_ani &
    + ((eta_ani - one) * lambda + two * eta_ani * mu * (one - one / xi_ani)) &
        * (deltaf(i,j) * s(k) * s(l) + deltaf(k,l) * s(i) * s(j))
    
c_ijkl_ani = c_ijkl_ani &
    + mu * (one / xi_ani - one) &
        * (deltaf(i,k) * s(j) * s(l) + deltaf(i,l) * s(j) * s(k) + &
           deltaf(j,k) * s(i) * s(l) + deltaf(j,l) * s(i) * s(k))

c_ijkl_ani = c_ijkl_ani &
    + ((one - two * eta_ani + phi_ani) * (lambda + two * mu) + (4. * eta_ani - 4.) * mu / xi_ani) &
        * (s(i) * s(j) * s(k) * s(l))

end function c_ijkl_ani
!=============================================================================

!-----------------------------------------------------------------------------
subroutine def_fluid_stiffness_terms(rho,massmat_kwts2)
!
! Fluid precomputed matrices definitions for all sources.
! Note that in this routine terms alpha etc. are scalars 
! (as opposed to the solid case of being elemental arrays).
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
include "mesh_params.h"
double precision, intent(in)  :: rho(0:npol,0:npol,nelem)
double precision, intent(in)  :: massmat_kwts2(0:npol,0:npol,nelem)

double precision, allocatable :: non_diag_fact(:,:)
double precision              :: local_crd_nodes(8,2)
double precision              :: alpha_wt_k,beta_wt_k,gamma_wt_k
double precision              :: delta_wt_k,epsil_wt_k,zeta_wt_k
integer                       :: iel,ipol,jpol,inode

  allocate(non_diag_fact(0:npol,1:nel_fluid))
  
  do iel=1,nel_fluid

     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1),&
             local_crd_nodes(inode,2),ielfluid(iel),inode)
     end do

     do jpol=0,npol
        do ipol=0, npol
! ::::::::::::::::non-axial elements::::::::::::::::
           if ( .not. axis_fluid(iel) ) then

              alpha_wt_k = alpha(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              beta_wt_k = beta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              gamma_wt_k = gamma1(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              delta_wt_k = delta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              epsil_wt_k = epsilon1(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              zeta_wt_k = zeta(eta(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)

 ! ::::::::::::::::axial elements::::::::::::::::
           elseif (axis_fluid(iel) ) then
              alpha_wt_k = alphak(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              beta_wt_k = betak(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel))&
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              gamma_wt_k = gammak(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              delta_wt_k = deltak(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              epsil_wt_k = epsilonk(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel))&
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              zeta_wt_k = zetak(xi_k(ipol),&
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
           endif

           M1chi_fl(ipol,jpol,iel)=(delta_wt_k + alpha_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
           M2chi_fl(ipol,jpol,iel)=(zeta_wt_k + gamma_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
           M4chi_fl(ipol,jpol,iel)=(epsil_wt_k + beta_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
        enddo
     enddo
  enddo

!  2nd order term for multipoles
  if (src_type(1)=='dipole' .or. src_type(1)=='quadpole') then

     non_diag_fact(:,:) = zero
     do iel=1,nel_fluid

! ::::::::::::::::all elements::::::::::::::::
        M_w_fl(:,:,iel)=massmat_kwts2(:,:,ielfluid(iel))/rho(:,:,ielfluid(iel))

! ::::::::::::::::axial elements::::::::::::::::
        if ( axis_fluid(iel) ) then
           do inode = 1, 8
              call compute_coordinates_mesh(local_crd_nodes(inode,1),&
                   local_crd_nodes(inode,2),ielfluid(iel),inode)
           end do
           do jpol = 0,npol 
              non_diag_fact(jpol,iel) = wt_axial_k(0)*wt(jpol)* &
                   jacobian(xi_k(0),eta(jpol),local_crd_nodes,ielfluid(iel))/&
                   s_over_oneplusxi_axis(xi_k(0), &
                   eta(jpol),local_crd_nodes,ielfluid(iel))
           end do

! axial masking of main term
           M_w_fl(0,0:npol,iel)=zero

! additional axial term
           M0_w_fl(:,iel)=non_diag_fact(:,iel)/rho(0,:,ielfluid(iel))
        endif ! axis
     enddo
     if (src_type(1)=='quadpole') then 
        M_w_fl=four*M_w_fl
        M0_w_fl=four*M0_w_fl
     endif
  endif ! multipole

  write(69,*)' '
  write(69,*)'Min/max M1chi_fl [m^4/kg]:',minval(M1chi_fl),maxval(M1chi_fl)

  deallocate(non_diag_fact)

end subroutine def_fluid_stiffness_terms
!=============================================================================

!-----------------------------------------------------------------------------
subroutine def_solid_fluid_boundary_terms
!
! Defines the 1-d vector-array bdry_matr which acts as the diagonal matrix
! to accomodate the exchange of fields across solid-fluid boundaries 
! in both directions. Take note of the sign conventions in accordance with 
! those used in the time loop.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use commun, only : psum
use data_io

include 'mesh_params.h'

!double precision, intent(in) :: rho(0:npol,0:npol,nelem)
double precision             :: local_crd_nodes(8,2)
double precision             :: s,z,r,theta,rf,thetaf
double precision             :: theta1,theta2,r1,r2,delta_th,bdry_sum
integer                      :: iel,ielglob,ipol,inode,idom
integer                      :: count_lower_disc,count_upper_disc

  bdry_sum = zero

! check if proc has boundary elements
  if (have_bdry_elem) then

  count_lower_disc=0; count_upper_disc=0

  do iel=1,nel_bdry

! Map from boundary to global indexing. Choice of the solid side is random...
     ielglob=ielsolid(bdry_solid_el(iel))

! closest to axis
     call compute_coordinates(s,z,r1,theta1,ielglob,0,bdry_jpol_solid(iel))

     call compute_coordinates(s,z,rf,thetaf,ielfluid(bdry_fluid_el(iel)),&
                              0,bdry_jpol_fluid(iel))

! test if the mapping of solid element & jpol numbers agrees for solid & fluid
! af rock-lgit
!    if (rf/=r1 .or. thetaf/=theta1) then 
     if (abs( (rf-r1) /r1 ) > 1.e-5 .or. abs((thetaf-theta1)) > 1.e-3) then 
        write(6,*)
        write(6,*)procstrg,'Problem with boundary term mapping near axis!'
        write(6,*)procstrg,'radius,theta solid index:',r1/1.d3,theta1/pi*180.
        write(6,*)procstrg,'radius,theta fluid index:',rf/1.d3,thetaf/pi*180.
        stop
     endif

! furthest away from axis
     call compute_coordinates(s,z,r2,theta2,ielglob,npol,bdry_jpol_solid(iel))
     call compute_coordinates(s,z,rf,thetaf,ielfluid(bdry_fluid_el(iel)),&
                             npol,bdry_jpol_fluid(iel))

! test if the mapping of solid element & jpol numbers agrees for solid & fluid
!af rock-lgit
!    if (rf/=r2 .or. thetaf/=theta2) then 
     if (abs( (rf-r2) /r2 ) > 1.e-5 .or. abs((thetaf-theta2)) > 1.e-3) then 
        write(6,*)
        write(6,*)procstrg,'Problem with boundary term mapping far axis!'
        write(6,*)procstrg,'radius,theta solid index:',r2/1.d3,theta2/pi*180.
        write(6,*)procstrg,'radius,theta fluid index:',rf/1.d3,thetaf/pi*180.
        stop
     endif

     if ( abs(r1-r2)>min_distance_dim) then 
        write(6,*)
        write(6,*)procstrg,'Problem with S/F boundary element',ielglob
        write(6,*)procstrg,'radii at min./max theta are not equal!'
        write(6,*)procstrg,'r1,r2 [km],theta [deg]:', &
                            r1/1000.,r2/1000.,theta1*180./pi
        stop
     endif
!1/2 comes from d_th=1/2 (th_2 - th_1)d_xi; abs() takes care of southern elems
     delta_th=half*abs(theta2-theta1)

! ::::::::::::::::axial elements::::::::::::::::
     if ( axis(ielglob) ) then

        if (abs(sin(theta1))*two*pi*r1>min_distance_dim) then
           write(6,*)
           write(6,*)procstrg,'Problem with axial S/F boundary element',ielglob
           write(6,*)procstrg,'Min theta is not exactly on the axis'
           write(6,*)procstrg,'r [km],theta [deg]:',r1/1000.,theta1*180./pi
           stop
        endif

        do inode = 1, 8
           call compute_coordinates_mesh(local_crd_nodes(inode,1),&
                local_crd_nodes(inode,2),ielglob,inode)
        end do

        ! I>0 (off the axis)
        do ipol=1,npol
           call compute_coordinates(s,z,r,theta,ielglob,ipol,&
                                    bdry_jpol_solid(iel))
           if(abs(r-r1)>min_distance_dim) then 
              write(6,*)
              write(6,*)procstrg,'Problem with axial S/F boundary element',&
                        ielglob
              write(6,*)procstrg,'radius at ipol=',ipol,'different from ipol=0'
              write(6,*)procstrg,'r,r1 [km],theta [deg]:',r/1000.,r1/1000., &
                                                           theta*180./pi 
              stop
           endif
           bdry_matr(ipol,iel,1)=delta_th*wt_axial_k(ipol)* &
                dsin(theta)/(one+xi_k(ipol))*dsin(theta)
           bdry_matr(ipol,iel,2)=delta_th*wt_axial_k(ipol)* &
                dsin(theta)/(one+xi_k(ipol))*dcos(theta)

! Test: Integrated over the whole boundary, the boundary term san sin/cos
! equals two: \int_0^pi sin(\theta) d\theta = two, i.e. 4 if 2 boundaries
        bdry_sum = bdry_sum + delta_th*wt_axial_k(ipol)*dsin(theta)/ &
                              (one+xi_k(ipol))
        enddo

! I=0 axis 
        bdry_sum = bdry_sum + 1/r * delta_th * wt_axial_k(0) * &
                           s_over_oneplusxi_axis(xi_k(0), &
                           eta(bdry_jpol_solid(iel)),local_crd_nodes,ielglob)

! I=0 axis itself
        bdry_matr(0,iel,1)=zero ! note sin(0) = 0    

! need factor 1/r to compensate for using s/(1+xi) rather than sin theta/(1+xi)
! note cos(0) = 1
        bdry_matr(0,iel,2)= 1/r * delta_th * wt_axial_k(0) * &
                           s_over_oneplusxi_axis(xi_k(0), &
                           eta(bdry_jpol_solid(iel)),local_crd_nodes,ielglob)
        write(69,*)
        write(69,11)'S/F axial r[km], theta[deg]   :',r1/1000.,theta1*180/pi
        write(69,10)'S/F axial element (bdry,glob.):',iel,ielglob
        write(69,9)'S/F deltath, sigma_0, s/(1+xi):',delta_th*180/pi, &
                                                     wt_axial_k(0),delta_th*r
10 format(a33,2(i8))
11 format(a33,2(1pe14.5))
9 format(a33,3(1pe14.5))

! testing some algebra...
        if (abs(s_over_oneplusxi_axis(xi_k(0),eta(bdry_jpol_solid(iel)), &
            local_crd_nodes,ielglob)-delta_th*r)> min_distance_dim) then 
           write(6,*)
           write(6,*)procstrg,&
                     'Problem with some axialgebra/definitions, elem:',ielglob
           write(6,*)procstrg,'r [km],theta [deg]:',r1/1000.,theta1*180/pi
           write(6,*)procstrg,'s_0 / (1+xi_0)  =',&
                                         s_over_oneplusxi_axis(xi_k(0),&
                                         eta(bdry_jpol_solid(iel)), &
                                         local_crd_nodes,ielglob)
           write(6,*)procstrg,'1/2 theta2 r_sf =',delta_th*r
           write(6,*)procstrg,'...are not the same :('
           write(6,*)procstrg,'xi,eta,eltype', &
                      xi_k(0),eta(bdry_jpol_solid(iel)),eltype(ielglob)
           write(6,*)procstrg,'theta1,theta2',theta1*180/pi,theta2*180/pi
           write(6,*)procstrg,'s,z min:', &
                     minval(local_crd_nodes(:,1)),minval(local_crd_nodes(:,2))
           stop
        endif

! ::::::::::::::::non-axial elements::::::::::::::::
        else
           do ipol=0,npol

           call compute_coordinates(s,z,r,theta,ielglob,ipol,&
                                    bdry_jpol_solid(iel))

!           write(6,12)'r,r1 [km],theta [deg]:',r/1000.,r1/1000., &
!                                              theta*180./pi 

           if (abs(r-r1)>min_distance_dim) then 
              write(6,*)
              write(6,*)procstrg,&
                        'Problem with non-axial S/F boundary element',ielglob
              write(6,*)procstrg,'radius at ipol=',ipol,'different from ipol=0'
              write(6,12)procstrg,'r,r1 [km],theta [deg]:',r/1000.,r1/1000., &
                                                           theta*180./pi 
              stop
           endif

           bdry_matr(ipol,iel,1)=delta_th*wt(ipol)*dsin(theta)*dsin(theta)
           bdry_matr(ipol,iel,2)=delta_th*wt(ipol)*dsin(theta)*dcos(theta)

! Test: Integrated over the whole boundary, the boundary term san sin/cos
! equals two: \int_0^pi sin(\theta) d\theta = two, i.e. 4 if 2 boundaries
        bdry_sum = bdry_sum + delta_th*wt(ipol)*dsin(theta)

     enddo
  endif ! ax/nonax

! Define the term such that B >(<) 0 of solid above(below) fluid
        do idom=1,ndisc
           if ( .not. solid_domain(idom) ) then
              ! run a check to make sure radius is either discontinuity
              if ( abs(r1-discont(idom))>min_distance_dim .and. &
                   abs(r1-discont(idom+1))>min_distance_dim ) then 
                 write(6,*)
                 write(6,*)procstrg, &
                           'Problem: S/F boundary radius is not one of the'
                 write(6,*)procstrg, &
                           '         two discontinuities bounding the fluid!!'
                 write(6,*)procstrg,'   r,elem(loc,glob):',r1,iel,ielglob
                 write(6,*)procstrg,'Upper/lower discont:',&
                                                  discont(idom),discont(idom+1)
                 stop
              endif
              ! if current radius=bottom radius of fluid layer, set negative
              if (abs(r1-discont(idom+1))<min_distance_dim) then 
                 bdry_matr(:,iel,:)=-bdry_matr(:,iel,:)
                 count_lower_disc=count_lower_disc+1
              else ! element is in upper radius of fluid layer, keep positive
                 count_upper_disc=count_upper_disc+1
              endif
           endif
        enddo ! idom

! Factor r^2 stemming from the integration over spherical domain
      bdry_matr(0:npol,iel,1:2)=bdry_matr(0:npol,iel,1:2)*r*r
      solflubdry_radius(iel) = r

! NEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEW
! Adopt fluid formulation of Chaljub & Valette 2004, Komatitsch et al 2006: 
! u = \nabla \chi (leaving out the 1/rho factor)
! This choice merely results in factors of rho for all boundary terms
! and (supposedly?) allows for density stratification in the fluid
!!$
!!$  write(6,*)'Fluid:',r/1000.,iel,bdry_jpol_fluid(iel)
!!$  write(6,*)'Fluid:',bdry_fluid_el(iel),ielfluid(bdry_fluid_el(iel))
!!$  write(6,*)'Fluid minrho:',minval(rho(0:npol,bdry_jpol_fluid(iel), &
!!$                            ielfluid(bdry_fluid_el(iel))))
!!$  write(6,*)'Fluid maxrho:',maxval(rho(0:npol,bdry_jpol_fluid(iel), &
!!$                            ielfluid(bdry_fluid_el(iel)))) 
!!$  write(6,*)
!!$  write(6,*)'Solid:',r/1000.,theta/pi*180.,iel,bdry_jpol_solid(iel)
!!$  write(6,*)'Solid:',bdry_solid_el(iel),ielsolid(bdry_solid_el(iel))
!!$  write(6,*)'Solid minrho:',minval(rho(0:npol,bdry_jpol_solid(iel), &
!!$                            ielsolid(bdry_solid_el(iel))))
!!$  write(6,*)'Solid maxrho:',maxval(rho(0:npol,bdry_jpol_solid(iel), &
!!$                            ielsolid(bdry_solid_el(iel)))) 
!!$
!!$  if (iel==1) write(6,*)'Chaljub/Komatitsch potential formulation w/o density' 
!!$  do i=1,2
!!$     bdry_matr_fluid(0:npol,iel,i)=bdry_matr(0:npol,iel,i) / &
!!$          rho(0:npol,bdry_jpol_fluid(iel),ielfluid(bdry_fluid_el(iel)))
!!$
!!$     bdry_matr_solid(0:npol,iel,i)=bdry_matr(0:npol,iel,i) * &
!!$          rho(0:npol,bdry_jpol_solid(iel),ielsolid(bdry_solid_el(iel)))
!!$  enddo
! NEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEWNEW

     enddo ! elements along solid-fluid boundary

  endif ! have_bdry_elem

  bdry_sum = psum(real(bdry_sum,kind=realkind))

! yet another check....see if # elements above fluid is multiple of # below
  if (mod(count_upper_disc,2*count_lower_disc)/=0) then
     write(6,*)procstrg,&
               'Problem: Number of elements found to be at discont above fluid'
     write(6,*)procstrg,&
               '   is not an even multiple of elements found to be below fluid'
     write(6,*)procstrg,'# elems above fluid:',count_upper_disc
     write(6,*)procstrg,'# elems below fluid:',count_lower_disc
     stop
  endif

12 format(a25,3(1pe14.6))

  write(69,*)
  write(69,*)'# bdry elems above fluid:',count_upper_disc
  write(69,*)'# bdry elems below fluid:',count_lower_disc
  write(69,*)'Min/max of boundary term [m^2]:', &
                                           minval(bdry_matr),maxval(bdry_matr)
  write(69,*)'Min location bdry term  :',minloc(bdry_matr)
  write(69,*)'Max location bdry term  :',maxloc(bdry_matr)
  write(69,*)'Integr. bdry term, diff :',bdry_sum,dbleabsreldiff(bdry_sum,four)
  write(69,*)

  if ( .not. dblreldiff_small(bdry_sum,four) ) then
     if (lpr) then
        write(6,*)'WARNING: boundary term not all that precise!'
        write(6,*)' Term should equal four for 2 boundaries (i.e. int (sin) )'
        write(6,*)' Actual numerical value:',bdry_sum
     endif
     stop
  endif
  
!\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\
! output boundary precomputable matrix with radius [k]m and colatitude [deg]


  open(unit=500+mynum,file=infopath(1:lfinfo)//'/boundary_term_sol'&
                           //appmynum//'.dat')
  open(unit=400+mynum,file=infopath(1:lfinfo)//'/boundary_term_flu'&
                           //appmynum//'.dat')
  write(69,*)
  write(69,*)'saving boundary matrix with solid radius/colatitude into ',&
              'boundary_term_sol_'//appmynum//'.dat'
  write(69,*)'saving boundary matrix with fluid radius/colatitude into ',&
              'boundary_term_flu_'//appmynum//'.dat'
  write(69,*)

  do iel=1,nel_bdry  
     ielglob=ielsolid(bdry_solid_el(iel))
     do ipol=0,npol
      write(500+mynum,15)rcoord(ipol,bdry_jpol_solid(iel),ielglob)/1.d3, &
                        thetacoord(ipol,bdry_jpol_solid(iel),ielglob)/pi*180.,&
                         bdry_matr(ipol,iel,1),bdry_matr(ipol,iel,2)

      write(400+mynum,15)rcoord(ipol,bdry_jpol_fluid(iel),&
                                ielfluid(bdry_fluid_el(iel)))/1.d3, &
                         thetacoord(ipol,bdry_jpol_fluid(iel),&
                                    ielfluid(bdry_fluid_el(iel)))/pi*180.,&
                         bdry_matr(ipol,iel,1),bdry_matr(ipol,iel,2)
     enddo
  enddo
  call flush(500+mynum); close(500+mynum)
  call flush(400+mynum); close(400+mynum)

15 format(4(1pe12.4))
!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

end subroutine def_solid_fluid_boundary_terms
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=========================
end module def_precomp_terms
!=========================
