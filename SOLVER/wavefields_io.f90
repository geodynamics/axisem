!=========================
module wavefields_io
!=========================
!
! Contains all routines that dump entire wavefields during the time loop. 
! Optimization of I/O therefore happens here and nowhere else.
! The corresponding meshes are dumped in meshes_io.
!
use global_parameters
use data_mesh
use data_proc
use data_io

implicit none

public 

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine glob_snapshot(f_sol,chi,ibeg,iend,jbeg,jend)
!
! Dumps the global displacement snapshots [m] in ASCII format
! When reading the fluid wavefield, one needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord
! as dumped by the corresponding routine dump_glob_grid!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f_sol(0:npol,0:npol,1:nel_solid,3)
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,iidim
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call define_io_appendix(appisnap,isnap)

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/snap_'&
                            //appmynum//'_'//appisnap//'.dat')

  if (have_fluid) then
     call axisym_laplacian_fluid(chi,usz_fluid)
     do iel=1,nel_fluid

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,0,dsdchi)
                 write(2500+mynum,*)usz_fluid(0,0,iel,1),prefac*dsdchi*&
                                  chi(0,0,iel),usz_fluid(0,0,iel,2)
              else
                 write(2500+mynum,*)usz_fluid(0,0,iel,1), &
                                    prefac*chi(0,0,iel), &
                                    usz_fluid(0,0,iel,2)
              endif

                 write(2500+mynum,*)usz_fluid(npol,0,iel,1), &
                                    prefac*chi(npol,0,iel), &
                                    usz_fluid(npol,0,iel,2)

                 write(2500+mynum,*)usz_fluid(npol,npol,iel,1), &
                                    prefac*chi(npol,npol,iel), &
                                    usz_fluid(npol,npol,iel,2)

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,npol,dsdchi)
                 write(2500+mynum,*)usz_fluid(0,npol,iel,1),prefac*dsdchi*&
                                  chi(0,npol,iel),usz_fluid(0,npol,iel,2)
              else
                 write(2500+mynum,*)usz_fluid(0,npol,iel,1), &
                                    prefac*chi(0,npol,iel), &
                                    usz_fluid(0,npol,iel,2)
              endif

     enddo
  endif ! have_fluid

  do iel=1,nel_solid
           write(2500+mynum,*)(f_sol(0,0,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(npol,0,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(npol,npol,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(0,npol,iel,iidim),iidim=1,3)
  enddo

  close(2500+mynum)


!  h_real=real(hmax/(period/(pts_wavelngth*real(npol))))
!  fname=trim(diagpath)//'/mesh_hmax'
!  call write_VTK_bin_scal(h_real,mesh2,neltot,fname)


end subroutine glob_snapshot
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_snapshot_midpoint(f_sol,chi,ibeg,iend,jbeg,jend)
!
! Dumps the global displacement snapshots [m] in ASCII format
! When reading the fluid wavefield, one needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord
! as dumped by the corresponding routine dump_glob_grid!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f_sol(0:npol,0:npol,1:nel_solid,3)
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,iidim
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call define_io_appendix(appisnap,isnap)

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/snap_'&
                            //appmynum//'_'//appisnap//'.dat' ,FORM="UNFORMATTED",STATUS="REPLACE")

  if (have_fluid) then
     call axisym_laplacian_fluid(chi,usz_fluid)
     do iel=1,nel_fluid
        do jpol=0,npol,npol/2
           do ipol=0,npol,npol/2

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,jpol,dsdchi)
                 write(2500+mynum)usz_fluid(ipol,jpol,iel,1),prefac*dsdchi*&
                                  chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
              else
                 write(2500+mynum)usz_fluid(ipol,jpol,iel,1), &
                                    prefac*chi(ipol,jpol,iel), &
                                    usz_fluid(ipol,jpol,iel,2)
              endif
              enddo
           enddo
     enddo
  endif ! have_fluid

  do iel=1,nel_solid
     do jpol=0,npol,npol/2
        do ipol=0,npol,npol/2
           write(2500+mynum)(f_sol(ipol,jpol,iel,iidim),iidim=1,3)
        enddo
     enddo
  enddo
  close(2500+mynum)

!  h_real=real(hmax/(period/(pts_wavelngth*real(npol))))
!  fname=trim(diagpath)//'/mesh_hmax'
!  call write_VTK_bin_scal(h_real,mesh2,neltot,fname)


end subroutine glob_snapshot_midpoint
!=============================================================================


!-----------------------------------------------------------------------------
subroutine solid_snapshot(f,ibeg,iend,jbeg,jend)
!
! Dumps the displacement snapshots [m] in the solid region in ASCII format
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f(0:npol,0:npol,1:nel_solid,3)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,idim

  call define_io_appendix(appisnap,isnap)

  open(unit=3500+mynum,file=datapath(1:lfdata)//'/snap_solid_'&
                            //appmynum//'_'//appisnap//'.dat')

  do iel=1,nel_solid
     do jpol=ibeg,iend
        do ipol=jbeg,jend
           write(3500+mynum,*)(f(ipol,jpol,iel,idim),idim=1,3)
        enddo
     enddo
  enddo
  close(3500+mynum)

end subroutine solid_snapshot
!=============================================================================

!-----------------------------------------------------------------------------

subroutine fluid_snapshot(chi,ibeg,iend,jbeg,jend)

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel,ipol,jpol
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call axisym_laplacian_fluid(chi,usz_fluid)

  call define_io_appendix(appisnap,isnap)

  open(unit=4500+mynum,file=datapath(1:lfdata)//'/snap_fluid_'&
                            //appmynum//'_'//appisnap//'.dat')

  do iel=1,nel_fluid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           if ( axis_fluid(iel) .and. ipol==0 ) then
              call dsdf_fluid_axis(chi(:,:,iel),iel,jpol,dsdchi)
              write(4500+mynum,*)usz_fluid(ipol,jpol,iel,1),prefac*dsdchi* &
                                chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
           else
              write(4500+mynum,*)usz_fluid(ipol,jpol,iel,1),prefac*&
                                 chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
           endif
        enddo
     enddo
  enddo

  close(4500+mynum)

end subroutine fluid_snapshot
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_field_1d(f,filename,appisnap,n)

use data_proc, ONLY : appmynum
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

integer, intent(in) :: n
real(kind=realkind),intent(inout) :: f(0:npol,0:npol,1:n)
character(len=16), intent(in)     :: filename
character(len=4), intent(in)      :: appisnap
real(kind=realkind), allocatable :: f1(:)
integer :: i,j,iel,ii

  open(unit=25000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="UNKNOWN",POSITION="REWIND")

  if (have_src .and. src_dump_type == 'mask' .and. n==nel_solid) &
       call eradicate_src_elem_values(f)

allocate(f1((iend-ibeg+1)**2*n))

  ii=0
  do iel=1,n
     do j=ibeg,iend
        do i=ibeg,iend
           ii=ii+1
           f1(ii)=f(i,j,iel)
        enddo
     enddo
  enddo

  write(25000+mynum) f1
    deallocate(f1)
  close(25000+mynum)

end subroutine dump_field_1d
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_field_over_s_solid_1d(f,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_solid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_solid)
integer                        :: iel,i
!real(kind=realkind) :: f1d((npol+1)**2*nel_solid)
real(kind=realkind)            :: g(0:npol,0:npol,nel_solid)


  open(unit=35000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

  g = inv_s_solid*f

  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(g)

  write(35000+mynum) g(ibeg:iend,ibeg:iend,:)

  close(35000+mynum)

end subroutine dump_field_over_s_solid_1d
!=============================================================================

!--------------------------------------------------------------------------
subroutine dump_field_over_s_solid_and_add(f,g,filename1,filename2,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the solid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_solid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_solid)
character(len=16), intent(in)  :: filename1,filename2
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_solid)
integer                        :: iel,i

  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

! construct masked f/s (e.g. Epp)
  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(inv_s_solid)

  open(unit=39000+mynum,file=datapath(1:lfdata)//filename1//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
  write(39000+mynum) inv_s_solid(ibeg:iend,ibeg:iend,:)* &
                     f(ibeg:iend,ibeg:iend,:)
  close(39000+mynum)

! construct sum of f/s and g (e.g. straintrace)
  g = inv_s_solid*f + g

  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(g)

  open(unit=35000+mynum,file=datapath(1:lfdata)//filename2//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
  write(35000+mynum) g(ibeg:iend,ibeg:iend,:)
  close(35000+mynum)

end subroutine dump_field_over_s_solid_and_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_half_field_over_s_solid_1d_add(f,g,filename,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the solid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in)    :: f(0:npol,0:npol,nel_solid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_solid)
character(len=16), intent(in)     :: filename
character(len=4), intent(in)      :: appisnap
real(kind=realkind)               :: dsdf(0:npol,naxel_solid)
integer                           :: iel

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

  g = real(.5,kind=realkind) * ( inv_s_solid * f + g )

  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(g)

  open(unit=35000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
  write(35000+mynum) g(ibeg:iend,ibeg:iend,:)
  close(35000+mynum)

end subroutine dump_half_field_over_s_solid_1d_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_field_over_s_fluid_and_add(f,g,filename1,filename2,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the fluid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid!,deviator
use pointwise_derivatives, ONLY: dsdf_fluid_allaxis
use data_source, ONLY : src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename1,filename2
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel,i

  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_fluid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_fluid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

! f/s (e.g. Epp)
  open(unit=39000+mynum,file=datapath(1:lfdata)//filename1//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
  write(39000+mynum) inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                     f(ibeg:iend,ibeg:iend,:)
  close(39000+mynum)

!  deviator(ibeg:iend,ibeg:iend,:,2) = inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
!                                      f(ibeg:iend,ibeg:iend,:)


! sum of f/s and g (e.g. straintrace)
  open(unit=35000+mynum,file=datapath(1:lfdata)//filename2//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
  write(35000+mynum) inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                     f(ibeg:iend,ibeg:iend,:) + g(ibeg:iend,ibeg:iend,:)
  close(35000+mynum)
  
!  deviator(ibeg:iend,ibeg:iend,:,3) = inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
!                                      f(ibeg:iend,ibeg:iend,:) + &
!                                      g(ibeg:iend,ibeg:iend,:)

end subroutine dump_field_over_s_fluid_and_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_half_f1_f2_over_s_fluid(f1,f2,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid
use pointwise_derivatives, ONLY : dsdf_fluid_allaxis
include 'mesh_params.h'

real(kind=realkind),intent(in) :: f1(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(in) :: f2(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel

  if (have_axis) then
     call dsdf_fluid_allaxis(f2,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

  open(unit=65000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  write(65000+mynum) 0.5*(f1(ibeg:iend,ibeg:iend,:) + &
                     inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                     f2(ibeg:iend,ibeg:iend,:))

  close(65000+mynum)

end subroutine dump_half_f1_f2_over_s_fluid
!=============================================================================

!--------------------------------------------------------------------------
subroutine dump_f1_f2_over_s_fluid(f1,f2,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid
use pointwise_derivatives, ONLY : dsdf_fluid_allaxis
include 'mesh_params.h'

real(kind=realkind),intent(in) :: f1(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(in) :: f2(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel

  if (have_axis) then
     call dsdf_fluid_allaxis(f2,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

  open(unit=65000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  write(65000+mynum) f1(ibeg:iend,ibeg:iend,:) + &
                     inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                     f2(ibeg:iend,ibeg:iend,:)

  close(65000+mynum)

end subroutine dump_f1_f2_over_s_fluid
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_disp(u,chi)

use data_source, ONLY : src_type,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: u(0:npol,0:npol,nel_solid,3)
real(kind=realkind),intent(in) :: chi(0:npol,0:npol,nel_fluid)
integer                        :: i
character(len=4)               :: appisnap
real(kind=realkind)            :: f(0:npol,0:npol,nel_solid,3)

  call define_io_appendix(appisnap,istrain)

! Dump solid displacement
  open(unit=75000+mynum,file=datapath(1:lfdata)//'/disp_sol_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  f = u

  if (src_dump_type == 'mask') &
       call eradicate_src_elem_vec_values(f)

  if (src_type(1)/='monopole') then
     write(75000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
  else
     write(75000+mynum) f(ibeg:iend,ibeg:iend,:,1:3:2)
  endif
  close(75000+mynum)

! Dump fluid potential 
  if (have_fluid) then 
     open(unit=76000+mynum,file=datapath(1:lfdata)//'/chi_flu_'&
                                //appmynum//'_'//appisnap//'.bindat',&
                                FORM="UNFORMATTED",STATUS="REPLACE")
  
     write(76000+mynum)chi
     close(76000+mynum)
  endif 

end subroutine dump_disp
!=============================================================================

!--------------------------------------------------------------------------
subroutine dump_velo_dchi(v,dchi)

use data_source, ONLY : src_type,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: v(0:npol,0:npol,nel_solid,3)
real(kind=realkind),intent(in) :: dchi(0:npol,0:npol,nel_fluid)
integer                        :: i
character(len=4)               :: appisnap
real(kind=realkind)            :: f(0:npol,0:npol,nel_solid,3)

  call define_io_appendix(appisnap,istrain)

! Dump solid velocity vector
  open(unit=85000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  f = v

  if (src_dump_type == 'mask') &
       call eradicate_src_elem_vec_values(f)

  if (src_type(1)/='monopole') then 
     write(85000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
  else
     write(85000+mynum) f(ibeg:iend,ibeg:iend,:,1),f(ibeg:iend,ibeg:iend,:,3)
  endif
  close(85000+mynum)

! Dump fluid potential 1st derivative
  if (have_fluid) then 
     open(unit=86000+mynum,file=datapath(1:lfdata)//'/dchi_flu_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
  
     write(86000+mynum)dchi
     close(86000+mynum)
  endif

end subroutine dump_velo_dchi
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_velo_global(v,dchi)

use data_pointwise, ONLY: inv_rho_fluid,inv_s_rho_fluid,usz_fluid
use data_source, ONLY : src_type,src_dump_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_allaxis
use unit_stride_colloc, ONLY : collocate0_1d

include 'mesh_params.h'

real(kind=realkind),intent(in) :: v(:,:,:,:)
real(kind=realkind),intent(in) :: dchi(:,:,:)

real(kind=realkind)            :: phicomp(0:npol,0:npol,nel_fluid)
integer                        :: i,iel
character(len=4)               :: appisnap
real(kind=realkind)            :: dsdchi(0:npol,naxel_fluid)
real(kind=realkind)            :: f(0:npol,0:npol,1:nel_solid,3)

  call define_io_appendix(appisnap,istrain)

! sssssssssssss dump velocity vector inside solid ssssssssssssssssssssssssssss

  open(unit=95000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")

  f=v
  if (src_dump_type == 'mask') &
       call eradicate_src_elem_vec_values(f)

  if (src_type(1)/='monopole') then
     write(95000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
  else
     write(95000+mynum) f(ibeg:iend,ibeg:iend,:,1), &
                        f(ibeg:iend,ibeg:iend,:,3)
  endif
  close(95000+mynum)

! ffffffff fluid region ffffffffffffffffffffffffffffffffffffffffffffffffffffff

  if (have_fluid) then 
! compute velocity vector inside fluid
     call axisym_laplacian_fluid(dchi,usz_fluid)

! phi component needs special care: m/(s rho) dchi
     call collocate0_1d(inv_s_rho_fluid,dchi,phicomp,npoint_fluid)

! Take care of axial singularity for phi component of fluid velocity
     if (have_axis) then
        call dsdf_fluid_allaxis(dchi,dsdchi)
        do iel=1,naxel_fluid
           phicomp(0,:,ax_el_fluid(iel))=dsdchi(:,iel)* &
                phicomp(0,:,ax_el_fluid(iel))
        enddo
     endif

  call define_io_appendix(appisnap,istrain)

! dump velocity vector inside fluid
     open(unit=960000+mynum,file=datapath(1:lfdata)//'/velo_flu_'&
                               //appmynum//'_'//appisnap//'.bindat',&
                               FORM="UNFORMATTED",STATUS="REPLACE")

     write(960000+mynum) inv_rho_fluid(ibeg:iend,ibeg:iend,:)* &
                         usz_fluid(ibeg:iend,ibeg:iend,:,1),phicomp, &
                         inv_rho_fluid(ibeg:iend,ibeg:iend,:)* &
                         usz_fluid(ibeg:iend,ibeg:iend,:,2)
     close(960000+mynum)
  endif ! have_fluid

end subroutine dump_velo_global
!=============================================================================


!-----------------------------------------------------------------------------
subroutine eradicate_src_elem_vec_values(u)
!
! Deletes all entries to vector field u on ALL GLL points inside
! elements that have a non-zero source term (i.e. including all 
! assembled neighboring elements)
! This is a preliminary test for the wavefield dumps.

use data_source, ONLY : nelsrc,ielsrc,have_src

include 'mesh_params.h'

real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel_solid,3)
integer :: iel

if (have_src) then
   do iel=1,nelsrc
      u(0:npol,0:npol,ielsrc(iel),1:3) = real(0.,kind=realkind)
   enddo
endif

end subroutine eradicate_src_elem_vec_values
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------------
subroutine eradicate_src_elem_values(u)
!
! Deletes all entries to scalar field u on ALL GLL points inside
! elements that have a non-zero source term (i.e. including all 
! assembled neighboring elements)
! This is a preliminary test for the wavefield dumps.

use data_source, ONLY : nelsrc,ielsrc,have_src

include 'mesh_params.h'

real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel_solid)
integer :: iel

if (have_src) then

   do iel=1,nelsrc
      u(0:npol,0:npol,ielsrc(iel)) = real(0.,kind=realkind)
   enddo
endif


end subroutine eradicate_src_elem_values
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal(u2,mesh,rows,filename)

 implicit none
 integer*4 :: i,rows
 real, dimension(1:rows), intent(in) :: u2
 real, dimension(1:rows) :: u1
real, dimension(1:rows,2), intent(in) :: mesh
 integer*4, dimension(1:rows*2) :: cell
 integer*4, dimension(1:rows) :: cell_type
  character (len=55) :: filename;
 character (len=50) :: ss; !stream
 
! write(6,*) size(W_vtk(:,1)),rows
!points structure

do i=2,rows*2,2
 cell(i-1)=1;
 cell(i)=(i/2)-1;
enddo
do i=1,rows
 cell_type(i)=1
enddo

u1=real(u2)
do i=1,rows
    if (abs(u1(i))<1.e-25) u1(i)=0.0
enddo
write(6789,*)size(u1),maxval(u1),minval(u1)

 write(6,*)'computing vtk file ',trim(filename),' ...'
open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')

write(100) '# vtk DataFile Version 4.0'//char(10)
write(100) 'mittico'//char(10)
write(100) 'BINARY'//char(10)
write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
write(100) ss//char(10)
!points
do i=1,rows
write(100) mesh(i,1),0.0,mesh(i,2)
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
end subroutine write_vtk_bin_scal
!-----------------------------------------------------------------------------


!================================
end module wavefields_io
!================================
