!=========================
module meshes_io
!=========================
!
! This module contains routines that compute and dump the respective meshes
! underlying the actual wavefields to be dumped in the time loop
! which is done in wavefields_io. This module needs pre-loop variables such
! as mesh coordinates and is therefore cut off from the dumping module.
!
use global_parameters
use data_mesh
use data_mesh_preloop, ONLY : ielsolid,ielfluid
use data_proc
use data_io

use utlity, ONLY : scoord, zcoord

implicit none

public

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------------
subroutine dump_test

use data_proc, ONLY : appmynum
use utlity, only : rcoord
include 'mesh_params.h'

double precision :: f(0:npol,0:npol,1:nel_fluid)
double precision :: g(0:npol,0:npol,1:nel_solid)
integer :: ipol,jpol,iel

if (lpr) write(6,*)'TEST:',nelem,nel_fluid,nel_solid

  do iel=1,nel_fluid
     do jpol=0,npol-1
        do ipol=0,npol-1
            f(ipol,jpol,iel) = &
                   2.d0*dsin(8.d0*pi*rcoord(ipol,jpol,ielfluid(iel))/router)
        enddo
     enddo
  enddo

  open(unit=25000+mynum,file=datapath(1:lfdata)//'dump_testflu'//'_'&
                             //appmynum//'.bindat',&
                            FORM="UNFORMATTED",STATUS="NEW")
  write(25000+mynum) f(ibeg:iend,ibeg:iend,:)

  close(25000+mynum)

  do iel=1,nel_solid
     do jpol=0,npol-1
        do ipol=0,npol-1
            g(ipol,jpol,iel) = &
                   2.d0*dsin(8.d0*pi*rcoord(ipol,jpol,ielsolid(iel))/router)
        enddo
     enddo
  enddo
  open(unit=25100+mynum,file=datapath(1:lfdata)//'dump_testsol'//'_'&
                            //appmynum//'.bindat',&
                            FORM="UNFORMATTED",STATUS="NEW")

  write(25100+mynum) g(ibeg:iend,ibeg:iend,:)    

  close(25100+mynum)

stop

end subroutine dump_test
!=============================================================================


!-----------------------------------------------------------------------------
subroutine dump_glob_grid(ibeg,iend,jbeg,jend)
!
! Dumps the mesh (s,z) [m] in ASCII format as needed to visualize global 
! snapshots, and additionally the constant factors preceding the displacement
! in the fluid, namely rho^{-1} and (rho s)^{-1}.
! When reading the fluid wavefield, one therefore needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : inv_rho_fluid

include 'mesh_params.h'

integer, intent(in) :: ibeg,iend,jbeg,jend 
integer             :: iel,ipol,jpol

  open(unit=25000+mynum,file=datapath(1:lfdata)//'/glob_grid_'&
                             //appmynum//'.dat')

  if (have_fluid) then 
     open(unit=26000+mynum,file=datapath(1:lfdata)// &
                                '/inv_rho_s_fluid_globsnaps_' &
                                //appmynum//'.dat')
     do iel=1,nel_fluid

           if ( axis_fluid(iel)  ) then
!       Axis s=0! write 1 instead of 1/s and then multiply
!       with the correct factor dsdchi, obtained by L'Hospital's rule
!       (see routine glob_snapshot in wavefields_io).
                 write(26000+mynum,*)inv_rho_fluid(0,0,iel),one
                 write(26000+mynum,*)inv_rho_fluid(npol,0,iel), &
                                     one/scoord(npol,0,ielfluid(iel))
                 write(26000+mynum,*)inv_rho_fluid(npol,npol,iel), &
                                     one/scoord(npol,npol,ielfluid(iel))
                 write(26000+mynum,*)inv_rho_fluid(0,npol,iel),one

              else
                 write(26000+mynum,*)inv_rho_fluid(0,0,iel), &
                                     one/scoord(0,0,ielfluid(iel))
                 write(26000+mynum,*)inv_rho_fluid(npol,0,iel), &
                                     one/scoord(npol,0,ielfluid(iel))
                 write(26000+mynum,*)inv_rho_fluid(npol,npol,iel), &
                                     one/scoord(npol,npol,ielfluid(iel))
                 write(26000+mynum,*)inv_rho_fluid(0,npol,iel), &
                                     one/scoord(0,npol,ielfluid(iel))
              endif

              write(25000+mynum,*)scoord(0,0,ielfluid(iel)), &
                                  zcoord(0,0,ielfluid(iel))
              write(25000+mynum,*)scoord(npol,0,ielfluid(iel)), &
                                  zcoord(npol,0,ielfluid(iel))
              write(25000+mynum,*)scoord(npol,npol,ielfluid(iel)), &
                                  zcoord(npol,npol,ielfluid(iel))
             write(25000+mynum,*)scoord(0,npol,ielfluid(iel)), &
                                  zcoord(0,npol,ielfluid(iel))

     enddo
     close(26000+mynum)
  endif ! have_fluid

  do iel=1,nel_solid

           write(25000+mynum,*)scoord(0,0,ielsolid(iel)), &
                               zcoord(0,0,ielsolid(iel))
           write(25000+mynum,*)scoord(npol,0,ielsolid(iel)), &
                               zcoord(npol,0,ielsolid(iel))
           write(25000+mynum,*)scoord(npol,npol,ielsolid(iel)), &
                               zcoord(npol,npol,ielsolid(iel))
           write(25000+mynum,*)scoord(0,npol,ielsolid(iel)), &
                               zcoord(0,npol,ielsolid(iel))
  enddo
  close(25000+mynum)

end subroutine dump_glob_grid
!=============================================================================



!-----------------------------------------------------------------------------
subroutine dump_glob_grid_midpoint(ibeg,iend,jbeg,jend)
!
! Dumps the mesh (s,z) [m] in ASCII format as needed to visualize global 
! snapshots, and additionally the constant factors preceding the displacement
! in the fluid, namely rho^{-1} and (rho s)^{-1}.
! When reading the fluid wavefield, one therefore needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : inv_rho_fluid

include 'mesh_params.h'

integer, intent(in) :: ibeg,iend,jbeg,jend 
integer             :: iel,ipol,jpol

  open(unit=25000+mynum,file=datapath(1:lfdata)//'/glob_grid_'&
                             //appmynum//'.dat')

  if (have_fluid) then 
     open(unit=26000+mynum,file=datapath(1:lfdata)// &
                                '/inv_rho_s_fluid_globsnaps_' &
                                //appmynum//'.dat')
     do iel=1,nel_fluid
        do jpol=0,npol,npol/2
           do ipol=0,npol,npol/2

           if ( axis_fluid(iel) .and. ipol==0 ) then
!       Axis s=0! write 1 instead of 1/s and then multiply
!       with the correct factor dsdchi, obtained by L'Hospital's rule
!       (see routine glob_snapshot in wavefields_io).
                 write(26000+mynum,*)inv_rho_fluid(ipol,jpol,iel),one
              else
                 write(26000+mynum,*)inv_rho_fluid(ipol,jpol,iel), &
                                     one/scoord(ipol,jpol,ielfluid(iel))
              endif

              write(25000+mynum,*)scoord(ipol,jpol,ielfluid(iel)), &
                                  zcoord(ipol,jpol,ielfluid(iel))
           enddo
        enddo
     enddo
     close(26000+mynum)
  endif ! have_fluid

  do iel=1,nel_solid
     do jpol=0,npol,npol/2
        do ipol=0,npol,npol/2
           write(25000+mynum,*)scoord(ipol,jpol,ielsolid(iel)), &
                               zcoord(ipol,jpol,ielsolid(iel))
           enddo
        enddo
  enddo
  close(25000+mynum)

end subroutine dump_glob_grid_midpoint
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_solid_grid(ibeg,iend,jbeg,jend)
!
! Dumps the mesh (s,z) [m] in ASCII format as needed to visualize snapshots 
! in the solid region only.
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include 'mesh_params.h'

integer, intent(in) :: ibeg,iend,jbeg,jend 
integer             :: iel, ipol,jpol

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/solid_grid_'&
                            //appmynum//'.dat')
  do iel=1,nel_solid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           write(2500+mynum,*)scoord(ipol,jpol,ielsolid(iel)), &
                              zcoord(ipol,jpol,ielsolid(iel))
        enddo
     enddo
  enddo
  close(2500+mynum)

end subroutine dump_solid_grid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_fluid_grid(ibeg,iend,jbeg,jend)
!
! Dumps the mesh (s,z) [m] in ASCII format as needed to visualize snapshots 
! in the fluid region only, and additionally the constant factors preceding 
! the displacement in the fluid, namely rho^{-1} and (rho s)^{-1}.
! When reading the fluid wavefield, one therefore needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : inv_rho_fluid

include 'mesh_params.h'

integer, intent(in) :: ibeg,iend,jbeg,jend
integer             :: iel, ipol,jpol

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  open(unit=2500+mynum,file=datapath(1:lfdata)//&
                            '/fluid_grid_'//appmynum//'.dat')
  open(unit=2600+mynum,file=datapath(1:lfdata)//&
                            '/inv_rho_scoord_fluid_flusnaps_'&
                            //appmynum//'.dat', STATUS="REPLACE")
  do iel=1,nel_fluid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           write(2500+mynum,*)scoord(ipol,jpol,ielfluid(iel)), &
                              zcoord(ipol,jpol,ielfluid(iel))
           if ( axis_fluid(iel) .and. ipol==0 ) then
!           Axis s=0! write 1 instead of 1/s and then multiply 
!           with the correct factor dsdchi, obtained by L'Hospital's rule 
!           (see routine fluid_snapshot below).
              write(2600+mynum,*)inv_rho_fluid(ipol,jpol,iel),one
           else  
              write(2600+mynum,*)inv_rho_fluid(ipol,jpol,iel), &
                                 one/scoord(ipol,jpol,ielfluid(iel))
           endif
        enddo
     enddo
  enddo
  close(2500+mynum)
  close(2600+mynum)

end subroutine dump_fluid_grid
!=============================================================================


!-----------------------------------------------------------------------------
subroutine dump_wavefields_mesh_1d
!
! Dumps the mesh (s,z) [m] and related constant fields in binary format as 
! needed to compute waveform kernels from the strain and velocity fields. 
! The distinction between different dumping methods is honored here, 
! and influences the amount of additional dumpsters (prefactors, derivatives). 
! In a nutshell, the end-member dumping methods constitute 
! 1) computing strain and velocity on-the-fly, i.e. only dumping the mesh here;
! 2) only dumping the already-known fields (displacement, potential) on-the-fly
!    and dump a number of constant fields here. 
! The latter choice is more memory- and CPU-efficient, but requires 
! significant post-processing AND dumping the entire SEM mesh. 
! See compute_strain in time_evol_wave.f90 for more info.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_io, ONLY : ibeg,iend,ndumppts_el
use data_spec, ONLY : G1T,G2T,G2
use data_pointwise

!double precision, dimension(:), allocatable :: ssol, zsol
!double precision, dimension(:), allocatable :: sflu, zflu
double precision, dimension(:,:,:), allocatable :: ssol, zsol
double precision, dimension(:,:,:), allocatable :: sflu, zflu

integer :: iel, ipol,jpol,i

! Dump entire (including duplicate) GLL point grid for displ_only
  if (dump_type=='displ_only') then
     
  else ! Free choice for other dumping method

     if (lpr) then
        write(6,*)'  set strain dumping GLL boundaries to:'
        write(6,*)'    ipol=',ibeg,iend

     endif
!     allocate(ssol(ndumppts_el*nel_solid))
!     allocate(zsol(ndumppts_el*nel_solid))
!     allocate(sflu(ndumppts_el*nel_fluid))
!     allocate(zflu(ndumppts_el*nel_fluid))

     allocate(ssol(ibeg:iend,ibeg:iend,nel_solid))
     allocate(zsol(ibeg:iend,ibeg:iend,nel_solid))
     allocate(sflu(ibeg:iend,ibeg:iend,nel_fluid))
     allocate(zflu(ibeg:iend,ibeg:iend,nel_fluid))

  endif

i=0

! compute solid grid
  do iel=1,nel_solid
      do jpol=ibeg,iend
        do ipol=ibeg,iend
!           i=i+1
!           ssol(i)=scoord(ipol,jpol,ielsolid(iel))
!           zsol(i)=zcoord(ipol,jpol,ielsolid(iel))
           
           ssol(ipol,jpol,iel) = scoord(ipol,jpol,ielsolid(iel))
           zsol(ipol,jpol,iel) = zcoord(ipol,jpol,ielsolid(iel))
        enddo
     enddo
  enddo


  if (lpr) &
  write(6,*)'  dumping solid submesh for kernel wavefields...'
  open(unit=2500+mynum,file=datapath(1:lfdata)//'/strain_mesh_sol_'&
                            //appmynum//'.dat', &
                            FORM="UNFORMATTED",STATUS="REPLACE")

  write(2500+mynum)ssol(ibeg:iend,ibeg:iend,:),zsol(ibeg:iend,ibeg:iend,:)
  close(2500+mynum)
  deallocate(ssol,zsol)

  if (have_fluid) then

i=0

! compute fluid grid
     do iel=1,nel_fluid
       do jpol=ibeg,iend
         do ipol=ibeg,iend
!              i=i+1
!              sflu(i)=scoord(ipol,jpol,ielfluid(iel))
!              zflu(i)=zcoord(ipol,jpol,ielfluid(iel))
            sflu(ipol,jpol,iel) = scoord(ipol,jpol,ielfluid(iel))
            zflu(ipol,jpol,iel) = zcoord(ipol,jpol,ielfluid(iel))
           enddo
        enddo
     enddo
     if (lpr) &
     write(6,*)'  dumping fluid submesh for kernel wavefields...'
     open(unit=2600+mynum,file=datapath(1:lfdata)//'/strain_mesh_flu_'&
          //appmynum//'.dat', &
          FORM="UNFORMATTED",STATUS="REPLACE")
     write(2600+mynum)sflu(ibeg:iend,ibeg:iend,:),zflu(ibeg:iend,ibeg:iend,:)
     close(2600+mynum)
     deallocate(sflu,zflu)
  endif ! have_fluid


! In the following: Only dumping additional arrays if displacements only 
! are dumped as wavefields to reconstruct the strains.

  select case (dump_type)
  case ('displ_only')
     if (lpr) then     
        write(6,*)'  strain dump: only displacement/velocity, potentials'
        write(6,*)'  ...now dumping global pointwise deriv. terms, etc....'
     endif

!    Dump pointwise derivative matrices in solid
     open(unit=2600+mynum,file=datapath(1:lfdata)//'/pointwise_deriv_sol_'&
                               //appmynum//'.dat', &
                               FORM="UNFORMATTED",STATUS="REPLACE")
     write(2600+mynum)DzDeta_over_J_sol,DzDxi_over_J_sol, &
                      DsDeta_over_J_sol,DsDxi_over_J_sol
     close(2600+mynum)

     if (have_fluid) then
!    Dump pointwise derivative matrices in fluid
        open(unit=2600+mynum,file=datapath(1:lfdata)//'/pointwise_deriv_flu_'&
             //appmynum//'.dat', &
             FORM="UNFORMATTED",STATUS="REPLACE")
        write(2600+mynum)DzDeta_over_J_flu,DzDxi_over_J_flu, &
             DsDeta_over_J_flu,DsDxi_over_J_flu
        close(2600+mynum)

!    Dump inverse density inside fluid
        open(unit=2600+mynum,file=datapath(1:lfdata)//'/inv_rho_fluid_'&
             //appmynum//'.dat', &
             FORM="UNFORMATTED",STATUS="REPLACE")
        write(2600+mynum)inv_rho_fluid
        close(2600+mynum)
     endif

!    Dump Lagrange interpolant derivatives
     open(unit=2600+mynum,file=datapath(1:lfdata)//'/lagrange_derivs_'&
                               //appmynum//'.dat', &
                               FORM="UNFORMATTED",STATUS="REPLACE")
     write(2600+mynum)G1T,G2T,G2
     close(2600+mynum)

     write(6,*)'  ...dumped it all.'

  case ('fullfields')
     if (lpr) then
        write(6,*)'  strain dump: Global strain tensor and velocity fields'
        write(6,*)'  ....no need to dump anything else.'
     endif
  case default
     if (lpr) then 
        write(6,*)'  wavefield dumping type',dump_type,' unknown!'
        write(6,*)'  select from 1) displ_only, 2) fullfields'
     endif
     stop
  end select

end subroutine dump_wavefields_mesh_1d
!=============================================================================

!-----------------------------------------------------------------------------
subroutine fldout_cyl2(fname,nel,f,ibeg,iend,jbeg,jend,flag_norm,domain)
! 
! Dumps the mesh (s,z) [m] along with corresponding field f in ASCII format.
! At this point only used in computing the valence. 
! Note that for higher-frequency meshes these files can be very large.
! flag_norm is to be set to one if the output is to be normalized.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use utlity, ONLY : compute_coordinates

include 'mesh_params.h'

character(len=80), intent(in)   :: fname
character(len=5), intent(in)    :: domain
integer, intent(in)             :: flag_norm,ibeg,iend,jbeg,jend,nel
real(kind=realkind), intent(in) :: f(ibeg:iend,jbeg:jend,nel)
integer                         :: lf,ielem, ipol,jpol,iel
double precision                :: r, theta, s, z
real(kind=realkind)             :: fnr,afnr
  
  lf=index(fname,' ')-1
  open(unit=10000+mynum,file=infopath(1:lfinfo)//'/'//fname(1:lf)//'_'&
                             //appmynum//'.dat')

  fnr = 1.
  if (flag_norm == 1) then
     fnr = zero
     do ielem = 1, nel
        do jpol = jbeg, jend
           do ipol = ibeg, iend
              fnr = max(fnr,abs(f(ielem,ipol,jpol)))
           end do
        end do
     end do   
  end if
!  if (nproc > 1) fnr=pmax(dble(fnr))
  afnr = one/fnr
  do ielem = 1, nel
     if (domain=='total') iel=ielem
     if (domain=='solid') iel=ielsolid(ielem)
     if (domain=='fluid') iel=ielfluid(ielem)
     do jpol = jbeg, jend !,npol/2
        do ipol = ibeg, iend !,npol/2
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
           write(10000+mynum,*) s/router,z/router,afnr*f(ipol,jpol,ielem)
        end do
     end do
  end do
  close(10000+mynum)

end subroutine fldout_cyl2
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!================================
end module meshes_io
!================================
