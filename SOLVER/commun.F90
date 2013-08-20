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

!> This is the communication module which loads/sorts data
!! to exchange/examine over the processors.
!!
!! APRIL 2007:
!! At this level, we only call parallel routines, but do not invoke MPI here. 
module commun
  
  use global_parameters
  use commpi
  use data_proc
  
  implicit none
  public :: comm2d ! the general assembly & communication routine
  public :: gather_elem_solid, scatter_elem_solid
  public :: gather_elem_fluid, scatter_elem_fluid
  public :: pdistsum_solid, pdistsum_fluid
  
  public :: assembmass_sum_solid, assembmass_sum_fluid ! assemble and sum massmat
  public :: broadcast_int, broadcast_int_arr
  public :: broadcast_dble, broadcast_char, broadcast_log
  public :: pinit, pend
  public :: pmin, pmax, pmax_int, psum, psum_int, psum_dble
  public :: barrier, pcheck
  public :: mpi_asynch_messaging_test_solid, mpi_asynch_messaging_test_fluid
  private

contains

!-----------------------------------------------------------------------------
!> This is a driver routine to call the assembly of field f of dimension nc
!! and either solid or fluid subdomains.
!! The global assembly is discarded as it is not necessary during the time loop
!! and therefore chose not to store any global numbering arrays.
!! If nproc>1, then internode message passing is applied where necessary.
subroutine comm2d(f, nel, nc, domainin)
  use data_mesh, only: npol 
  character(len=5), intent(in)       :: domainin
  integer, intent(in)                :: nc,nel
  real(kind=realkind), intent(inout) :: f(0:npol,0:npol,1:nel,1:nc)
  
    if (domainin=='total') then
       if (lpr) &
          write(6,'(a/a)') 'PROBLEM: Discarded this case since igloc is not',&
                           '         known in the solver any longer...'
       stop
    elseif (domainin=='solid') then
       call pdistsum_solid(f,nc)
    elseif (domainin=='fluid') then
       call pdistsum_fluid(f)
    else
       if (lpr) write(6,*) 'Assembly: Domain', domainin, ' non-existent!' 
       stop
    end if
  
end subroutine comm2d
!=============================================================================

!-----------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the solid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge 
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
!! If nproc>1, the asynchronous messaging scheme is invoked to additionally 
!! sum & exchange values on processor boundary points.
subroutine pdistsum_solid(vec, nc)
  
  use data_mesh,   only: igloc_solid
  use data_mesh,        only: gvec_solid, npol, nel_solid
  use data_time,        only: idmpi, iclockmpi
  use clocks_mod
  
  !include 'mesh_params.h' 
  
  integer, intent(in)                :: nc
  real(kind=realkind), intent(inout) :: vec(0:,0:,:,:)
  integer                            :: ic, iel, jpol, ipol, idest, ipt
  
  do ic = 1, nc
     ! Gather element boundaries
     gvec_solid(:) = 0.d0
     ipt = 1

     do iel = 1, nel_solid
         
        jpol = 0
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           gvec_solid(idest) = gvec_solid(idest) + vec(ipol,jpol,iel,ic)
           ipt = ipt + 1
        end do

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_solid(ipt)
           gvec_solid(idest) = gvec_solid(idest) + vec(ipol,jpol,iel,ic)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_solid(ipt)
           gvec_solid(idest) = gvec_solid(idest) + vec(ipol,jpol,iel,ic)
           ipt = ipt + 1
        end do

        jpol = npol
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           gvec_solid(idest) = gvec_solid(idest) + vec(ipol,jpol,iel,ic)
           ipt = ipt + 1
        end do

     end do

     ! Collect processor boundaries into buffer for each component
     iclockmpi = tick()
#ifndef serial
     if (nproc>1) call feed_buffer(ic)
#endif
     iclockmpi = tick(id=idmpi,since=iclockmpi)
  
     ! Scatter
     ipt = 1
     do iel = 1, nel_solid
        
        jpol = 0
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,ic) = gvec_solid(idest)
           ipt = ipt + 1
        end do

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,ic) = gvec_solid(idest)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,ic) = gvec_solid(idest)
           ipt = ipt + 1
        end do

        jpol = npol
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,ic) = gvec_solid(idest)
           ipt = ipt + 1
        end do

     end do
  end do
  
#ifndef serial
  iclockmpi = tick()
  if (nproc>1) then
     ! Do message-passing for all components at once
     call send_recv_buffers_solid(nc)
     ! Extract back into each component sequentially
     call extract_from_buffer(vec,nc)
  endif ! nproc>1
  iclockmpi = tick(id=idmpi,since=iclockmpi)
#endif
  
end subroutine pdistsum_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine gather_elem_solid(vec, ic, iel)
  
  use data_mesh,   only: igloc_solid
  use data_mesh,        only: gvec_solid, npol
  use data_time,        only: idmpi, iclockmpi
  use clocks_mod
  
  !include 'mesh_params.h' 
  
  integer, intent(in)                :: ic, iel
  real(kind=realkind), intent(inout) :: vec(0:,0:,:,:)
  integer                            :: jpol, ipol, idest, ipt
  
  do jpol = 0, npol
     do ipol = 0, npol
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        idest = igloc_solid(ipt)
        gvec_solid(idest) = gvec_solid(idest)+vec(ipol,jpol,iel,ic)
     end do
  end do
  
end subroutine gather_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine scatter_elem_solid(vec, ic, iel)
  
  use data_mesh,   only: igloc_solid
  use data_mesh,        only: gvec_solid, npol
  use data_time,        only: idmpi, iclockmpi
  use clocks_mod
  
  !include 'mesh_params.h' 
  
  integer, intent(in)                :: ic, iel
  real(kind=realkind), intent(inout) :: vec(0:,0:,:,:)
  integer                            :: jpol, ipol, idest, ipt
  
  do jpol = 0, npol
     do ipol = 0, npol
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        idest = igloc_solid(ipt)
        vec(ipol,jpol,iel,ic) = gvec_solid(idest)
     end do
  end do
  
end subroutine scatter_elem_solid
!=============================================================================

!-----------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the fluid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge 
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
subroutine pdistsum_fluid(vec)
  
  use data_mesh,   only: igloc_fluid
  use data_time,        only: idmpi, iclockmpi
  use data_mesh,        only: gvec_fluid, npol, nel_fluid
  use clocks_mod
  
  !include 'mesh_params.h' 
  
  real(kind=realkind), intent(inout) :: vec(0:npol,0:npol,nel_fluid)
  integer                            :: iel, jpol, ipol, idest, ipt

  gvec_fluid(:) = 0.d0

  ! Gather
  ipt = 1
  do iel = 1, nel_fluid

     jpol = 0
     do ipol = 0, npol
        idest = igloc_fluid(ipt)
        gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
        ipt = ipt + 1
     end do

     do jpol = 1, npol-1
        ipol = 0
        idest = igloc_fluid(ipt)
        gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
        ipt = ipt + npol

        ipol = npol
        idest = igloc_fluid(ipt)
        gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
        ipt = ipt + 1
     end do

     jpol = npol
     do ipol = 0, npol
        idest = igloc_fluid(ipt)
        gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
        ipt = ipt + 1
     end do

  end do

#ifndef serial
  iclockmpi = tick()
  if (nproc>1) call asynch_messaging_fluid
  iclockmpi = tick(id=idmpi, since=iclockmpi)
#endif

  ! Scatter
  ipt = 1
  do iel = 1, nel_fluid

     jpol = 0
     do ipol = 0, npol
        idest = igloc_fluid(ipt)
        vec(ipol,jpol,iel) = gvec_fluid(idest)
        ipt = ipt + 1
     end do

     do jpol = 1, npol-1
        ipol = 0
        idest = igloc_fluid(ipt)
        vec(ipol,jpol,iel) = gvec_fluid(idest)
        ipt = ipt + npol

        ipol = npol
        idest = igloc_fluid(ipt)
        vec(ipol,jpol,iel) = gvec_fluid(idest)
        ipt = ipt + 1
     end do

     jpol = npol
     do ipol = 0, npol
        idest = igloc_fluid(ipt)
        vec(ipol,jpol,iel) = gvec_fluid(idest)
        ipt = ipt + 1
     end do

  end do

end subroutine pdistsum_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine gather_elem_fluid(vec, iel)

  use data_mesh,   only: igloc_fluid
  use data_mesh,        only: gvec_fluid, npol
  !include 'mesh_params.h' 
  
  real(kind=realkind), intent(in)    :: vec(0:,0:,:)
  integer, intent(in)                :: iel
  integer                            :: jpol, ipol, idest, ipt

  do jpol = 0, npol
     do ipol = 0, npol
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        idest = igloc_fluid(ipt)
        !$omp atomic
        gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel) 
     end do
  end do
  
end subroutine gather_elem_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine scatter_elem_fluid(vec, iel)
     
  use data_mesh,   only: igloc_fluid
  use data_mesh,        only: gvec_fluid, npol, nel_fluid
  !include 'mesh_params.h' 
  
  real(kind=realkind), intent(out)   :: vec(0:npol,0:npol,nel_fluid)
  integer, intent(in)                :: iel
  integer                            :: jpol, ipol, idest, ipt

  do jpol = 0, npol
     do ipol = 0, npol
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        idest = igloc_fluid(ipt)
        vec(ipol,jpol,iel) = gvec_fluid(idest)
     end do
  end do
  
end subroutine scatter_elem_fluid
!=============================================================================

!-----------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the solid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge 
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
!! If nproc>1, the asynchronous messaging scheme is invoked to additionally 
!! sum & exchange values on processor boundary points.
!!
!! The local arrays are allocatable since this routine is only called before 
!! the time loop.
subroutine mpi_asynch_messaging_test_solid
  
  use data_mesh, only: igloc_solid, nglob_solid
  use data_mesh,      only: npol, nel_solid, nel_fluid
  !include 'mesh_params.h' 
  
  real(kind=realkind),allocatable   :: vec(:,:,:,:)
  real(kind=realkind),allocatable   :: gvec_solid2(:,:)
  integer                           :: ic, iel, jpol, ipol, idest, ipt

  allocate(vec(0:npol,0:npol,nel_solid,3))
  allocate(gvec_solid2(nglob_solid,3))

  gvec_solid2(:,:) = 0.d0

  do ic = 1, 3

     vec(:,:,:,ic)=real(ic,kind=realkind)

     ! Gather
     do iel = 1, nel_solid
        do jpol = 0, npol
           do ipol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_solid(ipt)
              gvec_solid2(idest,ic) = vec(ipol,jpol,iel,ic)
           end do
        end do
     end do
  end do

#ifndef serial
  if (nproc > 1) &
        call testing_asynch_messaging_solid(gvec_solid2,3)
#endif

  do ic = 1, 3
     ! Scatter
     do iel = 1, nel_solid
        do jpol = 0, npol
           do ipol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_solid(ipt)
              vec(ipol,jpol,iel,ic) = gvec_solid2(idest,ic)
           end do
        end do
     end do
  end do

  deallocate(gvec_solid2)
  deallocate(vec)

end subroutine mpi_asynch_messaging_test_solid
!=============================================================================

!-----------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the fluid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge 
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
!! If nproc>1, the asynchronous messaging scheme is invoked to additionally 
!! sum & exchange values on processor boundary points.
!!
!! The local arrays are allocatable since this routine is only called before 
!! the time loop.
subroutine mpi_asynch_messaging_test_fluid
  
  use data_mesh,   only: igloc_fluid
  use data_mesh,        only: gvec_fluid, npol, nel_fluid
  
  !include 'mesh_params.h' 
  
  real(kind=realkind),allocatable :: vec(:,:,:)
  integer                         :: iel,jpol,ipol,idest,ipt

  allocate(vec(0:npol,0:npol,nel_fluid))

  gvec_fluid(:) = 0.d0
  vec(:,:,:) = 1.d0

  ! Gather
  do iel = 1, nel_fluid
     do jpol = 0, npol
        do ipol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = vec(ipol,jpol,iel)
        end do
     end do
  end do

#ifndef serial
  if (nproc>1) call testing_asynch_messaging_fluid
#endif

  ! Scatter
  do iel = 1, nel_fluid
     do jpol = 0, npol
        do ipol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
           idest = igloc_fluid(ipt)
           vec(ipol,jpol,iel) = gvec_fluid(idest)
        end do
     end do
  end do

  deallocate(vec)

end subroutine mpi_asynch_messaging_test_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine assembmass_sum_solid(f1,res)

  use data_mesh,   only: igloc_solid
  use data_mesh,        only: gvec_solid, npol, nel_solid
  !include 'mesh_params.h'
  
  real(kind=realkind), intent(in)   :: f1(0:,0:,:)
  real(kind=dp)   , intent(out)     :: res
  integer                           :: ipt, idest, iel, ipol, jpol

  !!@TODO Optimise for npol = 4
  res = 0.d0 
  gvec_solid(:) = 0.d0
  do iel = 1, nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
           idest = igloc_solid(ipt)
           gvec_solid(idest) = gvec_solid(idest) + f1(ipol,jpol,iel)
        end do
     end do
  end do
  res = res + sum(gvec_solid(:))
#ifndef serial
  if (nproc>1) res=ppsum_dble(res)
#endif

end subroutine assembmass_sum_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine assembmass_sum_fluid(f1,res)

  use data_mesh,   only: igloc_fluid
  use data_mesh,        only: gvec_fluid
  use data_mesh,        only: gvec_solid, npol, nel_fluid
  !include 'mesh_params.h' 
  
  real(kind=realkind), intent(in)   :: f1(0:,0:,:)
  real(kind=dp)   , intent(out)     :: res
  integer ipt, idest
  integer iel, ipol, jpol

  res = 0.d0
  
  gvec_fluid(:) = 0.d0
  do iel = 1, nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + f1(ipol,jpol,iel)
        end do
     end do
  end do
  res = res + sum(gvec_fluid)

#ifndef serial
  if (nproc>1) res = ppsum_dble(res)
#endif

end subroutine assembmass_sum_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pinit

  !include 'mesh_params.h'
  use data_mesh, only: nproc_mesh
  
#ifndef serial
  ! Start message passing interface if parallel simulation
  if (nproc_mesh > 1) then 
     call ppinit
     if (nproc_mesh /= nproc) then        
        write(6,*) mynum, 'Problem with number of processors!'
        write(6,*) mynum, 'Mesh constructed for:', nproc_mesh
        write(6,*) mynum, 'Job submission for:', nproc
        stop
     endif
  else
     nproc = nproc_mesh
     mynum = 0
  endif
#endif

#ifdef serial
  if (nproc_mesh /= 1) &
        write(6,*) 'ERROR: Solver compiled with SERIAL flag, but mesh has nproc > 1: ', &
                    nproc_mesh
  nproc = 1
  mynum = 0
#endif

  lpr = .false.
  if (nproc>1) then
     if (mynum==nproc/2-1) lpr = .true.
  else
     lpr = .true.
  endif

  call define_io_appendix(appmynum, mynum)
 
  procstrg = 'Proc '//appmynum(3:4)//' '

  if (lpr) write(6,'(a,i5)') '    Initialized run for nproc =', nproc

end subroutine pinit
!=============================================================================

!-----------------------------------------------------------------------------
!! End message passing interface if parallel
subroutine pend

#ifndef serial
  if (nproc>1) call ppend
#endif

end subroutine pend
!=============================================================================

!-----------------------------------------------------------------------------
subroutine broadcast_char(input_char,input_proc)

  character(*), intent(inout)   :: input_char
  integer, intent(in)           :: input_proc

#ifndef serial
  if (nproc>1) call pbroadcast_char(input_char,input_proc)
#endif
   
end subroutine broadcast_char
!=============================================================================

!-----------------------------------------------------------------------------
subroutine broadcast_log(input_log,input_proc)

  integer, intent(in)    :: input_proc
  logical, intent(inout) :: input_log

#ifndef serial
  if (nproc>1) call pbroadcast_log(input_log,input_proc)
#endif
   
end subroutine broadcast_log
!=============================================================================

!-----------------------------------------------------------------------------
subroutine broadcast_int(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int

#ifndef serial
  if (nproc>1) call pbroadcast_int(input_int, input_proc)
#endif
   
end subroutine broadcast_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine broadcast_int_arr(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int(:)

#ifndef serial
  if (nproc>1) call pbroadcast_int_arr(input_int, input_proc)
#endif

end subroutine broadcast_int_arr
!=============================================================================

!-----------------------------------------------------------------------------
subroutine broadcast_dble(input_dble,input_proc)

  integer, intent(in)             :: input_proc
  real(kind=dp)   , intent(inout) :: input_dble

#ifndef serial
  if (nproc>1) call pbroadcast_dble(input_dble,input_proc)
#endif

end subroutine broadcast_dble
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function pmin(scal)

  real(kind=dp)    :: scal
  
  pmin = scal
#ifndef serial
  if (nproc>1) pmin = ppmin(scal)
#endif

end function pmin
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function pmax(scal)

  real(kind=dp)    :: scal

  pmax = scal
#ifndef serial
  if (nproc>1) pmax = ppmax(scal)
#endif

end function pmax
!=============================================================================

!-----------------------------------------------------------------------------
integer function pmax_int(scal)

  integer :: scal

  pmax_int=scal
#ifndef serial
  if (nproc>1) pmax_int = ppmax_int(scal)
#endif

end function pmax_int
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function psum(scal)

  real(kind=realkind) :: scal

  psum = scal
#ifndef serial
  if (nproc>1) psum = ppsum(scal)
#endif

end function psum
!=============================================================================

!-----------------------------------------------------------------------------
integer function psum_int(scal)

  integer :: scal

  psum_int = scal
#ifndef serial
  if (nproc>1) psum_int = ppsum_int(scal)
#endif

end function psum_int
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function psum_dble(scal)

  real(kind=dp)    :: scal

  psum_dble = scal
#ifndef serial
  if (nproc>1) psum_dble = ppsum_dble(scal)
#endif

end function psum_dble
!=============================================================================

!-----------------------------------------------------------------------------
subroutine barrier
 
#ifndef serial
  if (nproc>1) call pbarrier
#endif

end subroutine barrier
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pcheck(test, errmsg)
 
  logical, intent(in)            :: test
  character(len=*), intent(in)   :: errmsg
  
#ifndef serial
  if (nproc > 1) call ppcheck(test, errmsg)
#endif
  if (nproc == 1 .and. test) then
     print '(/,a,/,/,a,/)', 'ERROR:', trim(parse_nl(errmsg))
     stop
  endif

end subroutine pcheck
!=============================================================================

!====================
end module commun
!====================
