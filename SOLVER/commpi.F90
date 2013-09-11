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

!===============
module commpi
!===============
  
  ! Wrapper routines to invoke the MPI library. 
  ! This routine is the sole place for parallel interactions. 

  use global_parameters
  use data_proc
  use data_mesh,        only : gvec_solid,gvec_fluid
  use data_io,          only : verbose

  ! in case you have problems with the mpi module, you might try to use the
  ! include below, in which case you will have to specify the location in the 
  ! Makefile or copy to the build directory!
#ifndef serial
  use mpi
#endif
  implicit none
  !include 'mpif.h'
  
  public :: ppsum, ppsum_int, ppsum_dble
  public :: ppmin, ppmax, ppmax_int
  public :: ppinit, pbarrier, ppend
  !public :: testing_asynch_messaging_fluid
  !public :: testing_asynch_messaging_solid
  public :: asynch_messaging_fluid
  public :: feed_buffer, send_recv_buffers_solid, extract_from_buffer
  public :: pbroadcast_dble, pbroadcast_char, pbroadcast_log
  public :: pbroadcast_int_arr, pbroadcast_int
  public :: ppcheck, parse_nl
  private

contains

!----------------------------------------------------------------------------
subroutine ppcheck(test, errmsg)
  ! Routine that checks if an  error has occured at all ranks, at some ranks, or
  ! not at all. The message is only printed once if the error occured on all
  ! ranks, otherwise each processor spits its message stdout
  ! newlines in the error message can be achieved using '\n'
  ! IMPORTANT: As this routine contains a barrier, it needs to be allways called
  !            on ALL ranks 

  logical, intent(in)            :: test
  character(len=*), intent(in)   :: errmsg

  integer :: err, errsum, ierror
  
#ifndef serial
  if (test) then 
     err = 1
  else
     err = 0
  endif

  call mpi_allreduce(err, errsum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

  if (errsum == nproc) then
     if (mynum == 0) &
        print '(/,a,/,/,a,/)', 'ERROR on all ranks, error message:', &
                                trim(parse_nl(errmsg))
  elseif (errsum > 0 .and. err > 0) then
     print '(/,a,i4,a,/,/,a,/)', 'local ERROR on rank ', mynum, ', error message:', &
           trim(parse_nl(errmsg))
  endif

  call mpi_barrier(MPI_COMM_WORLD, ierror)

  if (errsum > 0) stop
#endif

end subroutine
!=============================================================================

!-----------------------------------------------------------------------------
function parse_nl(str)
  ! returns the input string with all '\n' in it converted to newlines

  character(len=*), intent(in)  :: str
  character(len=len(str))       :: parse_nl
  integer                       :: lenstr, i

  parse_nl = str

  lenstr = len_trim(parse_nl)

  do i=1, lenstr-1
     if (parse_nl(i:i+1) == '\n') parse_nl(i:i+1) = char(13) // char(10)
  enddo

end function
!=============================================================================

!----------------------------------------------------------------------------
subroutine ppinit
  ! Start message-passing interface, assigning the total number of processors 
  ! nproc and each processor with its local number mynum=0,...,nproc-1.

  use data_comm, only: mpi_realkind
  integer :: ierror
#ifndef serial
  
  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, mynum, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )
  
  if (nproc > 1) then
     if (realkind==4) mpi_realkind = MPI_REAL
     if (realkind==8) mpi_realkind = MPI_DOUBLE_PRECISION
  endif
#endif

end subroutine ppinit
!=============================================================================

!-----------------------------------------------------------------------------
subroutine ppend
  integer :: ierror

#ifndef serial
  call MPI_FINALIZE(ierror)
#endif

end subroutine ppend
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_char(input_char,input_proc)

  integer, intent(in)           :: input_proc
  character(*), intent(inout)   :: input_char
  integer                       :: ierror


#ifndef serial
  call mpi_bcast(input_char, len(input_char), MPI_CHARACTER, input_proc, &
                 MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbroadcast_char
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_log(input_log,input_proc)

  integer, intent(in)    :: input_proc
  logical, intent(inout) :: input_log
  integer                :: ierror

#ifndef serial
  call mpi_bcast(input_log, 1, MPI_LOGICAL, input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbroadcast_log
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_int(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int
  integer                :: ierror

#ifndef serial
  call mpi_bcast(input_int, 1, MPI_INTEGER, input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbroadcast_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_int_arr(input_int, input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int(:)
  integer                :: ierror

#ifndef serial
  call mpi_bcast(input_int, size(input_int), MPI_INTEGER, input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbroadcast_int_arr
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_dble(input_dble,input_proc)

  integer, intent(in)             :: input_proc
  real(kind=dp)   , intent(inout) :: input_dble
  integer                         :: ierror

#ifndef serial
  call mpi_bcast(input_dble, 1, MPI_DOUBLE_PRECISION, input_proc, &
                 MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbroadcast_dble
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function ppmin(scal)

  real(kind=dp)    :: scal
  real(kind=dp)    :: buff, buff2
  integer          :: ierror

  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, MPI_DOUBLE_PRECISION, MPI_MIN,  &
                     MPI_COMM_WORLD, ierror)
#endif

  ppmin = buff2

end function ppmin
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function ppmax(scal)

  real(kind=dp)    :: scal
  real(kind=dp)    :: buff, buff2
  integer          :: ierror

  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, MPI_DOUBLE_PRECISION, MPI_MAX,  &
                     MPI_COMM_WORLD, ierror)
#endif
  ppmax = buff2

end function ppmax
!=============================================================================

!-----------------------------------------------------------------------------
integer function ppmax_int(scal)

  integer :: scal
  integer :: buff, buff2
  integer :: ierror

  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, MPI_INTEGER, MPI_MAX,  &
                     MPI_COMM_WORLD, ierror)
#endif
  ppmax_int = buff2
  
end function ppmax_int
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function ppsum(scal)

  use data_comm, only: mpi_realkind
  real(kind=realkind) :: scal
  real(kind=realkind) :: buff, buff2
  integer             :: ierror

  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, mpi_realkind, MPI_SUM,  &
                     MPI_COMM_WORLD, ierror)
#endif

  ppsum = buff2

end function ppsum
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=dp) function ppsum_dble(scal)

  real(kind=dp)    :: scal
  real(kind=dp)    :: buff, buff2
  integer          :: ierror
  
  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  &
                     MPI_COMM_WORLD, ierror)
#endif

  ppsum_dble = buff2

end function ppsum_dble
!=============================================================================

!-----------------------------------------------------------------------------
integer function ppsum_int(scal)

  integer :: scal
  integer :: buff, buff2  
  integer :: ierror
  
  buff = scal
  buff2 = scal
#ifndef serial
  call MPI_ALLREDUCE(buff, buff2, 1, MPI_INTEGER, MPI_SUM,  &
                     MPI_COMM_WORLD, ierror)
#endif
  
  ppsum_int = buff2

end function ppsum_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbarrier
  integer :: ierror
 
#ifndef serial
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
#endif

end subroutine pbarrier
!=============================================================================

!-----------------------------------------------------------------------------
subroutine feed_buffer(nc)

  use data_comm
  use data_mesh, only: gvec_solid
  use linked_list
  
  integer, intent(in)            :: nc
  integer                        :: ic, imsg, ipg, ip
  integer                        :: sizemsg_solid
  class(link), pointer           :: buffs
  
#ifndef serial
  ! fill send buffer
  call buffs_all%resetcurrent()
  do imsg = 1, sizesend_solid
     buffs => buffs_all%getNext()
     sizemsg_solid = sizemsgsend_solid(imsg)
     do ip = 1, sizemsg_solid
        ipg = glocal_index_msg_send_solid(ip,imsg)
        do ic = 1, nc
           buffs%ldata(ip,ic) = gvec_solid(ipg,ic)
        enddo
     end do
  enddo
#endif

end subroutine feed_buffer
!=============================================================================

!-----------------------------------------------------------------------------
subroutine send_recv_buffers_solid(nc)
  ! Solid asynchronous communication pattern with one message per proc-proc pair.
  ! for a nc-component field gvec. The arrays to map global numbers along 
  ! processor-processor boundaries are determined in the mesher, routine pdb.f90
  ! (consulation thereof to be avoided if at all possible...)
  
  use data_comm
  use linked_list
  
  integer, intent(in) :: nc

#ifndef serial
  integer               :: imsg, sizeb, ipdes, ipsrc
  integer               :: ic, ip, ipg
  integer               :: msgnum, msgnum1, msgnum2
  integer               :: sizemsg_solid
  integer               :: recv_status(MPI_STATUS_SIZE, sizerecv_solid)
  integer               :: send_status(MPI_STATUS_SIZE, sizesend_solid)
  integer               :: ierror
  integer               :: recv_request(1:sizerecv_solid), send_request(1:sizesend_solid)
  class(link), pointer  :: buffs, buffr
  
  ! Send stuff around
  call buffs_all%resetcurrent()
  do imsg = 1, sizesend_solid
     buffs => buffs_all%getnext()
     sizemsg_solid = sizemsgsend_solid(imsg)
     sizeb  = nc * sizemsg_solid
     ipdes  = listsend_solid(imsg)
     msgnum = mynum * nproc + ipdes
     call MPI_ISEND(buffs%ldata, sizeb, mpi_realkind, ipdes, msgnum, &
                    MPI_COMM_WORLD, send_request(imsg), ierror)
  end do

  ! Receive data
  call buffr_all%resetcurrent()
  do imsg = 1, sizerecv_solid
     buffr => buffr_all%getnext()
     sizemsg_solid = sizemsgrecv_solid(imsg)
     sizeb = nc * sizemsg_solid
     ipsrc = listrecv_solid(imsg)
     msgnum1 = ipsrc * nproc + mynum
     call MPI_IRECV(buffr%ldata, sizeb, mpi_realkind, ipsrc, msgnum1, &
                    MPI_COMM_WORLD, recv_request(imsg), ierror)
  enddo
  call MPI_WAITALL(sizerecv_solid, recv_request, recv_status, ierror)
  call MPI_WAITALL(sizesend_solid, send_request, send_status, ierror)
#endif

end subroutine send_recv_buffers_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine extract_from_buffer(vec,nc)

  use data_mesh, only: npol, gvec_solid, igloc_solid
  use data_comm            
  use linked_list
  
  real(kind=realkind), intent(inout) :: vec(0:,0:,:,:)
  integer, intent(in)   :: nc
  integer               :: imsg, ipg, ip, ipol, jpol, iel, ipt, ic
  integer               :: sizemsg_solid
  class(link), pointer  :: buffr
  
#ifndef serial
  do ic = 1, nc
     ! Extract received from buffer
     call buffr_all%resetcurrent()
     do imsg = 1, sizerecv_solid
        buffr => buffr_all%getnext()
        sizemsg_solid = sizemsgrecv_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_recv_solid(ip,imsg)
           gvec_solid(ipg,ic) = gvec_solid(ipg,ic) + buffr%ldata(ip,ic)
        enddo
     enddo
     do ip = 1, num_recv_gll
        ipol = glob2el_recv(ip,1)
        jpol = glob2el_recv(ip,2)
        iel =  glob2el_recv(ip,3)
        ipt = (iel - 1) * (npol + 1)**2 + jpol * (npol + 1) + ipol + 1
        ipg = igloc_solid(ipt)
        vec(ipol,jpol,iel,ic) = gvec_solid(ipg,ic)
     enddo

  end do
#endif

end subroutine extract_from_buffer
!=============================================================================

!-----------------------------------------------------------------------------
subroutine asynch_messaging_fluid
  ! Fluid asynchronous communication pattern with one message per proc-proc pair
  ! for a single-component field gvec. The arrays to map global numbers along 
  ! processor-processor boundaries are determined in the mesher, routine pdb.f90
  ! (consulation thereof to be avoided if at all possible...)

  use data_comm
  use linked_list
  
#ifndef serial
  integer               :: imsg, ipg, ip, sizeb, ipdes, ipsrc
  integer               :: msgnum, msgnum1, msgnum2
  integer               :: sizemsg_fluid
  integer               :: ierror
  integer               :: recv_status(MPI_STATUS_SIZE, sizerecv_fluid)
  integer               :: send_status(MPI_STATUS_SIZE, sizesend_fluid)
  integer               :: recv_request(1:sizerecv_fluid), send_request(1:sizesend_fluid)
  class(link), pointer  :: buffs, buffr
  
  ! Prepare arrays to be sent.... MIGHT USE A POWER OF 2 STATEMENT THERE
  call buffs_all_fluid%resetcurrent()
  do imsg = 1, sizesend_fluid
     buffs => buffs_all_fluid%getnext()
     sizemsg_fluid = sizemsgsend_fluid(imsg)
     do ip = 1, sizemsg_fluid
        ipg = glocal_index_msg_send_fluid(ip,imsg)
        buffs%ldata(ip,1) = gvec_fluid(ipg)
     end do
     
     ! Send stuff around
     sizeb  = sizemsg_fluid
     ipdes  = listsend_fluid(imsg)
     msgnum = mynum*nproc + ipdes
     call MPI_ISEND(buffs%ldata, sizeb, mpi_realkind, ipdes, msgnum, &
                   MPI_COMM_WORLD, send_request(imsg), ierror)
  end do

  ! Receive data, sum things up
  call buffr_all_fluid%resetcurrent()
  do imsg = 1, sizerecv_fluid
     buffr => buffr_all_fluid%getnext()
     sizemsg_fluid = sizemsgrecv_fluid(imsg)
     sizeb = sizemsg_fluid
     ipsrc = listrecv_fluid(imsg)
     msgnum1 = ipsrc*nproc + mynum
     call MPI_IRECV(buffr%ldata, sizeb, mpi_realkind, ipsrc, msgnum1, &
                   MPI_COMM_WORLD, recv_request(imsg), ierror)
  end do

  call MPI_WAITALL(sizerecv_fluid, recv_request, recv_status, ierror)
  call MPI_WAITALL(sizesend_fluid, send_request, send_status, ierror)

  call buffr_all_fluid%resetcurrent()
  do imsg = 1, sizerecv_fluid
     buffr => buffr_all_fluid%getnext()
     sizemsg_fluid = sizemsgrecv_fluid(imsg)
     do ip = 1, sizemsg_fluid
        ipg = glocal_index_msg_recv_fluid(ip,imsg)
        ! add received buffer to own field at same global point
        gvec_fluid(ipg) = gvec_fluid(ipg) + buffr%ldata(ip,1)
     end do
  end do
#endif

end subroutine asynch_messaging_fluid
!=============================================================================

!===================
end module commpi
!===================
