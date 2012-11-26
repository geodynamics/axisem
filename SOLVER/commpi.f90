!===============
module commpi
!===============
!
! Wrapper routines to invoke the MPI library. 
! This routine is the sole place for parallel interactions. 
! In other words, the routines are called by wrappers that contain an 
! if (nproc>1) statement such that a serial version 
! (i.e. without MPI libraries) of this code *shall* run after merely taking 
! out this module commpi. Consequently, this module is 'public' 
! (every routine here needs to be called by external wrapper)
!
! WARNING: Need to make sure this is consistent for other issues like 
! source and receiver locations in global/local reference.
!
use global_parameters
use data_mesh, ONLY: gvec_solid,gvec_fluid
use data_proc
use data_comm

! in case you have problems with the mpi module, you might try to use the
! include below, in which case you will have to specify the location in the 
! Makefile or copy to the build directory!
use mpi
implicit none
!include 'mpif.h'

public 

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------------
subroutine ppinit
!
! Start message-passing interface, assigning the total number of processors 
! nproc and each processor with its local number mynum=0,...,nproc-1.
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, mynum, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )

end subroutine ppinit
!=============================================================================

!-----------------------------------------------------------------------------
subroutine ppend

  call MPI_FINALIZE(ierror)

end subroutine ppend
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_int(input_int,input_proc)

integer, intent(in)    :: input_proc
integer, intent(inout) :: input_int

  call mpi_bcast(input_int,1,mpi_integer,input_proc, mpi_comm_world,ierror)

end subroutine pbroadcast_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_dble(input_dble,input_proc)

integer, intent(in)             :: input_proc
double precision, intent(inout) :: input_dble

!write(69,*)'pbroadcast_dble1:',input_dble,input_proc;call flush(69)

  call mpi_bcast(input_dble,1,mpi_double_precision,input_proc,&
                 mpi_comm_world, ierror)
!write(69,*)'pbroadcast_dble2:',input_dble;call flush(69)

end subroutine pbroadcast_dble
!=============================================================================

!-----------------------------------------------------------------------------
double precision function ppmin(scal)

double precision :: scal
double precision :: buff,buff2

  buff = scal
  buff2 = scal
  call MPI_ALLREDUCE(buff,buff2,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                     MPI_COMM_WORLD,IERROR)

  ppmin = buff2

end function ppmin
!=============================================================================

!-----------------------------------------------------------------------------
double precision function ppmax(scal)

double precision :: scal
double precision :: buff,buff2

  buff = scal
  buff2 = scal
  write(69,*)'ppmax1:',scal
  call MPI_ALLREDUCE(buff,buff2,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                     MPI_COMM_WORLD,IERROR)
  ppmax=buff2
  write(69,*)'ppmax2:',buff2

end function ppmax
!=============================================================================

!-----------------------------------------------------------------------------
integer function ppmax_int(scal)

integer :: scal
integer :: buff,buff2

  buff = scal
  buff2 = scal
  call MPI_ALLREDUCE(buff,buff2,1,MPI_INTEGER,MPI_MAX, &
                     MPI_COMM_WORLD,IERROR)
  ppmax_int=buff2
  
end function ppmax_int
!=============================================================================

!-----------------------------------------------------------------------------
real(kind=realkind) function ppsum(scal)

real(kind=realkind) :: scal
real(kind=realkind) :: buff,buff2

  buff = scal
  buff2 = scal
 call MPI_ALLREDUCE(buff,buff2,1,mpi_realkind,MPI_SUM, &
                    MPI_COMM_WORLD,IERROR)

  ppsum = buff2

end function ppsum
!=============================================================================

!-----------------------------------------------------------------------------
double precision function ppsum_dble(scal)

double precision :: scal
double precision :: buff,buff2

  buff = scal
  buff2 = scal
  call MPI_ALLREDUCE(buff,buff2,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                     MPI_COMM_WORLD,IERROR)

  ppsum_dble = buff2

end function ppsum_dble
!=============================================================================

!-----------------------------------------------------------------------------
integer function ppsum_int(scal)

integer :: scal
integer :: buff,buff2
  
  buff = scal
  buff2 = scal
  call MPI_ALLREDUCE(buff,buff2,1,MPI_INTEGER,MPI_SUM, &
                     MPI_COMM_WORLD,IERROR)
  
  ppsum_int = buff2

end function ppsum_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbarrier
 
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine pbarrier
!=============================================================================

!-----------------------------------------------------------------------------
subroutine feed_buffer(ic)
!
!
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_comm, ONLY : buffs_all,buffr_all
use data_mesh, ONLY: gvec_solid

include 'mesh_params.h'
integer, intent(in)            :: ic
integer                        :: imsg,ipg,ip
integer                        :: sizemsg_solid

! Fill send buffer
  if (sizesend_solid>0) then
     buffs_all(:,ic,:) = 0.d0
     do imsg = 1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_send_solid(ip,imsg)
           buffs_all(ip,ic,imsg) = gvec_solid(ipg)
        end do
     enddo
  endif
  
! Fill receive buffer
  if (sizerecv_solid>0) then 
     buffr_all(:,ic,:) = 0.d0
     do imsg = 1, sizerecv_solid
        sizemsg_solid = sizemsgrecv_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_recv_solid(ip,imsg)
           buffr_all(ip,ic,imsg) = gvec_solid(ipg)
        end do
     end do
  endif

end subroutine feed_buffer
!=============================================================================

!-----------------------------------------------------------------------------
subroutine send_recv_buffers_solid(nc)
!
! Solid asynchronous communication pattern with one message per proc-proc pair.
! for a nc-component field gvec. The arrays to map global numbers along 
! processor-processor boundaries are determined in the mesher, routine pdb.f90
! (consulation thereof to be avoided if at all possible...)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
use data_comm, ONLY : buffr_solid,buffs_solid,buffs_all,buffr_all

include 'mesh_params.h'

integer, intent(in) :: nc
integer             :: imsg,sizeb,ipdes,ipsrc
integer             :: msgnum,msgnum1,msgnum2
integer             :: status(MPI_STATUS_SIZE),sizemsg_solid

! Send stuff around (phase 1)
  if (sizesend_solid>0) then
     do imsg = 1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        buffs_solid(1:sizemsg_solid,1:nc)=buffs_all(1:sizemsg_solid,1:nc,imsg)
        sizeb  = nc*sizemsg_solid
        ipdes  = listsend_solid(imsg)
        msgnum = mynum*nproc + ipdes
        call MPI_SEND(buffs_solid,sizeb,mpi_realkind,&
                      ipdes,msgnum,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive data, sum things up and send back to initial sender (phase 2)
  if (sizerecv_solid>0) then 
     do imsg = 1, sizerecv_solid
        sizemsg_solid = sizemsgrecv_solid(imsg)
        sizeb = nc*sizemsg_solid
        ipsrc = listrecv_solid(imsg)
        msgnum1 = ipsrc*nproc + mynum
        call MPI_RECV(buffr_solid,sizeb,mpi_realkind,&
                      ipsrc,msgnum1,MPI_COMM_WORLD,status,IERROR)

!       add received buffer to own field at same global point
        buffr_all(1:sizemsg_solid,1:nc,imsg) = &
                                  buffr_all(1:sizemsg_solid,1:nc,imsg) + &
                                  buffr_solid(1:sizemsg_solid,1:nc)

!       assuming that each global point is mapped one2one: each ip has one ipg
        buffr_solid(1:sizemsg_solid,1:nc)= buffr_all(1:sizemsg_solid,1:nc,imsg)

!       send joint data back, but stick into new envelope/msgnum
        msgnum2 = mynum*nproc + ipsrc
        call MPI_SEND(buffr_solid,sizeb,mpi_realkind,&
                      ipsrc,msgnum2,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive updated data back (phase 3)
  if (sizesend_solid>0) then
     do imsg =1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        sizeb = nc*sizemsg_solid
        ipsrc = listsend_solid(imsg)
        msgnum = ipsrc*nproc + mynum
        call MPI_RECV(buffs_solid,sizeb,mpi_realkind,&
                      ipsrc,msgnum,MPI_COMM_WORLD,status,IERROR)
        buffs_all(1:sizemsg_solid,1:nc,imsg)=buffs_solid(1:sizemsg_solid,1:nc)
     enddo
  endif

end subroutine send_recv_buffers_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine extract_from_buffer(vec,nc)
!
!
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_comm, ONLY : buffs_all,buffr_all,num_send_gll,num_recv_gll
use data_comm, ONLY : glob2el_send,glob2el_recv
use data_mesh, ONLY: gvec_solid
use data_numbering, ONLY: igloc_solid

include 'mesh_params.h'

integer, intent(in) :: nc
real(kind=realkind), intent(inout) :: vec(0:npol,0:npol,nel_solid,nc)
integer             :: imsg,ipg,ip,ipol,jpol,iel,ipt,ic
integer             :: sizemsg_solid


!----------
  do ic = 1, nc
!----------

! Extract back-received from buffer (phase 3)
  if (sizesend_solid>0) then
     do imsg =1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_send_solid(ip,imsg)
           gvec_solid(ipg) = buffs_all(ip,ic,imsg)
        enddo
     enddo
     do ip = 1, num_send_gll
        ipol = glob2el_send(ip,1)
        jpol = glob2el_send(ip,2)
        iel =  glob2el_send(ip,3)
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        ipg = igloc_solid(ipt)
        vec(ipol,jpol,iel,ic) = gvec_solid(ipg)
     enddo
  endif

! Extract received from buffer (phase2)
  if (sizerecv_solid>0) then
     do imsg =1, sizerecv_solid
        sizemsg_solid = sizemsgrecv_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_recv_solid(ip,imsg)
           gvec_solid(ipg) = buffr_all(ip,ic,imsg)
        enddo
     enddo
     do ip = 1, num_recv_gll
        ipol = glob2el_recv(ip,1)
        jpol = glob2el_recv(ip,2)
        iel =  glob2el_recv(ip,3)
        ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
        ipg = igloc_solid(ipt)
        vec(ipol,jpol,iel,ic) = gvec_solid(ipg)
     enddo
  endif

!----------
     end do
!----------

end subroutine extract_from_buffer
!=============================================================================

!-----------------------------------------------------------------------------
subroutine testing_asynch_messaging_solid(gvec_solid2,nc)
!
! Solid asynchronous communication pattern with one message per proc-proc pair.
! for a nc-component field gvec. The arrays to map global numbers along 
! processor-processor boundaries are determined in the mesher, routine pdb.f90
! (consulation thereof to be avoided if at all possible...)
! Same as above but with output, used solely for the mpi test preloop.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_comm, ONLY : buffr_solid,buffs_solid,sizesend_solid,sizerecv_solid

include 'mesh_params.h'
real(kind=realkind),intent(inout)  :: gvec_solid2(:,:)
integer, intent(in)                :: nc
integer                            :: imsg,ipg,ip,sizeb,ipdes,ipsrc
integer                            :: msgnum,msgnum1,msgnum2
integer                            :: status(MPI_STATUS_SIZE),sizemsg_solid

! Prepare arrays to be sent.... MIGHT USE A POWER OF 2 STATEMENT THERE
  write(69,*)' Asynchrounous solid communication test:'
  call flush(69)

  if (sizesend_solid>1) then
     write(6,*)
     write(6,*)'  PROBLEM:',procstrg,'sending more than one message!',&
               sizesend_solid
     stop
  endif
  if (sizerecv_solid>1) then 
     write(6,*)
     write(6,*)'  PROBLEM:',procstrg,'receiving more than one message!',&
               sizerecv_solid
     stop
  endif

  if (sizesend_solid>0) then
     do imsg = 1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_send_solid(ip,imsg)
           buffs_solid(ip,1:nc) = gvec_solid2(ipg,1:nc)
        end do
        
! Send stuff around
        sizeb  = nc*sizemsg_solid
        ipdes  = listsend_solid(imsg)
        msgnum = mynum*nproc + ipdes
        write(69,12)procstrg,'SENDING #',msgnum,sizemsg_solid,' to',ipdes
        call MPI_SEND(buffs_solid,sizeb,mpi_realkind,&
                      ipdes,msgnum,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive data, sum things up and send back to initial sender
  if (sizerecv_solid>0) then 
     do imsg = 1, sizerecv_solid
        sizemsg_solid = sizemsgrecv_solid(imsg)
        sizeb = nc*sizemsg_solid
        ipsrc = listrecv_solid(imsg)
        msgnum1 = ipsrc*nproc + mynum
        write(69,12)procstrg,'RECEIVING #',msgnum1,sizemsg_solid,' from',ipsrc
        call MPI_RECV(buffr_solid,sizeb,mpi_realkind,&
                      ipsrc,msgnum1,MPI_COMM_WORLD,status,IERROR)
!     Add received buffer to own field at same global point
           do ip = 1, sizemsg_solid
              ipg = glocal_index_msg_recv_solid(ip,imsg)
!           assuming here that each global point is mapped one2one from 
!           the buffer... i.e. each ip has one ipg
               gvec_solid2(ipg,1:nc) = gvec_solid2(ipg,1:nc)+buffr_solid(ip,1:nc)
               buffr_solid(ip,1:nc) = gvec_solid2(ipg,1:nc)
           end do
!     Send joint data back, but stick into new envelope/msgnum
        msgnum2 = mynum*nproc + ipsrc
        write(69,12)procstrg,'RETURNING #',msgnum2,sizemsg_solid,' to',ipsrc
        call MPI_SEND(buffr_solid,sizeb,mpi_realkind,&
                      ipsrc,msgnum2,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive updated data back
  if (sizesend_solid>0) then
     do imsg =1, sizesend_solid
        sizemsg_solid = sizemsgsend_solid(imsg)
        sizeb = nc*sizemsg_solid
        ipsrc = listsend_solid(imsg)
        msgnum = ipsrc*nproc + mynum
        write(69,12)procstrg,'RECV UPDATED #',msgnum,&
                    sizemsg_solid,' from',ipsrc
        call MPI_RECV(buffs_solid,sizeb,mpi_realkind,&
             ipsrc,msgnum,MPI_COMM_WORLD,status,IERROR)
        do ip = 1, sizemsg_solid
           ipg = glocal_index_msg_send_solid(ip,imsg)
           gvec_solid2(ipg,1:nc) = buffs_solid(ip,1:nc)
        enddo
     enddo
  endif

  write(69,13)procstrg; call flush(69)

12 format('   MPI: ',a8,a15,i3,', size=',i8,a6,i3)
13 format('   MPI: ',a8,' DONE.')

end subroutine testing_asynch_messaging_solid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine asynch_messaging_fluid
!
! Fluid asynchronous communication pattern with one message per proc-proc pair
! for a single-component field gvec. The arrays to map global numbers along 
! processor-processor boundaries are determined in the mesher, routine pdb.f90
! (consulation thereof to be avoided if at all possible...)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_comm, ONLY : buffr_fluid,buffs_fluid

include 'mesh_params.h'

integer                            :: imsg,ipg,ip,sizeb,ipdes,ipsrc
integer                            :: msgnum,msgnum1,msgnum2
integer                            :: status(MPI_STATUS_SIZE),sizemsg_fluid

! Prepare arrays to be sent.... MIGHT USE A POWER OF 2 STATEMENT THERE
  if (sizesend_fluid>0) then
     do imsg = 1, sizesend_fluid
        sizemsg_fluid = sizemsgsend_fluid(imsg)
        do ip = 1, sizemsg_fluid
           ipg = glocal_index_msg_send_fluid(ip,imsg)
           buffs_fluid(ip) = gvec_fluid(ipg)
        end do
        
! Send stuff around
        sizeb  = sizemsg_fluid
        ipdes  = listsend_fluid(imsg)
        msgnum = mynum*nproc + ipdes
        call MPI_SEND(buffs_fluid,sizeb,mpi_realkind,&
                      ipdes,msgnum,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive data, sum things up and send back to initial sender
  if (sizerecv_fluid>0) then 
     do imsg = 1, sizerecv_fluid
        sizemsg_fluid = sizemsgrecv_fluid(imsg)
        sizeb = sizemsg_fluid
        ipsrc = listrecv_fluid(imsg)
        msgnum1 = ipsrc*nproc + mynum
        call MPI_RECV(buffr_fluid,sizeb,mpi_realkind,&
                      ipsrc,msgnum1,MPI_COMM_WORLD,status,IERROR)
        do ip = 1, sizemsg_fluid
           ipg = glocal_index_msg_recv_fluid(ip,imsg)

!          add received buffer to own field at same global point
           gvec_fluid(ipg) = gvec_fluid(ipg) + buffr_fluid(ip)

!          assuming here that each global point is mapped one2one from 
!          the buffer... i.e. each ip has one ipg
           buffr_fluid(ip) = gvec_fluid(ipg)
        end do
!       send joint data back, but stick into new envelope/msgnum
        msgnum2 = mynum*nproc + ipsrc
        call MPI_SEND(buffr_fluid,sizeb,mpi_realkind,&
                      ipsrc,msgnum2,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive updated data back
  if (sizesend_fluid>0) then
     do imsg =1, sizesend_fluid
        sizemsg_fluid = sizemsgsend_fluid(imsg)
        sizeb = sizemsg_fluid
        ipsrc = listsend_fluid(imsg)
        msgnum = ipsrc*nproc + mynum
        call MPI_RECV(buffs_fluid,sizeb,mpi_realkind,&
                      ipsrc,msgnum,MPI_COMM_WORLD,status,IERROR)
        do ip = 1, sizemsg_fluid
           ipg = glocal_index_msg_send_fluid(ip,imsg)
           gvec_fluid(ipg) = buffs_fluid(ip)
        enddo
     enddo
  endif

end subroutine asynch_messaging_fluid
!=============================================================================

!-----------------------------------------------------------------------------
subroutine testing_asynch_messaging_fluid
!
! Fluid asynchronous communication pattern with one message per proc-proc pair.
! for a scalar field gvec. The arrays to map global numbers along 
! processor-processor boundaries are determined in the mesher, routine pdb.f90
! (consulation thereof to be avoided if at all possible...)
! Same as above but with output, used solely for the mpi test.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_comm, ONLY : buffr_fluid,buffs_fluid

include 'mesh_params.h'

integer                            :: imsg,ipg,ip,sizeb,ipdes,ipsrc
integer                            :: msgnum,msgnum1,msgnum2
integer                            :: status(MPI_STATUS_SIZE),sizemsg_fluid

! Prepare arrays to be sent.... MIGHT USE A POWER OF 2 STATEMENT THERE
  write(69,*)' Asynchrounous fluid communication test:'
  call flush(69)

  if (sizesend_fluid>0) then
     do imsg = 1, sizesend_fluid
        sizemsg_fluid = sizemsgsend_fluid(imsg)
        do ip = 1, sizemsg_fluid
           ipg = glocal_index_msg_send_fluid(ip,imsg)
           buffs_fluid(ip) = gvec_fluid(ipg)
        end do
        
! Send stuff around
        sizeb  = sizemsg_fluid
        ipdes  = listsend_fluid(imsg)
        msgnum = mynum*nproc + ipdes
        write(69,12)procstrg,'SENDING #',msgnum,sizemsg_fluid,' to',ipdes
        call MPI_SEND(buffs_fluid,sizeb,mpi_realkind,&
                      ipdes,msgnum,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive data, sum things up and send back to initial sender
  if (sizerecv_fluid>0) then 
     do imsg = 1, sizerecv_fluid
        sizemsg_fluid = sizemsgrecv_fluid(imsg)
        sizeb = sizemsg_fluid
        ipsrc = listrecv_fluid(imsg)
        msgnum1 = ipsrc*nproc + mynum
        write(69,12)procstrg,'RECEIVING #',msgnum1,sizemsg_fluid,' from',ipsrc
        call MPI_RECV(buffr_fluid,sizeb,mpi_realkind,&
                      ipsrc,msgnum1,MPI_COMM_WORLD,status,IERROR)
!     Add received buffer to own field at same global point
           do ip = 1, sizemsg_fluid
              ipg = glocal_index_msg_recv_fluid(ip,imsg)
!           assuming here that each global point is mapped one2one from 
!           the buffer... i.e. each ip has one ipg
               gvec_fluid(ipg) = gvec_fluid(ipg) + buffr_fluid(ip)
               buffr_fluid(ip) = gvec_fluid(ipg)
           end do
!     Send joint data back, but stick into new envelope/msgnum
        msgnum2 = mynum*nproc + ipsrc
        write(69,12)procstrg,'RETURNING #',msgnum2,sizemsg_fluid,' to',ipsrc
        call MPI_SEND(buffr_fluid,sizeb,mpi_realkind,&
                      ipsrc,msgnum2,MPI_COMM_WORLD,IERROR)
     end do
  endif

! Receive updated data back
  if (sizesend_fluid>0) then
     do imsg =1, sizesend_fluid
        sizemsg_fluid = sizemsgsend_fluid(imsg)
        sizeb = sizemsg_fluid
        ipsrc = listsend_fluid(imsg)
        msgnum = ipsrc*nproc + mynum
        write(69,12)procstrg,'RECV UPDATED #',msgnum,&
                    sizemsg_fluid,' from',ipsrc
        call MPI_RECV(buffs_fluid,sizeb,mpi_realkind,&
             ipsrc,msgnum,MPI_COMM_WORLD,status,IERROR)
        do ip = 1, sizemsg_fluid
           ipg = glocal_index_msg_send_fluid(ip,imsg)
           gvec_fluid(ipg) = buffs_fluid(ip)
        enddo
     enddo
  endif

  write(69,13)procstrg; call flush(69)

12 format('   MPI: ',a8,a15,i3,', size=',i8,a6,i3)
13 format('   MPI: ',a8,' DONE.')

end subroutine testing_asynch_messaging_fluid
!=============================================================================


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!===================
end module commpi
!===================
