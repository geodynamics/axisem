!===================
module data_proc
!===================
!
! General variables pertaining to process identification

  implicit none
  public 

  integer          :: nproc              ! Number of total processors
  integer          :: mynum              ! Local processor label, from 0 to nproc-1
  character(len=4) :: appnproc, appmynum ! processor-identifying file extension
  logical          :: lpr                ! last processor logical flag, for write stdout 
  character(len=8) :: procstrg           ! String containing mynum to include in writes

!=======================
end module data_proc
!=======================
