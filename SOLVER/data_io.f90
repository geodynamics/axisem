!=================
module data_io
!=================
!
! Miscellaneous variables relevant to any read/write process such as 
! paths, logicals describing what to save, sampling rate of dumps

implicit none
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  character(len=80) :: datapath,infopath
  integer           :: lfdata,lfinfo
  logical           :: save_large_tests 
  logical           :: dump_energy
  logical           :: dump_snaps_glob
  logical           :: dump_snaps_solflu
  logical           :: dump_wavefields
  logical           :: need_fluid_displ
  double precision  :: strain_samp
  integer           :: istrain,isnap
  character(len=12) :: dump_type
  character(len=8)  :: rec_file_type
  logical           :: correct_azi,sum_seis,sum_fields
  character(len=3)  :: rot_rec
  logical           :: srcvic,add_hetero,file_exists 

! indices to limit dumping to select contiguous range of GLL points:
! 0<=ibeg<=iend<=npol
! For the time being: dump the same in xeta and eta directions
  integer           :: ibeg,iend
! ndumppts_el=(iend-ibeg+1)**2
  integer           :: ndumppts_el

! rotations
  double precision :: rot_mat(3,3),trans_rot_mat(3,3)
  double precision, allocatable, dimension(:,:) :: recfac

!af
  character(len=80), dimension(:), allocatable :: fname_rec_seis
  character(len=80), dimension(:), allocatable :: fname_rec_velo
!end af
  

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=====================
end module data_io
!=====================
