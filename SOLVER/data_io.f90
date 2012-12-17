!> Miscellaneous variables relevant to any read/write process such as 
!! paths, logicals describing what to save, sampling rate of dumps
!=================
module data_io
!=================
!
! Miscellaneous variables relevant to any read/write process such as 
! paths, logicals describing what to save, sampling rate of dumps

implicit none
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  character(len=200) :: datapath,infopath
  integer           :: lfdata,lfinfo
  logical           :: save_large_tests 
  logical           :: dump_energy
  logical           :: dump_snaps_glob
  logical           :: dump_xdmf
  logical           :: dump_snaps_solflu
  !> N.B. This is not wavefield snapshots, but kernel wavefields. Belongs to nstrain and
  !! istrain
  logical           :: dump_wavefields 
  logical           :: need_fluid_displ
  double precision  :: strain_samp
  integer           :: iseismo  !< current seismogram sample
  integer           :: istrain  !< current kernel wavefield sample
  integer           :: isnap    !< current wavefield sample (for plotting)
  integer           :: nseismo  !< Number of seismogram samples
  integer           :: nstrain  !< Number of wavefield dumps for kernels
  integer           :: nsnap    !< Number of wavefield snapshots
  character(len=12) :: dump_type
  character(len=8)  :: rec_file_type
  logical           :: correct_azi,sum_seis,sum_fields
  character(len=3)  :: rot_rec
  logical           :: srcvic,add_hetero,file_exists,use_netcdf 
  character(len=6)  :: output_format  !< netcdf or binary
  logical           :: force_ani

  
! indices to limit dumping to select contiguous range of GLL points:
! 0<=ibeg<=iend<=npol
! For the time being: dump the same in xeta and eta directions
  integer           :: ibeg,iend
! ndumppts_el=(iend-ibeg+1)**2
  integer           :: ndumppts_el

! for xdmf dumps
  integer           :: i_n_xdmf, j_n_xdmf
  integer, allocatable :: i_arr_xdmf(:), j_arr_xdmf(:)

! rotations
  double precision :: rot_mat(3,3),trans_rot_mat(3,3)
  double precision, allocatable, dimension(:,:) :: recfac

  character(len=80), dimension(:), allocatable :: fname_rec_seis
  character(len=80), dimension(:), allocatable :: fname_rec_velo


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!=====================
end module data_io
!=====================
