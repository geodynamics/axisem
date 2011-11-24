!==================
module data_time
!==================

use global_parameters
implicit none
public 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  double precision    :: enforced_dt,enforced_period
  character(len=8)    :: time_scheme ! time extrapolation scheme, better be: 
                                     ! symplec4,sympqua4,symp_5_4,newmark2
  double precision    :: period      ! Dominant source period, given by mesher
  double precision    :: courant     ! Courant number, given by mesher
  double precision    :: t           ! current time
  double precision    :: deltat      ! Time step size
  double precision    :: seislength_t! seismogram length in seconds
  integer             :: niter       ! Number of iterations
  integer             :: iclockold,idold ! tick labels for timer
  integer             :: iclockcomm,idcomm ! tick labels for mpi timer
  integer             :: iclockmpi,idmpi ! tick labels for mpi timer
  integer             :: iclockstiff,idstiff ! tick labels for mpi timer
  integer             :: iclockdump,iddump ! tick labels for mpi timer
  double precision    :: seis_dt     ! seismogram sampling rate in seconds
  integer             :: seis_it     ! seismogram sampling rate in time steps
  double precision    :: snap_dt     ! time interval between snaps in seconds
  integer             :: snap_it     ! equivalent as snap_dt in time steps
  integer             :: strain_it   ! strain dump interval in time steps
  double precision    :: half_dt,half_dt_sq
  integer             :: nstages     ! number of substages in symplectic schemes
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!======================
end module data_time
!======================
