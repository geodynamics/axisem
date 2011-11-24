! Proc   1: Header for mesh information to run static solver
! created by the mesher on 10/15/2011, at 16h 47min
 
!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :          prem_light
!   Inner-core shear wave:         T
!   Dominant period [s]  :   40.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =      7448  !                   proc. els
 integer, parameter ::       npoint =    186200  !               proc. all pts
 integer, parameter ::    nel_solid =      6440  !             proc. solid els
 integer, parameter ::    nel_fluid =      1008  !             proc. fluid els
 integer, parameter :: npoint_solid =    161000  !             proc. solid pts
 integer, parameter :: npoint_fluid =     25200  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =     16521  !            proc. flocal pts
 integer, parameter ::     nel_bdry =       168  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =         9  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =         1  !        number of processors
 integer, parameter :: lfbkgrdmodel =        10  !   length of bkgrdmodel name
 
!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.4375
!   Min(h/vp),dt/courant :    2.9165    2.9165
!   max(h/vs),T0/wvlngth :   26.0919   26.6667
!   Inner core r_min [km]:  936.0116
!   Max(h) r/ns(icb) [km]:   68.5260
!   Max(h) precalc.  [km]:   71.3721
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
character(len=100) :: meshpath='/Users/tnm/CODE/AXISEM/SOLVER/MESHES/PREMLIGHT_40S_NP1'
