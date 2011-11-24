!========================
module parameters
!========================
!
! Read parameters for the general solver (i.e. NOT mesh, sources, receivers);
! compute other parameters for the simulation;
! write out summaries of all relevant simulation settings.
!

use global_parameters
use data_mesh 
use data_mesh_preloop 
use data_proc
use data_time 
use data_source
use data_io
use utlity
use commun

implicit none

public :: open_local_param_file,readin_parameters
public :: compute_numerical_parameters,write_parameters
private

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine open_local_param_file

  open(unit=69,file='output_proc'//appmynum//'.dat')
  write(69,*)
  write(69,*)'********** This is the OUTPUT for ',procstrg,&
             '***********************'
  write(69,*)
  write(69,*)' CONTAINS all information on the mesh, background model,'
  write(69,*)' source, receivers, precomputed matrices, and various tests '
  write(69,*)' such as surface areas, volume, valence, discontinuities,'
  write(69,*)' resolution test, axial masking, solid-fluid boundary copying,'
  write(69,*)' Lagrange interpolants, integration weights, message-passing.'
  write(69,*)
  write(69,*)'****************************************************************'
  write(69,*)

end subroutine open_local_param_file
!=============================================================================

!-----------------------------------------------------------------------------
subroutine readin_parameters_old
!
! Routine that reads in simulation parameters that are relevant at the 
! stage of the solver, i.e. number of time steps, integration scheme, 
! data paths, specification of wavefield dumping etc.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include 'mesh_params.h'
character(len=30) :: junk
integer :: i,ij

  open(5,file='inparam',POSITION='REWIND')
    read(5,10)datapath
    read(5,10)infopath
    read(5,*)num_simul
    read(5,*)seislength_t
    read(5,*)enforced_dt
    read(5,*)enforced_period
    read(5,*)src_file_type
    read(5,*)rec_file_type
    read(5,*)correct_azi
    read(5,*)sum_seis
    read(5,*)sum_fields
    read(5,*)rot_rec
    read(5,*)time_scheme
    read(5,*)seis_dt
    read(5,*)save_large_tests
    read(5,*)dump_energy
    read(5,*)dump_snaps_glob
    read(5,*)dump_snaps_solflu
    read(5,*)snap_dt
    read(5,*)dump_wavefields
    read(5,*)dump_type
    read(5,*)ibeg
    read(5,*)iend
    read(5,*)strain_samp
    read(5,*)src_dump_type
    read(5,*)make_homo
    read(5,*)vphomo,vshomo,rhohomo
    read(5,*)srcvic
    read(5,*)add_hetero
    read(5,*)do_mesh_tests
  close(5)
!af test
 vphomo = vphomo*1.e3
 vshomo = vshomo*1.e3
 rhohomo = rhohomo*1.e3
!
 iend = npol-iend

if (src_dump_type=='anal') then
   write(6,*)''
   write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(6,*)''
   write(6,*)'Analytical source wavefield dump not implemented YET!'
   write(6,*)'          DOING NOTHING INSTEAD......................'
   write(6,*)''
   write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(6,*)''
endif

10 format(a80)

  lfdata = index(datapath,' ')-1
  lfinfo = index(infopath,' ')-1

  call barrier
  if (lpr) then
     write(6,*)
     write(6,20)
     write(6,21)trim(datapath),trim(infopath),num_simul, seislength_t,enforced_dt, &
                enforced_period, &
                src_file_type,rec_file_type,correct_azi,sum_seis,sum_fields, &
                rot_rec,time_scheme,seis_dt,save_large_tests, &
                dump_energy,dump_snaps_glob,dump_snaps_solflu,dump_wavefields,&
                dump_type,ibeg,iend,strain_samp,src_dump_type,make_homo,srcvic, &
                add_hetero,do_mesh_tests

20 format(08x,&
       '///////////////////////////////////////////////////////////////',/&
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//                  A   X   I   S   E   M                    //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                  +-+      //',/  &
   08x,'//   Parallel spectral-element computation of    +-+---+-+   //',/  &
   08x,'//                                               | |   | |   //',/  &
   08x,'//     3-D seismic wave propagation for          | +---+ |   //',/  &
   08x,'//                                               |   _   |   //',/  &
   08x,'//    spherically symmetric background models    |_/_\\__|   //',/  &
   08x,'//                                               | \\_/  |   //',/  &
   08x,'//           in a global, 2-D domain             |       |   //',/  &
   08x,'//                                               | +---+ |   //',/  &
   08x,'//                                               | |   | |   //',/  &
   08x,'//                                               +-+---+-+   //',/  &
   08x,'//                                                  +-+      //',/  &
   08x,'//                                                           //',/  &
   08x,'//  Authors : Tarje Nissen-Meyer (tarje@princeton.edu)       //',/  &
   08x,'//            Alexandre Fournier (Grenoble)                  //',/  &
   08x,'//            Tony Dahlen (Princeton)                        //',/  &
   08x,'//                                                           //',/  &
   08x,'//       Comprehensive description of the underlying         //',/  &
   08x,'//           numerical analysis can be found in:             //',/  &
   08x,'//                                                           //')

21 format(08x,&
       '// (1) Tarje Nissen-Meyer, F. A. Dahlen, A Fournier (2007)   //',/&
   08x,'//     "Spherical-earth Frechet sensitivity kernels"         //',/  & 
   08x,'//     Geophysical Journal International 168(3),1051-1066.   //',/  & 
   08x,'//     doi:10.1111/j.1365-246X.2006.03123.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (2) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  & 
   08x,'//     "A two-dimensional spectral-element method for        //',/  &  
   08x,'//     spherical-earth seismograms-I. Moment-tensor source"  //',/  & 
   08x,'//     Geophysical Journal International 168(3), 1067-1092.  //',/  & 
   08x,'//     doi:10.1111/j.1365-246X.2006.03121.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (3) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  &
   08x,'//     "A two-dimensional spectral-element method for        //',/  &
   08x,'//     spherical-earth seismograms - II. Background models"  //',/  &
   08x,'//     submitted to Geophysical Journal International.       //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//  May 2007 : Version 1.1, includes a                       //',/  &
   08x,'//                                                           //',/  &
   08x,'///////////////////////////////////////////////////////////////',// &
   08x,'=============  I N P U T    P A R A M E T E R S ==============',/ &
   12x,'Data I/O path:                  ',a20,/                         &
   12x,'Info I/O path:                  ',a20,/                         &
   12x,'Number of source simulations:   ',i2,/                          &
   12x,'Simulation length [s]:          ',f9.3,/                          &
   12x,'Enforced time step [s]:         ',f7.3,/                        &
   12x,'Enforced source period [s]:     ',f7.3,/                        &
   12x,'Source file type:               ',a8,/                        &
   12x,'Receiver file type:             ',a8,/                        &
   12x,'Correct azimuth?                ',l2,/                          &
   12x,'Sum seismograms?                ',l2,/                          &
   12x,'Sum wavefields?                 ',l2,/                          &
   12x,'Receivers coordinates                 ',a3,/                          &
   12x,'Time extrapolation scheme:      ',a8,/                          &
   12x,'Seismogram sampling rate [s]:   ',f7.3,/                        &
   12x,'Save large tests?               ',l2,/                          &
   12x,'Dump kin./pot. energy?          ',l2,/                          &
   12x,'Dump global snaps?              ',l2,/                          &
   12x,'Dump solid/fluid snaps?         ',l2,/                          &
   12x,'Dump strain?                    ',l2,/                          &
   12x,'Wavefield dumping type:         ',a12,/                         &
   12x,'First GLL to save in strains:   ',i2,/                          &
   12x,'Last GLL to save in strains:    ',i2,/                          &
   12x,'Samples per period for strains: ',f7.3,/                        &
   12x,'Source dumping type:            ',a4,/                          &
   12x,'Homogenize background model?    ',l2,/                          &
   12x,'Analyt. homogen. radiation?     ',l2,/                          &
   12x,'Add heterogeneous region?       ',l2,/                          &
   12x,'Perform extensive mesh tests?       ',l2,/                          &
   08x,'==============================================================')
  write(6,*)
  write(6,*)'Processor-specific output is written to: output_proc<PROC ID>.dat'
  write(6,*)'All potential error messages will appear here...'
  endif !lpr

! Checking the consistency of some of the input parameters
  if ( mod(realkind,4)/=0 .or. realkind>8) then
     if (lpr) then
        write(6,*)
        write(6,*)'PROBLEM with REAL data kind!'
        write(6,*)'... can only handle real kinds 4 or 8.'
        write(6,*)'real kind here:', realkind
        write(6,*)'change parameter realkind in global_parameters.f90'
     endif
     stop
  endif

  if (strain_samp> 10) then
     if (lpr) then     
        write(6,*)
        write(6,*)"!!!!!! NOT GOING ANY FURTHER !!!!!!"
        write(6,*)"  It's just too much to save 10 frames of strain & velocity"
        write(6,*)"  per source period! Choose something reasonable."
     endif
     stop
  endif

  if (enforced_dt > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14)'maximal time step',enforced_dt
     endif
  endif

  if (enforced_period > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14)'min. source period',enforced_period
     endif
  endif

14 format('  WARNING: Overriding',a19,' with:',f8.3,' seconds')

  if (lpr) then
  if (dump_snaps_glob .and. dump_snaps_solflu) then 
     write(6,*)''
     write(6,*)" NOT dumping the same snapshots twice (global AND solid/fluid)"
     write(6,*)'...hence reverting to dumping global snaps only. Sorry.'
     dump_snaps_solflu=.false.
  endif

  if (srcvic) then 
     if (.not. make_homo ) then 
     if (lpr) then
     write(6,*)
     write(6,7)'0000000000000000 WARNING ABOUT PERFORMANCE 0000000000000000000'
     write(6,7)'00                                                          00'
     write(6,7)'00    Computing analyt. radiation for heterogeneous model?  00'
     write(6,7)'00         ...kinda silly, therefore turning it off...      00'
     write(6,7)'00                                                          00'
     write(6,7)'00000000000000000000000000000000000000000000000000000000000000'
     endif
     srcvic=.false.
     endif
  endif

7 format(04x,a62)

  if (realkind==4) then       
     if (lpr) then
     write(6,7)
     write(6,7)'44444444444444444444444444444444444444444444444444444444444444'
     write(6,7)'444   Running the solver time loop with SINGLE PRECISION   444'
     write(6,7)'44444444444444444444444444444444444444444444444444444444444444'
     endif
  elseif (realkind==8) then       
     if (lpr) then
     write(6,7)
     write(6,7)'88888888888888888888888888888888888888888888888888888888888888'
     write(6,7)'888   Running the solver time loop with DOUBLE PRECISION   888'
     write(6,7)'88888888888888888888888888888888888888888888888888888888888888'
     endif
  endif
  write(6,*)

  endif !mynum

! Need to decide here since this boolean is needed in def_precomp_terms
  need_fluid_displ = .false.
  if (dump_snaps_glob .or. dump_snaps_solflu .or. dump_energy .or. & 
       dump_wavefields .and. dump_type=='fullfields') then
! Need to add this for each new type of wavefield dumping method that 
! requires the fluid displacement/velocities
     need_fluid_displ = .true.
  endif
 

! define general small value

if (realkind==4) then
    smallval=smallval_sngl
elseif (realkind==8) then
   smallval=smallval_dble
endif
if (lpr) write(6,*)'  small value is:',smallval

!!$! derive number of simulations based on moment tensor entries
!!$if (src_file_type=='cmtsolut') then
!!$      open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk
!!$     read(20000,*)junk,Mij(1) !Mrr
!!$     read(20000,*)junk,Mij(2) !Mtt
!!$     read(20000,*)junk,Mij(3) !Mpp
!!$     read(20000,*)junk,Mij(4) !Mrt
!!$     read(20000,*)junk,Mij(5) !Mrp
!!$     read(20000,*)junk,Mij(6) !Mtp
!!$     close(20000)  
!!$
!!$! convert to [Nm]
!!$     Mij = Mij/1.E7
!!$     num_simul = 0
!!$     do i=1,6
!!$        if (Mij(i)>smallval*maxval(abs(Mij))) num_simul = num_simul + 1        
!!$     enddo
!!$     if (num_simul==6) then 
!!$        num_simul = 4
!!$     elseif (num_simul==1) then
!!$        ! nothing to do, will be dealt with in source.f90
!!$     else
!!$        write(6,*)'  havent done this case of less non-zero moment tensor elements than 6...'
!!$        write(6,*)'  ... just computing the full set of 4 simulations instead...'
!!$        num_simul = 4
!!$     endif
!!$endif

end subroutine readin_parameters_old
!=============================================================================

!-----------------------------------------------------------------------------
subroutine readin_parameters
!
! Routine that reads in simulation parameters that are relevant at the 
! stage of the solver, i.e. number of time steps, integration scheme, 
! data paths, specification of wavefield dumping etc.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include 'mesh_params.h'
character(len=100) :: junk
integer :: i

  open(5,file='inparam',POSITION='REWIND')
    read(5,*)junk
    read(5,*)num_simul;write(6,*)num_simul
    read(5,*)junk

    read(5,*)seislength_t
    read(5,*)enforced_dt
    read(5,*)time_scheme
    read(5,*)junk

    read(5,*)enforced_period
    read(5,*)src_file_type
    read(5,*)rec_file_type
    read(5,*)seis_dt
    read(5,*)junk

    read(5,10)datapath; datapath=trim(datapath)
    read(5,10)infopath; infopath=trim(infopath)
    read(5,*)dump_snaps_glob
    read(5,*)snap_dt
    read(5,*)junk

    read(5,*)dump_wavefields
    read(5,*)strain_samp
    read(5,*)src_dump_type
    read(5,*)ibeg,iend
    read(5,*)junk

    read(5,*)dump_energy
    read(5,*)make_homo
    read(5,*)vphomo,vshomo,rhohomo
    read(5,*)srcvic
    read(5,*)add_hetero
    read(5,*)do_mesh_tests
    read(5,*)save_large_tests
  close(5)

! now pre-set. Most of these are to be considered in the post processing stage now.
   correct_azi=.false.
   sum_seis=.false.
   sum_fields=.false.
   rot_rec='cyl'
   dump_snaps_solflu=.false.
   dump_type='fullfields'

!af test
 vphomo = vphomo*1.e3
 vshomo = vshomo*1.e3
 rhohomo = rhohomo*1.e3
!
 iend = npol-iend

if (src_dump_type=='anal') then
   write(6,*)''
   write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(6,*)''
   write(6,*)'Analytical source wavefield dump not implemented YET!'
   write(6,*)'          DOING NOTHING INSTEAD......................'
   write(6,*)''
   write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(6,*)''
endif

10 format(a80)

  lfdata = index(datapath,' ')-1
  lfinfo = index(infopath,' ')-1

  call barrier
  if (lpr) then
     write(6,*)
     write(6,20)
     write(6,21)datapath,infopath,num_simul, seislength_t,enforced_dt, &
                enforced_period, &
                src_file_type,rec_file_type,correct_azi,sum_seis,sum_fields, &
                rot_rec,time_scheme,seis_dt,save_large_tests, &
                dump_energy,dump_snaps_glob,dump_snaps_solflu,dump_wavefields,&
                dump_type,ibeg,iend,strain_samp,src_dump_type,make_homo,srcvic, &
                add_hetero,do_mesh_tests

20 format(08x,&
       '///////////////////////////////////////////////////////////////',/&
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//                  A   X   I   S   E   M                    //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                  +-+      //',/  &
   08x,'//   Parallel spectral-element computation of    +-+---+-+   //',/  &
   08x,'//                                               | |   | |   //',/  &
   08x,'//     3-D seismic wave propagation for          | +---+ |   //',/  &
   08x,'//                                               |   _   |   //',/  &
   08x,'//    spherically symmetric background models    |_/_\\__|   //',/  &
   08x,'//                                               | \\_/  |   //',/  &
   08x,'//           in a global, 2-D domain             |       |   //',/  &
   08x,'//                                               | +---+ |   //',/  &
   08x,'//                                               | |   | |   //',/  &
   08x,'//                                               +-+---+-+   //',/  &
   08x,'//                                                  +-+      //',/  &
   08x,'//                                                           //',/  &
   08x,'//  Authors : Tarje Nissen-Meyer (tarje@princeton.edu)       //',/  &
   08x,'//            Alexandre Fournier (Grenoble)                  //',/  &
   08x,'//            Tony Dahlen (Princeton)                        //',/  &
   08x,'//                                                           //',/  &
   08x,'//       Comprehensive description of the underlying         //',/  &
   08x,'//           numerical analysis can be found in:             //',/  &
   08x,'//                                                           //')

21 format(08x,&
       '// (1) Tarje Nissen-Meyer, F. A. Dahlen, A Fournier (2007)   //',/&
   08x,'//     "Spherical-earth Frechet sensitivity kernels"         //',/  & 
   08x,'//     Geophysical Journal International 168(3),1051-1066.   //',/  & 
   08x,'//     doi:10.1111/j.1365-246X.2006.03123.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (2) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  & 
   08x,'//     "A two-dimensional spectral-element method for        //',/  &  
   08x,'//     spherical-earth seismograms-I. Moment-tensor source"  //',/  & 
   08x,'//     Geophysical Journal International 168(3), 1067-1092.  //',/  & 
   08x,'//     doi:10.1111/j.1365-246X.2006.03121.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (3) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  &
   08x,'//     "A two-dimensional spectral-element method for        //',/  &
   08x,'//     spherical-earth seismograms - II. Background models"  //',/  &
   08x,'//     submitted to Geophysical Journal International.       //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//  May 2007 : Version 1.1, includes a                       //',/  &
   08x,'//                                                           //',/  &
   08x,'///////////////////////////////////////////////////////////////',// &
   08x,'=============  I N P U T    P A R A M E T E R S ==============',/ &
   12x,'Data I/O path:                  ',a20,/                         &
   12x,'Info I/O path:                  ',a20,/                         &
   12x,'Number of source simulations:   ',i2,/                          &
   12x,'Simulation length [s]:          ',f9.3,/                          &
   12x,'Enforced time step [s]:         ',f7.3,/                        &
   12x,'Enforced source period [s]:     ',f7.3,/                        &
   12x,'Source file type:               ',a7,/                        &
   12x,'Receiver file type:             ',a8,/                        &
   12x,'Correct azimuth?                ',l2,/                          &
   12x,'Sum seismograms?                ',l2,/                          &
   12x,'Sum wavefields?                 ',l2,/                          &
   12x,'Receivers coordinates                 ',a3,/                          &
   12x,'Time extrapolation scheme:      ',a8,/                          &
   12x,'Seismogram sampling rate [s]:   ',f7.3,/                        &
   12x,'Save large tests?               ',l2,/                          &
   12x,'Dump kin./pot. energy?          ',l2,/                          &
   12x,'Dump global snaps?              ',l2,/                          &
   12x,'Dump solid/fluid snaps?         ',l2,/                          &
   12x,'Dump strain?                    ',l2,/                          &
   12x,'Wavefield dumping type:         ',a12,/                         &
   12x,'First GLL to save in strains:   ',i2,/                          &
   12x,'Last GLL to save in strains:    ',i2,/                          &
   12x,'Samples per period for strains: ',f7.3,/                        &
   12x,'Source dumping type:            ',a4,/                          &
   12x,'Homogenize background model?    ',l2,/                          &
   12x,'Analyt. homogen. radiation?     ',l2,/                          &
   12x,'Add heterogeneous region?       ',l2,/                          &
   12x,'Perform extensive mesh tests?       ',l2,/                          &
   08x,'==============================================================')
  write(6,*)
  write(6,*)'Processor-specific output is written to: output_proc<PROC ID>.dat'
  write(6,*)'All potential error messages will appear here...'
  endif !lpr

! Checking the consistency of some of the input parameters
  if ( mod(realkind,4)/=0 .or. realkind>8) then
     if (lpr) then
        write(6,*)
        write(6,*)'PROBLEM with REAL data kind!'
        write(6,*)'... can only handle real kinds 4 or 8.'
        write(6,*)'real kind here:', realkind
        write(6,*)'change parameter realkind in global_parameters.f90'
     endif
     stop
  endif

  if (strain_samp> 10) then
     if (lpr) then     
        write(6,*)
        write(6,*)"!!!!!! NOT GOING ANY FURTHER !!!!!!"
        write(6,*)"  It's just too much to save 10 frames of strain & velocity"
        write(6,*)"  per source period! Choose something reasonable."
     endif
     stop
  endif

  if (enforced_dt > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14)'maximal time step',enforced_dt
     endif
  endif

  if (enforced_period > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14)'min. source period',enforced_period
     endif
  endif

14 format('  WARNING: Overriding',a19,' with:',f8.3,' seconds')

  if (lpr) then
  if (dump_snaps_glob .and. dump_snaps_solflu) then 
     write(6,*)''
     write(6,*)" NOT dumping the same snapshots twice (global AND solid/fluid)"
     write(6,*)'...hence reverting to dumping global snaps only. Sorry.'
     dump_snaps_solflu=.false.
  endif

  if (srcvic) then 
     if (.not. make_homo ) then 
     if (lpr) then
     write(6,*)
     write(6,7)'0000000000000000 WARNING ABOUT PERFORMANCE 0000000000000000000'
     write(6,7)'00                                                          00'
     write(6,7)'00    Computing analyt. radiation for heterogeneous model?  00'
     write(6,7)'00         ...kinda silly, therefore turning it off...      00'
     write(6,7)'00                                                          00'
     write(6,7)'00000000000000000000000000000000000000000000000000000000000000'
     endif
     srcvic=.false.
     endif
  endif

7 format(04x,a62)

  if (realkind==4) then       
     if (lpr) then
     write(6,7)
     write(6,7)'44444444444444444444444444444444444444444444444444444444444444'
     write(6,7)'444   Running the solver time loop with SINGLE PRECISION   444'
     write(6,7)'44444444444444444444444444444444444444444444444444444444444444'
     endif
  elseif (realkind==8) then       
     if (lpr) then
     write(6,7)
     write(6,7)'88888888888888888888888888888888888888888888888888888888888888'
     write(6,7)'888   Running the solver time loop with DOUBLE PRECISION   888'
     write(6,7)'88888888888888888888888888888888888888888888888888888888888888'
     endif
  endif
  write(6,*)

  endif !mynum

! Need to decide here since this boolean is needed in def_precomp_terms
  need_fluid_displ = .false.
  if (dump_snaps_glob .or. dump_snaps_solflu .or. dump_energy .or. & 
       dump_wavefields .and. dump_type=='fullfields') then
! Need to add this for each new type of wavefield dumping method that 
! requires the fluid displacement/velocities
     need_fluid_displ = .true.
  endif

! define general small value

if (realkind==4) then
    smallval=smallval_sngl
elseif (realkind==8) then
   smallval=smallval_dble
endif
if (lpr) write(6,*)'  small value is:',smallval

end subroutine readin_parameters
!=============================================================================

!-----------------------------------------------------------------------------
subroutine compute_numerical_parameters

include "mesh_params.h"

double precision :: s,z,r,theta,s_max
double precision :: dsaxis(0:npol-1,0:npol), dzaxis(0:npol-1) 
double precision :: minds(nelem),maxds(nelem),mindz(nelem),maxdz(nelem)
integer          :: ielem,ipol,jpol

  if (lpr) then
    write(6,*)
    write(6,*)'  Computing numerical parameters...'
  endif

! If source is not at the axis, then we just rotate receivers to spherical:
! This is done since then the rotation of the source-receiver system 
! is preserved in terms of components for the spherical-coordinate 
! based entries into the moment tensor, and the same for the receivers. 
! Were we to stay cylindrical, then a component rotation on the receiver 
! would be necessary upon rotating the source-receiver system to the north pole. 
! Is this clear? At least to TNM on Jan 30, 2011 ...
!  if (srccolat/=0.d0 .or. srclon/=0.d0 ) then 
!     if (lpr) then 
!        write(6,*)'  WARNING: Since source is not at the pole, we make sure that seismograms'
!        write(6,*)'                       are given in spherical coordinate components. Otherwise, this '
!        write(6,*)'                      may have necessitated additional rotations. See source code.'     
!     endif
!     rot_rec = 'sph'
!     if (lpr) write(6,*)'  .... changed receiver system to ',trim(rot_rec)
!  endif

! Override time step or source period if demanded by input
  if (enforced_dt>zero) then
     if (time_scheme=='newmark2' .and. enforced_dt>deltat .or. & 
          time_scheme=='symplec4' .and. enforced_dt>1.5*deltat .or. &
          time_scheme=='SS_35o10' .and. enforced_dt>3.0*deltat  ) then 
        if (lpr) then 
           write(6,*)
           write(6,*)'PROBLEM: Time step larger than allowed by mesh!'
           write(6,*)'Chosen value (in inparam file) [s] :',enforced_dt
           write(6,*)'Maximal time step for this mesh [s]:',deltat
           write(6,*)'Change time step in input file to lower than this max.'
           write(6,*)'or to zero to use precalculated (recommended)'
        endif
        stop
     else
        if (lpr) then 
           write(6,*)
           write(6,*)'    WARNING: Time step smaller than necessary by mesh!'
           write(6,20)enforced_dt,deltat
           write(6,19)100.-enforced_dt/deltat*100.
        endif
        deltat=enforced_dt
     endif
  else
     if (lpr) then 
        write(6,*)
        write(6,*)'    Using time step precalculated by the mesher:',deltat
     endif
  endif
20 format('     Chosen/maximal time step [s]:',2(f7.3))
19 format('     ...lengthens this simulation by',f6.2,' percent!')

! source period
  if (enforced_period>zero .and. stf_type/='dirac_0') then 
     if (enforced_period<period) then 
        if (lpr) then 
           write(6,*)
           write(6,*)'    WARNING: Period smaller than necessary by mesh!'
           write(6,21)'   Chosen value (in inparam file):',enforced_period
           write(6,21)'   Minimal period for this mesh  :',period
           write(6,*)'    Change period in input file to larger than this min.'
           write(6,*)'    or to zero to use precalculated (recommended)'
        endif
        t_0=enforced_period
     else
        if (lpr) then 
           write(6,*)
           write(6,*)'    WARNING: Using larger period than necessary by mesh!'
           write(6,23)enforced_period,period
        endif
        t_0=enforced_period
     endif
  else 
     if (stf_type/='dirac_0') then
        if (lpr) then
           write(6,*)
           write(6,*)'    Using period of the mesh:',period
        endif
        t_0=period
     else
        t_0=period ! Just for consistency -- is in fact never used for dirac.
     endif
  endif

21 format(a36,f8.3,' s')
23 format('     Chosen/minimal period   [s]:',2(f7.3))

! Compute number of iterations in time loop
  niter=floor((seislength_t+smallval_dble)/deltat)
  if (lpr) then 
     write(6,*)
     write(6,22)'    desired simulation length  :',seislength_t,' seconds'
     write(6,22)'    offered simulation length  :',niter*deltat,' seconds'
     write(6,11)'    number time loop iterations:',niter
  endif

! Compute seismogram sampling rate in time steps
  if (seis_dt > zero  .and. seis_dt >= deltat) then
     seis_it=floor((seis_dt+smallval_dble)/deltat)
  elseif (seis_dt < deltat) then
     if (lpr) write(6,*) 'seismogram sampling cannot be smaller than time step...'
     if (lpr) write(6,*) '...changing it to the time step'
     seis_dt = deltat
     seis_it = 1
  else
     seis_it = 1
  endif

  if (lpr) then
     write(6,*)
     write(6,22)'    desired seismogram sampling:',seis_dt,' seconds'
     write(6,22)'    offered seismogram sampling:',deltat*seis_it,' seconds'
     write(6,13)'    ...that is, every          :',seis_it,' timesteps'
     write(6,11)'    number of samples          :', &
                floor(real(niter)/real(seis_it))
  endif
22 format(a33,f9.2,a10)

! Initialize counters for I/O
  istrain = 0
  isnap = 0

! snapshot output, convert from interval given in seconds to 
! incremental time steps
   if (dump_snaps_glob .or. dump_snaps_solflu) then
     snap_it=floor(snap_dt/deltat)
     open(unit=2900+mynum,file=datapath(1:lfdata)//'/snap_info.dat'//appmynum)
     write(2900,*)floor(real(niter)/real(snap_it))
     do ielem=1,floor(real(niter)/real(snap_it))
        write(2900+mynum,*)real(ielem)*snap_dt,ielem*snap_it
     enddo
     close(2900+mynum)
     if (lpr) then
        write(6,*)
        write(6,11)'    Number of snapshots        :',&
                    floor(real(niter)/real(snap_it))
        write(6,12)'    ...approximately every     :',snap_dt,'seconds'
        write(6,13)'    ...that is, every          :',snap_it,'timesteps'
     endif
11 format(a33,i8)
12 format(a33,f8.2,a10)
13 format(a33,i8,a10)

  endif

! strain tensor output, convert from num of dumps per period into 
! incremental time steps
  if (dump_wavefields) then

     strain_it=floor(t_0/real(strain_samp)/deltat)
     open(unit=2900+mynum,file=datapath(1:lfdata)//'/strain_info.dat'//appmynum)
     write(2900,*)floor(real(niter)/real(strain_it))
     do ielem=1,floor(real(niter)/real(strain_it))
        write(2900+mynum,*)real(ielem)*t_0/real(strain_samp),ielem*strain_it
     enddo
     close(2900+mynum)
     if (lpr) then
        write(6,*)
        write(6,11)'    Number of wavefield dumps  :',&
                   floor(real(niter)/real(strain_it))
             
        write(6,12)'    ...approximately every     :', &
                   t_0/real(strain_samp),&
                   'seconds'
        write(6,13)'    ...that is, every          :',strain_it,'timestep'
     endif

     ndumppts_el=(iend-ibeg+1)**2
     if (lpr) then 
        write(6,*)'    Define limitation of GLL points in the dumped fields:'
        write(6,*)'      ibeg=',ibeg,'iend=',iend
        write(6,*)'      # points saved within an element:',ndumppts_el
     endif

  endif

! mesh info: coordinates of elements and collocation points               
  open(2222+mynum,file=infopath(1:lfinfo)//'/axial_points.dat'//appmynum)
  open(3333+mynum,file=infopath(1:lfinfo)//'/axial_ds_dz.dat'//appmynum)
  s_max=zero

! Set some parameters for faster access in time loop
  half_dt=half*deltat
  half_dt_sq=half*deltat**2

122 format(i8,i3,4(1pe15.5))

  do ielem = 1, nelem

!     write(6,*)'PARAMETERS EL LOOP before axial points:',ielem;call flush(6)

! write out axial points
     if (axis(ielem)) then
         call compute_coordinates(s,z,r,theta,ielem,0,npol)
        do jpol=0,npol
           call compute_coordinates(s,z,r,theta,ielem,0,jpol)
           write(2222+mynum,122)ielem,jpol,s,z,r,theta/pi*180.
        enddo

!     write(6,*)'PARAMETERS EL LOOP before grid spacing:',ielem;call flush(6)

! write out profile of grid spacing along Northern axis
        if (north(ielem)) then

         do jpol=0,npol
            do ipol=0,npol-1
               dsaxis(ipol,jpol) = dsqrt((scoord(ipol,jpol,ielem)-&
                    scoord(ipol+1,jpol,ielem))**2+&
                    (zcoord(ipol,jpol,ielem)-&
                    zcoord(ipol+1,jpol,ielem))**2 )
           enddo
         enddo
  
           do jpol=0,npol-1
           dzaxis(jpol) = dsqrt( (scoord(0,jpol,ielem) - &
                scoord(0,jpol+1,ielem) )**2 + &
                (zcoord(0,jpol,ielem) - &
                zcoord(0,jpol+1,ielem))**2 )
           enddo
           minds(naxel) = minval(dsaxis)
           maxds(naxel) = maxval(dsaxis)
           mindz(naxel) = minval(dzaxis)
           maxdz(naxel) = maxval(dzaxis)
           call compute_coordinates(s,z,r,theta,ielem,0,npol/2)        
           write(3333+mynum,123)ielem,r,minds(naxel),maxds(naxel), &
                mindz(naxel),maxdz(naxel),maxds(naxel)/minds(naxel)
        endif

     endif

123 format(i7,1pe14.3,5(1pe14.3))

  enddo

  close(2222+mynum)
  close(3333+mynum)

  if (lpr) write(6,*)


end subroutine compute_numerical_parameters
!=============================================================================

!-----------------------------------------------------------------------------
subroutine write_parameters

use data_comm
use data_numbering, ONLY : nglob,nglob_solid

include 'mesh_params.h'

integer          :: iel,curvel,linel,seminoel,semisoel,num_rec_glob
integer          :: curvel_solid,linel_solid,seminoel_solid,semisoel_solid
integer          :: curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid
integer          :: ipol,jpol,hmaxloc1(3),hminloc1(3),i,j
integer          :: maxprocssend_solid,maxprocsrecv_solid 
integer          :: maxprocssend_fluid,maxprocsrecv_fluid,nsim1
double precision :: dis1(0:npol-1,0:npol-1,nelem),dis2(0:npol-1,0:npol-1,nelem)
double precision :: s,z,r,theta,rminglob,thetaminglob,rmaxglob,thetamaxglob
double precision :: mysmin,myzmin,mysmax,myzmax
double precision :: myrmin,mythetamin,myrmax,mythetamax
double precision :: hmax,hmaxglob,hmin,hminglob
character(len=4) :: Mij_char(6)

  write(69,*)'  writing out all relevant simulation parameters...'
  call flush(69)

  write(69,*)'  number of respective element types...'; call flush(69)

  curvel=0; linel=0; seminoel=0; semisoel=0
  do iel=1,nelem
     if (eltype(iel)=='curved') curvel=curvel+1
     if (eltype(iel)=='linear') linel=linel+1
     if (eltype(iel)=='semino') seminoel=seminoel+1
     if (eltype(iel)=='semiso') semisoel=semisoel+1
  enddo

  curvel_solid=0; linel_solid=0; seminoel_solid=0; semisoel_solid=0
  do iel=1,nel_solid
     if (eltype(ielsolid(iel))=='curved') curvel_solid=curvel_solid+1
     if (eltype(ielsolid(iel))=='linear') linel_solid=linel_solid+1
     if (eltype(ielsolid(iel))=='semino') seminoel_solid=seminoel_solid+1
     if (eltype(ielsolid(iel))=='semiso') semisoel_solid=semisoel_solid+1
  enddo

  curvel_fluid=0; linel_fluid=0; seminoel_fluid=0; semisoel_fluid=0
  do iel=1,nel_fluid
     if (eltype(ielfluid(iel))=='curved') curvel_fluid=curvel_fluid+1
     if (eltype(ielfluid(iel))=='linear') linel_fluid=linel_fluid+1
     if (eltype(ielfluid(iel))=='semino') seminoel_fluid=seminoel_fluid+1
     if (eltype(ielfluid(iel))=='semiso') semisoel_fluid=semisoel_fluid+1
  enddo

  write(69,*)'  grid spacing min/max...'; call flush(69)
  do iel=1,nelem
     do ipol=0,npol-1
        do jpol=0,npol-1
           dis1(ipol,jpol,iel) = dsqrt(&
                (scoord(ipol,jpol,iel)-scoord(ipol+1,jpol,iel))**2&
                +(zcoord(ipol,jpol,iel)-zcoord(ipol+1,jpol,iel))**2) 
           dis2(ipol,jpol,iel) = dsqrt(&
                (scoord(ipol,jpol,iel)-scoord(ipol,jpol+1,iel))**2&
                +(zcoord(ipol,jpol,iel)-zcoord(ipol,jpol+1,iel))**2) 
        enddo
     enddo
  enddo

  write(69,*)'  calculating hmax...';call flush(69)
  hmax=max(maxval(dis1),maxval(dis2))
  write(69,*)'  hmaxstuff:',minval(dis1),minval(dis2),hmax;call flush(69)
  hmaxglob=pmax(hmax)
  write(69,*)'  hmaxglob:',hmaxglob;call flush(69)
  hmaxloc1=maxloc(dis1)
  if (maxval(dis2)<maxval(dis1)) hmaxloc1=maxloc(dis2)
  call compute_coordinates(s,z,rmaxglob,thetamaxglob,hmaxloc1(3), &
       hmaxloc1(1)-1,hmaxloc1(2)-1)
  write(69,*)' rmax,thetamax:',rmaxglob,thetamaxglob; call flush(69)

  write(69,*)'  calculating hmin...';call flush(69)
  hmin=min(minval(dis1),minval(dis2))
  write(69,*)'  hminstuff:',minval(dis1),minval(dis2),hmin;call flush(69)
  hminglob=pmin(hmin)
  write(69,*)'  hminglob:',hminglob;call flush(6)

  hminloc1=minloc(dis1)
  if (minval(dis2)<minval(dis1)) hminloc1=minloc(dis2)
  call compute_coordinates(s,z,rminglob,thetaminglob,hminloc1(3), &
          hminloc1(1)-1,hminloc1(2)-1)

! Checking potential issues with input parameter consistency
  if (lpr) write(6,*)
  if (lpr) write(6,*)'  checking input parameters for consistency...'
  call flush(6)
  call check_parameters(hmaxglob,hminglob,curvel,linel,seminoel,semisoel, &
                       curvel_solid,linel_solid,seminoel_solid,semisoel_solid,&
                       curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid)

  maxprocssend_solid = pmax_int(sizesend_solid)
  maxprocsrecv_solid = pmax_int(sizerecv_solid)
  maxprocssend_fluid = pmax_int(sizesend_fluid)
  maxprocsrecv_fluid = pmax_int(sizerecv_fluid)

! output to stdout, only by proc nproc-1
  if (lpr) then

     write(6,*)
     write(6,*)':::::::::::::::: SIMULATION PARAMETERS::::::::::::::::::::::::'

     write(6,*)'  Global mesh information______________________________'
     write(6,12)'     Background model  :',bkgrdmodel
     write(6,10)'     # discontinuities :',ndisc
     write(6,13)'     Have fluid region ?',have_fluid
     write(6,13)'     IC shear wave     ?',resolve_inner_shear
     write(6,11)'     Outer rad.     [m]:',router
     write(6,11)'     Inner rad.     [m]:',rmin
     write(6,10)'     Polynomial order  :',npol
     write(6,10)'     # control nodes   :',npoin
     write(6,10)'     Total elements    :',nelem    
     write(6,10)'     Total # points    :',npoint
     write(6,10)'     # global numbers  :',nglob
     write(6,10)'     # axial elements  :',naxel
     write(6,10)'     # curved elements :',curvel
     write(6,10)'     # linear elements :',linel
     write(6,10)'     # mixed elements  :',seminoel+semisoel
     write(6,11)'     Min. distance  [m]:',min_distance_dim
     write(6,11)'     Min. distance/r0  :',min_distance_nondim

     write(6,*)'  Grid spacing, velocities etc.________________________'
     write(6,17)'     Min. (pre,comp)[m]:',hmin_glob,hminglob
     write(6,17)'     Max. (pre,comp)[m]:',hmax_glob,hmaxglob
     write(6,17)'     Min. vp[m/s], r[m]:',vpmin,vpminr
     write(6,17)'     Min. vs[m/s], r[m]:',vsmin,vsminr
     write(6,17)'     Max. vp[m/s], r[m]:',vpmax,vpmaxr
     write(6,17)'     Max. vs[m/s], r[m]:',vsmax,vsmaxr
     write(6,11)'     Max. lead time [s]:',char_time_max
     write(6,17)'     r [m], theta [deg]:',char_time_max_rad*router,& 
                                           char_time_max_theta
     write(6,11)'     Min. lead time [s]:',char_time_min
     write(6,17)'     r [m], theta [deg]:',char_time_min_rad*router,& 
                                           char_time_min_theta

     write(6,*)'  Solid-Fluid configuration____________________________'
     write(6,15)'     S/F elements      :',nel_solid,nel_fluid  
     write(6,15)'     S/F # points      :',npoint_solid,npoint_fluid
     write(6,15)'     S/F global numbers:',nglob_solid,nglob_fluid
     write(6,15)'     S/F # axial elems :',naxel_solid,naxel_fluid
     write(6,10)'     # S/F boundary els:',nel_bdry
     write(6,15)'     S/F curved elems  :',curvel_solid,curvel_fluid
     write(6,15)'     S/F linear elems  :',linel_solid,linel_fluid
     write(6,15)'     S/F mixed elements:',seminoel_solid+semisoel_solid, &
                                           seminoel_fluid+semisoel_fluid

     write(6,*)'  Solid message passing_________________________________' 
     write(6,10)'     # processors      :',nproc
     write(6,10)'     max. sent messages:',maxprocssend_solid
     write(6,10)'     max. sent size    :',sizemsgsendmax_solid
     write(6,10)'     nax. recv messages:',maxprocsrecv_solid
     write(6,10)'     max. recv size    :',sizemsgrecvmax_solid

     if (have_fluid) then
        write(6,*)'  Fluid message passing_________________________________' 
        write(6,10)'     max. sent messages:',maxprocssend_fluid
        write(6,10)'     max. sent size    :',sizemsgsendmax_fluid
        write(6,10)'     nax. recv messages:',maxprocsrecv_fluid
        write(6,10)'     max. recv size    :',sizemsgrecvmax_fluid
     endif

     write(6,*)'  Source information___________________________________'
     write(6,16)'     Source type       :',src_type(1),src_type(2)
     write(6,11)'     Source depth   [m]:',zsrc
     write(6,11)'     Source colat [deg]:',srccolat*180./pi
     write(6,11)'     Source long  [deg]:',srclon*180./pi
     write(6,11)'     Magnitude    [N/m]:',magnitude
     write(6,12)'     Source time fct   :',stf_type
     write(6,11)'     Dom. period    [s]:',t_0
     write(6,*)'  Receiver information___________________________________'
     write(6,12)'     Receiver file type',rec_file_type
     write(6,19)'     Rotate to azimuth:',correct_azi
     write(6,19)'     Sum seismograms  :',sum_seis
     write(6,12)'     Components in:',rot_rec
     write(6,*)'  General numerical parameters_________________________'
     write(6,11)'     # elems/wavelength:',pts_wavelngth
     write(6,11)'     Courant number    :',courant
     write(6,11)'     Time step [s]     :',deltat
     write(6,10)'     # iterations      :',niter
     write(6,11)'     seismo length [s] :',niter*deltat
     write(6,12)'     time extrapolation:',time_scheme
     write(6,*)'  Input/Output information_____________________________'
     write(6,12)'     Output data path  :',datapath
     write(6,12)'     Output info path  :',infopath
     write(6,19)'     Save big testfiles:',save_large_tests
     write(6,19)'     Sum wavefields:', sum_fields
     write(6,19)'     Dump energy       :',dump_energy
     write(6,18)'     Glob/solflu snaps :',dump_snaps_glob,dump_snaps_solflu
     if (dump_snaps_glob .or. dump_snaps_solflu) then
        write(6,11)'     snap interval [s] :',snap_dt
        write(6,10)'     # snaps           :',snap_it
     endif
     write(6,19)'     Dump wavefields   :',dump_wavefields
     if (dump_wavefields) then 
        write(6,12)'     Dumping type      :',dump_type
        write(6,11)'     dump interval [s] :',period/real(strain_samp)
        write(6,10)'     # wavefield dumps :',strain_it
     endif
     write(6,19)'     Need fluid displ. :',need_fluid_displ
     write(6,*)
     write(6,*)':::::::::::::::: END SIMULATION PARAMETERS::::::::::::::::::::'
     write(6,*); call flush(6)

! additionally write a header for the kernel software
     if (dump_wavefields) call create_kernel_header

! write generic simulation info file
     open(unit=55,file='simulation.info')
     write(55,23)bkgrdmodel,'background model'
     write(55,21)deltat,'time step [s]'
     write(55,22)niter,'number of time steps'
     write(55,23)src_type(1),'source type'
     write(55,23)src_type(2),'source type'
     write(55,23)stf_type,'source time function'
     write(55,23)src_file_type,'source file type'
     write(55,21)period,'dominant source period'
     write(55,21)src_depth/1000.,'source depth [km]'
     write(55,25)magnitude,'scalar source magnitude'
     write(55,22)num_rec_tot,'number of receivers'
     write(55,22)floor(real(niter)/real(seis_it)),'length of seismogram [time samples]'
     write(55,21)deltat*seis_it,'seismogram sampling [s]'
     write(55,24)correct_azi,'compute seismograms at correct azimuth?'
     if (dump_wavefields) then
        write(55,22)floor(real(niter)/real(strain_it)),'number of strain dumps'
        write(55,21)period/real(strain_samp),'strain dump sampling rate [s]'
     else
        write(55,22)0,'number of strain dumps'       
        write(55,21)0.,'strain dump sampling rate [s]' 
     endif
   if (dump_snaps_glob .or. dump_snaps_solflu) then
      write(55,22)floor(real(niter)/real(snap_it)),'number of snapshot dumps'
      write(55,21)deltat*real(snap_it),'snapshot dump sampling rate [s]'      
   else
      write(55,22)0,'number of snapshot dumps'
      write(55,21)0.,'snapshot dump sampling rate [s]'      
   endif
   write(55,23)rot_rec,'receiver components '
   write(55,22)ibeg,'  ibeg: beginning gll index for wavefield dumps'
   write(55,22)iend,'iend: end gll index for wavefield dumps'
     close(55)

     write(6,*)
     write(6,*)'  wrote general simulation info into "simulation.info"'

21 format(f20.3,a45)
22 format(i20,a45)
23 format(a20,a45)
24 format(l20,a45)
25 format(1pe15.5,a45)

  endif ! lpr

! output for each processor==============================================

! extract processor location
  mysmin=router; myzmin=router; mysmax=zero; myzmax=zero
  myrmin=router; mythetamin=10.*pi; myrmax=zero; mythetamax=zero
  do iel=1,nelem
     do ipol=0,npol
        do jpol=0,npol
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
           if (s<mysmin) mysmin=s; if (z<myzmin) myzmin=z
           if (r<myrmin) myrmin=r; if (theta<mythetamin) mythetamin=theta
           if (s>mysmax) mysmax=s; if (z>myzmax) myzmax=z
           if (r>myrmax) myrmax=r; if (theta>mythetamax) mythetamax=theta
        enddo
     enddo
  enddo

  write(69,*)
  write(69,15)'My rank, total procs    :',mynum,nproc
  write(69,17)'Min./max. s [m]         :', mysmin,mysmax
  write(69,17)'Min./max. z [m]         :', myzmin,myzmax
  write(69,17)'Min./max. r [m]         :', myrmin,myrmax
  write(69,17)'Min./max. theta [deg]   :', mythetamin*180./pi, &
       mythetamax*180./pi
  write(69,13)'Have axis               ?',have_axis
  if (have_axis) then 
     write(69,10)'Axial total elems       :',naxel
     write(69,10)'Axial solid elems       :',naxel_solid
     write(69,10)'Axial fluid elems       :',naxel_fluid
  endif

  write(69,13)'Have source             ?',have_src
  if (have_src) then
     write(69,11)'Depth asked for      [m]:',src_depth
     write(69,11)'Computed depth       [m]:',router- &
                                    zcoord(ipol_src,jpol_src,ielsolid(iel_src))
  endif
  write(69,13)'Have boundary els       ?',have_bdry_elem
  if (have_bdry_elem) then
     write(69,10)'# boundary elements     :',nel_bdry
  endif
  write(69,*)
  write(69,*)'Solid message passing_____________________________'
  write(69,10)' # recv messages        :',sizerecv_solid
  write(69,10)' Max. size recv messages:',sizemsgrecvmax_solid
  write(69,10)' # sent messages        :',sizesend_solid
  write(69,10)' Max. size sent messages:',sizemsgsendmax_solid

  if (have_fluid) then
     write(69,*)'Fluid message passing_____________________________'
     write(69,10)' # recv messages        :',sizerecv_fluid
     write(69,10)' Max. size recv messages:',sizemsgrecvmax_fluid
     write(69,10)' # sent messages        :',sizesend_fluid
     write(69,10)' Max. size sent messages:',sizemsgsendmax_fluid
  endif !have_fluid

  call flush(69)

10 format(a25,i14)
11 format(a25,1pe14.5)
12 format(a25,'   ',a18)
13 format(a25,L14)
15 format(a25,i14,i9)
16 format(a25,' ',a12,a10)
17 format(a25,2(1pe13.3))
18 format(a25,2(L14))
19 format(a25,L14)

! write post processing file==============================================

if (lpr) then
   write(6,*)'  Writing post processing input file: param_post_processing'
   write(6,*)'  ... mainly based on guessing from the current simulation, make sure to edit!'

   open(unit=8,file='param_sum_seis') 
   read(8,*)nsim1
   close(8)

   open(unit=9,file="param_post_processing")
   if (rot_rec/='cyl') then
      write(9,222)'.false.','rotate receivers?'
      write(9,222)"'"//trim(rot_rec)//"'",'receiver components: enz,sph,cyl,xyz,src'
   else
      write(9,222)'.true.','rotate receivers?'
      write(9,222)"'sph'",'receiver components: enz,sph,cyl,xyz,src'
   endif
   if (src_file_type=='cmtsolut' .and. nsim1>1) then
      write(9,222)'.true.','sum to full Mij'
   elseif (nsim1>1) then 
      write(9,222)'.true.','sum to full Mij'
   else 
      write(9,222)'.true.','sum to full Mij'
   endif
   Mij_char = ['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp']
   do i=1,6
      write(9,223)Mij(i),Mij_char(i)
   enddo
   
   if (stf_type=='dirac_0' ) then!.or. stf_type=='quheavi' .or. stf_type=='heavis') then
      write(9,223)period,'convolve period ( 0. if not to be convolved)'
   else
      write(9,223)0.0,'convolve period (0. if not convolved)'
   endif
   write(9,222)"'gauss_0'",'source time function type for convolution'
   write(9,223)srccolat,'Source colatitude'
   write(9,223)srclon,'Source longitude'
   write(9,221)dump_snaps_glob,'plot global snaps?'
   write(9,224)'disp','disp or velo seismograms'
   write(9,224)"'Data_Postprocessing'",'Directory for post processed data'
   close(9)
221 format(l25,a50)
222 format(a25,a50)
223 format(1pe25.3,a50)
224 format(a25,a50)
   write(6,*)'    ... wrote file param_post_processing'
endif

! write param_snaps ==============================================
if (dump_snaps_glob) then 
   if (lpr) then
      write(6,*)'  Writing param_snaps for wavefield visualization'
      open(unit=9,file="param_snaps")

      write(9,222)"0.","starting phi (right cross section) [deg]"
      write(9,222)"85.", "increment phi (ending,left cross section) [deg]"
      write(9,223)router/1000.,"top radius [km]"
      write(9,222)"3190.","bottom radius [km]"
      write(9,222)"60. ", "meridional colatitude [deg]"
      write(9,225)1,floor(real(niter)/real(snap_it)),1, "snap starting number,end number, skipping factor"
      write(9,222)".false.","consider meridional cross section?"
      write(9,222)".true.","consider top surface?"
      write(9,222)".true. ","consider bottom surface?"
      close(9)
      write(6,*)'    ... wrote file param_snaps'; call flush(6)
   endif
225 format(i8,i8,i8,a50)
endif

end subroutine write_parameters
!=============================================================================

!----------------------------------------------------------------------------- 
subroutine create_kernel_header
!
! Static header with some info for subsequent kernel computations.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    character(len=8) :: mydate
    character(len=10) :: mytime
    character(len=80) :: dbname2
    integer :: lfdbname
    integer :: ndumps

    ndumps=floor(real(niter)/real(strain_it))

    call date_and_time(mydate,mytime)
    dbname2='mesh_params_kernel.h'
    lfdbname=index(dbname2,' ')-1
    open(97,file=dbname2(1:lfdbname))
    write(97,10) nproc
    write(97,11) mydate(5:6),mydate(7:8),mydate(1:4),mytime(1:2),mytime(3:4)
    write(97,*)''
    write(97,29)
    write(97,12)'Background model     :',bkgrdmodel
    write(97,13)'Inner-core shear wave:',resolve_inner_shear
    write(97,14)'Dominant period [s]  :',period
    write(97,14)'Elements/wavelength  :',pts_wavelngth
    write(97,14)'Courant number       :',courant
    write(97,30)
    write(97,*)''
    write(97,9)'nt',niter,'number of time steps'
    write(97,19)'deltat',deltat,'time step'
    write(97,9)'ndumps',ndumps, 'total wavefield dumps'
    write(97,9)'strain_samp',int(strain_samp),'dumps per period'
    write(97,18)"src_type",src_type(1),'source type'
    write(97,18)"src_type2",src_type(2),'source type'
    write(97,28)"bkgrdmodel",bkgrdmodel,&
                                                      'background model'
    write(97,9)'ibeg',ibeg,'dumped starting GLL within element'
    write(97,9)'iend',iend,'dumped ending GLL within element'
    
    if (have_fluid) then 
     write(97,31)
    else
     write(97,32)
    end if
    write(97,*)''
    write(97,30)
    write(97,*)''
    close(97)
    write(6,*)
    write(6,*)'wrote parameters for kerner into ',dbname2(1:lfdbname)
    write(6,*)

9 format(' integer, parameter :: ',A12,' =',i10,'  ! ',A27)
19 format(' real, parameter    :: ',A12,' =',f10.5,'  ! ',A27)
18 format(' character(len=10), parameter    :: ',A12," ='",A10,"'  ! ",A27)
28 format(' character(len=100), parameter    :: ',A12," ='",A10,"'  ! ",A27)
31 format(' logical, parameter    :: have_fluid=.true.')
32 format(' logical, parameter    :: have_fluid=.false.')
10 format('! Proc ',i3,': Header for kernel information to run static kerner')
11 format('! created by the solver on ', &
            A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
29 format('!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::')
12 format('!  ',A23,A20)
13 format('!  ',A23,L10)
14 format('!  ',A23,1f10.4)
30 format('!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

  end subroutine create_kernel_header
!------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine check_parameters(hmaxglob,hminglob,curvel,linel,seminoel,semisoel, &
                       curvel_solid,linel_solid,seminoel_solid,semisoel_solid,&
                       curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid)

use data_comm

include 'mesh_params.h'

double precision, intent(in) :: hmaxglob,hminglob
integer, intent(in) :: curvel,linel,seminoel,semisoel
integer, intent(in) :: curvel_solid,linel_solid,seminoel_solid,semisoel_solid
integer, intent(in) :: curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid

write(6,*)procstrg,'Checking solid message-passing...'
  if (nproc==1 .and. psum_int(sizesend_solid)>0 ) then 
     write(6,*)'Problem: Have only one proc but want to send messages..'
     stop
  endif

  if (nproc==1 .and. psum_int(sizerecv_solid)>0 ) then 
     write(6,*)'Problem: Have only one proc but want to receive messages...'
     stop
  endif

  if (nproc>1 .and. psum_int(sizesend_solid)==0 ) then 
     write(6,*)'Problem: No proc is willing to send anything....'
     stop
  endif

  if (nproc>1 .and. psum_int(sizesend_solid)==0 ) then 
     write(6,*)'Problem: No proc is willing to receive anything....'
     stop
  endif

  if (psum_int(sizesend_solid)< nproc-1 ) then 
     write(6,*)'Problem: Some proc(s) not willing to send anything...'
     stop
  endif

  if (psum_int(sizerecv_solid)< nproc-1 ) then 
     write(6,*)'Problem: Some proc(s) not willing to receive anything...'
     stop
  endif

  if (have_fluid) then
     write(6,*)procstrg,'Checking fluid message-passing...'
     if (nproc==1 .and. psum_int(sizesend_fluid)>0 ) then 
        write(6,*)'Problem: Have only one proc but want to send messages..'
        stop
     endif

     if (nproc==1 .and. psum_int(sizerecv_fluid)>0 ) then 
        write(6,*)'Problem: Have only one proc but want to receive messages...'
        stop
     endif

     if (nproc>1 .and. psum_int(sizesend_fluid)==0 ) then 
        write(6,*)'Problem: No proc is willing to send anything....'
        stop
     endif

     if (nproc>1 .and. psum_int(sizesend_fluid)==0 ) then 
        write(6,*)'Problem: No proc is willing to receive anything....'
        stop
     endif

     if (psum_int(sizesend_fluid)< nproc-1 ) then 
        write(6,*)'Problem: Some proc(s) not willing to send anything...'
        stop
     endif

     if (psum_int(sizerecv_fluid)< nproc-1 ) then 
        write(6,*)'Problem: Some proc(s) not willing to receive anything...'
        stop
     endif

  endif !have_fluid

! Even more tests.............
! stop if difference between loaded and on-the-fly mesh larger than 1 METER....
  if ( (hmin_glob-hminglob)>1.) then
     write(6,*)
     write(6,*)mynum,'Problem with minimal global grid spacing!'
     write(6,*)mynum,'Value from mesher        :',hmin_glob
     write(6,*)mynum,'Value computed on-the-fly:',hminglob
     stop
  endif

! stop if difference between loaded and on-the-fly mesh larger than 1 METER....
  if ((hmax_glob-hmaxglob)>1.) then
     write(6,*)
     write(6,*)mynum,'Problem with maximal global grid spacing!'
     write(6,*)mynum,'Value from mesher        :',hmax_glob
     write(6,*)mynum,'Value computed on-the-fly:',hmaxglob
     stop
  endif

! stop if sum of respective element types do not sum up to nelem
  if (curvel+linel+seminoel+semisoel/=nelem) then 
     write(6,*)
     write(6,*)mynum,'Problem with number of assigned global element types!'
     write(6,*)mynum,'curved,lin,semi(N),semi(S):',curvel,linel,seminoel, &
                                                   semisoel
     write(6,*)mynum,'Total # elements          :',nelem
     stop
  endif

  if (curvel_solid+linel_solid+seminoel_solid+semisoel_solid/=nel_solid) then 
     write(6,*)
     write(6,*)mynum,'Problem with number of assigned solid element types!'
     write(6,*)mynum,'curved,lin,semi(N),semi(S):',curvel_solid,linel_solid, &
                                                  seminoel_solid,semisoel_solid
     write(6,*)mynum,'Total # elements          :',nel_solid
     stop
  endif

  if (curvel_fluid+linel_fluid+seminoel_fluid+semisoel_fluid/=nel_fluid) then 
     write(6,*)
     write(6,*)mynum,'Problem with number of assigned fluid element types!'
     write(6,*)mynum,'curved,lin,semi(N),semi(S):',curvel_fluid,linel_fluid, &
                                                  seminoel_fluid,semisoel_fluid
     write(6,*)mynum,'Total # elements          :',nel_fluid
     stop
  endif

end subroutine check_parameters
!=============================================================================

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!========================
end module parameters
!========================

