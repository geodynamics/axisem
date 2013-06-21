!> Read parameters for the general solver (i.e. NOT mesh, sources, receivers);
!! compute other parameters for the simulation;
!! write out summaries of all relevant simulation settings.
module parameters

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

    character(len=100)  :: hostname, username, svn_version
    character(len=100)  :: compiler, compilerversion 

    public :: open_local_param_file, readin_parameters
    public :: compute_numerical_parameters, write_parameters
    private

contains

!-----------------------------------------------------------------------------
!> Open processor-specific output files
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
!> Routine that reads in simulation parameters that are relevant at the 
!! stage of the solver, i.e. number of time steps, integration scheme, 
!! data paths, specification of wavefield dumping etc.
subroutine readin_parameters

  include 'mesh_params.h'
  character(len=100)  :: junk
  integer             :: i

  call read_inparam_basic
  
  call read_inparam_advanced

  call get_runinfo


! now pre-set. Most of these are to be considered in the post processing stage now.
  sum_seis = .false.
  sum_fields = .false.
  rot_rec = 'cyl'
  dump_snaps_solflu = .false.
  dump_type = 'fullfields'
  num_simul = 1
   
  ! netcdf format
  output_format = 'binary'
  if (use_netcdf) output_format='netcdf'
  
  call barrier
  if (lpr) then
     write(6,20)
     write(6,21) datapath, infopath, num_simul,  seislength_t, enforced_dt,  &
                 enforced_period, trim(src_file_type), rec_file_type, &
                 sum_seis, sum_fields, rot_rec, time_scheme, seis_dt,  &
                 dump_energy, dump_vtk, dump_snaps_solflu, dump_wavefields, &
                 dump_type, ibeg, iend, strain_samp, src_dump_type, make_homo, &
                 add_hetero, do_mesh_tests, output_format

20 format(/&
   08x,'///////////////////////////////////////////////////////////////',/  &
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

21 format(&
   08x,'// (1) Tarje Nissen-Meyer, F. A. Dahlen, A Fournier (2007)   //',/  &
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
   08x,'=============  I N P U T    P A R A M E T E R S ===============',/  &
   12x,'Data I/O path:                      ',a20,/                         &
   12x,'Info I/O path:                      ',a20,/                         &
   12x,'Number of source simulations:       ',i2,/                          &
   12x,'Simulation length [s]:              ',f9.3,/                        &
   12x,'Enforced time step [s]:             ',f7.3,/                        &
   12x,'Enforced source period [s]:         ',f7.3,/                        &
   12x,'Source file type:                   ',a,/                           &
   12x,'Receiver file type:                 ',a8,/                          &
   12x,'Sum seismograms?                    ',l2,/                          &
   12x,'Sum wavefields?                     ',l2,/                          &
   12x,'Receivers coordinates               ',a3,/                          &
   12x,'Time extrapolation scheme:          ',a8,/                          &
   12x,'Seismogram sampling rate [s]:       ',f7.3,/                        &
   12x,'Dump kin./pot. energy?              ',l2,/                          &
   12x,'Dump global snaps?                  ',l2,/                          &
   12x,'Dump solid/fluid snaps?             ',l2,/                          &
   12x,'Dump strain?                        ',l2,/                          &
   12x,'Wavefield dumping type:             ',a12,/                         &
   12x,'First GLL to save in strains:       ',i2,/                          &
   12x,'Last GLL to save in strains:        ',i2,/                          &
   12x,'Samples per period for strains:     ',f7.3,/                        &
   12x,'Source dumping type:                ',a4,/                          &
   12x,'Homogenize background model?        ',l2,/                          &
   12x,'Add heterogeneous region?           ',l2,/                          &
   12x,'Perform extensive mesh tests?       ',l2,/                          &
   12x,'Output format (seism., wavefields): ',a6,/                          &
   08x,'===============================================================')

     write(6,'(a/a)') 'Processor-specific output is written to: output_proc<PROC ID>.dat', &
                      'All potential error messages will appear here...'
  endif !lpr

  call check_basic_parameters

  ! Need to decide here since this boolean is needed in def_precomp_terms
  need_fluid_displ = .false.
  if (dump_vtk .or. dump_xdmf .or. dump_snaps_solflu .or. dump_energy .or. & 
     dump_wavefields .and. dump_type=='fullfields') then
     ! Need to add this for each new type of wavefield dumping method that 
     ! requires the fluid displacement/velocities
     need_fluid_displ = .true.
  endif

  ! define general small value

  if (realkind==4) then
      smallval = smallval_sngl
  elseif (realkind==8) then
      smallval = smallval_dble
  endif

  if (lpr .and. verbose > 1) write(6,*)'     small value is:',smallval

end subroutine readin_parameters
!=============================================================================

!-----------------------------------------------------------------------------
!> Read file inparam_basic
subroutine read_inparam_basic
    use data_mesh,   only: meshname
    integer             :: iinparam_basic=500, ioerr, nval
    character(len=256)  :: line
    character(len=256)  :: keyword, keyvalue
    character(len=16)   :: simtype

    ! Default values
    enforced_period = 0.0
    enforced_dt = 0.0
    seislength_t = 1800.0
    seis_dt = 0.0
    add_hetero = .false.
    dump_vtk = .false.
    snap_dt = 20.
    dump_xdmf = .false.
    verbose = 1
    use_netcdf = .false.
    rec_file_type = 'stations'

    ! These values have to be set
    simtype = 'undefined'
    src_file_type = 'undefined'
    meshname = 'undefined'

    keyword = ' '
    keyvalue = ' '
    
    if (lpr .and. verbose > 1) write(6,'(A)', advance='no') '    Reading inparam_basic...'
    open(unit=iinparam_basic, file='inparam_basic', status='old', action='read',  iostat=ioerr)
    if (ioerr /= 0) stop 'Check input file ''inparam_basic''! Is it still there?' 
 
    do
        read(iinparam_basic, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

        read(line,*) keyword, keyvalue 
    
        parameter_to_read : select case(trim(keyword))
        
        case('SIMULATION_TYPE') 
            read(keyvalue, *) simtype
            select case(trim(simtype))
            case('single')
                src_file_type = 'sourceparams'
            case('forces')
                print *, 'FIXME: Forces needs to be checked!'
                stop 2
            case('moment')
                src_file_type = 'cmtsolut'
            case default
                if(lpr) then
                    print *, 'SIMULATION_TYPE in inparam_basic has invalid value ', simtype
                end if
                stop 2
            end select

!        case('SOURCE_TYPE')
!            read(keyvalue, *) src_file_type

        case('RECFILE_TYPE')
            rec_file_type = keyvalue

        case('SEISMOGRAM_LENGTH')
            read(keyvalue, *) seislength_t

        case('MESHNAME')
            meshname = keyvalue

        case('LAT_HETEROGENEITY')
            read(keyvalue, *) add_hetero

        case('ATTENUATION')
            read(keyvalue,*) do_anel 

        case('SAVE_SNAPSHOTS')
            read(keyvalue, *) dump_snaps_glob

        case('VERBOSITY')
            read(keyvalue, *) verbose


        end select parameter_to_read

    end do
    if (lpr .and. verbose > 1) print *, 'done'
    
end subroutine
!=============================================================================

!-----------------------------------------------------------------------------
!> Read file inparam_advanced
subroutine read_inparam_advanced
    include 'mesh_params.h'

    integer     :: iinparam_advanced=500, ioerr, nval
    character(len=256)                   :: line
    character(len=256) :: keyword, keyvalue

    ! Default values
    time_scheme = 'newmark2'
    datapath = './Data'
    infopath = './Info'
    do_mesh_tests = .false.
    dump_vtk = .false.
    dump_wavefields = .false.
    strain_samp = 8
    src_dump_type = 'mask'
    ibeg = 1
    iend = 1
    dump_energy = .false.
    make_homo = .false.
    force_ani = .true.
    do_anel = .false.
    deflate_level = 5

    keyword = ' '
    keyvalue = ' '

    if (lpr .and. verbose > 1) write(6, '(A)', advance='no') '    Reading inparam_advanced...'
    open(unit=iinparam_advanced, file='inparam_advanced', status='old', action='read', iostat=ioerr)
    if (ioerr /= 0) stop 'Check input file ''inparam_advanced''! Is it still there?' 
 
    do
        read(iinparam_advanced, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
        read(line,*) keyword, keyvalue 
    
        parameter_to_read : select case(trim(keyword))

        case('SAMPLING_RATE')
            read(keyvalue, *) seis_dt 
            
        case('SOURCE_PERIOD')
            read(keyvalue, *) enforced_period

        case('TIME_STEP')
            read(keyvalue, *) enforced_dt

        case('TIME_SCHEME')
            read(keyvalue,*) time_scheme

        case('DATA_DIR')
            datapath = keyvalue

        case('INFO_DIR')
            infopath = keyvalue

        case('MESH_TEST')
            read(keyvalue,*) do_mesh_tests

        case('KERNEL_WAVEFIELDS')
            read(keyvalue,*) dump_wavefields

        case('KERNEL_SPP')
            read(keyvalue,*) strain_samp

        case('KERNEL_SOURCE')
            read(keyvalue,*) src_dump_type

        case('KERNEL_IBEG')
            read(keyvalue,*) ibeg

        case('KERNEL_IEND')
            read(keyvalue,*) iend
            iend = npol - iend

        case('SAVE_ENERGY')
            read(keyvalue,*) dump_energy

        case('HOMO_MODEL')
            read(keyvalue,*) make_homo

        case('HOMO_VP')
            read(keyvalue,*) vphomo 
            vphomo = vphomo * 1.e3

        case('HOMO_VS')
            read(keyvalue,*) vshomo 
            vshomo = vshomo * 1.e3

        case('HOMO_RHO')
            read(keyvalue,*) rhohomo 
            rhohomo = rhohomo * 1.e3

        case('DEFLATE_LEVEL')
            read(keyvalue,*) deflate_level

        case('FORCE_ANISO')
            read(keyvalue,*) force_ani

        case('SNAPSHOT_DT')
            read(keyvalue, *) snap_dt

        case('SNAPSHOTS_FORMAT')
            if (dump_snaps_glob) then 
              select case (trim(keyvalue))
                case('xdmf') 
                    dump_xdmf = .true.
                case('vtk')
                    dump_vtk = .true.
                case default 
                    write(6,*)'invalid value for snapshots format!'; stop
                end select
            endif			

        case('USE_NETCDF')
            read(keyvalue, *) use_netcdf

        end select parameter_to_read

    end do
    lfdata = index(datapath,' ')-1
    lfinfo = index(infopath,' ')-1
    if (lpr .and. verbose > 1) print *, 'done'

end subroutine
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> Getting information like code revision, username and hostname
subroutine get_runinfo
    integer  :: iget_runinfo = 500, ioerr

    hostname = 'UNKNOWN'
    username = 'UNKNOWN'
    svn_version = 'UNKNOWN'

    if (lpr .and. verbose > 1) write(6, '(A)', advance='no') '    Reading runinfo... '
    open(unit=iget_runinfo, file='runinfo', status='old', action='read',  iostat=ioerr)
    if (ioerr /= 0) then
        if (lpr .and. verbose > 1) &
            write(6,*) 'No file ''runinfo'' found, continuing without.'
    else
        read(iget_runinfo,*) svn_version
        read(iget_runinfo,*) username
        read(iget_runinfo,*) hostname 
        if(lpr .and. verbose > 1) print *, 'done'
    end if
    close(iget_runinfo)

#if defined(__GFORTRAN__)
    compiler = 'gfortran'
#define gfortranversion __VERSION__
    compilerversion = gfortranversion
#undef gfortranversion
#endif

#if defined(__INTEL_COMPILER)
    compiler = 'ifort'
#define ifortversion __INTEL_COMPILER
    write(compilerversion, *) ifortversion
#undef ifortversion
#endif

end subroutine
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> Checking the consistency of some of the input parameters
subroutine check_basic_parameters

    if (trim(src_file_type) == 'undefined') then
        if (lpr) write(6,200) 'SIMULATION_TYPE', 'inparam_basic'
        stop
    end if
!    if (trim(rec_file_type).eq.'undefined') then
!        if (lpr) write(6,200) 'RECFILE_TYPE', 'inparam_basic'
!        stop
!    end if
    if (trim(meshname) == 'undefined') then
        if (lpr) write(6,200) 'MESHNAME', 'inparam_basic'
        stop
    end if

200 format(A,' is not defined in input file ',A, ', but has to')


#ifndef unc
    if (use_netcdf) then
       write(6,*) 'ERROR: trying to use netcdf IO but axisem was compiled without netcdf'
       stop 2
    endif
#endif

  if (src_dump_type == 'anal') then
     write(6,*) 'ERROR: Analytical source wavefield dump not implemented YET!'
     stop 2
  endif

  if (realkind /= 4 .and. realkind /= 8) then
     if (lpr) then
        write(6,'(a/a/ai4/a)') 'PROBLEM with REAL data kind!', &
                               '... can only handle real kinds 4 or 8.', &
                               'real kind here:', realkind, &
                               'change parameter realkind in global_parameters.f90'
     endif
     stop
  endif

  if (strain_samp > 15) then
     if (lpr) then     
        write(6,*)
        write(6,*) "!!!!!! NOT GOING ANY FURTHER !!!!!!"
        write(6,*) "  It's just too much to save 10 frames of strain & velocity"
        write(6,*) "  per source period! Choose something reasonable."
     endif
     stop
  endif

  if (enforced_dt > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14) 'maximal time step', enforced_dt
     endif
  endif

  if (enforced_period > zero) then
     if (lpr) then     
        write(6,*)
        write(6,14) 'min. source period', enforced_period
     endif
  endif

14 format('  WARNING: Overriding',a19,' with:',f8.3,' seconds')

  if (dump_vtk .and. dump_snaps_solflu) then 
      if (lpr) then
         write(6,*)''
         write(6,*)" NOT dumping the same snapshots twice (global AND solid/fluid)"
         write(6,*)'...hence reverting to dumping global snaps only. Sorry.'
      end if
      dump_snaps_solflu = .false.
  endif

7 format(04x,a62)

  if (verbose.ne.0) then
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
  endif
  write(6,*)


end subroutine check_basic_parameters
!=============================================================================

!-----------------------------------------------------------------------------
!> Compute numerical parameters, like time step, snapshot frequency
subroutine compute_numerical_parameters
  use attenuation, only: dump_memory_vars
  include "mesh_params.h"

  double precision :: s,z,r,theta,s_max,dshift
  double precision :: dsaxis(0:npol-1,0:npol), dzaxis(0:npol-1) 
  double precision :: minds(nelem),maxds(nelem),mindz(nelem),maxdz(nelem)
  integer          :: ielem,ipol,jpol,i
  logical          :: found_shift

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

  ! Overwrite time step or source period if demanded by input
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
  if (enforced_period>zero .and. (trim(stf_type)/='dirac_0' .or. trim(stf_type)/='quheavi') ) then 
     if (enforced_period<period) then 
        if (lpr) then 
           write(6,*)
           write(6,*)'    WARNING: Period smaller than necessary by mesh!'
           write(6,*)'    A pulse of this (short) chosen half width will produce numerical '
           write(6,*)'    noise on this (coarse) mesh'
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
     if (trim(stf_type)/='dirac_0' .or. trim(stf_type)/='quheavi') then
        if (lpr) then
           write(6,*)
           write(6,*)'    Using period of the mesh:',period
        endif
        t_0=period
     else
        t_0=period ! Just for consistency
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
  if (lpr) then
     write(6,*)
     write(6,22)'    desired seismogram sampling:',seis_dt,' seconds'
  endif
  if (seis_dt > zero  .and. seis_dt >= deltat ) then
     if (period/seis_dt<10.) then
        if (lpr) write(6,*) ' ! ! !   W A R N I N G   ! ! !'
        if (lpr) write(6,*) 'seismogram sampling too coarse w.r.t. mesh resolution:',period,seis_dt
        if (lpr) write(6,*) '... changing seismogram sampling rate to ',period/10.
        seis_dt=period/10.
     endif
     seis_it=floor((seis_dt+smallval_dble)/deltat)
     seis_dt=deltat*seis_it
  elseif (seis_dt < deltat) then
     if (lpr) write(6,*) 'seismogram sampling cannot be smaller than time step...'
     if (lpr) write(6,*) '...changing it to the time step'
     seis_dt = deltat
     seis_it = 1
  else
     seis_it = 1
     seis_dt = deltat
  endif
  deltat_coarse=seis_dt

  nseismo = floor(real(niter)/real(seis_it))

  if (lpr) then
     write(6,22)'    offered seismogram sampling:',deltat*seis_it,' seconds'
     write(6,13)'    ...that is, every          :',seis_it,' timesteps'
     write(6,11)'    number of samples          :',nseismo
  endif
22 format(a33,f9.2,a10)

  ! snapshot output, convert from interval given in seconds to 
  ! incremental time steps
  if (dump_vtk .or. dump_xdmf .or. dump_snaps_solflu .or. dump_memory_vars) then
     snap_it=floor(snap_dt/deltat)
     open(unit=2900+mynum,file=datapath(1:lfdata)//'/snap_info.dat'//appmynum)
     nsnap = floor(real(niter)/real(snap_it))
     write(2900,*) nsnap 
     do ielem=1,nsnap
        write(2900+mynum,*)real(ielem)*snap_dt,ielem*snap_it
     enddo
     close(2900+mynum)
     if (lpr) then
        write(6,*)
        write(6,11)'    Number of snapshots        :',nsnap
        write(6,12)'    ...approximately every     :',snap_dt,'seconds'
        write(6,13)'    ...that is, every          :',snap_it,'timesteps'
     endif
11   format(a33,i8)
12   format(a33,f8.2,a10)
13   format(a33,i8,a10)

  endif

  ! Source time function
  if (lpr) then
     write(6,*)''
     write(6,*)'  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
     write(6,*)'  SOURCE TIME FUNCTION: ',trim(stf_type)
     if (discrete_dirac) write(6,*)'    discrete Dirac type: ',discrete_choice
  endif

  period_vs_discrete_halfwidth = 8.
  sampling_per_a = 1
  decay=3.5d0

  if (.not. dump_wavefields & 
      .and. ( trim(stf_type)=='dirac_0' .or. trim(stf_type)=='queavi') &
      .and. deltat_coarse>1.9*deltat ) then 
     period_vs_discrete_halfwidth = period/(2.*deltat_coarse)
     if (period_vs_discrete_halfwidth<15.) period_vs_discrete_halfwidth=15.
     if (lpr) then 
        write(6,*)'    No wavefield dump, but discrete Dirac due to seismogram downsampling'
        write(6,*)'    Set discrete Dirac half width to [s]:',period/period_vs_discrete_halfwidth
        write(6,*)'    ... i.e. this part in the mesh period:',period_vs_discrete_halfwidth
     endif
  endif

  if (trim(stf_type)=='dirac_0' .or. trim(stf_type)=='quheavi')  then 
     discrete_dirac=.true.
     if (dump_wavefields .or. deltat_coarse> 1.9*deltat) then
        discrete_choice='gaussi'
     else
        discrete_choice = '1dirac'
     endif
     if (20.*seis_dt > period ) then 
        if (lpr) then 
           write(6,*)'   +++++++++++++++++++ W A R N I N G +++++++++++++++++++++ '
           write(6,*)'   The sampling rate of seismograms is quite coarse given the'
           write(6,*)'   Dirac delta source time function. We suggest to use at least'
           write(6,*)'   20 points per period to ensure accurate results, i.e. seis_dt<=',period/20.
           write(6,*)'   +++++++++++++++ E N D  o f  W A R N I N G +++++++++++++ '
        endif
     endif
  else
     discrete_dirac=.false.
  endif

  if (dump_wavefields) then 
     ! define coarse time step for strain wavefield dumps
     if (discrete_dirac) strain_samp = ceiling(period_vs_discrete_halfwidth*sampling_per_a)
     strain_it = floor(period/strain_samp / deltat)
     deltat_coarse=max(deltat_coarse,deltat*strain_it)
     strain_samp = period/deltat_coarse
     if (lpr) then
       write(6,*)'   dumping wavefields at sampling rate and deltat:',strain_samp,deltat_coarse
     end if
     strain_it=floor(t_0/real(strain_samp)/deltat)
     deltat_coarse=deltat*strain_it
  else
     strain_it=seis_it
  endif
  if (lpr) write(6,*)'   coarsest dump every', strain_it, 'th time step, dt:', deltat_coarse

  if (discrete_dirac) then
     discrete_dirac_halfwidth=period/period_vs_discrete_halfwidth
     t_0=discrete_dirac_halfwidth
     if (lpr) then 
        write(6,*)
        write(6,*)'   DISCRETE DIRAC DEFINITIONS:'
        write(6,*)'    Period discrete Dirac, mesh, simul. [s]:',&
                       real(discrete_dirac_halfwidth),real(period),real(t_0)
        write(6,*)'    period mesh/period discrete Dirac:',real(period_vs_discrete_halfwidth)
        write(6,*)"    deltat SEM, seis, coarse [s]:",real(deltat),real(seis_dt),real(deltat_coarse)
        write(6,*)'    # seismogram points per mesh,Dirac period:',&
                       real(period/seis_dt),real(discrete_dirac_halfwidth/seis_dt)
        if (dump_wavefields) then 
           write(6,*)'    # coarse points per mesh, Dirac period (int) :',strain_samp,sampling_per_a
           write(6,*)'    # coarse points per mesh, Dirac period (calc):',& 
                real(period/deltat_coarse),real(discrete_dirac_halfwidth/deltat_coarse)
        endif
        write(6,*)'    # SEM points per mesh, Dirac period:',&
                       real(period/deltat),real(discrete_dirac_halfwidth/deltat)
     endif

     ! parameters for source-time function time shift
     found_shift=.false.
     do i=1,ceiling(4.*discrete_dirac_halfwidth/deltat)
        dshift = deltat*ceiling(4.*discrete_dirac_halfwidth/deltat) + real(i)*deltat
        if ( .not. found_shift &
             .and. abs(nint(dshift/deltat_coarse)-dshift/deltat_coarse)<0.01*deltat &
             .and. abs(nint(dshift/deltat)-dshift/deltat)<0.01*deltat &
             .and. abs(nint(dshift/seis_dt)-dshift/seis_dt)<0.01*deltat) then 
           shift_fact_discrete_dirac = deltat_coarse*ceiling(dshift/deltat_coarse)
           found_shift=.true.
        endif
     enddo
     shift_fact=shift_fact_discrete_dirac

  else ! smooth source time function, e.g. Gauss
     shift_fact=deltat_coarse*ceiling(1.5*t_0/deltat_coarse)
  endif

  if (lpr) then 
     write(6,*)''
     write(6,*)'  SHIFT FACTOR of source time function [s]:',shift_fact
     write(6,*)'   # SEM, seis, coarse points per shift factor:',&
          real(shift_fact/deltat),real(shift_fact/seis_dt),real(shift_fact/deltat_coarse)
     write(6,*)'   # simul. half widths per shift factor:',real(shift_fact/t_0)
     if (discrete_dirac) write(6,*)'   # mesh halfwidths per shift fact',real(shift_fact/period)
     write(6,*)'  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
     write(6,*)''
  endif

! strain tensor output, convert from num of dumps per period into 
! incremental time steps
  if (dump_wavefields) then
     nstrain = floor(real(niter)/real(strain_it))

     open(unit=2900+mynum,file=datapath(1:lfdata)//'/strain_info.dat'//appmynum)
     write(2900,*) nstrain 
     do ielem= 1, nstrain 
        write(2900+mynum,*)real(ielem)*t_0/real(strain_samp),ielem*strain_it
     enddo
     close(2900+mynum)
     if (lpr) then
        write(6,*)
        write(6,11)'    Number of wavefield dumps  :', nstrain
             
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

  ! Initialize counters for I/O
  istrain = 0
  isnap = 0
  iseismo=0

  ! mesh info: coordinates of elements and collocation points               
  open(2222+mynum,file=infopath(1:lfinfo)//'/axial_points.dat'//appmynum)
  open(3333+mynum,file=infopath(1:lfinfo)//'/axial_ds_dz.dat'//appmynum)
  s_max=zero

  ! Set some parameters for faster access in time loop
  half_dt=half*deltat
  half_dt_sq=half*deltat**2

122 format(i8,i3,4(1pe15.5))

  do ielem = 1, nelem

     ! write out axial points
     if (axis(ielem)) then
        call compute_coordinates(s,z,r,theta,ielem,0,npol)
        do jpol=0,npol
           call compute_coordinates(s,z,r,theta,ielem,0,jpol)
           write(2222+mynum,122)ielem,jpol,s,z,r,theta/pi*180.
        enddo

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
  enddo
123 format(i7,1pe14.3,5(1pe14.3))
  close(2222+mynum)
  close(3333+mynum)

  if (lpr) write(6,*)

end subroutine compute_numerical_parameters
!=============================================================================

!-----------------------------------------------------------------------------
!> Writes out relevant simulation parameters, to simulation.info and the NetCDF
!! Output file.
subroutine write_parameters

    use data_comm
    use nc_routines
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
    character(len=7) :: clogic

    write(69,*)'  writing out all relevant simulation parameters...'

    write(69,*)'  number of respective element types...'

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

    write(69,*)'  grid spacing min/max...'
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

    write(69,*)'  calculating hmax...'
    hmax=max(maxval(dis1),maxval(dis2))
    write(69,*)'  hmaxstuff:',minval(dis1),minval(dis2),hmax
    hmaxglob=pmax(hmax)
    write(69,*)'  hmaxglob:',hmaxglob
    hmaxloc1=maxloc(dis1)
    if (maxval(dis2)<maxval(dis1)) hmaxloc1=maxloc(dis2)
    call compute_coordinates(s,z,rmaxglob,thetamaxglob,hmaxloc1(3), &
         hmaxloc1(1)-1,hmaxloc1(2)-1)
    write(69,*)' rmax,thetamax:',rmaxglob,thetamaxglob

    write(69,*)'  calculating hmin...'
    hmin=min(minval(dis1),minval(dis2))
    write(69,*)'  hminstuff:',minval(dis1),minval(dis2),hmin
    hminglob=pmin(hmin)
    write(69,*)'  hminglob:',hminglob

    hminloc1=minloc(dis1)
    if (minval(dis2)<minval(dis1)) hminloc1=minloc(dis2)
    call compute_coordinates(s,z,rminglob,thetaminglob,hminloc1(3), &
            hminloc1(1)-1,hminloc1(2)-1)

! Checking potential issues with input parameter consistency
    if (lpr) write(6,*)
    if (lpr) write(6,*)'  checking input parameters for consistency...'
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

        write(6,*)'  Code information_____________________________________'
        write(6,12)'     svn revision      :', svn_version
        write(6,12)'     username          :', username
        write(6,12)'     hostname          :', hostname
        write(6,12)'     compiler          :', compiler
        write(6,12)'     compilerversion   :', compilerversion
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
        write(6,12)'     Source time fct   :',trim(stf_type)
        write(6,11)'     Dom. period    [s]:',t_0
        write(6,*)'  Receiver information___________________________________'
        write(6,12)'     Receiver file type',rec_file_type
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
        write(6,12)'     Output data path  :',trim(datapath)
        write(6,12)'     Output info path  :',trim(infopath)
        write(6,19)'     Sum wavefields:', sum_fields
        write(6,19)'     Dump energy       :',dump_energy
        write(6,18)'     Glob/solflu snaps :',dump_vtk,dump_snaps_solflu
        write(6,18)'     XDMF VTK          :', dump_xdmf
        if (dump_vtk .or. dump_xdmf .or. dump_snaps_solflu) then
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
        write(6,*)
        call flush(6)

! additionally write a header for the kernel software
        if (dump_wavefields) call create_kernel_header



! write generic simulation info file
        open(unit=55,file='simulation.info')
        write(55,23)trim(bkgrdmodel),'background model'
        write(55,21)deltat,'time step [s]'
        write(55,22)niter,'number of time steps'
        write(55,23)trim(src_type(1)),'source type'
        write(55,23)trim(src_type(2)),'source type'
        write(55,23)trim(stf_type),'source time function'
        write(55,23)trim(src_file_type),'source file type'
        write(55,21)period,'dominant source period'
        write(55,21)src_depth/1000.,'source depth [km]'
        write(55,25)magnitude,'scalar source magnitude'
        write(55,22)num_rec_tot,'number of receivers'
        write(55,22)nseismo,'length of seismogram [time samples]'
        write(55,21)real(deltat)*real(seis_it),'seismogram sampling [s]'
        if (dump_wavefields) then
            write(55,22) nstrain,'number of strain dumps'
            write(55,21) real(period)/real(strain_samp),'strain dump sampling rate [s]'
        else
            write(55,22)0,'number of strain dumps'       
            write(55,21)0.,'strain dump sampling rate [s]' 
        endif
        if (dump_vtk .or. dump_xdmf .or. dump_snaps_solflu) then
            write(55,22) nsnap,'number of snapshot dumps'
            write(55,21)deltat*real(snap_it),'snapshot dump sampling rate [s]'      
        else
            write(55,22)0,'number of snapshot dumps'
            write(55,21)0.,'snapshot dump sampling rate [s]'      
        endif
        write(55,23)rot_rec,'receiver components '
        write(55,22)ibeg,'  ibeg: beginning gll index for wavefield dumps'
        write(55,22)iend,'iend: end gll index for wavefield dumps'
        write(55,21)shift_fact,'source shift factor [s]'
        write(55,22)int(shift_fact/deltat),'source shift factor for deltat'
        write(55,22)int(shift_fact/seis_dt),'source shift factor for seis_dt'
        write(55,22)int(shift_fact/deltat_coarse),'source shift factor for deltat_coarse'
        write(55,23)trim(rec_file_type),'receiver file type'
        write(55,21)dtheta_rec,'receiver spacing (0 if not even)'
        write(55,24)use_netcdf,'use netcdf for wavefield output?'
        close(55) ! simulation.info

        write(6,*)
        write(6,*)'  wrote general simulation info into "simulation.info"'

21      format(f22.7,a45)
22      format(i20,a45)
23      format(a20,a45)
24      format(l20,a45)
25      format(1pe15.5,a45)

    endif ! lpr

    if ((mynum.eq.0).and.(use_netcdf)) then !Only proc0 has the netcdf file open at that point
! write generic simulation info file
        write(6,*) ' Writing simulation info to netcdf file attributes' 
        call nc_write_att_char(trim(bkgrdmodel),'background model')
        call nc_write_att_char(trim(svn_version), 'SVN revision')
        call nc_write_att_char(trim(username),'user name')
        call nc_write_att_char(trim(hostname),'host name')
        call nc_write_att_char(trim(compiler),'compiler brand')
        call nc_write_att_char(trim(hostname),'compiler version')
        call nc_write_att_real(real(deltat),'time step in sec')
        call nc_write_att_int(niter,'number of time steps')
        call nc_write_att_char(trim(src_type(1)),'excitation type')
        call nc_write_att_char(trim(src_type(2)),'source type')
        call nc_write_att_char(trim(stf_type),'source time function')
        call nc_write_att_char(trim(src_file_type),'source file type')
        call nc_write_att_real(real(period),'dominant source period')
        call nc_write_att_real(real(src_depth/1000.),'source depth in km')
        call nc_write_att_real(real(magnitude),'scalar source magnitude')
        call nc_write_att_int(num_rec_tot,'number of receivers')
        call nc_write_att_int(nseismo,'length of seismogram  in time samples')
        call nc_write_att_real(real(deltat)*real(seis_it),'seismogram sampling in sec')
        if (dump_wavefields) then
           call nc_write_att_int(nstrain,'number of strain dumps')
           call nc_write_att_real(real(period)/real(strain_samp),'strain dump sampling rate in sec')
        else
           call nc_write_att_int(0,'number of strain dumps')       
           call nc_write_att_real(0.,'strain dump sampling rate in sec' )
        endif
        if (dump_vtk .or. dump_snaps_solflu) then
           call nc_write_att_int(nsnap,'number of snapshot dumps')
           call nc_write_att_real(real(deltat)*real(snap_it),'snapshot dump sampling rate in sec')      
        else
           call nc_write_att_int(0,'number of snapshot dumps')
           call nc_write_att_real(0.,'snapshot dump sampling rate in sec')
        endif
        call nc_write_att_char(rot_rec,'receiver components ')
        call nc_write_att_int(ibeg,'ibeg')
        call nc_write_att_int(iend,'iend')
        call nc_write_att_real(shift_fact,'source shift factor in sec')
        call nc_write_att_int(int(shift_fact/deltat),'source shift factor for deltat')
        call nc_write_att_int(int(shift_fact/seis_dt),'source shift factor for seis_dt')
        call nc_write_att_int(int(shift_fact/deltat_coarse),'source shift factor for deltat_coarse')
        call nc_write_att_char(trim(rec_file_type),'receiver file type')
        call nc_write_att_real(dtheta_rec,'receiver spacing (0 if not even)')
        write(clogic,*) use_netcdf
        call nc_write_att_char(clogic,'use netcdf for wavefield output?')
    end if


! output for each processor==============================================

! extract processor location
    mysmin=router; myzmin=router; mysmax=zero; myzmax=zero
    myrmin=router; mythetamin=10.*pi; myrmax=zero; mythetamax=zero
    do iel=1,nelem
        do ipol=0,npol
            do jpol=0,npol
                call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
                if (s<mysmin) mysmin=s
                if (s>mysmax) mysmax=s
                if (r<myrmin) myrmin=r
                if (z<myzmin) myzmin=z
                if (z>myzmax) myzmax=z
                if (r>myrmax) myrmax=r
                if (theta<mythetamin) mythetamin=theta
                if (theta>mythetamax) mythetamax=theta
            enddo
        enddo
    enddo
  
    write(69,*)
    write(69,15)'My rank, total procs    :', mynum,nproc
    write(69,17)'Min./max. s [m]         :', mysmin,mysmax
    write(69,17)'Min./max. z [m]         :', myzmin,myzmax
    write(69,17)'Min./max. r [m]         :', myrmin,myrmax
    write(69,17)'Min./max. theta [deg]   :', mythetamin*180./pi, mythetamax*180./pi
  
    write(69,10)'Axial total elems       :', naxel
    write(69,10)'Axial solid elems       :', naxel_solid
    write(69,10)'Axial fluid elems       :', naxel_fluid
  
    write(69,13)'Have source             ?', have_src
    if (have_src) then
        write(69,11)'Depth asked for      [m]:', src_depth
        write(69,11)'Computed depth       [m]:', router - &
                                      zcoord(ipol_src,jpol_src,ielsolid(iel_src))
    endif
    write(69,13)'Have boundary els       ?', have_bdry_elem
    if (have_bdry_elem) then
        write(69,10)'# boundary elements     :', nel_bdry
    end if

    write(69,*)
    write(69,*)'Solid message passing_____________________________'
    write(69,10)' # recv messages        :', sizerecv_solid
    write(69,10)' Max. size recv messages:', sizemsgrecvmax_solid
    write(69,10)' # sent messages        :', sizesend_solid
    write(69,10)' Max. size sent messages:', sizemsgsendmax_solid
  
    if (have_fluid) then
        write(69,*)'Fluid message passing_____________________________'
        write(69,10)' # recv messages        :', sizerecv_fluid
        write(69,10)' Max. size recv messages:', sizemsgrecvmax_fluid
        write(69,10)' # sent messages        :', sizesend_fluid
        write(69,10)' Max. size sent messages:', sizemsgsendmax_fluid
    endif !have_fluid

   call flush(69)

10  format(a25,i14)
11  format(a25,1pe14.5)
12  format(a25,'   ',a18)
13  format(a25,L14)
15  format(a25,i14,i9)
16  format(a25,' ',a12,a10)
17  format(a25,2(1pe13.3))
18  format(a25,2(L14))
19  format(a25,L14)

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
        write(9,222)"'enz'",'receiver components: enz,sph,cyl,xyz,src'
    endif
    if (trim(src_file_type)=='cmtsolut' .and. nsim1>1) then
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
    
    if (stf_type=='dirac_0' .or. stf_type=='quheavi' ) then
        write(9,223)period,'convolve period ( 0. if not to be convolved)'
    else
        write(9,223)0.0,'convolve period (0. if not convolved)'
    endif
    write(9,222)"'gauss_0'",'source time function type for convolution'
    write(9,223)srccolat,'Source colatitude'
    write(9,223)srclon,'Source longitude'
    write(9,221)dump_vtk,'plot global snaps?'
    write(9,224)'disp','disp or velo seismograms'
    write(9,224)"'Data_Postprocessing'",'Directory for post processed data'
    write(9,221).true.,'seismograms at negative time (0 at max. of stf)'
    close(9)
221 format(l25,a50)
222 format(a25,a50)
223 format(1pe25.6,a50)
224 format(a25,a50)
    write(6,*)'    ... wrote file param_post_processing'
endif

! write param_snaps ==============================================
if (dump_vtk) then 
    if (lpr) then
        write(6,*)'  Writing param_snaps for wavefield visualization'
        open(unit=9,file="param_snaps")

        write(9,222)"0.","starting phi (right cross section) [deg]"
        write(9,222)"85.", "increment phi (ending,left cross section) [deg]"
        write(9,223)router/1000.,"top radius [km]"
        write(9,222)"3190.","bottom radius [km]"
        write(9,222)"60. ", "meridional colatitude [deg]"
        write(9,225)1,nsnap,1, "snap starting number,end number, skipping factor"
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
!> Static header with some info for subsequent kernel computations.
!! Produces file mesh_params_kernel.h
subroutine create_kernel_header

    use data_mesh, ONLY: hmax_glob
    use data_io, ONLY: nstrain
    character(len=8)  :: mydate
    character(len=10) :: mytime
    character(len=80) :: dbname2
    integer           :: lfdbname

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
    write(97,17)'hmax', hmax_glob,'maximum element size'
    write(97,17)'vpmax', vpmax,'maximum p-velocity'
    write(97,17)'rmax', router, 'maximum radius'
    write(97,9)'ndumps',nstrain, 'total wavefield dumps'
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

9   format(' integer, parameter :: ',A12,' =',i11,'  ! ',A27)
17  format(' real, parameter    :: ',A12,' =',f11.2,'  ! ',A27)
19  format(' real, parameter    :: ',A12,' =',f11.5,'  ! ',A27)
18  format(' character(len=10), parameter    :: ',A12," ='",A10,"'  ! ",A27)
28  format(' character(len=100), parameter    :: ',A12," ='",A10,"'  ! ",A27)
31  format(' logical, parameter    :: have_fluid=.true.')
32  format(' logical, parameter    :: have_fluid=.false.')
10  format('! Proc ',i3,': Header for kernel information to run static kerner')
11  format('! created by the solver on ', &
             A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
29  format('!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::')
12  format('!  ',A23,A20)
13  format('!  ',A23,L10)
14  format('!  ',A23,1f10.4)
30  format('!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

end subroutine create_kernel_header
!------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!< Checking some mesh parameters and message parsing
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

!========================
end module parameters
!========================

