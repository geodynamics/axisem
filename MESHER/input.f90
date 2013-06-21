!===================
module input
!===================

    use data_grid
    use data_diag
    use data_coarse
    use data_bkgrdmodel
  
    implicit none
  
    public :: read_params
    private
  
    contains

!----------------------------------------------------------------------------
subroutine read_params
  
    use global_parameters
    use data_mesh
    use data_spec
  
    character(len=100)    :: keyword, keyvalue, line
    integer               :: iinparam_mesh = 500, ioerr
  
    keyword = ' '
    keyvalue = ' '
    line = ' ' 
 
    ! Default values
    bkgrdmodel = 'UNDEFINED'
    nproc_target = -1
    period = -1
    dump_mesh_vtk = .true.
    nc_init = 3
    resolve_inner_shear = .true.
    npol = 4
    pts_wavelngth = 1.5
    courant = 0.6
    router = 6.371e6
    dump_mesh_info_files = .false.
    dump_mesh_info_screen = .false.
    diagpath = 'Diags'
    lfdiag = index(diagpath,' ') - 1 


    write(6, '(A)', advance='no') 'Reading inparam_mesh...'
    open(unit=iinparam_mesh, file='./inparam_mesh', status='old', action='read',  iostat=ioerr)
    if (ioerr.ne.0) stop 'Check input file ''inparam_mesh''! Is it still there?' 
 
    do
        read(iinparam_mesh,fmt='(a100)',iostat=ioerr) line
        if (ioerr.lt.0) exit
        if (len(trim(line)).lt.1.or.line(1:1).eq.'#') cycle
       
        read(line,*) keyword, keyvalue 
        parameter_to_read : select case(trim(keyword))
        
        case('BACKGROUND_MODEL') 
            bkgrdmodel = keyvalue
            lfbkgrdmodel = index(bkgrdmodel,' ') - 1 

        case('EXT_MODEL')
            fnam_ext_model = keyvalue

        case('DOMINANT_PERIOD')
            read(keyvalue, *) period

        case('NCPU')
            read(keyvalue, *) nproc_target

        case('WRITE_VTK')
            read(keyvalue, *) dump_mesh_vtk

        case('COARSENING_LAYERS')
            read(keyvalue, *) nc_init

        case('IC_SHEAR_WAVE')
            read(keyvalue, *) resolve_inner_shear

        case('NPOL')
            read(keyvalue, *) npol

        case('EL_PER_LAMBDA')
            read(keyvalue, *) pts_wavelngth

        case('COURANT_NR')
            read(keyvalue, *) courant

        case('RADIUS') 
            read(keyvalue, *) router

        case('SAVE_MESH')
            read(keyvalue, *) dump_mesh_info_files

        case('VERBOSE')
            read(keyvalue, *) dump_mesh_info_screen

        end select parameter_to_read
    end do

    if (trim(bkgrdmodel).eq.'UNDEFINED') then
        write(6,20) 'BACKGROUND_MODEL' 
        stop
    end if
    if (nproc_target==-1) then
        write(6,20) 'NCPU'
        stop
    end if
    if (period==-1) then
        write(6,20) 'DOMINANT_PERIOD'
        stop
    end if
20  format('ERROR: Parameter ', A, ' not set in inparam_mesh')
    write(6,*) 'done'

  
    write(6,*) ''
    write(6,*) 'PREDEFINED MODEL/SIMULATION PARAMETERS'
    write(6,*) 'Background model                 : ',bkgrdmodel(1:lfbkgrdmodel)
    write(6,*) 'Resolve inner core shear wave    : ',resolve_inner_shear
    write(6,*) 'Dominant period [s]              : ',period
    write(6,*) 'Elements per dominant wavelength : ',pts_wavelngth
    write(6,*) 'Courant number                   : ',courant
    write(6,*) 'coarsening levels                : ',nc_init
    write(6,*) 'processors used in solver        : ',nproc_target
    write(6,*) 'outer radius [m]                 : ',router
    write(6,*) 'save mesh info files?            : ',dump_mesh_info_files
    write(6,*) 'print mesh info to screen?       : ',dump_mesh_info_screen
    write(6,*) 'path to dump output files        : ',trim(diagpath)
    write(6,*) 
    call flush(6)
  
  !MvD: This did not make to much sense
  !if (realkind==4) then 
  !!  smallval = smallval_sngl
  !  smallval = smallval_dble
  !elseif (realkind==8) then 
  !  smallval = smallval_dble
  !endif
  smallval = 1.E-10

end subroutine read_params
!----------------------------------------------------------------------------

!===================
end module input
!===================
