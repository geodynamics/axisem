!=========================
module wavefields_io
!=========================
!
! Contains all routines that dump entire wavefields during the time loop. 
! Optimization of I/O therefore happens here and nowhere else.
! The corresponding meshes are dumped in meshes_io.
!
use global_parameters
use data_mesh
use data_proc
use data_io
use nc_routines

implicit none

public 

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
subroutine glob_snapshot(f_sol,chi,ibeg,iend,jbeg,jend)
!
! Dumps the global displacement snapshots [m] in ASCII format
! When reading the fluid wavefield, one needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord
! as dumped by the corresponding routine dump_glob_grid!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f_sol(0:npol,0:npol,1:nel_solid,3)
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,iidim
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call define_io_appendix(appisnap,isnap)

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/snap_'&
                            //appmynum//'_'//appisnap//'.dat')

  if (have_fluid) then
     call axisym_laplacian_fluid(chi,usz_fluid)
     do iel=1,nel_fluid

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,0,dsdchi)
                 write(2500+mynum,*)usz_fluid(0,0,iel,1),prefac*dsdchi*&
                                  chi(0,0,iel),usz_fluid(0,0,iel,2)
              else
                 write(2500+mynum,*)usz_fluid(0,0,iel,1), &
                                    prefac*chi(0,0,iel), &
                                    usz_fluid(0,0,iel,2)
              endif

                 write(2500+mynum,*)usz_fluid(npol,0,iel,1), &
                                    prefac*chi(npol,0,iel), &
                                    usz_fluid(npol,0,iel,2)

                 write(2500+mynum,*)usz_fluid(npol,npol,iel,1), &
                                    prefac*chi(npol,npol,iel), &
                                    usz_fluid(npol,npol,iel,2)

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,npol,dsdchi)
                 write(2500+mynum,*)usz_fluid(0,npol,iel,1),prefac*dsdchi*&
                                  chi(0,npol,iel),usz_fluid(0,npol,iel,2)
              else
                 write(2500+mynum,*)usz_fluid(0,npol,iel,1), &
                                    prefac*chi(0,npol,iel), &
                                    usz_fluid(0,npol,iel,2)
              endif

     enddo
  endif ! have_fluid

  do iel=1,nel_solid
           write(2500+mynum,*)(f_sol(0,0,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(npol,0,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(npol,npol,iel,iidim),iidim=1,3)
           write(2500+mynum,*)(f_sol(0,npol,iel,iidim),iidim=1,3)
  enddo

  close(2500+mynum)


!  h_real=real(hmax/(period/(pts_wavelngth*real(npol))))
!  fname=trim(diagpath)//'/mesh_hmax'
!  call write_VTK_bin_scal(h_real,mesh2,neltot,fname)


end subroutine glob_snapshot
!=============================================================================

!-----------------------------------------------------------------------------
subroutine glob_snapshot_midpoint(f_sol,chi,ibeg,iend,jbeg,jend)
!
! Dumps the global displacement snapshots [m] in ASCII format
! When reading the fluid wavefield, one needs to multiply all 
! components with inv_rho_fluid and the phi component with one/scoord
! as dumped by the corresponding routine dump_glob_grid!
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f_sol(0:npol,0:npol,1:nel_solid,3)
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,iidim
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call define_io_appendix(appisnap,isnap)

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/snap_'&
                            //appmynum//'_'//appisnap//'.dat' ,FORM="UNFORMATTED",STATUS="REPLACE")

  if (have_fluid) then
     call axisym_laplacian_fluid(chi,usz_fluid)
     do iel=1,nel_fluid
        do jpol=0,npol,npol/2
           do ipol=0,npol,npol/2

              if ( axis_fluid(iel)) then
                 call dsdf_fluid_axis(chi(:,:,iel),iel,jpol,dsdchi)
                 write(2500+mynum)usz_fluid(ipol,jpol,iel,1),prefac*dsdchi*&
                                  chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
              else
                 write(2500+mynum)usz_fluid(ipol,jpol,iel,1), &
                                    prefac*chi(ipol,jpol,iel), &
                                    usz_fluid(ipol,jpol,iel,2)
              endif
              enddo
           enddo
     enddo
  endif ! have_fluid

  do iel=1,nel_solid
     do jpol=0,npol,npol/2
        do ipol=0,npol,npol/2
           write(2500+mynum)(f_sol(ipol,jpol,iel,iidim),iidim=1,3)
        enddo
     enddo
  enddo
  close(2500+mynum)

!  h_real=real(hmax/(period/(pts_wavelngth*real(npol))))
!  fname=trim(diagpath)//'/mesh_hmax'
!  call write_VTK_bin_scal(h_real,mesh2,neltot,fname)


end subroutine glob_snapshot_midpoint
!=============================================================================


!-----------------------------------------------------------------------------
subroutine glob_snapshot_xdmf(f_sol, chi)

    use data_source, ONLY : src_type
    use data_pointwise, ONLY : inv_rho_fluid, inv_s_rho_fluid
    use pointwise_derivatives, ONLY: axisym_laplacian_fluid, dsdf_fluid_axis
    use data_time, only : t
    
    include 'mesh_params.h'
    
    real(kind=realkind), intent(in) :: f_sol(0:npol,0:npol,1:nel_solid,3)
    real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
    character(len=4)     :: appisnap
    integer              :: iel, iidim, ct, ipol, jpol, ipol1, jpol1, i, j
    real(kind=realkind)  :: dsdchi, prefac
    real*4, allocatable  :: u(:,:), u_new(:,:), usz_fl(:,:,:,:), up_fl(:,:,:)
    character(len=120)   :: fname

    allocate(usz_fl(0:npol,0:npol,1:nel_fluid,2))
    allocate(up_fl(0:npol,0:npol,1:nel_fluid))
    
    allocate(u(1:3,1:npoint_plot))

    if (src_type(1)=='monopole') prefac = 0.
    if (src_type(1)=='dipole')   prefac = 1.
    if (src_type(1)=='quadpole') prefac = 2.
 
    call define_io_appendix(appisnap,isnap)
 
    if (have_fluid) then
       call axisym_laplacian_fluid(chi,usz_fl)
       usz_fl(:,:,:,1) = usz_fl(:,:,:,1) * inv_rho_fluid
       usz_fl(:,:,:,2) = usz_fl(:,:,:,2) * inv_rho_fluid

       up_fl(:,:,:) = prefac * chi * inv_s_rho_fluid

       do iel=1, nel_fluid
           do i=1, i_n_xdmf - 1
               ipol = i_arr_xdmf(i)
               ipol1 = i_arr_xdmf(i+1)

               do j=1, j_n_xdmf - 1
                   jpol = j_arr_xdmf(j)
                   jpol1 = j_arr_xdmf(j+1)
       
                   if (plotting_mask(i,j,iel)) then
                       ct = mapping_ijel_iplot(i,j,iel)
                       u(1, ct) = usz_fl(ipol,jpol,iel,1)
                       u(2, ct) =  up_fl(ipol,jpol,iel)
                       u(3, ct) = usz_fl(ipol,jpol,iel,2)
                   endif

                   if (plotting_mask(i+1,j,iel)) then
                       ct = mapping_ijel_iplot(i+1,j,iel)
                       u(1, ct) = usz_fl(ipol1,jpol,iel,1)
                       u(2, ct) =  up_fl(ipol1,jpol,iel)
                       u(3, ct) = usz_fl(ipol1,jpol,iel,2)
                   endif

                   if (plotting_mask(i+1,j+1,iel)) then
                       ct = mapping_ijel_iplot(i+1,j+1,iel)
                       u(1, ct) = usz_fl(ipol1,jpol1,iel,1)
                       u(2, ct) =  up_fl(ipol1,jpol1,iel)
                       u(3, ct) = usz_fl(ipol1,jpol1,iel,2)
                   endif

                   if (plotting_mask(i,j+1,iel)) then
                       ct = mapping_ijel_iplot(i,j+1,iel)
                       u(1, ct) = usz_fl(ipol,jpol1,iel,1)
                       u(2, ct) =  up_fl(ipol,jpol1,iel)
                       u(3, ct) = usz_fl(ipol,jpol1,iel,2)
                   endif
               enddo
           enddo
       enddo
    endif
    
    deallocate(usz_fl, up_fl)
    
    do iel=1, nel_solid
        do i=1, i_n_xdmf - 1
            ipol = i_arr_xdmf(i)
            ipol1 = i_arr_xdmf(i+1)

            do j=1, j_n_xdmf - 1
                jpol = j_arr_xdmf(j)
                jpol1 = j_arr_xdmf(j+1)
    
                if (plotting_mask(i,j,iel + nel_fluid)) then
                    ct = mapping_ijel_iplot(i,j,iel + nel_fluid)
                    u(1, ct) = f_sol(ipol,jpol,iel,1)
                    u(2, ct) = f_sol(ipol,jpol,iel,2)
                    u(3, ct) = f_sol(ipol,jpol,iel,3)
                endif
                
                if (plotting_mask(i+1,j,iel + nel_fluid)) then
                    ct = mapping_ijel_iplot(i+1,j,iel + nel_fluid)
                    u(1, ct) = f_sol(ipol1,jpol,iel,1)
                    u(2, ct) = f_sol(ipol1,jpol,iel,2)
                    u(3, ct) = f_sol(ipol1,jpol,iel,3)
                endif
                
                if (plotting_mask(i+1,j+1,iel + nel_fluid)) then
                    ct = mapping_ijel_iplot(i+1,j+1,iel + nel_fluid)
                    u(1, ct) = f_sol(ipol1,jpol1,iel,1)
                    u(2, ct) = f_sol(ipol1,jpol1,iel,2)
                    u(3, ct) = f_sol(ipol1,jpol1,iel,3)
                endif
                
                if (plotting_mask(i,j+1,iel + nel_fluid)) then
                    ct = mapping_ijel_iplot(i,j+1,iel + nel_fluid)
                    u(1, ct) = f_sol(ipol,jpol1,iel,1)
                    u(2, ct) = f_sol(ipol,jpol1,iel,2)
                    u(3, ct) = f_sol(ipol,jpol1,iel,3)
                endif
            enddo
        enddo
    enddo
      
    fname = datapath(1:lfdata)//'/xdmf_snap_s_' //appmynum//'.dat'
    open(100, file=trim(fname), access='stream', status='unknown', &
        convert='little_endian', position='append')
    write(100) u(1,:)
    close(100)

    if (.not. src_type(1)=='monopole') then
        fname = datapath(1:lfdata)//'/xdmf_snap_p_' //appmynum//'.dat'
        open(100, file=trim(fname), access='stream', status='unknown', &
            convert='little_endian', position='append')
        write(100) u(2,:)
        close(100)
    endif

    fname = datapath(1:lfdata)//'/xdmf_snap_z_' //appmynum//'.dat'
    open(100, file=trim(fname), access='stream', status='unknown', &
        convert='little_endian', position='append')
    write(100) u(3,:)
    close(100)

    deallocate(u)
    
    fname = datapath(1:lfdata) // '/xdmf_xml_' // appmynum // '.xdmf'
    open(100, file=trim(fname), access='append')

    if (src_type(1)=='monopole') then
        write(100, 734) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                    npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, 'xdmf_snap_s_'//appmynum//'.dat', &
                    npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, 'xdmf_snap_z_'//appmynum//'.dat', &
                    npoint_plot, appisnap, appisnap
    else
        write(100, 735) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                    npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, 'xdmf_snap_s_'//appmynum//'.dat', &
                    npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, 'xdmf_snap_p_'//appmynum//'.dat', &
                    npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, 'xdmf_snap_z_'//appmynum//'.dat', &
                    npoint_plot, appisnap, appisnap, appisnap
    endif

734 format(&    
    '    <Grid Name="', A,'" GridType="Uniform">',/&
    '        <Time Value="',F8.2,'" />',/&
    '        <Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '        </Topology>',/&
    '        <Geometry GeometryType="XY">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '        </Geometry>',/&
    '        <Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '        <Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '        <Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="Function" Function="sqrt($0 * $0 + $1 * $1)" Dimensions="', I10,'">',/&
    '                <DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '                </DataItem>',/&
    '                <DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '    </Grid>',/)

735 format(&    
    '    <Grid Name="', A,'" GridType="Uniform">',/&
    '        <Time Value="',F8.2,'" />',/&
    '        <Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '        </Topology>',/&
    '        <Geometry GeometryType="XY">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '        </Geometry>',/&
    '        <Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '        <Attribute Name="u_p" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '        <Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '        <Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="Function" Function="sqrt($0 * $0 + $1 * $1 + $2 * $2)" Dimensions="', I10,'">',/&
    '                <DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '                </DataItem>',/&
    '                <DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_p"]/DataItem[1]',/&
    '                </DataItem>',/&
    '                <DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>',/&
    '    </Grid>',/)
    
    close(100)

end subroutine glob_snapshot_xdmf
!=============================================================================


!-----------------------------------------------------------------------------
subroutine solid_snapshot(f,ibeg,iend,jbeg,jend)
!
! Dumps the displacement snapshots [m] in the solid region in ASCII format
! Convention for order in the file: First the fluid, then the solid domain.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: f(0:npol,0:npol,1:nel_solid,3)
character(len=4)                :: appisnap
integer                         :: iel, ipol,jpol,idim

  call define_io_appendix(appisnap,isnap)

  open(unit=3500+mynum,file=datapath(1:lfdata)//'/snap_solid_'&
                            //appmynum//'_'//appisnap//'.dat')

  do iel=1,nel_solid
     do jpol=ibeg,iend
        do ipol=jbeg,jend
           write(3500+mynum,*)(f(ipol,jpol,iel,idim),idim=1,3)
        enddo
     enddo
  enddo
  close(3500+mynum)

end subroutine solid_snapshot
!=============================================================================

!-----------------------------------------------------------------------------

subroutine fluid_snapshot(chi,ibeg,iend,jbeg,jend)

use data_pointwise, ONLY : usz_fluid
use data_source, ONLY : src_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_axis

include 'mesh_params.h'

integer, intent(in)             :: ibeg,iend,jbeg,jend
real(kind=realkind), intent(in) :: chi(0:npol,0:npol,1:nel_fluid)
character(len=4)                :: appisnap
integer                         :: iel,ipol,jpol
real(kind=realkind)             :: dsdchi,prefac

! When reading the fluid wavefield, one needs to multiply all components 
! with inv_rho_fluid and the phi component with one/scoord!!

  if (src_type(1)=='monopole') prefac=zero
  if (src_type(1)=='dipole')   prefac=one
  if (src_type(1)=='quadpole') prefac=two

  call axisym_laplacian_fluid(chi,usz_fluid)

  call define_io_appendix(appisnap,isnap)

  open(unit=4500+mynum,file=datapath(1:lfdata)//'/snap_fluid_'&
                            //appmynum//'_'//appisnap//'.dat')

  do iel=1,nel_fluid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           if ( axis_fluid(iel) .and. ipol==0 ) then
              call dsdf_fluid_axis(chi(:,:,iel),iel,jpol,dsdchi)
              write(4500+mynum,*)usz_fluid(ipol,jpol,iel,1),prefac*dsdchi* &
                                chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
           else
              write(4500+mynum,*)usz_fluid(ipol,jpol,iel,1),prefac*&
                                 chi(ipol,jpol,iel),usz_fluid(ipol,jpol,iel,2)
           endif
        enddo
     enddo
  enddo

  close(4500+mynum)

end subroutine fluid_snapshot
!=============================================================================



!--------------------------------------------------------------------------
subroutine dump_field_1d(f,filename,appisnap,n)

use data_proc, ONLY : appmynum
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

integer, intent(in) :: n
real(kind=realkind),intent(inout) :: f(0:npol,0:npol,1:n)
character(len=16), intent(in)     :: filename
character(len=4), intent(in)      :: appisnap
real(kind=realkind), allocatable :: f1(:),f2(:)
integer :: i,j,iel,ii, f1len

  if (have_src .and. src_dump_type == 'mask' .and. n==nel_solid) &
       call eradicate_src_elem_values(f)

  f1len = (iend-ibeg+1)**2*n
  allocate(f1(f1len))

!  ii=0
!  do iel=1,n
!     do j=ibeg,iend
!        do i=ibeg,iend
!           ii=ii+1
!           f1(ii)=f(i,j,iel)
!        enddo
!     enddo
!  enddo
  
!  f1 = reshape(f(ibeg:iend,ibeg:iend,1:n),(/f1len/))
  f1 = pack(f(ibeg:iend,ibeg:iend,1:n),.true.)
  
  if (use_netcdf) then
    call nc_dump_field_1d(f1,f1len,filename(2:),appisnap)
  
  else ! Binary
    open(unit=25000+mynum,file=datapath(1:lfdata)//filename//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="UNKNOWN",POSITION="REWIND")
    write(25000+mynum) f1
    close(25000+mynum)
  end if

  deallocate(f1)


end subroutine dump_field_1d
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_field_over_s_solid_1d(f,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_solid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_solid)
integer                        :: iel,i,glen
!real(kind=realkind) :: f1d((npol+1)**2*nel_solid)
real(kind=realkind)            :: g(0:npol,0:npol,nel_solid)



  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

  g = inv_s_solid*f

  if (have_src .and. src_dump_type == 'mask') then 
    call eradicate_src_elem_values(g)
  end if
  
  glen = size(g)
  if (use_netcdf) then
    call nc_dump_field_1d(g,glen,filename(2:),appisnap)
  else
    open(unit=35000+mynum,file=datapath(1:lfdata)//filename//'_'&
                            //appmynum//'_'//appisnap//'.bindat',&
                            FORM="UNFORMATTED",STATUS="REPLACE")
    write(35000+mynum) g(ibeg:iend,ibeg:iend,:)
    close(35000+mynum)
  end if

end subroutine dump_field_over_s_solid_1d
!=============================================================================




!--------------------------------------------------------------------------
subroutine dump_field_over_s_solid_and_add(f,g,filename1,filename2,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the solid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_solid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_solid)
character(len=16), intent(in)  :: filename1,filename2
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_solid)
integer                        :: iel,i,glen

  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

! construct masked f/s (e.g. Epp)
  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(inv_s_solid)

! construct sum of f/s and g (e.g. straintrace)
  g = inv_s_solid*f + g

  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(g)


  glen = size(g(ibeg:iend,ibeg:iend,:))
  if (use_netcdf) then
    call nc_dump_field_1d(&
          &  (inv_s_solid(ibeg:iend,ibeg:iend,:)*f(ibeg:iend,ibeg:iend,:)), &
          &  glen, filename1(2:), appisnap)
    call nc_dump_field_1d(g(ibeg:iend,ibeg:iend,:), glen, filename2(2:), appisnap)
  else

    open(unit=39000+mynum,file=datapath(1:lfdata)//filename1//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    write(39000+mynum) inv_s_solid(ibeg:iend,ibeg:iend,:)* &
                       f(ibeg:iend,ibeg:iend,:)
    close(39000+mynum)

    open(unit=35000+mynum,file=datapath(1:lfdata)//filename2//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    write(35000+mynum) g(ibeg:iend,ibeg:iend,:)
    close(35000+mynum)

  end if

end subroutine dump_field_over_s_solid_and_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_half_field_over_s_solid_1d_add(f,g,filename,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the solid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_solid
use pointwise_derivatives, ONLY: dsdf_solid_allaxis
use data_source, ONLY : have_src,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in)    :: f(0:npol,0:npol,nel_solid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_solid)
character(len=16), intent(in)     :: filename
character(len=4), intent(in)      :: appisnap
real(kind=realkind)               :: dsdf(0:npol,naxel_solid)
integer                           :: iel, glen

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_solid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_solid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_solid
        inv_s_solid(0,:,ax_el_solid(iel))=dsdf(:,iel)
     enddo
  endif

  g = real(.5,kind=realkind) * ( inv_s_solid * f + g )

  if (have_src .and. src_dump_type == 'mask') &
       call eradicate_src_elem_values(g)

  glen = size(g(ibeg:iend,ibeg:iend,:))
  if (use_netcdf) then
    call nc_dump_field_1d(g(ibeg:iend,ibeg:iend,:), glen, filename(2:), appisnap)
  else
    open(unit=35000+mynum,file=datapath(1:lfdata)//filename//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    write(35000+mynum) g(ibeg:iend,ibeg:iend,:)
    close(35000+mynum)
  end if

end subroutine dump_half_field_over_s_solid_1d_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_field_over_s_fluid_and_add(f,g,filename1,filename2,appisnap)
!
! This routine acts like the one above, calculating the term f/s in the fluid,
! but additionally adds field g to the dump. This is convenient for the 
! strain trace, where (dsus+dzuz) has been computed beforehand.
!
use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid!,deviator
use pointwise_derivatives, ONLY: dsdf_fluid_allaxis
use data_source, ONLY : src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: f(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(inout) :: g(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename1,filename2
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel,i,glen

  i=0

! TARJE JAN 14 2009: Confusion. Methinks this needs to be such that 
! inv_s_fluid is always multiplied into the wavefield... not just at 
! the axis! Changed it to this case.... not warranted though!

  if (have_axis) then 
     call dsdf_fluid_allaxis(f,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

  glen = size(f(ibeg:iend,ibeg:iend,:))
  if (use_netcdf) then
    call nc_dump_field_1d(inv_s_fluid(ibeg:iend,ibeg:iend,:)* f(ibeg:iend,ibeg:iend,:),&
      &  glen, filename1(2:), appisnap)
    call nc_dump_field_1d(inv_s_fluid(ibeg:iend,ibeg:iend,:)* f(ibeg:iend,ibeg:iend,:) &
      & + g(ibeg:iend,ibeg:iend,:), glen, filename2(2:), appisnap)
  else !Binary
  ! f/s (e.g. Epp)
    open(unit=39000+mynum,file=datapath(1:lfdata)//filename1//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    write(39000+mynum) inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       f(ibeg:iend,ibeg:iend,:)
    close(39000+mynum)



  ! sum of f/s and g (e.g. straintrace)
    open(unit=35000+mynum,file=datapath(1:lfdata)//filename2//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    write(35000+mynum) inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       f(ibeg:iend,ibeg:iend,:) + g(ibeg:iend,ibeg:iend,:)
    close(35000+mynum)
  end if 

end subroutine dump_field_over_s_fluid_and_add
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_half_f1_f2_over_s_fluid(f1,f2,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid
use pointwise_derivatives, ONLY : dsdf_fluid_allaxis
include 'mesh_params.h'

real(kind=realkind),intent(in) :: f1(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(in) :: f2(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel,glen

  if (have_axis) then
     call dsdf_fluid_allaxis(f2,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

  glen = size(f1(ibeg:iend,ibeg:iend,:))
  if (use_netcdf) then
    call nc_dump_field_1d(0.5*(f1(ibeg:iend,ibeg:iend,:) + &
                       &   inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       &   f2(ibeg:iend,ibeg:iend,:)),&
                       &   glen, filename(2:), appisnap)
  else
    open(unit=65000+mynum,file=datapath(1:lfdata)//filename//'_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")

    write(65000+mynum) 0.5*(f1(ibeg:iend,ibeg:iend,:) + &
                       inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       f2(ibeg:iend,ibeg:iend,:))

    close(65000+mynum)
  end if

end subroutine dump_half_f1_f2_over_s_fluid
!=============================================================================

!--------------------------------------------------------------------------
subroutine dump_f1_f2_over_s_fluid(f1,f2,filename,appisnap)

use data_proc, ONLY : appmynum
use data_pointwise, ONLY: inv_s_fluid
use pointwise_derivatives, ONLY : dsdf_fluid_allaxis
include 'mesh_params.h'

real(kind=realkind),intent(in) :: f1(0:npol,0:npol,nel_fluid)
real(kind=realkind),intent(in) :: f2(0:npol,0:npol,nel_fluid)
character(len=16), intent(in)  :: filename
character(len=4), intent(in)   :: appisnap
real(kind=realkind)            :: dsdf(0:npol,naxel_fluid)
integer                        :: iel,glen

  if (have_axis) then
     call dsdf_fluid_allaxis(f2,dsdf) ! axial f/s
     do iel=1,naxel_fluid
        inv_s_fluid(0,:,ax_el_fluid(iel))=dsdf(:,iel)
     enddo
  endif

  glen = size(f1(ibeg:iend,ibeg:iend,:))
  if (use_netcdf) then
    call nc_dump_field_1d((f1(ibeg:iend,ibeg:iend,:) + &
                       &   inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       &   f2(ibeg:iend,ibeg:iend,:)),&
                       &   glen, filename(2:), appisnap)
  
   else 
     open(unit=65000+mynum,file=datapath(1:lfdata)//filename//'_'&
          //appmynum//'_'//appisnap//'.bindat',&
          FORM="UNFORMATTED",STATUS="REPLACE")

     write(65000+mynum) f1(ibeg:iend,ibeg:iend,:) + &
                       inv_s_fluid(ibeg:iend,ibeg:iend,:)* &
                       f2(ibeg:iend,ibeg:iend,:)

     close(65000+mynum)
  end if

end subroutine dump_f1_f2_over_s_fluid
!=============================================================================


!--------------------------------------------------------------------------
subroutine dump_disp(u,chi)

use data_source, ONLY : src_type,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: u(0:npol,0:npol,nel_solid,3)
real(kind=realkind),intent(in) :: chi(0:npol,0:npol,nel_fluid)
integer                        :: i, glen
character(len=4)               :: appisnap
real(kind=realkind)            :: f(0:npol,0:npol,nel_solid,3)

  call define_io_appendix(appisnap,istrain)

  f = u

  if (src_dump_type == 'mask') then
    call eradicate_src_elem_vec_values(f)
  end if


! Dump solid displacement
  glen = size(f(ibeg:iend,ibeg:iend,:,1))

  if (use_netcdf) then
    if (src_type(1)/='monopole') then
      call nc_dump_field_1d((f(ibeg:iend,ibeg:iend,:,:)), &
                       &   glen*3, 'disp_sol', appisnap)
    else
      call nc_dump_field_1d(f(ibeg:iend,ibeg:iend,:,(/1,3/)), &
                       &   glen*2, 'disp_sol', appisnap)
    end if
  else
    open(unit=75000+mynum,file=datapath(1:lfdata)//'/disp_sol_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")

    if (src_type(1)/='monopole') then
       write(75000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
    else
       write(75000+mynum) f(ibeg:iend,ibeg:iend,:,1:3:2)
    endif
    close(75000+mynum)
  end if

! Dump fluid potential 
  if (have_fluid) then 
    if (use_netcdf) then
      glen = size(chi)
      call nc_dump_field_1d( chi, glen, 'chi_flu', appisnap)
    else
      open(unit=76000+mynum,file=datapath(1:lfdata)//'/chi_flu_'&
                                //appmynum//'_'//appisnap//'.bindat',&
                                FORM="UNFORMATTED",STATUS="REPLACE")
  
      write(76000+mynum)chi
      close(76000+mynum)
    end if
  endif 

end subroutine dump_disp
!=============================================================================

!--------------------------------------------------------------------------
subroutine dump_velo_dchi(v,dchi)

use data_source, ONLY : src_type,src_dump_type

include 'mesh_params.h'

real(kind=realkind),intent(in) :: v(0:npol,0:npol,nel_solid,3)
real(kind=realkind),intent(in) :: dchi(0:npol,0:npol,nel_fluid)
integer                        :: i, glen
character(len=4)               :: appisnap
real(kind=realkind)            :: f(0:npol,0:npol,nel_solid,3)

  call define_io_appendix(appisnap,istrain)

  f = v

  if (src_dump_type == 'mask') then
       call eradicate_src_elem_vec_values(f)
  end if
  
  glen = size(f(ibeg:iend,ibeg:iend,:,1))

! Dump solid velocity vector
  if (use_netcdf) then
    if (src_type(1)/='monopole') then
      call nc_dump_field_1d((f(ibeg:iend,ibeg:iend,:,:)), &
                       &   glen*3, 'velo_sol', appisnap)
    else
      call nc_dump_field_1d(f(ibeg:iend,ibeg:iend,:,(/1,3/)), &
                       &   glen*2, 'velo_sol', appisnap)
    end if
  
  else ! Binary
    open(unit=85000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")

    if (src_type(1)/='monopole') then 
       write(85000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
    else
       write(85000+mynum) f(ibeg:iend,ibeg:iend,:,1),f(ibeg:iend,ibeg:iend,:,3)
    endif
    close(85000+mynum)
  end if

! Dump fluid potential 1st derivative
  if (have_fluid) then 
    if (use_netcdf) then
      glen = size(dchi)
      call nc_dump_field_1d( dchi, glen, 'dchi_flu', appisnap)
    else
      open(unit=86000+mynum,file=datapath(1:lfdata)//'/dchi_flu_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
  
      write(86000+mynum)dchi
      close(86000+mynum)
    end if
  endif


end subroutine dump_velo_dchi
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_velo_global(v,dchi)

use data_pointwise, ONLY: inv_rho_fluid,inv_s_rho_fluid,usz_fluid
use data_source, ONLY : src_type,src_dump_type
use pointwise_derivatives, ONLY: axisym_laplacian_fluid,dsdf_fluid_allaxis
use unit_stride_colloc, ONLY : collocate0_1d

include 'mesh_params.h'

real(kind=realkind),intent(in) :: v(:,:,:,:)
real(kind=realkind),intent(in) :: dchi(:,:,:)

real(kind=realkind)            :: phicomp(0:npol,0:npol,nel_fluid)
integer                        :: i,iel,glen
character(len=4)               :: appisnap
real(kind=realkind)            :: dsdchi(0:npol,naxel_fluid)
real(kind=realkind)            :: f(0:npol,0:npol,1:nel_solid,3)
real(kind=realkind)            :: fflu(0:npol,0:npol,1:nel_fluid,3)

  call define_io_appendix(appisnap,istrain)

! sssssssssssss dump velocity vector inside solid ssssssssssssssssssssssssssss


  f=v
  if (src_dump_type == 'mask') then
    call eradicate_src_elem_vec_values(f)
  end if

  glen = size(f(ibeg:iend,ibeg:iend,:,1))

  if (use_netcdf) then
    if (src_type(1)/='monopole') then
      call nc_dump_field_1d((f(ibeg:iend,ibeg:iend,:,:)), &
                       &   glen*3, 'velo_sol', appisnap)
    else
      call nc_dump_field_1d(f(ibeg:iend,ibeg:iend,:,(/1,3/)), &
                       &   glen*2, 'velo_sol', appisnap)
    end if
  else
    open(unit=95000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                              //appmynum//'_'//appisnap//'.bindat',&
                              FORM="UNFORMATTED",STATUS="REPLACE")
    if (src_type(1)/='monopole') then
       write(95000+mynum) (f(ibeg:iend,ibeg:iend,:,i),i=1,3)
    else
       write(95000+mynum) f(ibeg:iend,ibeg:iend,:,1), &
                          f(ibeg:iend,ibeg:iend,:,3)
    end if
    close(95000+mynum)
  end if

! ffffffff fluid region ffffffffffffffffffffffffffffffffffffffffffffffffffffff

  if (have_fluid) then 
! compute velocity vector inside fluid
    call axisym_laplacian_fluid(dchi,usz_fluid)

! phi component needs special care: m/(s rho) dchi
    call collocate0_1d(inv_s_rho_fluid,dchi,phicomp,npoint_fluid)

! Take care of axial singularity for phi component of fluid velocity
    if (have_axis) then
       call dsdf_fluid_allaxis(dchi,dsdchi)
       do iel=1,naxel_fluid
          phicomp(0,:,ax_el_fluid(iel))=dsdchi(:,iel)* &
               phicomp(0,:,ax_el_fluid(iel))
       enddo
    endif

    call define_io_appendix(appisnap,istrain)

    fflu(ibeg:iend,ibeg:iend,:,1) = inv_rho_fluid(ibeg:iend,ibeg:iend,:) * &
                        &        usz_fluid(ibeg:iend,ibeg:iend,:,1)
    fflu(ibeg:iend,ibeg:iend,:,2) = phicomp
    fflu(ibeg:iend,ibeg:iend,:,3) = inv_rho_fluid(ibeg:iend,ibeg:iend,:) * &
                        &        usz_fluid(ibeg:iend,ibeg:iend,:,2)      
  ! dump velocity vector inside fluid
    if (use_netcdf) then
      call nc_dump_field_1d(fflu(ibeg:iend,ibeg:iend,:,:), &
                       &   glen*2, 'velo_flu', appisnap)
    else
      open(unit=960000+mynum,file=datapath(1:lfdata)//'/velo_flu_'&
                                //appmynum//'_'//appisnap//'.bindat',&
                                 FORM="UNFORMATTED",STATUS="REPLACE")

      write(960000+mynum) (fflu(ibeg:iend,ibeg:iend,:,i),i=1,3)
 !                         nv_rho_fluid(ibeg:iend,ibeg:iend,:)* &
 !                         usz_fluid(ibeg:iend,ibeg:iend,:,1),phicomp, &
 !                         inv_rho_fluid(ibeg:iend,ibeg:iend,:)* &
 !                         usz_fluid(ibeg:iend,ibeg:iend,:,2)
      close(960000+mynum)
    end if ! netcdf
  endif ! have_fluid

end subroutine dump_velo_global
!=============================================================================


!-----------------------------------------------------------------------------
subroutine eradicate_src_elem_vec_values(u)
!
! Deletes all entries to vector field u on ALL GLL points inside
! elements that have a non-zero source term (i.e. including all 
! assembled neighboring elements)
! This is a preliminary test for the wavefield dumps.

use data_source, ONLY : nelsrc,ielsrc,have_src

include 'mesh_params.h'

real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel_solid,3)
integer :: iel

if (have_src) then
   do iel=1,nelsrc
      u(0:npol,0:npol,ielsrc(iel),1:3) = real(0.,kind=realkind)
   enddo
endif

end subroutine eradicate_src_elem_vec_values
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------------
subroutine eradicate_src_elem_values(u)
!
! Deletes all entries to scalar field u on ALL GLL points inside
! elements that have a non-zero source term (i.e. including all 
! assembled neighboring elements)
! This is a preliminary test for the wavefield dumps.

use data_source, ONLY : nelsrc,ielsrc,have_src

include 'mesh_params.h'

real(kind=realkind),intent(inout) :: u(0:npol,0:npol,nel_solid)
integer :: iel

if (have_src) then

   do iel=1,nelsrc
      u(0:npol,0:npol,ielsrc(iel)) = real(0.,kind=realkind)
   enddo
endif

end subroutine eradicate_src_elem_values

!================================
end module wavefields_io
!================================
