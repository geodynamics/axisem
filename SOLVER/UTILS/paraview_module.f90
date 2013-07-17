!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stephanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!==================
 module paraview
!==================

 use data_mesh
 use data_spec
 implicit none
 public :: paraview_test
 public :: paraview_one_chunk
 private
 contains
!//////////////////////////////////////////////////////////////////////////
!
!dk paraview_test----------------------------------------------------------
 subroutine paraview_test(scalar_field)
!
!    This routine outputs velocity and temperature
! in a format (The VISUAL TOOLKIT format) compliant with
! the use of the Paraview (http://www.paraview.org) visualization software.
!
 real, dimension(0:npol,0:npol,0:nrpol,1:nelt) :: scalar_field
 integer :: ielem, ipol, jpol, iphi,j, kpol
 integer :: ipt,npt,ncell,nsize_cell,icell,nelem_axial
 integer,dimension(:,:,:,:), allocatable :: num
 integer,dimension(:),allocatable :: elem_invnum
 integer,dimension(:),allocatable :: ipol_invnum
 integer,dimension(:),allocatable :: jpol_invnum
 integer,dimension(:),allocatable :: kpol_invnum
 integer :: k
 real :: s,z,phi,x,y,r,theta,cp,sp
 character(len=4) :: appcount
!
 real :: t,vx,vy,vz,vs,vp
!
!----------------------------------------------------------------------------------
 nelem_axial = 0 ! sum(axis,dim=1) ! Number of axial elements
 npt = nelt*(nrpol+1)*(npol+1)**2
                               ! Total number of grid points
! allocate(num(0:npol,0:npol,nelem,0:N_phi-1))
 allocate(num(0:npol,0:npol,0:nrpol,nelt))
 num(:,:,:,:) = 0
 allocate(elem_invnum(0:npt-1),ipol_invnum(0:npt-1),&
          jpol_invnum(0:npt-1),kpol_invnum(0:npt-1))
 elem_invnum(:) = 0 ; ipol_invnum(:) = 0
 jpol_invnum(:) = 0 ; kpol_invnum(:) = 0
!----------------------------------------------------------------------------------
! Open Visual ToolKit file
! open(10,file=datapath(1:lfdata)//'/testpara'//appmynum//'.vtk')
 open(10,file='./scalarfield_cs.vtk')
!
! write header
 write(10,20)  
! write points information
 write(10,30) npt
!
 ipt = -1

 do ielem = 1, nelt
       do kpol = 0, nrpol
       do jpol = 0, npol
       do ipol = 0, npol
         ipt = ipt + 1 ; num(ipol,jpol,kpol,ielem) = ipt
         elem_invnum(ipt) = ielem
         ipol_invnum(ipt) = ipol
         jpol_invnum(ipt) = jpol
         kpol_invnum(ipt) = kpol
         x = xcol(ipol,jpol,kpol,ielem)
         y = ycol(ipol,jpol,kpol,ielem)
         z = zcol(ipol,jpol,kpol,ielem)
         write(10,40) x,y,z
       end do
       end do
       end do
 end do
!
!-------------------------------------------------------------------------------------
! Compute total number of cells
!
! ncell = nelem_axial*(npol*N_phi + (npol-1)*npol*N_phi) &
!       +(nelem-nelem_axial)*(N_phi*npol**2)
 ncell = nelt*nrpol*npol**2
! nsize_cell = nelem_axial*(npol*(N_phi)*7 + (N_phi)*(npol-1)*npol*9)&
!            +(nelem-nelem_axial)*(9*N_phi*npol**2)
 nsize_cell = nelt*(9*nrpol*npol**2)
! Write this info into VTK file
 write(10,50) ncell, nsize_cell
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
! Now write cells topology
 do ielem = 1, nelt
       do ipol = 0, npol-1 ! HEXAHEDRA ONLY
       do jpol = 0, npol-1
       do kpol = 0, nrpol-1
             write(10,60)                              8,&
                         num(ipol  ,jpol  ,kpol  ,ielem),&
                         num(ipol+1,jpol  ,kpol  ,ielem),&
                         num(ipol+1,jpol  ,kpol+1,ielem),&
                         num(ipol  ,jpol  ,kpol+1,ielem),&
                         num(ipol  ,jpol+1,kpol  ,ielem),&
                         num(ipol+1,jpol+1,kpol  ,ielem),&
                         num(ipol+1,jpol+1,kpol+1,ielem),&
                         num(ipol  ,jpol+1,kpol+1,ielem)
       end do
       end do
       end do
 end do
!
!---------------------------------------------------------------------
! Now write down cell types list
!
 write(10,70) ncell
 do ielem = 1, nelt
    do icell = 1, nrpol*(npol**2)
       write(10,80) 12
    end do
 end do
!---------------------------------------------------------------------
 write(10,90) npt
 write(10,100)
 do ipt = 0, npt-1
    ipol = ipol_invnum(ipt)
    jpol = jpol_invnum(ipt)
    kpol = kpol_invnum(ipt)
    ielem = elem_invnum(ipt)
    write(10,110) scalar_field(ipol,jpol,kpol,ielem)
 end do
 close(10)
!
!-------------------------------------------------------------------------------
!
! Deallocate arrays
 deallocate(elem_invnum,ipol_invnum,jpol_invnum,kpol_invnum)
 deallocate(num)

20 format('# vtk DataFile Version 2.0',/&
         'THIS IS A  TEST FILE',/&
         'ASCII               ',//,&
         'DATASET UNSTRUCTURED_GRID')
30 format('POINTS',x,i7,x,'float')
40 format(3(1pe12.5,2x))
50 format(/,'CELLS',1x,i7,1x,i7)
60 format(i1,10(i7,x))
70 format(/,'CELL_TYPES',1x,i7)
80 format(i2)
90 format(/,'POINT_DATA',1x,i7)
100 format('SCALARS Temperature float',/&
          'LOOKUP_TABLE default')
101 format('SCALARS KE float',/&
          'LOOKUP_TABLE default')
102 format('SCALARS Helicity float',/&
          'LOOKUP_TABLE default')
103 format('SCALARS dTdr float',/&
          'LOOKUP_TABLE default')
110 format(1pe12.5)
120 format(/,'VECTORS velocity float')
130 format(1pe12.5,2x,1pe12.5,2x,1pe12.5)

 end subroutine paraview_test
!------------------------------------------------------------------------------
!
!dk paraview_onechunk----------------------------------------------------------
 subroutine paraview_one_chunk(scalar_field,ichunk)
!
!    This routine outputs velocity and temperature
! in a format (The VISUAL TOOLKIT format) compliant with
! the use of the Paraview (http://www.paraview.org) visualization software.
!
 real, dimension(0:npol,0:npol,0:nrpol,1:nelt) :: scalar_field
 integer, intent(in) :: ichunk
 integer :: ielem, ipol, jpol, iphi,j, kpol
 integer :: ipt,npt,ncell,nsize_cell,icell,nelem_axial
 integer,dimension(:,:,:,:), allocatable :: num
 integer,dimension(:),allocatable :: elem_invnum
 integer,dimension(:),allocatable :: ipol_invnum
 integer,dimension(:),allocatable :: jpol_invnum
 integer,dimension(:),allocatable :: kpol_invnum
 integer :: k
 real :: s,z,phi,x,y,r,theta,cp,sp
 character(len=4) :: appchunk
 integer :: numerator
 integer :: denominator
 integer :: ibegelem
!
 real :: t,vx,vy,vz,vs,vp
!
 call define_io_appendix(appchunk,ichunk)
!
!----------------------------------------------------------------------------------
 nelem_axial = 0 ! sum(axis,dim=1) ! Number of axial elements
 npt = nelt*(nrpol+1)*(npol+1)**2
                               ! Total number of grid points
! allocate(num(0:npol,0:npol,nelem,0:N_phi-1))
 allocate(num(0:npol,0:npol,0:nrpol,nelt))
 num(:,:,:,:) = 0
 allocate(elem_invnum(0:npt-1),ipol_invnum(0:npt-1),&
          jpol_invnum(0:npt-1),kpol_invnum(0:npt-1))
 elem_invnum(:) = 0 ; ipol_invnum(:) = 0
 jpol_invnum(:) = 0 ; kpol_invnum(:) = 0
!----------------------------------------------------------------------------------
! Open Visual ToolKit file
! open(10,file=datapath(1:lfdata)//'/testpara'//appmynum//'.vtk')
 open(10,file='./scalarfield_chunk'//appchunk//'.vtk')
!
! write header
 write(10,20)  
! write points information
! write(10,30) npt/6
! write(10,30) npt/2
 numerator = 1
 denominator = 6
 write(10,30) numerator*npt/denominator
!
 ipt = -1
!
 ibegelem = (ichunk-1)*(nelt/6)
!
 do ielem = ibegelem+1, ibegelem+numerator*nelt/denominator
       do kpol = 0, nrpol
       do jpol = 0, npol
       do ipol = 0, npol
         ipt = ipt + 1 ; num(ipol,jpol,kpol,ielem) = ipt
         elem_invnum(ipt) = ielem
         ipol_invnum(ipt) = ipol
         jpol_invnum(ipt) = jpol
         kpol_invnum(ipt) = kpol
         x = xcol(ipol,jpol,kpol,ielem)
         y = ycol(ipol,jpol,kpol,ielem)
         z = zcol(ipol,jpol,kpol,ielem)
         write(10,40) x,y,z
       end do
       end do
       end do
 end do
!
!-------------------------------------------------------------------------------------
! Compute total number of cells
 ncell = numerator*nelt*nrpol*npol**2/denominator
 nsize_cell = numerator*nelt*(9*nrpol*npol**2)/denominator
! Write this info into VTK file
 write(10,50) ncell, nsize_cell
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
! Now write cells topology
 do ielem = ibegelem+1, ibegelem+numerator*nelt/denominator
       do ipol = 0, npol-1 ! HEXAHEDRA ONLY
       do jpol = 0, npol-1
       do kpol = 0, nrpol-1
             write(10,60)                              8,&
                         num(ipol  ,jpol  ,kpol  ,ielem),&
                         num(ipol+1,jpol  ,kpol  ,ielem),&
                         num(ipol+1,jpol  ,kpol+1,ielem),&
                         num(ipol  ,jpol  ,kpol+1,ielem),&
                         num(ipol  ,jpol+1,kpol  ,ielem),&
                         num(ipol+1,jpol+1,kpol  ,ielem),&
                         num(ipol+1,jpol+1,kpol+1,ielem),&
                         num(ipol  ,jpol+1,kpol+1,ielem)
       end do
       end do
       end do
 end do
!
!---------------------------------------------------------------------
! Now write down cell types list
!
 write(10,70) ncell
 do ielem = ibegelem+1, ibegelem+nelt*numerator/denominator
    do icell = 1, nrpol*(npol**2)
       write(10,80) 12
    end do
 end do
!---------------------------------------------------------------------
 write(10,90) npt*numerator/denominator
 write(10,100)
 do ipt = 0, npt*numerator/denominator-1
    ipol = ipol_invnum(ipt)
    jpol = jpol_invnum(ipt)
    kpol = kpol_invnum(ipt)
    ielem = elem_invnum(ipt)
    write(10,110) scalar_field(ipol,jpol,kpol,ielem)
 end do
 close(10)
!
!-------------------------------------------------------------------------------
!
! Deallocate arrays
 deallocate(elem_invnum,ipol_invnum,jpol_invnum,kpol_invnum)
 deallocate(num)

20 format('# vtk DataFile Version 2.0',/&
         'THIS IS A  TEST FILE',/&
         'ASCII               ',//,&
         'DATASET UNSTRUCTURED_GRID')
30 format('POINTS',x,i7,x,'float')
40 format(3(1pe12.5,2x))
50 format(/,'CELLS',1x,i7,1x,i7)
60 format(i1,10(i7,x))
70 format(/,'CELL_TYPES',1x,i7)
80 format(i2)
90 format(/,'POINT_DATA',1x,i7)
100 format('SCALARS Br float',/&
          'LOOKUP_TABLE default')
101 format('SCALARS KE float',/&
          'LOOKUP_TABLE default')
102 format('SCALARS Helicity float',/&
          'LOOKUP_TABLE default')
103 format('SCALARS dTdr float',/&
          'LOOKUP_TABLE default')
110 format(1pe12.5)
120 format(/,'VECTORS velocity float')
130 format(1pe12.5,2x,1pe12.5,2x,1pe12.5)

 end subroutine paraview_one_chunk
!------------------------------------------------------------------------------
!dk define_io_appendix
 subroutine define_io_appendix(app,iproc)
!
! Defines the 4 digit character string appended to any
! data or io file related to process myid.
!
 integer :: iproc
 character(len=4) :: app
 character(len=1) :: milp,cenp,dizp,unip
 milp = char(48+    iproc/1000)
 cenp = char(48+mod(iproc/100,10))
 dizp = char(48+mod(iproc/10,10))
 unip = char(48+mod(iproc,10))
 app = milp//cenp//dizp//unip
 end subroutine define_io_appendix
!--------------------------------------------------------------------------

!
!//////////////////////////////////////////////////////////////////////////////
!
!======================
 end module paraview
!======================
