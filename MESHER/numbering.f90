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

  module numbering

  use data_gllmesh
  use data_numbering
  use data_spec
  use data_mesh
  use data_grid
  use data_diag

  implicit none
  public :: define_global_global_numbering
  public :: define_global_flobal_numbering
  public :: define_global_slobal_numbering
  public :: get_global
  private
  contains

!--------------------------------------------------------------------------
! dk define_global_global_numbering----------------------------------------
subroutine define_global_global_numbering

  integer npointot
  double precision, dimension(:), allocatable :: sgtmp,zgtmp
  logical, dimension(:), allocatable ::   ifseg
  integer, dimension(:), allocatable :: loc
  integer :: iel, jpol,ipol, ipt
!
  ngllcube = (npol+1)**2 
  npointot = neltot * (npol+1)**2

  if (dump_mesh_info_screen) then
   write(6,*) 
   write(6,*) 'NPOINTOT GLOBAL IS ' , npointot
  end if
!  
if (dump_mesh_info_files) then
  open(2,file=diagpath(1:lfdiag)//'/crds',form="UNFORMATTED")
   write(2) sgll
   write(2) zgll
  close(2)
endif 

  allocate(sgtmp(npointot)) ; sgtmp(:) = 0. 
  do iel = 1, neltot
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     sgtmp(ipt) = sgll(ipol,jpol,iel)
    end do
   end do
  end do
  if (dump_mesh_info_files)  deallocate(sgll)
  allocate(zgtmp(npointot)) ; zgtmp(:) = 0.d0 
  do iel = 1, neltot
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     zgtmp(ipt) = zgll(ipol,jpol,iel)
    end do
   end do
  end do
  if (dump_mesh_info_files)  deallocate(zgll)

  allocate(iglob(npointot)); iglob(:) = 0
  allocate(loc(npointot)); loc(:) = 0
  allocate(ifseg(npointot))
  call get_global(neltot,sgtmp,zgtmp,iglob,loc,ifseg,nglobglob,npointot,ngllcube,NDIM)
  deallocate(ifseg)
  deallocate(loc)
  deallocate(sgtmp)
  deallocate(zgtmp)

  if (dump_mesh_info_files) then
  allocate(zgll(0:npol,0:npol,neltot))
  allocate(sgll(0:npol,0:npol,neltot))
  open(2,file=diagpath(1:lfdiag)//'/crds',form="unformatted")
   read(2) sgll
   read(2) zgll
  close(2)
endif

  if (dump_mesh_info_screen) write(6,*) 'NGLOBGLOB IS ' , NGLOBGLOB

end subroutine define_global_global_numbering
!--------------------------------------------------------------------------
!
!dk define_global_flobal_numbering-----------------------------------------
  subroutine define_global_flobal_numbering
  integer npointot
  double precision, dimension(:), allocatable :: sgtmp,zgtmp
  integer, dimension(:), allocatable :: loc_fluid
  logical, dimension(:), allocatable ::   ifseg
  integer :: iel, jpol,ipol, ipt
!

  npointot = neltot_fluid * (npol+1)**2
!
  if (dump_mesh_info_screen) then 
   write(6,*) 
   write(6,*) 'NPOINTOT FLOBAL IS ' , npointot
  end if
!
if (dump_mesh_info_files) then
  open(2,file=diagpath(1:lfdiag)//'/crds',form="UNFORMATTED")
   write(2) sgll_fluid
   write(2) zgll_fluid
  close(2)
endif

  allocate(sgtmp(npointot))
  do iel = 1, neltot_fluid
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     sgtmp(ipt) = sgll_fluid(ipol,jpol,iel)
    end do
   end do
  end do

  deallocate(sgll_fluid)
  allocate(zgtmp(npointot))
  do iel = 1, neltot_fluid
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     zgtmp(ipt) = zgll_fluid(ipol,jpol,iel)
    end do
   end do
  end do
  deallocate(zgll_fluid)
!
  allocate(iglob_fluid(npointot)) ; iglob_fluid(:) = 0
  allocate(loc_fluid(npointot)) ;   loc_fluid(:) = 0
  allocate(ifseg(npointot))
!
  call get_global(neltot_fluid,sgtmp,zgtmp,iglob_fluid,loc_fluid,ifseg,nglobflob,npointot,NGLLcube,NDIM)
!
  deallocate(ifseg)
  deallocate(loc_fluid)
  deallocate(zgtmp)
  deallocate(sgtmp)
!
! allocate(zgll_fluid(0:npol,0:npol,neltot_fluid))
! allocate(sgll_fluid(0:npol,0:npol,neltot_fluid))
! open(2,file='crds',form="UNFORMATTED")
!  read(2) sgll_fluid
!  read(2) zgll_fluid
! close(2)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBFLOB IS ' , NGLOBFLOB

  end subroutine define_global_flobal_numbering
!
!-------------------------------------------------------------------------
! dk define_global_slobal_numbering---------------------------------------
  subroutine define_global_slobal_numbering
  integer npointot
  double precision, dimension(:), allocatable :: sgtmp,zgtmp
  integer, dimension(:), allocatable :: loc_solid
  logical, dimension(:), allocatable ::   ifseg
!
  integer :: iel, jpol,ipol, ipt
!
! test 
!  double precision, dimension(:), allocatable :: utest, uglob

  npointot = neltot_solid * (npol+1)**2
!
  if (dump_mesh_info_screen) then 
   write(6,*) 
   write(6,*) 'NPOINTOT SLOBAL IS ' , npointot
  end if
! To save some memory 
if (dump_mesh_info_files) then
  open(2,file=diagpath(1:lfdiag)//'/crds',form="UNFORMATTED")
   write(2) sgll
   write(2) zgll
  close(2)
  deallocate(sgll,zgll)
endif
!
  allocate(sgtmp(npointot))
  do iel = 1, neltot_solid
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     sgtmp(ipt) = sgll_solid(ipol,jpol,iel)
    end do
   end do
  end do
  deallocate(sgll_solid) ! not needed anymore 
  allocate(zgtmp(npointot))
  do iel = 1, neltot_solid
   do jpol = 0, npol
    do ipol = 0, npol
     ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
     zgtmp(ipt) = zgll_solid(ipol,jpol,iel)
    end do
   end do
  end do
  deallocate(zgll_solid) ! not needed anymore
!
  allocate(iglob_solid(npointot)) ; iglob_solid(:) = 0
  allocate(loc_solid(npointot)) ;   loc_solid(:) = 0
  allocate(ifseg(npointot))

  call get_global(neltot_solid,sgtmp,zgtmp,iglob_solid,loc_solid,ifseg,nglobslob,npointot,NGLLcube,NDIM)

  deallocate(ifseg)
  deallocate(loc_solid)
  deallocate(zgtmp)
  deallocate(sgtmp)
! now load global coordinate arrays back in
if (dump_mesh_info_files) then
  allocate(zgll(0:npol,0:npol,neltot))
  allocate(sgll(0:npol,0:npol,neltot))
  open(2,file=diagpath(1:lfdiag)//'/crds',form="UNFORMATTED")
   read(2) sgll
   read(2) zgll
  close(2)
endif
!
  if (dump_mesh_info_screen) write(6,*) 'NGLOBSLOB IS ' , NGLOBSLOB 
!
  end subroutine define_global_slobal_numbering
!

!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

subroutine get_global(nspec2,xp,yp,iglob2,loc2,ifseg2,nglob2,npointot2,NGLLCUBE2,NDIM2)

  ! this routine MUST be in double precision to avoid sensitivity
  ! to roundoff errors in the coordinates of the points

  ! non-structured global numbering software provided by Paul F. Fischer

  ! leave sorting subroutines in same source file to allow for inlining

  implicit none

  !  include "constants.h"

  integer, intent(in) ::  nspec2,npointot2,NGLLCUBE2,NDIM2
  double precision, intent(in) ::  xp(npointot2),yp(npointot2)
  integer, intent(out) :: iglob2(npointot2),loc2(npointot2),nglob2
  logical, intent(out) :: ifseg2(npointot2)

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! TNM: that's what I had
! double precision, parameter :: SMALLVALTOL = 1.d-15
  double precision, parameter :: SMALLVALTOL = 1.d-08

! write(6,*)'GLOBAL NUMBERING npointot2,nspec2,NGLLCUBE2:',npointot2,nspec2,NGLLCUBE2
! write(6,*)'GLOBAL NUMBERING xp yp max:', maxval(abs(xp)),maxval(abs(yp))

! establish initial pointers
  do ispec=1,nspec2
     ieoff=NGLLCUBE2*(ispec-1)
     do ilocnum=1,NGLLCUBE2
        loc2(ilocnum+ieoff)=ilocnum+ieoff
     enddo
  enddo

  ifseg2(:)=.false.

! dynamically allocate arrays
  allocate(ind(npointot2))
  allocate(ninseg(npointot2))
  allocate(iwork(npointot2))
  allocate(work(npointot2))

  nseg=1
  ifseg2(1)=.true.
  ninseg(1)=npointot2

!==========================================
  do j=1,NDIM2
!==========================================

! sort within each segment
     ioff=1
     do iseg=1,nseg
        if(j == 1) then
           call rank(xp(ioff),ind,ninseg(iseg))
        else
           call rank(yp(ioff),ind,ninseg(iseg))
        endif
!af
        call swap_all(loc2(ioff),xp(ioff),yp(ioff),iwork,work,ind,ninseg(iseg))
!end af
        ioff=ioff+ninseg(iseg)
     enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
     if(j == 1) then
        do i=2,npointot2
           if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg2(i)=.true.
!          if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) write(6666,*)'DISTANCE X:',i,loc2(i),dabs(xp(i)-xp(i-1))
!          if(dabs(xp(i)-xp(i-1)) < SMALLVALTOL) write(6667,*)'DISTANCE X:',i,loc2(i),dabs(xp(i)-xp(i-1))
        enddo
     else
        do i=2,npointot2
           if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg2(i)=.true.
!          if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) write(6666,*)'DISTANCE Y:',i,loc2(i),dabs(yp(i)-yp(i-1))
!          if(dabs(yp(i)-yp(i-1)) < SMALLVALTOL) write(6667,*)'DISTANCE Y:',i,loc2(i),dabs(yp(i)-yp(i-1))
        enddo

     endif

! count up number of different segments
     nseg=0
     do i=1,npointot2
        if(ifseg2(i)) then
           nseg=nseg+1
           ninseg(nseg)=1
        else
           ninseg(nseg)=ninseg(nseg)+1
        endif

     enddo

!==========================================
  enddo ! NDIM2 loop
!==========================================

! deallocate arrays
  deallocate(ind)
  deallocate(iwork)
  deallocate(work)
  deallocate(ninseg)

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot2
     if(ifseg2(i)) ig=ig+1
     iglob2(loc2(i))=ig
  enddo
  nglob2=ig

end subroutine get_global
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! sorting routines put in same file to allow for inlining

subroutine rank(A,IND,N)
  !
  ! Use Heap Sort (Numerical Recipes)
  !
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
     IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
100 CONTINUE
  IF (l>1) THEN
     l=l-1
     indx=ind(l)
     q=a(indx)
  ELSE
     indx=ind(ir)
     q=a(indx)
     ind(ir)=ind(1)
     ir=ir-1
     if (ir == 1) then
        ind(1)=indx

        return
     endif
  ENDIF
  i=l
  j=l+l
200 CONTINUE
  IF (J <= IR) THEN
     IF (J<IR) THEN
        IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
     ENDIF
     IF (q<A(IND(j))) THEN
        IND(I)=IND(J)
        I=J
        J=J+J
     ELSE
        J=IR+1
     ENDIF
     goto 200
  ENDIF
  IND(I)=INDX
  goto 100

end subroutine rank

! ------------------------------------------------------------------

subroutine swap_all(IA,A,B,IW,W,ind,n)
  !
  ! swap arrays IA, A, B and C according to addressing in array IND
  !
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
     IA(i)=IW(ind(i))
     A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
     B(i)=W(ind(i))
  enddo

end subroutine swap_all


!=========================
  end module numbering
!=========================
