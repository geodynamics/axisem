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
subroutine define_global_global_numbering

  double precision, dimension(:), allocatable   :: sgtmp, zgtmp
  logical, dimension(:), allocatable            :: ifseg
  integer, dimension(:), allocatable            :: loc
  integer   :: iel, jpol,ipol, ipt
  integer   :: npointot

  ngllcube = (npol + 1)**2 
  npointot = neltot * (npol+1)**2

  if (dump_mesh_info_screen) then
   write(6,*) 
   write(6,*) 'NPOINTOT GLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll, .true.)

  allocate(zgtmp(npointot))
  sgtmp = pack(zgll, .true.)

  allocate(iglob(npointot))
  iglob(:) = 0
  allocate(loc(npointot))
  loc(:) = 0
  allocate(ifseg(npointot))

  call get_global(neltot, sgtmp, zgtmp, iglob, loc, ifseg, nglobglob, npointot, ngllcube, NDIM)

  deallocate(ifseg)
  deallocate(loc)
  deallocate(sgtmp)
  deallocate(zgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBGLOB IS ' , NGLOBGLOB

end subroutine define_global_global_numbering
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine define_global_flobal_numbering
  integer npointot
  double precision, dimension(:), allocatable :: sgtmp,zgtmp
  integer, dimension(:), allocatable :: loc_fluid
  logical, dimension(:), allocatable ::   ifseg
  integer :: iel, jpol,ipol, ipt

  npointot = neltot_fluid * (npol+1)**2

  if (dump_mesh_info_screen) then 
     write(6,*) 
     write(6,*) 'NPOINTOT FLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll_fluid, .true.)
  deallocate(sgll_fluid)

  allocate(zgtmp(npointot))
  zgtmp = pack(zgll_fluid, .true.)
  deallocate(zgll_fluid)

  allocate(iglob_fluid(npointot))
  iglob_fluid(:) = 0
  allocate(loc_fluid(npointot))
  loc_fluid(:) = 0
  allocate(ifseg(npointot))

  call get_global(neltot_fluid, sgtmp, zgtmp, iglob_fluid, loc_fluid, ifseg, &
                  nglobflob, npointot, NGLLcube, NDIM)

  deallocate(ifseg)
  deallocate(loc_fluid)
  deallocate(zgtmp)
  deallocate(sgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBFLOB IS ' , NGLOBFLOB

end subroutine define_global_flobal_numbering
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine define_global_slobal_numbering

  double precision, dimension(:), allocatable   :: sgtmp, zgtmp
  integer, dimension(:), allocatable            :: loc_solid
  logical, dimension(:), allocatable            :: ifseg

  integer :: npointot
  integer :: iel, jpol,ipol, ipt

  npointot = neltot_solid * (npol+1)**2

  if (dump_mesh_info_screen) then 
     write(6,*) 
     write(6,*) 'NPOINTOT SLOBAL IS ' , npointot
  end if

  allocate(sgtmp(npointot))
  sgtmp = pack(sgll_solid, .true.)
  deallocate(sgll_solid) ! not needed anymore 
  
  allocate(zgtmp(npointot))
  zgtmp = pack(zgll_solid, .true.)
  deallocate(zgll_solid) ! not needed anymore

  allocate(iglob_solid(npointot))
  iglob_solid(:) = 0
  allocate(loc_solid(npointot))
  loc_solid(:) = 0
  allocate(ifseg(npointot))

  call get_global(neltot_solid, sgtmp, zgtmp, iglob_solid, loc_solid, ifseg, nglobslob, &
                  npointot, NGLLcube, NDIM)

  deallocate(ifseg)
  deallocate(loc_solid)
  deallocate(zgtmp)
  deallocate(sgtmp)

  if (dump_mesh_info_screen) write(6,*) 'NGLOBSLOB IS ' , NGLOBSLOB 

end subroutine define_global_slobal_numbering
!-------------------------------------------------------------------------


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

  !$ use omp_lib     
  implicit none

  integer, intent(in)               :: nspec2, npointot2, NGLLCUBE2, NDIM2
  double precision, intent(inout)   ::  xp(npointot2), yp(npointot2)
  integer, intent(out)              :: iglob2(npointot2), loc2(npointot2), nglob2
  logical, intent(out)              :: ifseg2(npointot2)

  integer :: ioffs(npointot2)

  integer :: ispec, i, j,nthreads, inttemp
  integer :: ieoff, ilocnum, nseg, ioff, iseg, ig
  double precision :: realtemp

  integer, dimension(:), allocatable :: ind, ninseg

  double precision, parameter :: SMALLVALTOL = 1.d-08

  ! establish initial pointers
  do ispec=1, nspec2
     ieoff = NGLLCUBE2 * (ispec - 1)
     do ilocnum=1, NGLLCUBE2
        loc2(ilocnum + ieoff) = ilocnum + ieoff
     enddo
  enddo

  ifseg2(:)=.false.

  allocate(ind(npointot2))
  allocate(ninseg(npointot2))

  ninseg = 0
  ind = 0

  nseg = 1
  ifseg2(1) = .true.
  ninseg(1) = npointot2

  do j=1, NDIM2

    ! sort within each segment
     ioff = 1

     if(j == 1) then
          call mergesort(xp, npointot2, yp, loc2)
     else
        !$ nthreads = min(OMP_get_max_threads(),8)
        !!$ print *, 'Using ', nthreads, ' threads!'
        ioffs(1) = 1
        do iseg=2,nseg
           ioffs(iseg) = ioffs(iseg-1) + ninseg(iseg)
           !write(23106,*) ninseg(iseg)
        end do
        !!!$omp parallel do shared(xp,yp,loc2,ninseg) private(ind,ioff)
        do iseg=1, nseg
           ioff = ioffs(iseg)
           ! First ordered by xp (j==1), then by yp (j==1)
           !call mergesort_serial(xp(ioffs(iseg):ioffs(iseg)+ninseg(iseg)-1), ninseg(iseg), &
           !               yp(ioffs(iseg):ioffs(iseg)+ninseg(iseg)-1), &
           !               loc2(ioffs(iseg):ioffs(iseg)+ninseg(iseg)-1))
           if (ninseg(iseg) == 2) then
               if (yp(ioff).gt.yp(ioff+1)) then
                   realtemp   = yp(ioff)
                   yp(ioff)   = yp(ioff+1)
                   yp(ioff+1) = realtemp
                   realtemp   = xp(ioff)
                   xp(ioff)   = xp(ioff+1)
                   xp(ioff+1) = realtemp
                   !yp(ioff:ioff+1) = yp([ioff+1, ioff])
                   !xp(ioff:ioff+1) = xp([ioff+1, ioff])
                   inttemp    = loc2(ioff)
                   loc2(ioff) = loc2(ioff+1)
                   loc2(ioff+1) = inttemp
                   !loc2(ioff:ioff+1) = loc2([ioff+1, ioff])
               end if
           elseif (ninseg(iseg) == 4) then
               call rank_y(yp(ioff), ind, 4)
               call swapall(loc2(ioff), xp(ioff), yp(ioff), ind, ninseg(iseg))
               !loc2(ioff:ioff+ninseg(iseg)-1) = loc2(ioff - 1 + ind(1:ninseg(iseg)))
               !xp(ioff:ioff+ninseg(iseg)-1) = xp(ioff - 1 + ind(1:ninseg(iseg)))
               !yp(ioff:ioff+ninseg(iseg)-1) = yp(ioff - 1 + ind(1:ninseg(iseg)))
           else
               call rank_y(yp(ioff), ind, ninseg(iseg))
               call swapall(loc2(ioff), xp(ioff), yp(ioff), ind, ninseg(iseg))
               !loc2(ioff:ioff+ninseg(iseg)-1) = loc2(ioff - 1 + ind(1:ninseg(iseg)))
               !xp(ioff:ioff+ninseg(iseg)-1) = xp(ioff - 1 + ind(1:ninseg(iseg)))
               !yp(ioff:ioff+ninseg(iseg)-1) = yp(ioff - 1 + ind(1:ninseg(iseg)))
           end if
        enddo
        !!!$omp end parallel do        
     endif
        

    ! check for jumps in current coordinate
    ! compare the coordinates of the points within a small tolerance
     if(j == 1) then
        do i=2, npointot2
           if (dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg2(i) = .true.
        enddo
     else
        do i=2,npointot2
           if (dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg2(i) = .true.
        enddo

     endif

    ! count up number of different segments
     nseg = 0
     do i=1, npointot2
        if(ifseg2(i)) then
           nseg = nseg + 1
           ninseg(nseg) = 1
        else
           ninseg(nseg) = ninseg(nseg) + 1
        endif

     enddo
  enddo ! NDIM2 loop

  !close(23106)
  deallocate(ind)
  deallocate(ninseg)

  ! assign global node numbers (now sorted lexicographically)
  ig = 0
  do i=1, npointot2
     if(ifseg2(i)) ig = ig + 1
     iglob2(loc2(i)) = ig
  enddo
  nglob2 = ig

end subroutine get_global
!-------------------------------------------------------------------------

! sorting routines put in same file to allow for inlining
!-------------------------------------------------------------------------
subroutine rank_y(A,IND,N)
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

end subroutine rank_y
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine swapall(IA,A,B,ind,n)
  !
  ! swap arrays IA, A, B and C according to addressing in array IND
  !
  implicit none

  integer, intent(in)             :: n
  integer, intent(in)             :: IND(n)
  integer, intent(inout)          :: IA(n)
  double precision,intent(inout)  :: A(n),B(n)
  double precision                :: W(n)
  integer                         :: IW(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)
  W(:) = B(:)

  do i=1,n
     IA(i)=IW(ind(i))
     A(i)=W(ind(i))
     B(i)=W(ind(i))
  enddo

end subroutine swapall
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine mergesort(A, N, D, E)
!$  use omp_lib     
    integer, intent(in)    :: N
    double precision, intent(inout) :: A(N), D(N)
    integer, intent(inout) :: E(N)
    integer                :: T(N)
    integer                :: nthreads

!$  call omp_set_nested(.true.)
!$  nthreads = min(OMP_get_max_threads(), 8)
!!$  print *, 'Using ', nthreads, ' threads!'

!$  call MergeSort_parallel(A,N,D,E,nthreads)
!$  if(.false.) then
    call MergeSort_serial(A, N, D, E)
!$  endif

end subroutine mergesort
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
recursive subroutine MergeSort_parallel(A, N, D, E, Threads)
 
   integer, intent(in)                  :: N, Threads
   double precision, intent(inout)      :: A(N), D(N)
   integer, intent(inout)               :: E(N)
   real(8), dimension(N)                :: Dtemp
   integer, dimension(N)                :: Etemp
   integer, dimension((N+1)/2)          :: INDA, INDB
   double precision, dimension((N+1)/2) :: T
 
   integer                              :: NA,NB, Vint, i
   real(8)                              :: V

!   print *, 'In MergeSort_parallel, N=', N, ', Threads=', Threads
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
         V = D(1)
         D(1) = D(2)
         D(2) = V
         Vint = E(1)
         E(1) = E(2)
         E(2) = Vint
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA

   if (Threads==1) then
       call MergeSort_serial(A(1)  , NA, D(1),   E(1)  )
       call MergeSort_serial(A(NA+1), NB, D(NA+1), E(NA+1))
   elseif(Threads>1) then
!$omp parallel sections shared(A, D, E) 
!$omp section
       call MergeSort_parallel(A(1),   NA, D(1),   E(1),   Threads/2)
!$omp section
       call MergeSort_parallel(A(NA+1), NB, D(NA+1), E(NA+1), Threads/2)
!$omp end parallel sections
   end if
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T, NA, A(NA+1), NB, A(1), INDA, INDB)
   
      Dtemp = D
      !D(INDA(1:NA)) = Dtemp(1:NA)
      !D(INDB(1:NB)) = Dtemp(NA+1:N)
      Etemp = E
      !Etemp(INDA(1:NA)) = Etemp(1:NA)
      !Etemp(INDB(1:NB)) = Etemp(NA+1:N)
      do i=1,NA
          D(inda(i)) = Dtemp(i)
          E(inda(i)) = Etemp(i)
      end do
      do i=1,NB
          D(indb(i)) = Dtemp(NA+i)
          E(indb(i)) = Etemp(NA+i)
      end do
      
   endif
   return
 
end subroutine MergeSort_parallel
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
recursive subroutine MergeSort_serial(A,N,D,E)
 
   integer, intent(in) :: N
   double precision, intent(inout) :: A(N), D(N)
   integer, intent(inout) :: E(N)
   real(8), dimension(N)                 :: Dtemp
   integer, dimension(N)                 :: Etemp
   integer, dimension((N+1)/2)           :: INDA,INDB
   double precision, dimension((N+1)/2)  :: T
 
   integer :: NA,NB,Vint,i
   double precision                      :: V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
         V = D(1)
         D(1) = D(2)
         D(2) = V
         Vint = E(1)
         E(1) = E(2)
         E(2) = Vint
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA

   call MergeSort_serial(A(1)  , NA, D(1),   E(1)  )
   call MergeSort_serial(A(NA+1), NB, D(NA+1), E(NA+1))
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T, NA, A(NA+1), NB, A, INDA, INDB)
   
      Dtemp = D
      Etemp = E
      do i=1,NA
          D(inda(i)) = Dtemp(i)
          E(inda(i)) = Etemp(i)
      end do
      do i=1,NB
          D(indb(i)) = Dtemp(NA+i)
          E(indb(i)) = Etemp(NA+i)
      end do
      !D(INDA(1:NA)) = Dtemp(1:NA)
      !D(INDB(1:NB)) = Dtemp(NA+1:N)
      !Etemp(INDA(1:NA)) = Etemp(1:NA)
      !Etemp(INDB(1:NB)) = Etemp(NA+1:N)
   endif
   return
 
end subroutine MergeSort_serial
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!@TODO merge is a bad name
subroutine Merge(A,NA,B,NB,C,INDA,INDB)
 
   integer, intent(in)           :: NA,NB         ! Normal usage: NA+NB = NC
   double precision, intent(in)  :: A(NA)        ! B overlays C(NA+1:NC)
   double precision, intent(in)  :: B(NB)
   double precision, intent(out) :: C(NA+NB)
   integer, intent(out)          :: INDA(NA)
   integer, intent(out)          :: INDB(NB)
 
   integer                       :: I,J,K

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         INDA(I) = K
         I = I+1
      else
         C(K) = B(J)
         INDB(J) = K
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      INDA(I) = K
      I = I + 1
      K = K + 1
   enddo
   do while (J <= NB)
      C(K) = B(J)
      INDB(J) = K
      J = J + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
!-------------------------------------------------------------------------

!=========================
end module numbering
!=========================
