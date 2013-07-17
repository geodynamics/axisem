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

!==============
  module splib
!==============
!
! Core of the spectral method. 
!
  use global_parameters
!
  public :: zemngr,valepo,inlegl,zelegl,zemngl2,               &
            delegl,hn_j,hn_jprime,lag_interp, lag_interp_deriv,&
            lag_interp_deriv_wgl,                              &
            get_welegl,get_welegl_axial,vamnpo
  private
!
  contains
!
!dk order---------------------------------------------------------------
  subroutine order(vin,vout,n)
!
!       This routine reorders array vin(n) in increasing order and
! outputs array vout(n).
!
  integer :: n
  double precision, dimension(n) :: vin,vout
  integer :: rankmax
  integer , dimension (n) :: rank
  integer :: i, j

  rankmax = 1

  do i = 1, n
!
     rank(i) = 1
!
     do j = 1, n
        if((vin(i) > vin(j)) .and. (i /= j) ) rank(i) = rank(i) + 1
     end do
!
     rankmax = max(rank(i),rankmax)
     vout(rank(i)) = vin(i)
!
  end do
!
  end subroutine order
!-------------------------------------------------------------------------
!
!--------------------------------------------
  double precision function hn_j(zeta,xi,j,N)
!--------------------------------------------
! Compute the value of the j-th Lagrange polynomial 
! of order N defined by the N+1 GLL points xi at a
! point zeta in [-1,1].
!
  implicit none
!
  integer :: j, N
  integer :: i 
  double precision :: zeta, DX,D2X
  double precision :: VN (0:N), QN(0:N)
  double precision :: xi(0:N)
!
  hn_j = 0.d0
  VN(:)= 0.d0
  QN(:)= 0.d0
!
! 
 do i = 0, N
    call VALEPO(N,xi(i), VN(i), DX,D2X)
    if (i == j) QN(i) = 1.d0
 end do
!
 call INLEGL(N,xi,VN,QN,zeta,hn_j)
! 
 end function hn_j
!------------------
!
!dk lag_interp------------------------------------------------ 
  double precision function lag_interp(zeta,xi,j,N,ibeg,iend)   
!
! Computes the value of the j-th Lagrange polynomial
! of order N defined by the N+1 points xi at a
! point zeta in [-1,1].
!
  implicit none
!
  integer :: j, N,ibeg,iend
  integer :: i
  double precision :: zeta, adist
  double precision :: xi(0:N)
!
  lag_interp = 1.d0
!
!
  do i = ibeg,iend
     if ( i /= j) then 
        adist = 1./(xi(j)-xi(i))
        lag_interp = lag_interp*(zeta-xi(i))*adist
     end if
  end do
!
  end function lag_interp
!-------------------------
!
!dk lag_interp_deriv---------------------------------------------------
  subroutine lag_interp_deriv(df,xi,j,N)
!
! Computes the value of the derivative of the j-th Lagrange polynomial
! of order N defined by the N+1 points xi at each of these collocation
! points xi. 
!
  implicit none
!
  integer          :: j, N
  double precision :: xi(0:N), df(0:N)
  integer          :: i,ixi,l
  double precision :: dist,fact,prod(0:N,0:N)
!
  df(:) = 0.
  fact  = 1.
!
  do i = 0,N
     if ( i /= j) then
        dist = (xi(j)-xi(i))
        fact = fact*dist 
     end if
  end do

  do ixi = 0, N
     do i = 0, N
        prod(i,ixi) = 1.
        do l = 0, N
           if ( l /= i .and. l /= j) prod(i,ixi) = prod(i,ixi)*(xi(ixi)-xi(l))
        end do
     end do
  end do

  do  ixi = 0, N
     do i = 0,N
        if( i /= j) df(ixi) = df(ixi) + prod(i,ixi)
     end do
  end do

  df(:) = df(:)/fact
!
!---------------------------------
  end subroutine lag_interp_deriv
!---------------------------------
!
!dk lag_interp_deriv_wgl---------------------------------
  subroutine lag_interp_deriv_wgl(dl,xi,i,N)
!
!	Applies more robust formula to return
! value of the derivative of the i-th Lagrangian interpolant
! defined over the weighted GLL points computed at these
! weighted GLL points.
!
! Alex 08282003: New (robust) routine. 
! Apply formula (A11) of appendix A of JCP draft.
! 
!
  implicit none
  integer :: N, j
  double precision :: dl(0:N),xi(0:N)
  double precision :: mn_xi_i,mnprime_xi_i,mnprimeprime_xi_i
  double precision :: mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j
  double precision :: DN
  integer :: i
!
  if ( i > N ) stop
  DN = dble(N)
  call vamnpo(N,xi(i),mn_xi_i,mnprime_xi_i,mnprimeprime_xi_i)
!
  if ( i == 0) then

     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)

        if (j == 0) &
                     dl(j) = -DN*(DN+two)/6.d0 
        if (j > 0 .and. j < N) &
                     dl(j) = two*((-1.)**N)*mn_xi_j/((1.+xi(j))*(DN+1.))
        if (j == N) &
                     dl(j) = ((-1.)**N)/(DN+1.) 
     end do
!
  elseif (i == N) then
!
     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)
        if (j == 0) &
                     dl(j) = ((-1.)**(N+1))*(DN+1.)/4.d0
        if (j > 0 .and. j <  N) & 
                     dl(j) = -mn_xi_j/(1.-xi(j))
        if (j == N) &
                     dl(j) = (DN*(DN+two)-1.)/4.d0 
     end do
!
  else
!
     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)
        if (j == 0) &
                     dl(j) = ( ((-1.)**(N+1)) * (DN+1.) )&
                            /(two * mn_xi_i * (1. + xi(i)))
        if (j > 0 .and. j < N .and. j /= i) &
                     dl(j) = ((xi(j)-xi(i))**(-1)) * mn_xi_j/mn_xi_i
        if (j > 0 .and. j < N .and. j == i) &
                     dl(j) = - half/(1. + xi(j))
        if (j == N) &
                     dl(j) = (mn_xi_i*(1.-xi(i)))**(-1)
     end do
!
  end if

  end subroutine lag_interp_deriv_wgl
!--------------------------------------------------------
!
!--------------------------------
 subroutine hn_jprime(xi,j,N,dhj)
!--------------------------------
!Compute the value of the derivative of the j-th Lagrange polynomial
!of order N defined by the N+1 GLL points xi evaluated at these very
!same N+1 GLL points. 
 implicit none
!
 integer :: j, N
 integer :: i
 double precision :: DX,D2X
 double precision :: VN (0:N), QN(0:N)
 double precision :: xi(0:N), dhj(0:N)
!
 dhj(:) = 0.d0
 VN(:)= 0.d0
 QN(:)= 0.d0
!
!
 do i = 0, N
    call VALEPO(N,xi(i), VN(i), DX,D2X)
    if (i == j) QN(i) = 1.d0
 end do
!
 call DELEGL(N,xi,VN,QN,dhj)
!
!------------------------
 end subroutine hn_jprime
!------------------------
!
!dk zelegl
      SUBROUTINE ZELEGL(N,ET,VN)
!********************************************************************
!   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
!   N  = ORDER OF THE FORMULA
!   ET = VECTOR OF THE NODES, ET(I), I=0,N
!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
!********************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION ET(0:*), VN(0:*)
  IF (N  ==  0) RETURN

     N2 = (N-1)/2
     SN = DFLOAT(2*N-4*N2-3)
     ET(0) = -1.D0
     ET(N) = 1.D0
     VN(0) = SN
     VN(N) = 1.D0
  IF (N  ==  1) RETURN

     ET(N2+1) = 0.D0
     X = 0.D0
  CALL VALEPO(N,X,Y,DY,D2Y)
     VN(N2+1) = Y
  IF(N  ==  2) RETURN

     C  = PI/DFLOAT(N)
  DO 1 I=1,N2
     ETX = DCOS(C*DFLOAT(I))
  DO 2 IT=1,8
  CALL VALEPO(N,ETX,Y,DY,D2Y)
     ETX = ETX-DY/D2Y
2 CONTINUE   
     ET(I) = -ETX
     ET(N-I) = ETX
     VN(I) = Y*SN
     VN(N-I) = Y
1 CONTINUE   

  RETURN     
!-----------------------
  end subroutine zelegl
!-----------------------
!
!dk zelegl2------------------------------------------------------------------------------------
      SUBROUTINE ZELEGL2(N,ET)
!********************************************************************
!   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
!   N  = ORDER OF THE FORMULA
!   ET = VECTOR OF THE NODES, ET(I), I=0,N
!********************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION ET(0:*)
  IF (N  ==  0) RETURN

     N2 = (N-1)/2
     SN = DFLOAT(2*N-4*N2-3)
     ET(0) = -1.D0
     ET(N) = 1.D0
  IF (N  ==  1) RETURN

     ET(N2+1) = 0.D0
     X = 0.D0
  CALL VALEPO(N,X,Y,DY,D2Y)
  IF(N  ==  2) RETURN

     C  = PI/DFLOAT(N)
  DO 1 I=1,N2
     ETX = DCOS(C*DFLOAT(I))
  DO 2 IT=1,8
  CALL VALEPO(N,ETX,Y,DY,D2Y)
     ETX = ETX-DY/D2Y
2 CONTINUE
     ET(I) = -ETX
     ET(N-I) = ETX
1 CONTINUE

  RETURN
!
  end subroutine zelegl2
!---------------------------------------------------------------------------------
!
!dk zemngl2----------------------------------------------------------
      SUBROUTINE ZEmnGL2(N,ET)
!********************************************************************
!
!   COMPUTES THE NODES RELATIVE TO THE modified LEGENDRE GAUSS-LOBATTO
!   FORMULA along the s-axis
!   N  = ORDER OF THE FORMULA
!   ET = VECTOR OF THE NODES, ET(I), I=0,N. 
!   Relies on computing the eigenvalues of tridiagonal matrix. 
!   The nodes correspond to the second quadrature formula proposed
!   by Azaiez et al.  
!
!********************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION ET(0:*)
  double precision, dimension(n-1) :: d, e
  integer :: i

  IF (N  ==  0) RETURN

     N2 = (N-1)/2
     ET(0) = -1.D0
     ET(N) = 1.D0
  IF (N  ==  1) RETURN

     ET(N2+1) = 2D-1
     X = 2d-1
  IF(N  ==  2) RETURN
!
! Form the matrix diagonals and subdiagonals according to
! formulae page 109 of Azaiez, Bernardi, Dauge and Maday.
!
  do i = 1, n-1
     d(i) = three/(four*(dble(i)+half)*(dble(i)+three*half))
  end do

  do i = 1, n-2
     e(i+1) =   dsqrt(dble(i)*(dble(i)+three)) &
                      /(two*(dble(i)+three*half))
  end do
!
! Compute eigenvalues
  call tqli(d,e,n-1)
!
! Sort them in increasing order
  call order(d,e,n-1)
!
  ET(1:n-1) = e(1:n-1)
!
!
  end subroutine zemngl2
!------------------------
!
!dk zemngr----------------------------------------------------------
      SUBROUTINE ZEmnGr(N,ET)
!********************************************************************
!
!   COMPUTES THE NODES RELATIVE TO THE pseudo-Gauss Radau
!   FORMULA along the s-axis
!   N  = ORDER OF THE FORMULA
!   ET = VECTOR OF THE NODES, ET(I), I=0,N.
!   Relies on computing the eigenvalues of tridiagonal matrix.
!   The nodes correspond to the third quadrature formula proposed
!   by Azaiez et al.  
!
!********************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION ET(0:*)
! double precision, dimension(n-1) :: d, e  ! ALEX CHANGE 10102001
  double precision, dimension(n) :: d, e
  integer :: i

  IF (N  ==  0) RETURN

  ET(N) = 1.D0
!
! Form the matrix diagonals and subdiagonals according to
! formulae page 109 of Azaiez, Bernardi, Dauge and Maday.
!
  do i = 1, n
     d(i) = -1.d0/(four*(dble(i)-half)*(dble(i)+half))
  end do

  do i = 1, n-1
     e(i+1) =   dsqrt(dble(i)*(dble(i)+1.)) &
                      /(two*(dble(i)+half))
  end do
!
! Compute eigenvalues
  call tqli(d,e,n)
!
! Sort them in increasing order
  call order(d,e,n)
!
  ET(0:n-1) = e(1:n)
!
  end subroutine zemngr
!-------------------------------------------------------------------------
!
!dk tqli------------------------------------------------------------------
  subroutine tqli(d,e,n)
!
!  	This routines returns the eigenvalues of the tridiagonal matrix 
!which diagonal and subdiagonal coefficients are contained in d(1:n) and
!e(2:n) respectively. e(1) is free. The eigenvalues are returned in array
!d. 
!
  implicit none
!
  integer :: n
  double precision :: d(n),e(n)
  integer ::  i,iter,l,m
  double precision :: b,c,dd,f,g,p,r,s,one_dble
!
  one_dble=one

  do i = 2, n
    e(i-1) = e(i)
  end do

  e(n)=0. 
  do l=1,n
     iter=0
     iterate: do
     do m = l, n-1
       dd = dabs(d(m))+dabs(d(m+1))
       if (dabs(e(m))+dd.eq.dd) exit
     end do
     if( m == l ) exit iterate
     if( iter == 30 ) stop 'too many iterations in tqli'
     iter=iter+1
     g = (d(l+1)-d(l))/(2.d0*e(l))
     r = pythag(g,one_dble)
     g = d(m)-d(l)+e(l)/(g+sign(r,g))
     s = one
     c = one
     p = zero
     do i = m-1,l,-1
        f      = s*e(i)
        b      = c*e(i)
        r      = pythag(f,g)
        e(i+1) = r
        if(r == 0.d0 )then
           d(i+1) = d(i+1)-p
           e(m)   = zero
           cycle iterate
        endif
        s      = f/r
        c      = g/r
        g      = d(i+1)-p
        r      = (d(i)-g)*s+2.d0*c*b
        p      = s*r
        d(i+1) = g+p
        g      = c*r-b
     end do
     d(l) = d(l)-p
     e(l) = g
     e(m) = zero
     end do iterate
  end do

  end subroutine tqli
!------------------------------------------------------------------------
!dk pythag---------------------------------------------------------------
  double precision FUNCTION pythag(a,b)
!
  implicit none  
!
  double precision :: a,b
  double precision :: absa,absb
!
  absa=dabs(a)
  absb=dabs(b)
!
  if(absa.gt.absb)then

     pythag=absa*dsqrt(1.+(absb/absa)**2)

  else
     if(absb.eq.0.)then
        pythag=0.
     else
        pythag=absb*dsqrt(1.+(absa/absb)**2)
     endif

  endif

  end function pythag
!-------------------------------------------------------------------------
  SUBROUTINE DELEGL(N,ET,VN,QN,DQN)
!***********************************************************************
!  COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LEGENDRE GAUSS-LOBATTO
!  NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
!   N   = THE DEGREE OF THE POLYNOMIAL
!   ET  = VECTOR OF THE NODES, ET(I), I=0,N
!   VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
!   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
!   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N
!***********************************************************************
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   DIMENSION ET(0:*), VN(0:*), QN(0:*), DQN(0:*)
       DQN(0) = 0.D0
   IF (N .EQ. 0) RETURN

   DO 1 I=0,N
       SU = 0.D0
       VI = VN(I)
       EI = ET(I)
   DO 2 J=0,N 
   IF (I .EQ. J) GOTO 2
       VJ = VN(J)
       EJ = ET(J)
       SU = SU+QN(J)/(VJ*(EI-EJ))
2  CONTINUE   
       DQN(I) = VI*SU
1  CONTINUE   

       DN = DFLOAT(N)
       C  = .25D0*DN*(DN+1.D0)
       DQN(0) = DQN(0)-C*QN(0)
       DQN(N) = DQN(N)+C*QN(N)

   RETURN     
!
 end subroutine delegl
!--------------------------------
!
!dk valepo------------------------------------------------------
  SUBROUTINE VALEPO(N,X,Y,DY,D2Y)
  !*************************************************************
  !   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N
  !   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
  !   N  = DEGREE OF THE POLYNOMIAL
  !   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
  !   Y  = VALUE OF THE POLYNOMIAL IN X
  !   DY = VALUE OF THE FIRST DERIVATIVE IN X
  !   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
  !*************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)

   Y   = 1.D0
   DY  = 0.D0
   D2Y = 0.D0
  IF (N  ==  0) RETURN

   Y   = X
   DY  = 1.D0
   D2Y = 0.D0
  IF(N  ==  1) RETURN

   YP   = 1.D0
   DYP  = 0.D0
   D2YP = 0.D0
  DO 1 I=2,N
   C1 = DFLOAT(I)
   C2 = 2.D0*C1-1.D0
   C4 = C1-1.D0
   YM = Y
   Y  = (C2*X*Y-C4*YP)/C1
   YP = YM
   DYM  = DY
   DY   = (C2*X*DY-C4*DYP+C2*YP)/C1
   DYP  = DYM
   D2YM = D2Y
   D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
   D2YP = D2YM
  1     CONTINUE

  RETURN

!        
  end subroutine valepo
!-----------------------
!
  SUBROUTINE INLEGL(N,ET,VN,QN,X,QX)
  !*********************************************************************
  !   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
  !   BY THE VALUES ATTAINED AT THE LEGENDRE GAUSS-LOBATTO NODES
  !   N  = THE DEGREE OF THE POLYNOMIAL
  !   ET = VECTOR OF THE NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
  !   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
  !   QX = VALUE OF THE POLYNOMIAL IN X
  !*********************************************************************
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   DIMENSION ET(0:*), VN(0:*), QN(0:*)
   IF (N == 0) RETURN

       EPS = 1.D-14
   CALL VALEPO(N,X,Y,DY,D2Y)
       DN = DFLOAT(N)
       C  = 1.D0/(DN*(DN+1.D0))
       QX = QN(0)*C*DY*(X-1.D0)/VN(0)
       QX = QX+QN(N)*C*DY*(X+1.D0)/VN(N)
   IF (N .EQ. 1) RETURN

   DO 1 J=1,N-1
       ED = X-ET(J)
   IF(DABS(ED) < EPS) THEN
       QX = QN(J)
   RETURN     
   ELSE
       QX = QX+QN(J)*C*DY*(X*X-1.D0)/(VN(J)*ED)
   ENDIF      
1  CONTINUE   

   RETURN     
!----------------------  
  end subroutine inlegl
!----------------------
!
!dk get_welegl------------------------------------------------------------------
  subroutine get_welegl(N,xi,wt)
!	
!	This routine computes the N+1 weights associated with the
!Gauss-Lobatto-Legendre quadrature formula of order N.
!
  implicit none
!
  integer :: N
  double precision ::  xi(0:N),wt(0:N)
  integer :: j
  double precision :: y,dy,d2y,fact 
!
  fact = two/(dble(N)*dble(N+1))
!
  wt(:) = 0.d0
!
  do j = 0, N
     call VALEPO(N,xi(j),y,dy,d2y)
     wt(j) =  fact*y**(-2)
  end do
!
  end subroutine get_welegl
!-------------------------------------------------------------------------------
!
!dk get_welegl_axial------------------------------------------------------------
  subroutine get_welegl_axial(N,xi,wt,iflag)
!
!       This routine computes the N+1 weights associated with the
!Gauss-Lobatto-Legendre quadrature formula of order N that one 
!to apply for elements having a non-zero intersection with the
!axis of symmetry of the Earth.
!
!
! iflag = 2 : Second  quadrature formula proposed by Bernardi et al.
!             Formula : (VI.1.12), page 104             
! iflag = 3 : Third   quadrature formula proposed by Bernardi et al.
!             Formula : (VI.1.20), page 107            
!
  implicit none
!
  integer :: N,iflag
  double precision ::  xi(0:N),wt(0:N)
  integer :: j
  double precision :: y,dy,d2y,fact
!
  wt(:) = 0.d0
!
       if (iflag == 2 ) then 

          fact = four/(dble(N)*dble(N+2))
          do j = 0, N
             call VAmnPO(N,xi(j),y,dy,d2y)
             wt(j) =  fact*y**(-2)
             if (j == 0) wt(j) = two*wt(j)
          end do

       elseif ( iflag == 3 ) then 

          fact = 1./(dble(N+1)*dble(N+1))
          do j = 0, N
             call valepo(N,xi(j),y,dy,d2y)
             wt(j) = (fact*(1+xi(j))**2)*y**(-2)
          end do  

       end if
!
  end subroutine get_welegl_axial
!-------------------------------------------------------------------------------
!
!dk vamnpo---------------------------------------------------------
  SUBROUTINE VAMNPO(N,X,Y,DY,D2Y)
  !*************************************************************
  !   COMPUTES THE VALUE OF the "cylindrical" polynomial
  !   M_n = (L_n + L_{n+1})/(1+x) OF DEGREE N
  !   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
  !   N  = DEGREE OF THE POLYNOMIAL 
  !   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
  !   Y  = VALUE OF THE POLYNOMIAL IN X
  !   DY = VALUE OF THE FIRST DERIVATIVE IN X
  !   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
  !
  !   Implemented after Bernardi et al., page 57, eq. (III.1.10)
  ! 
  !*************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
   Y   = 1.D0
   DY  = 0.D0
   D2Y = 0.D0
  IF (N  ==  0) RETURN

   Y   = half*(three*X-1.)
   DY  = half*three
   D2Y = 0.D0
  IF(N  ==  1) RETURN

   YP   = 1.D0
   DYP  = 0.D0
   D2YP = 0.D0
  do i=2,N
      C1 = dble(I-1)
      YM = Y
       Y = (X-1./((2*C1+1.)*(2*C1+three)) ) * Y &
          - (C1/(two*C1+1.))*YP
       Y = (two*C1+three)*Y/(c1+two)
      YP = YM
     DYM = DY
      DY = (X-1./((2*C1+1.)*(2*C1+three)) ) * DY &
           +YP - (C1/(two*C1+1.))*DYP
      DY = (two*C1+three)*DY/(c1+two)
     DYP = DYM
    D2YM = D2Y
    D2y  = two*dyp + (X-1./((2*C1+1.)*(2*C1+three)) ) * D2Y &
           - (C1/(two*C1+1.))*D2YP
    D2Y  = (two*C1+three)*D2Y/(c1+two)
    D2YP = D2YM
  end do
!
  RETURN
!
!-----------------------
  end subroutine vamnpo
!-----------------------
!
!dk welegl----------------------------------------
      SUBROUTINE WELEGL(N,ET,VN,WT)
!**********************************************************************
!   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
!   N  = ORDER OF THE FORMULA
!   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
!   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
!**********************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION ET(0:*), VN(0:*), WT(0:*)
  IF (N .EQ. 0) RETURN

      N2 = (N-1)/2
      DN = DFLOAT(N)
      C  = 2.D0/(DN*(DN+1.D0))
  DO 1 I=0,N2
      X = ET(I)
      Y = VN(I)
      WTX = C/(Y*Y)
      WT(I) = WTX
      WT(N-I) = WTX
1 CONTINUE

  IF(N-1 .EQ. 2*N2) RETURN
      X = 0.D0
      Y = VN(N2+1)
      WT(N2+1) = C/(Y*Y)

  RETURN
 end subroutine welegl
!-------------------------------------------------
!
!dk polint-----------------------------------------------------------------
  subroutine polint(xa,ya,n,x,y,dy)
!
! 	This routine computes the Lagrange interpolated value y at point x
! associated to the function defined by the n values ya at n distinct points
! xa. dy is the estimate of the error made on the interpolation.
!

  implicit none
!
  integer ::  n
  double precision ::  dy,x,y,xa(n),ya(n)
  integer ::  i,m,ns
  double precision :: den,dif,dift,ho,hp,w,c(n),d(n)
!
!
  ns=1
  dif=dabs(x-xa(1))
!
  do i=1,n
!
     dift=dabs(x-xa(i))
     if (dift < dif) then
       ns=i    
       dif=dift
     endif    
     c(i)=ya(i)       
     d(i)=ya(i)       
! 
  end do
 
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den == 0.) stop 'failure in polint'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     end do
     if (2*ns < n-m)then
       dy=c(ns+1)
     else
       dy=d(ns)
       ns=ns-1
     endif
     y=y+dy

  end do
!
  end subroutine polint
!----------------------------------------------------------------------------
!
!========================
 end module splib
!========================
