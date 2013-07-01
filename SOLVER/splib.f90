!> Core of the spectral method. 
module splib
 
  use global_parameters
  
  public :: zelegl, zemngl2,                  &
            hn_jprime, lag_interp_deriv_wgl,  &
            get_welegl, get_welegl_axial
  private

  contains
 
!-----------------------------------------------------------------------------
!> This routine reorders array vin(n) in increasing order and
!! outputs array vout(n).
subroutine order(vin,vout,n)

  implicit none
  integer, intent(in)            :: n
  double precision, intent(in)   :: vin(n)
  double precision, intent(out)  :: vout(n)
  integer                        :: rankmax
  integer , dimension (n)        :: rank
  integer                        :: i, j

  rankmax = 1

  do i = 1, n

     rank(i) = 1

     do j = 1, n
        if((vin(i) > vin(j)) .and. (i /= j) ) rank(i) = rank(i) + 1
     end do

     rankmax = max(rank(i),rankmax)
     vout(rank(i)) = vin(i)

  end do

end subroutine order
!=============================================================================

!-----------------------------------------------------------------------------
!> Applies more robust formula to return
!! value of the derivative of the i-th Lagrangian interpolant
!! defined over the weighted GLL points computed at these
!! weighted GLL points.
subroutine lag_interp_deriv_wgl(dl,xi,i,N)
  
  implicit none
  integer, intent(in)           :: N, i
  double precision, intent(in)  :: xi(0:N)
  double precision, intent(out) :: dl(0:N)
  double precision              :: mn_xi_i, mnprime_xi_i
  double precision              :: mnprimeprime_xi_i
  double precision              :: mn_xi_j, mnprime_xi_j 
  double precision              :: mnprimeprime_xi_j
  double precision              :: DN
  integer                       :: j

  if ( i > N ) stop
  DN = dble(N)
  call vamnpo(N,xi(i),mn_xi_i,mnprime_xi_i,mnprimeprime_xi_i)

  if ( i == 0) then

     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)

        if (j == 0) &
                     dl(j) = -DN*(DN+two)/6.d0 
        if (j > 0 .and. j < N) &
                     dl(j) = two*((-one)**N)*mn_xi_j/((one+xi(j))*(DN+one))
        if (j == N) &
                     dl(j) = ((-one)**N)/(DN+one) 
     end do

  elseif (i == N) then

     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)
        if (j == 0) &
                     dl(j) = ((-one)**(N+1))*(DN+one)/4.d0
        if (j > 0 .and. j <  N) & 
                     dl(j) = -mn_xi_j/(one-xi(j))
        if (j == N) &
                     dl(j) = (DN*(DN+two)-one)/4.d0 
     end do

  else

     do j = 0, N
        call vamnpo(N,xi(j),mn_xi_j,mnprime_xi_j,mnprimeprime_xi_j)
        if (j == 0) &
                     dl(j) = ( ((-one)**(N+1)) * (DN+one) )&
                            /(two * mn_xi_i * (one + xi(i)))
        if (j > 0 .and. j < N .and. j /= i) &
                     dl(j) = ((xi(j)-xi(i))**(-1)) * mn_xi_j/mn_xi_i
        if (j > 0 .and. j < N .and. j == i) &
                     dl(j) = - half/(one + xi(j))
        if (j == N) &
                     dl(j) = (mn_xi_i*(one-xi(i)))**(-1)
     end do

  end if

end subroutine lag_interp_deriv_wgl
!=============================================================================

!-----------------------------------------------------------------------------
!> Compute the value of the derivative of the j-th Lagrange polynomial
!! of order N defined by the N+1 GLL points xi evaluated at these very
!! same N+1 GLL points. 
subroutine hn_jprime(xi,j,N,dhj)
 
  implicit none
  integer,intent(in)            :: N
  integer,intent(in)            :: j
  integer                       :: i
  double precision, intent(out) :: dhj(0:N)
  double precision              :: DX,D2X
  double precision              :: VN (0:N), QN(0:N)
  double precision              :: xi(0:N)
 
  dhj(:) = 0d0
  VN(:)= 0d0
  QN(:)= 0d0
  
  
  do i = 0, N
     call valepo(N, xi(i), VN(i), DX, D2X)
     if (i == j) QN(i) = 1d0
  end do
  
  call delegl(N, xi, VN, QN, dhj)
 
end subroutine hn_jprime
!=============================================================================

!-----------------------------------------------------------------------------
!> computes the nodes relative to the legendre gauss-lobatto formula
subroutine zelegl(n,et,vn)

  implicit double precision (a-h,o-z)
  integer, intent(in)            :: n !< Order of the formula
  double precision, intent(out)  :: ET(0:n) ! Vector of the nodes
  double precision, intent(out)  :: VN(0:n) ! Values of the Legendre polynomial at the nodes
  
  if (n  ==  0) return

  n2 = (n-1)/2
  sn = dfloat(2*n-4*n2-3)
  et(0) = -1.d0
  et(n) = 1.d0
  vn(0) = sn
  vn(n) = 1.d0
  
  if (n  ==  1) return

  et(n2+1) = 0.d0
  x = 0.d0
  call valepo(n,x,y,dy,d2y)
  
  vn(n2+1) = y
  
  if(n  ==  2) return

  c  = pi/dfloat(n)

  do i=1, n2
     etx = dcos(c*dfloat(i))
     do it=1, 8
        call valepo(n,etx,y,dy,d2y)
        etx = etx-dy/d2y
     end do
     et(i) = -etx
     et(n-i) = etx
     vn(i) = y*sn
     vn(n-i) = y
  end do

end subroutine zelegl
!=============================================================================

!-----------------------------------------------------------------------------
!> computes the nodes relative to the legendre gauss-lobatto formula
subroutine zelegl2(n,et)

  implicit double precision (a-h,o-z)
  integer, intent(in)            :: n !< Order of the formula
  double precision, intent(out)  :: ET(0:n) ! Vector of the nodes
  
  if (n  ==  0) return

  n2 = (n-1)/2
  sn = dfloat(2*n-4*n2-3)
  et(0) = -1.d0
  et(n) = 1.d0

  if (n  ==  1) return

  et(n2+1) = 0.d0
  x = 0.d0
  call valepo(n,x,y,dy,d2y)

  if (n .gt. 2) then

     c  = pi/dfloat(n)
      do i=1,n2
         etx = dcos(c*dfloat(i))
         do it=1,8
            call valepo(n,etx,y,dy,d2y)
            etx = etx-dy/d2y
         enddo   
         et(i) = -etx
         et(n-i) = etx
      enddo
  endif

  return

end subroutine zelegl2
!=============================================================================

!-----------------------------------------------------------------------------
!>   Computes the nodes relative to the modified legendre gauss-lobatto
!!   FORMULA along the s-axis
!!   Relies on computing the eigenvalues of tridiagonal matrix. 
!!   The nodes correspond to the second quadrature formula proposed
!!   by Azaiez et al.  
subroutine zemngl2(n,et)
  
  implicit double precision (a-h,o-z)
  integer, intent(in)              :: n       !< Order of the formula
  double precision, intent(out)    :: et(0:n) !< vector of the nodes, et(i), i=0,n.
  double precision, dimension(n-1) :: d, e
  integer                          :: i

  if (n  ==  0) return

     n2 = (n-1)/2
     et(0) = -1.d0
     et(n) = 1.d0
  if (n  ==  1) return

     et(n2+1) = 2d-1
     x = 2d-1
  if(n  ==  2) return

  ! Form the matrix diagonals and subdiagonals according to
  ! formulae page 109 of Azaiez, Bernardi, Dauge and Maday.

  do i = 1, n-1
     d(i) = three/(four*(dble(i)+half)*(dble(i)+three*half))
  end do

  do i = 1, n-2
     e(i+1) =   dsqrt(dble(i)*(dble(i)+three)) &
                      /(two*(dble(i)+three*half))
  end do

  ! Compute eigenvalues
  call tqli(d,e,n-1)

  ! Sort them in increasing order
  call order(d,e,n-1)

  ET(1:n-1) = e(1:n-1)

end subroutine zemngl2
!=============================================================================

!-----------------------------------------------------------------------------
!> This routines returns the eigenvalues of the tridiagonal matrix 
!! which diagonal and subdiagonal coefficients are contained in d(1:n) and
!! e(2:n) respectively. e(1) is free. The eigenvalues are returned in array d
subroutine tqli(d,e,n)

  implicit none
  integer, intent(in)             :: n
  double precision, intent(inout) :: d(n)
  double precision, intent(inout) :: e(n)
  integer                         :: i,iter,l,m
  double precision                :: b, c, dd, f, g, p, r, s

  do i = 2, n
    e(i-1) = e(i)
  end do

  e(n)=zero 
  do l=1,n
     iter=0
     iterate: do
     do m = l, n-1
       dd = abs(d(m))+abs(d(m+1))
       if (abs(e(m))+dd.eq.dd) exit
     end do
     if( m == l ) exit iterate
     if( iter == 30 ) stop 'too many iterations in tqli'
     iter=iter+1
     g = (d(l+1)-d(l))/(2.*e(l))
     r = pythag(g,one)
     g = d(m)-d(l)+e(l)/(g+sign(r,g))
     s = one
     c = one
     p = zero
     do i = m-1,l,-1
        f      = s*e(i)
        b      = c*e(i)
        r      = pythag(f,g)
        e(i+1) = r
        if(r == zero )then
           d(i+1) = d(i+1)-p
           e(m)   = zero 
           cycle iterate
        endif
        s      = f/r
        c      = g/r
        g      = d(i+1)-p
        r      = (d(i)-g)*s+2.*c*b
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
!=============================================================================

!-----------------------------------------------------------------------------
!> L2 norm of a and b  
double precision function pythag(a,b)

  implicit none
  double precision, intent(in) :: a, b
  double precision             :: absa,absb

  absa=dabs(a)
  absb=dabs(b)

  if(absa.gt.absb)then

     pythag=absa*sqrt(1.+(absb/absa)**2)

  else
     if(absb.eq.zero)then
        pythag=zero
     else
        pythag=absb*sqrt(1.+(absa/absb)**2)
     endif

  endif

end function pythag
!=============================================================================

!-------------------------------------------------------------------------
!>  computes the derivative of a polynomial at the legendre gauss-lobatto
!!  nodes from the values of the polynomial attained at the same points
subroutine delegl(n,et,vn,qn,dqn)

   implicit double precision (a-h,o-z)
   integer, intent(in)           ::  n        !< the degree of the polynomial
   double precision, intent(in)  ::  et(0:n)  !< vector of the nodes, et(i), i=0,n
   double precision, intent(out) ::  vn(0:n)  !< values of the legendre polynomial at the nodes, vn(i), i=0,n
   double precision, intent(out) ::  qn(0:n)  !< values of the polynomial at the nodes, qn(i), i=0,n
   double precision, intent(out) ::  dqn(0:n) !< derivatives of the polynomial at the nodes, dqz(i), i=0,n
       
   dqn(0) = 0.d0
   if (n .eq. 0) return

   do i=0,n
       su = 0.d0
       vi = vn(i)
       ei = et(i)
       do j=0,n 
           if (i .eq. j) cycle !goto 2
           vj = vn(j)
           ej = et(j)
           su = su+qn(j)/(vj*(ei-ej))
       enddo !2  continue   
       dqn(i) = vi*su
    enddo !1  continue   

    dn = dfloat(n)
    c  = .25d0 * dn * (dn+1.d0)
    dqn(0) = dqn(0) - c * qn(0)
    dqn(n) = dqn(n) + c * qn(n)

end subroutine delegl
!=============================================================================

!-----------------------------------------------------------------------------
!> computes the value of the legendre polynomial of degree n
!! and its first and second derivatives at a given point
subroutine valepo(n,x,y,dy,d2y)

  implicit double precision (a-h,o-z)

  integer, intent(in)           ::  n   !< degree of the polynomial
  double precision, intent(in)  ::  x   !< point in which the computation is performed
  double precision, intent(out) ::  y   !< value of the polynomial in x
  double precision, intent(out) ::  dy  !< value of the first derivative in x
  double precision, intent(out) ::  d2y !< value of the second derivative in x

  y   = 1.d0
  dy  = 0.d0
  d2y = 0.d0
  if(n == 0) return

  y   = x
  dy  = 1.d0
  d2y = 0.d0
  if(n == 1) return

  yp   = 1.d0
  dyp  = 0.d0
  d2yp = 0.d0
  do i = 2, n
     c1   = dfloat(i)
     c2   = 2.d0*c1-1.d0
     c4   = c1-1.d0
     ym   = y
     y    = (c2*x*y-c4*yp)/c1
     yp   = ym
     dym  = dy
     dy   = (c2*x*dy-c4*dyp+c2*yp)/c1
     dyp  = dym
     d2ym = d2y
     d2y  = (c2*x*d2y-c4*d2yp+2.d0*c2*dyp)/c1
     d2yp = d2ym
  enddo

end subroutine valepo
!=============================================================================

!-----------------------------------------------------------------------------
!> This routine computes the N+1 weights associated with the
!! Gauss-Lobatto-Legendre quadrature formula of order N.
subroutine get_welegl(N,xi,wt)

  implicit none
  integer           :: N
  double precision  :: xi(0:N), wt(0:N)
  integer           :: j
  double precision  :: y,dy,d2y,fact 

  fact = 2.0d0/(dble(N*(N+1)))

  wt(:) = 0.0

  do j = 0, N
     call valepo(N,xi(j),y,dy,d2y)
     wt(j) =  fact*y**(-2)
  end do

end subroutine get_welegl
!=============================================================================

!-----------------------------------------------------------------------------
!> This routine computes the N+1 weights associated with the
!! Gauss-Lobatto-Legendre quadrature formula of order N that one 
!! to apply for elements having a non-zero intersection with the
!! axis of symmetry of the Earth.
subroutine get_welegl_axial(N,xi,wt,iflag)

  implicit none
  integer, intent(in)           :: N       !< order of GLL quadrature formula
  integer, intent(in)           :: iflag   !< Selector for quadrature formulae proposed 
                                           !! by Bernardi et al.
                                           !! iflag = 2 : Second formula 
                                           !!             Formula : (VI.1.12), page 104             
                                           !! iflag = 3 : Third formula
                                           !!             Formula : (VI.1.20), page 107            
  double precision, intent(in)  :: xi(0:N) !< Support points
  double precision, intent(out) :: wt(0:N) !< Weighting factor at support points
  integer                       :: j
  double precision              :: y, dy, d2y, fact

  wt(:) = 0.0 

  if (iflag == 2 ) then 

     fact = 4d0 / dble(N*(N+2)) !four/(dble(N)*dble(N+2))
     do j = 0, N
        call vamnpo(N, xi(j), y, dy, d2y)
        wt(j) =  fact / (y*y)
        if (j == 0) wt(j) = 2.0 * wt(j)
     end do

  elseif ( iflag == 3 ) then 

     fact = 1.0 / dble((N+1)*(N+1))
     do j = 0, N
        call valepo(N, xi(j), y, dy, d2y)
        wt(j) = (fact * (1+xi(j))**2) / (y*y) 
     end do  

  end if

end subroutine get_welegl_axial
!=============================================================================

!-----------------------------------------------------------------------------
!>   Computes the value of the "cylindrical" polynomial
!!   m_n = (l_n + l_{n+1})/(1+x) of degree n
!!   and its first and second derivatives at a given point
!!
!!   implemented after bernardi et al., page 57, eq. (iii.1.10)
subroutine vamnpo(n,x,y,dy,d2y)
  
  implicit double precision (a-h,o-z)
  integer, intent(in)           :: n   !< degree of the polynomial 
  double precision, intent(in)  :: x   !< point in which the computation is performed
  double precision, intent(out) :: y   !< value of the polynomial in x
  double precision, intent(out) :: dy  !< value of the first derivative in x
  double precision, intent(out) :: d2y !< value of the second derivative in x
  
  
   y   = 1.d0
   dy  = 0.d0
   d2y = 0.d0
  if (n  ==  0) return

   y   = 1.5d0*x - 0.5d0 !half*(three*x-one)
   dy  = 1.5d0 !half*three
   d2y = 0.d0
  if(n  ==  1) return

   yp   = 1.d0
   dyp  = 0.d0
   d2yp = 0.d0
  do i=2,n
      c1 = dble(i-1)
      ym = y
       y = (x-one/((2*c1+one)*(2*c1+three)) ) * y &
          - (c1/(two*c1+one))*yp
       y = (two*c1+three)*y/(c1+two)
      yp = ym
     dym = dy
      dy = (x-one/((2*c1+one)*(2*c1+three)) ) * dy &
           +yp - (c1/(two*c1+one))*dyp
      dy = (two*c1+three)*dy/(c1+two)
     dyp = dym
    d2ym = d2y
    d2y  = two*dyp + (x-one/((2*c1+one)*(2*c1+three)) ) * d2y &
           - (c1/(two*c1+one))*d2yp
    d2y  = (two*c1+three)*d2y/(c1+two)
    d2yp = d2ym
  end do
  
end subroutine vamnpo
!=============================================================================

!-----------------------------------------------------------------------------
!> This routine computes the Lagrange interpolated value y at point x
!! associated to the function defined by the n values ya at n distinct points
!! xa. dy is the estimate of the error made on the interpolation.
subroutine polint(xa,ya,n,x,y,dy)

  integer ::  n
  double precision ::  dy,x,y,xa(n),ya(n)
  integer ::  i,m,ns
  double precision :: den,dif,dift,ho,hp,w,c(n),d(n)


  ns=1
  dif=abs(x-xa(1))

  do i=1,n

     dift=abs(x-xa(i))
     if (dift < dif) then
       ns=i    
       dif=dift
     endif    
     c(i)=ya(i)       
     d(i)=ya(i)       
 
  end do
 
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den == zero) stop 'failure in polint'
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

end subroutine polint
!----------------------------------------------------------------------------

end module splib
