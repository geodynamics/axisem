!========================
 module unit_stride_colloc
!========================

use global_parameters

implicit none
public
contains

!--------------------------------------------------------------------------
subroutine  collocate12s_sum_1d(a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6, &
                               a7,a8,a9,a10,b7,b8,b9,s1,s2,n)
          
integer, intent(in) :: n
real(kind=realkind), intent(in) :: a1(n),b1(n),a2(n),b2(n)
real(kind=realkind), intent(in) :: a3(n),b3(n),a4(n),b4(n),a5(n),b5(n)
real(kind=realkind), intent(in) :: a6(n),b6(n)
real(kind=realkind), intent(in) :: a7(n),a8(n),a9(n),a10(n),b7(n),b8(n),b9(n)
real(kind=realkind), intent(out) :: s1(n),s2(n)
real(kind=realkind) :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12
integer :: i     
  
  do i = 1, n
    c1 = a1(i) * b1(i)
    c2 = a2(i) * b2(i)
    c3 = a3(i) * b3(i)
    c4 = a4(i) * b4(i)
    c5 = a5(i) * b5(i)
    c6 = a6(i) * b6(i)

    c7 = a7(i) * b7(i)
    c8 = a8(i) * b8(i)
    c9 = a9(i) * b9(i)

    c10 =  a8(i) * b7(i)
    c11 = a10(i) * b8(i)
    c12 =  a9(i) * b9(i)

    s1(i) = c1 + c2 + two * (c3 + c4) + c5 + c6 + c7 + c8 + two * c9
    s2(i) = -two*(c1 + c2 + c5 + c6) - (c3 + c4) + c10 + c11 - c12
  end do 

end subroutine collocate12s_sum_1d
!==========================================================================

!--------------------------------------------------------------------------
subroutine collocate28s_sum_1d(m11,m21,m41,m12,m22,m42,m13,m23,m33,m43,&
                              X1,X2,X3,X4,X5,X6,u2,u3, &
                              z_eta_p,s_eta_mu,z_eta_m, &
                              z_xi_p,s_xi_mu,z_xi_m, &
                              s1p,s2p,s1m,s2m,n)

integer, intent(in) :: n
real(kind=realkind), intent(in), dimension(n)  :: m11,m21,m41,m12,m22,m42
real(kind=realkind), intent(in), dimension(n)  :: m13,m23,m33,m43
real(kind=realkind), intent(in), dimension(n)  :: X1,X2,X3,X4,X5,X6,u2,u3
real(kind=realkind), intent(in), dimension(n)  :: z_eta_p,s_eta_mu,z_eta_m
real(kind=realkind), intent(in), dimension(n)  :: z_xi_p,s_xi_mu,z_xi_m
real(kind=realkind), intent(out), dimension(n) :: s1p,s2p,s1m,s2m
real(kind=realkind) :: pm11,pm21,pm31,pm12,pm22,pm32
real(kind=realkind) :: p11,p21,p31,p41,p51,p12,p22,p32,p42,p52
real(kind=realkind) :: mi11,mi21,mi31,mi41,mi51,mi12,mi22,mi32,mi42,mi52
integer :: i
                                                                    
  do i = 1, n

    pm11 = m13(i) * X6(i)
    pm21 = m23(i) * X3(i)
    pm12 = m33(i) * X3(i)
    pm22 = m43(i) * X6(i)
    pm31 = s_eta_mu(i) * u3(i)
    pm32 = s_xi_mu(i) * u3(i)
    
    p11 = m11(i) * X4(i)
    p21 = m21(i) * X1(i)
    p31 = m12(i) * X5(i)
    p41 = m22(i) * X2(i)
    p51 = z_eta_p(i) * u2(i)

    mi11 = m11(i) * X5(i)
    mi21 = m21(i) * X2(i)
    mi31 = m12(i) * X4(i)
    mi41 = m22(i) * X1(i)
    mi51 = z_eta_m(i) * u2(i)
    
    p12 = m11(i) * X1(i)
    p22 = m41(i) * X4(i)
    p32 = m12(i) * X2(i)
    p42 = m42(i) * X5(i)
    p52 = z_xi_p(i) * u2(i)

    mi12 = m11(i) * X2(i)
    mi22 = m41(i) * X5(i)
    mi32 = m12(i) * X1(i)
    mi42 = m42(i) * X4(i)
    mi52 = z_xi_m(i)*u2(i)

    s1p(i) = pm11 + pm21 + pm31 +  p11 +  p21 +  p31 +  p41 +  p51
    s1m(i) = pm11 + pm21 - pm31 + mi11 + mi21 + mi31 + mi41 + mi51

    s2p(i) = pm12 + pm22 + pm32 + p12  + p22  + p32  + p42  + p52
    s2m(i) = pm12 + pm22 - pm32 + mi12 + mi22 + mi32 + mi42 + mi52

  end do 

end subroutine collocate28s_sum_1d
!==========================================================================

!========================
end module unit_stride_colloc
!========================
