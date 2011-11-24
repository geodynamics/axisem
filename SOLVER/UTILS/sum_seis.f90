subroutine sum_seis(nsim,simdir,src_type,appmynum,mij,magnitude,nt,colat,lon,time,seis_sglcomp,seis)
    
integer, intent(in) :: nt,nsim
real, intent(in) :: lon,colat
real, dimension(nt),intent(in) :: time
character(len=34), intent(in) :: simdir(nsim)
character(len=34), intent(in) :: src_type(nsim,2)
character(len=4), intent(in) :: appmynum
real, intent(in) :: mij(6),magnitude(nsim)
real, intent(in), dimension(nt,3,nsim) :: seis_sglcomp
real, dimension(nt,3) :: seisxx,seiszz,seisxz,seisxy
real, dimension(nt,3),intent(out) :: seis

seis(:,1) = seis(:,1)+ prefact*mij *seisin(:,1) 


seis(:,1) = mij(2)*(seisxx(:,1) + seisxy(:,1)*cos(2.*lon) ) +& !Mxx
                   mij(3)*(seisxx(:,1) - seisxy(:,1)*cos(2.*lon) ) + & !Myy
                   mij(1)*seiszz(:,1) + & !Mzz
                   mij(4)*seisxz(:,1)*cos(lon) +& ! Mxz
                   mij(5)*seisxz(:,1)*sin(lon) +& !Myz
                   mij(6)*seisxy(:,1)*sin(2.*lon) ! Mxy

seis(:,2) = mij(2)*(seisxx(:,2) + seisxy(:,2)*sin(2.*lon) ) +& !Mxx
                   mij(3)*(seisxx(:,2) - seisxy(:,2)*sin(2.*lon) ) + & !Myy
                   mij(1)*seiszz(:,2) + & !Mzz
                   mij(4)*seisxz(:,2)*sin(lon) +& ! Mxz
                   mij(5)*seisxz(:,2)*cos(lon) +& !Myz
                   mij(6)*seisxy(:,2)*cos(2.*lon) ! Mxy

seis(:,3) = mij(2)*(seisxx(:,3) + seisxy(:,3)*cos(2.*lon) ) +& !Mxx
                   mij(3)*(seisxx(:,3) - seisxy(:,3)*cos(2.*lon) ) + & !Myy
                   mij(1)*seiszz(:,3) + & !Mzz
                   mij(4)*seisxz(:,3)*cos(lon) +& ! Mxz
                   mij(5)*seisxz(:,3)*sin(lon) +& !Myz
                   mij(6)*seisxy(:,3)*sin(2.*lon) ! Mxy

end subroutine sum_seis
!=============================================================================
