!========================
module background_models
!========================
!
! This module is identical in the mesher and solver.
! Function "velocity" retrieves the velocity (or density) of the given 
! background model (see different cases and respective routines below)
! for radius r0 given its subdomain identification number idom 
! (determined respectively in mesher and solver). 
!
! When adding new background models, one needs to define them both in terms 
! of subdomains for this module (e.g. as a polynomial), and in terms of the 
! discontinuities and their above/below elastic values for the mesher only 
! (see module model_discontinuities).

implicit none

public :: velocity, model_is_ani
private
contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------------
double precision function velocity(r0,param,idom,bkgrdmodel2,lfbkgrdmodel2)
!
! Wrapper function to call velocities upon different background models 
! for a given radius r0 [m], parameter type param (rho,vs,vp) and idom
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in)   :: r0
integer, intent(in)            :: idom
character(len=100), intent(in) :: bkgrdmodel2
integer, intent(in)            :: lfbkgrdmodel2
character(len=3)               :: param !rho, vs,vp

  select case(bkgrdmodel2(1:lfbkgrdmodel2))
  case('ak135')
     velocity=ak135(r0,param,idom)
  case('prem')
     velocity=prem_sub(r0,param,idom)
  case('prem_ani')
     velocity=prem_ani_sub(r0,param,idom)
  case('prem_solid')
     velocity=prem_solid_sub(r0,param,idom)
  case('prem_light')
     velocity=prem_light_sub(r0,param,idom)
  case('prem_light_ani')
     velocity=prem_light_ani_sub(r0,param,idom)
  case('prem_onecrust')
     velocity=prem_onecrust_sub(r0,param,idom)
  case('prem_onecrust_ani')
     velocity=prem_onecrust_ani_sub(r0,param,idom)
  case('prem_solid_light')
     velocity=prem_solid_light_sub(r0,param,idom)
  case('iasp91')
     velocity=iasp91_sub(r0,param,idom)
  case('solar')
     velocity=arbitr_sub_solar(r0,param,idom,bkgrdmodel2)
  case default
     velocity=arbitr_sub(param,idom,bkgrdmodel2)
  end select

end function velocity
!=============================================================================

!-----------------------------------------------------------------------------
logical function model_is_ani(bkgrdmodel2)
!
! returns true if the model is radially anisotrpoic
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

character(len=100), intent(in) :: bkgrdmodel2

  select case(trim(bkgrdmodel2))
  case('prem_ani')
    model_is_ani = .true.
  case('prem_onecrust_ani')
    model_is_ani = .true.
  case('prem_light_ani')
    model_is_ani = .true.
  case default
    model_is_ani = .false.
  end select

end function model_is_ani
!=============================================================================

!-----------------------------------------------------------------------------
double precision function ak135(r0,param,idom)
! from Kennett, Engdahl and Buland, 1995
! interpolated between discontinuities using matlab's polyfit, use radii!!!
! use routine axisem_ak135_fitting.m and ak135/iasp whatever as nd-files
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_ak
double precision             :: ro_ak,vp_ak,vs_ak
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  x_ak=r/6371. ! normalized


  IF(idom==1) THEN
     ro_ak=2.72
     vp_ak=5.8
     vs_ak=3.46
  ELSEIF(idom==2) THEN
     ro_ak=2.92
     vp_ak=6.5
     vs_ak=3.85
  ELSEIF(idom==3) THEN 
  ! moho -> 210
     ro_ak=7.1576-3.859*x_ak
     vp_ak=17.4734-9.5332*x_ak
     vs_ak=5.8556-1.3825*x_ak
  ELSEIF(idom==4) THEN
  ! 210 -> 410
     ro_ak=7.1594-3.8608*x_ak
     vp_ak=30.7877-23.2542*x_ak
     vs_ak=15.2181-11.0601*x_ak
  ELSEIF(idom==5) THEN
  ! 410 -> 660
     ro_ak=11.1204-7.8713*x_ak
     vp_ak=29.389-21.4066*x_ak
     vs_ak=17.7173-13.5065*x_ak
  ELSEIF(idom==6) THEN
  ! 660 -> D''
     ro_ak=6.8294-1.7227*x_ak-1.1064*x_ak*x_ak-0.034409*x_ak*x_ak*x_ak
     vp_ak=26.8598-48.9644*x_ak+63.7326*x_ak*x_ak-32.4155*x_ak*x_ak*x_ak
     vs_ak=18.0019-43.6346*x_ak+60.4205*x_ak*x_ak-29.689*x_ak*x_ak*x_ak
  ELSEIF(idom==7) THEN
  ! D'' -> CMB
     ro_ak=-65.8145+386.221*x_ak-691.6551*x_ak*x_ak+409.6742*x_ak*x_ak*x_ak
     vp_ak=3.4872+55.1872*x_ak-99.0089*x_ak*x_ak+58.7141*x_ak*x_ak*x_ak
     vs_ak=-22.9553+164.0287*x_ak-294.2766*x_ak*x_ak+174.5113*x_ak*x_ak*x_ak
  ELSEIF(idom==8) THEN
  ! CMB -> ICB
     ro_ak=12.592-1.778*x_ak-1.6964*x_ak*x_ak-7.3524*x_ak*x_ak*x_ak
     vp_ak=10.7738-2.4831*x_ak+3.2584*x_ak*x_ak-14.9171*x_ak*x_ak*x_ak
     vs_ak=0
  ELSEIF(idom==9) THEN
  ! inner core
     ro_ak=13.0122-0.0011863*x_ak-8.4449*x_ak*x_ak
     vp_ak=11.2641-0.090247*x_ak-5.7431*x_ak*x_ak
     vs_ak=3.6677+0.0049932*x_ak-4.4808*x_ak*x_ak
  ENDIF

  if (param=='rho') then
     ak135=ro_ak*1000.
  elseif (param=='v_p') then
     ak135=vp_ak*1000.
  elseif (param=='v_s') then
     ak135=vs_ak*1000.
  elseif (param=='vpv') then
     ak135=vp_ak*1000.
  elseif (param=='vsv') then
     ak135=vs_ak*1000.
  elseif (param=='vph') then
     ak135=vp_ak*1000.
  elseif (param=='vsh') then
     ak135=vs_ak*1000.
  elseif (param=='eta') then
     ak135=1.
  else
     write(6,*)'ERROR IN AK135 FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function ak135
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(idom==1)THEN        ! upper crustal layer
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  ELSEIF(idom==2)THEN
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  ELSEIF(idom==3)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(idom==4)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(idom==5)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(idom==6)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(idom==7)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(idom==9)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(idom==10)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  ELSEIF(idom==11)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  if (param=='rho') then
     prem_sub=ro_prem*1000.
  elseif (param=='v_p') then
     prem_sub=vp_prem*1000.
  elseif (param=='v_s') then
     prem_sub=vs_prem*1000.
  elseif (param=='vpv') then
     prem_sub=vp_prem*1000.
  elseif (param=='vsv') then
     prem_sub=vs_prem*1000.
  elseif (param=='vph') then
     prem_sub=vp_prem*1000.
  elseif (param=='vsh') then
     prem_sub=vs_prem*1000.
  elseif (param=='eta') then
     prem_sub=1.
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_ani_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem, vpv_prem, vsv_prem, vph_prem 
double precision             :: vsh_prem, eta_aniso, Qmu, Qkappa
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  eta_aniso=1.

  IF(idom==1)THEN       ! upper crustal layer
     ro_prem=2.6
     vpv_prem=5.8
     vsv_prem=3.2
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==2)THEN   ! lower crustal layer
     ro_prem=2.9
     vpv_prem=6.8
     vsv_prem=3.9
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==3)THEN   ! upper mantle
     ro_prem=2.6910+0.6924*x_prem
     vpv_prem=0.8317+7.2180*x_prem
     vph_prem=3.5908+4.6172*x_prem
     vsv_prem=5.8582-1.4678*x_prem
     vsh_prem=-1.0839+5.7176*x_prem
     eta_aniso=3.3687-2.4778*x_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==4)THEN
     ro_prem=7.1089-3.8045*x_prem
     vpv_prem=20.3926-12.2569*x_prem
     vsv_prem=8.9496-4.4597*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==5)THEN
     ro_prem=11.2494-8.0298*x_prem
     vpv_prem=39.7027-32.6166*x_prem
     vsv_prem=22.3512-18.5856*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==6)THEN
     ro_prem=5.3197-1.4836*x_prem
     vpv_prem=19.0957-9.8672*x_prem
     vsv_prem=9.9839-4.9324*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==7)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vsv_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==9)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==10)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vpv_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vsv_prem=0.0
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=0.0
     Qkappa=57827.0
  ELSEIF(idom==11)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vpv_prem=11.2622-6.3640*x_prem**2
     vsv_prem=3.6678-4.4475*x_prem**2
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=84.6
     Qkappa=1327.7
  ENDIF

  if (param=='rho') then
     prem_ani_sub=ro_prem*1000.
  elseif (param=='vpv') then
     prem_ani_sub=vpv_prem*1000.
  elseif (param=='vsv') then
     prem_ani_sub=vsv_prem*1000.
  elseif (param=='vph') then
     prem_ani_sub=vph_prem*1000.
  elseif (param=='vsh') then
     prem_ani_sub=vsh_prem*1000.
  elseif (param=='eta') then
     prem_ani_sub=eta_aniso
  elseif (param=='Qmu') then
     prem_ani_sub=Qmu
  elseif (param=='Qka') then
     prem_ani_sub=Qkappa
  !min/max velocities needed for the mesher:
  elseif (param=='v_p') then
     prem_ani_sub=max(vpv_prem, vph_prem)*1000.
  elseif (param=='v_s') then
     prem_ani_sub=min(vsv_prem, vsh_prem)*1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:',param,' NOT AN OPTION'
     stop
  endif

end function prem_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_solid_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
! No fluid outer core, but instead vs=vp/sqrt(3)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(idom==1)THEN        ! upper crustal layer
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  ELSEIF(idom==2)THEN
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  ELSEIF(idom==3)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(idom==4)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(idom==5)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(idom==6)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(idom==7)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(idom==9)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(idom==10)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=vp_prem/sqrt(3.)
  ELSEIF(idom==11)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  if (param=='rho') then
     prem_solid_sub=ro_prem*1000.
  elseif (param=='v_p') then
     prem_solid_sub=vp_prem*1000.
  elseif (param=='v_s') then
     prem_solid_sub=vs_prem*1000.
  elseif (param=='vpv') then
     prem_solid_sub=vp_prem*1000.
  elseif (param=='vsv') then
     prem_solid_sub=vs_prem*1000.
  elseif (param=='vph') then
     prem_solid_sub=vp_prem*1000.
  elseif (param=='vsh') then
     prem_solid_sub=vs_prem*1000.
  elseif (param=='eta') then
     prem_solid_sub=1.
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_solid_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_onecrust_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
! but with lower crust extended to the surface
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(idom==1)THEN        ! one crustal layer
     ro_prem=2.6                    
     vp_prem=5.8
     vs_prem=3.2
  ELSEIF(idom==2)THEN      ! upper mantle
     ro_prem=2.691+.6924*x_prem             
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(idom==3)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(idom==4)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(idom==5)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(idom==6)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(idom==7)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(idom==9)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  ELSEIF(idom==10)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  if (param=='rho') then
     prem_onecrust_sub=ro_prem*1000.
  elseif (param=='v_p') then
     prem_onecrust_sub=vp_prem*1000.
  elseif (param=='v_s') then
     prem_onecrust_sub=vs_prem*1000.
  elseif (param=='vpv') then
     prem_onecrust_sub=vp_prem*1000.
  elseif (param=='vsv') then
     prem_onecrust_sub=vs_prem*1000.
  elseif (param=='vph') then
     prem_onecrust_sub=vp_prem*1000.
  elseif (param=='vsh') then
     prem_onecrust_sub=vs_prem*1000.
  elseif (param=='eta') then
     prem_onecrust_sub=1.
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_onecrust_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_onecrust_ani_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem, vpv_prem, vsv_prem, vph_prem 
double precision             :: vsh_prem, eta_aniso, Qmu, Qkappa
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  eta_aniso=1.

  IF(idom==1)THEN       ! upper crustal layer
     ro_prem=2.6
     vpv_prem=5.8
     vsv_prem=3.2
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==2)THEN   ! upper mantle
     ro_prem=2.6910+0.6924*x_prem
     vpv_prem=0.8317+7.2180*x_prem
     vph_prem=3.5908+4.6172*x_prem
     vsv_prem=5.8582-1.4678*x_prem
     vsh_prem=-1.0839+5.7176*x_prem
     eta_aniso=3.3687-2.4778*x_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==3)THEN
     ro_prem=7.1089-3.8045*x_prem
     vpv_prem=20.3926-12.2569*x_prem
     vsv_prem=8.9496-4.4597*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==4)THEN
     ro_prem=11.2494-8.0298*x_prem
     vpv_prem=39.7027-32.6166*x_prem
     vsv_prem=22.3512-18.5856*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==5)THEN
     ro_prem=5.3197-1.4836*x_prem
     vpv_prem=19.0957-9.8672*x_prem
     vsv_prem=9.9839-4.9324*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==6)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==7)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vsv_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==8)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==9)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vpv_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vsv_prem=0.0
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=0.0
     Qkappa=57827.0
  ELSEIF(idom==10)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vpv_prem=11.2622-6.3640*x_prem**2
     vsv_prem=3.6678-4.4475*x_prem**2
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=84.6
     Qkappa=1327.7
  ENDIF

  if (param=='rho') then
     prem_onecrust_ani_sub=ro_prem*1000.
  elseif (param=='vpv') then
     prem_onecrust_ani_sub=vpv_prem*1000.
  elseif (param=='vsv') then
     prem_onecrust_ani_sub=vsv_prem*1000.
  elseif (param=='vph') then
     prem_onecrust_ani_sub=vph_prem*1000.
  elseif (param=='vsh') then
     prem_onecrust_ani_sub=vsh_prem*1000.
  elseif (param=='eta') then
     prem_onecrust_ani_sub=eta_aniso
  elseif (param=='Qmu') then
     prem_onecrust_ani_sub=Qmu
  elseif (param=='Qka') then
     prem_onecrust_ani_sub=Qkappa
  !min/max velocities needed for the mesher:
  elseif (param=='v_p') then
     prem_onecrust_ani_sub=max(vpv_prem, vph_prem)*1000.
  elseif (param=='v_s') then
     prem_onecrust_ani_sub=min(vsv_prem, vsh_prem)*1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:',param,' NOT AN OPTION'
     stop
  endif

end function prem_onecrust_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_light_sub(r0,param,idom)
!
! prem_light model (crust removed) in terms of domains separated by disconts.
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(idom==1)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(idom==2)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(idom==3)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(idom==4)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(idom==5)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(idom==6)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(idom==7)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(idom==8)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.0
  ELSEIF(idom==9)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  if (param=='rho') then
     prem_light_sub=ro_prem*1000.
  elseif (param=='v_p') then
     prem_light_sub=vp_prem*1000.
  elseif (param=='v_s') then
     prem_light_sub=vs_prem*1000.
  elseif (param=='vpv') then
     prem_light_sub=vp_prem*1000.
  elseif (param=='vsv') then
     prem_light_sub=vs_prem*1000.
  elseif (param=='vph') then
     prem_light_sub=vp_prem*1000.
  elseif (param=='vsh') then
     prem_light_sub=vs_prem*1000.
  elseif (param=='eta') then
     prem_light_sub=1.
  else
     write(6,*)'ERROR IN PREM_LIGHT_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_light_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_light_ani_sub(r0,param,idom)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem, vpv_prem, vsv_prem, vph_prem 
double precision             :: vsh_prem, eta_aniso, Qmu, Qkappa
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  eta_aniso=1.

  IF(idom==1)THEN
     ro_prem=2.6910+0.6924*x_prem
     vpv_prem=0.8317+7.2180*x_prem
     vph_prem=3.5908+4.6172*x_prem
     vsv_prem=5.8582-1.4678*x_prem
     vsh_prem=-1.0839+5.7176*x_prem
     eta_aniso=3.3687-2.4778*x_prem
     Qmu=600.0
     Qkappa=57827.0
  ELSEIF(idom==2)THEN
     ro_prem=7.1089-3.8045*x_prem
     vpv_prem=20.3926-12.2569*x_prem
     vsv_prem=8.9496-4.4597*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==3)THEN
     ro_prem=11.2494-8.0298*x_prem
     vpv_prem=39.7027-32.6166*x_prem
     vsv_prem=22.3512-18.5856*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==4)THEN
     ro_prem=5.3197-1.4836*x_prem
     vpv_prem=19.0957-9.8672*x_prem
     vsv_prem=9.9839-4.9324*x_prem
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=143.0
     Qkappa=57827.0
  ELSEIF(idom==5)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==6)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vsv_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==7)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vpv_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vsv_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=312.0
     Qkappa=57827.0
  ELSEIF(idom==8)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vpv_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vsv_prem=0.0
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=0.0
     Qkappa=57827.0
  ELSEIF(idom==9)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vpv_prem=11.2622-6.3640*x_prem**2
     vsv_prem=3.6678-4.4475*x_prem**2
     vph_prem=vpv_prem
     vsh_prem=vsv_prem
     Qmu=84.6
     Qkappa=1327.7
  ENDIF

  if (param=='rho') then
     prem_light_ani_sub=ro_prem*1000.
  elseif (param=='vpv') then
     prem_light_ani_sub=vpv_prem*1000.
  elseif (param=='vsv') then
     prem_light_ani_sub=vsv_prem*1000.
  elseif (param=='vph') then
     prem_light_ani_sub=vph_prem*1000.
  elseif (param=='vsh') then
     prem_light_ani_sub=vsh_prem*1000.
  elseif (param=='eta') then
     prem_light_ani_sub=eta_aniso
  elseif (param=='Qmu') then
     prem_light_ani_sub=Qmu
  elseif (param=='Qka') then
     prem_light_ani_sub=Qkappa
  !min/max velocities needed for the mesher:
  elseif (param=='v_p') then
     prem_light_ani_sub=max(vpv_prem, vph_prem)*1000.
  elseif (param=='v_s') then
     prem_light_ani_sub=min(vsv_prem, vsh_prem)*1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:',param,' NOT AN OPTION'
     stop
  endif

end function prem_light_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
double precision function prem_solid_light_sub(r0,param,idom)
!
! prem_light model (crust removed) in terms of domains separated by disconts.
! No fluid outer core, but instead vs=vp/sqrt(3)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x_prem
double precision             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(idom==1)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(idom==2)THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(idom==3)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(idom==4)THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(idom==5)THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(idom==6)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(idom==7)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(idom==8)THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=vp_prem/sqrt(3.)
  ELSEIF(idom==9)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ENDIF

  if (param=='rho') then
     prem_solid_light_sub=ro_prem*1000.
  elseif (param=='v_p') then
     prem_solid_light_sub=vp_prem*1000.
  elseif (param=='v_s') then
     prem_solid_light_sub=vs_prem*1000.
  elseif (param=='vpv') then
     prem_solid_light_sub=vp_prem*1000.
  elseif (param=='vsv') then
     prem_solid_light_sub=vs_prem*1000.
  elseif (param=='vph') then
     prem_solid_light_sub=vp_prem*1000.
  elseif (param=='vsh') then
     prem_solid_light_sub=vs_prem*1000.
  elseif (param=='eta') then
     prem_solid_light_sub=1.
  else
     write(6,*)'ERROR IN PREM_LIGHT_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_solid_light_sub
!=============================================================================



!-----------------------------------------------------------------------------
double precision function iasp91_sub(r0,param,idom)
!
! iasp91 model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double precision, intent(in) :: r0
integer, intent(in)          :: idom
double precision             :: r,x
double precision             :: rho,vp,vs
character(len=3), intent(in) :: param !rho, vs,vp

double precision :: RICB,RCMB,RTOPDDOUBLEPRIME,R771,R670,R400,R220
double precision :: R120,RMOHO,RMIDDLE_CRUST,ROCEAN

double precision x1,x2

! compute real physical radius in meters

  r=r0
  x=r/6371000.     ! Radius (normalized to x(surface)=1 )


! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0
    R120 = 6251000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
! there is no d600 discontinuity in IASP91 therefore this value is useless
! but it needs to be there for compatibility with other subroutines
!    R600 = R_EARTH - 600000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    x1 = R120 / ROCEAN
    x2 = RMOHO / ROCEAN

  IF(idom==1)THEN ! upper crust 
     vp = 5.8d0
     vs = 3.36d0
     rho = 2.72d0
  ELSEIF(idom==2)THEN ! lower crust
     vp = 6.5d0
     vs = 3.75d0
     rho = 2.92d0
  ELSEIF(idom==3)THEN ! R120 < r <= RMOHO
     vp = 8.78541d0-0.74953d0*x
     vs = 6.706231d0-2.248585d0*x
     rho = 3.3713d0 + (3.3198d0-3.3713d0)*(x-x1)/(x2-x1)
     if(rho < 3.30d0 .or. rho > 3.38d0) then 
        write(6,*) 'incorrect density computed for IASP91',rho
        stop
     endif
  ELSEIF(idom==4)THEN ! R220 < r <= R120
     rho=2.6910d0+0.6924d0*x
     vp=25.41389-17.69722*x
     vs=5.75020-1.27420*x
  ELSEIF(idom==5)THEN ! R400 < r <= R220
     rho=7.1089d0-3.8045d0*x
     vp=30.78765-23.25415*x
     vs=15.24213-11.08552*x
  ELSEIF(idom==6)THEN ! R670 < r <= R400
     rho=5.3197d0-1.4836d0*x
     vp=29.38896-21.40656*x
     vs=17.70732-13.50652*x
  ELSEIF(idom==7)THEN ! R771 < r <= R670
     rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
     vp=25.96984-16.93412*x
     vs=20.76890-16.53147*x
  ELSEIF(idom==8)THEN ! RTOPDDOUBLEPRIME < r <= R771
     rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
     vp=25.1486-41.1538*x+51.9932*x**2-26.6083*x**3
     vs=12.9303-21.2590*x+27.8988*x**2-14.1080*x**3
  ELSEIF(idom==9)THEN ! RCMB < r <= RTOPDDOUBLEPRIME
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=14.49470-1.47089*x
    vs=8.16616-1.58206*x
  ELSEIF(idom==10)THEN ! RICB < r <= RCMB
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=10.03904+3.75665*x-13.67046*x**2
    vs=0.0d0
  ELSEIF(idom==11)THEN ! 0.d0 < r <= RICB
    rho=13.0885d0-8.8381d0*x*x
    vp=11.24094-4.09689*x**2
    vs=3.56454-3.45241*x**2
  ELSE
     write(6,*) 'iasp91_sub: error with domain idom=',idom
     stop
  ENDIF

  if (param=='rho') then
     iasp91_sub=rho*1000.
  elseif (param=='v_p') then
     iasp91_sub=vp*1000.
  elseif (param=='v_s') then
     iasp91_sub=vs*1000.
  elseif (param=='vpv') then
     iasp91_sub=vp*1000.
  elseif (param=='vsv') then
     iasp91_sub=vs*1000.
  elseif (param=='vph') then
     iasp91_sub=vp*1000.
  elseif (param=='vsh') then
     iasp91_sub=vs*1000.
  elseif (param=='eta') then
     iasp91_sub=1.
  else
     write(6,*)'ERROR IN IASP91_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function iasp91_sub
!=============================================================================



!-----------------------------------------------------------------------------
double precision function arbitr_sub(param,idom,bkgrdmodel2)
!
! file-based, step-wise model in terms of domains separated by disconts.
! format:
! ndisc
! r vp vs rho
! ...
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

integer, intent(in) :: idom
character(len=100), intent(in) :: bkgrdmodel2
character(len=3), intent(in) :: param !rho, vs,vp
double precision, allocatable, dimension(:) :: disconttmp,rhotmp,vstmp,vptmp
integer :: ndisctmp,i,ndisctmp2
logical :: bkgrdmodelfile_exists

! Does the file bkgrdmodel".bm" exist?
    inquire(file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm', &
                exist=bkgrdmodelfile_exists)
  
    if (bkgrdmodelfile_exists) then
        open(unit=77,file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm')
        read(77,*)ndisctmp
  
        ndisctmp2=ndisctmp
        if (ndisctmp==1) ndisctmp2 = 2
        allocate(disconttmp(1:ndisctmp2))
        allocate(vptmp(1:ndisctmp2),vstmp(1:ndisctmp2),rhotmp(1:ndisctmp2))
        do i=1, ndisctmp
            read(77,*)disconttmp(i),rhotmp(i),vptmp(i),vstmp(i)
        enddo
        close(77)
        if (ndisctmp==1) then
            disconttmp(2)=disconttmp(1)/4.
            vptmp(2) = vptmp(1)
            vstmp(2) = vstmp(1)
            rhotmp(2)=rhotmp(1)
        end if
  
        if (param=='rho') arbitr_sub = rhotmp(idom)
        if (param=='v_p') arbitr_sub = vptmp(idom)
        if (param=='v_s') arbitr_sub = vstmp(idom)
        if (param=='vpv') arbitr_sub = vptmp(idom)
        if (param=='vsv') arbitr_sub = vstmp(idom)
        if (param=='vph') arbitr_sub = vptmp(idom)
        if (param=='vsh') arbitr_sub = vstmp(idom)
        if (param=='eta') arbitr_sub = 1.
        deallocate(disconttmp,vstmp,vptmp,rhotmp)
    else 
        write(6,*)'Background model file', &
            bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm','does not exist!!!'
        stop
    endif

end function arbitr_sub
!=============================================================================



!-----------------------------------------------------------------------------
double precision function arbitr_sub_solar(r0,param,idom,bkgrdmodel2)
!
! file-based, step-wise model in terms of domains separated by disconts.
! format:
! ndisc
! r vp vs rho
! ...
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double precision, intent(in) :: r0
integer, intent(in) :: idom
character(len=100), intent(in) :: bkgrdmodel2
character(len=3), intent(in) :: param !rho, vs,vp
double precision, allocatable, dimension(:) :: disconttmp,rhotmp,vstmp,vptmp
integer :: ndisctmp,i,ndisctmp2,ind(2)
logical :: bkgrdmodelfile_exists
double precision :: w(2),wsum

! Does the file bkgrdmodel".bm" exist?
  inquire(file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm', &
          exist=bkgrdmodelfile_exists)
  if (bkgrdmodelfile_exists) then
     open(unit=77,file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm')
     read(77,*)ndisctmp
!     write(6,*)'num discont:',ndisctmp
!     write(6,*)'radius:',r0
     allocate(disconttmp(1:ndisctmp))
     allocate(vptmp(1:ndisctmp),vstmp(1:ndisctmp),rhotmp(1:ndisctmp))
     do i=1, ndisctmp
        read(77,*)disconttmp(i),rhotmp(i),vptmp(i),vstmp(i)
     enddo
     close(77)
     
     call interp_vel(r0,disconttmp(1:ndisctmp),ndisctmp,ind,w,wsum)

     if (param=='rho') arbitr_sub_solar=sum(w*rhotmp(ind))*wsum
     if (param=='v_p') &!arbitr_sub_solar=sum(w*vptmp(ind))*wsum
          arbitr_sub_solar=(w(1)*vptmp(ind(1))+w(2)*vptmp(ind(2)))*wsum
     if (param=='v_s') arbitr_sub_solar=sum(w*vstmp(ind))*wsum
     deallocate(disconttmp,vstmp,vptmp,rhotmp)
  else 
     write(6,*)'Background model file', &
          bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm','does not exist!!!'
     stop
  endif

end function arbitr_sub_solar
!=============================================================================


!=============================================================================
subroutine interp_vel(r0,r,n,ind,w,wsum)

integer, intent(in) :: n
double precision, intent(in) :: r0,r(1:n)
integer, intent(out) :: ind(2)
double precision, intent(out) :: w(2),wsum
integer :: i,p
double precision :: dr1,dr2

p=1

     i=minloc(abs(r-r0),1)
!     write(6,*)'INTERP;',r0,i
!     write(6,*)'INTERP:',r(i)

     if (r0>0.d0) then
        if ((r(i)-r0)/r0> 1.d-8) then ! closest discont. at larger radius
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))
        elseif ((r0-r(i))/r0> 1.d-8) then  ! closest discont. at smaller radius
           if (r0>maxval(r)) then ! for round-off errors where mesh is above surface
              ind(1)=i
              ind(2)=i
              dr1=1.d0
              dr2=1.d0
           else
              ind(1)=i-1
              ind(2)=i
              dr1=r(ind(1))-r0
              dr2=r0-r(ind(2))
            endif
        elseif (abs((r(i)-r0)/r0)< 1.d-8) then ! closest discont identical
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           write(6,*)'problem with round-off errors in interpolating......'
           write(6,*)'r0,r(i),i',r0,r(i),abs((r(i)-r0)/r0),i
           stop
        endif
     else !r0=0
        if (r(i)==0.d0) then ! center of the sun
           ind(1)=i
           ind(2)=i
           dr1=1.d0
           dr2=1.d0
        else
           ind(1)=i
           ind(2)=i+1
           dr1=r(ind(1))-r0
           dr2=r0-r(ind(2))        
        endif
     endif

! inverse distance weighting
     w(1)=(dr1)**(-p)
     w(2)=(dr2)**(-p)
     wsum=1./sum(w)

end subroutine interp_vel
!=============================================================================



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!========================
end module background_models
!========================
