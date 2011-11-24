!========================
module analyt_homog_radiation
!========================

  use global_parameters
  use data_mesh
  use data_time
  use data_source, ONLY: t_0,zsrc,src_type,stf_type,magnitude,iel_src,jpol_src
  use data_io

  implicit none

  public :: prepare_src_vicinity,src_vicinity
  private

  integer                       :: nrdist
  double precision, allocatable :: rdist(:),s_srcvic(:),z_srcvic(:)
  integer, allocatable          :: iel_srcvic(:),ipol_srcvic(:),jpol_srcvic(:)
  logical                       :: surface_reflection

  contains

!//////////////////////////////////////////////////

!--------------------------------------------------------------
!dk prepare_src_vicinity
!--------------------------------------------------------------
subroutine prepare_src_vicinity

  use data_proc, ONLY : appmynum

  include 'mesh_params.h'

  integer :: iel
  character(len=4)  :: appiter

  surface_reflection=.true.

!     call srcvic_midmantle_recs
     call srcvic_corenuke_recs
     write(6,*)'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(6,*)'oooooooooooo Analytical P-wave radiation ooooooooooooooooo'
     write(6,*)'  Vp [km/s]         :',vphomo/1.d3
     write(6,*)'  Vs [km/s]         :',vshomo/1.d3
     write(6,*)'  Density [g/cm^3]  :',rhohomo/1.d3
     write(6,*)'  Surface reflection?',surface_reflection 
     write(6,*)
     write(6,*)'Number of receivers in source vicinity :',nrdist

     allocate(rdist(1:nrdist))
     allocate(s_srcvic(1:nrdist))
     allocate(z_srcvic(1:nrdist))
     allocate(iel_srcvic(1:nrdist))
     allocate(ipol_srcvic(1:nrdist))
     allocate(jpol_srcvic(1:nrdist))

     write(6,*)'Number of axial solid elements:',naxel_solid
     write(6,*)'P-wavelength [km]:',10.d3*t_0/1000.
     write(6,*)'Closest distance (r [km], wavelengths) :', &
                minval(rdist(1:nrdist))/1000.,minval(rdist(1:nrdist))/10.d3/t_0
     write(6,*)'Furthest distance (r [km], wavelengths):', &
                maxval(rdist(1:nrdist))/1000.,maxval(rdist(1:nrdist))/10.d3/t_0
     write(6,*)'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(6,*)'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(6,*)

     open(unit=1111,file=infopath(1:lfinfo)//'src_vicinity_radii.dat'//appmynum)
     do iel = 1,nrdist
        write(1111,12)iel_srcvic(iel),rdist(iel), & 
                      s_srcvic(iel)*rdist(iel),z_srcvic(iel)*rdist(iel)+zsrc
        call define_io_appendix(appiter,iel)
        open(unit=11111+iel,file=datapath(1:lfdata)//&
                                 '/srcvic_seis.dat'//appiter)
        open(unit=12111+iel,file=datapath(1:lfdata)//&
                                 '/srcvic_seisnearfar.dat'//appiter)
     enddo
     close(1111)
12 format(i10,3(1pe15.5))


end subroutine prepare_src_vicinity
!--------------------------------------------------------------------------

!dk src_vicinity
!--------------------------------------------------------------------------
subroutine src_vicinity(iter,disp)

  use utlity, ONLY: eps2zero

  integer, intent(in)                      :: iter
  real(kind=realkind), intent(in)          :: disp(0:npol,0:npol,1:nel_solid,3)
  real(kind=realkind),dimension(nel_solid) :: ur_anal_near,ur_anal_far
  integer                                  :: iel

  if ( mod(iter,5)==0 ) then
     
     if (surface_reflection) then 
        call corenukesurf_expl_analyt(t,ur_anal_near,ur_anal_far)
     else
        call prerefl_expl_analyt(t,ur_anal_near,ur_anal_far)  
     endif

     do iel = 1,nrdist  
        
        write(11111+iel,*)t,ur_anal_far(iel)+ur_anal_near(iel),s_srcvic(iel)* &
                 disp(ipol_srcvic(iel),jpol_srcvic(iel),iel_srcvic(iel),1) + &
                 z_srcvic(iel) * &
                 disp(ipol_srcvic(iel),jpol_srcvic(iel),iel_srcvic(iel),3)

        write(12111+iel,12)t,eps2zero(ur_anal_near(iel)), &
                             eps2zero(ur_anal_far(iel)), &
                             eps2zero(ur_anal_near(iel)) + &
                             eps2zero(ur_anal_far(iel))
     enddo
  endif ! mod(iter,10)

12 format(1pe15.4,3(1pe17.9))

end subroutine src_vicinity
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
subroutine srcvic_midmantle_recs

  use data_mesh_preloop
  use utlity

  double precision  :: s,z,r,theta
  integer           :: iel

  write(6,*)
  write(6,*)'SOURCE VICINITY: DEFINING MIDMANTLE RECEIVERS...'

     nrdist = 0
     do iel = 1, naxel_solid
        if ( north(ielsolid(ax_el_solid(iel))) ) then
           call compute_coordinates(s,z,r,theta,ielsolid(ax_el_solid(iel)), &
                                    0,npol/2)
           nrdist = nrdist + 1
           rdist(nrdist) = abs(zsrc-z)
           iel_srcvic(nrdist) = ax_el_solid(iel)
           ipol_srcvic(nrdist) = 0
           jpol_srcvic(nrdist) = npol/2
           s_srcvic(nrdist) = s/rdist(nrdist)
           z_srcvic(nrdist) = (z-zsrc)/rdist(nrdist)
        endif ! north
     enddo

     iel = 1
     theta=0
     do while ( theta < 85./180.*pi ) 
        iel = iel + 1
        call compute_coordinates(s,z,r,theta,ielsolid(iel_src+iel), &
                                 npol/2,jpol_src)
        nrdist = nrdist + 1
        rdist(nrdist) = sqrt( (zsrc-z)**2  + s**2 )
        iel_srcvic(nrdist) = iel_src+iel
        ipol_srcvic(nrdist) = npol/2
        jpol_srcvic(nrdist) = jpol_src
        s_srcvic(nrdist) = s/rdist(nrdist)
        z_srcvic(nrdist) = (z-zsrc)/rdist(nrdist)
     enddo

end subroutine srcvic_midmantle_recs
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine srcvic_corenuke_recs

  use data_mesh_preloop
  use utlity

  double precision  :: s,z,r,theta
  integer           :: iel

  write(6,*)
  write(6,*)'SOURCE VICINITY: DEFINING CORE-NUKE RECEIVERS....'

  nrdist = 0
   
! axial elements
   do iel = 1, naxel_solid

   ! north/south pole
      call compute_coordinates(s,z,r,theta,ielsolid(ax_el_solid(iel)),0,npol)

      if (dblreldiff_small(r,router)) then
         write(6,*)'...found a pole! r,z [km]:',r/1.d3,z/1.d3
         nrdist = nrdist + 1
         rdist(nrdist) = abs(z)
         iel_srcvic(nrdist) = ax_el_solid(iel)
         ipol_srcvic(nrdist) = 0
         jpol_srcvic(nrdist) = npol/2
         s_srcvic(nrdist) = s/rdist(nrdist)
         z_srcvic(nrdist) = z/rdist(nrdist)         
      endif

      call compute_coordinates(s,z,r,theta,ielsolid(ax_el_solid(iel)),0,npol/2)

!    mid-element
      if (r /=zero) then
         nrdist = nrdist + 1
         rdist(nrdist) = abs(z)
         iel_srcvic(nrdist) = ax_el_solid(iel)
         ipol_srcvic(nrdist) = 0
         jpol_srcvic(nrdist) = npol/2
         s_srcvic(nrdist) = s/rdist(nrdist)
         z_srcvic(nrdist) = z/rdist(nrdist)
      endif
   enddo

! equatorial elements
     do iel = 1, nel_solid
        if ( north(ielsolid(iel)) ) then
           call compute_coordinates(s,z,r,theta,ielsolid(iel),npol,npol/2)
           if (z<min_distance_dim .and. r/=zero) then
              nrdist = nrdist + 1
              rdist(nrdist) = abs(s)
              iel_srcvic(nrdist) = iel
              ipol_srcvic(nrdist) = npol
              jpol_srcvic(nrdist) = npol/2
              s_srcvic(nrdist) = s/rdist(nrdist)
              z_srcvic(nrdist) = z/rdist(nrdist)
           endif
        endif
     enddo

end subroutine srcvic_corenuke_recs
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine prerefl_expl_analyt(t_in,ur_anal_near,ur_anal_far)

  include 'mesh_params.h'

  double precision,intent(in)   :: t_in
  real(kind=realkind), intent(out) :: ur_anal_far(nrdist),ur_anal_near(nrdist)
  integer          :: ir
  double precision :: tt,ttp,tts,decay,shift_fact,stf,deriv_stf

  decay=3.5d0
  shift_fact=1.5d0

  if (src_type(2) /= 'explosion') then 
    write(6,*)''
    write(6,*)'Problem in analytical explosion radiation!'
    write(6,*)'Analytical expressions not implemented for ',src_type(2)
    stop
  endif

  do ir=1,nrdist

     if (rdist(ir)==zero) then
        write(6,*)'Analytical expression not defined for r=0',ir,rdist(ir)
        ur_anal_near(ir)=half*stf
        ur_anal_far(ir)=half*stf
     else

      ttp=t_in-rdist(ir)/vphomo
      tts=t_in-rdist(ir)/vshomo
      tt=t_in-rdist(ir)/vphomo

! Kennett, Vol 1, p.204, eq. (11.4.23)

      if (stf_type=='gauss_0') then
        stf = exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        deriv_stf= -two*(decay/t_0)**2*(tt-shift_fact*t_0) * &
                     exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )

      elseif (stf_type=='gauss_1') then
        stf = -(decay/t_0)**2*(tt-shift_fact*t_0) * &
                exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        deriv_stf=(decay/t_0)**2 * ( two*(decay/t_0)**2 * &
                  (tt-shift_fact*t_0)**2 - 1 ) *&
                   exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        stf=stf/( decay/t_0*sqrt(half)*exp(-half) )
        deriv_stf=deriv_stf/( decay/t_0*sqrt(half)*exp(-half) )

!!$     elseif (stf_type=='gauss_2') then
!!$
!!$        integrp=(decay/t_0*(ttp-shift_fact*t_0))**2
!!$        integrs=(decay/t_0*(tts-shift_fact*t_0))**2
!!$
!!$        prefactp=(decay/t_0)**2 * ( (decay/t_0)**2 *(ttp-shift_fact*t_0)**2 -1)
!!$        prefacts=(decay/t_0)**2 * ( (decay/t_0)**2 *(tts-shift_fact*t_0)**2 -1)
!!$        deriv_prefactp=two*(decay/t_0)**4*(ttp-shift_fact*t_0)*( three-integrp)
!!$
!!$        scalefac= two*((decay/t_0)**2)*exp(-three/two)
!!$
!!$        ur_anal(ir)=magnitude/(two*pi*rho*scalefac) * &
!!$       (three/rdist(ir)**4*((half*rdist(ir)/vs*integrs-0.25d0)*exp(-integrs)- &
!!$                           (half*rdist(ir)/vp*integrp-0.25d0)*exp(-integrp))+ &
!!$           three/(two*vp**2*rdist(ir)**2)*exp(-integrp)- &
!!$           exp(-integrs)/(vs**2*rdist(ir)**2)+ &
!!$           deriv_prefactp/(two*vp**2*rdist(ir))*exp(-integrp))

     else
        write(6,*)''
        write(6,*)'Problem in analytical explosion radiation!'
        write(6,*)'Analytical expressions not implemented for ',stf_type
        stop
      endif

!      ur_anal(ir)= magnitude*(stf/rdist(ir) + deriv_stf/vp)/&
!                   (four*pi*rho*rdist(ir)*vp**2)

      ur_anal_near(ir)= stf/rdist(ir)**2
      ur_anal_far(ir)= deriv_stf/(vphomo*rdist(ir))

      endif ! r=0 or not

  enddo

  ur_anal_far = ur_anal_far*magnitude/(four*pi*rhohomo*vphomo**2)
  ur_anal_near = ur_anal_near*magnitude/(four*pi*rhohomo*vphomo**2)

end subroutine prerefl_expl_analyt
!---------------------------------------------------

!--------------------------------------------------------------------------
subroutine corenukesurf_expl_analyt(t_in,ur_anal_near,ur_anal_far)

  include 'mesh_params.h'

  double precision,intent(in)   :: t_in
  real(kind=realkind), intent(out) :: ur_anal_far(nrdist),ur_anal_near(nrdist)
  integer                       :: ir,i
  double precision :: tt,ttp,decay,shift_fact,stf,deriv_stf

  decay=3.5d0
  shift_fact=1.5d0
  ur_anal_near = zero; ur_anal_far = zero

  if (src_type(2) /= 'explosion') then 
    write(6,*)''
    write(6,*)'Problem in analytical explosion radiation!'
    write(6,*)'Analytical expressions not implemented for ',src_type(2)
    stop
  endif

  do ir=1,nrdist ! receiver locations
     do i = 1,5  ! upgoing, downgoing, upgoing waves

     if (rdist(ir)==zero) then
        write(6,*)'Analytical expression not defined for r=0',ir,rdist(ir)
        ur_anal_near(ir)=half*stf
        ur_anal_far(ir)=half*stf
     else      
        if (i==1) then ! upgoing, direct
           ttp=t_in-rdist(ir)/vphomo
           tt =t_in-rdist(ir)/vphomo
        elseif (i==2) then ! downgoing, reflected
           ttp=t_in-rdist(ir)/vphomo
           tt =t_in+(rdist(ir)-two*router)/vphomo
        elseif (i==3) then ! upgoing 2nd time
           ttp=t_in-(rdist(ir)+two*router)/vphomo
           tt =t_in-(rdist(ir)+two*router)/vphomo
        elseif (i==4) then ! downgoing 2nd time
           ttp=t_in-rdist(ir)/vphomo
           tt =t_in+(rdist(ir)-four*router)/vphomo
        elseif (i==5) then ! upgoing 3rd time
           ttp=t_in-(rdist(ir)+four*router)/vphomo
           tt =t_in-(rdist(ir)+four*router)/vphomo
        endif

! Kennett, Vol 1, p.204, eq. (11.4.23)
      if (stf_type=='gauss_0') then
        stf = exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        deriv_stf= -two*(decay/t_0)**2*(tt-shift_fact*t_0) * &
                     exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
!        stf = exp(-( (decay/t_0*(ttp-shift_fact*t_0))**2) )
!        deriv_stf= -two*(decay/t_0)**2*(tt-shift_fact*t_0) * &
!                     exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )

      elseif (stf_type=='gauss_1') then
        stf = -(decay/t_0)**2*(tt-shift_fact*t_0) * &
                exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        deriv_stf=(decay/t_0)**2 * ( two*(decay/t_0)**2 * &
                  (tt-shift_fact*t_0)**2 - 1 ) *&
                   exp(-( (decay/t_0*(tt-shift_fact*t_0))**2) )
        stf=stf/( decay/t_0*sqrt(half)*exp(-half) )
        deriv_stf=deriv_stf/( decay/t_0*sqrt(half)*exp(-half) )

     else
        write(6,*)''
        write(6,*)'Problem in analytical explosion radiation!'
        write(6,*)'Analytical expressions not implemented for ',stf_type
        stop
      endif

      if (mod(i,2)==0) then 
         ur_anal_near(ir) = ur_anal_near(ir) + stf/rdist(ir)**2 
!         ur_anal_near(ir) = ur_anal_near(ir) + 0.5*stf/rdist(ir)**2 - &
!                            one/(8.d3*pi**vp*rdist(ir))*(tt-shift_fact*t_0)

      else
        ur_anal_near(ir) = ur_anal_near(ir) + stf/rdist(ir)**2
      endif
      ur_anal_far(ir)  = ur_anal_far(ir) + deriv_stf/(vphomo*rdist(ir))

      endif ! r=0 or not
   enddo !up,down waves
enddo

  ur_anal_far = ur_anal_far*magnitude/(four*pi*rhohomo*vphomo**2)
  ur_anal_near = ur_anal_near*magnitude/(four*pi*rhohomo*vphomo**2)

end subroutine corenukesurf_expl_analyt
!---------------------------------------------------


!========================
end module analyt_homog_radiation
!========================
