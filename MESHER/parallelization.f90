!=========================
  module parallelization
!=========================
  use data_grid
  use data_gllmesh
  use data_mesh
  use data_spec
  use data_pdb
  use data_diag
  use data_bkgrdmodel
  implicit none
  public :: create_domain_decomposition

  logical, dimension(:), allocatable :: attributed

  private 
  contains

!dk create_domain_decomposition_theta--------------------------------------
  subroutine create_domain_decomposition
!
! nel : number of glocal elements, i.e. total number of a processor's elements
! nelmax: maximal number of glocal elements
! neltot: global total number of elements
! procel: 

  integer :: iproc,iel,nelmax,nelmax_fluid,nelmax_solid

  write(6,*)'     creating domain decomposition....'

! DD parameters
  nproc = nproc_target
  call check_nproc(nproc)

  allocate(attributed(1:neltot))
  attributed(:) = .false.
!
  allocate(nel(0:nproc-1))
  allocate(nel_fluid(0:nproc-1))
  allocate(nel_solid(0:nproc-1))
  if (mod(neltot_fluid,nproc)/=0) then
   write(6,*) ' neltot_fluid is not a multiple of nproc '
   stop
  end if
  if (mod(neltot_solid,nproc)/=0) then
   write(6,*)'NELTOT_SOLID=',neltot_solid,'NPROC=',nproc
   write(6,*) ' neltot_solid is not a multiple of nproc '
   stop
  end if
!
  if ( mod(neltot,nproc) == 0) then
   do iproc = 0, nproc-1
    nel(iproc) = neltot / nproc
    nel_fluid(iproc) = neltot_fluid / nproc
    nel_solid(iproc) = neltot_solid / nproc
   end do
  else
   do iproc = 0, nproc -1
    if (iproc < nproc -1 )  nel(iproc) = neltot / (nproc-1)
    if (iproc == (nproc-1)) nel(iproc) = neltot - SUM (nel(0:nproc-2))
   end do
  end if
  nelmax = maxval(nel)
  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)
!==================================  
  if (dump_mesh_info_screen) then 
!==================================  

  write(6,*)'************** VARIOUS NUMBERS OF ELEMENTS**********************'
  do iproc=0,nproc-1
     write(6,13)iproc,'has tot/sol/flu number      =',&
                   nel(iproc),nel_solid(iproc),nel_fluid(iproc)
  enddo
  write(6,12)'Maximal glocal total number nelmax      =',&
              nelmax,nelmax_solid,nelmax_fluid
  write(6,*)
  write(6,14)'Sum over total elements for all procs   =',sum(nel)
  write(6,14)'Global, total number neltot             =',neltot
  write(6,14)'Sum over solid elements for all procs   =',sum(nel_solid)
  write(6,14)'Global, solid-domain number neltot_solid=',neltot_solid
  write(6,14)'Sum over fluid elements for all procs   =',sum(nel_fluid)
  write(6,14)'Global, fluid-domain number neltot_fluid=',neltot_fluid

13 format(i4,a37,3(i8))
12 format(a41,3(i8))
14 format(a41,i8)
  write(6,*)'*****************************************************************'

!==================================  
  end if
!==================================  

  allocate(procel(nelmax,0:nproc-1))
  allocate(procel_fluid(nelmax_fluid,0:nproc-1))
  allocate(procel_solid(nelmax_solid,0:nproc-1))
  allocate(el2proc(neltot)); el2proc(:)=-1

! TNM: NOT ENTIRELY SURE IF THIS IS CORRECT!!!!! seems so far... but be aware..
!  allocate(inv_procel(nelmax,0:nproc-1)); inv_procel(:,:)=-1
  allocate(inv_procel(neltot,0:nproc-1)); inv_procel(:,:)=-1

! Decompose such that each processor owns a cake piece in the theta direction,
! i.e. same amount of elements in solid and fluid respectively. 
! The inner cube is done such that each processor maximally has 2 neighbors.
  call domain_decomposition_theta

! write out procel arrays
if (dump_mesh_info_files) then
open(unit=666,file=diagpath(1:lfdiag)//'/inv_procel.dat')
do iproc=0,nproc-1
do iel=1,neltot
   if (el2proc(iel)==iproc) then 
      write(666,*)iproc,iel,inv_procel(iel,iproc),&
                  procel(inv_procel(iel,iproc),iproc)
   endif
enddo
enddo
close(666)

open(unit=666,file=diagpath(1:lfdiag)//'/procel.dat')
do iproc=0,nproc-1
do iel=1,nel(iproc)
   write(666,*)iproc,iel,procel(iel,iproc),inv_procel(procel(iel,iproc),iproc)
enddo
enddo
close(666)

end if

! check that every element has been assigned
    if (minval(el2proc)==-1) then
       write(6,*) ' ' 
       write(6,*) 'Element(s) not assigned to any processor:', minloc(el2proc)
       stop
    endif

   do iel = 1, neltot
    if (.not.attributed(iel) ) then 
       write(6,*) ' NOT ATTRIBUTED '
       write(6,*) iel, thetacom(iel),solid(iel),fluid(iel)
    endif
   end do

   do iel = 1, neltot
    if (.not.attributed(iel) ) stop
   end do

  if (dump_mesh_info_screen) then
  write(6,*)
  write(6,*)'NUMBER OF ELEMENTS IN EACH SUBDOMAIN:'
  do iproc=0,nproc-1
     write(6,11)iproc, nel_solid(iproc),nel_fluid(iproc),nel(iproc) 
  enddo
  call flush(6)
  write(6,*)
11 format('Proc ',i3, ' has ',i8, ' solid,',i6,' fluid,',i9,' total elements')
  end if

! end DD parameters

deallocate(attributed)

call plot_dd_vtk
  
!  call out_dd ! to plot dd

  end subroutine create_domain_decomposition
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine plot_dd_vtk
use test_bkgrdmodel, only : write_VTK_bin_scal_old

real, allocatable, dimension(:,:) :: mesh1
integer :: iel
character(len=80) :: fname
real(kind=realkind), dimension(:), allocatable :: wel2proc



allocate(mesh1(neltot,2))
do iel=1,neltot
   mesh1(iel,1)=real(sgll(npol/2,npol/2,iel))
   mesh1(iel,2)=real(zgll(npol/2,npol/2,iel))
enddo
fname=trim(diagpath)//'/mesh_domaindecomposition'
allocate(wel2proc(1:neltot))
wel2proc=real(el2proc)
call write_VTK_bin_scal_old(wel2proc,mesh1,neltot,fname)
deallocate(wel2proc)
deallocate(mesh1)
end subroutine plot_dd_vtk
!--------------------------------------------------------------------------



!dk domain_decomposition_theta--------------------------------------
  subroutine domain_decomposition_theta
!
! nel : number of glocal elements, i.e. total number of a processor's elements
! nelmax: maximal number of glocal elements
! neltot: global total number of elements
! procel: 

  integer :: iproc,iiproc,iel
  integer :: iel0_solid,iel0_fluid,mycount
  double precision :: deltatheta
  integer, allocatable :: central_count(:)
  double precision :: pi2
  pi2 = two*dasin(one)

  write(6,*)'     THETA-SLICING as domain decomposition....'

  allocate(central_count(0:nproc-1))

! Create colatitude bounds array for outer shell
  allocate(theta_min_proc(0:nproc-1),theta_max_proc(0:nproc-1))
  theta_min_proc(:) = 0.d0 ; theta_max_proc(:) = 0.d0
  deltatheta = pi2/dble(nproc)
  do iproc = 0, nproc-1
     theta_min_proc(iproc) = dble(iproc  )*deltatheta
     theta_max_proc(iproc) = dble(iproc+1)*deltatheta
  end do

! **************** INNER CUBE **********************
!   call decompose_inner_cube(iproc,theta_min,theta_max)

! define quadratic functions to delineate processor boundaries. 
! Works for nproc=2,4,8,16 at this point
   call decompose_inner_cube_quadratic_fcts(central_count)

! **************** END OF INNER CUBE****************

  iel0_solid = 1; iel0_fluid = 1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iproc = 0, nproc -1
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   mycount = 0

!><><><><><><><><><
   do iel = 1,neltot
!><><><><><><><><><

!   I add the extra requirement that element iel to be in appropriate 
!   theta slice

!   if (fluid(iel)) then
    if (fluid(iel) .and. (thetacom(iel) > theta_min_proc(iproc)) .and.  &
        (thetacom(iel) < theta_max_proc(iproc)) ) then
     mycount = mycount+1
     procel_fluid(mycount,iproc) = iel
     attributed(iel) = .true.
    end if
    if ( mycount == nel_fluid(iproc) ) then
!    I have to change this as well
!    iel0_fluid = iel + 1
     iel0_fluid = 1
     exit
    end if
!><><><><><><><><><
   end do ! iel
!><><><><><><><><><

!!!!!!
   mycount = central_count(iproc)
!!!!!


!  At this stage we have assigned nel_fluid fluid elements
!  to each processor, stored in a procel_fluid(1:nel_fluid,iproc)
!  do iel = iel0_solid,neltot
!  Here we start the loop over solid elements and try to assign them to iproc

!><><><><><><><><><
   do iel = 1,neltot
!><><><><><><><><><
      if (.not.fluid(iel) .and. .not.(attributed(iel)) .and. &
           (thetacom(iel) >= theta_min_proc(iproc))     .and. &
           (thetacom(iel) <= theta_max_proc(iproc))     ) then
         mycount = mycount+1
         procel_solid(mycount,iproc) = iel
         attributed(iel) = .true.
      end if
      if ( mycount == nel_solid(iproc) ) then
!      iel0_solid = iel + 1
         iel0_solid = 1
         if (dump_mesh_info_screen) write(6,*) ' PROC ', iproc ,' has everybody it needs ' ; call flush(6)
         exit
      end if
      if ( mycount == nel_solid(iproc) ) exit

!><><><><><><><><><
   end do
!><><><><><><><><><

   if (mycount< nel_solid(iproc)) then 
      write(6,*)
      write(6,*)'Problem: not all solid elements attributed for proc',iproc,&
                 mycount
      do iiproc=0,nproc-1
         write(6,*)'nel_solid(iproc),centralcount:',&
                   iiproc,nel_solid(iiproc),central_count(iiproc)
      enddo
      stop
   endif

! procel contains
! procel(1:nel_fluid) : the nel_fluid element numbers pertaining to iproc
! procel(nel_fluid+1:nel(iproc)) : the nel_solid solid element numbers 
! belonging to iproc
! Element numbers are defined in a global sense (solid+fluid whole mesh)
   do iel = 1, nel(iproc)
      if (iel <= nel_fluid(iproc) ) then
         procel(iel,iproc) = procel_fluid(iel,iproc)
      else
         procel(iel,iproc) = procel_solid(iel-nel_fluid(iproc),iproc)
         if (procel(iel,iproc)<=0) then
            write(6,*)'PROCEL ZERO!',iproc,eltypeg(iel),iel,nel_fluid(iproc)
            stop
         endif
      end if
      el2proc(procel(iel,iproc))=iproc
      inv_procel(procel(iel,iproc),iproc)=iel
   end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  end do !nproc-1
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! _/_/_/_/_/_/_OUTPUT OF DOMAIN DECOMPOSITION_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  if (dump_mesh_info_files) then
 
  write(6,*)'Writing out the domain decomposition...'; call flush(6)

! central dd only, but checking through whole grid
  open(unit=647,file=diagpath(1:lfdiag)//'/dd_central_sz_ielglob_iproc.dat')
  do iel=1,neltot
     if (eltypeg(iel)=='linear') write(647,14)scom(iel),zcom(iel),iel,&
                                              el2proc(iel)
  enddo
  close(647)
14 format(2(1pe13.3),2(i8))

! THIS ONE WILL BE KINDA LARGE (but still only one per element)
  open(unit=648,file=diagpath(1:lfdiag)//'/dd_sz_ielglob_iproc.dat')
  do iel=1,neltot
     write(648,14)scom(iel),zcom(iel),iel,el2proc(iel)
  enddo
  close(648)

  end if
! _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  end subroutine domain_decomposition_theta
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
subroutine decompose_inner_cube_old(iproc,theta_min,theta_max)

!  Trying to domain decompose the inner region
!  Jan 6th 06
!
  integer, intent(in) :: iproc
  double precision, intent(out) :: theta_min, theta_max
  double precision :: pi2
  pi2 = two*dasin(one)

   if (iproc == 0) then
    theta_min = 0.
!   if (nproc >1) theta_max = 1.25*theta_max_proc(iproc)
    if (nproc >1) theta_max = one*theta_max_proc(iproc)
    if (nproc==1) theta_max = pi2
   elseif (iproc > 0 .and. iproc <= (nproc/2-1) ) then
    if (iproc == ( nproc/2 -1 ) .and. (southern) ) then
     theta_min = zero
     theta_max = half*pi2
    else
     theta_min =     zero
     if (1.15*theta_max_proc(iproc) < half * pi2 ) then
!     theta_max =     1.15*theta_max_proc(iproc)
      theta_max =     one*theta_max_proc(iproc)
     else
      theta_max = half*pi2
     end if
    end if
   elseif (iproc > 0 .and. iproc > (nproc/2-1)  .and. iproc < (nproc-1) ) then
    if (southern) then
     theta_min =      half*pi2
!    theta_max =   1.05*theta_max_proc(iproc)
     theta_max =   one*theta_max_proc(iproc)
    else
     theta_min = theta_min_proc(iproc-1)
     theta_max = theta_max_proc(iproc+1)
    end if
   elseif (iproc == (nproc-1) ) then
    theta_max = pi2
!   if (nproc >1) theta_min = theta_min_proc(iproc-1)
!   if (nproc==1) theta_min = 0.
    theta_min = half*pi2
   end if

end subroutine decompose_inner_cube_old
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
subroutine decompose_inner_cube_quadratic_fcts(central_count)

integer, intent(out) :: central_count(0:nproc-1)
integer :: iproc,is,iz,nproc2
integer :: icount,i2count,iicount,missing
integer :: arclngth,area,CapA,proccount,quadels
integer,allocatable :: proc_central(:,:),num_columns(:),upper_boundary_el(:)
integer,allocatable :: num_columns_hi(:),num_columns_lo(:),num_el(:)
integer,allocatable :: count_assi(:)
double precision :: a,b

if (dump_mesh_info_screen) then 
    write(6,*)
    write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><>'
    write(6,*)'CENTRAL LINEAR DOMAIN: quadratic function decomposition!'
    write(6,*)'ndivs,nproc:',ndivs,nproc
    write(6,*)'==> each processor should have el=',ndivs**2/nproc*2
    write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><>'
    write(6,*)
end if

nproc2=nproc/2-1
if (nproc==1) nproc2=0

allocate(proc_central(1:ndivs,1:ndivs))
proc_central(1:ndivs,1:ndivs)=-1

allocate(num_columns(1:ndivs))
allocate(upper_boundary_el(1:ndivs))
allocate(num_columns_hi(1:ndivs),num_columns_lo(1:ndivs))
allocate(num_el(0:nproc2))

!========================
if (nproc==1) then 
!========================
    proc_central(1:ndivs,1:ndivs)=0
    num_el(0)=ndivs**2

!========================
elseif (nproc==2) then 
!========================
! This is the same as nproc=1 for the North.
! The distinctive southern processor will be taken care of later on...
    proc_central(1:ndivs,1:ndivs)=0
    num_el(0)=ndivs**2

!========================
else ! domain decomposition needed
!========================

    arclngth=int(4*ndivs/nproc)
    area=2*ndivs**2/nproc ! per proc, e.g. north

    if (mod(2*ndivs**2,nproc)/=0) then 
       write(6,*)'Central cube area not divisible into equal areas!'
       stop
    endif

    ! arclngth = processor's section along the boundary of the central region
    if (mod(4*ndivs,nproc)/=0) then 
       write(6,*)'PROBLEM with number of central-region elements and nproc:'
       write(6,*)'ndivs,nproc:',ndivs,nproc
       write(6,*)'ndivs (number of northern elements in one direction)'
       write(6,*)'needs to be multiple of nproc/4...'
       stop
    endif

    ! area well defined?
    if (mod(2*ndivs**2,nproc)/=0) then 
       write(6,*)'PROBLEM with number of central-region elements and nproc:'
       write(6,*)'els,nproc:',ndivs**2,nproc
       write(6,*)'number of elements needs to be multiple of nproc/2...'
       stop
    endif

    !%%%%%%%%%%%% equator mynum=nproc/2-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! define quadratic function with area = ndivs**2/nproc
    ! starting at arclength = 4 ndivs/nproc for s=ndivs
    ! 
    a = (3.d0/2.d0)*(dble(arclngth)-dble(area)/dble(ndivs))/dble(ndivs)**2
    b = ((3.d0/2.d0)*dble(area)/dble(ndivs)-dble(arclngth)/2.d0) 
    icount=0;

    if (nproc>4) then
    ! make sure the column at ndivs equals arclength
       proc_central(ndivs,arclngth)=nproc2
       icount=icount+1
    endif


    do is=ndivs,1,-1
       do iz=1,ndivs
          if (iz < a*is**2 + b) then
             proc_central(is,iz)=nproc2
             icount=icount+1
          endif

    ! keep the first two elements off the spherical/linear boundary
    ! belonging to the same proc
          if (is==ndivs-1 .and. proc_central(ndivs,iz)==nproc2 .and. &
               proc_central(is,iz)/=nproc2 ) then 
             proc_central(ndivs-1,iz)=nproc2
             icount=icount+1
          endif
       enddo
    enddo

    missing=area-icount

    if (missing<0) then
       write(6,*)'Problem: included more elements than needed in equat. element!'
       write(6,*)'included/needed:',icount,area
       stop

    elseif (missing>0) then
       ! fill the missing elements
       if (dump_mesh_info_screen)& 
       write(6,*)'Equatorial processor nproc/2-1: Filling missing els:',missing
       icount=0
          do while (icount<missing)
          is=ndivs
          do while (is>2)
             is=is-1
             do iz=1,ndivs
             if (proc_central(is,iz)==nproc2) then
                if (proc_central(is-1,iz)/=nproc2) then
                   proc_central(is-1,iz)=nproc2
                   icount=icount+1
                   is=is-1
                   if (icount==missing .or. is<2) exit
                endif
                if (icount==missing .or. is<2) exit
             endif
             if (icount==missing .or. is<2) exit
          enddo
          if (icount==missing ) exit
       enddo
       if (icount==missing ) exit
       enddo
    endif ! if missing

    ! duplicate this distribution to next processors as # columns per is
    num_columns(1:ndivs)=0
    do is=1,ndivs
       do iz=1,ndivs-1
          if (proc_central(is,iz)==nproc2) then
             num_columns(is)=num_columns(is)+1
             if (proc_central(is,iz+1)/=nproc2) upper_boundary_el(is)=iz
          endif
       enddo
    enddo

    ! Check if proc has right amount of elements
    call check_my_els(nproc2,proc_central)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (nproc==4) then 
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       
       do is=1,ndivs
          do iz=1,ndivs
             if (proc_central(is,iz)==-1) then
                proc_central(is,iz)=0
             end if
          end do
       end do

    ! Check if proc has right amount of elements
    call check_my_els(0,proc_central)

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else ! nproc>4
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        !%%%%%%%%% processors along the s=ndivs vertical edge of the central cube
        !%%%%%%%%% nproc/2-1 down to corner/45-deg proc nproc/4

        proccount=0
        num_columns_hi(1:ndivs)=num_columns(1:ndivs)

        !------------------------------
        do iproc=nproc/2-2,nproc/4,-1
        !------------------------------

           proccount=proccount+1

           num_columns_lo = num_columns_hi+1
           num_columns_hi =  num_columns_hi  + num_columns

        ! try to avoid jumps growing with each new processor
        !   do is  = 2,ndivs-2
           do is  = 2,ndivs-1
              if (num_columns_hi(is-1) < num_columns_hi(is)-1) then

                 num_columns_hi(is-1) = num_columns_hi(is-1) + 1
                 num_columns_hi(is) = num_columns_hi(is) - 1

        ! avoid diagonally touching processor boundaries, i.e. elements known to 
        ! three processors 
                 if (num_columns_hi(is)<=num_columns_lo(is+1)+1) then 
                    num_columns_hi(is)=num_columns_hi(is)+1
                    num_columns_hi(is+1)=num_columns_hi(is+1)-1
                 endif

              endif
                 
        ! try to avoid length plateaus followed by abrupt steps
        !      if (num_columns_hi(is-1) == num_columns_hi(is) .and. &
        !          num_columns_hi(is+1) == num_columns_hi(is) ) then
        !
        !         num_columns_hi(is-1) = num_columns_hi(is-1) - 1
        !         num_columns_hi(is+1) = num_columns_hi(is+1) + 1
        !      endif

        ! try to avoid jumps over two
              if (num_columns_hi(is-1) == num_columns_hi(is) .and. &
                  num_columns_hi(is+1) >= num_columns_hi(is)+2 ) then

                 num_columns_hi(is) = num_columns_hi(is) + 1
                 num_columns_hi(is+1) = num_columns_hi(is+1) - 1

              endif

        ! try to avoid drops in the positive is direction
              if (num_columns_hi(is+1) < num_columns_hi(is)) then

                 num_columns_hi(is) = num_columns_hi(is) - 1
                 num_columns_hi(is+1) = num_columns_hi(is+1) + 1

              endif

        ! try to avoid jumps over two (again...)
              if (num_columns_hi(is-1) == num_columns_hi(is) .and. &
                  num_columns_hi(is+1) >= num_columns_hi(is)+2 ) then

                 num_columns_hi(is) = num_columns_hi(is) + 1
                 num_columns_hi(is+1) = num_columns_hi(is+1) - 1

              endif

           enddo

        ! fill with belongings...
           do is=1,ndivs
              proc_central(is,num_columns_lo(is):num_columns_hi(is))=iproc
           enddo

        ! Check if proc has right amount of elements
        call check_my_els(iproc,proc_central)

        !------------------------------
        enddo ! iproc
        !------------------------------

        !%%%%%%%%%%%% mynum 0 (northernmost) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        proc_central(1:arclngth,ndivs/2+1:ndivs)=0

        ! Check if proc has right amount of elements
        call check_my_els(0,proc_central)

    ! >>>>>>>>>>>>>>>>
    endif ! if nproc>4
    ! >>>>>>>>>>>>>>>>

    ! Now start working on the procs on the upper half%%%%%%%%%%%%%%%%%%%%%%%%%%

    !<<<<<<<<<<<<<<<<<<<<<<
    if (nproc==8) then  
    !<<<<<<<<<<<<<<<<<<<<<<

        do is=1,ndivs
           do iz=1,ndivs
              if (proc_central(is,iz)==-1) then
                 proc_central(is,iz)=1
              endif
           enddo
        enddo

        ! Check if proc has right amount of elements
        call check_my_els(1,proc_central)

    !<<<<<<<<<<<<<<<<<<<<<<
    elseif (nproc==16) then
    !<<<<<<<<<<<<<<<<<<<<<<

        !%%%%%% mynum=1, second processor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! upper rectangle
        proc_central(arclngth+1:floor(1.5d0*arclngth),ndivs/2+1:ndivs)=1
        i2count=floor(dble(arclngth)/2.d0*dble(ndivs)/2.d0)

        ! right slope along rectangle
        do is=floor(1.5d0*arclngth)+1,2*arclngth
           do iz=ndivs/2+1,ndivs
              if (iz>ndivs/arclngth*is-ndivs) then
                 proc_central(is,iz)=1
                 i2count=i2count+1
              endif
           enddo
        enddo
        proccount=0

        do is=1,floor(3.d0*arclngth/2.d0)
           do iz=floor(ndivs/3.d0),ndivs/2
              if (iz>ceiling(dble(ndivs)/9.d0/dble(arclngth)*is)+dble(ndivs)/3.d0 .and. &
                  iz<=ndivs .and. proc_central(is,iz)==-1 ) then
                 proc_central(is,iz)=1
                 i2count=i2count+1
              endif
           enddo
        enddo

        missing=area-i2count

        if (missing<0) then 
           write(6,*)'Problem with second processor!'
           write(6,*)'included too many elements:',i2count,'instead of',area
           stop
        endif

        i2count=0
        if (missing>0) then
           if (dump_mesh_info_screen)&
           write(6,*)'Second processor 1: Filling missing els:',missing
           do while (i2count <missing)
              do is=floor(1.5d0*arclngth)+1,2*arclngth
                 do iz=ndivs/2+1,ndivs
                    if (proc_central(is,iz)==1 .and. proc_central(is,iz-1)/=1) then
                       proc_central(is,iz-1)=1
                       i2count=i2count+1
                       if (i2count==missing) exit
                    endif
                    if (i2count==missing) exit
                 enddo
              if (i2count==missing) exit
              enddo
              if (i2count==missing) exit
           enddo
        endif !missing

        ! Check if proc has right amount of elements
        call check_my_els(1,proc_central)

        ! proc in the 45-deg corner nproc/4-1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iproc=nproc/4-1
        icount=0

        ! first fill in the procs do s> ndivs-arclngth (the part in the "corner")
        CapA = ndivs - arclngth
        do is=CapA,ndivs
           do iz=1,ndivs
              if  (proc_central(is,iz)==-1) then
                 proc_central(is,iz)=iproc
                 icount=icount+1
              endif
           enddo
        enddo

        ! quadratic part s< ndivs-arclngth

        ! first fill one row above previous processor to make sure it goes to axis
        do is=1,ndivs-arclngth+1
           iz=1
           do while (iz<ndivs)
              iz=iz+1
        !   do iz=2,ndivs
              if (proc_central(is,iz)==-1) then

                 if (proc_central(is,iz-1)==iproc+1 .or. &
                     proc_central(is+1,iz)==iproc+1 .or. &
                     proc_central(is+1,iz-1)==iproc+1) then 

        !            proc_central(is,iz)=iproc
        !            proc_central(is,iz+1)=iproc
        !            iz=ndivs
        !            icount=icount+2

        ! NEED TO MAKE SURE THAT THERE IS MORE THAN ONE ELEMENT DIAGONALLY!
                    proc_central(is,iz)=iproc
                    proc_central(is,iz+1)=iproc
                    icount=icount+2

        ! fill diagonal element 
                    if (is>1 ) then
                       if ( proc_central(is-1,iz)/=-1 .and. &
                            proc_central(is-1,iz+1)==-1) then
                          proc_central(is-1,iz+1)=iproc
                          if (dump_mesh_info_screen) &
                          write(6,*)'Filled diagonal element in corner proc:',is,iz
                          icount=icount+1
                       endif
                    endif

                    if (proc_central(is+1,iz+1)==iproc+1 ) then 
                       proc_central(is,iz+2)=iproc
                       if (dump_mesh_info_screen) &
                       write(6,*)'Filled diagonal element in corner proc:',is,iz
                       icount=icount+1
                    endif
                    iz=ndivs

                 end if
              endif
           enddo
        enddo

        quadels=area-icount

        ! now check by quadratic function to fill up the rest until quadels
        a = 3.d0/2.d0*ndivs**2/CapA**3 * (1.d0/2.d0-2.d0/dble(nproc) )
        b= 3.d0/2.d0*ndivs**2/CapA*(2.d0/dble(nproc) - 0.5d0) + ndivs

        iicount=0
        do while (iicount<quadels)
            do is=ndivs-arclngth,1,-1
               do iz=1,ndivs
                  if (iz<a*is**2 + b .and. proc_central(is,iz)==-1) then
                     proc_central(is,iz)=iproc
                     iicount=iicount+1
                     if (iicount>=quadels) exit
                  endif
                  if (iicount>=quadels) exit
               enddo
               if (iicount>=quadels) exit
            enddo
        enddo

        ! Check if proc has right amount of elements
        call check_my_els(iproc,proc_central)

        !<<<<<<<<<<<<<<<<<<<<<<
        if (nproc==16) then
        !<<<<<<<<<<<<<<<<<<<<<<

            ! Filling the rest with values to the last proc 
               iproc=2
               do is=1,ndivs
                  do iz=1,ndivs
                     if (proc_central(is,iz)==-1) then
                        proc_central(is,iz)=iproc
                     endif
                  enddo
               enddo

            ! Check if proc has right amount of elements
            call check_my_els(iproc,proc_central)

        !<<<<<<<<<<<<<<<<<<<<<<
        endif ! nproc=16
        !<<<<<<<<<<<<<<<<<<<<<<

    else

       write(6,*)
       write(6,*)'PROBLEM: central cube decomposition not yet done for nproc',nproc
       stop

    !<<<<<<<<<<<<<<<<<<<<<<
    endif ! nproc >= 16
    !<<<<<<<<<<<<<<<<<<<<<<

!========================
endif
!========================

allocate(count_assi(0:nproc2))
count_assi=0
! %%%%%% %%%%%%count respective processors elements %%%%%%%%%%%
num_el(0:nproc2)=0
do is=1,ndivs
   do iz=1,ndivs
      do iproc=0,nproc2
         if (proc_central(is,iz)==iproc) then 
            num_el(iproc)=num_el(iproc)+1
            count_assi(iproc)=count_assi(iproc)+1
         endif
      enddo
      if (proc_central(is,iz)<0 .or. proc_central(is,iz)> nproc2) then
         write(6,*)'Problem:',is,iz,'has no processor!'
         stop
      endif
   enddo
enddo

do iproc=0,nproc2
if (count_assi(iproc)/= ndivs**2/(nproc2+1)) then 
   write(6,*)'Problem: Not every element is assigned to processor',iproc
   write(6,*)'Counted assigned els/total els:',count_assi(iproc), &
                                               ndivs**2/(nproc2+1)
   if (iproc<nproc2) write(6,*)'els for other procs:',&
                     count_assi(iproc+1:nproc2)
   stop
endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(dump_mesh_info_screen) then 
 write(6,*)
  do iproc=0,nproc2
    write(6,12)iproc,num_el(iproc)
  enddo
12 format('Central cube: proc',i3' has',i6' elements')
end if

! connect these processor dependencies to the global element numbering scheme
  central_count(0:nproc-1)=0
  do iz = 1, ndivs
   do is = 1, ndivs
      if (proc_central(is,iz)/=-1) then
         central_count(proc_central(is,iz))=&
                      central_count(proc_central(is,iz))+1

         attributed(central_is_iz_to_globiel(is,iz))=.true.
         attributed(central_is_iz_to_globiel(is,iz)+neltot/2)=.true.
      else
         write(6,*)'Unassigned element in the central cube!'
         write(6,*)'is,iz:',is,iz
         stop
      endif

     if (neltot_solid>0 ) then 
     procel_solid(central_count(proc_central(is,iz)),proc_central(is,iz)) = &
                                                central_is_iz_to_globiel(is,iz)
     else
        procel_fluid(central_count(proc_central(is,iz)),proc_central(is,iz)) = &
                                                      central_is_iz_to_globiel(is,iz)
     endif

! South: inverted copy
     if (nproc>1) then 
        if (neltot_solid>0 ) then 
           procel_solid(central_count(proc_central(is,iz)), &
                nproc-1-proc_central(is,iz)) = &
                central_is_iz_to_globiel(is,iz)+neltot/2
        else
           procel_fluid(central_count(proc_central(is,iz)), &
                nproc-1-proc_central(is,iz)) = &
                central_is_iz_to_globiel(is,iz)+neltot/2
        endif
     endif
    enddo
  enddo

! South: 
if (nproc>1) then
   do iproc=0,nproc/2-1
      central_count(nproc-iproc-1) = central_count(iproc)
   enddo
endif
     
! special case one processor... still needs to count the south!
if (nproc==1) then 
  do iz = 1, ndivs
   do is = 1, ndivs
        central_count(proc_central(is,iz))=&
                      central_count(proc_central(is,iz))+1
        if (neltot_solid>0 ) then         
           procel_solid(central_count(proc_central(is,iz)),0) = &
                central_is_iz_to_globiel(is,iz)+neltot/2
        else
           procel_fluid(central_count(proc_central(is,iz)),0) = &
                central_is_iz_to_globiel(is,iz)+neltot/2
        endif
    enddo
  enddo
endif

! check if all central-cube elements are assigned
  do is=1,neltot
     if (eltypeg(is)=='linear' ) then 
        if (.not. attributed(is)) then 
           write(6,*)
           write(6,*)'Problem: Central cube element not assigned!',is
           stop
        endif
     endif
  enddo

! _/_/_/_/_write out the central cube decomposition_/_/_/_/_/_/_/_/_/_/_/_/_/_/

if(dump_mesh_info_files) then 
open(unit=10008,file=diagpath(1:lfdiag)//'/central_locind_locsz_iproc.dat')
do iz = 1, ndivs
   do is = 1, ndivs
      write(10008,13)is,iz,s_arr(is,iz),z_arr(is,iz),proc_central(is,iz)
   enddo
enddo
do iz = 1, ndivs
   do is = 1, ndivs
      write(10008,13)is,iz,s_arr(is,iz),-z_arr(is,iz),nproc-1-proc_central(is,iz)
   enddo
enddo
close(10008)
13 format(2(i4),2(f9.3),i4)
end if

! _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  if (dump_mesh_info_screen) then
   write(6,*)
   write(6,*)'<><><><> Finished the central domain decomposition!<><><><>'
   write(6,*)
   call flush(6)
  end if

end subroutine decompose_inner_cube_quadratic_fcts
!--------------------------------------------------------------------------
!
!dk check_my_els-----------------------------------------------------------
subroutine check_my_els(iproc,proc_central)

integer, intent(in) :: iproc,proc_central(1:ndivs,1:ndivs)
integer :: procelcount,is,iz

! Check if proc has right amount of elements
  procelcount=0
  do is=1,ndivs
     do iz=1,ndivs
        if (proc_central(is,iz)==iproc) procelcount=procelcount+1
     enddo
  enddo
  if (procelcount/=2*ndivs**2/nproc) then 
     write(6,*)
     write(6,12)iproc,procelcount
     write(6,*)'Needed:',2*ndivs**2/nproc
     stop
  else 
     if (dump_mesh_info_screen) then 
      write(6,*)
      write(6,*)'>>>', iproc,' has right amount of elements:',2*ndivs**2/nproc
      call flush(6)
     end if
  endif

12 format('Problem: Processor',i4,' has ',i6,' elements!')

end subroutine check_my_els
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!dk create_domain_decomposition_naive--------------------------------------
  subroutine create_domain_decomposition_naive
  integer :: iel,iproc,nelmax
 
  write(6,*) ' GLOBAL NUMBER OF ELEMENTS IS ', neltot

! DD parameters
  nproc = nproc_target 
  allocate(nel(0:nproc-1))

! Simple division a la iproc <-- 1:neltot/nproc
  if ( mod(neltot,nproc) == 0) then
   do iproc = 0, nproc-1
    nel(iproc) = neltot / nproc
   end do
  else
   do iproc = 0, nproc -1
     if (iproc < nproc -1 )  nel(iproc) = neltot / (nproc-1)
     if (iproc == (nproc-1)) nel(iproc) = neltot - SUM (nel(0:nproc-2))
   end do
  end if
  nelmax = maxval(nel)
  allocate(procel(nelmax,0:nproc-1))
  do iproc = 0, nproc -1
   do iel = 1, nel(iproc)
    procel(iel,iproc) =   iel + SUM(nel(0:iproc-1))
   end do
  end do
  write(6,*) ' Number of elements per proc '
  do iproc = 0, nproc-1
   write(6,*) ' PROC ', iproc, nel(iproc)
!  write(6,*) '   PROC ', iproc
!  write(6,*) (procel(iel,iproc),iel=1,nel(iproc))
  end do
! end DD parameters
!!

  call out_dd

  end subroutine create_domain_decomposition_naive
!--------------------------------------------------------------------------


!dk check_nproc------------------------------------------------------------
  subroutine check_nproc(np)
  use data_coarse
  integer :: np,nb,itest,ritest
  nb = ns_ib ! Number of elements at the inner spherical boundary (defined while
             ! generating the skeleton)
  ritest = dble(2*nb)/dble(np)
   itest = 2*nb/np
  if (ritest /= dble(itest) ) then 
   write(6,*) ritest,itest,2*nb,np
   write(6,*) 2*nb,np
   write(6,*) ' Number of processes not compliant with mesh topology '
   stop
  end if
  end subroutine check_nproc
!--------------------------------------------------------------------------
!
!  
!dk out_dd---------------------------------------------------------------
  subroutine out_dd
  use data_gllmesh
  integer :: iel,ipol,jpol,iproc,ielg
  do iproc = 0, nproc -1
   do iel = 1, nel(iproc)
    ielg = procel(iel,iproc)
    do jpol = 0, npol
     do ipol = 0, npol
!   do jpol = 1, npol-1
!    do ipol = 1, npol-1
      write(55,*) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg),iproc
     end do
    end do
    do jpol = 0, npol, npol
     do ipol = 0, npol
      write(55,*) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg),-1
     end do
    end do
    do ipol = 0, npol, npol
     do jpol = 0, npol
      write(55,*) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg),-1
     end do
    end do
    do ipol = 0, npol-1
     do jpol = 0, npol-1
      write(56,20) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg)
      write(56,20) sgll(ipol+1,jpol,ielg),zgll(ipol+1,jpol,ielg)
      write(56,20) sgll(ipol+1,jpol+1,ielg),zgll(ipol+1,jpol+1,ielg)
      write(56,20) sgll(ipol,jpol+1,ielg),zgll(ipol,jpol+1,ielg)
!     if ((jpol == (npol-1)) .and. (ielg < neltot) ) then
       write(56,20) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg)
      write(56,22)
!     else
!      write(56,20) sgll(ipol,jpol,ielg),zgll(ipol,jpol,ielg)
!     end if
     end do
    end do
   end do
  end do
20 format(1pe12.5,2x,1pe12.5)
!21 format(1pe12.5,2x,1pe12.5,' > ')
22 format('>')
  end subroutine out_dd
!------------------------------------------------------------------------
!
  end module parallelization 
