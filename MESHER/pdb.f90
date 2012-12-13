!============
  module pdb
!============
  use data_mesh
  use data_numbering
  use data_spec
  use data_diag
  use data_pdb
  use data_bkgrdmodel
  use data_grid

  implicit none
  public :: create_pdb
  private
  contains

!dk create_pdb-----------------------------------------------------------
subroutine create_pdb
  ! Wrapper routine to define everything in the parallel global and solid/fluid!
  ! worlds and dump the database for each processor. 
  ! Additionally creating the header mesh_params.h containing all major mesh 
  ! sizes for the solver and the solver routine unrolled_loops.f90 that depends 
  ! on the polynomial order given in the mesher.
  ! Whenever DATABASE occurs as a comment below, those parameters will be saved 
  ! into the solver database.
  ! 
  ! Terminology: 
  ! global - counting through the entire solid+fluid domain across processors
  ! slobal - counting through the entire solid domain across processors
  ! glocal - counting through a single processor's entire solid+fluid domain
  ! slocal - counting through a single processor's solid domain
  ! accordingly for fluid: flobal, flocal, and colloquially for sflobal, sflocal.
  
  use data_gllmesh, only : sgll, zgll

  integer   :: nelmax
  if (dump_mesh_info_screen) then
    write(6,*)
    write(6,*)' ||||||||||||||| CREATING THE PARALLEL DATABASE ||||||||||||||||'
    write(6,*)
  end if

  if (allocated(iglob)) deallocate(iglob)

  write(6,*) '  define glocal numbering....'; call flush(6)
  call define_glocal_numbering ! needs sgll,zgll, creates igloc

  ! Solid-fluid distinction
  write(6,*) '  define solflu coordinates....'; call flush(6)
  call define_sflocal_coordinates ! needs sgll, zgll 
                                  ! creates procel_solidp     procel_fluidp 
                                  !         inv_procel_solidp inv_procel_fluidp  
  nelmax = maxval(nel)

  write(6,*) '  define axial elems....'; call flush(6)
  call define_axial_elem ! needs sgll, sgll_solid, sgll_fluid
  

  write(6,*) '  define solflu numbering....'; call flush(6)
  call define_sflocal_numbering   ! needs sgll, zgll
                                  ! creates igloc_solid, igloc_fluid
  deallocate(sgll,zgll)

  write(6,*) '  define flobal2flocal etc....'; call flush(6)
  call define_sflobal2sflocal     ! needs iglob_solid, iglob_fluid
                                  ! needs igloc_fluid, igloc_fluid

  write(6,*) '  define search sflobal index....'; call flush(6)
  call define_search_sflobal_index ! needs iglob_solid, iglob_fluid

  write(6,*) '  partition sflobal index....'; call flush(6)
  call partition_sflobal_index
 
  write(6,*) '  define local bdry elems....'; call flush(6)
  call define_local_bdry_elem


  ! For compliance with solver, we need the control points info for each process
  write(6,*) '  generate processor serendipity....'; call flush(6)
  call generate_serendipity_per_proc(sg,zg) ! needs sgp, zgp

  write(6,*) '  define element type....'; call flush(6)
  call define_element_type ! 

  ! Write out mesh database
  write(6,*) '  write database....'; call flush(6) 
  call write_db

  ! Check that database is written out correctly
  write(6,*) '  read database....'; call flush(6) 
  ! af That test is not necessary
  call read_db 

  write(6,*) ' create static header mesh_params.h ....';call flush(6) 
  call create_static_header

  write(6,*) ' create solver routine unrolled_loops.f90 ....';call flush(6) 
  call create_unrolled_loop_routines

end subroutine create_pdb
!------------------------------------------------------------------------
  
!dk define_sflocal_coordinates-------------------------------------------
subroutine define_sflocal_coordinates
!
! NOTE: procel_fluid and procel_solid map the GLOBAL element numbers!!!
! i.e.: procel(1:nel_fluid(iproc),iproc)=procel_fluid(:,iproc)
!       procel(nel_fluid(iproc)+1:nel(iproc),iproc)= procel_solid(:,iproc)
!
! Also re-defining procel_solid here, and defining new mapping
! from solf/flu parallel elem to global parallel elem number
! This will be the output to the solver to become ielsolid and ielfluid!!
!
! procel_solidp: mapping from slocal to glocal elem number
! inv_procel_solidp: mapping from glocal to slocal elem number
!
! DATABASE: procel_solidp,procel_fluidp
!
  use data_gllmesh
  integer           :: iproc, iel, ielg, jpol, ipol
  integer           :: nelmax_solid, nelmax_fluid
  !double precision,allocatable :: rsol(:,:,:,:),rflu(:,:,:,:)
  double precision  :: rsol_min, rsol_max, rflu_min, rflu_max

  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)

  ! initialize to crazy values
  rsol_max = -1e30
  rflu_max = -1e30
  rsol_min = 1e30
  rflu_min = 1e30

! Solid
  do iproc = 0, nproc-1
   do iel = 1, nel_solid(iproc)
    ielg = procel_solid(iel,iproc)
    do jpol = 0, npol
     do ipol = 0, npol
      rsol_max = max(rsol_max, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
      rsol_min = min(rsol_min, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
     end do
    end do
   end do
  end do


  write(6,*)'HAVE FLUID?',have_fluid

! Fluid
  if (have_fluid) then
     do iproc = 0, nproc-1
      do iel = 1, nel_fluid(iproc)
       ielg = procel_fluid(iel,iproc)
       do jpol = 0, npol
        do ipol = 0, npol
         rflu_max = max(rflu_max, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
         rflu_min = min(rflu_min, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
        end do
       end do
      end do
     end do
  endif ! have_fluid

  rsol_max = sqrt(rsol_max)
  rsol_min = sqrt(rsol_min)
  rflu_max = sqrt(rflu_max)
  rflu_min = sqrt(rflu_min)

! Min/max values
  if (dump_mesh_info_screen) then
   write(6,*)'Solid-fluid coordinates:'
   if (have_solid) &
      write(6,*)'Min/max radius solid:', rsol_max, rsol_min
   if (have_fluid) &
      write(6,*)'Min/max radius fluid:', rflu_max, rflu_min
  end if

! Define processor-specific mapping: 
! procel_solidp: given slocal el, return glocal el number
! inv_procel_solidp, given glocal, return slocal el number

  allocate(procel_solidp(nelmax_solid,0:nproc-1))
  allocate(procel_fluidp(nelmax_fluid,0:nproc-1))
  allocate(inv_procel_solidp(maxval(nel),0:nproc-1))
  allocate(inv_procel_fluidp(maxval(nel),0:nproc-1))

  do iproc=0,nproc-1
     do iel=1,nel_solid(iproc)
        procel_solidp(iel,iproc)=inv_procel(procel_solid(iel,iproc),iproc)
        inv_procel_solidp(procel_solidp(iel,iproc),iproc)=iel
     enddo
     do iel=1,nel_fluid(iproc)
        procel_fluidp(iel,iproc)=inv_procel(procel_fluid(iel,iproc),iproc)
        inv_procel_fluidp(procel_fluidp(iel,iproc),iproc)=iel
     enddo
  enddo

! check procel_solidp and inverses
  do iproc=0,nproc-1
     do iel=1,nel_solid(iproc)
        if (iel .ne. inv_procel_solidp(procel_solidp(iel,iproc),iproc) ) then 
           write(6,*)'PROBLEM in procel_solidp vs. inv_procel_solidp!'
           write(6,*)'iproc,iel,inv(procel_solidp(iel)):',iproc,iel, &
                      inv_procel_solidp(procel_solidp(iel,iproc),iproc)
           stop
        endif

        if ( .not. solid(procel( procel_solidp(iel,iproc),iproc) ) ) then
           write(6,*)'PROBLEM with procel or procel_solidp or both!'
           write(6,*)'slocal,glocal,global el:',iel,procel_solidp(iel,iproc), &
                     solid(procel(procel_solidp(iel,iproc),iproc))
           stop
        endif
     enddo
  enddo

! check procel_fluidp and inverses
  do iproc=0,nproc-1
     do iel=1,nel_fluid(iproc)
        if (iel .ne. inv_procel_fluidp(procel_fluidp(iel,iproc),iproc) ) then
           write(6,*)'PROBLEM in procel_fluidp vs. inv_procel_fluidp!'
           write(6,*)'iproc,iel,inv(procel_fluidp(iel)):',iproc,iel, &
                      inv_procel_fluidp(procel_fluidp(iel,iproc),iproc)
           stop
        endif

        if ( .not. fluid(procel(procel_fluidp(iel,iproc),iproc)) ) then
           write(6,*)'PROBLEM with procel or procel_fluidp or both!'
           write(6,*)'flocal,glocal,global el:',iel,procel_fluidp(iel,iproc), &
                     fluid(procel(procel_fluidp(iel,iproc),iproc))
           stop
        endif
     enddo
  enddo

! check procel_solidp/procel_fluidp and inverses per global search 
  do iproc=0,nproc-1
     do iel=1,nel(iproc)
        if (solid(procel(iel,iproc))) then
           if (iel .ne. procel_solidp(inv_procel_solidp(iel,iproc),iproc) ) then 
              write(6,*)'PROBLEM in procel_solidp vs. inv_procel_solidp!'
              write(6,*)'iproc,iel,procel_solidp(inv(iel)):',iproc,iel, &
                   procel_solidp(inv_procel_solidp(iel,iproc),iproc)
              stop
           endif

        elseif (fluid(procel(iel,iproc))) then
!          write(1009,*)iproc,iel,inv_procel_fluidp(iel,iproc)
           if (iel .ne. procel_fluidp(inv_procel_fluidp(iel,iproc),iproc) ) then 
              write(6,*)'PROBLEM in procel_fluidp vs. inv_procel_fluidp!'
              write(6,*)'iproc,iel,procel_fluidp(inv(iel)):',iproc,iel, &
                         procel_fluidp(inv_procel_fluidp(iel,iproc),iproc)
              stop
           endif
        else
           write(6,*)'If not solid nor fluid, what then??' 
           write(6,*)'iel (glocal),iproc:',iel,iproc
           stop
        endif
     enddo
  enddo

  !deallocate(rsol,rflu)

end subroutine define_sflocal_coordinates
!------------------------------------------------------------------------

!dk define_glocal_numbering----------------------------------------------
subroutine define_glocal_numbering
!
! The glocal numbering for each processor's entire domain.
! The glocal number bookkeeping is not passed to the database, 
! as we only need solid and fluid global numbers.
!
  use numbering 
  use data_gllmesh, only : sgll,zgll
  integer :: nelmax,nelp
  integer :: npointotp,wnglob
  integer :: iproc,ipol,jpol,iel,ipt,ielg
!  real, dimension(:), allocatable :: wsgll,wzgll
  double precision, dimension(:), allocatable :: wsgll,wzgll
  integer, dimension(:), allocatable :: wigloc,wloc
  logical, dimension(:), allocatable :: wifseg

! valence test
  double precision, dimension(:), allocatable :: uglob2
  double precision, dimension(:), allocatable :: val(:,:,:)
  integer :: idest,i
  integer :: valnum_cent(6),totvalnum_cent
  integer :: valnum_semi(6),totvalnum_semi
  character(len=4) :: appiproc

  nelmax = maxval(nel)
  allocate (igloc(nelmax*(npol+1)**2,0:nproc-1))
  allocate (nglobp(0:nproc-1))

!===================
  do iproc = 0, nproc-1
!===================

   nelp = nel(iproc)  
   npointotp = nelp*(npol+1)**2
   allocate(wsgll(npointotp),wzgll(npointotp))
   do iel = 1, nelp
    ielg = procel(iel,iproc)
    do jpol = 0, npol
     do ipol = 0, npol
      ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1 
      wsgll(ipt) = sgll(ipol,jpol,ielg)
      wzgll(ipt) = zgll(ipol,jpol,ielg)
     end do
    end do
   end do

   allocate(wigloc(npointotp))
   allocate(wifseg(npointotp))
   allocate(wloc(npointotp))
   call get_global(nelp,wsgll,wzgll,wigloc,wloc,wifseg,wnglob,&
                   npointotp,NGLLCUBE,NDIM)
   do ipt = 1, npointotp
    igloc(ipt,iproc) = wigloc(ipt)
   end do
   nglobp(iproc) = wnglob

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! valence: test global numbering
!
  allocate(uglob2(nglobp(iproc))) 
  allocate(val(0:npol,0:npol,nel(iproc)))

! valence test, equivalent to how assembly is used in the solver
! use script plot_proc_valence.csh to generate GMT valence grids for each proc
! glocally and zoomed into r<0.2 ( denoted as *_central )
 val(:,:,:) = 1.0

 valnum_cent=0
 valnum_semi=0

 uglob2(:) = 0.d0
 do iel = 1, nel(iproc)
  do ipol = 0, npol
   do jpol = 0, npol
    ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
    idest = wigloc(ipt)
    uglob2(idest) = uglob2(idest) + val(ipol,jpol,iel)
   end do
  end do
 end do
!QT
 do iel = 1, nel(iproc)
  do ipol = 0, npol
   do jpol = 0, npol
    ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
    idest = wigloc(ipt)
    val(ipol,jpol,iel) = uglob2(idest)

! Check valence/global number inside central region
     if (eltypeg(procel(iel,iproc))=='linear') then 
        do i=1,6 !possible valences
           if (val(ipol,jpol,iel)==i) valnum_cent(i)=valnum_cent(i)+1
        enddo
     endif

     if (eltypeg(procel(iel,iproc))=='semino' .or. &
         eltypeg(procel(iel,iproc))=='semiso') then 
        do i=1,6 !possible valences
           if (val(ipol,jpol,iel)==i) valnum_semi(i)=valnum_semi(i)+1
        enddo
     endif

   end do
  end do
 end do
 deallocate(uglob2,val)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 deallocate(wloc)
 deallocate(wifseg)
 deallocate(wigloc)
 deallocate(wzgll,wsgll)

 totvalnum_cent = 0
 totvalnum_semi = 0
 do i=1,6
    totvalnum_cent = totvalnum_cent + valnum_cent(i)/i
    totvalnum_semi = totvalnum_semi + valnum_semi(i)/i
 enddo

 if (dump_mesh_info_screen) then 
  write(6,*)
  write(6,*)iproc,'glocal number      :',nglobp(iproc)
  write(6,*)iproc,'central region only:',totvalnum_cent
  write(6,*)iproc,'semicurved els only:',totvalnum_semi
  write(6,*)iproc,'everywhere else    :',nglobp(iproc)-&
                                         totvalnum_cent-totvalnum_semi
 end if
!===================
end do
!===================

if (dump_mesh_info_screen) then
  do iproc = 0, nproc-1
   write(6,*) 'proc.,nglob:', iproc, nglobp(iproc)
  end do
  write(6,*) 'SUM = ',SUM(nglobp(0:nproc-1))
  write(6,*)
end if

  end subroutine define_glocal_numbering
!------------------------------------------------------------------------

!dk define_sflocal_numbering----------------------------------------------
  subroutine define_sflocal_numbering
!
! The sflocal numbering for each processor's solid & fluid domains.
! This will be passed to the solver as the defining bookkeeping arrays 
! for the assembly procedure.
! DATABASE: igloc_solid,igloc_solid
!
  use numbering 
  use data_gllmesh, only : sgll, zgll
  integer :: nelmax_solid,nelmax_fluid,nelp_solid,nelp_fluid
  integer :: npointotp_solid,npointotp_fluid,wnglob_solid,wnglob_fluid
  integer :: iproc,ipol,jpol,iel,ipt, ielg
  double precision, dimension(:), allocatable :: wsgll_solid,wzgll_solid
  double precision, dimension(:), allocatable :: wsgll_fluid,wzgll_fluid

  integer, dimension(:), allocatable :: wigloc_solid,wigloc_fluid
  integer, dimension(:), allocatable :: wloc_solid,wloc_fluid
  logical, dimension(:), allocatable :: wifseg_solid,wifseg_fluid

! valence test
  integer :: idest,i
  character(len=4) :: appiproc
  double precision, dimension(:), allocatable :: uglob2_solid
  double precision, dimension(:), allocatable :: val_solid(:,:,:)
  integer :: valnum_cent_solid(6),totvalnum_cent_solid
  integer :: valnum_semi_solid(6),totvalnum_semi_solid

  double precision, dimension(:), allocatable :: uglob2_fluid
  double precision, dimension(:), allocatable :: val_fluid(:,:,:)
  integer :: valnum_fluid(6),totvalnum_fluid

  nelmax_solid = maxval(nel_solid);   nelmax_fluid = maxval(nel_fluid)

  allocate (igloc_solid(nelmax_solid*(npol+1)**2,0:nproc-1))
  allocate (igloc_fluid(nelmax_fluid*(npol+1)**2,0:nproc-1))

  allocate (nglobp_solid(0:nproc-1),nglobp_fluid(0:nproc-1))
  nglobp_solid(:) = 0;  nglobp_fluid(:) = 0

!-----------------------
  do iproc = 0, nproc-1
!-----------------------

   nelp_solid = nel_solid(iproc) 
   nelp_fluid = nel_fluid(iproc) 
   npointotp_solid = nelp_solid*(npol+1)**2
   npointotp_fluid = nelp_fluid*(npol+1)**2
   if (have_solid) &
        allocate(wsgll_solid(npointotp_solid),wzgll_solid(npointotp_solid))

   if (have_fluid) &
        allocate(wsgll_fluid(npointotp_fluid),wzgll_fluid(npointotp_fluid))

! Solid
   if (have_solid) then
   do iel = 1, nelp_solid
      ielg = procel_solid(iel,iproc)
      do jpol = 0, npol
         do ipol = 0, npol
            ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
            wsgll_solid(ipt) = sgll(ipol,jpol,ielg)
            wzgll_solid(ipt) = zgll(ipol,jpol,ielg)
         end do
      end do
   end do

   allocate(wigloc_solid(npointotp_solid))
   allocate(wifseg_solid(npointotp_solid))
   allocate(wloc_solid(npointotp_solid))
   call get_global(nelp_solid,wsgll_solid,wzgll_solid,wigloc_solid, &
        wloc_solid,wifseg_solid,wnglob_solid,npointotp_solid,NGLLCUBE,NDIM)
   do ipt = 1, npointotp_solid
      igloc_solid(ipt,iproc) = wigloc_solid(ipt)
   end do
   nglobp_solid(iproc) = wnglob_solid

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SOLID valence: test slocal numbering
  allocate(uglob2_solid(nglobp_solid(iproc))) 
  allocate(val_solid(0:npol,0:npol,nel_solid(iproc)))

! valence test, equivalent to how assembly is used in the solver
! use script plot_proc_valence.csh to generate GMT valence grids for each proc
! glocally and zoomed into r<0.2 ( denoted as *_central )
 val_solid(:,:,:) = 1.0d0

 valnum_cent_solid=0
 valnum_semi_solid=0

 uglob2_solid(:) = 0.d0
 do iel = 1, nel_solid(iproc)
    do ipol = 0, npol
       do jpol = 0, npol
          ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
          idest = wigloc_solid(ipt)
          uglob2_solid(idest) = uglob2_solid(idest) + val_solid(ipol,jpol,iel)
       end do
    end do
 end do
 !QT
 do iel = 1, nel_solid(iproc)
    do ipol = 0, npol
       do jpol = 0, npol
          ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
          idest = wigloc_solid(ipt)
          val_solid(ipol,jpol,iel) = uglob2_solid(idest)
          
! Check valence/global number inside central region
          if (eltypeg(procel_solid(iel,iproc))=='linear') then 
             do i=1,6 !possible valences
                if (val_solid(ipol,jpol,iel)==i) &
                     valnum_cent_solid(i)=valnum_cent_solid(i)+1
             enddo
          endif
          
          if (eltypeg(procel_solid(iel,iproc))=='semino' .or. &
               eltypeg(procel_solid(iel,iproc))=='semiso') then 
             do i=1,6 !possible valences
                if (val_solid(ipol,jpol,iel)==i) &
                     valnum_semi_solid(i)=valnum_semi_solid(i)+1
             enddo
          endif
          
       end do
    end do
 end do
 
 deallocate(uglob2_solid,val_solid)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   deallocate(wloc_solid)
   deallocate(wifseg_solid)
   deallocate(wigloc_solid)
   deallocate(wzgll_solid,wsgll_solid)
   
   endif ! have_solid

! Fluid
   if (have_fluid) then
      do iel = 1, nelp_fluid
         ielg = procel_fluid(iel,iproc)
         do jpol = 0, npol
            do ipol = 0, npol
               ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1 
               wsgll_fluid(ipt) = sgll(ipol,jpol,ielg)
               wzgll_fluid(ipt) = zgll(ipol,jpol,ielg)
            end do
         end do
      end do

      allocate(wigloc_fluid(npointotp_fluid))
      allocate(wifseg_fluid(npointotp_fluid))
      allocate(wloc_fluid(npointotp_fluid))
      call get_global(nelp_fluid,wsgll_fluid,wzgll_fluid,wigloc_fluid, &
           wloc_fluid,wifseg_fluid,wnglob_fluid,npointotp_fluid,NGLLCUBE,NDIM)
      do ipt = 1, npointotp_fluid
         igloc_fluid(ipt,iproc) = wigloc_fluid(ipt)
      end do
      nglobp_fluid(iproc) = wnglob_fluid
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! FLUID valence: test flocal numbering
  allocate(uglob2_fluid(nglobp_fluid(iproc))) 
  allocate(val_fluid(0:npol,0:npol,nel_fluid(iproc)))

! valence test, equivalent to how assembly is used in the solver
! use script plot_proc_valence.csh to generate GMT valence grids for each proc
! glocally and zoomed into r<0.2 ( denoted as *_central )
 val_fluid(:,:,:) = 1.0d0

 valnum_fluid=0

 uglob2_fluid(:) = 0.d0
 do iel = 1, nel_fluid(iproc)      
  do ipol = 0, npol
   do jpol = 0, npol
    ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
    idest = wigloc_fluid(ipt)
    uglob2_fluid(idest) = uglob2_fluid(idest) + val_fluid(ipol,jpol,iel)
   end do
  end do
 end do
!QT
 do iel = 1, nel_fluid(iproc)
  do ipol = 0, npol
   do jpol = 0, npol
    ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
    idest = wigloc_fluid(ipt)
    val_fluid(ipol,jpol,iel) = uglob2_fluid(idest)

! Check valence/global number inside central region
        do i=1,6 !possible valences
           if (val_fluid(ipol,jpol,iel)==i) &
                valnum_fluid(i)=valnum_fluid(i)+1
        enddo

   end do
  end do

 end do

 deallocate(uglob2_fluid,val_fluid)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     deallocate(wloc_fluid)
     deallocate(wifseg_fluid)
     deallocate(wigloc_fluid)
     deallocate(wzgll_fluid,wsgll_fluid)     

   endif ! have_fluid

   if (have_solid) then
      totvalnum_cent_solid = 0
      totvalnum_semi_solid = 0
      do i=1,6
         totvalnum_cent_solid = totvalnum_cent_solid + valnum_cent_solid(i)/i
         totvalnum_semi_solid = totvalnum_semi_solid + valnum_semi_solid(i)/i
      enddo

      if (dump_mesh_info_screen) then 
         write(6,*)
         write(6,*)' SOLID VALENCE & GLOBAL NUMBERING:'
         write(6,*)iproc,'slocal number      :',nglobp_solid(iproc)
         write(6,*)iproc,'scentral region only:',totvalnum_cent_solid
         write(6,*)iproc,'ssemicurved els only:',totvalnum_semi_solid
         write(6,*)iproc,'severywhere else    :',nglobp_solid(iproc)-&
              totvalnum_cent_solid-totvalnum_semi_solid
      end if      
   endif ! have_solid

   if (have_fluid) then
 totvalnum_fluid = 0
 do i=1,6
    totvalnum_fluid = totvalnum_fluid + valnum_fluid(i)/i

 enddo

 if (dump_mesh_info_screen) then 
  write(6,*)' FLUID VALENCE & GLOBAL NUMBERING:'
  write(6,*)iproc,'flocal number      :',nglobp_fluid(iproc)
  write(6,*)iproc,'fluid region only  :',totvalnum_fluid
  write(6,*)iproc,'feverywhere else   :',nglobp_fluid(iproc)-&
                                          totvalnum_fluid
 end if
endif ! have_fluid

!-----------------------
  end do ! iproc
!-----------------------

  if (dump_mesh_info_screen) then 
! af 
! here are used some arrays such as nglobp that are not needed
! nglobp will not be defined in subsequent versions of the mesher
! this whole printout is bound to be erased
  write(6,*)
  write(6,*)'Amount of unique grid points in each processor:'
  do iproc = 0, nproc-1
   write(6,21)iproc,nglobp_solid(iproc),nglobp_fluid(iproc),nglobp(iproc),&
              nglobp_solid(iproc)+nglobp_fluid(iproc)-nglobp(iproc)
  end do

  write(6,*)
  write(6,*) 'SUM # grid points Solid:       ',SUM(nglobp_solid(0:nproc-1))
  write(6,*) 'SUM # grid points Fluid:       ',SUM(nglobp_fluid(0:nproc-1))
  write(6,*) 'SUM # grid points Total:       ',SUM(nglobp(0:nproc-1))
  write(6,*) '    # total grid points:       ',nglobglob
  write(6,*) 'SUM # grid points S/F Boundary:',SUM(nglobp_fluid(0:nproc-1))+ &
                            SUM(nglobp_solid(0:nproc-1))-SUM(nglobp(0:nproc-1))
  write(6,*) 'Predicted S/F Boundary points :', &
                              2*(npol+1)+(sum(nbelem)/2-nbcnd)*(npol)

 end if

! This test is not valid for nproc>1 since the proc boundaries are singular 
! global points for each proc and hence one would need to subtract the 
! sum of the proc boundary points here. 
! Could be done ;)
!
!  if (have_fluid .and. 2*(npol+1)+(sum(nbelem)/2-nbcnd)*(npol) /= &
!          SUM(nglobp_fluid(0:nproc-1))+ &
!          SUM(nglobp_solid(0:nproc-1))-SUM(nglobp(0:nproc-1)) ) then
!     write(6,*)'...something wrong with global number of S/F boundary points..'
!     call flush(6)
!     stop
!  endif

  if (dump_mesh_info_screen) write(6,*)

21 format('Proc',i3, ' has',i9, ' solid ',i9,' fluid',i9,&
          ' total and',i6,' S/F boundary pts')
! END ALEX
end subroutine define_sflocal_numbering
!------------------------------------------------------------------------

!dk define_global2glocal-------------------------------------------------
subroutine define_global2glocal
!
! Define bookkeeping array that for a given global number, returns 
! the glocal number. Used later on to partition global index into message bins,
! and deallocated there, in partition_global_index
! Both this and partition_global_index might well be redundant after all....
!
integer :: ipt,iel,ipol,jpol,iproc,ielg,iptg,igg,igp
allocate(glob2gloc(nglobglob,0:nproc-1)) ; glob2gloc(:,:) = 0

  do iproc = 0, nproc-1
     do iel = 1, nel(iproc)
        ielg = procel(iel,iproc)
        do jpol = 0, npol
           do ipol = 0, npol
              iptg = (ielg-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1 
              igg = iglob(iptg)      ! global #
              igp = igloc(ipt,iproc) ! glocal # 
              glob2gloc(igg,iproc) = igp 
           end do
        end do
     end do
  end do

  deallocate(igloc)

end subroutine define_global2glocal
!------------------------------------------------------------------------

!dk define_sflobal2sflocal-------------------------------------------------
subroutine define_sflobal2sflocal
!
! Define bookkeeping array that for a given sflobal number, returns 
! the sflocal number. Used later on to partition global index into message 
! bins for the solid/fluid subdomains, and should be deallocated there, 
! in partition_sflobal_index
!
integer :: ipt,iel,ipol,jpol,iproc,ielg,iptg,igg,igp

  if (dump_mesh_info_screen) then 
   write(6,*)'Info on global numbering array sizes/values:'
   write(6,*)'MAX INV_SOL/FLU     :',maxval(inv_ielem_solid),&
                                    maxval(inv_ielem_fluid)
!   write(6,*)'MAX/SIZE IGLOB      :',maxval(iglob),size(iglob)
    if (have_solid) write(6,*)'MAX/SIZE IGLOB_SOLID:',maxval(iglob_solid),size(iglob_solid)
  if (have_fluid) & 
   write(6,*)'MAX/SIZ IGLOB_FLUID :',maxval(iglob_fluid), size(iglob_fluid)

    if (have_solid) write(6,*)'MAX/SIZE SOL IGLOC  :',maxval(igloc_solid),size(igloc_solid)
   if (have_fluid) & 
   write(6,*)'MAX/SIZE FLU IGLOC  :',maxval(igloc_fluid),size(igloc_fluid)

!  write(6,*)'MAX LOC/SOL/FLU    :',maxval(loc),maxval(loc_solid), &
!                                   maxval(loc_fluid)
!  write(6,*)'SIZ LOC/SOL/FLU    :',size(loc),size(loc_solid), &
!                                   size(loc_fluid)

   write(6,*)'MAX PROCEL/SOL/FLU  :',maxval(procel),maxval(procel_solid),&
                                    maxval(procel_fluid)
   write(6,*)'SIZ PROCEL/SOL/FLU  :',size(procel),size(procel_solid),&
                                    size(procel_fluid)
   write(6,*)
  call flush(6)
  end if 

! Solid
  allocate(slob2sloc(nglobslob,0:nproc-1)) ; slob2sloc(:,:) = 0
  do iproc = 0, nproc-1
   do iel = 1, nel_solid(iproc)
    ielg = procel_solid(iel,iproc) 
    do jpol = 0, npol
     do ipol = 0, npol
      iptg = (inv_ielem_solid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
      ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
      igg = iglob_solid(iptg)      ! slobal #
      igp = igloc_solid(ipt,iproc) ! slocal #
      slob2sloc(igg,iproc) = igp 
     end do
    end do
   end do
  end do

! Fluid
  if (have_fluid) then
  allocate(flob2floc(nglobflob,0:nproc-1)) ; flob2floc(:,:) = 0
  do iproc = 0, nproc-1
   do iel = 1, nel_fluid(iproc)
    ielg = procel_fluid(iel,iproc)
    do jpol = 0, npol
     do ipol = 0, npol
      iptg = (inv_ielem_fluid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
      ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
      igg = iglob_fluid(iptg)
      igp = igloc_fluid(ipt,iproc) 
      flob2floc(igg,iproc) = igp 
     end do
    end do
   end do
  end do
  endif

end subroutine define_sflobal2sflocal
!------------------------------------------------------------------------

!dk define_search_global_index--------------------------------------------
subroutine define_search_global_index
!
! Define a global index that given GLL coords for a proc returns global number.
! Also defined are
! nprocb: Given global number returns number of procs it belongs to (valence)
! lprocb: Given global number and index of valence returns processor ID
!
integer :: iproc,ielg,iel,ipol,jpol,ipt,nelmax
integer :: il,nprocbmax
integer, allocatable :: nbelong2(:,:)

  write(6,*) ' NEW GLOBAL INDEX SEARCH... ', nglobglob ; call flush(6)

  allocate(nbelong(nglobglob))
  allocate(nbelong2(0:nproc-1,nglobglob))
  allocate(nprocb(nglobglob))

  nbelong(:)=0; nbelong2(:,:) = 0; nprocbmax = 1
  nelmax = maxval(nel)
  allocate(global_index(0:npol,0:npol,nelmax,0:nproc-1))

  do iproc = 0, nproc -1
     do iel = 1, nel(iproc)
        ielg = procel(iel,iproc) ! global element number
        do jpol = 0, npol
           do ipol = 0, npol
              ipt = (ielg-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
              global_index(ipol,jpol,iel,iproc) = iglob(ipt)
              
              nbelong(iglob(ipt)) = nbelong(iglob(ipt)) + 1 ! global valence
              nbelong2(iproc,iglob(ipt)) = 1 ! this glob number is in iproc
              if (nprocbmax < sum(nbelong2(:,iglob(ipt)))) then
                 nprocbmax = sum(nbelong2(:,iglob(ipt)))
              endif
           end do
        end do
     end do
  end do !iproc
  
  write(6,*)'nprocbmax=',nprocbmax
  allocate(lprocb(nprocbmax,nglobglob))

  do ipt=1,nglobglob
     il=0
     do iproc=0,nproc-1
        if (nbelong2(iproc,ipt)==1) then
           il = il + 1
           lprocb(il,ipt)=iproc
           nprocb(ipt)=il
        endif
     enddo
  enddo

  deallocate(nbelong2)

!  do iproc=0,nproc-1
!  do ipt=1,nglobglob
!     write(88888,*)lprocb(iproc+1,ipt)
!  enddo
!  enddo

  end subroutine define_search_global_index
!------------------------------------------------------------------------

!dk define_search_sflobal_index----------------------------------------------
subroutine define_search_sflobal_index 
!
! Define a sflobal index: given GLL coords and proc returns sflobal number.
! Also defined are
! nprocb_solid: Given slobal number returns # of procs it belongs to (valence)
! lprocb_solid: Given slobal number and index of valence returns processor ID
! .... and correspondingly for fluid arrays.
! These are crucial for the partitioning, i.e. the heart of message passing.
!
integer :: iproc,ielg,iel,ipol,jpol,ipt,nelmax_solid,nelmax_fluid
integer :: il,nprocbmax_solid,nprocbmax_fluid
integer, allocatable :: nbelong2_solid(:,:),nbelong2_fluid(:,:)

  nelmax_solid = maxval(nel_solid); nelmax_fluid = maxval(nel_fluid)
  allocate(slobal_index(0:npol,0:npol,nelmax_solid,0:nproc-1))
  allocate(flobal_index(0:npol,0:npol,nelmax_fluid,0:nproc-1))
  allocate(nprocb_solid(nglobslob),nprocb_fluid(nglobflob))

  allocate(nbelong_solid(nglobslob))
  allocate(nbelong2_solid(0:nproc-1,nglobslob))
  nbelong_solid(:)=0; nbelong2_solid(:,:) = 0; nprocbmax_solid = 1


if (have_fluid) then
  allocate(nbelong_fluid(nglobflob))
  allocate(nbelong2_fluid(0:nproc-1,nglobflob))
  nbelong_fluid(:)=0; nbelong2_fluid(:,:) = 0; nprocbmax_fluid = 1
endif

  do iproc = 0, nproc -1 

! solid
   do iel = 1, nel_solid(iproc)
    ielg = procel_solid(iel,iproc) ! global element number
    do jpol = 0, npol
     do ipol = 0, npol
      ipt = (inv_ielem_solid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
      slobal_index(ipol,jpol,iel,iproc) = iglob_solid(ipt)

      nbelong_solid(iglob_solid(ipt)) = nbelong_solid(iglob_solid(ipt)) + 1 
      nbelong2_solid(iproc,iglob_solid(ipt)) = 1 ! this glob number is in iproc
      if (nprocbmax_solid < sum(nbelong2_solid(:,iglob_solid(ipt)))) then
        nprocbmax_solid = sum(nbelong2_solid(:,iglob_solid(ipt)))
      endif  
     
     end do
    end do
   end do

! fluid
if (have_fluid) then
   do iel = 1, nel_fluid(iproc)
    ielg = procel_fluid(iel,iproc) ! global element number
    do jpol = 0, npol
     do ipol = 0, npol
      ipt = (inv_ielem_fluid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
      flobal_index(ipol,jpol,iel,iproc) = iglob_fluid(ipt)

      nbelong_fluid(iglob_fluid(ipt)) = nbelong_fluid(iglob_fluid(ipt)) + 1 
      nbelong2_fluid(iproc,iglob_fluid(ipt)) = 1 ! this glob number is in iproc
      if (nprocbmax_fluid < sum(nbelong2_fluid(:,iglob_fluid(ipt)))) then
        nprocbmax_fluid = sum(nbelong2_fluid(:,iglob_fluid(ipt)))
      endif  

     end do
    end do
   end do 

end if


  end do !iproc

  allocate(lprocb_solid(nprocbmax_solid,nglobslob))

  do ipt=1,nglobslob
     il=0
     do iproc=0,nproc-1
        if (nbelong2_solid(iproc,ipt)==1) then
           il=il+1
           lprocb_solid(il,ipt)=iproc
           nprocb_solid(ipt) = il
        endif
     enddo
  enddo

  deallocate(nbelong2_solid)

!  do iproc=0,nproc-1
!  do ipt=1,nglobslob
!    write(88888,*)lprocb_solid(iproc+1,ipt)
!  enddo
!  enddo

if (have_fluid) then
allocate(lprocb_fluid(nprocbmax_fluid,nglobflob))

  do ipt=1,nglobflob
     il=0
     do iproc=0,nproc-1
        if (nbelong2_fluid(iproc,ipt)==1) then
           il = il + 1
           lprocb_fluid(il,ipt)=iproc
           nprocb_fluid(ipt) = il
        endif
     enddo
  enddo

  deallocate(nbelong2_fluid)

endif

!  do iproc=0,nproc-1
!  do ipt=1,nglobslob
!    write(88888,*)lprocb_fluid(iproc+1,ipt)
!  enddo
!  enddo

  end subroutine define_search_sflobal_index
!------------------------------------------------------------------------

!dk partition_global_index-----------------------------------------------
subroutine partition_global_index
!
! Partitioning the index defined in define_search_global_index, i.e. 
! creating bins for the message passing across processors. These inherently 
! depend on the glocal number as the message passing is done during the 
! assembly, i.e. glocal stage in the solver. 
!
! SOME TERM EXPLANATIONS:
! nglobglob: global number of points
! sizebin: size (bin (iproc) ) 
! binp: array of global indices (i.e. glob point #) attributed to a given iproc
!
! This will be redundant as we only need bins and message passing for solid 
! and fluid subdomains! TO BE DISCARDED!
!
integer :: ipt,iproct,sizebinmax,ibp,ibel,ig,ip
integer :: ipdes,ipsrc,imsg
integer, dimension(:), allocatable :: ibin 
integer, dimension(:,:), allocatable :: sizemsg
integer, dimension(:,:), allocatable :: index_msg
integer, dimension(:,:), allocatable :: binp
integer, dimension(:), allocatable :: sizebin
integer, dimension(:,:,:), allocatable :: global_index_msg
integer :: sizemsgmax,sizerecvpmax,sizesendpmax

  write(6,*) ' Creating bins ' 
  allocate (sizebin(0:nproc-1))
  sizebin(0:nproc-1) = 0 
  do ipt = 1, nglobglob
     iproct = lprocb(1,ipt)! completely ad hoc choice
     sizebin(iproct)=sizebin(iproct) + 1
  end do

  write(6,*) ' BIN SIZES ' 
  do iproct = 0, nproc - 1
     write(6,*) iproct, sizebin(iproct)
  end do

  write(6,*) SUM(sizebin),nglobglob
  write(6,*) 'end check'
  sizebinmax = maxval (sizebin(:)) 

  write(6,*) ' sizebinmax = ' , sizebinmax

  allocate(binp(sizebinmax,0:nproc-1)) 
  binp(1:sizebinmax,0:nproc-1) = 0 

  allocate(ibin(0:nproc-1))
  ibin(0:nproc-1) = 0  

  do ipt = 1, nglobglob
     iproct =  lprocb(1,ipt)! same choice 
     ibin(iproct) = ibin(iproct) + 1
     ibp = ibin(iproct)
     binp(ibp,iproct) = ipt
  end do
  deallocate(ibin)

  ! do iproct = 0, nproc-1
  !  write(6,20) (binp(ibp,iproct),ibp=1,sizebin(iproct))
  ! end do
  allocate(sizemsg(0:nproc-1,0:nproc-1))
  sizemsg(0:nproc-1,0:nproc-1) = 0 

  write(6,*) ' DEPENDENCIES ' 
  do iproct = 0, nproc - 1
     !  write(6,*) ' PROC IS ' , iproct
     do ipt = 1, sizebin(iproct)
        ig = binp(ipt,iproct)
        !   write(6,*) nprocb(ig)
        do ibel = 2, nprocb(ig)
           !    write(6,*) 'iglob =',ig, 'nprocb = ',nprocb(ig), &
           ! (lprocb(ibel,ig),ibel=2,nprocb(ig))
           ipdes = lprocb(1,ig) ! same choice 
           ipsrc = lprocb(ibel,ig)
           sizemsg(ipsrc,ipdes) = sizemsg(ipsrc,ipdes) + 1 
        end do
     end do
  end do

  do iproct = 0,  nproc-1
     write(6,21) (sizemsg(iproct,ipdes),ipdes=0,nproc-1)
  end do
  write(6,*) SUM(SUM(sizemsg,DIM=1))

  sizemsgmax = maxval(maxval(sizemsg,DIM=1))
  write(6,*) 'size msg max is ' , sizemsgmax

  allocate(index_msg(0:nproc-1,0:nproc-1))
  index_msg(0:nproc-1,0:nproc-1) = 0

  allocate(global_index_msg(sizemsgmax,0:nproc-1,0:nproc-1))
  do iproct = 0, nproc-1
     do ipt = 1, sizebin(iproct)
        ig = binp(ipt,iproct)
        do ibel = 2, nprocb(ig)
           ipdes = lprocb(1,ig) ! same choice again
           ipsrc = lprocb(ibel,ig)
           index_msg(ipsrc,ipdes) = index_msg(ipsrc,ipdes) + 1
           imsg = index_msg(ipsrc,ipdes)
           global_index_msg(imsg,ipsrc,ipdes) = ig
        end do
     end do
  end do

! TNM: Why = 5 ??? changed this...
! ipsrc = 5
  ipsrc = nproc-1
  do ipdes = 0, nproc-1 
     write(6,20) ipdes,sizemsg(ipsrc,ipdes),&
          (global_index_msg(imsg,ipsrc,ipdes),imsg=1,sizemsg(ipsrc,ipdes))
  end do
  deallocate(index_msg)

! How many messages am I sending / receiving ? 
  allocate(sizerecvp(0:nproc-1),sizesendp(0:nproc-1))
  sizerecvp(0:nproc-1) = 0
  sizesendp(0:nproc-1) = 0 
  do iproct = 0, nproc -1
     do ipsrc = 0, nproc-1
        if ( sizemsg(ipsrc,iproct) > 0 ) sizerecvp(iproct) = sizerecvp(iproct) + 1 
     end do
     do ipdes = 0, nproc-1 
        if ( sizemsg(iproct,ipdes) > 0 ) sizesendp(iproct) = sizesendp(iproct) + 1 
     end do
  end do

  do iproct = 0, nproc -1
     write(6,*) ' processor ', iproct, 'will receive stuff from '  , &
          sizerecvp(iproct), ' procs'
     write(6,*) ' processor ', iproct, 'will send    stuff to   '  , &
          sizesendp(iproct), ' procs'
  end do
  sizerecvpmax = maxval(sizerecvp(:))
  sizesendpmax = maxval(sizesendp(:))

! To which processor ? 
  allocate(listrecvp(sizerecvpmax,0:nproc-1))
  allocate(listsendp(sizesendpmax,0:nproc-1))
  listrecvp(1:sizerecvpmax,0:nproc-1) = -1 
  listsendp(1:sizesendpmax,0:nproc-1) = -1

  do iproct = 0, nproc -1 
     ip = 0
     do ipsrc = 0, nproc -1
        if ( sizemsg(ipsrc,iproct) > 0 ) then
           ip = ip + 1
           listrecvp(ip,iproct) = ipsrc
        end if
     end do
     ip = 0
     do ipdes = 0, nproc -1
        if ( sizemsg(iproct,ipdes) > 0 ) then
           ip = ip + 1
           listsendp(ip,iproct) = ipdes
        end if
     end do
  end do

  do iproct = 0, nproc-1
     write(6,22) iproct, sizerecvp(iproct), (listrecvp(ip,iproct),ip=1,sizerecvp(iproct))
     write(6,23) iproct, sizesendp(iproct), (listsendp(ip,iproct),ip=1,sizesendp(iproct))
  end do
22 format('Proc', i2, ' will receive ', i2, ' messages from procs ', 20(i3,1x))
23 format('Proc', i2, ' will send    ', i2, ' messages to   procs ', 20(i3,1x))

! What size ? 
  allocate(sizemsgrecvp(sizerecvpmax,0:nproc-1))
  allocate(sizemsgsendp(sizesendpmax,0:nproc-1))
  do iproct = 0, nproc-1
     do ip = 1, sizerecvp(iproct)
        ipsrc = listrecvp(ip,iproct)
        sizemsgrecvp(ip,iproct) = sizemsg(ipsrc,iproct)
     end do
     do ip = 1, sizesendp(iproct) 
        ipdes = listsendp(ip,iproct)
        sizemsgsendp(ip,iproct) = sizemsg(iproct,ipdes)
     end do
  end do
  do iproct = 0, nproc -1
     write(6,*) 'PROC ', iproct, ' RECEIVING INFO '
     do ip = 1, sizerecvp(iproct)
        write(6,*) listrecvp(ip,iproct),sizemsgrecvp(ip,iproct)
     end do
     write(6,*) ' SENDING INFO '
     do ip = 1, sizesendp(iproct)
        write(6,*) listsendp(ip,iproct),sizemsgsendp(ip,iproct)
     end do
  end do
! NOW CREATE GLOCAL INDEX FOR MESSAGES 
! Which glocal indices?
!!  
! ALEX TUESDAY Pm : must change the second indices everywhere !!!!
! TO sizemsgrecvpmax AND sizemsgsenpmax RESPECTIVELY
  allocate(glocal_index_msg_recvp(sizemsgmax,sizerecvpmax,0:nproc-1))
  allocate(glocal_index_msg_sendp(sizemsgmax,sizesendpmax,0:nproc-1))
! END ALEX TUESDAY Pm
!!
  glocal_index_msg_recvp(1:sizemsgmax,1:sizerecvpmax,0:nproc-1) = 0
  glocal_index_msg_sendp(1:sizemsgmax,1:sizesendpmax,0:nproc-1) = 0
  do iproct = 0, nproc-1
!! do iproct=0,0
!     write(6,*)
!     write(6,*) '--------------- PROC ',iproct,'------------------------'
!     write(6,*) iproct,'  SHALL RECEIVE: '
!!   write(6,*) iproct,sizerecvp(iproct)
     do ip = 1, sizerecvp(iproct)
        ipsrc = listrecvp(ip,iproct)
!        write(6,*) ipsrc, sizemsgrecvp(ip,iproct)
        do ipt = 1, sizemsgrecvp(ip,iproct)  
           ig = global_index_msg(ipt,ipsrc,iproct)
           glocal_index_msg_recvp(ipt,ip,iproct) = glob2gloc(ig,iproct)
           !!     write(6,*) ipsrc,ipt,ig,glocal_index_msg_recvp(ipt,ip,iproct)
        end do
     end do
!     write(6,*) iproct,'  SHALL SEND   : '
     do ip = 1, sizesendp(iproct)
        ipdes = listsendp(ip,iproct)
!        write(6,*) ipdes, sizemsgsendp(ip,iproct)
        do ipt = 1, sizemsgsendp(ip,iproct)
           ig = global_index_msg(ipt,iproct,ipdes)
           glocal_index_msg_sendp(ipt,ip,iproct) = glob2gloc(ig,iproct)
!           write(6,*) ipdes,ipt,ig,glocal_index_msg_sendp(ipt,ip,iproct)
        end do
     end do
  end do
20 format(200(i5,1x))
21 format(80(i3,1x))

!  deallocate(glob2gloc)
  deallocate(binp,sizebin,global_index_msg)

  end subroutine partition_global_index
!------------------------------------------------------------------------

!dk partition_sflobal_index-----------------------------------------------
subroutine partition_sflobal_index
!
! Partitioning the index defined in define_search_sflobal_index, i.e. 
! creating bins for the message passing across processors. These inherently 
! depend on the sflocal number as the message passing is done during the 
! assembly, i.e. sflocal stage in the solver. 
!
! SOME TERM EXPLANATIONS (All have _solid appended, left out for brevity):
!
! nglobglob: global number of points
! nprocb: Given global numnber, lists the number of procs it belongs to
! lprocb: given global number, indices 1,..,nprocb list proc it belongs to
! sizebin: given a proc, lists glocal number
! binp: given glocal number and processor, lists global number
! sizemsg: for each proc-proc pair, lists number of global points shared
! global_index_msg: for each proc-proc pair and shared points 1,...,sizemsgmax
!                   lists global number
! THE FOLLOWING ENTER THE DATABASE FOR THE SOLVER:
! sizerecvp,sizesendp: given a proc, lists # other procs to communicate with
! listrecvp,listsendp: given a proc and all its shared points, lists proc ID 
!                      to send/receive that point with
! sizemsgrecvp,sizemsgsendp: given proc and each communicating proc, list size 
!                            of respective messages sent/received
! glocal_index_msg_recvp,glocal_index_msg_sendp: same as glocal_ind..for glocal
!
! DATABASE: sizerecvp_solid,listrecvp_solid
! DATABASE: sizemsgrecvp_solid
! DATABASE: glocal_index_msg_recvp_solid
! DATABASE: sizesendp_solid,listsendp_solid
! DATABASE: sizemsgsendpg_solid
! DATABASE: glocal_index_msg_sendp_solid
!
! DATABASE: sizerecvp_fluid,listrecvp_fluid
! DATABASE: sizemsgrecvp_fluid
! DATABASE: glocal_index_msg_recvp_fluid
! DATABASE: sizesendp_fluid,listsendp_fluid
! DATABASE: sizemsgsendpg_fluid
! DATABASE: glocal_index_msg_sendp_fluid

integer :: ipt,iproct,ibp,ibel,ig,ip
integer :: ipdes,ipsrc,imsg

integer, dimension(:), allocatable :: ibin_solid
integer, dimension(:,:), allocatable :: sizemsg_solid
integer, dimension(:,:), allocatable :: index_msg_solid
integer, dimension(:,:,:), allocatable :: global_index_msg_solid
integer, dimension(:,:), allocatable :: binp_solid
integer, dimension(:), allocatable :: sizebin_solid
integer :: sizebinmax_solid,sizemsgmax_solid
integer :: sizerecvpmax_solid,sizesendpmax_solid

integer, dimension(:), allocatable :: ibin_fluid
integer, dimension(:,:), allocatable :: sizemsg_fluid
integer, dimension(:,:), allocatable :: index_msg_fluid
integer, dimension(:,:,:), allocatable :: global_index_msg_fluid
integer, dimension(:,:), allocatable :: binp_fluid
integer, dimension(:), allocatable :: sizebin_fluid
integer :: sizebinmax_fluid,sizemsgmax_fluid
integer :: sizerecvpmax_fluid,sizesendpmax_fluid

! valence test
  integer :: idest,iel,ipol,jpol
  character(len=4) :: appiproc
  integer, dimension(:), allocatable :: uglob2_solid
  integer, dimension(:), allocatable :: val_solid(:,:,:)


  if (dump_mesh_info_screen) then 
   write(6,*)
   write(6,*)'****************************************************************'
   write(6,*)'******************* CREATING MESSAGING ARRAYS ******************'
   write(6,*)'****************************************************************'
  end if
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SOLID SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS 

  if (dump_mesh_info_screen) write(6,*) 'Creating solid bins ' 
  allocate (sizebin_solid(0:nproc-1))
  sizebin_solid(0:nproc-1) = 0
  do ipt = 1, nglobslob
     iproct = lprocb_solid(1,ipt)! completely ad hoc choice
     sizebin_solid(iproct)=sizebin_solid(iproct) + 1
  end do

  if (dump_mesh_info_screen) then 
   do iproct = 0, nproc - 1
      write(6,*) 'proc:',iproct,'size solid bin:',sizebin_solid(iproct)
   end do
  end if
 
  if (dump_mesh_info_screen) write(6,*) 'sum sizebin,slobal:',SUM(sizebin_solid),nglobslob
  if (SUM(sizebin_solid)/=nglobslob) then 
     write(6,*)'PROBLEM: sum of solid bins not equal to slobal num!'
     stop
  endif
  sizebinmax_solid = maxval (sizebin_solid(:)) 

  if (dump_mesh_info_screen) write(6,*) ' sizebinmax_solid = ' , sizebinmax_solid

  allocate(binp_solid(sizebinmax_solid,0:nproc-1)) 
  binp_solid(1:sizebinmax_solid,0:nproc-1) = 0 

  allocate(ibin_solid(0:nproc-1))
  ibin_solid(0:nproc-1) = 0  

  do ipt = 1, nglobslob
     iproct =  lprocb_solid(1,ipt)! same choice 
     ibin_solid(iproct) = ibin_solid(iproct) + 1
     ibp = ibin_solid(iproct)
     binp_solid(ibp,iproct) = ipt
  end do
  deallocate(ibin_solid)

  ! do iproct = 0, nproc-1
  !  write(6,20) (binp(ibp,iproct),ibp=1,sizebin(iproct))
  ! end do
  allocate(sizemsg_solid(0:nproc-1,0:nproc-1))
  sizemsg_solid(0:nproc-1,0:nproc-1) = 0 

!!$!  write(6,*) ' SOLID DEPENDENCIES ' 
!!$  do iproct = 0, nproc - 1
!!$!  write(6,*) ' PROC IS ' , iproct
!!$     do ipt = 1, sizebin_solid(iproct)
!!$        ig = binp_solid(ipt,iproct)
!!$        !   write(6,*) nprocb(ig)
!!$        do ibel = 2, nprocb_solid(ig)
!!$           !    write(6,*) 'iglob =',ig, 'nprocb = ',nprocb(ig), &
!!$           ! (lprocb(ibel,ig),ibel=2,nprocb(ig))       
!!$           ipdes = lprocb_solid(1,ig) ! same choice 
!!$           ipsrc = lprocb_solid(ibel,ig)
!!$           sizemsg_solid(ipsrc,ipdes) = sizemsg_solid(ipsrc,ipdes) + 1 
!!$        end do
!!$     end do
!!$  end do

! TNM NEW ATTEMPT.....
!  write(6,*) ' SOLID DEPENDENCIES ' 
  do iproct = 0, nproc - 1
!  write(6,*) ' PROC IS ' , iproct
     do ipt = 1, sizebin_solid(iproct)
        ig = binp_solid(ipt,iproct)
        !   write(6,*) nprocb(ig)
        do ibel = 2, nprocb_solid(ig)
           !    write(6,*) 'iglob =',ig, 'nprocb = ',nprocb(ig), &
           ! (lprocb(ibel,ig),ibel=2,nprocb(ig))       
           ipdes = lprocb_solid(1,ig) ! same choice           
           ipsrc = lprocb_solid(ibel,ig)
           sizemsg_solid(ipsrc,ipdes) = sizemsg_solid(ipsrc,ipdes) + 1 
        end do
     end do
  end do

! TNM MAY07: MAKE THIS SYMMETRIC; send/receive need to happen in both ways...
! NOTE THE ORDER OF THE LOOP! Other way around is ZERO throughout............
!  do ipdes = 0,nproc-1
!     do ipsrc = 0,nproc-1
!        sizemsg_solid(ipdes,ipsrc) = sizemsg_solid(ipsrc,ipdes)
!     enddo
!  enddo

  if (dump_mesh_info_screen) then 
   write(6,*)'Size of solid messages for each proc-proc pair:'
   write(6,*)'---> destination proc, down: my proc'
   do iproct = 0,  nproc-1
      write(6,21) (sizemsg_solid(iproct,ipdes),ipdes=0,nproc-1)
   end do
   write(6,*) 'Total solid messages size:',SUM(SUM(sizemsg_solid,DIM=1))
  end if

! TNM: DIM=1 ????! 
! af : means max along first dimension
  sizemsgmax_solid = maxval(maxval(sizemsg_solid,DIM=1))
  if (dump_mesh_info_screen) write(6,*) 'size msg max solid is ' , sizemsgmax_solid

  allocate(index_msg_solid(0:nproc-1,0:nproc-1))
  index_msg_solid(0:nproc-1,0:nproc-1) = 0

  allocate(global_index_msg_solid(sizemsgmax_solid,0:nproc-1,0:nproc-1))
  do iproct = 0, nproc-1
     do ipt = 1, sizebin_solid(iproct)
        ig = binp_solid(ipt,iproct)
        do ibel = 2, nprocb_solid(ig)
           ipdes = lprocb_solid(1,ig) ! same choice again
           ipsrc = lprocb_solid(ibel,ig)
           index_msg_solid(ipsrc,ipdes) = index_msg_solid(ipsrc,ipdes) + 1
           imsg = index_msg_solid(ipsrc,ipdes)
           global_index_msg_solid(imsg,ipsrc,ipdes) = ig
        end do
     end do
  end do

! TNM MAY07: NEEDS TO BE SYMMETRIC TOO then.... 
!  index_msg_solid(0:nproc-1,0:nproc-1) = 0
!  do ipdes = 0,nproc-1
!     do ipsrc = 0,nproc-1
!        index_msg_solid(ipsrc,ipdes) = index_msg_solid(ipsrc,ipdes) + 1
!        global_index_msg_solid(:,ipdes,ipsrc) = &
!                     global_index_msg_solid(:,ipsrc,ipdes)
!     enddo
!  enddo

! TNM: Why = 5 ??? changed this...
! ipsrc = 5
!  ipsrc = nproc-1
!  do ipdes = 0, nproc-1 
!     write(6,20) ipdes,sizemsg_solid(ipsrc,ipdes),&
!          (global_index_msg_solid(imsg,ipsrc,ipdes),&
!          imsg=1,sizemsg_solid(ipsrc,ipdes))
!  end do

  deallocate(index_msg_solid)

! How many messages am I sending / receiving ? 
  allocate(sizerecvp_solid(0:nproc-1),sizesendp_solid(0:nproc-1))
  sizerecvp_solid(0:nproc-1) = 0
  sizesendp_solid(0:nproc-1) = 0 
  do iproct = 0, nproc -1
     do ipsrc = 0, nproc-1
        if ( sizemsg_solid(ipsrc,iproct) > 0 ) &
             sizerecvp_solid(iproct) = sizerecvp_solid(iproct) + 1 
     end do
     do ipdes = 0, nproc-1 
        if ( sizemsg_solid(iproct,ipdes) > 0 ) &
             sizesendp_solid(iproct) = sizesendp_solid(iproct) + 1 
     end do
  end do


  if (dump_mesh_info_screen) then 
   do iproct = 0, nproc -1
     write(6,12)iproct,sizerecvp_solid(iproct),sizesendp_solid(iproct)
   end do
   12 format('Proc', i3, ' receives solid stuff from',i3, &
             ' procs and sends to',i3,' procs')
  end if

  sizerecvpmax_solid = maxval(sizerecvp_solid(:))
  sizesendpmax_solid = maxval(sizesendp_solid(:))

  if (dump_mesh_info_screen) then 
   write(6,*)'max size recv solid:',sizerecvpmax_solid
   write(6,*)'max size send solid:',sizerecvpmax_solid
  end if

! To which processors ? 
  allocate(listrecvp_solid(sizerecvpmax_solid,0:nproc-1))
  allocate(listsendp_solid(sizesendpmax_solid,0:nproc-1))
  listrecvp_solid(1:sizerecvpmax_solid,0:nproc-1) = -1 
  listsendp_solid(1:sizesendpmax_solid,0:nproc-1) = -1

  do iproct = 0, nproc -1
     ip = 0
     do ipsrc = 0, nproc -1
        if ( sizemsg_solid(ipsrc,iproct) > 0 ) then
           ip = ip + 1
           listrecvp_solid(ip,iproct) = ipsrc
        end if
     end do
     ip = 0
     do ipdes = 0, nproc -1
        if ( sizemsg_solid(iproct,ipdes) > 0 ) then
           ip = ip + 1
           listsendp_solid(ip,iproct) = ipdes
        end if
     end do
  end do
!
!af
  if (dump_mesh_info_screen .AND. nproc>1) then 
   do iproct = 0, nproc-1
      write(6,22) iproct, sizerecvp_solid(iproct), &
                  (listrecvp_solid(ip,iproct),ip=1,sizerecvp_solid(iproct))
      write(6,23) iproct, sizesendp_solid(iproct), &
                  (listsendp_solid(ip,iproct),ip=1,sizesendp_solid(iproct))
   end do
22 format('Proc', i3, ' will receive ', i2, ' solid messages from procs ', &
          20(i3,1x))
23 format('Proc', i3, ' will send    ', i2, ' solid messages to   procs ', &
           20(i3,1x))
  end if
  
! What size ? 
  allocate(sizemsgrecvp_solid(sizerecvpmax_solid,0:nproc-1))
  allocate(sizemsgsendp_solid(sizesendpmax_solid,0:nproc-1))

  do iproct = 0, nproc-1
     if (sizerecvp_solid(iproct)>0) then 
     do ip = 1, sizerecvp_solid(iproct)
        ipsrc = listrecvp_solid(ip,iproct)
        sizemsgrecvp_solid(ip,iproct) = sizemsg_solid(ipsrc,iproct)
     end do
     endif 

     if (sizesendp_solid(iproct)>0) then 
     do ip = 1, sizesendp_solid(iproct) 
        ipdes = listsendp_solid(ip,iproct)
        sizemsgsendp_solid(ip,iproct) = sizemsg_solid(iproct,ipdes)
     end do
     endif
  end do

! OUTPUT message size
  if (dump_mesh_info_screen) then 
   write(6,*)
   do iproct = 0, nproc -1
      if (sizerecvp_solid(iproct)>0) then 
      do ip = 1, sizerecvp_solid(iproct)
         write(6,17)iproct,listrecvp_solid(ip,iproct),&
                    sizemsgrecvp_solid(ip,iproct)
      end do
      end if

      if (sizesendp_solid(iproct)>0) then 
      do ip = 1, sizesendp_solid(iproct)
         write(6,18)iproct,listsendp_solid(ip,iproct),&
                    sizemsgsendp_solid(ip,iproct)
      end do
      endif
   end do
17 format('Proc',i3,' receiving solid message from',i3,' sized',i6)
18 format('Proc',i3,' sending   solid message to  ',i3,' sized',i6)
  end if

! NOW CREATE GLOCAL INDEX FOR MESSAGES
! Which glocal indices?
!!  
! ALEX TUESDAY Pm : must change the second indices everywhere !!!!
! TO sizemsgrecvpmax AND sizemsgsenpmax RESPECTIVELY
  allocate(glocal_index_msg_recvp_solid&
              (sizemsgmax_solid,sizerecvpmax_solid,0:nproc-1))
  allocate(glocal_index_msg_sendp_solid&
              (sizemsgmax_solid,sizesendpmax_solid,0:nproc-1))
!  allocate(glocal_index_msg_recvp_solid&
!              (sizemsgmax_solid,maxval(sizemsgrecvp_solid),0:nproc-1))
!  allocate(glocal_index_msg_sendp_solid&
!              (sizemsgmax_solid,maxval(sizemsgsendp_solid),0:nproc-1))

! END ALEX TUESDAY Pm
!!

  write(6,*)

  glocal_index_msg_recvp_solid = 0
  glocal_index_msg_sendp_solid = 0
  do iproct = 0, nproc-1
!! do iproct=0,0
!     write(6,*)
!     write(6,*) '--------------- PROC ',iproct,'------------------------'
     if (sizerecvp_solid(iproct)>0) then 
     do ip = 1, sizerecvp_solid(iproct)
        ipsrc = listrecvp_solid(ip,iproct)
        if (dump_mesh_info_screen) &
        write(6,*) ip,sizerecvp_solid(iproct),ipsrc, sizemsgrecvp_solid(ip,iproct)
        call flush(6)
        do ipt = 1, sizemsgrecvp_solid(ip,iproct)  
           ig = global_index_msg_solid(ipt,ipsrc,iproct)
           if (ig<1) then 
              write(6,*)ipsrc,ipt,ip,ig
           endif
!           write(6,*)'ipt,ip,ig',ipt,ip,ig; call flush(6)
           glocal_index_msg_recvp_solid(ipt,ip,iproct) = slob2sloc(ig,iproct)
        end do
     end do
     endif
     if (sizesendp_solid(iproct)>0) then 
     do ip = 1, sizesendp_solid(iproct)
        ipdes = listsendp_solid(ip,iproct)
        if (dump_mesh_info_screen) &
        write(6,*) ipdes, sizemsgsendp_solid(ip,iproct); call flush(6)
        do ipt = 1, sizemsgsendp_solid(ip,iproct)
           ig = global_index_msg_solid(ipt,iproct,ipdes)
           glocal_index_msg_sendp_solid(ipt,ip,iproct) = slob2sloc(ig,iproct)
        end do
     end do
     endif
  end do

21 format(80(i5,1x))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Write out domain for each processor and an array that is iproc on its 
! non-communicating points, and assumes the processor number of its partner
! for send/receive on these communicating global points.
! Use GMT script plot_proc_messaging.csh to plot processor domains 
!
do iproct=0,nproc-1
  allocate(uglob2_solid(nglobp_solid(iproct))) 
  allocate(val_solid(0:npol,0:npol,nel_solid(iproct)))
  call define_io_appendix(appiproc,iproct)

! send
  uglob2_solid(:) = iproct
  val_solid(:,:,:)=0
  if (sizesendp_solid(iproct)>0) then 
     do ip = 1, sizesendp_solid(iproct)
        ipdes = listsendp_solid(ip,iproct)
        do ipt = 1, sizemsgsendp_solid(ip,iproct)
           ig = global_index_msg_solid(ipt,iproct,ipdes)
           uglob2_solid(glocal_index_msg_sendp_solid(ipt,ip,iproct)) = ipdes
        enddo
     enddo

     do iel = 1, nel_solid(iproct)
        do ipol = 0, npol
           do jpol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
              idest = igloc_solid(ipt,iproct)
              val_solid(ipol,jpol,iel) = uglob2_solid(idest)
           enddo
        enddo
     enddo
  endif

  call define_io_appendix(appiproc,iproct)

! receive
  uglob2_solid(:) = iproct
  val_solid(:,:,:)=0

  if (sizerecvp_solid(iproct)>0) then 
     do ip = 1, sizerecvp_solid(iproct)
        ipdes = listrecvp_solid(ip,iproct)
        do ipt = 1, sizemsgrecvp_solid(ip,iproct)
           ig = global_index_msg_solid(ipt,iproct,ipdes)
           uglob2_solid(glocal_index_msg_recvp_solid(ipt,ip,iproct)) = ipdes
        enddo
     enddo

     do iel = 1, nel_solid(iproct)
        do ipol = 0, npol
           do jpol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
              idest = igloc_solid(ipt,iproct)
              val_solid(ipol,jpol,iel) = uglob2_solid(idest)
           enddo
        enddo
     enddo
  endif

  deallocate(uglob2_solid,val_solid)

enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  deallocate(sizebin_solid,binp_solid)
  deallocate(global_index_msg_solid)

  write(6,*)'End of solid messaging'; call flush(6) 
 
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FLUID FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF


if (have_fluid) then
  if (dump_mesh_info_screen) write(6,*) ' Creating fluid bins ' 
  allocate (sizebin_fluid(0:nproc-1))
  sizebin_fluid(0:nproc-1) = 0 
  do ipt = 1, nglobflob
     iproct = lprocb_fluid(1,ipt)! completely ad hoc choice
     sizebin_fluid(iproct)=sizebin_fluid(iproct) + 1
  end do

  if (dump_mesh_info_screen) then 
   do iproct = 0, nproc - 1
     write(6,*) 'proc:',iproct,'size fluid bin:',sizebin_fluid(iproct)
   end do
  end if

  if (dump_mesh_info_screen)  write(6,*) SUM(sizebin_fluid),nglobflob
  if (SUM(sizebin_fluid)/=nglobflob) then 
     write(6,*)'PROBLEM: sum of fluid bins not equal to flobal num!'
     stop
  endif
  sizebinmax_fluid = maxval (sizebin_fluid(:)) 

  if (dump_mesh_info_screen) &
  write(6,*) ' sizebinmax_fluid = ' , sizebinmax_fluid

  allocate(binp_fluid(sizebinmax_fluid,0:nproc-1))
  binp_fluid(1:sizebinmax_fluid,0:nproc-1) = 0 

  allocate(ibin_fluid(0:nproc-1))
  ibin_fluid(0:nproc-1) = 0

  do ipt = 1, nglobflob
     iproct =  lprocb_fluid(1,ipt)! same choice 
     ibin_fluid(iproct) = ibin_fluid(iproct) + 1
     ibp = ibin_fluid(iproct)
     binp_fluid(ibp,iproct) = ipt
  end do
  deallocate(ibin_fluid)

  ! do iproct = 0, nproc-1
  !  write(6,20) (binp(ibp,iproct),ibp=1,sizebin(iproct))
  ! end do
  allocate(sizemsg_fluid(0:nproc-1,0:nproc-1))
  sizemsg_fluid(0:nproc-1,0:nproc-1) = 0 

!  write(6,*) ' FLUID DEPENDENCIES ' 
  do iproct = 0, nproc - 1
     !  write(6,*) ' PROC IS ' , iproct
     do ipt = 1, sizebin_fluid(iproct)
        ig = binp_fluid(ipt,iproct)
        !   write(6,*) nprocb(ig)
        do ibel = 2, nprocb_fluid(ig)
           !    write(6,*) 'iglob =',ig, 'nprocb = ',nprocb(ig), &
           ! (lprocb(ibel,ig),ibel=2,nprocb(ig))
           ipdes = lprocb_fluid(1,ig) ! same choice 
           ipsrc = lprocb_fluid(ibel,ig)
           sizemsg_fluid(ipsrc,ipdes) = sizemsg_fluid(ipsrc,ipdes) + 1 
        end do
     end do
  end do

! TNM MAY07: MAKE THIS SYMMETRIC; send/receive need to happen in both ways...
! NOTE THE ORDER OF THE LOOP! Other way around is ZERO throughout............
!  do ipdes = 0,nproc-1
!     do ipsrc = 0,nproc-1
!        sizemsg_fluid(ipdes,ipsrc) = sizemsg_fluid(ipsrc,ipdes)
!     enddo
!  enddo

  if (dump_mesh_info_screen) then 
   write(6,*)'Size of solid messages for each proc-proc pair:'
   write(6,*)'---> destination proc, down: my proc'
   do iproct = 0,  nproc-1
      write(6,21) (sizemsg_fluid(iproct,ipdes),ipdes=0,nproc-1)
   end do
   write(6,*) 'Total fluid messages size:',SUM(SUM(sizemsg_fluid,DIM=1))
  end if
  sizemsgmax_fluid = maxval(maxval(sizemsg_fluid,DIM=1))
  if (dump_mesh_info_screen) write(6,*) 'size msg max fluid is ' , sizemsgmax_fluid

  allocate(index_msg_fluid(0:nproc-1,0:nproc-1))
  index_msg_fluid(0:nproc-1,0:nproc-1) = 0

  allocate(global_index_msg_fluid(sizemsgmax_fluid,0:nproc-1,0:nproc-1))
  do iproct = 0, nproc-1
     do ipt = 1, sizebin_fluid(iproct)
        ig = binp_fluid(ipt,iproct)
        do ibel = 2, nprocb_fluid(ig)
           ipdes = lprocb_fluid(1,ig) ! same choice again
           ipsrc = lprocb_fluid(ibel,ig)
           index_msg_fluid(ipsrc,ipdes) = index_msg_fluid(ipsrc,ipdes) + 1
           imsg = index_msg_fluid(ipsrc,ipdes)
           global_index_msg_fluid(imsg,ipsrc,ipdes) = ig
        end do
     end do
  end do

! TNM MAY07: NEEDS TO BE SYMMETRIC TOO then.... 
!  index_msg_fluid(0:nproc-1,0:nproc-1) = 0
!  do ipdes = 0,nproc-1
!     do ipsrc = 0,nproc-1
!        index_msg_fluid(ipsrc,ipdes) = index_msg_fluid(ipsrc,ipdes) + 1
!        global_index_msg_fluid(:,ipdes,ipsrc) = &
!                     global_index_msg_fluid(:,ipsrc,ipdes)
!     enddo
!  enddo

! TNM: Why = 5 ??? changed this...
! ipsrc = 5
!  ipsrc = nproc-1
!  do ipdes = 0, nproc-1 
!     write(6,20) ipdes,sizemsg_fluid(ipsrc,ipdes),&
!          (global_index_msg_fluid(imsg,ipsrc,ipdes),&
!          imsg=1,sizemsg_fluid(ipsrc,ipdes))
!  end do

  deallocate(index_msg_fluid)

! How many messages am I sending / receiving ? 
  allocate(sizerecvp_fluid(0:nproc-1),sizesendp_fluid(0:nproc-1))
  sizerecvp_fluid(0:nproc-1) = 0
  sizesendp_fluid(0:nproc-1) = 0 
  do iproct = 0, nproc -1
     do ipsrc = 0, nproc-1
        if ( sizemsg_fluid(ipsrc,iproct) > 0 ) &
             sizerecvp_fluid(iproct) = sizerecvp_fluid(iproct) + 1 
     end do
     do ipdes = 0, nproc-1 
        if ( sizemsg_fluid(iproct,ipdes) > 0 ) &
             sizesendp_fluid(iproct) = sizesendp_fluid(iproct) + 1 
     end do
  end do

  if (dump_mesh_info_screen) then 
   do iproct = 0, nproc -1
      write(6,13) iproct,sizerecvp_fluid(iproct),sizesendp_fluid(iproct)
   end do
  13 format('Proc', i3, ' receives fluid stuff from',i3, &
             ' procs and sends to',i3,' procs')
  end if

  sizerecvpmax_fluid = maxval(sizerecvp_fluid(:))
  sizesendpmax_fluid = maxval(sizesendp_fluid(:))
 
  if (dump_mesh_info_screen) then 
   write(6,*)'max size recv fluid:',sizerecvpmax_fluid
   write(6,*)'max size send fluid:',sizerecvpmax_fluid
  end if

! To which processor ? 
  allocate(listrecvp_fluid(sizerecvpmax_fluid,0:nproc-1))
  allocate(listsendp_fluid(sizesendpmax_fluid,0:nproc-1))
  listrecvp_fluid(1:sizerecvpmax_fluid,0:nproc-1) = -1
  listsendp_fluid(1:sizesendpmax_fluid,0:nproc-1) = -1

  do iproct = 0, nproc -1
     ip = 0
     do ipsrc = 0, nproc -1
        if ( sizemsg_fluid(ipsrc,iproct) > 0 ) then
           ip = ip + 1
           listrecvp_fluid(ip,iproct) = ipsrc
        end if
     end do
     ip = 0
     do ipdes = 0, nproc -1
        if ( sizemsg_fluid(iproct,ipdes) > 0 ) then
           ip = ip + 1
           listsendp_fluid(ip,iproct) = ipdes
        end if
     end do
  end do
!
!af
  if (dump_mesh_info_screen .and. nproc>1) then 
   do iproct = 0, nproc-1
      write(6,32) iproct, sizerecvp_fluid(iproct), &
                  (listrecvp_fluid(ip,iproct),ip=1,sizerecvp_fluid(iproct))
      write(6,33) iproct, sizesendp_fluid(iproct), &
                  (listsendp_fluid(ip,iproct),ip=1,sizesendp_fluid(iproct))
   end do
32 format('Proc', i2, ' will receive ', i2, ' fluid messages from procs ', &
          20(i3,1x))
33 format('Proc', i2, ' will send    ', i2, ' fluid messages to   procs ', &
           20(i3,1x))
  end if

! What size ? 
  allocate(sizemsgrecvp_fluid(sizerecvpmax_fluid,0:nproc-1))
  allocate(sizemsgsendp_fluid(sizesendpmax_fluid,0:nproc-1))
  do iproct = 0, nproc-1
     if (sizerecvp_fluid(iproct)>0) then 
     do ip = 1, sizerecvp_fluid(iproct)
        ipsrc = listrecvp_fluid(ip,iproct)
        sizemsgrecvp_fluid(ip,iproct) = sizemsg_fluid(ipsrc,iproct)
     end do
     endif

     if (sizesendp_fluid(iproct)>0) then 
     do ip = 1, sizesendp_fluid(iproct) 
        ipdes = listsendp_fluid(ip,iproct)
        sizemsgsendp_fluid(ip,iproct) = sizemsg_fluid(iproct,ipdes)
     end do
     endif
  end do

! OUTPUT message size
  if (dump_mesh_info_screen .and. nproc>1) then
   write(6,*)
   do iproct = 0, nproc -1
      if (sizerecvp_fluid(iproct)>0) then 
      do ip = 1, sizerecvp_fluid(iproct)
         write(6,37)iproct,listrecvp_fluid(ip,iproct),&
                    sizemsgrecvp_fluid(ip,iproct)
      end do
      endif

      if (sizesendp_fluid(iproct)>0) then
      do ip = 1, sizesendp_fluid(iproct)
         write(6,38)iproct,listsendp_fluid(ip,iproct),&
                    sizemsgsendp_fluid(ip,iproct)
      end do
      endif
   end do
37 format('Proc',i3,' receiving fluid message from',i3,' sized',i6)
38 format('Proc',i3,' sending   fluid message to  ',i3,' sized',i6)
  end if

! NOW CREATE GLOCAL INDEX FOR MESSAGES 
! Which glocal indices?
!!  
! ALEX TUESDAY Pm : must change the second indices everywhere !!!!
! TO sizemsgrecvpmax AND sizemsgsenpmax RESPECTIVELY
  allocate(glocal_index_msg_recvp_fluid&
              (sizemsgmax_fluid,sizerecvpmax_fluid,0:nproc-1))
  allocate(glocal_index_msg_sendp_fluid&
              (sizemsgmax_fluid,sizesendpmax_fluid,0:nproc-1))
!  allocate(glocal_index_msg_recvp_fluid&
!              (sizemsgmax_fluid,maxval(sizerecvp_fluid),0:nproc-1))
!  allocate(glocal_index_msg_sendp_fluid&
!              (sizemsgmax_fluid,maxval(sizesendp_fluid),0:nproc-1))

! END ALEX TUESDAY Pm
!!
  glocal_index_msg_recvp_fluid = 0
  glocal_index_msg_sendp_fluid = 0
  do iproct = 0, nproc-1
!! do iproct=0,0
!     write(6,*)
!     write(6,*) '--------------- PROC ',iproct,'------------------------'
!     write(6,*) iproct,'  SHALL RECEIVE: '
!!   write(6,*) iproct,sizerecvp_fluid(iproct)
     if (sizerecvp_fluid(iproct)>0) then 
     do ip = 1, sizerecvp_fluid(iproct)
        ipsrc = listrecvp_fluid(ip,iproct)
!        write(6,*) ipsrc, sizemsgrecvp_fluid(ip,iproct)
        do ipt = 1, sizemsgrecvp_fluid(ip,iproct)  
           ig = global_index_msg_fluid(ipt,ipsrc,iproct)
           glocal_index_msg_recvp_fluid(ipt,ip,iproct) = flob2floc(ig,iproct)
       !! write(6,*) ipsrc,ipt,ig,glocal_index_msg_recvp_fluid(ipt,ip,iproct)
        end do
     end do
     endif
!     write(6,*) iproct,'  SHALL SEND   : '

     if (sizesendp_fluid(iproct)>0) then 
     do ip = 1, sizesendp_fluid(iproct)
        ipdes = listsendp_fluid(ip,iproct)
!        write(6,*) ipdes, sizemsgsendp_fluid(ip,iproct)
        do ipt = 1, sizemsgsendp_fluid(ip,iproct)
           ig = global_index_msg_fluid(ipt,iproct,ipdes)
           glocal_index_msg_sendp_fluid(ipt,ip,iproct) = flob2floc(ig,iproct)
!           write(6,*) ipdes,ipt,ig,glocal_index_msg_sendp_fluid(ipt,ip,iproct)
        end do
     end do
     endif
  end do

end if ! have_fluid


  if (dump_mesh_info_screen .and. nproc>1) then 
   write(6,*)
   write(6,*)'----------------------------------------------------------------'
   write(6,*)'Sum over all solid message sizes  :',SUM(SUM(sizemsg_solid,DIM=1))
   write(6,*)'Total global points in solid,ratio:',nglobslob,&
                            real(SUM(SUM(sizemsg_solid,DIM=1)))/real(nglobslob)
   write(6,*)'----------------------------------------------------------------'

if (have_fluid) then
   write(6,*)'----------------------------------------------------------------'
   write(6,*)'Sum over all fluid message sizes  :',SUM(SUM(sizemsg_fluid,DIM=1))
   write(6,*)'Total global points in fluid,ratio:',nglobflob,&
                            real(SUM(SUM(sizemsg_fluid,DIM=1)))/real(nglobflob)
   write(6,*)'----------------------------------------------------------------'
   write(6,*)'----------------------------------------------------------------'
   write(6,*)'Sum over all s/f message sizes:',SUM(SUM(sizemsg_fluid,DIM=1))+&
                                               SUM(SUM(sizemsg_solid,DIM=1)) 
   write(6,*)'Total global points, ratio    :',nglobglob,&
                real(SUM(SUM(sizemsg_fluid,DIM=1))+&
                     SUM(SUM(sizemsg_solid,DIM=1)))/real(nglobglob)
   write(6,*)'----------------------------------------------------------------'
 end if !have_fluid 

   write(6,*)
   write(6,*)'****************************************************************'
   write(6,*)'******************* END OF MESSAGING ARRAYS ********************'
   write(6,*)'****************************************************************'
  end if


if (have_fluid) then

  deallocate(sizemsg_fluid)
  deallocate(sizebin_fluid,binp_fluid,global_index_msg_fluid)

end if ! have_fluid

  deallocate(sizemsg_solid)

end subroutine partition_sflobal_index
!------------------------------------------------------------------------

!dk define_local_bdry_elem---------------------------------------------------
subroutine define_local_bdry_elem

  use global_parameters
  use data_gllmesh, only : sgll, zgll

integer ::iel,ielg,iproc,j,solid_count(0:nproc-1),fluid_count(0:nproc-1)
integer, allocatable, dimension(:,:) :: tmpsolid,tmpfluid,jpolsol,jpolflu
double precision :: rbound
integer :: myproc,herproc,myielglob,herielglob

! test parameters
integer count,ipol,jpol
double precision, allocatable :: tmpsolfield(:,:,:,:), tmpflufield(:,:,:,:)
double precision :: rtmp

  allocate(nbdry_el(0:nproc-1),have_bdry_elemp(0:nproc-1))
  allocate(belemp(sum(nbelem),0:nproc-1))
  allocate(tmpsolid(sum(nbelem)/2,0:nproc-1),tmpfluid(sum(nbelem)/2,0:nproc-1))
  allocate(jpolsol(sum(nbelem)/2,0:nproc-1),jpolflu(sum(nbelem)/2,0:nproc-1))

  tmpsolid = -1; tmpfluid = -1

  solid_count(:)=0; fluid_count(:)=0
  nbdry_el(:)=0
  do j=1,nbcnd

   if (mod(j,2)/=0) then ! upper boundary of fluid region
      rbound = discont(idom_fluid(j))/router
   else  ! lower boundary of fluid region
      rbound = discont(idom_fluid(j-1)+1)/router
   endif

     do iel=1,nbelem(j)
        nbdry_el(el2proc(belem(iel,j)))=nbdry_el(el2proc(belem(iel,j)))+1

! belemp: given element counter on boundary for a proc, return glocal element number
        belemp(nbdry_el(el2proc(belem(iel,j))),el2proc(belem(iel,j)))=&
             inv_procel(belem(iel,j),el2proc(belem(iel,j)))

!       write(1010,*)iel,nbdry_el(el2proc(belem(iel,j))),el2proc(belem(iel,j)), &
!                    belem(iel,j),&
!                    belemp(nbdry_el(el2proc(belem(iel,j))),el2proc(belem(iel,j)))

!       call flush(1010)

        if (solid(belem(iel,j))) then

           myielglob = belem(iel,j)
           herielglob = belem(my_neighbour(iel,j),j)

           myproc =  el2proc(myielglob)
           herproc = el2proc(herielglob)

           if (herproc /= myproc) then 
              write(6,*)
              write(6,*)'PROBLEM: changing processors across solid/fluid boundary!'
              write(6,*)'This case is not implemented at this point as it requires'
              write(6,*)'message passing when adding the boundary term on both sides'
              write(6,*)'and therefore seems highly ineffective. Sorry...'
              write(6,*)'solid domain proc:',myproc
              write(6,*)'fluid domain proc:',herproc
              write(6,*)'boundary,global element num:',j,iel,myielglob
              stop
           endif

           solid_count(myproc)=solid_count(myproc)+1
                
! Determine jpol at the boundary: 0,npol depending on above/below, North/South
           if (rcom(myielglob) > rbound) then !solid above solid-fluid boundary
              if (zcom(myielglob)>=0.d0) then ! North
                 jpolsol(solid_count(myproc),myproc)=0
                 jpolflu(solid_count(myproc),myproc)=npol
              else ! South
                 jpolsol(solid_count(myproc),myproc)=npol
                 jpolflu(solid_count(myproc),myproc)=0
              endif
           else ! solid below solid-fluid boundary
              if (zcom(myielglob)>=0.d0) then ! North
                 jpolsol(solid_count(myproc),myproc)=npol
                 jpolflu(solid_count(myproc),myproc)=0
              else ! South
                 jpolsol(solid_count(myproc),myproc)=0
                 jpolflu(solid_count(myproc),myproc)=npol
              endif                  
! crude way to accomodate the case of having the buffer layer below the ICB.
! ...i.e., the jpol indices are switched for 45 < theta < 135 deg
              if (eltypeg(myielglob)/='curved') then 
                 if ( scom(myielglob) > abs(zcom(myielglob)) ) then 
                    jpolsol(solid_count(myproc),myproc)=&
                         abs(jpolsol(solid_count(myproc),myproc)-npol)
                 endif
              endif

           endif

           tmpsolid ( solid_count(myproc), myproc) = &
              inv_procel_solidp( inv_procel(myielglob,myproc), myproc)

           tmpfluid(solid_count(myproc),herproc) = &
              inv_procel_fluidp( inv_procel(herielglob,herproc), herproc)

        else
           fluid_count(myproc)=fluid_count(myproc)+1
        endif
     enddo
  enddo

  if (dump_mesh_info_screen)&
   write(6,*)'ended solid-fluid boundary loop'; call flush(6)

  if (sum(solid_count) .ne. sum(fluid_count) ) then
     write(6,*)'Something wrong with the # solid/fluid bdry elements:'
     write(6,*)'Solid count:',sum(solid_count)
     write(6,*)'Fluid count:',sum(fluid_count)
     stop
  endif
  
  nbdry_el=solid_count

  allocate(bdry_solid_elp(maxval(nbdry_el),0:nproc-1))
  allocate(bdry_fluid_elp(maxval(nbdry_el),0:nproc-1))

  do iproc=0,nproc-1
     if (minval(tmpsolid(1:nbdry_el(iproc),iproc)) <1) then 
        write(6,*)'Problem with bdry_solid count!'
        write(6,*)minval(tmpsolid(1:nbdry_el(iproc),iproc))
        write(6,*)'Unassigned elements...'
        stop
     endif
     if (minval(tmpfluid(1:nbdry_el(iproc),iproc)) <1) then 
        write(6,*)'Problem with bdry_fluid count!'
        write(6,*)minval(tmpfluid(1:nbdry_el(iproc),iproc))
        write(6,*)'Unassigned elements...'
        stop
     endif
     bdry_solid_elp(1:nbdry_el(iproc),iproc)=tmpsolid(1:nbdry_el(iproc),iproc)
     bdry_fluid_elp(1:nbdry_el(iproc),iproc)=tmpfluid(1:nbdry_el(iproc),iproc)
     if (dump_mesh_info_screen) then
      write(6,101)iproc,nbdry_el(iproc)
      if (maxval(nbdry_el)>=1) write(6,*)bdry_solid_elp(1,iproc),bdry_fluid_elp(1,iproc)
     end if
  enddo

101 format('Proc', i3, ' has ', i5, ' S/F boundary elements (on one side)')

  if (dump_mesh_info_screen) then 
   write(6,*)'max bdry_solid_elp,max loc el:',&
              maxval(bdry_solid_elp),maxval(nel_solid)
   if (have_fluid) &
   write(6,*)'max bdry_fluid_elp,max loc el:',&
              maxval(bdry_fluid_elp),maxval(nel_fluid)
  end if

  if (have_fluid) then
     allocate(bdry_jpol_solidp(maxval(nbdry_el),0:nproc-1))
     bdry_jpol_solidp(:,:)=jpolsol(1:maxval(nbdry_el),:)

     allocate(bdry_jpol_fluidp(maxval(nbdry_el),0:nproc-1))
     bdry_jpol_fluidp(:,:)=jpolflu(1:maxval(nbdry_el),:) 
  endif

! boolean array to check whether proc has solid-fluid boundary
  have_bdry_elemp(:)=.false.
  do iproc=0,nproc-1
     if (nbdry_el(iproc)>0 ) have_bdry_elemp(iproc)=.true.
  enddo

  
  if (dump_mesh_info_screen .and. have_fluid) then 
   write(6,*)'BOUNDARY TERMS:'
   do iproc=0,nproc-1
      write(6,*)iproc,'minmax bdry sol:',minval(bdry_solid_elp(1:nbdry_el(iproc),iproc)), &
                maxval(bdry_solid_elp(1:nbdry_el(iproc),iproc))
      write(6,*)iproc,'minmax bdry flu:',minval(bdry_fluid_elp(1:nbdry_el(iproc),iproc)), &
                maxval(bdry_fluid_elp(1:nbdry_el(iproc),iproc))
   enddo
  end if


 return







! ===========  T E S  T   F R O M   H E R E   O N  ===========================
! MvD 12/2012: this test is never executed (return statement above), causes a
! segfault and uses loads of memory. 


  !+++++++++++++
  if (have_fluid) then
  !++++++++++++

     if (dump_mesh_info_screen) then 
        write(6,*)'Solid-Fluid boundary tests'; call flush(6)
     end if

     !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
     ! New test: define field on fluid side, copy to solid, check difference
     !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     allocate(tmpflufield(maxval(nel_fluid),0:npol,0:npol,0:nproc-1))
     do iproc=0,nproc-1
        do iel=1,nel_fluid(iproc)
           ielg = procel_fluid(iel,iproc)
           do ipol=0,npol
              do jpol=0,npol
                 ! make tmpflufield equal to theta
                 rtmp=sqrt(sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
                 tmpflufield(iel,ipol,jpol,iproc)= &
                    asin(sgll(ipol,jpol,ielg)/rtmp)*180./pi
              enddo
           enddo
        enddo
     enddo
     
     if (dump_mesh_info_screen)&
               write(6,*)'Test bdry: maxval/loc fluid [deg]:', &
               maxval(tmpflufield),maxloc(tmpflufield); call flush(6)

     ! Now copy the S/F boundary values to the solid domain
     allocate(tmpsolfield(maxval(nel_solid),0:npol,0:npol,0:nproc-1))
     tmpsolfield(1:maxval(nel_solid),0:npol,0:npol,0:nproc-1)=-45.

     do iproc=0,nproc-1
        do iel=1,nbdry_el(iproc)
           tmpsolfield(bdry_solid_elp(iel,iproc),0:npol, &
                       bdry_jpol_solidp(iel,iproc),iproc)= &
                       tmpflufield(bdry_fluid_elp(iel,iproc), &
                       0:npol,bdry_jpol_fluidp(iel,iproc),iproc)
        enddo
     enddo

     if (dump_mesh_info_screen) &
                write(6,*)'Test bdry: maxval/loc solid [deg/indices]:', &
                maxval(tmpsolfield),maxloc(tmpsolfield); call flush(6)

        ! Write out the radii at which the solid field takes values >1 
        ! (i.e., if other than S/F boundary radii, something's wrong....)

     if (dump_mesh_info_screen .and. have_fluid) then 
        write(6,*)'Testing S/F boundary copying...'
        write(6,*)'File: bdrytest_solflubdry.dat contains:'
        write(6,*)'iproc,iel_solid, rad,latitude, value at which things have been copied:'
     end if
     count=0
     do iproc=0,nproc-1
        do iel=1,nel_solid(iproc)
           ielg = procel_solid(iel,iproc)
           do ipol=0,npol
              do jpol=0,npol
                 rtmp=sqrt(sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
                 if ( tmpsolfield(iel,ipol,jpol,iproc) >= 0.) then 
                    count=count+1
                 endif
              enddo
           enddo
        enddo
     enddo

110 format(i3,i9,3(1pe14.3))


  if (count .ne. sum(nbdry_el)*(npol+1) ) then
     write(6,*)
     write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(6,*)'E R R O R at S/F boundary copying!'
     write(6,*)'expected # bdry points        :',sum(nbdry_el)*(npol+1)
     write(6,*)'actually copied # bdry points :',count
     write(6,*)'...see file bdrytest_solflubdry.dat for details...'; call flush(6)
     stop
  endif

  !+++++++++++++
  endif ! have_fluid
  !++++++++++++




!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! End New test: define field on fluid side, copy to solid, check difference
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

end subroutine define_local_bdry_elem
!------------------------------------------------------------------------


!dk define_axial_elem----------------------------------------------------
subroutine define_axial_elem

use data_gllmesh

integer :: iel, iproc, ielg
!integer, dimension(maxval(nel),0:nproc-1) :: dummyax_elp
!integer, dimension(maxval(nel_solid),0:nproc-1) :: dummyax_el_solidp
!integer, dimension(maxval(nel_fluid),0:nproc-1) :: dummyax_el_fluidp
integer, dimension(:,:), allocatable :: dummyax_elp
integer, dimension(:,:), allocatable :: dummyax_el_solidp
integer, dimension(:,:), allocatable :: dummyax_el_fluidp
!af
allocate(dummyax_elp(maxval(nel),0:nproc-1))
allocate(dummyax_el_solidp(maxval(nel_solid),0:nproc-1))
allocate(dummyax_el_fluidp(maxval(nel_fluid),0:nproc-1))

allocate(naxelp(0:nproc-1),naxel_solidp(0:nproc-1),naxel_fluidp(0:nproc-1))
allocate(have_axis(0:nproc-1))
allocate(axis(maxval(nel),0:nproc-1))
allocate(axis_solid(maxval(nel_solid),0:nproc-1))
allocate(axis_fluid(maxval(nel_fluid),0:nproc-1))

axis=0; axis_solid=0; axis_fluid=0
dummyax_elp=0; dummyax_el_solidp=0; dummyax_el_fluidp=0
have_axis=.false.
naxelp=0; naxel_fluidp=0; naxel_solidp=0;

if (dump_mesh_info_screen) then 
 write(6,*)
 write(6,*)'Defining axial elements:'
end if

do iproc=0,nproc-1

! glocal domain
   do iel=1,nel(iproc)
      ! jpol=npol is random choice
      ielg = procel(iel,iproc)
      
      if ( sgll(0,npol,ielg) < min_distance_nondim ) then 
         axis(iel,iproc)=1
         naxelp(iproc)=naxelp(iproc)+1
         dummyax_elp(naxelp(iproc),iproc)=iel
      endif
   enddo ! iel

   if (naxelp(iproc)>0) have_axis(iproc)=.true.

!slocal domain
   do iel=1,nel_solid(iproc)
      ielg = procel_solid(iel,iproc)
      ! jpol=npol is random choice

      if ( sgll(0,npol,ielg) < min_distance_nondim ) then 
         axis_solid(iel,iproc)=1
         naxel_solidp(iproc)=naxel_solidp(iproc)+1
         dummyax_el_solidp(naxel_solidp(iproc),iproc)=iel
      endif
   enddo ! iel

!flocal domain
   do iel=1,nel_fluid(iproc)
      ielg = procel_fluid(iel,iproc)
      ! jpol=npol is random choice
      if ( sgll(0,npol,ielg) < min_distance_nondim ) then 
         axis_fluid(iel,iproc)=1
         naxel_fluidp(iproc)=naxel_fluidp(iproc)+1
         dummyax_el_fluidp(naxel_fluidp(iproc),iproc)=iel
      endif
   enddo ! iel

enddo !iproc

write(6,*)'NAXEL:',naxelp(0),naxel_fluidp(0),min_distance_nondim

allocate(ax_elp(maxval(naxelp),0:nproc-1))
allocate(ax_el_solidp(maxval(naxel_solidp),0:nproc-1))
allocate(ax_el_fluidp(maxval(naxel_fluidp),0:nproc-1))

do iproc=0,nproc-1
   if (dump_mesh_info_screen) &
   write(6,20)iproc,naxelp(iproc),naxel_solidp(iproc),naxel_fluidp(iproc)
20 format(' Proc', i3, ' has ', i4,' total, ',i4,&
          ' solid, and ',i4,' fluid axial elements')

ax_elp(1:naxelp(iproc),iproc)=dummyax_elp(1:naxelp(iproc),iproc)
ax_el_solidp(1:naxel_solidp(iproc),iproc) = &
            dummyax_el_solidp(1:naxel_solidp(iproc),iproc)
ax_el_fluidp(1:naxel_fluidp(iproc),iproc) = &
            dummyax_el_fluidp(1:naxel_fluidp(iproc),iproc)

if (naxelp(iproc) /=naxel_solidp(iproc)+naxel_fluidp(iproc)) then
   write(6,*)'Processor',iproc,': PROBLEM in counting axial elements:'
   write(6,*)'Solid axial elements:',naxel_solidp(iproc)
   write(6,*)'Fluid axial elements:',naxel_fluidp(iproc)
   write(6,*)'Total axial elements:',naxelp(iproc)
   stop
endif

enddo

deallocate(dummyax_elp)
deallocate(dummyax_el_solidp)
deallocate(dummyax_el_fluidp)

end subroutine define_axial_elem
!------------------------------------------------------------------------

!dk generate_serendipity_per_proc----------------------------------------
  subroutine generate_serendipity_per_proc(sg,zg)
  use numbering
  double precision, dimension(4,neltot), intent(in) :: sg,zg

  integer :: iproc,nelp,iel,ielg
  integer :: npointotp,nelpmax,nglobmeshpmax
  integer :: ncp
!  real, dimension(:,:,:), allocatable :: sgp, zgp
!  real, dimension(:,:), allocatable :: sgpw,zgpw
  double precision, dimension(:,:,:), allocatable :: sgp, zgp
  double precision, dimension(:,:), allocatable :: sgpw,zgpw
  integer, dimension(:), allocatable :: wiglob,wloc
  logical, dimension(:), allocatable :: wifseg
  integer :: wnglob,iptcp,inode
  character(len=4) :: appiproc

  if (dump_mesh_info_screen) then 
   write(6,*) 
   write(6,*) ' SERENDIPITY DB '
   write(6,*) 
  end if

  ncp = 8 ! exclusively using the serendipity quadrilateral element topology

! FIRST DETERMINE NUMBER OF CONTROL POINTS FOR A GIVEN PROCESSOR
  allocate(nglobmeshp(0:nproc-1))
  nelpmax = maxval(nel(:))

! elemental physical coordinates of serendipity control nodes, parallel
  allocate(sgp(8,nelpmax,0:nproc-1),zgp(8,nelpmax,0:nproc-1))
  sgp(8,nelpmax,0:nproc-1) = 0.d0
  zgp(8,nelpmax,0:nproc-1) = 0.d0

! Global number of control nodes, parallel
  allocate(iglobcp(nelpmax*8,0:nproc-1))
  iglobcp(nelpmax*8,0:nproc-1) = 0
  do iproc = 0, nproc-1
     nelp = nel(iproc)
     npointotp = nelp*8
     do iel = 1, nelp
        ielg = procel(iel,iproc)
        sgp(1,iel,iproc) = sg(1,ielg) ; zgp(1,iel,iproc) = zg(1,ielg)
        sgp(3,iel,iproc) = sg(2,ielg) ; zgp(3,iel,iproc) = zg(2,ielg)
        sgp(5,iel,iproc) = sg(3,ielg) ; zgp(5,iel,iproc) = zg(3,ielg)
        sgp(7,iel,iproc) = sg(4,ielg) ; zgp(7,iel,iproc) = zg(4,ielg)
        sgp(2,iel,iproc) = .5d0*(sgp(1,iel,iproc)+sgp(3,iel,iproc))
        sgp(4,iel,iproc) = .5d0*(sgp(3,iel,iproc)+sgp(5,iel,iproc))
        sgp(6,iel,iproc) = .5d0*(sgp(5,iel,iproc)+sgp(7,iel,iproc))
        sgp(8,iel,iproc) = .5d0*(sgp(7,iel,iproc)+sgp(1,iel,iproc))
        zgp(2,iel,iproc) = .5d0*(zgp(1,iel,iproc)+zgp(3,iel,iproc))
        zgp(4,iel,iproc) = .5d0*(zgp(3,iel,iproc)+zgp(5,iel,iproc))
        zgp(6,iel,iproc) = .5d0*(zgp(5,iel,iproc)+zgp(7,iel,iproc))
        zgp(8,iel,iproc) = .5d0*(zgp(7,iel,iproc)+zgp(1,iel,iproc))
     end do

     allocate(sgpw(8,nelp),zgpw(8,nelp))
     do iel = 1, nelp
        sgpw(:,iel) = sgp(:,iel,iproc)
        zgpw(:,iel) = zgp(:,iel,iproc)
     end do

     allocate(wiglob(npointotp),wloc(npointotp),wifseg(npointotp))
     wiglob(1:npointotp) = 0 
     wloc(1:npointotp) = 0

     call get_global(nelp,sgpw,zgpw,wiglob,wloc,wifseg,wnglob,npointotp,ncp,ndim)

     nglobmeshp(iproc) = wnglob
     if (dump_mesh_info_screen) write(6,*) ' iproc = ', iproc, ' nglobmesh = ', nglobmeshp(iproc)
     do iel = 1, nelp
        do inode = 1, 8
           iptcp = (iel-1)*8 + inode
           iglobcp(iptcp,iproc) = wiglob(iptcp)
        end do
     end do
     deallocate(wifseg,wloc,wiglob)
     deallocate(zgpw,sgpw)
  end do

! COORDINATES OF THESE CONTROL POINTS (GLOCALLY NUMBERED)
! TOPOLOGY OF THE ELEMENTS
  nglobmeshpmax = maxval(nglobmeshp)
  allocate(scpp(nglobmeshpmax,0:nproc-1),zcpp(nglobmeshpmax,0:nproc-1))
  scpp(nglobmeshpmax,0:nproc-1) = 0.d0
  zcpp(nglobmeshpmax,0:nproc-1) = 0.d0 
  allocate(lnodescp(nelpmax,8,0:nproc-1))

  do iproc = 0, nproc - 1
     nelp = nel(iproc)
     do iel = 1, nelp
!        ielg = procel(iel,iproc) ! NO NEED???
        do inode = 1, 8
           iptcp = (iel-1)*8 + inode

! glocally numbered physical coordinates of serendipity control nodes, parallel
           scpp(iglobcp(iptcp,iproc),iproc) = sgp(inode,iel,iproc)
           zcpp(iglobcp(iptcp,iproc),iproc) = zgp(inode,iel,iproc)

! global number of control nodes stored by element, parallel
           lnodescp(iel,inode,iproc) = iglobcp(iptcp,iproc)
        end do
     end do
  end do

  end subroutine generate_serendipity_per_proc
!------------------------------------------------------------------------

!dk define_element_type--------------------------------------------------
  subroutine define_element_type

  integer :: nelp,nelpmax
  integer :: iproc,iel,ielg

  if (dump_mesh_info_screen) write(6,*) ' DEFINING ELEMENT TYPE ARRAY FOR EACH PROCESSOR '
  nelpmax=maxval(nel)

  allocate(eltypep(nelpmax,0:nproc-1),coarsingp(nelpmax,0:nproc-1))
  do iproc = 0, nproc-1
   nelp = nel(iproc)
   do iel = 1, nelp
    ielg = procel(iel,iproc)
    eltypep(iel,iproc) = eltypeg(ielg)
    coarsingp(iel,iproc)=coarsing(ielg)
   end do
  end do

! Fluid
  nelpmax=maxval(nel_fluid)
  allocate(eltypep_fluid(nelpmax,0:nproc-1))
  do iproc = 0, nproc-1
   do iel = 1, nel_fluid(iproc)
    ielg = procel_fluid(iel,iproc)
    eltypep_fluid(iel,iproc) = eltypeg(ielg)
   end do
  end do

! Solid
  nelpmax=maxval(nel_solid)
  allocate(eltypep_solid(nelpmax,0:nproc-1))
  do iproc = 0, nproc-1
   do iel = 1, nel_solid(iproc)
    ielg = procel_solid(iel,iproc)
    eltypep_solid(iel,iproc) = eltypeg(ielg)
   end do
  end do

  end subroutine define_element_type
!--------------------------------------------------------------------------
!
!dk out_dd-----------------------------------------------------------------
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
22 format('>')
  end subroutine out_dd
!------------------------------------------------------------------------ 

!dk write_db-------------------------------------------------------------
subroutine write_db
!
! Writes out a database file to be read by the solver for each processor.
!
!
use data_gllmesh

integer :: iproc,iptp,npointotp,ipsrc,ipdes,imsg,iptcp,iel,inode,ielg,idom,i
character(len=4) :: appiproc
character(len=80) :: dbname
integer :: lfdbname

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do iproc=0, nproc-1
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   call define_io_appendix(appiproc,iproc)
   dbname='meshdb.dat'//appiproc
   lfdbname=index(dbname,' ')-1
   if (dump_mesh_info_screen) write(6,*) 'WRITING OUT DATABASE TO ',dbname(1:lfdbname)
   open(10,file=dbname(1:lfdbname), FORM="UNFORMATTED", STATUS="REPLACE")

    npointotp = nel(iproc)*(npol+1)**2 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Coordinates of control points
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing coordinates/control points...',iproc

! global number of control nodes (slightly differs for each processor!)
   write(10) nglobmeshp(iproc)

   do iptcp = 1, nglobmeshp(iproc) 
      write(10) iptcp, router*scpp(iptcp,iproc),router*zcpp(iptcp,iproc)
   end do

!  Topology of control points
   do iel = 1, nel(iproc)
    ielg = procel(iel,iproc)
!   if (zgll(npol/2,npol/2,ielg)>= 0.) then ! NORTH
! af
    if (zcom(ielg)>= 0.) then ! NORTH
       write(10) (lnodescp(iel,inode,iproc),inode=1,8)
   else ! SOUTH
       write(10) (lnodescp(iel,inode,iproc),inode=7,1,-1),lnodescp(iel,8,iproc)
    endif
   end do

! Number of global distinct points (slightly differs for each processor!)
    write(10) nglobp(iproc)

!  Element types
   do iel = 1, nel(iproc)
    write(10) eltypep(iel,iproc),coarsingp(iel,iproc)
   end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Solid-Fluid distinction
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing solid/fluid domain info...',iproc

! mapping from sol/flu to global element numbers
    do iel=1,nel_solid(iproc)
       write(10) procel_solidp(iel,iproc)
    enddo

    do iel=1,nel_fluid(iproc)
       write(10) procel_fluidp(iel,iproc)
    enddo

! Number of distinct points in solid (slightly differs for each processor!)
    write(10) nglobp_solid(iproc)

!  slocal numbering 
    npointotp = nel_solid(iproc)*(npol+1)**2 
    do iptp = 1, npointotp
       write(10) igloc_solid(iptp,iproc)
    end do

!  flocal numbering
    npointotp = nel_fluid(iproc)*(npol+1)**2 
    do iptp = 1, npointotp
       write(10) igloc_fluid(iptp,iproc)
    end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Solid-Fluid boundary
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing solid/fluid boundary info...',iproc
   write(10)have_bdry_elemp(iproc)
   do iel=1,nbdry_el(iproc)
      write(10)bdry_solid_elp(iel,iproc),bdry_fluid_elp(iel,iproc)
   enddo
   do iel=1,nbdry_el(iproc)
      write(10)bdry_jpol_solidp(iel,iproc),bdry_jpol_fluidp(iel,iproc)
   enddo

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! General numerical input/output parameters
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing numerical parameters...',iproc
   write(10)pts_wavelngth,period,courant,dt

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Background model
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing background model info...',iproc
   write(10) bkgrdmodel(1:lfbkgrdmodel)
   write(10) router,resolve_inner_shear,have_fluid
   do idom=1,ndisc
      write(10) discont(idom),solid_domain(idom),idom_fluid(idom)
   enddo
   write(10)rmin,minh_ic,maxh_ic,maxh_icb

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Min/max grid spacing
   write(10)hmin_glob,hmax_glob
   write(10)min_distance_dim,min_distance_nondim

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! critical ratios h/v min/max and locations
   write(10)char_time_max,char_time_max_globel
   write(10)char_time_max_rad,char_time_max_theta
   write(10)char_time_min,char_time_min_globel
   write(10)char_time_min_rad,char_time_min_theta

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Axial element arrays
   write(10)have_axis(iproc)
   write(10)naxelp(iproc),naxel_solidp(iproc),naxel_fluidp(iproc)

   do i=1,naxelp(iproc)
      write(10)ax_elp(i,iproc)
   enddo
   do i=1,naxel_solidp(iproc)
      write(10)ax_el_solidp(i,iproc)
   enddo
   do i=1,naxel_fluidp(iproc)
      write(10)ax_el_fluidp(i,iproc)
   enddo

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (dump_mesh_info_screen) write(6,*)'PARALLEL DATABASE: writing communication info...',iproc

! SSSSSSSSSSSS SOLID MESSAGING SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!  number of processors to receive data from
   write(10) sizerecvp_solid(iproc)

   if (sizerecvp_solid(iproc) > 0 ) then
!     list of processors to receive data from
      write(10) listrecvp_solid(1:sizerecvp_solid(iproc),iproc)

!     size of messages coming from each of these processors
      write(10) sizemsgrecvp_solid(1:sizerecvp_solid(iproc),iproc)

!     glocal_index corresponding to each of these messages
      do imsg = 1, sizerecvp_solid(iproc)
         ipsrc = listrecvp_solid(imsg,iproc) ! To be deleted
         do iptp = 1, sizemsgrecvp_solid(imsg,iproc)
! TO DOUBLE CHECK
            write(10) glocal_index_msg_recvp_solid(iptp,imsg,iproc)
         end do
      end do
   end if

!  number of processors to send data to
   write(10) sizesendp_solid(iproc)

!  list of processors to send data to
   if (sizesendp_solid(iproc) > 0 ) then
      write(10) listsendp_solid(1:sizesendp_solid(iproc),iproc)

!     size of messages going to each of these processors
      write(10) sizemsgsendp_solid(1:sizesendp_solid(iproc),iproc)

!     glocal_index corresponding to each of these messages
      do imsg = 1, sizesendp_solid(iproc)
         ipdes = listsendp_solid(imsg,iproc) ! to be deleted
         do iptp = 1, sizemsgsendp_solid(imsg,iproc)
!!!      TO DOUBLE CHECK
            write(10) glocal_index_msg_sendp_solid(iptp,imsg,iproc)
         end do
      end do
   end if

! FFFFFFFFFFFFF FLUID MESSAGING FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Only written if there is a fluid!
   if (have_fluid) then
      if (dump_mesh_info_screen) write(6,*)'writing fluid messaging'
!     number of processors to receive data from
      write(10) sizerecvp_fluid(iproc)

      if (sizerecvp_fluid(iproc) > 0 ) then
!        list of processors to receive data from
         write(10) listrecvp_fluid(1:sizerecvp_fluid(iproc),iproc)

!        size of messages coming from each of these processors
         write(10) sizemsgrecvp_fluid(1:sizerecvp_fluid(iproc),iproc)

!        glocal_index corresponding to each of these messages
         do imsg = 1, sizerecvp_fluid(iproc)
            ipsrc = listrecvp_fluid(imsg,iproc) ! To be deleted
            do iptp = 1, sizemsgrecvp_fluid(imsg,iproc)
            ! TO DOUBLE CHECK
               write(10) glocal_index_msg_recvp_fluid(iptp,imsg,iproc)
            end do
         end do
      end if

!     number of processors to send data to
      write(10) sizesendp_fluid(iproc)

!     list of processors to send data to
      if (sizesendp_fluid(iproc) > 0 ) then
         write(10) listsendp_fluid(1:sizesendp_fluid(iproc),iproc)

!        size of messages going to each of these processors
         write(10) sizemsgsendp_fluid(1:sizesendp_fluid(iproc),iproc)

!        glocal_index corresponding to each of these messages
         do imsg = 1, sizesendp_fluid(iproc)
            ipdes = listsendp_fluid(imsg,iproc) ! to be deleted
            do iptp = 1, sizemsgsendp_fluid(imsg,iproc)
!!!      TO DOUBLE CHECK
               write(10) glocal_index_msg_sendp_fluid(iptp,imsg,iproc)
            end do
         end do
      end if
   endif ! have_fluid

   close(10)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   write(6,*)'....Wrote database for processor',iproc; call flush(6)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  end do
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine write_db
!------------------------------------------------------------------------

!dk read_db--------------------------------------------------------------
  subroutine read_db
    
    use data_gllmesh

  integer :: iproc,iptp,npointotp,ipsrc,ipdes,imsg,ipt,inode,iptcp,iel

  character(len=4) :: appiproc
  character(len=80) :: dbname
  integer :: lfdbname

  integer :: nelem,npoint,npoin

  integer :: tmpsizemsgrecvmax_solid,tmpsizemsgsendmax_solid
  integer :: tmpsizerecv_solid,tmpsizesend_solid
  integer, dimension(:), allocatable :: tmplistrecv_solid,tmpsizemsgrecv_solid
  integer, dimension(:), allocatable :: tmplistsend_solid,tmpsizemsgsend_solid
  integer, dimension(:,:), allocatable :: tmpglocal_index_msg_recv_solid
  integer, dimension(:,:), allocatable :: tmpglocal_index_msg_send_solid

  integer :: tmpsizemsgrecvmax_fluid,tmpsizemsgsendmax_fluid
  integer :: tmpsizerecv_fluid,tmpsizesend_fluid
  integer, dimension(:), allocatable :: tmplistrecv_fluid,tmpsizemsgrecv_fluid
  integer, dimension(:), allocatable :: tmplistsend_fluid,tmpsizemsgsend_fluid
  integer, dimension(:,:), allocatable :: tmpglocal_index_msg_recv_fluid
  integer, dimension(:,:), allocatable :: tmpglocal_index_msg_send_fluid

!  real, dimension(:,:), allocatable :: crd_nodes
  double precision, dimension(:,:), allocatable :: crd_nodes
  integer, dimension(:,:), allocatable :: lnods
  character(len=6), dimension(:), allocatable :: eltype
  logical, dimension(:), allocatable          :: coarsingtmp
  integer :: tmpnglob,tmpnglob_solid

! Solid fluid variables
    integer tmpnel_solid,tmpnel_fluid
    integer, allocatable :: tmpigloc_solid(:)
    integer, allocatable :: tmpigloc_fluid(:)
    integer :: tmpnbdry_el, tmphave_bdry_elem
    integer, allocatable :: tmpbdry_solid_elp(:),tmpbdry_fluid_elp(:)
    integer, allocatable :: tmpbdry_jpol_solidp(:),tmpbdry_jpol_fluidp(:)
    integer, allocatable :: tmpprocel_solid(:),tmpprocel_fluid(:)
  
! General numerical input/output parameters
    integer :: tmppts_wavelngth
    real :: tmpperiod,tmpcourant,tmpdt

! Background model
   character(len=100) :: tmpbkgrdmodel
   double precision :: tmprouter
   logical :: tmpresolve_inner_shear,tmphave_fluid
   integer :: tmpndisc,idom,i
   double precision, allocatable ::   tmpdiscont(:)
   logical, allocatable :: tmpsolid_domain(:)
   integer, allocatable :: tmpidom_fluid(:)
   double precision :: tmprmin,tmpminh_ic,tmpmaxh_ic,tmpmaxh_icb

! Min/max grid spacing
   double precision :: tmphmin_glob,tmphmax_glob
   double precision :: tmpmin_distance_dim,tmpmin_distance_nondim

! critical ratios h/v min/max and locations
   integer :: tmpchar_time_max_globel
   double precision :: tmpchar_time_max
   double precision :: tmpchar_time_max_rad,tmpchar_time_max_theta
   integer :: tmpchar_time_min_globel
   double precision :: tmpchar_time_min
   double precision :: tmpchar_time_min_rad,tmpchar_time_min_theta

! Axial element arrays
   logical :: tmphave_axis
   integer :: tmpnaxelp,tmpnaxel_solidp,tmpnaxel_fluidp
   integer, allocatable :: tmpax_elp(:)
   integer, allocatable :: tmpax_el_solidp(:)
   integer, allocatable :: tmpax_el_fluidp(:)

   write(6,*)'TEST: read in new database...'
  do iproc=0, nproc-1
   call define_io_appendix(appiproc,iproc)
   dbname='meshdb.dat'//appiproc
   lfdbname=index(dbname,' ')-1
   if (dump_mesh_info_screen) &
    write(6,*) iproc,'READING IN DATABASE FROM ',dbname(1:lfdbname)
   open(10,file=dbname(1:lfdbname), FORM="UNFORMATTED", STATUS="OLD")

   nelem=nel(iproc)
   npoint = nelem*(npol+1)**2

   read(10) npoin

   allocate(crd_nodes(npoin,2))

   do iptcp = 1, npoin 
     read(10) ipt, crd_nodes(iptcp,1),crd_nodes(iptcp,2)
   end do

   allocate(lnods(nelem,1:8))
   do iel = 1, nelem
    read(10) (lnods(iel,inode), inode=1,8)
   end do

! Number of global distinct points (slightly differs for each processor!)
   read(10) tmpnglob

!  Element type
   allocate(eltype(nelem),coarsingtmp(nelem))
   do iel = 1, nelem
    read(10) eltype(iel),coarsingtmp(iel)
   end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Solid-Fluid distinction

!  number of elements for proc iproc
!    read(10) tmpnel_solid,tmpnel_fluid
   tmpnel_solid=nel_solid(iproc)
   tmpnel_fluid=nel_fluid(iproc)

   if (dump_mesh_info_screen) write(6,*)'READ nel_solid:',tmpnel_solid

! mapping from sol/flu to global element numbers
    allocate(tmpprocel_solid(tmpnel_solid))
    allocate(tmpprocel_fluid(tmpnel_fluid))
    do iel=1,tmpnel_solid
       read(10) tmpprocel_solid(iel)
    enddo
    do iel=1,tmpnel_fluid
       read(10) tmpprocel_fluid(iel)
    enddo

! Number of distinct points in solid (slightly differs for each processor!)
    read(10) tmpnglob_solid
   if (dump_mesh_info_screen) write(6,*)'READ nglob_solid:',tmpnglob_solid

!  glocal number of points for proc iproc
!    read(10) tmpnglobp_solid,tmpnglobp_fluid

    npointotp = tmpnel_solid*(npol+1)**2 
    allocate(tmpigloc_solid(npointotp))
!  slocal numbering 
    do iptp = 1, npointotp
     read(10) tmpigloc_solid(iptp)
    end do    

    npointotp = tmpnel_fluid*(npol+1)**2 
    allocate(tmpigloc_fluid(npointotp))
!  flocal numbering 
    do iptp = 1, npointotp
     read(10) tmpigloc_fluid(iptp)
    end do    

! Assuming here that we don't need the solid & fluid 
! control node definitions, therefore left out

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Solid-Fluid boundary
!   read(10)tmpnbdry_el
   tmpnbdry_el=nbdry_el(iproc)
   read(10)tmphave_bdry_elem
   allocate(tmpbdry_solid_elp(tmpnbdry_el))
   allocate(tmpbdry_fluid_elp(tmpnbdry_el))
   do iel=1,tmpnbdry_el
      read(10)tmpbdry_solid_elp(iel),tmpbdry_fluid_elp(iel)
   enddo

   allocate(tmpbdry_jpol_solidp(tmpnbdry_el))
   allocate(tmpbdry_jpol_fluidp(tmpnbdry_el))
   do iel=1,tmpnbdry_el
      read(10)tmpbdry_jpol_solidp(iel),tmpbdry_jpol_fluidp(iel)
   enddo

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! General numerical input/output parameters
!   read(10) tmpnpol,tmppts_wavelngth
   read(10) tmppts_wavelngth,tmpperiod,tmpcourant,tmpdt

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Background model
!   read(10) tmplfbkgrdmodel
   read(10) tmpbkgrdmodel(1:lfbkgrdmodel)
 ! af
   read(10) tmprouter,tmpresolve_inner_shear,tmphave_fluid
!   read(10) tmpndisc
   tmpndisc=ndisc
   allocate(tmpdiscont(tmpndisc),tmpsolid_domain(tmpndisc)) 
   allocate(tmpidom_fluid(tmpndisc))
   do idom=1,tmpndisc
      read(10) tmpdiscont(idom),tmpsolid_domain(idom),tmpidom_fluid(idom)
   enddo
!   read(10)tmpnz_glob,tmpns_glob,tmpnc_glob
!   allocate(tmpdz_glob(tmpnz_glob))
!   do i=1,tmpnz_glob
!      read(10)tmpdz_glob(i)
!   enddo
   read(10)tmprmin,tmpminh_ic,tmpmaxh_ic,tmpmaxh_icb

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Min/max grid spacing
   read(10)tmphmin_glob,tmphmax_glob
   read(10)tmpmin_distance_dim,tmpmin_distance_nondim

! critical ratios h/v min/max and locations
   read(10)tmpchar_time_max,tmpchar_time_max_globel
   read(10)tmpchar_time_max_rad,tmpchar_time_max_theta
   read(10)tmpchar_time_min,tmpchar_time_min_globel
   read(10)tmpchar_time_min_rad,tmpchar_time_min_theta

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Axial element arrays
   read(10)tmphave_axis
   read(10)tmpnaxelp,tmpnaxel_solidp,tmpnaxel_fluidp

   allocate(tmpax_elp(tmpnaxelp),tmpax_el_solidp(tmpnaxel_solidp))
   allocate(tmpax_el_fluidp(tmpnaxel_fluidp))
   do i=1,tmpnaxelp
      read(10)tmpax_elp(i)
   enddo
   do i=1,tmpnaxel_solidp
      read(10)tmpax_el_solidp(i)
   enddo
   do i=1,tmpnaxel_fluidp
      read(10)tmpax_el_fluidp(i)
   enddo

   read(10) tmpsizerecv_solid
   if (dump_mesh_info_screen) write(6,*)'READ: size solid recv:',tmpsizerecv_solid; call flush(6)
   if ( tmpsizerecv_solid > 0) then 
    allocate(tmplistrecv_solid(tmpsizerecv_solid)) ; tmplistrecv_solid(:) = -1
    allocate(tmpsizemsgrecv_solid(tmpsizerecv_solid)) ; tmpsizemsgrecv_solid(:) = 0
    if (dump_mesh_info_screen) then
     write(6,*)'read listrecv_solid'; call flush(6)
    end if
    read(10) tmplistrecv_solid(:)
    if (dump_mesh_info_screen) then
     write(6,*)'read sizemsgrecv_solid'; call flush(6)
    end if
    read(10) tmpsizemsgrecv_solid(:)
    tmpsizemsgrecvmax_solid = maxval(tmpsizemsgrecv_solid(:))
    if (dump_mesh_info_screen) then
     write(6,*)'READ: sizemsgrecvmax_solid:',tmpsizemsgrecvmax_solid; call flush(6)
    end if
    allocate(tmpglocal_index_msg_recv_solid(tmpsizemsgrecvmax_solid,tmpsizerecv_solid))
    tmpglocal_index_msg_recv_solid(:,:) = 0
    do imsg = 1, tmpsizerecv_solid
!      write(6,*)'imsg=',imsg; call flush(6)
     ipsrc = tmplistrecv_solid(imsg) 
     do iptp = 1, tmpsizemsgrecv_solid(imsg)
!      write(6,*)'iptp=',iptp; call flush(6)
       read(10) tmpglocal_index_msg_recv_solid(iptp,imsg)
     end do
    end do
   end if
   read(10) tmpsizesend_solid
   if (dump_mesh_info_screen) write(6,*)'READ: size solid send:',tmpsizesend_solid
   if ( tmpsizesend_solid > 0) then 
    allocate(tmplistsend_solid(tmpsizesend_solid)) ; tmplistsend_solid(:) = -1
    allocate(tmpsizemsgsend_solid(tmpsizesend_solid)) ; tmpsizemsgsend_solid(:) = 0
    read(10) tmplistsend_solid(:)
    read(10) tmpsizemsgsend_solid(:) 
    tmpsizemsgsendmax_solid = maxval(tmpsizemsgsend_solid(:))
    allocate(tmpglocal_index_msg_send_solid(tmpsizemsgsendmax_solid,tmpsizesend_solid))
    tmpglocal_index_msg_send_solid(:,:) = 0
    do imsg = 1, tmpsizesend_solid
     ipdes = tmplistsend_solid(imsg)
      do iptp = 1, tmpsizemsgsend_solid(imsg)
       read(10) tmpglocal_index_msg_send_solid(iptp,imsg)
      end do
    end do
   end if

   if (tmphave_fluid) then
      if (dump_mesh_info_screen) write(6,*)'reading fluid messaging'
      read(10) tmpsizerecv_fluid
      if (dump_mesh_info_screen) write(6,*)'READ: size fluid recv:',tmpsizerecv_fluid
      if ( tmpsizerecv_fluid > 0) then 
         allocate(tmplistrecv_fluid(tmpsizerecv_fluid)) ; tmplistrecv_fluid(:) = -1
         allocate(tmpsizemsgrecv_fluid(tmpsizerecv_fluid)) ; tmpsizemsgrecv_fluid(:) = 0
         read(10) tmplistrecv_fluid(:)
         read(10) tmpsizemsgrecv_fluid(:) 
         tmpsizemsgrecvmax_fluid = maxval(tmpsizemsgrecv_fluid(:))
         allocate(tmpglocal_index_msg_recv_fluid&
              (tmpsizemsgrecvmax_fluid,tmpsizerecv_fluid))
         tmpglocal_index_msg_recv_fluid(:,:) = 0
         do imsg = 1, tmpsizerecv_fluid
            ipsrc = tmplistrecv_fluid(imsg) 
            do iptp = 1, tmpsizemsgrecv_fluid(imsg)
               read(10) tmpglocal_index_msg_recv_fluid(iptp,imsg)
            end do
         end do
      end if
      read(10) tmpsizesend_fluid
      if (dump_mesh_info_screen) write(6,*)'READ: size fluid send:',tmpsizesend_fluid
      if ( tmpsizesend_fluid > 0) then 
         allocate(tmplistsend_fluid(tmpsizesend_fluid)) ; tmplistsend_fluid(:) = -1
         allocate(tmpsizemsgsend_fluid(tmpsizesend_fluid)) ; tmpsizemsgsend_fluid(:) = 0
         read(10) tmplistsend_fluid(:)
         read(10) tmpsizemsgsend_fluid(:) 
         tmpsizemsgsendmax_fluid = maxval(tmpsizemsgsend_fluid(:))
         allocate(tmpglocal_index_msg_send_fluid&
                  (tmpsizemsgsendmax_fluid,tmpsizesend_fluid))
         tmpglocal_index_msg_send_fluid(:,:) = 0
         do imsg = 1, tmpsizesend_fluid
            ipdes = tmplistsend_fluid(imsg)
            do iptp = 1, tmpsizemsgsend_fluid(imsg)
               read(10) tmpglocal_index_msg_send_fluid(iptp,imsg)
            end do
         end do
      end if
   endif ! have_fluid
   if (dump_mesh_info_screen) write(6,*)iproc,'read its database'; call flush(6)

   close(10)

!  Deallocation phase for serial parallel test
   deallocate(eltype,coarsingtmp)
   deallocate(lnods)
   deallocate(crd_nodes)

! Solid fluid
    deallocate(tmpprocel_solid,tmpprocel_fluid)
    deallocate(tmpigloc_solid)
    deallocate(tmpigloc_fluid)
    deallocate(tmpbdry_solid_elp,tmpbdry_fluid_elp)
    deallocate(tmpbdry_jpol_solidp,tmpbdry_jpol_fluidp)

    deallocate(tmpdiscont,tmpsolid_domain) 
    deallocate(tmpidom_fluid)

   deallocate(tmpax_elp,tmpax_el_solidp)
   deallocate(tmpax_el_fluidp)

   if (allocated(tmpglocal_index_msg_send_solid)) &
        deallocate(tmpglocal_index_msg_send_solid)
   if (allocated(tmpsizemsgsend_solid)) deallocate (tmpsizemsgsend_solid)
   if (allocated(tmplistsend_solid)) deallocate (tmplistsend_solid)

   if (allocated(tmpglocal_index_msg_recv_solid)) &
        deallocate(tmpglocal_index_msg_recv_solid)
   if (allocated(tmpsizemsgrecv_solid)) deallocate (tmpsizemsgrecv_solid)
   if (allocated(tmplistrecv_solid)) deallocate (tmplistrecv_solid)

   if (allocated(tmpglocal_index_msg_send_fluid)) &
        deallocate(tmpglocal_index_msg_send_fluid)
   if (allocated(tmpsizemsgsend_fluid)) deallocate (tmpsizemsgsend_fluid)
   if (allocated(tmplistsend_fluid)) deallocate (tmplistsend_fluid)

   if (allocated(tmpglocal_index_msg_recv_fluid)) &
        deallocate(tmpglocal_index_msg_recv_fluid)
   if (allocated(tmpsizemsgrecv_fluid)) deallocate (tmpsizemsgrecv_fluid)
   if (allocated(tmplistrecv_fluid)) deallocate (tmplistrecv_fluid)

  end do
  write(6,*)'TEST: Successful.'; call flush(6)

  end subroutine read_db
!------------------------------------------------------------------------

!dk create_static_header------------------------------------------------
  subroutine create_static_header

    integer :: iproc
    character(len=8) :: mydate
    character(len=10) :: mytime
    character(len=4) :: appiproc
    character(len=80) :: dbname,dbname2
    integer :: lfdbname

    call date_and_time(mydate,mytime)
    do iproc=0,nproc-1
    call define_io_appendix(appiproc,iproc)
    dbname='mesh_params.h'//appiproc
    dbname2='mesh_params.h'
    lfdbname=index(dbname,' ')-1
    open(98,file=dbname(1:lfdbname), STATUS="REPLACE")
    write(98,10) nproc
    write(98,11) mydate(5:6),mydate(7:8),mydate(1:4),mytime(1:2),mytime(3:4)
    write(98,*)''
    write(98,29)
    write(98,12)'Background model     :',bkgrdmodel(1:lfbkgrdmodel) ! af
    write(98,13)'Inner-core shear wave:',resolve_inner_shear
    write(98,14)'Dominant period [s]  :',period
    write(98,14)'Elements/wavelength  :',pts_wavelngth
    write(98,14)'Courant number       :',courant
    write(98,15)'Coarsening levels    :',nc_init
    write(98,30)
    write(98,*)''
    write(98,9)'npol',npol,'polynomial order'
    write(98,9)'nelem',nel(iproc),'proc. els'
    write(98,9)'npoint',nel(iproc)*(npol+1)**2,'proc. all pts'
    write(98,9)'nel_solid',nel_solid(iproc),'proc. solid els'
    write(98,9)'nel_fluid',nel_fluid(iproc),'proc. fluid els'
    write(98,9)'npoint_solid',nel_solid(iproc)*(npol+1)**2,'proc. solid pts'
    write(98,9)'npoint_fluid',nel_fluid(iproc)*(npol+1)**2,'proc. fluid pts'
    write(98,9)'nglob',nglobp(iproc),'proc. glocal pts'
    write(98,9)'nglob_solid',nglobp_solid(iproc),'proc. slocal pts'
    write(98,9)'nglob_fluid',nglobp_fluid(iproc),'proc. flocal pts'
    write(98,9)'nel_bdry',nbdry_el(iproc),'proc. solid-fluid bndry els'
    write(98,9)'npoin',nglobmeshp(iproc),'proc. control nodes'
    write(98,9)'ndisc',ndisc,'# disconts in bkgrd model'
    write(98,9)'nproc_mesh',nproc,'number of processors'
    write(98,9)'lfbkgrdmodel',lfbkgrdmodel,'length of bkgrdmodel name'
    write(98,*)''
    write(98,31)
    write(98,14)'Time step [s]        :',dt
    write(98,16)'Min(h/vp),dt/courant :',minhvp,dt/courant*real(npol)
    write(98,16)'max(h/vs),T0/wvlngth :',maxhvs,period/pts_wavelngth
    write(98,14)'Inner core r_min [km]:',rmin/1.d3
    write(98,14)'Max(h) r/ns(icb) [km]:',maxhnsicb/1.d3
    write(98,14)'Max(h) precalc.  [km]:',maxh_icb/1.d3
    write(98,30)
    write(98,*)''
    close(98)
    write(6,*)'wrote parameters for static solver into ',dbname(1:lfdbname)
    enddo

    lfdbname=index(dbname2,' ')-1
    open(97,file=dbname2(1:lfdbname), STATUS="REPLACE")
    iproc=0
    write(97,10) nproc
    write(97,11) mydate(5:6),mydate(7:8),mydate(1:4),mytime(1:2),mytime(3:4)
    write(97,*)''
    write(97,29)
    write(97,12)'Background model     :',bkgrdmodel(1:lfbkgrdmodel)
    write(97,13)'Inner-core shear wave:',resolve_inner_shear
    write(97,14)'Dominant period [s]  :',period
    write(97,14)'Elements/wavelength  :',pts_wavelngth
    write(97,14)'Courant number       :',courant
    write(97,15)'Coarsening levels    :',nc_init
    write(97,30)
    write(97,*)''
    write(97,9)'npol',npol,'polynomial order'
    write(97,9)'nelem',nel(iproc),'proc. els'
    write(97,9)'npoint',nel(iproc)*(npol+1)**2,'proc. all pts'
    write(97,9)'nel_solid',nel_solid(iproc),'proc. solid els'
    write(97,9)'nel_fluid',nel_fluid(iproc),'proc. fluid els'
    write(97,9)'npoint_solid',nel_solid(iproc)*(npol+1)**2,'proc. solid pts'
    write(97,9)'npoint_fluid',nel_fluid(iproc)*(npol+1)**2,'proc. fluid pts'
!    write(97,9)'nglob',nglobp(iproc),'proc. glocal pts'
!    write(97,9)'nglob_solid',nglobp_solid(iproc),'proc. slocal pts'
    write(97,9)'nglob_fluid',nglobp_fluid(iproc),'proc. flocal pts'
    write(97,9)'nel_bdry',nbdry_el(iproc),'proc. solid-fluid bndry els'
!    write(97,9)'npoin',nglobmeshp(iproc),'proc. control nodes'
    write(97,9)'ndisc',ndisc,'# disconts in bkgrd model'
    write(97,9)'nproc_mesh',nproc,'number of processors'
    write(97,9)'lfbkgrdmodel',lfbkgrdmodel,'length of bkgrdmodel name'

    write(97,*)''
    write(97,31)
    write(97,14)'Time step [s]        :',dt
    write(97,16)'Min(h/vp),dt/courant :',minhvp,dt/courant*real(npol)
    write(97,16)'max(h/vs),T0/wvlngth :',maxhvs,period/pts_wavelngth
    write(97,14)'Inner core r_min [km]:',rmin/1.d3
    write(97,14)'Max(h) r/ns(icb) [km]:',maxhnsicb/1.d3
    write(97,14)'Max(h) precalc.  [km]:',maxh_icb/1.d3
    write(97,30)
    write(97,*)''
    close(97)
    write(6,*)'wrote parameters for static solver into ',dbname2(1:lfdbname)

9 format(' integer, parameter :: ',A12,' =',i10,'  ! ',A27)
10 format('! Proc ',i3,': Header for mesh information to run static solver')
11 format('! created by the mesher on ', &
            A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
29 format('!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::')
12 format('!  ',A23,A20)
13 format('!  ',A23,L10)
14 format('!  ',A23,1f10.4)
16 format('!  ',A23,2(f10.4))
15 format('!  ',A23,I10)
30 format('!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

31 format('!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::')


  end subroutine create_static_header
!------------------------------------------------------------------------

!dk create_unrolled_loop_routines----------------------------------------------
subroutine create_unrolled_loop_routines

    character(len=8) :: mydate
    character(len=10) :: mytime
    integer :: ipol

! Creates the unrolled-loop routines for the solver that specifically 
! depend upon the polynomial order npol
    
call date_and_time(mydate,mytime)

open(unit=101,file='unrolled_loops.f90')
write(101,*)'!========================'
write(101,*)'module unrolled_loops'
write(101,*)'!========================'
write(101,*)
write(101,*)'! This module has been created by the mesher'
write(101,5) npol
!write(101,11) mydate(5:6),mydate(7:8),mydate(1:4),mytime(1:2),mytime(3:4)
write(101,*)
write(101,*)'  use global_parameters'
write(101,*)
write(101,*)'  implicit none'
write(101,*)
write(101,*)'  public :: mxm,vxm'
write(101,*)'  private'
write(101,*)
write(101,*)'  contains' 
write(101,*)
write(101,*)'!============================================================'
write(101,*)
write(101,*)'subroutine mxm(a,b,c)'
write(101,*)
write(101,*)'  include "mesh_params.h" '
write(101,*)
write(101,9)  npol,npol,npol,npol
write(101,10) npol,npol
write(101,*)'  integer i,j'
write(101,*)
write(101,6) npol
write(101,61)
write(101,62)
write(101,63) npol
write(101,64)
write(101,65)
write(101,*)
write(101,7) npol
write(101,8) npol
write(101,*)'      c(i,j) = & '
do ipol=0,npol-1
   write(101,12)ipol,ipol
enddo
write(101,13)npol,npol
write(101,*)'    end do'
write(101,*)'  end do '
write(101,*)'  return'
write(101,*)
write(101,*)'end subroutine mxm'
write(101,*)'!-------------------------------------------------------------'
write(101,*)
write(101,*)
write(101,*)'!-------------------------------------------------------------'
write(101,*)'subroutine vxm(a,b,c)'
write(101,*)
write(101,*)'  include "mesh_params.h" '
write(101,*)
write(101,16) npol,npol,npol
write(101,17) npol
write(101,*)'  integer j'
write(101,*)
write(101,6) npol
write(101,61)
write(101,62)
write(101,63) npol
write(101,64)
write(101,65)
write(101,*)
write(101,7) npol
write(101,*)'      c(j) = & '
do ipol=0,npol-1
   write(101,14)ipol,ipol
enddo
write(101,15)npol,npol
write(101,*)'  end do '
write(101,*)'  return'
write(101,*)
write(101,*)'end subroutine vxm'
write(101,*)'!-------------------------------------------------------------'
write(101,*)
write(101,*)'!========================'
write(101,*)'end module unrolled_loops'
write(101,*)'!========================'
close(101)

write(6,*)'wrote module unrolled_loops into unrolled_loops.f90'

5   format(' ! for polynomial order ',i2)
6   format('   if ( npol /= ',i2,' ) then')
61  format('      write(6,*)"Problem: unrolled_loops.f90 has different" ')
62  format('      write(6,*)"         polynomial order than mesh_params.h:" ')
63  format('      write(6,*)"        ',i2,'",npol')
64  format('      stop')
65  format('   endif')
7   format('   do j = 0,',i2)
8   format('     do i = 0,',i2)
9   format('   real(kind=realkind), intent(in)  :: a(0:',i2, &
           ',0:',i2,'),b(0:',i2,',0:',i2,')')
10  format('   real(kind=realkind), intent(out) :: c(0:',i2,',0:',i2,')')
11  format(' ! on ', A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
12  format('          + a(i,',i2,') * b(',i2,',j) &')
13  format('          + a(i,',i2,') * b(',i2,',j)')
14  format('          + a(',i2,') * b(',i2,',j) &')
15  format('          + a(',i2,') * b(',i2,',j)')
16  format('   real(kind=realkind), intent(in)  :: a(0:',i2, &
           '),b(0:',i2,',0:',i2,')')
17  format('   real(kind=realkind), intent(out) :: c(0:',i2,')')

end subroutine create_unrolled_loop_routines
!------------------------------------------------------------------------


!dk define_io_appendix---------------------------------------------------
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
!------------------------------------------------------------------------

!////////////////////////////////////////////////////////////////////////
!=================
  end module pdb
!=================
