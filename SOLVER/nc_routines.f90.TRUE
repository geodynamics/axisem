module nc_routines
contains
!!!!!!!!!!AUXUILIAR----------------------------------------NETCDF----------!!!!!!!
!subroutine dump_smtg_ncdf(valore,ishot,n,filename)
!implicit none
!character(len=55), intent(in)              :: filename
!integer,intent(in) :: n,ishot
!real(kind=realkind),dimension(n),intent(in) :: valore

!!write(6,*)'computing nc file ',trim(filename),' ...'

!if(ishot==1) then
! 		call init_ncdf(ishot,trim(filename),real(valore,4),n)
!else
! 		call write_ncdf(ishot,trim(filename),real(valore,4),n)
!endif
!end subroutine dump_smtg_ncdf

!!!!!!!!!!AUXUILIAR----------------------------------------NETCDF----------!!!!!!!
!subroutine init_ncdf(ishot,filename,u,n)
!use netcdf
!implicit none

!integer*4 :: ishot,i
!character (len = *) :: filename
!integer, parameter :: NDIMS = 2;
!integer :: n;
!character (len = *), parameter :: U_NAME = "data"
!character (len = *), parameter :: REC_NAME = "time"
!character (len = *), parameter :: UNITS = "index"
!integer :: ncid, varid, dimids(NDIMS)
!integer :: u_dimid,t_dimid
!integer :: start(NDIMS), count(NDIMS)
!real*4, dimension(n):: u

! call check(nf90_create(filename,NF90_CLOBBER, ncid))
! call check(nf90_def_dim(ncid,UNITS,n,u_dimid))
! call check(nf90_def_dim(ncid,REC_NAME,NF90_UNLIMITED,t_dimid))
! dimids=(/u_dimid, t_dimid/)
! call check( nf90_def_var(ncid,U_NAME,NF90_REAL,dimids,varid))
! call check( nf90_enddef(ncid))
! count = (/n,1/)
! start = (/1,1/)
! start(2)=ishot;
! call check(nf90_put_var(ncid,varid,real(u,4),start=start,count=count));
! call check(nf90_close(ncid)) 

!contains
!  subroutine check(status)
!    integer, intent ( in) :: status
!    
!    if(status /= nf90_noerr) then 
!      print *, trim(nf90_strerror(status))
!      stop 2
!    end if
!  end subroutine check  
!end subroutine init_ncdf

!!!!---------------------------------------------------------------------!!!!!
!subroutine write_ncdf(ishot,filename,u,n)
!use NETCDF
!implicit none

!integer*4 :: ishot,i
!character (len = *):: filename
!character (len = *), parameter :: U_NAME = "data"
!integer, parameter :: NDIMS = 2;
!integer :: n,ncid,varid;
!integer :: start(NDIMS), count(NDIMS)
!real*4, dimension(n)::u

! call check(nf90_open(filename,NF90_WRITE, ncid))
!	call check(nf90_inq_varid(ncid,U_NAME,varid))
! count = (/n,1/)
! start = (/1,ishot/)
! call check(nf90_put_var(ncid,varid,real(u,4),start=start,count=count));
! call check(nf90_close(ncid))

!!	write(6,*) 'shot:',ishot,'ncid:',ncid,'varid',varid

!contains
!  subroutine check(status)
!    integer, intent ( in) :: status
!    
!    if(status /= nf90_noerr) then 
!      print *, trim(nf90_strerror(status))
!      stop 2
!    end if
!  end subroutine check  
!end subroutine write_ncdf
!!!!!!!!!!!!!!-----------------------------------------___!!!!!!!!!!
!subroutine read_ncdf(filename,u,n,ishot)
!use netcdf
!implicit none
!integer :: ncid,ishot,i
!character (len=*) :: filename
!character(4) :: appdump
!character (len = *), parameter :: U_NAME = "data"
!integer,parameter :: NDIMS=2
!integer :: n,varid;
!integer :: start(NDIMS),count(NDIMS);
!real*4 :: dump(n),u(n);dump=0.0;
!call check(nf90_open(filename,nf90_nowrite,ncid))
!call check(nf90_inq_varid(ncid,U_NAME,varid))
!count = (/n,1/)
!	start = (/1,ishot/)
!	call check(nf90_get_var(ncid,varid,dump,start=start,count=count))

!!	write(6,*) u(2)-dump(2),dump(2)
!	do i=1,n
!		if(dump(i)/=u(i))stop 2
!	enddo
!	call check(nf90_close(ncid))
!	write(6,*) 'checked'
!contains
!  subroutine check(status)
!    integer, intent ( in) :: status
!    if(status /= nf90_noerr) then 
!      print *, trim(nf90_strerror(status))
!      stop 2
!    end if
!  end subroutine check 
!end subroutine read_ncdf

!!subroutine load_snap_bdr2(u,rows,cols,ishot,folder)
!!implicit none
!!integer :: rows,cols,ishot
!!real*4,dimension(rows,cols),intent(out) :: u
!!	character(50) :: folder
!!write(6,*)'reading nc file ',trim(folder),' ...'
!!		call read_ncdf_matrix(trim(folder)//'.nc',u,rows,cols,ishot)

!!end subroutine load_snap_bdr2

!!!!!!!!!AUXUILIAR----------------------------------------NETCDF----------!!!!!!!
subroutine dump_matrix_ncdf(valore,ishot,rows,cols,filename)
implicit none
character(len=50), intent(in)              :: filename
integer,intent(in) :: rows,cols,ishot
real*4,dimension(rows,cols),intent(in) :: valore

write(6,*)'computing nc file ',trim(filename),' ...'
if(ishot==1) then
 		call init_ncdf_matrix(ishot,trim(filename),valore,rows,cols)
else
 		call write_ncdf_matrix(ishot,trim(filename),valore,rows,cols)
endif
end subroutine dump_matrix_ncdf

!!!!!!!!!AUXUILIAR----------------------------------------NETCDF----------!!!!!!!
subroutine init_ncdf_matrix(ishot,filename,u,rows,cols)
use netcdf
implicit none

integer :: ishot,i
character (len = *) :: filename
integer, parameter :: NDIMS = 3;
integer :: rows,cols
character (len = *), parameter :: U_NAME = "data"
character (len = *), parameter :: REC_NAME = "time"
character (len = *), parameter :: UNITS1 = "values"
character (len = *), parameter :: UNITS2 = "component"
integer :: ncid, varid, dimids(NDIMS)
integer :: u1_dimid,u2_dimid,t_dimid
integer :: start(NDIMS), count(NDIMS)
real*4, dimension(rows,cols):: u

 call check(nf90_create(filename,NF90_SHARE, ncid))
 call check(nf90_def_dim(ncid,UNITS1,rows,u1_dimid))
 call check(nf90_def_dim(ncid,UNITS2,cols,u2_dimid))
 call check(nf90_def_dim(ncid,REC_NAME,NF90_UNLIMITED,t_dimid))
 dimids=(/u1_dimid, u2_dimid, t_dimid/)
 call check(nf90_def_var(ncid,U_NAME,NF90_REAL,dimids,varid))
 call check(nf90_enddef(ncid))
 count = (/rows,cols,1/)
 start = (/1,1,1/)
 call check(nf90_put_var(ncid,varid,u,start=start,count=count));
 call check(nf90_close(ncid)) 

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check  
end subroutine init_ncdf_matrix

!!!---------------------------------------------------------------------!!!!!
subroutine write_ncdf_matrix(ishot,filename,u,rows,cols)
use netcdf
implicit none
integer :: ishot,i
character (len = *):: filename
character (len = *), parameter :: U_NAME = "data"
integer, parameter :: NDIMS = 3;
integer :: rows,cols,ncid,varid
integer :: start(NDIMS), count(NDIMS)
real*4, dimension(rows,cols)::u

 call check(nf90_open(filename,NF90_WRITE, ncid))
	call check(nf90_inq_varid(ncid,U_NAME,varid))
 count = (/rows,cols,1/)
 start = (/1,1,ishot/)
 call check(nf90_put_var(ncid,varid,u,start=start,count=count));
 call check(nf90_close(ncid))

!	write(6,*) 'shot:',ishot,'ncid:',ncid,'varid',varid

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check  
end subroutine write_ncdf_matrix


end module nc_routines
