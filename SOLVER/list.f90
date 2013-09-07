!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
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

!=================
module linked_list
!=================

  implicit none

  private
  public :: list

  type :: link
     private
     type(link), pointer        :: next => null()   ! next link in list
     type(link), pointer        :: prev => null()   ! next link in list
     integer                    :: ldata
     contains
     procedure, pass :: getData      ! return value pointer
     procedure, pass :: getNextLink  ! return next pointer
     procedure, pass :: setNextLink  ! set next pointer
     procedure, pass :: getPrevLink  ! return next pointer
     procedure, pass :: setPrevLink  ! set next pointer
  end type link

  interface link
     procedure constructor ! construct/initialize a link
  end interface

  type :: list
     private
     class(link), pointer :: firstLink => null()    ! first link in list
     class(link), pointer :: lastLink => null()     ! last link in list
     class(link), pointer :: currentLink => null()  ! iterator
     integer              :: length = 0
     contains
     procedure, pass :: free            ! empty the list and free the memory
     procedure, pass :: append          ! append an element to the end of the list
     procedure, pass :: insert          ! insert an element to beginning of the list
     procedure, pass :: getFirst        ! return first element
     procedure, pass :: getLast         ! return last element
     procedure, pass :: getCurrent      ! iterator, can be moved with getNext
                                        ! if not set, returns first element
     procedure, pass :: resetCurrent    ! reset the iterator
     procedure, pass :: getNext         ! get the next element
                                        ! if current not set, returns first element
     procedure, pass :: getPrev         ! get the previous element
                                        ! if current not set, returns last element
     procedure, pass :: getLength       ! return length of the list
  end type list

contains

!-----------------------------------------------------------------------------------------
function getNextLink(this)
  class(link)           :: this
  type(link), pointer   :: getNextLink
  getNextLink => this%next
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setNextLink(this,next)
  class(link)           :: this
  type(link), pointer   :: next
  this%next => next
end subroutine setNextLink
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getPrevLink(this)
  class(link)           :: this
  type(link), pointer   :: getPrevLink
  getPrevLink => this%prev
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setPrevLink(this,prev)
  class(link)           :: this
  type(link), pointer   :: prev
  this%prev => prev
end subroutine setPrevLink
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getData(this)
  class(link)           :: this
  getData = this%ldata
end function getData
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function constructor(ldata, next, prev)
  class(link), pointer    :: constructor
  integer                 :: ldata
  type(link), pointer     :: next, prev

  allocate(constructor)
  constructor%next => next
  constructor%prev => prev
  constructor%ldata = ldata
end function constructor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine append(this, ldata)
  class(list)             :: this
  integer                 :: ldata
  class(link), pointer    :: newLink

  if (.not. associated(this%firstLink)) then
     this%firstLink => link(ldata, null(), null())
     this%lastLink => this%firstLink
  else
     newLink => link(ldata, null(), this%lastLink)
     call this%lastLink%setNextLink(newLink)
     this%lastLink => newLink
  end if
   
  this%length = this%length + 1
end subroutine append
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine insert(this, ldata)
  class(list)             :: this
  integer                 :: ldata
  class(link), pointer    :: newLink

  if (.not. associated(this%firstLink)) then
     this%firstLink => link(ldata, null(), null())
     this%lastLink => this%firstLink
  else
     newLink => link(ldata, this%firstLink, null())
     call this%firstLink%setPrevLink(newLink)
     this%firstLink => newLink
  end if

  this%length = this%length + 1
end subroutine insert
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getFirst(this)
  class(list)           :: this

  if(associated(this%firstLink)) then
     getFirst = this%firstLink%getData()
  else
     error stop 'trying to go access data, but list is empty'
  endif
end function getFirst
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getLast(this)
  class(list)           :: this

  if(associated(this%lastLink)) then
     getLast = this%lastLink%getData()
  else
     error stop 'trying to go access data, but list is empty'
  endif
end function getLast
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getCurrent(this)
  class(list)           :: this
  if (.not. associated(this%currentLink)) then
     if(associated(this%firstLink)) then
        this%currentLink => this%firstLink
     else
        error stop 'trying to go access data, but list is empty'
     endif
  endif
  getCurrent = this%currentLink%getData()
end function getCurrent
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine resetCurrent(this)
  class(list)           :: this
  this%currentLink => null()
end subroutine resetCurrent
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getNext(this)
  class(list)           :: this
  if (.not. associated(this%currentLink)) then
     this%currentLink => this%firstLink
     getNext = this%currentLink%getData()
  elseif (associated(this%currentLink%getNextLink())) then
     this%currentLink => this%currentLink%getNextLink()
     getNext = this%currentLink%getData()
  else
     error stop 'trying to go beyond last element in list'
  end if 
end function getNext
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getPrev(this)
  class(list)           :: this
  if (.not. associated(this%currentLink)) then
     this%currentLink => this%lastLink
     getPrev = this%currentLink%getData()
  elseif (associated(this%currentLink%getPrevLink())) then
     this%currentLink => this%currentLink%getPrevLink()
     getPrev = this%currentLink%getData()
  else
     error stop 'trying to go beyond first element in list'
  end if 
end function getPrev
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getLength(this)
  class(list)           :: this
  getLength = this%length
end function getLength
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine free(this)
  class(list)               :: this
  class(link), pointer      :: current, next

  if(associated(this%firstLink)) then
     next => this%firstLink
     do while ( associated(next) )
        current => next
        write(6,*) 'free', current%getData()
        next => current%getNextLink()
        deallocate( current )
     enddo
     this%firstLink => null()
  endif

  this%length = 0
end subroutine free
!-----------------------------------------------------------------------------------------

!=====================
end module linked_list
!=====================


program test_list
  use linked_list
  implicit none

  type(list)                :: l

  call l%append(1)
  write(6,*) 'bla'
  write(6,*) l%getFirst()
  write(6,*) l%getCurrent()
  write(6,*) l%getLast()
  call l%insert(0)
  call l%append(2)
  write(6,*) 'bla'
  write(6,*) l%getFirst()
  write(6,*) l%getCurrent()
  write(6,*) l%getLast()
  call l%append(3)
  write(6,*) 'bla'
  write(6,*) 'length', l%getLength()
  write(6,*) l%getFirst()
  write(6,*) l%getCurrent()
  write(6,*) l%getLast()
  write(6,*) 'bla'
  call l%resetCurrent()
  write(6,*) l%getPrev()
  write(6,*) l%getPrev()

  write(6,*) 'bla'
  call l%free()
  write(6,*) 'length', l%getLength()
  write(6,*) 'bla'
  write(6,*) l%getFirst()
end program
