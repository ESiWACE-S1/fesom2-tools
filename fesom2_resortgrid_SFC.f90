program resort_grid

  implicit none

  integer              :: nod2d, elem2d, nlvl

  real(kind=8), allocatable :: lat(:), lon(:)
  character(len=20), allocatable :: depth(:) , lvl(:) 
  integer, allocatable      :: tri(:,:), sfc(:), idx(:)
  integer, allocatable      :: sfc_perm(:,:)
  

  character(len=100)   :: MeshPath 
  character(len=100)   :: File_nod, File_elem, File_depth
  character(len=100)   :: File_nod_sfc, File_elem_sfc, File_depth_sfc

  integer :: n, i, a, b, c
  real    :: x,y

  if (iargc() /=1) then
     print *,'Usage: ./resortgrid /path/to/mesh'
     stop
  endif
  call getarg(1,MeshPath)
  
  File_nod   = trim(MeshPath)//"nod2d.out"
  File_elem  = trim(MeshPath)//"elem2d.out"
  File_depth = trim(MeshPath)//"aux3d.out"

  File_nod_sfc   = trim(File_nod)//"_sfc"
  File_elem_sfc  = trim(File_elem)//"_sfc"
  File_depth_sfc = trim(File_depth)//"_sfc"

  
! Reorder the nodes along the SFC
  open(51,file=trim(File_nod),status='old')
  read(51,*) nod2d
  allocate(lat(nod2d),lon(nod2d),idx(nod2d),sfc(nod2d))
  do n=1,nod2d
     read(51,*) i,x,y,b
     lat(i) = x
     lon(i) = y
     idx(i) = b
  enddo
  close(51)

  call compute_sfc_index(nod2d,lat,lon,sfc)

  allocate(sfc_perm(2,nod2d))
  
  do n=1,nod2d
     sfc_perm(1,n) = sfc(n)
     sfc_perm(2,n) = n
  enddo

  call quick_sort(sfc_perm,2,nod2d)

    
  
  open(61,file=trim(File_nod_sfc),status='new')
  write(61,*) nod2d
  do n=1,nod2d
     write(61,'(i8,f10.4,f10.4,i3)') n, lat(sfc_perm(2,n)), lon(sfc_perm(2,n)), idx(sfc_perm(2,n))
  enddo
  close(61)

  deallocate(lat,lon,idx,sfc)

! Reorder the depth in the same way
! Strings are used here to retain the accuracy
  allocate(depth(nod2d))

! Adjusted for FESOM2 aux3d
  open(52,file=trim(File_depth),status='old')
  read (52,*) nlvl
  allocate(lvl(nlvl))
  do n=1,nlvl
     read(52,*) lvl(n)
  end do
 do n=1,nod2d
     read(52,*) depth(n)
  enddo
  close(52)
  open(62,file=trim(File_depth_sfc),status='new')
  write(62,*) nlvl
  do n=1,nlvl
     write(62,*) trim(lvl(n))
  enddo
  do n=1,nod2d
     write(62,*) trim(depth(sfc_perm(2,n)))
  enddo
  close(62)


! And the elements. We need the inverse permutation - reuse sfc_perm(1,1:nod2d)
  do n=1,nod2d
     sfc_perm(1,sfc_perm(2,n)) = n
  enddo

  open(53,file=trim(File_elem),status='old')
  read(53,*) elem2d
  allocate(tri(3,elem2d))
  do n=1,elem2d
     read(53,*) tri(1:3,n)
  enddo
  close(53)

  do n=1,elem2d
! new node number (inverse SFC perm)
     a = sfc_perm(1,tri(1,n))
     b = sfc_perm(1,tri(2,n))
     c = sfc_perm(1,tri(3,n))

! re-order each elements: smallest number first. Do not change orientation!
     if (a < min (b,c)) then
        tri(1:3,n) = (/ a, b, c /)
     elseif (b < min(a,c)) then
        tri(1:3,n) = (/ b, c, a /)
     else
        tri(1:3,n) = (/ c, a, b /)
     endif
  enddo
  
! resort the elements numerically
  call quick_sort(tri,3,elem2d)


  open(63,file=trim(File_elem_sfc),status='new')
  write(63,*) elem2d
  do n=1,elem2d
     write(63,'(3i8)') tri(1:3,n)
  enddo
  close(63)

  deallocate(sfc_perm)

end program resort_grid

!==============================

subroutine compute_sfc_index(dim,x,y,sfc)

  implicit none

  integer, intent(in)          :: dim
  double precision, intent(in) :: x(dim), y(dim)
  integer, intent(out)         :: sfc(dim)

  integer            :: i, j, k
  integer, parameter :: maxlevel = 14
  integer            :: sfc_add
  integer            :: halve
  integer            :: lc
  double precision   :: xmin, xmax, xscale, ymin, ymax, yscale, mx, ny, xhalve, yhalve
  integer            :: nmax
  integer            :: q
  integer, parameter :: A=1, D=2, C=3, U=4
  
  nmax = 2**maxlevel

  xmin = minval(x)
  xmax = maxval(x)
  ymin = minval(y)
  ymax = maxval(y)

  xscale = real(nmax,8)/(xmax-xmin)
  yscale = real(nmax,8)/(ymax-ymin)
    
  do i=1,dim

    sfc(i) = 0

    mx = x(i) - xmin
    ny = y(i) - ymin

    xhalve = xmax-xmin
    yhalve = ymax-ymin
    
    q = A
    
    sfc_add = nmax*nmax

    do j =1,maxlevel 

       sfc_add = sfc_add / 4

       xhalve = .5*xhalve
       yhalve = .5*yhalve

       if (mx < xhalve) then
          lc = 0
       else
          lc = 2
          mx = mx - xhalve
       end if
       if (ny >= yhalve) then
          lc = lc+1
          ny = ny - yhalve
       endif
      
       select case (q)
	
          case(A)

             select case(lc)
	  
             case (0)
                q = D
             case(2)
                q = C
                sfc(i) = sfc(i) + 3*sfc_add
             case(3)
                sfc(i) = sfc(i) + 2*sfc_add
             case(1)
                sfc(i) = sfc(i) +   sfc_add   	  
             end select

          case(U)
	
             select case(lc)

             case(3)
                q = C
             case(1)
                q = D
                sfc(i) = sfc(i) + 3*sfc_add	  
             case(0)
                sfc(i) = sfc(i) + 2*sfc_add	  
             case(2)
                sfc(i) = sfc(i) + sfc_add
             end select
	
          case(C)
	
             select case(lc)
	  
             case(3) 
                q = U
             case(2)
                q = A
                sfc(i) = sfc(i) + 3*sfc_add	  
             case(0)
                sfc(i) = sfc(i) + 2*sfc_add	  
             case(1)
                sfc(i) = sfc(i) + sfc_add
             end select

          case(D)
	
             select case(lc)
	  
             case(0)
                q = A
             case(1)
                q = U
                sfc(i) = sfc(i) + 3*sfc_add	  
             case(3)
                sfc(i) = sfc(i) + 2*sfc_add	  
             case(2)
                sfc(i) = sfc(i) + sfc_add
             end select
          
          end select
       end do
    end do
  
 
end subroutine compute_sfc_index
!======================================================================
! Fast in-line QSORT+INSERTION SORT for Fortran.
! Author: Joseph M. Krahn
! FILE: qsort_inline.inc
!
! Rewritten to subroutine for one specific task (sort integer array along
!   the first dimension), Natalja Rakowsky, AWI, 2015
!
! PURPOSE:
! Generate a custom array sort procedure for a specific type,
! without the comparison-callback overhead of a generic sort procedure.
! This is essentially the same as an in-line optimization, which generally
! is not feasible for a library-based generic sort procedure.
!
! This implementation is as generic as possible, while avoiding the need
! for a code pre-processor. The success of this approach assumes that
! internal procedures are always in-lined with optimized Fortran compilation.
!
! USAGE:
! This file contains the sort subroutine body. You must supply
! an integer parameter QSORT_THRESHOLD, and internal procedures:
!    subroutine SWAP(a,b)
!    subroutine RSHIFT(left,right)
!
! Descriptions:
!
!
! SUBROUTINE SWAP(A,B)
!   Swap array members 'a' and 'b'
!
! SUBROUTINE RSHIFT(LEFT,RIGHT)
!   Perform a circular shift of array members LEFT through RIGHT,
!   shifting the element at RIGHT back to the position at LEFT.
!
! QSORT_THRESHOLD:
!   The QSORT is used down to the QSORT_THRESHOLD size sorted blocks.
!   Then insertion sort is used for the remainder, because it is faster
!   for small sort ranges. The optimal size is not critical. Most of
!   the benefit is in blocks of 8 or less, and values of 16 to 128
!   are generally about equal speed. However, the optimal value
!   depends a lot on the hardware and the data being sorted, so this
!   is left as a tunable parameter for cases where ther is an
!   effect on performance.
!
!---------------------------------------------------------------------
! NOTES:
! The procedure uses a optimized combination of QSORT and INSERTION
! sorting. The algorithm is based on code used in GLIBC. 
! A stack is used in place of recursive calls. The stack size must
! be at least as big as the number of bits in the largest array index.
!
! Sorting vectors of a multidimensional allocatable array can be
! VERY slow. In this case, or with large derived types, it is better
! to sort a simple derived type of key/index pairs, then reorder
! tha actual data using the sorted indices.
!
!---------------------------------------------------------------------
subroutine quick_sort(array,dim,array_size)

  implicit none

  integer, intent(in)    :: dim, array_size
  integer, intent(inout) :: array(dim,array_size)

  integer, parameter :: QSORT_THRESHOLD=16

  integer :: stack_top, right_size, left_size
  integer :: mid, left, right, low, high

! A stack of 32 can handle the entire extent of a 32-bit
! index, so this value is fixed. If you have 64-bit indexed
! arrays, which might contain more thant 2^32 elements, this
! should be set to 64.
  integer, parameter :: QSORT_STACK_SIZE = 32
  type qsort_stack 
     integer :: low, high; 
  end type qsort_stack
  type(qsort_stack) :: stack(QSORT_STACK_SIZE)

  integer :: n



  if (array_size > QSORT_THRESHOLD) then
    low = 1
    high = array_size
    stack_top = 0

    QSORT_LOOP: &
    do
      mid = (low + high)/2
      if (array(1,mid) < array(1,low)) then
        call SWAP(mid,low)
      end if
      if (array(1,high) < array(1,mid)) then
        call SWAP(high,mid)
        if (array(1,mid) < array(1,low)) then
          call SWAP(mid,low)
        end if
      end if
      left  = low + 1
      right = high - 1

      COLLAPSE_WALLS: &
      do
        do while (array(1,left) < array(1,mid))
          left=left+1
        end do
        do while (array(1,mid) < array(1,right))
          right=right-1
        end do
        if (left < right) then
          call SWAP(left,right)
          if (mid == left) then
            mid = right
          else if (mid == right) then
            mid = left
          end if
          left=left+1
          right=right-1
        else
          if (left == right) then
            left=left+1
            right=right-1
          end if
          exit COLLAPSE_WALLS
        end if
      end do COLLAPSE_WALLS

! Set up indices for the next iteration.
! Determine left and right partition sizes.
! Defer partitions smaller than the QSORT_THRESHOLD.
! If both partitions are significant,
! push the larger one onto the stack.
      right_size = right - low
      left_size = high - left
      if (right_size <= QSORT_THRESHOLD) then
        if (left_size <= QSORT_THRESHOLD) then
          ! Ignore both small partitions: Pop a partition or exit.
          if (stack_top<1) exit QSORT_LOOP
          low=stack(stack_top)%low; high=stack(stack_top)%high
          stack_top=stack_top-1
        else
          ! Ignore small left partition.
          low = left
        end if
      else if (left_size <= QSORT_THRESHOLD) then
        ! Ignore small right partition.
        high = right
      else if (right_size > left_size) then
        ! Push larger left partition indices.
        stack_top=stack_top+1
        stack(stack_top)=qsort_stack(low,right)
        low = left
      else
        ! Push larger right partition indices.
        stack_top=stack_top+1
        stack(stack_top)=qsort_stack(left,high)
        high = right
      end if
    end do QSORT_LOOP
  end if ! (array_size > QSORT_THRESHOLD)

! Sort the remaining small partitions using insertion sort,
! which should be faster for partitions smaller than the
! appropriate QSORT_THRESHOLD.

! First, find smallest element in first QSORT_THRESHOLD and
! place it at the array's beginning. This places a lower
! bound 'guard' position, and speeds up the inner loop
! below, because it will not need a lower-bound test.
  low = 1
  high = array_size

! left is the MIN_LOC index here:
  left=low
  do right = low+1, min(low+QSORT_THRESHOLD,high)
    if (array(1,right) < array(1,left)) left=right
  end do
  if (left/=low) call SWAP(left,low)

! Insertion sort, from left to right.
! (assuming that the left is the lowest numbered index)
  INSERTION_SORT: &
  do right = low+2,high
    left=right-1
    if (array(1,right) < array(1,left)) then
      do while (array(1,right) < array(1,left-1))
        left=left-1
      end do
      call RSHIFT(left,right)
    end if
  end do INSERTION_SORT

contains

  subroutine swap(m,n)
    integer, intent(in) :: m, n
    integer :: hold(dim)

    hold(1:dim)   = array(1:dim,m)
    array(1:dim,m) = array(1:dim,n)
    array(1:dim,n) = hold(1:dim)

  end subroutine swap

  subroutine rshift(left,right)
    implicit none
    integer, intent(in) :: left, right
    integer :: hold(dim), i

    hold(1:dim)=array(1:dim,right)
    do i=right,left+1,-1
      array(1:dim,i)=array(1:dim,i-1)
    end do
    array(1:dim,left)=hold
    
  end subroutine rshift

end subroutine quick_sort
!--------------------------------------------------------------
