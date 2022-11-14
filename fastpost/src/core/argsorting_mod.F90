module argsorting

   ! Module for sorting arrays.
   ! Based on code written by John E. Pask, LLNL.

   implicit none

   private

   ! overload argsort
   interface argsort
      module procedure iargsort, rargsort
   end interface

   public argsort

contains

   function iargsort(a) result(b)
      ! Returns the indices that would sort an array.
      !
      ! Arguments
      ! ---------
      !
      integer, intent(in):: a(:)    ! array of numbers
      integer :: b(size(a))         ! indices into the array 'a' that sort it
      !
      ! Example
      ! -------
      !
      ! iargsort([10, 9, 8, 7, 6])   ! Returns [5, 4, 3, 2, 1]

      integer :: N                           ! number of numbers/vectors
      integer :: i, imin                      ! indices: i, i of smallest
      integer :: temp                        ! temporary
      integer :: a2(size(a))
      a2 = a
      N = size(a)
      do i = 1, N
         b(i) = i
      end do
      do i = 1, N - 1
         ! find ith smallest in 'a'
         imin = minloc(a2(i:), 1) + i - 1

         ! swap to position i in 'a' and 'b', if not already there
         if (imin /= i) then
            temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
            temp = b(i); b(i) = b(imin); b(imin) = temp
         end if
      end do
   end function

   function rargsort(a) result(b)
      ! Returns the indices that would sort an array.
      !
      ! Arguments
      ! ---------
      !
      real*8, intent(in):: a(:)   ! array of numbers
      integer :: b(size(a))         ! indices into the array 'a' that sort it
      !
      ! Example
      ! -------
      !
      ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

      integer :: N                           ! number of numbers/vectors
      integer :: i, imin                      ! indices: i, i of smallest
      integer :: temp1                       ! temporary
      real*8   :: temp2
      real*8   :: a2(size(a))
      a2 = a
      N = size(a)
      do i = 1, N
         b(i) = i
      end do
      do i = 1, N - 1
         ! find ith smallest in 'a'
         imin = minloc(a2(i:), 1) + i - 1
         ! swap to position i in 'a' and 'b', if not already there
         if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
         end if
      end do
   end function

end module argsorting
