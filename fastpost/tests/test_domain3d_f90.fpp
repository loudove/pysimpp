#include "f90tw_test.h"

      TESTMODULE(test_domain3d, domain3d)

#ifdef TESTF90
      real(8), parameter :: DEPS = epsilon(1.d0) ;
      real(8), dimension(3), parameter :: dv = [ DEPS, DEPS, DEPS] ;
      real(8), dimension(3), parameter :: zero = [ 0.D0, 0.D0, 0.D0] ;
      real(8), dimension(3), parameter :: one = [ 1.D0, 1.D0, 1.D0] ;
#endif

TESTCONTAINS

#ifdef TESTF90
      subroutine init_random()
            integer(4) :: n=1, iseed(1)=(/ 9234128/);
            call random_seed(size=n) ;
            call random_seed(put=iseed(1:n)) ;            
      end subroutine
#endif

TESTCODE(TEST, domain3d, simulationbox, test_periodic,
      use vector3d ;
      type(SimulationBox) :: box ;
      real(8) :: minbox ;
      integer(4) i, j, k, cnt ;
      real(8) :: a, b, c, av(3), bv(3), cv(3), v(3), w(3) ;
      logical(KIND=C_BOOL) cond; logical(4) :: periodic = .True. ;
      do cnt = 1, 1000 ;
            call random_number(minbox); minbox = minbox + 0.5D0 ;
            call random_number(a) ; call random_number(b) ; call random_number(c) ;
            a = (1.D0 + a) * minbox ; b = ( 1.D0 + b ) * minbox ; c = ( 1.D0 + c ) * minbox ;
            av = [ a, 0.D0, 0.D0] ; bv = [ 0.D0, b, 0.D0] ; cv = [ 0.D0, 0.D0, c ] ;
            call box%initialize( zero, av, bv, cv, periodic) ;
            call F90_ASSERT_EQ ( box%volume, a*b*c, "box volume" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_EQ ( length( box%minImg(av)), 0.D0, "pbc along a" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_EQ ( length( box%minImg(bv)), 0.D0, "pbc along b" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_EQ ( length( box%minImg(cv)), 0.D0, "pbc along c" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(av,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,0,0]), KIND=C_BOOL ), "pbc indexes along a" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(bv,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,1,0]), KIND=C_BOOL ), "pbc indexes along b" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(cv,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,0,1]), KIND=C_BOOL ), "pbc indexes along c" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(av+bv,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,1,0]), KIND=C_BOOL ), "pbc indexes along c" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(av+cv,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,0,1]), KIND=C_BOOL ), "pbc indexes along c" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet(bv+cv,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,1,1]), KIND=C_BOOL ), "pbc indexes along c" F90CONCAT C_NULL_CHAR ) ;
            v = box%minImgGet([a,b,c]*one,i,j,k) ;
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,1,1]), KIND=C_BOOL ), "pbc indexes along c" F90CONCAT C_NULL_CHAR ) ;
            do i = 1, 1000 ;
                  call random_number(v) ; v = 0.5D0 * [ a, b, c] - [ a, b, c] * v ;
                  call F90_ASSERT_DOUBLE_EQ ( length( box%minImg(v)), length(v), "pbc not applied" F90CONCAT C_NULL_CHAR ) ;
                  call random_number(v) ; v = 0.5D0 * [ a, b, c] * ( one + v ) ;
                  call F90_ASSERT_NE ( length( box%minImg(v)), length(v), "pbc applied" F90CONCAT C_NULL_CHAR ) ;
            end do ;
      end do
)

ENDTESTMODULE(test_domain3d)
