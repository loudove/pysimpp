#include "f90tw_test.h"

      TESTMODULE(test_vector3d, vector3d)

#ifdef TESTF90
      real(8), parameter :: DEPS = epsilon(1.d0) ;
      real(8), dimension(3), parameter :: zero = [ 0.d0,0.d0,0.d0] ;
      real(8), dimension(3), parameter :: one = [ 1.d0,1.d0,1.d0] ;
      real(8), dimension(3), parameter :: x = [ 1.d0,0.d0,0.d0] ;
      real(8), dimension(3), parameter :: y = [ 0.d0,1.d0,0.d0] ;
      real(8), dimension(3), parameter :: z = [ 0.d0,0.d0,1.d0] ;
      real(8), dimension(3), parameter :: dv = [ DEPS, DEPS, DEPS] ;
#endif

TESTCONTAINS

#ifdef TESTF90
      subroutine init_random()
            integer(4) :: n=1, iseed(1)=(/ 9234128/);
            call random_seed(size=n) ;
            call random_seed(put=iseed(1:n)) ;            
      end subroutine
#endif

TESTCODE(TEST, vector3d, simple, test_simple,
      real(8), dimension(3) :: v ;
      logical(KIND=C_BOOL) cond;
      call F90_ASSERT_EQ ( square(zero), 0.d0 ) ;
      call F90_ASSERT_EQ ( square(one), 3.d0 ) ;
      call F90_ASSERT_EQ ( length(zero), 0.d0 ) ;
      call F90_ASSERT_EQ ( length(one), dsqrt(3.d0) ) ;
      call F90_ASSERT_EQ ( dot(zero, one), 0.d0 ) ;
      call F90_ASSERT_EQ ( dot(x, y), 0.d0 ) ;
      call F90_ASSERT_EQ ( dot(one, one), 3.d0 ) ;
      call F90_ASSERT_EQ ( angle(zero, zero), 0.d0 ) ; ! zero angle
      call F90_ASSERT_EQ ( angle(zero, one), 0.d0 ) ;  ! zero angle
      call F90_ASSERT_EQ ( angle(one, one), 0.d0 ) ;
      call F90_ASSERT_TRUE ( logical( ALL ( cross( x, y) == z ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(cross(y,z) == x ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(cross(z,x) == y ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(cross(y,x) == -z ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(cross(z,y) == -x ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(cross(x,z) == -y ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(project(x,y) == zero ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(project(one,one) == one ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(project(one,zero) == zero ), kind=C_BOOL) ) ;
      call F90_ASSERT_TRUE ( logical( ALL(project(zero,one) == zero ), kind=C_BOOL) ) ; 
      call F90_ASSERT_EQ ( length( unit(one)), 1.d0 ) ;
      call F90_ASSERT_EQ ( length( unit(x)), 1.d0  ) ;
      call F90_ASSERT_EQ ( length( unit(y)), 1.d0  ) ;
      call F90_ASSERT_EQ ( length( unit(z)), 1.d0  ) ;
      call F90_ASSERT_EQ ( length( unit(zero)), 0.d0 ) ;
      v = one ; call toUnit(v) ;
      call F90_ASSERT_EQ ( length( v), 1.d0) ;
      v = x ; call toUnit(v) ;
      call F90_ASSERT_EQ ( length( v), 1.d0) ;
)

TESTCODE(TEST, vector3d, square, test_square,
      call F90_ASSERT_DOUBLE_EQ( square(one+dv), square(one),
            "quare(one+dv)=square(one)" F90CONCAT C_NULL_CHAR ) ;
      call F90_ASSERT_DOUBLE_EQ( length(1.D04*one+dv), length(1.D04*one),
            "quare(1.D04*one+dv)=square(1.D04*one)" F90CONCAT C_NULL_CHAR ) ;
)

TESTCODE(TEST, vector3d, linear, test_linear,
      real(8) :: a ; integer(4) :: i;
      real(8), dimension(3) :: v, w ;
      real(8) :: myeps=DEPS*1.d5 ;
      call init_random() ;
      do i = 1, 100000 ;
            call random_number(a) ; call random_number(v) ; call random_number(w) ;
            call F90_ASSERT_DOUBLE_EQ( length(a*v), a*length(v), 
                              "length(a*v)=a*lengt(v)" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_DOUBLE_EQ( angle(v,w), angle(-v,-w), 
                              "angle(v,w)=angle(-v,-w)" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_NEAR( length( project(a*v,w)), a*length( project(v,w)), myeps, 
                              "project(a*v,w)=a*project(a*v,w)" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_NEAR( angle(v,w), dacos(-1.d0)-angle(-1.d0*v,w), myeps, 
                              "angle(v,w)=pi-angle(v,w)" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_NEAR( length( cross(v,w)), length(v)*length(w)*dsin(angle(v,w)), myeps,
                              "length(cross(v,w))=length(v)*length(w)*sin(angle(v,w)" F90CONCAT C_NULL_CHAR ) ;
            call F90_ASSERT_TRUE ( logical( ALL(cross(v,w) == -cross(w,v) ), kind=C_BOOL),
                              "cross(v,w)=-cross(w,v)" F90CONCAT C_NULL_CHAR ) ; 

      enddo
)

ENDTESTMODULE(test_vector3d)
