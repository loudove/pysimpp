      module test_domain3d

      use assertions_gtest
      use domain3d
      use vector3d

      implicit none

      real(8), parameter :: DEPS = epsilon(1.d0)
      real(8), parameter :: DPI = acos(-1.D0)
      real(8), dimension(3), parameter :: dv = [ DEPS, DEPS, DEPS]
      real(8), dimension(3), parameter :: zero = [ 0.D0, 0.D0, 0.D0]
      real(8), dimension(3), parameter :: one = [ 1.D0, 1.D0, 1.D0]
      real(8), dimension(3), parameter :: two = [ 2.D0, 2.D0, 2.D0]

      contains

      !> initialize random number generator.
      subroutine init_random()
            integer(4) :: n=1, iseed(1)=(/ 9234128 /)
            call random_seed(size=n)
            call random_seed(put=iseed(1:n))
      end subroutine

      !> initialize a simulation box using the given periodic and orthogonal flags.
      subroutine init_random_box(box, periodic, ortho)
            type(SimulationBox), intent(inout) :: box
            logical, intent(in) :: ortho, periodic
            real(8) :: xx, yy, zz, xy, xz, yz, tmp(6)
            call random_number(tmp)
            xx = ( 5.D0 + tmp(1) * 100.d0 )
            yy = ( 5.D0 + tmp(2) * 100.d0 )
            zz = ( 5.D0 + tmp(3) * 100.d0 )
            if ( .not. ortho ) then
                  xy = 0.3D0 * xx * ( 1.D0 - 2.D0 * tmp(4))
                  xz = 0.3D0 * xx * ( 1.D0 - 2.D0 * tmp(5))
                  yz = 0.3D0 * yy * ( 1.D0 - 2.D0 * tmp(6))
            else
                  xy=0.D0
                  xz=0.D0
                  yz=0.D0
            endif
            call box%initialize( zero, [ xx, 0.D0, 0.D0], [xy, yy, 0.D0], [xz, yz, zz], periodic)
      end subroutine init_random_box

      !> initialize a simulation box using the given periodic and orthogonal flags.
      function random_unit() result(u)
            real(8), dimension(3) :: u
            real(8) th, s, c, rnd(2)
            call random_number(rnd)
            th = 2.D0 * DPI * rnd(1)
            c = 1.D0 - 2.D0 * rnd(2)
            s = dsqrt(1.D0 - c * c)
            u = [ dcos(th)*s, dsin(th)*s, c ]
            return
      end function random_unit

      !$f90tw TESTCODE(TEST, domain3d, box_periodic, test_periodic, DUMMY)
      subroutine test_periodic() BIND(C,name="test_periodic")
      type(SimulationBox) :: box
      integer(4) i, j, k, cnt
      real(8) :: v(3), w(3), px, py, pz, p
      logical(KIND=C_BOOL) cond
      logical(4) :: periodic = .True., ortho = .False.
      call init_random()
      do cnt = 1, 100
            call init_random_box( box, periodic, ortho)
            call F90_ASSERT_EQ ( box%volume, dot(cross(box%vx,box%vy),box%vz), "box volume" // C_NULL_CHAR )
            call F90_ASSERT_EQ ( length( box%minImg(box%vx)), 0.D0, "pbc along a" // C_NULL_CHAR )
            call F90_ASSERT_EQ ( length( box%minImg(box%vy)), 0.D0, "pbc along b" // C_NULL_CHAR )
            call F90_ASSERT_EQ ( length( box%minImg(box%vz)), 0.D0, "pbc along c" // C_NULL_CHAR )
            v = box%minImgGet(box%vx,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,0,0]), KIND=C_BOOL ), "pbc indexes along a" // C_NULL_CHAR )
            v = box%minImgGet(box%vy,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,1,0]), KIND=C_BOOL ), "pbc indexes along b" // C_NULL_CHAR )
            v = box%minImgGet(box%vz,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,0,1]), KIND=C_BOOL ), "pbc indexes along c" // C_NULL_CHAR )
            v = box%minImgGet(box%vx+box%vy,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,1,0]), KIND=C_BOOL ), "pbc indexes along ab" // C_NULL_CHAR )
            v = box%minImgGet(box%vx+box%vz,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,0,1]), KIND=C_BOOL ), "pbc indexes along ac" // C_NULL_CHAR )
            v = box%minImgGet(box%vy+box%vz,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[0,1,1]), KIND=C_BOOL ), "pbc indexes along bc" // C_NULL_CHAR )
            v = box%minImgGet(box%vx+box%vy+box%vz,i,j,k)
            call F90_ASSERT_TRUE( logical( ALL([i,j,k]==[1,1,1]), KIND=C_BOOL ), "pbc indexes along abc" // C_NULL_CHAR )
            v =  0.4999999D0*(box%vx+box%vy+box%vz)
            call F90_ASSERT_DOUBLE_EQ ( length( box%minImg(v)), length(v), "pbc not applied" // C_NULL_CHAR )
            w = 0.25D0*(box%vx+box%vy+box%vz)
            call box%projections(px, py, pz)
            p = 0.25D0 * min( px, py, pz)
            do i = 1, 2000
                  call random_number(v) 
                  v =  w + random_unit() * p
                  call F90_ASSERT_DOUBLE_EQ ( length( box%minImg(v)), length(v), "pbc not applied" // C_NULL_CHAR )
            end do
            w = box%vx+box%vy+box%vz 
            p = 0.5D0 * min( px, py, pz)
            do i = 1, 2000
                  call random_number(v)
                  v =  w + random_unit() * p
                  call F90_ASSERT_NE ( length( box%minImg(v)), length(v), "pbc applied" // C_NULL_CHAR )
            end do
      end do
      end subroutine test_periodic

      !$f90tw TESTCODE(TEST, domain3d, box_matrix, test_matrix, DUMMY)
      subroutine test_matrix() BIND(C,name="test_matrix")

      real(8) :: myeps=DEPS*10000.d0
      logical, parameter :: periodic = .True., ortho = .True.
      integer(4) :: cnt
      type(SimulationBox) :: box
      do cnt = 1, 1000
            call init_random_box( box, periodic, ortho)
            call box%setMatrix()
            call F90_ASSERT_NEAR ( box%matrix(1), box%a, myeps, "m(1,1) matrix" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( box%matrix(2), 0.D0, myeps, "m(1,2) matrix" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( box%matrix(3), 0.D0, myeps, "m(1,3) matrix" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( box%matrix(4), box%b, myeps, "m(2,2) matrix" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( box%matrix(5), 0.D0, myeps, "m(2,3) matrix" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( box%matrix(6), box%c, myeps, "m(3,3) matrix" // C_NULL_CHAR )
      enddo
      end subroutine test_matrix

      !$f90tw TESTCODE(TEST, domain3d, box_fractional, test_fractional, DUMMY)
      subroutine test_fractional() BIND(C,name="test_fractional")
      real(8) :: myeps=DEPS*100.d0
      type(SimulationBox) :: box
      integer(4) :: i, j, k, cnt
      real(8) :: v(3), w(3,3), px, py, pz, p
      logical(KIND=C_BOOL) cond
      logical(4) :: periodic = .True., ortho = .False.
      do cnt = 1, 1000
            call init_random_box( box, periodic, ortho)
            w(:,1) = box%vx
            call box%toFractional(1,w)
            call F90_ASSERT_NEAR ( w(1,1), 1.D0, myeps, "a fractional" // C_NULL_CHAR )
            call box%toCartesian(1,w)
            call F90_ASSERT_NEAR ( w(1,1), box%vx(1), myeps, "a(1) cartesian" // C_NULL_CHAR )
            w(:,1) = box%vy 
            call box%toFractional(1,w)
            call F90_ASSERT_NEAR ( w(2,1), 1.D0, myeps, "b fractional" // C_NULL_CHAR )
            call box%toCartesian(1,w)
            call F90_ASSERT_NEAR ( w(1,1), box%vy(1), myeps, "b(1) cartesian" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( w(2,1), box%vy(2), myeps, "b(2) cartesian" // C_NULL_CHAR )
            w(:,1) = box%vz 
            call box%toFractional(1,w)
            call F90_ASSERT_NEAR ( w(3,1), 1.D0, myeps, "f fractional" // C_NULL_CHAR )
            call box%toCartesian(1,w)
            call F90_ASSERT_NEAR ( w(1,1), box%vz(1), myeps, "c(1) cartesian" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( w(2,1), box%vz(2), myeps, "c(2) cartesian" // C_NULL_CHAR )
            call F90_ASSERT_NEAR ( w(3,1), box%vz(3), myeps, "c(3) cartesian" // C_NULL_CHAR )
            do i = 1, 2000
                  call random_number(w(:,1)) 
                  v = w(:,1)
                  call box%toCartesian(1,w) 
                  call box%toFractional(1,w)
                  call F90_ASSERT_NEAR ( w(1,1), v(1), myeps, "f(1) fractional" // C_NULL_CHAR )
                  call F90_ASSERT_NEAR ( w(2,1), v(2), myeps, "f(2) fractional" // C_NULL_CHAR )
                  call F90_ASSERT_NEAR ( w(3,1), v(3), myeps, "f(3) fractional" // C_NULL_CHAR )
                  call random_number(w)
                  w(:,1) = w(:,1) * box%vx 
                  w(:,2) = w(:,2) * box%vy 
                  w(:,3) = w(:,3) * box%vz
                  v = w(:,1) + w(:,2) + w(:,3) 
                  w(:,1) = v
                  call box%toFractional(1,w) 
                  call box%toCartesian(1,w)
                  call F90_ASSERT_NEAR ( w(1,1), v(1), myeps*10.D0, "c(1) cartesian" // C_NULL_CHAR )
                  call F90_ASSERT_NEAR ( w(2,1), v(2), myeps*10.D0, "c(2) cartesian" // C_NULL_CHAR )
                  call F90_ASSERT_NEAR ( w(3,1), v(3), myeps*10.D0, "c(3) cartesian" // C_NULL_CHAR )
            enddo
      end do
      end subroutine test_fractional

      !$f90tw TESTCODE(TEST, domain3d, cells_fill, test_fill, DUMMY)
      subroutine test_fill() BIND(C,name="test_fill")
      real(8) :: myeps=DEPS*100.d0
      type(SimulationBox) :: box
      type(BoxCells) :: celllist
      integer(4) :: i, j, k, n(3), nsites, isite, cnt, cnat
      real(8) :: px, py, pz, dx(3), hdx(3), p(3), w(3)
      real(8), allocatable, dimension(:,:) :: x
      logical(KIND=C_BOOL) cond 
      logical(4) :: periodic = .True., ortho = .False.
      integer(4), pointer, dimension (:) :: cat
      integer(4) :: ncc 
      integer(4), dimension(27) :: cc
      do cnt = 1, 100
            call init_random_box( box, periodic, ortho)
            call random_number(p) 
            n = 3 + int(50.d0*p,kind=4)
            nsites = n(1) * n(2) * n(3)
            if ( allocated(x) .and. size(x) < 3*nsites ) deallocate(x)
            if ( .not. allocated(x) ) allocate( x(3, nsites))
            isite = 1 
            dx = [ 1.D0/n(1), 1.D0/n(2), 1.D0/n(3)] 
            hdx = 0.5D0 * dx
            do k = 1, n(3) 
                  pz = (k-1)*dx(3)+hdx(3) 
                  do j = 1, n(2)  
                        py = (j-1)*dx(2)+hdx(2) 
                        do i = 1, n(1)
                              x(:,isite) = [ (i-1)*dx(1)+hdx(1), py, pz] 
                              isite = isite + 1
                        enddo
                  enddo
            enddo
            call box%toCartesian( nsites, x(:,1:nsites))
            ! call random_number(w) 
            ! w = (1.D0 - 2.D0*w) * 10.D0
            ! do i = 1, nsites 
            !      x(:,i) = x(:,i) + w 
            ! enddo
            call celllist%splitAndfill(box, n(1), n(2), n(3), nsites, x )
            do i = 1, nsites
                  call celllist%getCellAtoms(i,cnat,cat)
                  call F90_ASSERT_EQ( cnat, 1, "incorrect number of sites" // C_NULL_CHAR )
                  call F90_ASSERT_EQ( cat(1), i, "incorrect cell-site bind" // C_NULL_CHAR )
                  call celllist%getNeighbourCells(i,ncc,cc)
                  call F90_ASSERT_EQ( ncc, 27, "incorrect number of neighbors" // C_NULL_CHAR )
                  k = 0 
                  do j = 1, 27 
                        if ( cc(j) .ne. -1) k = k + 1
                  enddo
                  call F90_ASSERT_EQ( k, 27, "incorrect number neighbor cells" // C_NULL_CHAR )
            enddo
            do i = 1, nsites
                  call random_number(p) 
                  p = (one - two * p) * hdx
                  x(:,i) = x(:,i) - w + p
            enddo
            do i = 1, nsites
                  call celllist%getCellAtoms(i,cnat,cat)
                  call F90_ASSERT_EQ( cnat, 1, "incorrect number of sites" // C_NULL_CHAR )
                  call celllist%getNeighbourCells(i,ncc,cc)
                  call F90_ASSERT_EQ( ncc, 27, "incorrect number of neighbors" // C_NULL_CHAR )
                  k = 0 
                  do j = 1, 27 
                        if ( cc(j) .ne. -1) k = k + 1 
                  enddo
                  call F90_ASSERT_EQ( k, 27, "incorrect number of neighbor cells" // C_NULL_CHAR )
            enddo
            call celllist%destroy()
      end do
      end subroutine test_fill

end module test_domain3d
