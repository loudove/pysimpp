module vector3d

    use iso_c_binding, only : c_int
    use iso_fortran_env, only : stderr=>error_unit, stdout=>output_unit

    implicit none

    private

    real(8), parameter :: EPS = epsilon(1.0)
    real(8), parameter :: PI = acos(-1.d0)
    real(8), parameter :: TWOPI = 2.0*PI
    integer(c_int) :: FASTPOSTDEBUG = 0

    interface tostr
        module procedure m_tostr
    end interface

    interface square
        module procedure m_square
    end interface

    interface length
        module procedure m_length
    end interface

    interface angle
        module procedure m_angle_vectors, m_angle_points
    end interface

    interface dihedral
        module procedure m_dihedral_vectors, m_dihedral_points
    end interface

    interface cross
        module procedure m_cross
    end interface

    interface dot
        module procedure m_dot
    end interface

    interface project
        module procedure m_project
    end interface

    interface normal
        module procedure m_normal
    end interface

    interface decompose
        module procedure m_decompose
    end interface

    interface rotate
        module procedure m_rotate
    end interface

    interface random
        module procedure m_random_on_sphere
        module procedure m_random_vertical
    end interface

    interface unit
        module procedure m_getUnit
    end interface

    interface toUnit
        module procedure m_setToUnit
    end interface  

    interface debug
        module procedure m_debug_set
        module procedure m_debug_env
    end interface

    interface build_tetrahedral
        module procedure m_fast_build_tetrahedral
    end interface

    interface build_dihedral
        module procedure m_build_dihedral
    end interface

    public debug, tostr
    public square, length, angle, dihedral
    public dot, cross, unit, toUnit, random
    public project, normal, decompose, rotate
    public build_tetrahedral, build_dihedral

contains

    !> set the debug state of the module
    subroutine m_debug_set(ivar) BIND(C,name="m_debug_set")
        integer(c_int), intent(in) :: ivar

        FASTPOSTDEBUG = int(ivar, kind=4)
    end subroutine m_debug_set

    !> read environment variable ${FASTPOSTDEBUG} and set 
    !> the debug state of the module
    subroutine m_debug_env() BIND(C,name="m_debug_env")
        character(LEN=255) :: var
        integer(c_int) :: ivar
        integer(4) :: ierr

        call get_environment_variable("FASTPOSTDEBUG", var)
        read( var, *, IOSTAT=ierr) ivar
        if ( ierr .ne. 0 ) FASTPOSTDEBUG = ivar
    end subroutine m_debug_env

    !> return the square length
    real(8) function m_square(v)
        real(8), dimension(3), intent(in) :: v

        m_square = v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
    end function m_square

    !> return the length
    real(8) function m_length(v)
        real(8), dimension(3) :: v

        m_length = dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
    end function m_length

    !> return the dot product v.u
    function m_dot(v, u) result(d)
        real(8), dimension(3), intent(in) :: v, u
        real(8) :: d

        d = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)
    end function m_dot

    !> return the cross product v x u
    function m_cross(v, u) result(vxu)
        real(8), dimension(3), intent(in) :: v, u
        real(8), dimension(3) :: vxu

        vxu(1) = v(2)*u(3)-v(3)*u(2)
        vxu(2) = v(3)*u(1)-v(1)*u(3)
        vxu(3) = v(1)*u(2)-v(2)*u(1)
    end function m_cross

    !> return the projection of v on u vector
    function m_project(v, u) result(v_u)
        real(8), dimension(3), intent(in) :: v, u
        real(8), dimension(3) :: v_u
        real(8) :: sq

        sq = m_square(u)
        if ( sq .gt. 0.d0 ) then
            v_u = u * m_dot(v, u) / sq
        else
            v_u = (/ 0.d0, 0.d0, 0.d0 /)
            write(stderr,'("vector3d:ZeroProject")')
            write(stdout,'("vector3d:ZeroProject")')
        endif
    end function m_project

    !> return the unit v
    function m_getUnit(v) result(u)
        real(8), dimension(3), intent(in) :: v
        real(8), dimension(3) :: u
        u = v
        call m_setToUnit(u)
        return
    end function m_getUnit

    !> set v to unit
    subroutine m_setToUnit(v)
        real(8), dimension(3), intent(inout) :: v
        real(8) sq
        sq = m_square(v)
        if ( sq .gt. 0.d0 ) then
            v = v / dsqrt(sq)
        else
            write(stderr,'("vector3d:ZeroUnit")')
            write(stdout,'("vector3d:ZeroUnit")')           
        endif
        return
    end subroutine m_setToUnit

    !> return the angles between the vectors v and u
    function m_angle_vectors(v,u) result(theta)
        real(8), dimension(3), intent(in) :: v, u
        real(8) :: sqv, squ
        real(8) :: theta

        sqv = m_square(v)
        squ = m_square(u)
        if ( sqv .gt. 0.d0 .and. squ .gt. 0.d0 ) then
            theta = (v(1)*u(1)+v(2)*u(2)+v(3)*u(3))/dsqrt( sqv * squ )
            if ( dabs(theta) > 1.d0 ) then
                if ( dabs(theta) - 1.d0 > EPS ) stop
                theta = dsign(1.d0, theta)
            end if
            theta = dacos(theta)
        else
            theta = 0.d0
            write(stderr,'("vector3d:ZeroAngle")')
            write(stdout,'("vector3d:ZeroAngle")')           
        endif
    end function m_angle_vectors

    ! return the angle defined by the given sequence of points
    function m_angle_points(p1, p2, p3) result(theta)
        real(8), dimension(3), intent(in) :: p1, p2, p3
        real(8) :: theta
        
        theta = m_angle_vectors(p1-p2, p3-p2)
    end function m_angle_points

    !> return the dihedral form by the central vector u
    !> and the adjacent vectors v and w using trans-zero, 
    !> right-hand (anti-clockwise) convention
    function m_dihedral_vectors(v, u, w) result(phi)
        real(8), dimension(3), intent(in) :: v, u, w
        real(8) :: phi
        real(8), dimension(3) :: vxu, uxw

        vxu = m_cross(u,v)
        uxw = m_cross(u,w)
        phi = sign( m_angle_vectors( vxu, uxw), m_dot(w, vxu))
    end function m_dihedral_vectors

    !> return the dihedral defined by the given sequence of points
    !> using trans-zero, right-hand (anti-clockwise) convention
    function m_dihedral_points(p1, p2, p3, p4) result(phi)
        real(8), dimension(3), intent(in) :: p1, p2, p3, p4
        real(8) :: phi

        phi = m_dihedral_vectors(p1-p2, p3-p2, p4-p2)
    end function m_dihedral_points

    !> return the normal to u, defined by v
    function m_normal(v, u) result( n)
        real(8), dimension(3), intent(in) :: v, u
        real(8), dimension(3) :: n

        n = v -  m_project(v, u)
    end function m_normal
    
    !> decompose vector v to its projection p and normal n with
    !> respect to verctir u
    subroutine m_decompose(v, u, n, p) 
        real(8), dimension(3), intent(in) :: v, u
        real(8), dimension(3), intent(out) :: n, p

        p =  m_project(v, u)
        n = v - p
    end subroutine m_decompose

    !> rotate unit vector v around unit vector u by 
    !> phi=acos(c)=asin(s), right-handed
    function m_rotate(v, u, s, c) result( w)
        real(8), dimension(3), intent(in) :: v, u
        real(8), intent(in) :: s, c
        real(8), dimension(3) :: w  
        real(8) :: fact
        
        fact = m_dot( v, u) * ( 1.d0 - c )
        w = c * v + fact * u + s * m_cross(u,v)
    end function m_rotate

    !> return a vector with its end uniformaly distributed
    !> on the the surface of a unit shphere. it is assumed 
    !> that the intrinsic Random Number Generator of fortran
    !> compiler has been initialized properly by the caller
    function m_random_on_sphere() result (u)
        real(8), dimension(3) :: u
        real(8) :: theta, phi, s, rnd
        call random_number( rnd)
        theta = PI * rnd
        call random_number( rnd)
        phi = acos(1.d0-2.d0*rnd);
        s = sin(phi);
        u = (/ cos(theta)*s, sin(theta)*s, cos(phi) /)
    end function m_random_on_sphere

    !> return a random unit vector verticle to vector v. 
    !> it is assumed that the intrinsic Random Number
    !> Generator of fortran compiler has been initialized
    !> properly by the caller
    function m_random_vertical(v) result (u)
        real(8), dimension(3), intent(in) :: v
        real(8), dimension(3) :: u, uv
        real(8) lv, phi, s, c, rnd;

        u = (/ -v(2), v(1)+v(3), -v(2) /)
        call m_setToUnit(u)
        lv = m_length(v)
        uv = v / lv
        call random_number( rnd)
        phi = TWOPI * rnd
        s = sin(phi)
        c = cos(phi)
        u = m_rotate(u,uv,s,c)
     end function m_random_vertical

    !> given three points [p0, p1, p2] build the fourth [p3] and
    !> the fifth [p4] points forming a tetrahedral structure with
    !> p1 positioned at its center. the length of the bonds [p4-p1]
    !> and [p5-p1] will set to b while the angle (p4,p1,p5) will set
    !> to theta. the construction is symmetrical, i.e., points p3 and
    !> p4 are equidistant from the plane [p0, p1, p3], and the plane
    !> [p3, p1, p4] bisects angle (p0,p1,p2).
    subroutine m_fast_build_tetrahedral(p0, p1, p2, b, theta, p3, p4)
        real(8), dimension(3), intent(in) :: p0, p1, p2
        real(8), intent(in) :: b, theta
        real(8), dimension(3), intent(out) :: p3, p4

        real(8) :: halftheta
        real(8), dimension(3) :: v0, v1
        real(8), dimension(3) :: uy, uz

        !> uy lies in plane [p0,p1,p2] and bisects the angle (p0,p1,p2)
        v0 = p0 - p1
        v1 = p2 - p1
        uy = m_getUnit( -(v0+v1))
        !> v0 is normal to uy defined by v1
        v0 = m_normal( v1, uy)
        !> uz is normal to the plane [p0,p1,p2]
        uz = m_cross( uy, v0)
        call m_setToUnit( uz)
        halftheta = 0.5 * theta
        v0 = p1 + uy * b * cos( halftheta)
        v1 = uz * b * sin( halftheta)
        p3 = v0 + v1
        p4 = v0 - v1
    end subroutine m_fast_build_tetrahedral
     
    !> given three points [p0, p1, p2] build the fourth point [p3]
    !> with using its internal dofs (bond [b], angle [theta], and
    !> dihedral [phi]). theta is defined as the angle between the
    !> vectors (p0-p1, p3-p2) and the dihedrla using the trans zero
    !> right-hand (anti-clockwise) convention
    function m_build_dihedral( p0, p1, p2, b, theta, phi) result (p3)
        real(8), dimension(3), intent(in) :: p0, p1, p2
        real(8), intent(in) :: b, theta, phi
        real(8), dimension(3) :: p3

        real(8), dimension(3) :: ux, uy, uz
        real(8) :: s, c

        !> construct local frame [ux, uy, uz]
        ux = m_getUnit( p2-p1)
        uy = m_normal( p1-p0, ux)
        call m_setToUnit( uy)
        uz = m_cross( ux, uy)
        !> find the new positions
        s = sin(theta)
        c = cos(theta)
        p3 = p2 - b * (ux * c - &
                       uy * s * cos(phi) - &
                       uz * s * sin(phi))
    end function m_build_dihedral

    function m_tostr(v) result(str)
        real(8), dimension(3) :: v
        character(LEN=46) str
        write( str, '("(",E14.7,",",E14.7,",",E14.7")")') v
        return
    end function m_tostr

endmodule 
