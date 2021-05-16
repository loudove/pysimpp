module vector3d

    use iso_fortran_env, only : stderr=>error_unit, stdout=>output_unit

    implicit none

    private

    real(8), parameter :: EPS = epsilon(1.0)
    integer(4) :: FASTPOSTDEBUG = 0

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
        module procedure m_angle
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

    public debug
    public square, length, dot, cross, project, unit, toUnit, angle, tostr

contains

    subroutine m_debug_set(ivar) BIND(C,name="m_debug_set")
        integer(4), intent(in) :: ivar
            FASTPOSTDEBUG = ivar
    end subroutine m_debug_set

    subroutine m_debug_env() BIND(C,name="m_debug_env")
        character(LEN=255) :: var
        integer(4) :: ivar, ierr
        call get_environment_variable("FASTPOSTDEBUG", var)
        read( var, *, IOSTAT=ierr) ivar
        if ( ierr .ne. 0 ) FASTPOSTDEBUG = ivar
    end subroutine m_debug_env

    real(8) function m_square(v)
        real(8), dimension(3), intent(in) :: v
        m_square = v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
        return
    end function m_square

    real(8) function m_length(v)
        real(8), dimension(3) :: v
        m_length = dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
        return
    end function m_length

    function m_dot(v, u) result(dot)
        real(8), dimension(3), intent(in) :: v, u
        real(8) :: dot

        dot = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)

        return
    end function m_dot

    function m_cross(v, u) result(vxu)
        real(8), dimension(3), intent(in) :: v, u
        real(8), dimension(3) :: vxu

        vxu(1) = v(2)*u(3)-v(3)*u(2)
        vxu(2) = v(3)*u(1)-v(1)*u(3)
        vxu(3) = v(1)*u(2)-v(2)*u(1)

        return
    end function m_cross

    !> project v to u vector.
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
        return
    end function m_project

    function m_getUnit(v) result(u)
        real(8), dimension(3), intent(in) :: v
        real(8), dimension(3) :: u
        u = v
        call m_setToUnit(u)
        return
    end function m_getUnit

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

    real(8) function m_angle(v,u)
        real(8), dimension(3), intent(in) :: v, u
        real(8) sqv, squ, angle
        sqv = m_square(v)
        squ = m_square(u)
        if ( sqv .gt. 0.d0 .and. squ .gt. 0.d0 ) then
            angle = (v(1)*u(1)+v(2)*u(2)+v(3)*u(3))/dsqrt( sqv * squ )
            if ( dabs(angle) > 1.d0 ) then
                if ( dabs(angle) - 1.d0 > EPS ) stop
                angle = dsign(1.d0, angle)
            end if
            m_angle = dacos(angle)
        else
            m_angle = 0.d0
            write(stderr,'("vector3d:ZeroAngle")')
            write(stdout,'("vector3d:ZeroAngle")')           
        endif
        return
    end function m_angle

    function m_tostr(v) result(str)
        real(8), dimension(3) :: v
        character(LEN=46) str
        write( str, '("(",E14.7,",",E14.7,",",E14.7")")') v
        return
    end function m_tostr

endmodule 
