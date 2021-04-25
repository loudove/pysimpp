module vector3d

    implicit none

    private

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

    public square, length, dot, cross, project, unit, toUnit, angle

contains

    real(8) function m_square(v)
        real(8), dimension(3) :: v
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

        v_u = u * m_dot(v, u) / m_dot(u, u)

        return
    end function m_project

    function m_getUnit(v) result(u)
        real(8), dimension(3), intent(in) :: v
        real(8), dimension(3) :: u
        real(8) rlength
        rlength = 1.d0 / m_length(v)
        u = v * rlength
        return
    end function m_getUnit

    subroutine m_setToUnit(v)
        real(8), dimension(3), intent(inout) :: v
        real(8) rlength
        rlength = 1.d0 / m_length(v)
        v = v * rlength
        return
    end subroutine m_setToUnit

    real(8) function m_angle(v,u)
        real(8), dimension(3), intent(in) :: v, u
        real(8) vl, ul, rl, angle
        vl = m_length(v)
        ul = m_length(u)
        rl = 1.d0 / ( vl * ul )
         ! if we want to be safe
         angle = (v(1)*u(1)+v(2)*u(2)+v(3)*u(3))*rl
         if ( dabs(angle) > 1.d0 ) then
           if ( dabs(angle) - 1.d0 > 1.0d-06 ) stop
           angle = dsign(1.d0, angle)
         end if
        m_angle = dacos(angle)
        !m_angle = dacos((v(1)*u(1)+v(2)*u(2)+v(3)*u(3))*rl)
        return
    end function m_angle

endmodule 
