!> Accelerate bond lengths calculateion
subroutine fastbonds(n, r, va, vb, vc, nbnd, bnd, lengths)

    use vector3d
    use domain3d

    implicit none

    integer, intent(in) :: n                            !< number of atoms
    real*8, dimension(0:2,0:n-1), intent(in)  :: r      !< atoms unwrapped coordinates
    real*8, dimension(0:2), intent(in) :: va, vb, vc    !< simulation box spanning vectors
    integer, intent(in) :: nbnd                         !< number of bonds
    integer, dimension(0:1,0:nbnd-1), intent(in) :: bnd !< bond info b(i,:) = (atom1, atom2) for ith bond
    real*8, dimension(0:nbnd-1), intent(out) :: lengths !< bond lengths

    integer :: i
    real*8, dimension(0:2) :: dr
    type(SimulationBox) :: box
    
    lengths = 0.D0
    call box%initialize( [0.D0, 0.D0, 0.D0], va, vb, vc, .true.)

!$omp parallel do private(i, dr)
    do i = 0, nbnd-1
      dr = r( :, bnd(0,i)) - r( :, bnd(1,i))
      call box%minImgTo(dr)
      lengths(i) = length(dr)
    end do
!$omp end parallel do

end subroutine fastbonds

subroutine fastbondstest(n, r, va, vb, vc, nbnd, bnd, lengths)

    use vector3d
    use domain3d

    implicit none

    integer, intent(in) :: n                            !< number of atoms
    real*8, dimension(0:n-1,0:2), intent(in)  :: r      !< atoms unwrapped coordinates
    real*8, dimension(0:2), intent(in) :: va, vb, vc    !< simulation box spanning vectors
    integer, intent(in) :: nbnd                         !< number of bonds
    integer, dimension(0:nbnd-1,0:1), intent(in) :: bnd !< bond info b(i,:) = (type, atom1, atom2) for ith bond
    real*8, dimension(0:nbnd-1), intent(out) :: lengths !< bond lengths

    integer :: i
    real*8, dimension(0:2) :: dr
    type(SimulationBox) :: box
    
    lengths = 0.D0
    call box%initialize( [0.D0, 0.D0, 0.D0], va, vb, vc, .true.)

!$omp parallel do private(i, dr)
    do i = 0, nbnd-1
      dr = r( bnd(i,0), :) - r( bnd(i,1), :)
      call box%minImgTo(dr)
      lengths(i) = length(dr)
    end do
!$omp end parallel do

end subroutine fastbondstest

!> Accelerate bend angles calculation
subroutine fastangles(n, r, va, vb, vc, nang, ang, angles)

    use vector3d
    use domain3d  

    implicit none

    integer, intent(in) :: n                            !< number of atoms
    real*8, dimension(0:2,0:n-1), intent(in) :: r       !< atoms unwrapped coordinates
    real*8, dimension(0:2), intent(in) :: va, vb, vc    !< simulation box spanning vectors
    integer, intent(in) :: nang                         !< number of angles
    integer, dimension(0:2,0:nang-1), intent(in) :: ang !< angle info b(i,:) = (type, atom1, atom2, atom3) for ith angle
    real*8, dimension(0:nang-1), intent(out) :: angles  !< angle value

    integer :: i
    real*8, dimension(0:2) :: dr1, dr2
    type(SimulationBox) :: box

    ! real*8, parameter :: pi = dacos(-1.d0)
    ! real*8, parameter :: RAD2DEG = pi / 180.d0
    ! real*8, parameter :: EPS = 1.0D-12
    
    angles = 0.D0
    call box%initialize( [0.d0, 0.d0, 0.d0 ], va, vb, vc, .true.)

!$omp parallel do private(i, dr1, dr2)
    do i = 0, nang-1
      dr1 = r( :, ang(0,i)) - r( :, ang(1, i))
      dr2 = r( :, ang(2,i)) - r( :, ang(1, i))
      call box%minImgTo(dr1)
      call box%minImgTo(dr2)
      angles( i) = angle( dr1, dr2)
    end do
!$omp end parallel do

end subroutine fastangles

!> Accelerate bend angles calculation
subroutine fastdihedrals(n, r, va, vb, vc ,ndih, dih, angles)

    use vector3d
    use domain3d    

    implicit none

    integer, intent(in) :: n                            !< number of atoms
    real*8, dimension(0:2,0:n-1), intent(in) :: r       !< atoms unwrapped coordinates
    real*8, dimension(0:2), intent(in) :: va, vb, vc    !< simulation box spanning vectors
    integer, intent(in) :: ndih                         !< number of angles
    integer, dimension(0:3,0:ndih-1), intent(in) :: dih !< angle info b(i,:) = (type, atom1, atom2, atom3) for ith angle
    real*8, dimension(0:ndih-1), intent(out) :: angles  !< angle value

    integer :: i
    real*8, dimension(0:2) :: dr1, dr2, dr3
    real*8, dimension(0:2) :: dr1xdr2, dr2xdr3
    real*8 phi, s
    type(SimulationBox) :: box

    ! real*8, parameter :: pi = dacos(-1.d0)
    ! real*8, parameter :: RAD2DEG = pi / 180.d0
    ! real*8, parameter :: EPS = 1.0D-12

    angles = 0.D0
    call box%initialize( [0.d0, 0.d0, 0.d0 ], va, vb, vc, .true.)

!$omp parallel do private(i, dr1, dr2, dr3, dr1xdr2, dr2xdr3, phi, s)
    do i = 0, ndih-1
        ! bond vectors
        dr1 = r( :, dih(1, i)) - r( :, dih(0, i))
        dr2 = r( :, dih(2, i)) - r( :, dih(1, i))
        dr3 = r( :, dih(3, i)) - r( :, dih(2, i))
        call box%minImgTo(dr1)
        call box%minImgTo(dr2)
        call box%minImgTo(dr3)

        ! cis zero righ-handed
        dr1xdr2 = unit( cross(dr1, dr2))
        dr2xdr3 = unit( cross(dr2, dr3))
    
        phi = angle(dr1xdr2, dr2xdr3)
        s = dot(dr1xdr2, dr3)
    
        angles(i) = dsign(phi, s)

    end do
!$omp end parallel do

end subroutine fastdihedrals
