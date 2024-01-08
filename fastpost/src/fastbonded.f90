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

!> Find the bonds in the system
subroutine fastfindbonds(n, atid, r, v0, va, vb, vc, &
        npairid, pairid, pairrc, maxnbd, nbd, bdid, bd, bdln, ierr)

    use vector3d
    use domain3d
    use argsorting

    implicit none

    integer(4), intent(in) :: n !< number of atoms
    integer(4), dimension(0:n-1), intent(in) :: atid !< atoms id
    real(8), dimension(0:2,0:n-1), intent(in) :: r !< atoms wrapped coordinates
    real(8), dimension(0:2), intent(in) :: v0, va, vb, vc !< box edges vectors
    integer(4), intent(in) :: npairid !< number of pairs
    integer(4), dimension(0:npairid-1,0:npairid-1), intent(in) :: pairid !< pair id
    real(8), dimension(0:npairid-1,0:npairid-1), intent(in) :: pairrc !< pair critical distnace
    integer(4), intent(in) :: maxnbd !< max # bonds
    integer(4), intent(out) :: nbd !< # bonds
    integer(4), dimension(0:maxnbd-1), intent(out) :: bdid   !< bond id
    integer(4), dimension(0:1,0:maxnbd-1), intent(out) :: bd !< bond atoms
    real(8), dimension(0:maxnbd-1), intent(out) :: bdln !< bond atoms
    integer(4), intent(out) :: ierr !< error flag (0: no error, <0 atom where maxatnb isexceed )

    type(SimulationBox) :: box                ! simulation cell
    type(BoxCells) :: cells                   ! grid sub-cells
    integer(4), dimension(0:2) :: dn          ! grid sub-division
    integer(4) :: ncc                         ! number of neighbor cells
    integer(4), dimension(27) :: cc           ! neighbor cells indexes
    integer(4) :: cnat                        ! cell number of atoms
    integer(4), pointer, dimension (:) :: cat ! cell atoms

    integer(4) :: iat, icell, iatid, jcell, j, jat, jatid
    integer(4) :: ic, ijbdid
    real(8) :: drsq, rcmax
    real(8), dimension(3) :: ri, dr
    real(8), dimension(0:npairid-1,0:npairid-1) :: pairrcsq

    ! initialize
    rcmax = maxval(pairrc)
    where(pairid .gt. 0) 
        pairrcsq = pairrc * pairrc
    elsewhere
        pairrcsq = 1.d0
    endwhere
    nbd = 0
    bdid = -1
    
    ! fill cells
    call box%initialize( v0, va, vb, vc, .true.)
    dn = [ max(3,int(1.d0/box%px/rcmax)), &
           max(3,int(1.d0/box%py/rcmax)), &
           max(3,int(1.d0/box%pz/rcmax)) ]
    call cells%splitAndfill( box, dn(0), dn(1), dn(2), n, r)
    if ( .not. cells%initialized) then
        ierr = 1
        return
    endif
    ! find bonds
    !$omp parallel do private(iat, ri, iatid, icell, &
    !$omp   ncc, cc, ic, jcell, cnat, cat, j, jat, jatid, &
    !$omp   ijbdid, dr, drsq)
    ATOMLOOP: do iat = 0, n-1
        ri = r(:,iat)
        iatid = atid(iat)
        icell = cells%findcell( ri) !< r(:,iat) no need to be wrapped
        call cells%getNeighbourCells( icell, ncc, cc) !< get the neighbour cells
        NBCELLSLOOP: do ic = 1, 27
            jcell = cc(ic)
            if ( jcell .eq. -1) cycle
            ! loop over cell's atoms
            call cells%getCellAtoms(jcell, cnat, cat)
            NBATOMSLOOP: do j = 1, cnat
                jat = cat(j) - 1
                if ( iat .ge. jat) cycle NBATOMSLOOP
                jatid = atid(jat)
                ijbdid = pairid(iatid,jatid)
                if ( ijbdid .lt. 0 ) cycle NBATOMSLOOP
                dr = r(:,jat) - ri
                call box%minImgTo( dr)
                drsq = square(dr)
                if ( drsq .le. pairrcsq(iatid,jatid)) then
                    !$omp critical
                    bd(0,nbd) = iat
                    bd(1,nbd) = jat
                    bdid(nbd) = ijbdid
                    bdln(nbd) = sqrt( drsq )
                    nbd = nbd + 1
                    if ( nbd .ge. maxnbd ) then
                        ierr = -2
                        return
                    endif
                    !$omp end critical
                endif
            end do NBATOMSLOOP
        end do NBCELLSLOOP
    end do ATOMLOOP
    !$omp end parallel do

end subroutine fastfindbonds
