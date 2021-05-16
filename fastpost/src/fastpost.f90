!> Accelerate the msd calculation using fft. Calculate the
!> s1 term for the msd formula for the one dimensional case
!> (http://dx.doi.org/10.1051/sfn/201112010)
subroutine fasts1(n, d, s1)

      implicit none

      integer, intent(in) :: n
      real*8, dimension(0:n-1), intent(in) :: d
      real*8, dimension(0:n-1), intent(out) :: s1

      integer i
      real*8 :: q

      q = 2.d0*sum( d)
      s1 = 0.d0

      s1(0) = q / n
      do i = 1, n-1
        q = q - d(i-1) - d(n-i)
        s1(i) = q / (n-i)
      end do

end subroutine fasts1

!> Accelerate the msd calculation using fft. Calculate the
!> s1 term for the msd formula for the m dimensional case.
!> Boost fasts1 even more.
!> (http://dx.doi.org/10.1051/sfn/201112010)
subroutine fasts1x(n, m, d, s1)

      implicit none

      integer, intent(in) :: n, m
      real*8, dimension(0:n-1,0:m-1), intent(in) :: d
      real*8, dimension(0:m-1,0:n-1), intent(out) :: s1

      integer i, j
      !real*8, dimension(0:m-1) :: q
      real*8 :: q

      s1 = 0.d0
      do j = 0, m-1
        q = 2.d0*sum( d(:,j))
        s1(j,0) = q / n
        do i = 1, n-1
          q = q - d(i-1,j) - d(n-i,j)
          s1(j,i) = q / (n-i)
        end do
      end do

      !q = 2.d0*sum( d, dim=1)
      !s1(:,0) = q / n
      !do i = 1, n-1
      !   q = q - d(i-1,:) - d(n-i,:)
      !   s1(:,i) = q / (n-i)
      !enddo


end subroutine fasts1x

! Accelerate placing coordinates in (n,3) array and wrapping them (per
! configuration) using the periodic indexes
subroutine fastunwrapv(n, rw, ip, a, b, c, r)
    implicit none
      integer, intent(in) :: n
      real*8, dimension(0:n-1,0:2), intent(in) :: rw
      integer, dimension(0:n-1,0:2), intent(in) :: ip
      real*8, dimension( 0:2), intent(in) :: a, b, c
      real*8, dimension(0:n-1,0:2), intent(out) :: r

      integer :: i

      r = 0.d0
      do i = 0, n-1
          r(i,0:2) = rw(i,0:2) + ip(i,0) * a  + ip(i,1) * b + ip(i,2) * c
      end do

end subroutine fastunwrapv

!> Accelerate placing coordinates in (n,3) array and wrapping them
!> (per configuration) using the periodic indexes
subroutine fastunwrapf(n, x, y, z, ix, iy, iz, va, vb, vc, r)

    implicit none

    integer, intent(in) :: n                            !< # of atoms
    real*8, dimension( 0:n-1), intent(in) :: x, y, z    !< atoms wrapped coordinates
    integer, dimension( 0:n-1), intent(in) :: ix, iy, iz !< atoms periodic indexes
    real*8, dimension( 0:2), intent(in) :: va, vb, vc   !< box spaning vectors
    real*8, dimension( 0:2, 0:n-1), intent(out) :: r    !< atoms unwrapped coordinates

    integer :: i
    real*8, dimension( 0:2) :: v

    r = 0.d0
!$omp parallel do private(i, v)
    do i = 0, n-1
        v =  (/ x(i), y(i), z(i) /)
        r(0:2,i) = v + ix(i) * va  + iy(i) * vb + iz(i) * vc
    end do
!$omp end parallel do

end subroutine fastunwrapf

!> Accelerate placing coordinates in (n,3) array and wrapping them
!> (per configuration) using the periodic indexes
subroutine fastunwrap(n, rw, ip, va, vb, vc, r)

    implicit none

    integer, intent(in) :: n                            !< # of atoms
    real*8, dimension( 0:2, 0:n-1), intent(in) :: rw    !< atoms wrapped coordinates
    integer, dimension( 0:2, 0:n-1), intent(in) :: ip   !< atoms periodic indexes
    real*8, dimension( 0:2), intent(in) :: va, vb, vc   !< box spaning vectors
    real*8, dimension( 0:2, 0:n-1), intent(out) :: r    !< atoms unwrapped coordinates

    integer :: i
    real*8, dimension( 0:2) :: v

    r = 0.d0
!$omp parallel do private(i)
    do i = 0, n-1
        r(:,i) = rw(:,i) + ip(0,i) * va  + ip(1,i) * vb + ip(2,i) * vc
    end do
!$omp end parallel do

end subroutine fastunwrap

!> Accelerate making molecules whole again (ONLY cubic boxes)
subroutine fastwhole(n, rw, molecule, nbonds, bonds, a, b, c, r)
    implicit none
    integer, intent(in) :: n                        !< number of atoms
    real*8, dimension(0:n-1,0:2), intent(in) :: rw  !< atoms wrapped coordinates
    !> atoms molecule (zero indexed array)
    integer, dimension(0:n-1), intent(in) :: molecule
    integer, intent(in) :: nbonds                   !< number of bonds
    !> bond atoms b(i,:) = (atom1, atom2) for ith bond
    integer, dimension(0:nbonds-1,0:1), intent(in) :: bonds
    real*8, intent(in) :: a, b, c                   !< box edges
    real*8, dimension(0:n-1,0:2), intent(out) :: r  !< atoms whole coordinates

    integer :: nmolecules, i, ii, j, jj, k, kk, it, im, ib, ib1, ib2, maxn

    integer, allocatable, dimension(:) :: molnat    !< # of atoms per molecule
    !> starting index for the atoms of the molecules in the molat array
    integer, allocatable, dimension(:) :: moliat
    integer, dimension(0:n-1) :: molat  !< atoms per molecule array
    integer, dimension(0:n-1) :: iatmol !< atom local numbering in the molecule

    integer, dimension(0:n-1) :: atnbnd !< # of bonds per atom
    !> starting index for the bonds of the atoms in the atbnd array
    integer, dimension(0:n-1) :: atibnd
    integer, allocatable, dimension(:) :: atbnd     !< atoms bonds

    integer :: ntrack, itrack, nmask, maxmolnat, nbuff
    integer, allocatable, dimension(:) :: track, buff
    logical, allocatable, dimension(:) :: mask
    real(8), dimension(0:2) :: box, dr

    box = (/ a, b, c /) ! box array
    r = rw              ! initialize whole coordinates

    ! find the number of atoms per molecules
    nmolecules = maxval( molecule) + 1
    allocate( molnat(0:nmolecules-1), moliat(0:nmolecules-1))
    molnat = 0
    do i=0,n-1
        im = molecule(i)
        molnat( im) = molnat( im) + 1
    enddo
    ! construct the atoms per molecule (molat) array and
    ! the corresponding indexing array (moiat)
    ! construct the atoms local numbering array (iatmol)
    moliat(0) = 0
    do im = 1, nmolecules-1
        moliat(im) = moliat(im-1)+molnat(im-1)
    enddo
    molnat = 0
    do i=0,n-1
        im = molecule(i)
        ii = moliat(im)+molnat(im)
        molat( ii) = i
        iatmol(i) = ii
        molnat( im) = molnat( im) + 1
    enddo

    ! find the number of bonds per atom
    atnbnd = 0
    do i=0,nbonds-1
        ib1 = bonds(i,0)
        ib2 = bonds(i,1)
        atnbnd(ib1) = atnbnd(ib1) + 1
        atnbnd(ib2) = atnbnd(ib2) + 1
    enddo
    ! construct the bonds per atom array (atbnd) and the
    ! corresponding index array (atibnd)
    atibnd(0) = 0
    do i = 1, n-1
        atibnd(i) = atibnd(i-1)+atnbnd(i-1)
    enddo
    allocate( atbnd( int( sum(atnbnd))))
    atnbnd = 0
    do i=0,nbonds-1
        ib1 = bonds(i,0)
        ib2 = bonds(i,1)
        atbnd( atibnd(ib1)+atnbnd(ib1)) = ib2
        atbnd( atibnd(ib2)+atnbnd(ib2)) = ib1
        atnbnd(ib1) = atnbnd(ib1) + 1
        atnbnd(ib2) = atnbnd(ib2) + 1
    enddo

    ! allocate utility arrays and make the molecules whole
    maxmolnat = int( maxval(molnat))
    allocate ( track( 0:maxmolnat-1), mask( 0:maxmolnat-1))
    allocate ( buff(maxmolnat))
    do im = 0, nmolecules-1 ! loop over the molecules
        ! track the fragments in the molecule
        ntrack= 0
        track = -1
        do while ( .not. any(track .eq. -1) ) ! do while every atom is assigned to a track/fragment
            ! find the first atom of the current track/fragment
            nbuff = 0
            do ii = 0, molnat(im)-1
                if ( track(ii) .eq. -1 ) then
                    track(ii) = ntrack
                    nbuff = nbuff + 1
                    buff(nbuff) = ii
                    exit
                endif
            enddo
            ! find the rest of the atoms of the current track/fragment
            do while ( nbuff .gt. 0)
                ii = buff(nbuff)
                nbuff = nbuff - 1
                i = molat( moliat(im)+ii)
                ! loop over atom bonded neighbors and if they are not
                ! tracked already and their corresponding bonds do not
                ! cross box boundaries assign them in the same fragment
                do ib = atibnd(i), atibnd(i)+atnbnd(i)-1
                    j = atbnd(ib)
                    jj = iatmol(j)
                    if ( .not. track(jj) .eq. -1 ) continue
                    dr = rw(j,:) - rw(i,:)
                    if ( all( anint(dr/box) .eq. 0 ) ) then
                        track(jj) = ntrack
                        nbuff = nbuff + 1
                        buff(nbuff) = jj
                    endif
                enddo
            enddo
            ! increase the current number of fragments by one
            ntrack = ntrack + 1
        enddo

        ! if only one fragment has been identified for the molecule do nothing
        ! (either the whole molecule is in the box or it is infinit periodic)
        if ( ntrack .eq. 1 ) continue

        ! make it whole
        ! trace the larger fragment and start from there
        itrack = 0
        maxn = count(track .eq. itrack)
        do it = 1, ntrack-1
            mask = track .eq. it
            nmask = count(mask)
            if ( nmask .gt. maxn ) then
                itrack = it
                maxn = nmask
            endif
        enddo

        ! find the connections of other fragments
        ! with itrack and join them.
        nbuff = 1
        buff(nbuff) = itrack
        mask = .false. ! use mask as atoms state bit
        where (track .eq. itrack) mask = .true.
        do while ( nbuff .gt. 0)
            it = buff(nbuff)
            nbuff = nbuff - 1
            ! loop over the atoms fo the current fragment
            do ii = 0, molnat(im)-1
                if (track(ii) .eq. it) then
                    i = molat( moliat(im)+ii)
                    ! loop over the neighbors and check the ones that
                    ! are not placed to find neighboring fragments not
                    ! processed so far
                    do ib = atibnd(i), atibnd(i)+atnbnd(i)-1
                        j = atbnd(ib)
                        jj = iatmol(j)
                        ! if its allready placed continue
                        if ( mask(jj)) continue
                        ! append the track/fragment of the atom into the buffer
                        nbuff = nbuff + 1
                        buff(nbuff) = track(jj)
                        ! translate all the atoms of the track/fragmet and mark
                        ! them as placed
                        dr = rw(j,:) - rw(i,:)
                        dr = - anint(dr/box)
                        do kk = 0, molnat(im)-1
                            if ( track(kk) .eq. track(jj) ) then
                                k = molat( moliat(im)+ii)
                                r(k,:) = r(i,:) + dr
                                mask(kk) = .true.
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
    enddo

    deallocate( track, mask, buff)
    deallocate( atbnd)
    deallocate( molnat)
    deallocate( moliat)

end subroutine fastwhole

!> Accelerate wrapping coordinates per configuration for orthogonal boxes
subroutine fastwrapo(n, r, v0, a, b, c, rw)

    implicit none

    integer, intent(in) :: n                          !< # of atoms
    real*8, dimension( 0:2, 0:n-1), intent(in) :: r   !< atoms coordinates
    real*8, dimension( 0:2), intent(in) :: v0         !< box origin
    real*8, intent(in) :: a, b, c                     !< box edges
    real*8, dimension( 0:2, 0:n-1), intent(out) :: rw !< atoms wrapped coordinates

    integer :: i, j, k
    real*8, dimension( 0:2) :: v, d

    rw = 0.d0

    d = (/ a, b, c /)
!$omp parallel do private(i,j,k,v)
    do i = 0, n-1
        v =  r(:,i) - v0
        do j = 0, 2
            k = int( v(j) / d(j) )
            if ( k /= 0       ) v(j) = v(j) - k * d(j)
            if ( v(j) < 0.d0  ) v(j) = v(j) + d(j)
            if ( v(j) >= d(j) ) v(j) = v(j) - d(j)
        end do
        rw(:,i) = v0 + v
    end do
!$omp end parallel do

end subroutine fastwrapo

!> Accelerate wrapping coordinates per configuration
subroutine fastwrap(n, r, v0, va, vb, vc, rw)

    use domain3d

    implicit none

    integer, intent(in) :: n                          !< # of atoms
    real*8, dimension( 0:2, 0:n-1), intent(in) :: r   !< atoms coordinates
    real*8, dimension( 0:2), intent(in) :: v0         !< box origin
    real*8, dimension( 0:2), intent(in) :: va, vb, vc !< box spaning vectors
    real*8, dimension( 0:2, 0:n-1), intent(out) :: rw !< atoms wrapped coordinates

    integer :: i, j, k
    real*8, dimension( 0:2) :: v, d
    type(SimulationBox) :: box

    rw = r
    call box%initialize( v0, va, vb, vc, .true.)
    call box%toFractional( n, rw)
    call fastwrapo(n, rw, v0, 1.D0, 1.D0, 1.D0, rw)
    call box%toCartesian(n, rw)

end subroutine fastwrap

!> Accelerate coordinateas unwrapping for a polyethylene model
!> (only ortho boxes are supported for now)
subroutine pefastunwrap(n, rw, nch, maxnatch, natch, chat , a, b, c, r)
    implicit none
    !> # of atoms, # of chains, max atoms per chain
    integer, intent(in) :: n, nch, maxnatch
    !> unwrapped coordinates
    real*8, dimension(0:n-1,0:2), intent(in) :: rw
    !> # of atoms per chain
    integer, dimension(0:nch-1), intent(in) :: natch
    !> atoms for each chain
    integer, dimension(0:nch-1, 0:maxnatch-1), intent(in) :: chat
    !> box spanning vectors
    real*8, dimension(0:2), intent(in) :: a, b, c
    !> unwrapped coordinates
    real*8, dimension(0:n-1,0:2), intent(out) :: r

    integer :: i, j, k, l
    real*8, dimension( 0:2) :: dr
    real*8, dimension(0:2) :: box, rbox

    box = (/ a(0), b(1), c(2) /)
    rbox = 1.d0 / box
    r = 0.d0

    do i = 0, nch-1
        k = chat(i,0)
        r( k, 0:2) = rw( k, 0:2)
        do j = 1, natch(i)-1
            l = chat(i,j)
            k = chat(i,j-1)
            dr = rw(l,:) - rw(k,:)
            dr = dr - box * anint(dr*rbox)
            r(l,0:2) = r(k,0:2) + dr
        end do
    end do

end subroutine pefastunwrap

!> calculate energy and virial pressure (molecular) for pe
!> ua-model using trappe force field (CH2 for all the atoms)
!> (only ortho boxes are supported for now)
subroutine peenrgpress(n, r, rw, nch, rcm, atch , a, b, c, temperature,     &
                       rcut, enrg, einter, eintra, eimage, etrunc,                  &
                       press, apress, pid, ptrunc, mstten, astten)

!enrg,einter,eintra,eimage,press,apress,mstens[:,:],astens[:,:]
!    = peenrgpress(r,rw,rcm,molecules,box.va, box.vb, box.vc, T, rcut)
    use vector3d
    use domain3d
    use lj_interactions

    implicit none

    !> # of atom, # of chains
    integer, intent(in) :: n, nch
    !> unwrapped and wrapped coordinates
    real*8, dimension(0:n-1,0:2), intent(in) :: r, rw
    !> unwrapped center of mass coordinates
    real*8, dimension(0:nch-1,0:2), intent(in) :: rcm
    !> chain for each atom
    integer, dimension(0:n-1), intent(in) :: atch
    !>  cutoff distance
    real*8, intent(inout) :: rcut
    !> temperature
    real*8, intent(in) :: temperature
    !> box spaning vectors
    real*8, dimension(0:2), intent(in) :: a, b, c
    !> energy term
    real*8, intent(out) :: enrg, einter, eintra, eimage, etrunc
    !> pressure
    real*8,  intent(out) :: press, apress, pid, ptrunc
    !> atomic and molecular stress tensor
    real*8, dimension(0:2,0:2), intent(out) :: mstten, astten
    ! real*8, dimension(0:2,0:2) :: stcm

    call d_peenrgpress(n, r, rw, nch, rcm, atch , a, b, c, temperature,     &
    rcut, enrg, einter, eintra, eimage, etrunc,                  &
    press, apress, pid, ptrunc, mstten, astten)
!
    ! !> box stuff
    ! type(SimulationBox) :: box
    ! real*8, dimension(0:2) :: r0 = (/ 0.d0, 0.d0, 0.d0/)
    ! logical(1) :: isperiodic = .true.
    ! !> energy and pressure calculation stuff
    ! logical(4) :: isIntra
    ! integer*4, dimension(0:2) :: dn
    ! integer*4 :: i, j, k, iat, jat, ich, jch, icell, jcell, iiat, jjat
    ! integer*4 :: hi, hj, hk, ii, jj
    ! real*8, dimension(0:2) :: rij, f, dx
    ! real*8, dimension(0:2) :: cmi, cmj
    ! real*8, dimension(0:2,0:2) :: astrij, mstrij
    ! real*8 :: dijsq, dij, uij, duij
    ! real*8, allocatable, dimension (:,:) :: fij, fcm
    ! !> cell neighbor stuff
    ! integer*4 :: ncc                            ! number of neighbor cells
    ! integer*4, dimension(27) :: cc              ! neighbor cells indexes
    ! integer*4 :: cnat                           ! cell number of atoms
    ! integer*4, pointer, dimension (:) :: cat    ! cell atoms
    ! real*8 :: utrunc, ftrunc

    ! real*8, dimension (0:2, 0:n-1) :: rtmp

    ! real*8 :: pid, ptrunc, pidtrunc

    ! allocate( fij(3,n), fcm(3,nch))

    ! ! initialize interactions and box partitioning for neighbor cells
    ! call lj_initialize(utrunc, ftrunc)
    ! rcut = lj_rcut()
    ! dn = int( (/ a(0), b(1), c(2) /) / rcut)
    ! where(dn < 3) dn = 3 ! at least a 3x3x3 partition
    ! call box%initialize( r0, a, b, c, isperiodic)

    ! rtmp = transpose(rw)
    ! call box%splitAndfill( dn(0), dn(1), dn(2), n, rtmp)

    ! enrg = 0.d0
    ! einter = 0.d0
    ! eintra = 0.d0
    ! eimage = 0.d0
    ! fij = 0.d0
    ! fcm = 0.d0
    ! astrij = 0.d0
    ! mstrij = 0.d0

    ! do iat = 1, n ! NOTE box lists are indexed from 1 (std fortran)
    !     icell =  box%getAtomCell(iat)
    !     iiat = iat-1
    !     ich = atch(iiat)
    !     call box%getNeighbourCells(icell,ncc,cc)
    !     do j = 1, 27
    !         jcell = cc(j)
    !         if ( jcell .eq. -1) cycle
    !         call box%getCellAtoms(jcell, cnat, cat)
    !         do k = 1, cnat
    !             jat = cat(k)
    !             if ( jat .ge. iat) cycle
    !             jjat = jat-1
    !             jch = atch(jjat)
    !             if ( ich .eq. jch .and. (iat - jat) .lt. 4) cycle

    !             rij = r(iiat,:) - r(jjat,:)
    !             call box%minImgToGet( rij, hi, hj, hk)
    !             dijsq = square(rij)
    !             dij = dsqrt(dijsq)
    !             rij = rij / dij

    !             if (dij .gt. rcut ) cycle

    !             isIntra = (ich .eq. jch .and. hi .eq. 0 .and. hj .eq. 0 .and. hk .eq. 0 )

    !             uij = lj_energy(dijsq, dij, duij)
    !             f = - rij * duij
    !             fij(1:3,iat) = fij(1:3,iat) + f(0:2)
    !             fij(1:3,jat) = fij(1:3,jat) - f(0:2)

    !             enrg = enrg + uij
    !             if ( isIntra) then
    !                 eintra = eintra + uij
    !             else
    !                 einter = einter + uij
    !                 if ( ich .eq. jch ) eimage = eimage + uij
    !                 ! molecular virial
    !                 cmi = rcm(ich,:) + rw(iiat,:) - r(iiat,:)
    !                 cmj = rcm(jch,:) + rw(iiat,:) - rij(:) - r(jjat,:)
    !                 do ii = 0, 2
    !                     mstrij(ii,ii) = mstrij(ii,ii) + (cmi(ii) - cmj(ii))*f(ii)
    !                     do jj = ii+1, 2
    !                         mstrij(ii,jj) = mstrij(ii,jj) + (cmi(ii) - cmj(ii))*f(jj)
    !                         mstrij(jj,ii) = mstrij(ii,jj)
    !                     end do
    !                 end do
    !                 ! atomic virial
    !                 do jj = 0, 2
    !                     do ii = 0, 2
    !                         astrij(ii,jj) = astrij(ii,jj) + rij(ii) * f(jj)
    !                     end do
    !                 end do
    !             endif

    !         end do
    !     end do
    ! end do

    ! pid = - nch / box%volume * temperature
    ! etrunc = utrunc * dble(n)*dble(n) / 2.d0 / box%volume
    ! ptrunc = ftrunc * dble(n)*dble(n) / 2.d0 / box%volume / box%volume
    ! pidtrunc  = pid + ptrunc

    ! mstten = 0.d0
    ! astten = 0.d0
    ! stcm = 0.d0

    ! do i = 0, 2
    !     mstten(i,i) = - mstrij(i,i) / box%volume + pidtrunc
    !     do j = i+1, 2
    !         mstten(i,j) = - mstrij(i,j) / box%volume
    !         mstten(j,i) = mstten(i,j)
    !     end do
    ! end do

    ! do iat = 0, n-1
    !     dx(:) = r(iat,:) - rcm(atch(iat),:)
    !     do j = 0, 2
    !         do i  =0, 2
    !             stcm(i,j) = stcm(i,j) + ( dx(i) * fij(j,iat) )
    !         end do
    !     end do
    ! end do

    ! do j = 0, 2
    !     do i = 0, 2
    !         astten(i,j) = ( stcm(i,j) - astrij(i,j) ) / box%volume
    !     end do
    ! end do

    ! do i = 0, 2
    !   astten(i,i) = astten(i,i) + pidtrunc
    ! end do

    ! press = (mstten(0,0)+mstten(1,1)+mstten(2,2))/3.d0
    ! apress = (astten(0,0)+astten(1,1)+astten(2,2))/3.d0

    ! deallocate( fij, fcm)

end subroutine peenrgpress

!> calculate conformational tensor for pe ua-model
!> (only ortho boxes are supported for now)
subroutine peconftens(n, r, nch, maxnatch, natch, chat , ct)
    use vector3d, only : length
    implicit none
    !> # of atoms, # of chains, maximum number of atoms per chain
    integer, intent(in) :: n, nch, maxnatch
    !> unwrapped coordinates
    real*8, dimension(0:n-1,0:2), intent(in) :: r
    !> # of atoms per chain
    integer, dimension(0:nch-1), intent(in) :: natch
    !> atoms per chain
    integer, dimension(0:nch-1, 0:maxnatch-1), intent(in) :: chat
    !> conformation tensor (mean value of molecular conformation
    !> tensors for the given configuration).
    real*8, dimension(0:2,0:2), intent(out) :: ct

    ! use the same reference as in polybead code
    real*8, parameter :: alpha0 =    9.13127122d0
    real*8, parameter :: alpha1 =  -75.1865516d0
    real*8, parameter :: alpha2 =  315.74216d0
    real*8, parameter :: alpha3 = -500.351884d0
    real*8, parameter :: bondlsq = 1.54d0*1.54d0

    integer :: ich, iat, jat
    real*8 :: ichnat
    real*8, dimension(0:2) :: ree
    real*8 :: sqree0

    ct(:,:) = 0.d0
    do ich = 0, nch-1
        iat = chat(ich,0)
        jat = chat(ich,natch(ich)-1)
        ree = r(iat,:) - r(jat,:) ! end-to-end vector using unwrpped coordinatess
        ichnat = dble( natch(ich)-1)

        sqree0 = (  alpha0 +                   &
     &              alpha1 / ichnat +          &
     &              alpha2 / ichnat / ichnat + &
     &              alpha3 / ichnat / ichnat / ichnat ) *  ichnat * bondlsq

        ct(0,0) = ct(0,0) + 3.d0*ree(0)*ree(0)/sqree0
        ct(0,1) = ct(0,1) + 3.d0*ree(0)*ree(1)/sqree0
        ct(0,2) = ct(0,2) + 3.d0*ree(0)*ree(2)/sqree0
        ct(1,1) = ct(1,1) + 3.d0*ree(1)*ree(1)/sqree0
        ct(1,2) = ct(1,2) + 3.d0*ree(1)*ree(2)/sqree0
        ct(2,2) = ct(2,2) + 3.d0*ree(2)*ree(2)/sqree0
    enddo
    ct(1,0) = ct(0,1)
    ct(2,0) = ct(0,2)
    ct(2,1) = ct(1,2)
    ct = ct /dble(nch)

    return
end subroutine peconftens

!> Accelerate wrapped coordinates check (ONLY cubic boxes)
!> This call fixes wrapped coordinates which are slightly out of the box
subroutine chkwrap(n, r, r0, a, b, c, rw)
    implicit none
      integer, intent(in) :: n
      real*8, dimension(0:n-1,0:2), intent(in) :: r
      real*8, dimension(0:2), intent(in) :: r0
      real*8, intent(in) :: a, b, c
      real*8, dimension(0:n-1,0:2), intent(out) :: rw

      integer :: i, j
      real*8, dimension(0:2) :: r1, d

      d = (/ a, b, c /)
      r1 = r0 + d

      ! TODO improve check accuracy

      do i = 0, n-1
        do j = 0, 2
            if ( r(i,j) < r0(j)) then
                rw(i,j) = r(i,j) + d(j)
            elseif ( r(i,j) > r1(j)) then
                rw(i,j) = r(i,j) - d(j)
            else
                rw(i,j) = r(i,j)
            endif
        end do
      end do

end subroutine chkwrap

!> Accelerate fast center of masses calculation
!> for the whole collection of molecules
subroutine fastcom(n, r, mass, molecule, m, cm)
    implicit none
    !> number of atoms
    integer, intent(in) :: n
    !> unwrapped coordinates
    real*8, dimension(0:n-1,0:2), intent(in) :: r
    !> atom mass
    real*8, dimension(0:n-1), intent(in) :: mass
    !> atoms molecule (zero indexed array)
    integer, dimension(0:n-1), intent(in) :: molecule
    !> number of molecules
    integer, intent(in) :: m
    !> center of mass array
    real*8, dimension(0:m-1,0:2), intent(out) :: cm

    integer :: i, im
    integer :: nmolecules
    real*8, dimension(0:m) :: mmass

    cm = 0.d0

    nmolecules = maxval( molecule) + 1
    if ( m /= nmolecules) return

    mmass = 0.d0
    do i = 0, n-1
        im = molecule(i)
        cm( im, :) = cm( im, :) + mass(i)*r(i,:)
        mmass( im) = mmass( im) + mass(i)
    end do

    do i = 0, m-1
        cm( i, :) = cm(i, :) / mmass(i)
    end do

end subroutine fastcom

!> Accelerate fast radious of gyratio calculation
!> (the given coordinates should be unwrapped)
subroutine fastrg(n, r, mass, molecule, m, excluded, cm, rg)
    implicit none
    integer, intent(in) :: n
    !> the coordinates should be unwrapped
    real*8, dimension(0:n-1,0:2), intent(in) :: r
    real*8, dimension(0:n-1), intent(in) :: mass
    integer, dimension(0:n-1), intent(in) :: molecule ! provided with zero starting indexes
    integer, intent(in) :: m
    logical, dimension(0:m-1), intent(in) :: excluded
    real*8, dimension(0:m-1,0:2), intent(out) :: cm
    real*8, dimension(0:m-1), intent(out) :: rg

    integer :: i, im
    integer :: nmolecules
    real*8, dimension(0:m) :: mmass
    real*8, dimension(0:2) :: dx
    real*8 :: sqdx

    ! LDP: q&d optimize it (exercise)
    cm = 0.d0
    rg = 0.d0

    nmolecules = maxval( molecule) + 1
    if ( m /= nmolecules) return

    mmass = 0.d0
    do i = 0, n-1
        im = molecule(i)
        if ( excluded(im)) cycle
        cm( im, :) = cm( im, :) + mass(i)*r(i,:)
        mmass( im) = mmass( im) + mass(i)
    end do

    do im = 0, m-1
        if ( excluded(im)) cycle
        cm( im, :) = cm(im, :) / mmass(im)
    end do

    do i = 0, n-1
        im = molecule(i)
        if ( excluded(im)) cycle
        dx = r(i,:) - cm(im,:)
        sqdx = dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2)
        rg( im) = rg( im) + mass(i)*sqdx
    end do

    do im = 0, m-1
        if ( excluded(im)) cycle
        rg( im) = rg( im) / mmass( im)
    end do

end subroutine fastrg

!> Accelerate fast center of masses calculation
!> for the whole collection of atoms
subroutine fastcom_total(n, r, mass, cm, mtot)
    implicit none
    integer, intent(in) :: n
    !> the coordinates should be unwrapped
    real*8, dimension(0:n-1,0:2), intent(in) :: r
    real*8, dimension(0:n-1), intent(in) :: mass
    real*8, dimension(0:2), intent(out) :: cm
    real*8, intent(out) :: mtot

    integer i

    cm = 0.d0
    mtot = 0.d0
    do i = 0, n-1
        cm = cm + mass(i)*r(i,:)
        mtot = mtot + mass(i)
    end do

    cm = cm / mtot

end subroutine fastcom_total

!> Accelerate fast center of masses calculation
subroutine fastcom_exclude(n, r, mass, molecule, exclude, m, cm)
    implicit none
    integer, intent(in) :: n
    !> the coordinates should be unwrapped
    real*8, dimension(0:n-1,0:2), intent(in) :: r
    real*8, dimension(0:n-1), intent(in) :: mass
    integer, dimension(0:n-1), intent(in) :: molecule ! provided with zero starting indexes
    logical, dimension(0:n-1), intent(in) :: exclude  ! atoms to be excluded
    integer, intent(in) :: m
    real*8, dimension(0:m-1,0:2), intent(out) :: cm

    integer :: i, im
    integer :: nmolecules
    real*8, dimension(0:m) :: mmass

    cm = 0.d0

    nmolecules = maxval( molecule) + 1
    if ( m /= nmolecules) return

    mmass = 0.d0
    do i = 0, n-1
        im = molecule(i)
        if ( exclude(i)) cycle
        cm( im, :) = cm( im, :) + mass(i)*r(i,:)
        mmass( im) = mmass( im) + mass(i)
    end do

    do i = 0, m-1
        !write(*,'(i5,3(x,f15.8)," - ",f15.8)') i, cm(i,0), cm(i,1), cm(i,2), mmass(i)
        cm( i, :) = cm(i, :) / mmass(i)
    end do

end subroutine fastcom_exclude

!> Accelerate bond lengths calculateion
! subroutine fastbonds(n, r, nb, b, lengths)
!     implicit none
!     integer, intent(in) :: n                          !< number of atoms
!     real*8, dimension(0:n-1,0:2), intent(in)  :: r    !< atoms unwrapped coordinates
!     integer, intent(in) :: nb                         !< number of bonds
!     integer, dimension(0:nb-1,0:1), intent(in) :: b   !< bond info b(i,:) = (type, atom1, atom2) for ith bond
!     real*8, dimension(0:nb-1), intent(out) :: lengths !< bond lengths

!     integer :: i
!     real*8, dimension(0:2) :: dr

!     do i = 0, nb-1
!       dr = r( b(i,0), :) - r( b(i,1), :)
!       lengths(i) = dsqrt( dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2))
!     end do

! end subroutine fastbonds

!> Accelerate bend angles calculation
! subroutine fastangles(n, r, na, a, angles)
!     implicit none
!     integer, intent(in) :: n                          !< number of atoms
!     real*8, dimension(0:n-1,0:2), intent(in) :: r     !< atoms unwrapped coordinates
!     integer, intent(in) :: na                         !< number of angles
!     integer, dimension(0:na-1,0:3), intent(in) :: a   !< angle info b(i,:) = (type, atom1, atom2, atom3) for ith angle
!     real*8, dimension(0:na-1), intent(out) :: angles  !< angle value

!     integer :: i
!     real*8, dimension(0:2) :: dr1, dr2
!     real*8 :: length1, length2
!     real*8 :: c
!     real*8 :: pi, RAD2DEG

!     real*8, parameter :: EPS = 1.0D-12

!     pi = dacos(-1.d0)
!     RAD2DEG = pi / 180.d0
!     angles = 0.D0

!     do i = 0, na-1
!       dr1 = r( a(i,0), :) - r( a(i,1), :)
!       dr2 = r( a(i,2), :) - r( a(i,1), :)
!       length1 = dsqrt( dr1(0)*dr1(0) + dr1(1)*dr1(1) + dr1(2)*dr1(2))
!       length2 = dsqrt( dr2(0)*dr2(0) + dr2(1)*dr2(1) + dr2(2)*dr2(2))
!       c = ( dr1(0)*dr2(0) + dr1(1)*dr2(1) + dr1(2)*dr2(2)) / length1 / length2
!       c = dabs( c)
!       if ( c > 1.D0) then
!           if ( c - 1.D0 < EPS) c = dsign( 1.D0, c)
!       endif
!       angles( i) = dacos( c)
!     end do

! end subroutine fastangles

!> Accelerate bend angles calculation
! subroutine fastdihedrals(n, r, vo, va, vb, vc ,nd, d, angles)

!     use vector3d
!     use domain3d

!     implicit none
!     integer, intent(in) :: n                          !< number of atoms
!     real*8, dimension(0:n-1,0:2), intent(in) :: r     !< atoms unwrapped coordinates
!     real*8, dimension(0:2), intent(in) :: vo, va, vb, vc
!     integer, intent(in) :: nd                         !< number of angles
!     integer, dimension(0:nd-1,0:3), intent(in) :: d   !< angle info b(i,:) = (type, atom1, atom2, atom3) for ith angle
!     real*8, dimension(0:nd-1), intent(out) :: angles  !< angle value

!     integer :: i
!     real*8, dimension(0:2) :: dr1, dr2, dr3
!     real*8, dimension(0:2) :: dr1xdr2, dr2xdr3
!     real*8 dih, s
!     real*8 :: pi, RAD2DEG

!     real*8, parameter :: EPS = 1.0D-12

!     pi = dacos(-1.d0)
!     RAD2DEG = pi / 180.d0
!     angles = 0.D0

!     do i = 0, nd-1
!         ! bond vectors
!         dr1 = r( d(i,1), :) - r( d(i,0), :)
!         dr2 = r( d(i,2), :) - r( d(i,1), :)
!         dr3 = r( d(i,3), :) - r( d(i,2), :)

!         ! cis zero righ-handed
!         dr1xdr2 = unit( cross(dr1, dr2))
!         dr2xdr3 = unit( cross(dr2, dr3))

!         dih = angle(dr1xdr2, dr2xdr3)
!         s = dot(dr1xdr2, dr3)

!         angles(i) = dsign(dih, s)

!     end do

! end subroutine fastdihedrals

!> get the size for the seed array (used to initialize the random_number fortran intrinsic
subroutine random_seed_size( n)
    implicit none
    !> the size of the integer array needed to initialize the random_number
    !> fortran intrinsic
    integer, intent(out) :: n
    call random_seed( size = n )
end subroutine random_seed_size

!> Accelerate the calculation of the local density using randomly inserted cell
!> (ONLY cubic boxes)
subroutine localdensity(n, r, mass, r0, a0, b0, c0, a, b, c, nseed, seed, np, d)

      implicit none

      !> number of atoms
      integer, intent(in) :: n
      !> atom wrapped coordinates
      real*8, dimension(0:n-1,0:2), intent(in) :: r
      !> atom mass
      real*8, dimension(0:n-1), intent(in) :: mass
      !> simulation box origin
      real*8, dimension(0:2), intent(in) :: r0
      !> simulation box dimensions (cubic)
      real*8, intent(in) :: a0, b0, c0
      !> box dimension to be used for the calculation of the local density
      !> (cubic)
      real*8, intent(in) :: a, b, c
      !> size of the randomo number generator seed. Should be obtained by calling
      !> the random_seed_size subroutine
      integer, intent(in) :: nseed
      !> random number generator seed it should be initialized prior to the call
      integer, dimension(0:nseed-1), intent(inout) :: seed
      !> number of local density samples
      integer, intent(in) :: np
      !> the sequence of local densities calculated
      real*8, dimension(0:np-1), intent(out) :: d

      integer :: i, j
      real*8 :: vol
      real*8, dimension(0:2) :: x0, x1, x2, x3, x4
      logical, dimension(0:n-1) :: chk

      real*8, parameter :: EPS = 1.D-12

      ! check
      if ( a > a0 .or. b > b0 .or. c > c0) then
          d = sum( mass) / ( a0 * b0 * c0)
          return
      endif

      ! initiailize random seed
      call random_seed( put = seed(0:nseed-1))

      x2 = (/ a, b, c /)
      x3 = (/ a0-a+EPS, b0-b+EPS, c0-c+EPS /)
      vol = a * b * c

      do i = 0, np-1
         ! get the random point x1 inside the cube [r0,r0+x3] = [r0,r0+(/a0-a,b0-b,c0-c)]
         ! to be used as the origing of the cube (a,b,c)
         call random_number( x4)
         x0 = r0 + x3 * x4
         x1 = x0 + x2
         !> check if the points are inside the cube
         forall(j=0:n-1) chk(j) = ALL( r(j,:) >= x0) .and. ALL(r(j,:) <= x1)
         d(i) = sum( mass, chk) / vol
      enddo

      call random_seed( get = seed(0:nseed-1))

end subroutine localdensity

!> Accelerate the calculation of the local density using a grid (ONLY cubic boxes)
subroutine grd_localdensity(n, r, mass, r0, a0, b0, c0, f, na, nb, nc, d)

      implicit none

      !> number of atoms
      integer, intent(in) :: n
      !> atom wrapped coordinates
      real*8, dimension(0:n-1,0:2), intent(in) :: r
      !> atom mass
      real*8, dimension(0:n-1), intent(in) :: mass
      !> simulation box origin
      real*8, dimension(0:2), intent(in) :: r0
      !> simulation box dimensions (cubic)
      real*8, intent(in) :: a0, b0, c0
      !> fraction of the original volume to be considered
      real*8, intent(in) :: f
      !> number grid points for each edge (cubic)
      integer, intent(in) :: na, nb, nc
      !> the sequence of local densities calculated
      real*8, dimension(0:na*nb*nc-1), intent(out) :: d

      integer :: i, ia, ib, ic, indx, nab, np
      real*8 ::vol, tmp
      real*8, dimension(0:2) :: x0, x1, fdelta, delta

      real*8, parameter :: EPS = 1.D-12

      ! check
      np = na * nb * nc

      ! calculate cubic box corners
      tmp = 0.5d0 * f**(1.d0/3.d0)
      ! calculate the grid on the cubic box
      x0 = r0 + (/ a0, b0, c0 /) * tmp
      x1 = r0 + (/ a0, b0, c0 /) * (1.d0 - tmp)
      fdelta = x1 - x0
      delta = (/ fdelta(0) / na, fdelta(1) / nb, fdelta(2) / nc /)
      vol = delta(0)*delta(1)*delta(2)

      ! initialize and calculate cells density
      d = 0.d0
      nab = na * nb
      do i = 0, n-1
         if ( ALL( r(i,:) >= x0) .and. ALL(r(i,:) <= x1) ) then
            ia = int( (r(i,0)-x0(0)) / delta(0))
            ib = int( (r(i,1)-x0(1)) / delta(1))
            ic = int( (r(i,2)-x0(2)) / delta(2))
            if ( ia == na) ia = na - 1
            if ( ib == nb) ib = nb - 1
            if ( ic == nc) ic = nc - 1
            indx = ic * nab + ib * na + ia ! cell index
            d(indx) = d(indx) + mass(i)
         end if
      enddo
      d = d / vol

end subroutine grd_localdensity

!> Accelerate gyration tensors calculations
subroutine gyration(n, r, mass, molecule, m, mexclude, rp, rg, eigval, eigvec, ierr)

    use vector3d

    implicit none

    integer, intent(in) :: n                          ! # atoms
    real*8, dimension(0:n-1,0:2), intent(in) :: r     ! unwrapped atom coordinates
    real*8, dimension(0:n-1), intent(in) :: mass      ! atom mass
    integer, dimension(0:n-1), intent(in) :: molecule ! atom molecule (zero starting indexes)
    integer, intent(in) :: m                          ! # molecules
    logical, dimension(0:m-1), intent(in) :: mexclude ! molecules to be excluded
    real*8, dimension(0:n-1,0:2), intent(out) :: rp   ! coordinates alligned to the invariant, primary (gyration tensor) frame
    real*8, dimension(0:m-1,0:5), intent(out) :: rg
    real*8, dimension(0:m-1,0:2), intent(out) :: eigval
    real*8, dimension(0:m-1,0:8), intent(out) :: eigvec
    integer, intent(out) :: ierr                       ! error flag (0 on normal exit)

    integer i, im, i1, i2, i3
    real*8, dimension(0:2) :: dr, dr_, x_, y_, z_
    real*8, dimension(0:m-1) :: mmass
    real*8, dimension(0:m-1,0:2) :: cm
    real*8, dimension(0:m-1,3,3) :: frame

    ! eigenvectors stuff
    character(LEN=1), parameter :: UPLO_ = 'U'
    integer, parameter :: N_ = 3, LDA_ = 3
    integer :: LWORK_, NB_, INFO_
    real*8 :: A_(LDA_,N_), W_(3)
    real*8, allocatable, dimension (:) :: WORK_
    integer, external :: ILAENV

    logical, dimension(0:n-1) :: exclude_ ! atoms to be excluded

    ! initialize
    ierr = 0
    rp = 0.d0
    rg = 0.d0
    eigval = 0.d0
    eigvec = 0.d0
    frame = 0.d0
    ! mexclude = .FALSE.

    ! calculate cm
    do i = 0, n-1
        exclude_(i) = mexclude( molecule(i))
    enddo
    call fastcom_exclude(n, r, mass, molecule, exclude_, m, cm)

    ! calculate gyration tensors
    mmass = 0.D0
    rg = 0.D0
    do i = 0, n-1
        im = molecule(i)
        if ( mexclude(im) ) cycle
        dr = r(i,:) - cm(im,:)
        rg(im,0) = rg(im,0) + mass(i)*dr(0)*dr(0) ! xx
        rg(im,1) = rg(im,1) + mass(i)*dr(1)*dr(1) ! yy
        rg(im,2) = rg(im,2) + mass(i)*dr(2)*dr(2) ! zz
        rg(im,3) = rg(im,3) + mass(i)*dr(0)*dr(1) ! xy
        rg(im,4) = rg(im,4) + mass(i)*dr(0)*dr(2) ! xz
        rg(im,5) = rg(im,5) + mass(i)*dr(1)*dr(2) ! yz
        mmass(im) = mmass(im) + mass(i)
    end do
    do im = 0, m-1
      if ( mexclude(im) ) cycle
      rg(im,:) = rg(im,:) / mmass(im)
    end do

    ! calculate eigen and specify the primary axes for each molecule
    NB_ = ILAENV( 1, 'DSYTRD',  UPLO_, N_, -1, -1, -1 )
    LWORK_ = MAX( 1, ( NB_+2 ) * N_ )
    allocate ( WORK_(LWORK_))

    do im = 0, m-1
        if ( mexclude(im) ) cycle
        A_(1,1) = rg(im,0)
        A_(2,2) = rg(im,1)
        A_(3,3) = rg(im,2)
        A_(1,2) = rg(im,3)
        A_(1,3) = rg(im,4)
        A_(2,3) = rg(im,5)
        CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
        if (INFO_ .ne. 0) then
            print*, 'fastpost.f90:gyration: problem in eigenvectors calculation.'
            ! mexclude(im) = .TRUE.
            ! or establish a local frame of convience
            ! find three atoms and use
            i1 = -1
            i2 = -1
            i3 = -1
            do i = 0, n-1, 3
                ! find three atoms form the molecule
                if (molecule(i) == im ) then
                    if ( i1 == -1 ) then
                        i1 = i
                    else if ( i2 == -1 ) then
                        i2 = i
                    else if ( i3 == -1 ) then
                        i3 = i
                    end if
                endif
                ! if found construct the orthorormal frame
                if (i1*i2*i3 .gt. 0) then
                    x_ = cm(im,:) - r(i1,:) ! x' -> (cm-i1)
                    call toUnit(x_)
                    A_(:,1) = x_

                    dr = r(i3,:) - r(i1,:)  ! y' -> (i2-i1) normal to x'
                    dr_ = project( dr, x_)
                    y_ = dr - dr_
                    call toUnit(y_)
                    A_(:,2) = y_

                    z_ = cross(x_, x_)       ! z' -> cross(x',y')
                    call toUnit(z_)
                    A_(:,3) = z_

                    exit
                endif
            end do
            ierr = 1
            if (i1*i2*i3 .gt. 0) then
                write(*,'(i6,x,f15.8)') im, mmass(im)
                write(*,'(3(x,f15.8))') rg(im,0), rg(im,3), rg(im,4)
                write(*,'(3(x,f15.8))') rg(im,3), rg(im,1), rg(im,5)
                write(*,'(3(x,f15.8))') rg(im,4), rg(im,5), rg(im,2)
                print*, "-------"
                A_(:,1) = (/1.,0.,0./)
                A_(:,2) = (/0.,1.,0./)
                A_(:,3) = (/0.,0.,1./)
            endif
        else
            frame(im,:,:) = A_
            eigval(im,0) = W_(3)
            eigval(im,1) = W_(2)
            eigval(im,2) = W_(1)
            eigvec(im,0:2) = A_(:,3)
            eigvec(im,3:5) = A_(:,2)
            eigvec(im,6:8) = A_(:,1)
        endif
    enddo
    deallocate( WORK_)

    ! express molecules coordinates in the local molecular frames
    ! defined from the gyration tensor
    do i = 0, n-1
        im = molecule(i)
        if ( mexclude(im) ) cycle
        dr = r(i,:) - cm(im,:)
        rp(i,2) = dr(0) * frame(im,1,1) + dr(1) * frame(im,2,1) + dr(2) * frame(im,3,1)
        rp(i,1) = dr(0) * frame(im,1,2) + dr(1) * frame(im,2,2) + dr(2) * frame(im,3,2)
        rp(i,0) = dr(0) * frame(im,1,3) + dr(1) * frame(im,2,3) + dr(2) * frame(im,3,3)
    enddo

end subroutine gyration

!> Accelerate inertia tensors calculations
subroutine inertia(n, r, mass, molecule, m, mexclude, inert, eigval, eigvec, ierr)

    use vector3d

    implicit none

    integer, intent(in) :: n                          ! # atoms
    real*8, dimension(0:n-1,0:2), intent(in) :: r     ! unwrapped atom coordinates
    real*8, dimension(0:n-1), intent(in) :: mass      ! atom mass
    integer, dimension(0:n-1), intent(in) :: molecule ! atom molecule (zero starting indexes)
    integer, intent(in) :: m                          ! # molecules
    logical, dimension(0:m-1), intent(in) :: mexclude ! molecules to be excluded
    real*8, dimension(0:m-1,0:5), intent(out) :: inert ! inertia tensor
    real*8, dimension(0:m-1,0:2), intent(out) :: eigval
    real*8, dimension(0:m-1,0:8), intent(out) :: eigvec
    integer, intent(out) :: ierr                       ! error flag (0 on normal exit)

    integer i, im
    real*8, dimension(0:2) :: dr
    real*8, dimension(0:m-1,0:2) :: cm

    ! eigenvectors stuff
    character(LEN=1), parameter :: UPLO_ = 'U'
    integer, parameter :: N_ = 3, LDA_ = 3
    integer :: LWORK_, NB_, INFO_
    real*8 :: A_(LDA_,N_), W_(3)
    real*8, allocatable, dimension (:) :: WORK_
    integer, external :: ILAENV

    logical, dimension(0:n-1) :: exclude_ ! atoms to be excluded

    ! initialize
    ierr = 0
    inert = 0.d0
    eigval = 0.d0
    eigvec = 0.d0
    ! mexclude = .FALSE.

    ! calculate cm
    do i = 0, n-1
        exclude_(i) = mexclude( molecule(i))
    enddo
    call fastcom_exclude(n, r, mass, molecule, exclude_, m, cm)

    ! calculate inertia tensors
    inert = 0.D0
    do i = 0, n-1
        im = molecule(i)
        if ( mexclude(im) ) cycle
        dr = r(i,:) - cm(im,:)
        inert(im,0) = inert(im,0) + mass(i)*(dr(1)*dr(1)+dr(2)*dr(2)) ! xx
        inert(im,1) = inert(im,1) + mass(i)*(dr(0)*dr(0)+dr(2)*dr(2)) ! yy
        inert(im,2) = inert(im,2) + mass(i)*(dr(0)*dr(0)+dr(1)*dr(1)) ! zz
        inert(im,3) = inert(im,3) - mass(i)*dr(0)*dr(1) ! xy
        inert(im,4) = inert(im,4) - mass(i)*dr(0)*dr(2) ! xz
        inert(im,5) = inert(im,5) - mass(i)*dr(1)*dr(2) ! yz
    end do

    ! calculate eigen and specify the primary axes for each molecule
    NB_ = ILAENV( 1, 'DSYTRD',  UPLO_, N_, -1, -1, -1 )
    LWORK_ = MAX( 1, ( NB_+2 ) * N_ )
    allocate ( WORK_(LWORK_))

    do im = 0, m-1
        if ( mexclude(im) ) cycle
        A_(1,1) = inert(im,0)
        A_(2,2) = inert(im,1)
        A_(3,3) = inert(im,2)
        A_(1,2) = inert(im,3)
        A_(1,3) = inert(im,4)
        A_(2,3) = inert(im,5)
        CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
        if (INFO_ .ne. 0) then
            print*, 'fastpost.f90:inertia: problem in eigenvectors calculation.'
            IERR = 1
            return
        else
            eigval(im,0) = W_(3)
            eigval(im,1) = W_(2)
            eigval(im,2) = W_(1)
            eigvec(im,0:2) = A_(:,3)
            eigvec(im,3:5) = A_(:,2)
            eigvec(im,6:8) = A_(:,1)
        endif
    enddo
    deallocate( WORK_)

end subroutine inertia

!> Accelerate order parameter calculations
subroutine order_parameter( v, exclude, m, Q, eigval, eigvec, IERR)

    use vector3d

    implicit none

    real*8, dimension(0:m-1,0:2), intent(in) :: v       ! vectors to be considered
    logical, dimension(0:m-1), intent(in) :: exclude  ! vectors to be excluded
    integer, intent(in) :: m                          ! # vectors
    ! logical, dimension(0:m-1), intent(in) :: excludem
    real*8, dimension(0:8), intent(out) :: Q        ! order tensor
    real*8, dimension(0:2), intent(out) :: eigval   ! Q eigenvalues
    real*8, dimension(0:8), intent(out) :: eigvec   ! Q eigenvectors
    integer, intent(out) :: IERR

    integer im
    real*8, dimension(0:2) :: u

    ! eigenvectors stuff
    character(LEN=1), parameter :: UPLO_ = 'U'
    integer, parameter :: N_ = 3, LDA_ = 3
    integer :: LWORK_, NB_, INFO_
    real*8 :: A_(LDA_,N_), W_(3)
    real*8, allocatable, dimension (:) :: WORK_
    integer, external :: ILAENV

    IERR = 0    ! initialize

    ! calculate eigen and specify the primary axes for each molecule
    NB_ = ILAENV( 1, 'DSYTRD',  UPLO_, N_, -1, -1, -1 )
    LWORK_ = MAX( 1, ( NB_+2 ) * N_ )
    allocate ( WORK_(LWORK_))

    do im = 0, m-1
        if ( exclude(im)) continue
        u = unit( v(im,:))
        Q(0) = Q(0) + u(0)*u(0)
        Q(1) = Q(1) + u(1)*u(1)
        Q(2) = Q(2) + u(2)*u(2)
        Q(3) = Q(3) + u(0)*u(1)
        Q(4) = Q(4) + u(0)*u(2)
        Q(5) = Q(5) + u(1)*u(2)
    enddo

    Q =  3.d0 * Q / 2.d0 / m
    Q(0) = Q(0) - 0.5d0
    Q(1) = Q(1) - 0.5d0
    Q(2) = Q(2) - 0.5d0

    A_(1,1) = Q(0)
    A_(2,2) = Q(1)
    A_(3,3) = Q(2)
    A_(1,2) = Q(3)
    A_(1,3) = Q(4)
    A_(2,3) = Q(5)
    ! print*, A_(1,:)
    ! print*, A_(2,:)
    ! print*, A_(3,:)
    CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
    if (INFO_ .ne. 0) then
        IERR = 1
        print*, 'fastpost.f90:order_parameter: problem in order tensor eigenvectors calculation.'
        return
    else
        eigval(0) = W_(3)
        eigval(1) = W_(2)
        eigval(2) = W_(1)
        eigvec(0:2) = A_(:,3)
        eigvec(3:5) = A_(:,2)
        eigvec(6:8) = A_(:,1)
    endif

    deallocate( WORK_)

end subroutine order_parameter

subroutine order_parameter_local( v, nneighbors, neighbors, m, n, bin, nbins, q, h, IERR)

    use vector3d

    implicit none

    real*8, dimension(0:m-1,0:2), intent(in) :: v       ! vectors to be considered
    integer, dimension(0:m-1), intent(in) :: nneighbors
    integer, dimension(0:n-1), intent(in) :: neighbors
    integer, intent(in) :: m                          ! # of vectors
    integer, intent(in) :: n                          ! # of neighbors
    real*8, intent(in) :: bin                         ! local order histogram bin length
    real*8, dimension(0:nbins-1), intent(out) :: h      ! local order histogram
    integer, intent(in) :: nbins                      ! local order histogram number of bins
    ! logical, dimension(0:m-1), intent(in) :: excludem
    real*8, intent(out) :: q                        ! local order parameter
    integer, intent(out) :: IERR

    integer im, j, jm, n0, ibin
    real*8 :: c, ord, qm
    real*8, dimension(0:2,0:m-1) :: u
    real*8, dimension(0:2) :: iu
    integer, dimension(0:m-1) :: nneighbors0

    IERR = 0    ! initialize

    ! check thet total number of neighbors
    if ( sum(nneighbors) .ne. n ) then
        IERR = 1
        return
    endif
    ! find starting index of each molecule in the neighbors array
    nneighbors0 = 0
    u(:,0) = unit(v(0,:))
    do im = 1, m-1
        nneighbors0(im) = nneighbors0(im-1)+nneighbors(im-1)
        u(:,im) = unit(v(im,:))
    end do

    q = 0.d0
    do im = 0, m-1
        iu = u(:,im)
        n0 = nneighbors0(im)
        qm = 0.d0
        if ( nneighbors(im) .eq. 0 ) cycle
        do j = 0, nneighbors(im)-1
            jm = neighbors(n0+j)
            c = dot(iu,u(:,jm))
            ord = 0.5 * ( 3.d0 * c * c - 1.d0)
            qm = qm + ord
            ! update histogram
            ibin = int((ord + 0.5d0) / bin)
            if (ibin .gt. nbins-1) ibin = nbins-1
            h(ibin) = h(ibin)+1
        enddo
        q = q + qm / nneighbors(im)
    enddo
    q = q / dble(m)

end subroutine order_parameter_local

!> Accelerate order parameter calculations
subroutine order_parameter_inertia( n, r, mass, molecule, mexclude, m, Q, eigval, eigvec, IERR)

    use vector3d

    implicit none

    integer, intent(in) :: n                          ! # atoms
    real*8, dimension(0:n-1,0:2), intent(in) :: r     ! unwrapped atom coordinates
    real*8, dimension(0:n-1), intent(in) :: mass      ! atom mass
    integer, dimension(0:n-1), intent(in) :: molecule ! atom molecule (zero starting indexes)
    logical, dimension(0:m-1), intent(in) :: mexclude ! molecules to be excluded
    integer, intent(in) :: m                          ! # molecules
    ! logical, dimension(0:m-1), intent(in) :: excludem
    real*8, dimension(0:8), intent(out) :: Q        ! order tensor
    real*8, dimension(0:2), intent(out) :: eigval   ! Q eigenvalues
    real*8, dimension(0:8), intent(out) :: eigvec   ! Q eigenvectors
    integer, intent(out) :: IERR

    integer i, im
    logical, dimension(0:n-1) :: exclude_ ! atoms to be excluded
    real*8, dimension(0:2) :: dr
    real*8, dimension(0:m-1,0:2) :: cm      ! center of mass
    real*8, dimension(0:m-1,0:5) :: inert   ! inertia tensor
    real*8, dimension(0:2) :: u

    ! eigenvectors stuff
    character(LEN=1), parameter :: UPLO_ = 'U'
    integer, parameter :: N_ = 3, LDA_ = 3
    integer :: LWORK_, NB_, INFO_
    real*8 :: A_(LDA_,N_), W_(3)
    real*8, allocatable, dimension (:) :: WORK_
    integer, external :: ILAENV

    IERR = 0    ! initialize

    ! calculate cm
    do i = 0, n-1
        exclude_(i) = mexclude( molecule(i))
    enddo
    call fastcom_exclude(n, r, mass, molecule, exclude_, m, cm)

    ! calculate inertia tensors
    inert = 0.D0
    do i = 0, n-1
        im = molecule(i)
        if ( mexclude(im) ) cycle
        dr = r(i,:) - cm(im,:)
        inert(im,0) = inert(im,0) + mass(i)*(dr(1)*dr(1)+dr(2)*dr(2)) ! xx
        inert(im,1) = inert(im,1) + mass(i)*(dr(0)*dr(0)+dr(2)*dr(2)) ! yy
        inert(im,2) = inert(im,2) + mass(i)*(dr(0)*dr(0)+dr(1)*dr(1)) ! zz
        inert(im,3) = inert(im,3) - mass(i)*dr(0)*dr(1) ! xy
        inert(im,4) = inert(im,4) - mass(i)*dr(0)*dr(2) ! xz
        inert(im,5) = inert(im,5) - mass(i)*dr(1)*dr(2) ! yz
    end do

    ! calculate eigen and specify the primary axes for each molecule
    NB_ = ILAENV( 1, 'DSYTRD',  UPLO_, N_, -1, -1, -1 )
    LWORK_ = MAX( 1, ( NB_+2 ) * N_ )
    allocate ( WORK_(LWORK_))

    do im = 0, m-1
        ! if ( excludem(im) ) cycle
        A_(1,1) = inert(im,0)
        A_(2,2) = inert(im,1)
        A_(3,3) = inert(im,2)
        A_(1,2) = inert(im,3)
        A_(1,3) = inert(im,4)
        A_(2,3) = inert(im,5)
        CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
        if (INFO_ .ne. 0) then
            print*, 'fastpost.f90:order_parameter: problem in inertia tensor eigenvectors calculation.'
            IERR = 1
            return
        else
            u(:) = A_(:,1)       ! get the eigenvector consponds to the lower eigenvalue
            Q(0) = Q(0) + u(0)*u(0)
            Q(1) = Q(1) + u(1)*u(1)
            Q(2) = Q(2) + u(2)*u(2)
            Q(3) = Q(3) + u(0)*u(1)
            Q(4) = Q(4) + u(0)*u(2)
            Q(5) = Q(5) + u(1)*u(2)
        endif
    enddo

    Q =  3.d0 * Q / 2.d0 / m
    Q(0) = Q(0) - 0.5d0
    Q(1) = Q(1) - 0.5d0
    Q(2) = Q(2) - 0.5d0

    A_(1,1) = Q(0)
    A_(2,2) = Q(1)
    A_(3,3) = Q(2)
    A_(1,2) = Q(3)
    A_(1,3) = Q(4)
    A_(2,3) = Q(5)
    ! print*, A_(1,:)
    ! print*, A_(2,:)
    ! print*, A_(3,:)
    CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
    if (INFO_ .ne. 0) then
        IERR = 1
        print*, 'fastpost.f90:order_parameter: problem in order tensor eigenvectors calculation.'
        return
    else
        eigval(0) = W_(3)
        eigval(1) = W_(2)
        eigval(2) = W_(1)
        eigvec(0:2) = A_(:,3)
        eigvec(3:5) = A_(:,2)
        eigvec(6:8) = A_(:,1)
    endif

    deallocate( WORK_)

end subroutine order_parameter_inertia

!> Accelerate structure factor calculation
subroutine order_parameter_vectors(n, r,  m, P, director, mode, IERR)

    use vector3d

    implicit none

    integer, intent(in) :: n                          ! # atoms
    real*8, dimension(0:n-1,0:2), intent(in) :: r     ! unwrapped atom coordinates
    integer, intent(in) :: m                          ! # molecules
    real*8, intent(out) :: P                          ! order parameter
    real*8, dimension(0:2), intent(out) :: director   ! director
    character(LEN=2), intent(in), optional :: mode    ! 'ee' end-to-end, 'sq' backbone sequence
    integer, intent(out) :: IERR

    character(LEN=2) :: mode_
    integer :: nm, i, im, j
    real*8, dimension(0:5) :: Q     ! ordering tensor - symmetric
    real*8 :: drsq
    real*8, dimension(0:2) :: u, dr
    ! real*8, dimension(0:2) :: eigval
    ! real*8, dimension(0:2,0:2) :: eigvec

    ! eigenvectors stuff
    character(LEN=1), parameter :: UPLO_ = 'U'
    integer, parameter :: N_ = 3, LDA_ = 3
    integer :: LWORK_, NB_, INFO_
    real*8 :: A_(LDA_,N_), W_(3)
    real*8, allocatable, dimension (:) :: WORK_
    integer, external :: ILAENV

    IERR = 0    ! initialize
    mode_ = 'sq'
    if ( present(mode) ) mode_ = mode
    nm = n / m  ! number of atoms per molecule

    Q = 0.d0    ! calculate ordering tensor
    if (mode_ == 'ee') then
        do im = 0, m-1
                j=im*nm
                i=j+nm-1
                dr = r(i,:) - r(j,:)
                drsq = dsqrt( dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2))
                u = dr / drsq
                Q(0) = Q(0) + u(0)*u(0)
                Q(1) = Q(1) + u(1)*u(1)
                Q(2) = Q(2) + u(2)*u(2)
                Q(3) = Q(3) + u(0)*u(1)
                Q(4) = Q(4) + u(0)*u(2)
                Q(5) = Q(5) + u(1)*u(2)
        end do
    else if ( mode_ == 'sq' ) then
        do im = 0, m-1
            do i = 1, nm-1
                j = im * nm + i
                dr = r(j,:) - r(j-1,:)
                drsq = dsqrt( dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2))
                u = dr / drsq
                Q(0) = Q(0) + u(0)*u(0)
                Q(1) = Q(1) + u(1)*u(1)
                Q(2) = Q(2) + u(2)*u(2)
                Q(3) = Q(3) + u(0)*u(1)
                Q(4) = Q(4) + u(0)*u(2)
                Q(5) = Q(5) + u(1)*u(2)
            end do
        end do
    else
        print*, 'fastpost.f90:order_parameter_vector: argument mode is invalid (should be either "sq" or "ee").'
    endif
    Q =  3.d0 * Q
    Q(0) = Q(0) - 1.d0
    Q(1) = Q(1) - 1.d0
    Q(2) = Q(2) - 1.d0
    Q = Q / 2.d0 / ( m * ( nm-1))

    ! calculate eigen and specify the primary axes for each molecule
    NB_ = ILAENV( 1, 'DSYTRD',  UPLO_, N_, -1, -1, -1 )
    LWORK_ = MAX( 1, ( NB_+2 ) * N_ )
    allocate ( WORK_(LWORK_))

    A_ = 0.d0
    A_(1,1) = Q(0)
    A_(2,2) = Q(1)
    A_(3,3) = Q(2)
    A_(1,2) = Q(3)
    A_(1,3) = Q(4)
    A_(2,3) = Q(5)
    ! print*, A_(1,:)
    ! print*, A_(2,:)
    ! print*, A_(3,:)
    CALL DSYEV( 'V', 'U', N_, A_, LDA_, W_, WORK_, LWORK_, INFO_ )
    if (INFO_ .ne. 0) then
        IERR = 1
        print*, 'fastpost.f90:order_parameter_vector: problem in eigenvectors calculation.'
    else
        ! eigvec(0:2,0:2) = A_ ! eigenvectors are stored columnwise
        ! eigval(0) = W_(3)
        ! eigval(1) = W_(2)
        ! eigval(2) = W_(1)
        P = W_(3)
        director(0:2) = A_(:,3)
        ! print*, P
        ! print*, director
    endif

    deallocate( WORK_)

end subroutine order_parameter_vectors

!> Accelerate percistance length calculation
! subroutine percistance_length(n, r,  m, lp_acf, lp_flory, lp_curvature, IERR)

!     use vector3d
!     use, intrinsic :: iso_c_binding

!     implicit none

!     include 'fftw3.f03'
!     type(C_PTR) :: FPLAN, BPLAN

!     real(C_DOUBLE), pointer, dimension(:) :: FFTIN
!     complex(C_DOUBLE_COMPLEX), pointer, dimension(:) :: FFTOUT
!     type(C_PTR) :: p1, p2

!     integer, intent(in) :: n                                ! # points
!     real*8, dimension(0:n-1,0:2), intent(in) :: r           ! unwrapped points coordinates
!     integer, intent(in) :: m                                ! # batches
!     integer, intent(in) :: nm                               ! # points per batches
!     real*8, dimension(0:nm-1), intent(out) :: lp_acf        ! percistance length using tangent vectors auto-correlation function
!     real*8, dimension(0:nm-1), intent(out) :: lp_flory      ! percistance length usnig flory definition
!     real*8, dimension(0:nm-1), intent(out) :: lp_curvature  ! percistance length usnig curvature
!     integer, intent(out) :: IERR

!     real*8, dimension(3) :: re
!     real*8, dimension(1:n-1,3) :: u
!     real*8, dimension(1:n-1) :: l
!     integer :: i, j, k, kk, m

!     ! fft stuff use dynamic allocation
!     type(C_PTR) :: PLANF, PLANB
!     real(C_DOUBLE), pointer, dimension(:) :: FFTIN
!     complex(C_DOUBLE_COMPLEX), pointer, dimension(:) :: FFTOUT
!     type(C_PTR) :: p1, p2
!     real(C_DOUBLE) :: dimension(NM) :: FFTNORM

!     ! initialize fft
!     p1 = fftw_alloc_real(int(2*NM, C_SIZE_T))
!     call c_f_pointer(p1, IN, [NM])
!     p2 = fftw_alloc_complex(int(2*NM,C_SIZE_T))
!     call c_f_pointer(p2, OUT, [NM])
!     forall (i=0:NM-1) FFTNORM(i) = NM - i

!     ! create the tangent vectors (unit)
!     ! store the end-to-end vector to the first slot of each batch
!     do i = 0, m-1
!         k = i * nm
!         u(k,:) = r(k+nm-1,:) - r(k,:)
!         do j = 1, nm
!             kk = k+j
!             u(kk,:) = r(kk,:) - r(kk-1,:)
!             l = length( u(kk,:))
!             u(kk,:) = u(kk,:) / l
!         enddo
!     enddo

!     ! calculate bond/bond acf using fft
!     call dfftw_plan_dft_r2c_1d( PLANF, NM, FFTIN, FFTOUT, FFTW_ESTIMATE) ! forward
!     call dfftw_plan_dft_c2r_1d( PLANB, NM, FFTOUT, FFTIN, FFTW_ESTIMATE) ! backward

!     do i = 0, m-1
!         j = i * m + 1
!         k = j + m - 1
!         do kk = 1, 3
!             FFTIN=u(j:k,kk)
!             call dfftw_execute( PLANF)
!             FFTOUT(1:NM) = FFTOUT(1:NM)*conjg(FFTOUT(1:NM))
!             call dfftw_execute( PLANB)
!             FFTOUT = FFTOUT / N    !<- normalization missing from fftw (but not in np.fft.ifft)
!             FFTOUT = FFTOUT / NORM
!             lp_acf = lp_acf + FFTOUT
!         enddo
!         lp_flory(i) = lp_flory(i) + dot(u(j,:),u(j-1,:))
!     enddo
!     lp_acf = lp_acf / m
!     lp_flory = lp_flory / m

!     ! free fftw stuff
!     call fftw_free(p1)
!     call fftw_free(p2)
!     call fftw_destroy_plan(PLANF)
!     call fftw_destroy_plan(PLANB)

! end subroutine
