module lj_interactions

    implicit none

    private

    real*8, parameter :: epsreal = 0.091402d0, sigreal = 3.95d0
    real*8, parameter :: sigsq = sigreal * sigreal
    real*8, parameter :: frtrunc = 1.45d0
    real*8, parameter :: frcut   = frtrunc + 0.85d0
    real*8, parameter :: sig6    = sigreal ** 6
    real*8, parameter :: sig12   = sigreal ** 12    
    real*8, parameter :: anb     = 4.d0 * epsreal * sig12
    real*8, parameter :: cnb     = 4.d0 * epsreal * sig6
    real*8, parameter :: danb    = -12.d0 * anb
    real*8, parameter :: dcnb    =   6.d0 * cnb

    real*8, parameter :: r1nb    = 1.45d0 * sigreal
    real*8, parameter :: r1nbsq  = r1nb * r1nb
    
    real*8, parameter :: r2nb    = (1.45d0 + 0.85d0) * sigreal
    real*8, parameter :: r2nbsq  = r2nb * r2nb
    real*8, parameter :: r12nb   = r2nb - r1nb

    real*8, parameter :: shell     = 1.3d0
    real*8, parameter :: r2nbp     = (r2nb + shell)**2
    real*8, parameter :: fdel      = (r2nb - r1nb) / sigreal
    real*8, parameter :: U1        = (anb/r1nb**12) - (cnb/r1nb**6)
    real*8, parameter :: U1p       = -12.d0 * (anb/r1nb**13) +  6.d0 * (cnb/r1nb**7)
    real*8, parameter :: U1pp      = 156.d0 * (anb/r1nb**14) - 42.d0 * (cnb/r1nb**8)
    real*8, parameter :: c1        = U1
    real*8, parameter :: c2        = fdel * U1p * sigreal
    real*8, parameter :: c3        = 0.5d0 * (fdel*fdel) * U1pp * (sigreal*sigreal)
    real*8, parameter :: pi = 3.14159265358979323846D0
    real*8, dimension(0:5) :: coeff
    real*8 :: ILJ, Ispl, fILJ, fIspl
    real*8 :: truncscal, ftruncscal, small

    interface lj_initialize
        module procedure m_initialize
    end interface

    interface lj_energy
        module procedure m_energy
    end interface

    interface lj_rcut
        module procedure m_getRcut
    end interface

    interface d_peenrgpress
        module procedure m_peenrgpress
    end interface

    public lj_initialize, lj_energy, lj_rcut, d_peenrgpress

    contains

    function m_getRcut() result(rcut)
        real*8 :: rcut 
        rcut = r2nb
    end function m_getRcut

    subroutine m_initialize(utrunc, ftrunc)
        real*8, intent(out) :: utrunc, ftrunc
 
        coeff(0) = c1
        coeff(1) = c2
        coeff(2) = c3
        coeff(3) = -(10.d0 * c1) - (6.d0 * c2) - (3.d0 * c3)
        coeff(4) =  (15.d0 * c1) + (8.d0 * c2) + (3.d0 * c3)
        coeff(5) = -(6.d0  * c1) - (3.d0 * c2) -         c3
        ! Integral of the Lennard-Jones from Rtrunc to rmaxnb:
        ILJ = (- (sig12/9.d0) * (1.d0/r2nb**9 - 1.d0/r1nb**9)  &
               + (sig6 /3.d0) * (1.d0/r2nb**3 - 1.d0/r1nb**3)) &
               * 4.d0 * epsreal
        ! Integral of the quintic spline from Rtrunc to rmaxnb:
        Ispl =  ( coeff(0) + coeff(1)/2.d0 + coeff(2)/3.d0             &
                + coeff(3)/4.d0 + coeff(4)/5.d0 + coeff(5)/6.d0) *     &
                    (r2nb - r1nb)     * r1nb**2 +                      &
                ( coeff(0)/2.d0 + coeff(1)/3.d0 + coeff(2)/4.d0        &
                + coeff(3)/5.d0 + coeff(4)/6.d0 + coeff(5)/7.d0) *     &
                    (r2nb - r1nb)**2  * r1nb * 2.d0 +                  &
                ( coeff(0)/3.d0 + coeff(1)/4.d0 + coeff(2)/5.d0        &
                + coeff(3)/6.d0 + coeff(4)/7.d0 + coeff(5)/8.d0) *     &
                    (r2nb - r1nb)**3
        ! Integral of the Lennard-Jones from Rtrunc to rmaxnb:
        fILJ =  12.d0 * anb * (1.d0/r2nb**9 - 1.d0/r1nb**9) / 9.d0     &   
              -  2.d0 * cnb * (1.d0/r2nb**3 - 1.d0/r1nb**3)
        ! Integral of the force quintic spline from Rtrunc to rmaxnb:
        fIspl = (      coeff(1)/4.d0 + 2.d0*coeff(2)/5.d0 +              &
                  3.d0*coeff(3)/6.d0 + 4.d0*coeff(4)/7.d0 +              &
                  5.d0*coeff(5)/8.d0 ) * r12nb**3                        &
            +  (       coeff(1)/3.d0 + 2.d0*coeff(2)/4.d0 +              &
                   3.d0*coeff(3)/5.d0 + 4.d0*coeff(4)/6.d0 +             &
                   5.d0*coeff(5)/7.d0 ) * r12nb**2 * r1nb * 3.d0         &
            +  (        coeff(1)/2.d0 + 2.d0*coeff(2)/3.d0 +             &
                    3.d0*coeff(3)/4.d0 + 4.d0*coeff(4)/5.d0 +            &
                    5.d0*coeff(5)/6.d0 ) * r12nb * r1nbsq * 3.d0         &
            +  (         coeff(1)      +      coeff(2)      +            &
                         coeff(3)      +      coeff(4)      +            &
                         coeff(5)      ) * r1nb**3
     
        truncscal  = 4.d0 * pi * ( (ILJ - Ispl) - (cnb/3.d0) / (r2nb**3) )
        ftruncscal = 4.d0 * pi * ( (fILJ - fIspl) + cnb * 2.d0 / r2nb**3 ) / 3.d0
        small      = sigsq / 8.d0

        utrunc = truncscal
        ftrunc = ftruncscal

    end subroutine m_initialize

    function m_energy(rijsq, rij, duij) result(uij)
        real*8, intent(in) :: rijsq, rij
        real*8 :: uij, duij
        real*8 :: r12, r6
        real*8 :: ratiosq, fpsi

        if      (rijsq.lt.small) then
            uij = 1.d5
            duij = 1.d5
        else if (rijsq.lt.r1nbsq) then
            ratiosq  = sigsq / rijsq
            r6 =  ratiosq **3
            r12 = r6*r6
            uij = 4.d0 * epsreal * (r12 - r6)
            duij = 4.d0 * epsreal * (-12.d0*r12+6.0*r6)/rij
        else if (rijsq.lt.r2nbsq) then
            fpsi = (rij - r1nb) / r12nb
            uij = coeff(0) + fpsi * (  &
                  coeff(1) + fpsi * (  &
                  coeff(2) + fpsi * (  &
                  coeff(3) + fpsi * (  &
                  coeff(4) + fpsi * (  &
                  coeff(5))))))
            duij = (      coeff(1) + fpsi * (     &
                   2.d0 * coeff(2) + fpsi * (     &
                   3.d0 * coeff(3) + fpsi * (     &
                   4.d0 * coeff(4) + fpsi * (     &
                   5.d0 * coeff(5)))))) / r12nb                 
        else
            uij = 0.d0
            duij = 0.d0
        end if
    end function m_energy


    subroutine m_peenrgpress(n, r, rw, nch, rcm, atch , a, b, c, temperature,     &
            rcut, enrg, einter, eintra, eimage, etrunc,                  &
            press, apress, pid, ptrunc, mstten, astten)
    !enrg,einter,eintra,eimage,press,apress,mstens[:,:],astens[:,:] 
    !    = peenrgpress(r,rw,rcm,molecules,box.va, box.vb, box.vc, T, rcut)
    use vector3d
    use domain3d
    ! use lj_interactions

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
    real*8, dimension(0:2,0:2) :: stcm

    !> box stuff
    type(SimulationBox) :: box
    type(BoxCells) :: cells
    real*8, dimension(0:2) :: r0 = (/ 0.d0, 0.d0, 0.d0/)
    logical(4) :: isperiodic = .true.
    !> energy and pressure calculation stuff
    logical(4) :: isIntra
    integer*4, dimension(0:2) :: dn
    integer*4 :: i, j, k, iat, jat, ich, jch, icell, jcell, iiat, jjat
    integer*4 :: hi, hj, hk, ii, jj
    real*8, dimension(0:2) :: rij, f, dx
    real*8, dimension(0:2) :: cmi, cmj
    real*8, dimension(0:2,0:2) :: astrij, mstrij    
    real*8 :: dijsq, dij, uij, duij
    real*8, allocatable, dimension (:,:) :: fij, fcm
    !> cell neighbor stuff
    integer*4 :: ncc                            ! number of neighbor cells
    integer*4, dimension(27) :: cc              ! neighbor cells indexes
    integer*4 :: cnat                           ! cell number of atoms
    integer*4, pointer, dimension (:) :: cat    ! cell atoms
    real*8 :: utrunc, ftrunc

    !> convert ideal pressure term (nRT/V) to simulation units kcal/mol/A3 (Kb*NA*4.184/1000.0)
    real*8, parameter :: pidcnv = 1.38064852*6.02214076/4.184d0/1000.d0
    real*8, dimension (0:2, 0:n-1) :: rtmp

    real*8 :: pidtrunc

    allocate( fij(3,n), fcm(3,nch))

    ! initialize interactions and box partitioning for neighbor cells
    call lj_initialize(utrunc, ftrunc)
    rcut = lj_rcut()
    dn = int( (/ a(0), b(1), c(2) /) / rcut)
    where(dn < 3) dn = 3 ! at least a 3x3x3 partition
    call box%initialize( r0, a, b, c, isperiodic)

    rtmp = transpose(rw)
    call cells%splitAndfill( box, dn(0), dn(1), dn(2), n, rtmp)

    enrg = 0.d0
    einter = 0.d0
    eintra = 0.d0
    eimage = 0.d0
    fij = 0.d0
    fcm = 0.d0
    astrij = 0.d0
    mstrij = 0.d0

    do iat = 1, n ! NOTE box lists are indexed from 1 (std fortran)
        icell =  cells%getAtomCell(iat)
        iiat = iat-1
        ich = atch(iiat)
        call cells%getNeighbourCells(icell,ncc,cc)
        do j = 1, 27
            jcell = cc(j)
            if ( jcell .eq. -1) cycle
            call cells%getCellAtoms(jcell, cnat, cat)
            do k = 1, cnat
                jat = cat(k)
                if ( jat .ge. iat) cycle
                jjat = jat-1
                jch = atch(jjat)
                if ( ich .eq. jch .and. (iat - jat) .lt. 4) cycle

                rij = r(iiat,:) - r(jjat,:) 
                call box%minImgToGet( rij, hi, hj, hk)
                dijsq = square(rij)
                dij = dsqrt(dijsq)

                if (dij .gt. rcut ) cycle

                isIntra = (ich .eq. jch .and. hi .eq. 0 .and. hj .eq. 0 .and. hk .eq. 0 )

                uij = lj_energy(dijsq, dij, duij)
                f = - rij * duij / dij
                fij(1:3,iat) = fij(1:3,iat) + f(0:2)
                fij(1:3,jat) = fij(1:3,jat) - f(0:2)

                enrg = enrg + uij                 
                if ( isIntra) then
                    eintra = eintra + uij
                else
                    einter = einter + uij
                    if ( ich .eq. jch ) eimage = eimage + uij
                    ! molecular virial
                    cmi = rcm(ich,:) + rw(iiat,:) - r(iiat,:)
                    cmj = rcm(jch,:) + rw(iiat,:) - rij(:) - r(jjat,:)
                    do ii = 0, 2
                        mstrij(ii,ii) = mstrij(ii,ii) + (cmi(ii) - cmj(ii))*f(ii)
                        do jj = ii+1, 2
                            mstrij(ii,jj) = mstrij(ii,jj) + (cmi(ii) - cmj(ii))*f(jj)
                            mstrij(jj,ii) = mstrij(ii,jj)
                        end do
                    end do
                    ! atomic virial
                    do ii = 0, 2
                        do jj = 0, 2
                            astrij(ii,jj) = astrij(ii,jj) + rij(ii) * f(jj)
                        end do
                    end do
                endif
            end do
        end do
    end do

    pid = nch / box%volume * temperature * pidcnv
    etrunc = utrunc * dble(n)*dble(n) / 2.d0 / box%volume
    ptrunc = ftrunc * dble(n)*dble(n) / 2.d0 / box%volume / box%volume
    pidtrunc  = pid + ptrunc

    mstten = 0.d0
    astten = 0.d0
    stcm = 0.d0

    do i = 0, 2
        mstten(i,i) = - mstrij(i,i) / box%volume - pid
        do j = i+1, 2
            mstten(i,j) = - mstrij(i,j) / box%volume
            mstten(j,i) = mstten(i,j)
        end do
    end do

    do iat = 1, n
        iiat = iat -1
        dx(:) = r(iiat,:) - rcm(atch(iiat),:)
        do i = 0, 2
            do j = 0, 2
                stcm(i,j) = stcm(i,j) + ( dx(i) * fij(j,iat) )
            end do
        end do
    end do

    do i = 0, 2
        do j = 0, 2
            ! LDP: should be ( stcm(i,j) - astrij(i,j) ) / box%volume + ptrunc
            astten(i,j) = - ( stcm(i,j) - astrij(i,j) ) / box%volume
        end do
    end do

    ! do i = 0, 2
    !     astten(i,i) = astten(i,i) + ptrucn
    ! end do

    press = (mstten(0,0)+mstten(1,1)+mstten(2,2))/3.d0
    apress = (astten(0,0)+astten(1,1)+astten(2,2))/3.d0

    deallocate( fij, fcm)

    end subroutine m_peenrgpress

end module lj_interactions