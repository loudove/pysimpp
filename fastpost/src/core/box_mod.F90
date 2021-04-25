module simulation_box

    implicit none

    private

    real*8, parameter :: EPS = 1.d-12    ! precision used

    type CellAtomsData
        integer :: n = 0
        integer :: start = 0
    end Type CellAtomsData

    type, public :: SimulationBox
        !> true if cell is periodic
        logical(4) :: periodic
        !> flag indicating that all the arrays needed for
        !> the calculations are allocated correctly.
        logical :: initialized=.false.
        !> cell lengths
        real(8) :: a, b, c
        !> cell angles
        real(8) :: alpha, beta, gama
        !> cell origin
        real(8), dimension (3) :: r0
        !> cell vectors
        real(8), dimension (3) :: ux, uy, uz
        !> cell normal vectors (unit)
        real(8), dimension (3) :: nXY, nYZ, nZX
        !> cell normal vector "projection" lengths
        real(8) :: px, py, pz
        !> cell surfaces
        real(8) :: sXY, sYZ, sZX
        !> cell volume
        real(8) :: volume
        !> cell transformation matrix from fractional to Cartesian coordinates.
        real(8), dimension (6) :: matrix
        !> # cells in x dimension
        integer, private :: nxCells = 0
        !> # cells in y dimension
        integer, private :: nyCells = 0
        !> # cells in z dimension
        integer, private :: nzCells = 0
        !> # cells in xy plane
        integer, private :: nxyCells = 0
        !> # cells
        integer, private :: nCells = 0
        !> cells length and reverse length in x dimension
        real(8), private :: dx = 0.d0, rdx = 0.d0
        !> cells length and reverse length in y dimension
        real(8), private :: dy = 0.d0, rdy = 0.d0
        !> cells length and reverse length in z dimension
        real(8), private :: dz = 0.d0, rdz = 0.d0
        !> number of atoms
        integer, private :: nAtoms = 0
        !> cell that atom belongs
        integer, private, allocatable, dimension (:)   :: atomCell
        !> number of points
        integer, private :: nPoints = 0
        !> cell that point belong
        integer, private, allocatable, dimension (:)   :: pointCell

        !< when the box grid is filled the points are placed in each cell in
        !< ordered fashion {xmin->xmax, ymin->ymax, zmin->zmax} i.e. at each
        !< cell the first atoms will be closer to X0, Y0, and Z0 planes. This is
        !< convenient if you want to calculate periodic indexes conforming to
        !< specific direction which is useful if for example, you analyze Voronoi
        !< cells vertices and you want to find the periodic indexes of each vertex
        !< in order to use it for cluster connectivity checks.
        logical, private :: placeordered = .False.

        !> an array contains atoms indexes. Helps in avoiding the
        !> allocation/deallocation of memory. The slab:
        !> atoms(cellAtoms(i)%start:cellAtoms(i)%start+cellAtoms(i)%n)
        !> contains the atoms indexes in the ith cell
        integer, private, allocatable, dimension (:) :: atoms
        !> a array contains point indexes. Helps in avoiding the allocation
        !> deallocation of memory
        integer, private, allocatable, dimension (:) :: points

        type(CellAtomsData), private, allocatable, dimension(:) :: cellAtoms
        type(CellAtomsData), private, allocatable, dimension(:) :: cellPoints

        contains
        procedure :: setOrdered => m_setOrdered
        procedure :: initialize => m_initialize
        procedure :: minImg => m_MinimumImage
        procedure :: minImgTo => m_MinimumImageTo
        procedure :: minImgGet => m_iMinimumImage
        procedure :: minImgToGet => m_iMinimumImageTo
        procedure :: setMatrix => m_setMatrix
        procedure :: split => m_split
        procedure :: fill => m_fill
        procedure :: fillPoints => m_fillPoints
        procedure :: fillCells => m_fillCells
        procedure :: splitAndfill => m_splitAndfill
        procedure :: getAtomCell => m_getAtomCell
        procedure :: getCellAtoms => m_getCellAtoms
        procedure :: getCellPoints => m_getCellPoints
        procedure :: getNeighbourCells => m_getNeighbourCells
        procedure :: toCartesian => m_toCartesian
        procedure :: toFractional => m_toFractional
        procedure :: getProjections => m_getProjections
        procedure :: getVectors => m_getVectors
        procedure :: destroy => m_destroy

    end Type SimulationBox

    contains

    ! initialize the placeordered flag
    subroutine m_setOrdered( this, on)
        class(SimulationBox), intent(inout) :: this
        logical, intent(in) :: on
        this%placeordered = on
    end subroutine m_setOrdered

    ! initialize the box
    subroutine m_initialize( this, r, x, y, z, isPeriodic)
        use vector3d, only : length, angle, cross, dot
        class(SimulationBox), intent(inout) :: this
        real(8), dimension(3), intent(in) :: r, x, y, z
        logical(4), intent(in) :: isperiodic

        ! initialize variables
        this%nAtoms=0

        ! copy the values.
        this%r0 = r
        this%ux = x
        this%uy = y
        this%uz = z
        this%periodic = isPeriodic

        ! find lengths.
        this%a = length(x)
        this%b = length(y)
        this%c = length(z)

        ! find angles.
        this%gama  = angle(x,y)
        this%alpha = angle(y,z)
        this%beta = angle(z,x)

        ! find normal vectors.
        this%nXY = cross(x,y)
        this%nYZ = cross(y,z)
        this%nZX = cross(z,x)

        ! find surfaces.
        this%sXY = length(this%nXY)
        this%sYZ = length(this%nYZ)
        this%sZX = length(this%nZX)

        ! find volume.
        this%volume = dot(z,this%nXY)

        ! set normal vectors to unit.
        this%nXY = this%nXY * (1.d0/this%sXY)
        this%nYZ = this%nYZ * (1.d0/this%sYZ)
        this%nZX = this%nZX * (1.d0/this%sZX)

        ! find normal vectors projection lengths.
        this%px = this%sYZ / this%volume
        this%py = this%sZX / this%volume
        this%pz = this%sXY / this%volume

        return
    end subroutine m_initialize

    function m_MinimumImage(this, v) result(u)
      class(SimulationBox), intent(in) :: this
      real(8), dimension(3), intent(in) :: v
      real(8), dimension(3) :: u

      u = v
      call m_MinimumImageTo(this, u)

      return
    end function m_MinimumImage

    subroutine m_MinimumImageTo(this, v)
        use vector3d, only : dot
        class(SimulationBox), intent(in) :: this
        real(8), dimension(3), intent(inout) :: v
        real(8) project
        real(8) d

        if ( .NOT. this%periodic ) return
        project = dot(v,this%nYZ)
        d = dnint(project*this%px)
        v = v - d * this%ux
        project = dot(v,this%nZX)
        d = dnint(project*this%py)
        v = v - d * this%uy
        project = dot(v,this%nXY)
        d = dnint(project*this%pz)
        v = v - d * this%uz

        return
    end subroutine m_MinimumImageTo

    function m_iMinimumImage(this, v, i, j, k) result(u)
        class(SimulationBox), intent(in) :: this
        real(8), dimension(3), intent(in) :: v
        real(8), dimension(3) :: u
        integer :: i, j, k

        u=v
        call m_iMinimumImageTo(this,u,i,j,k)

        return
    end function m_iMinimumImage

    subroutine m_iMinimumImageTo(this, v, i, j, k)
        use vector3d, only : dot
        class(SimulationBox), intent(in) :: this
        real(8), dimension(3), intent(inout) :: v
        integer, intent(out) :: i, j, k
        real(8) project
        real(8) d

        if ( .NOT. this%periodic ) return
        project = dot(v,this%nYZ)
        d = dnint(project*this%px)
        i = int(d)
        v = v - d * this%ux
        project = dot(v,this%nZX)
        d = dnint(project*this%py)
        j = int(d)
        v = v - d * this%uy
        project = dot(v,this%nXY)
        d = dnint(project*this%pz)
        k = int(d)
        v = v - d * this%uz

        return
    end subroutine m_iMinimumImageTo

    subroutine m_setMatrix(this)
        class(SimulationBox), intent(inout) :: this

        this%matrix(1)=this%a;
        this%matrix(2)=this%b*dcos(this%gama);
        this%matrix(3)=this%c*dcos(this%beta);
        this%matrix(4)=dsqrt(this%b*this%b-this%matrix(2)*this%matrix(2));
        this%matrix(5)=(this%b*this%c*dcos(this%alpha)-this%matrix(2)*this%matrix(3))/this%matrix(4);
        this%matrix(6)=dsqrt(this%c*this%c-this%matrix(3)*this%matrix(3)-this%matrix(5)*this%matrix(5));

        return
    end subroutine m_setMatrix

    subroutine m_split(this, nx, ny, nz)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz

        if ( nx < 3 .AND. ny < 3 .AND. nz < 3 ) return

        this%nxCells = nx
        this%nyCells = ny
        this%nzCells = nz
        this%dx = 1.d0 / this%nxCells
        this%dy = 1.d0 / this%nyCells
        this%dz = 1.d0 / this%nzCells
        this%rdx = dble(this%nxCells)
        this%rdy = dble(this%nyCells)
        this%rdz = dble(this%nzCells)
        this%nxyCells = this%nxCells*this%nyCells
        this%nCells = this%nxyCells*this%nzCells

        if ( allocated(this%cellAtoms) ) deallocate(this%cellAtoms)
        allocate(this%cellAtoms(this%nCells))

        return
    end subroutine m_split

    subroutine m_fill(this, n, x)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(in) :: x

        ! handle memory stuff
        if ( this%nAtoms /= n ) then
            this%nAtoms = n
            if ( allocated(this%atomCell) ) deallocate(this%atomCell)
            allocate(this%atomCell(this%nAtoms))
            if ( allocated(this%atoms) ) deallocate(this%atoms)
            allocate(this%atoms(this%nAtoms))
        end if

        if ( .NOT. this%initialized ) return

        call m_fillcells(this, n, x, this%atomCell, this%cellAtoms, this%atoms)

    end subroutine m_fill

    subroutine m_fillPoints(this, n, x)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(in) :: x

        ! handle memory stuff
        if ( this%nPoints /= n ) then
            this%nPoints = n
            if ( allocated(this%pointCell) ) deallocate(this%pointCell)
            allocate(this%pointCell(this%nPoints))
            if ( allocated(this%points) ) deallocate(this%points)
            allocate(this%points(this%nPoints))
        end if

        if ( .NOT. this%initialized ) return

        ! check cellPoints allocation
        if ( allocated(this%cellPoints) ) then
            if ( size(this%cellPoints) /= this%nCells )  deallocate(this%cellPoints)
        endif
        if ( .not. allocated(this%cellPoints) ) allocate( this%cellPoints(this%nCells))

        call m_fillcells(this, n, x, this%pointCell, this%cellPoints, this%points)

    end subroutine m_fillPoints

    subroutine m_fillCells(this, n, xin, xCell, Cellx, xarray)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(in) :: xin
        integer, dimension(n) :: xCell
        type(CellAtomsData), dimension(n) :: Cellx
        integer, dimension(n) :: xarray
        real(8), dimension(3) :: r
        integer i, j, ii
        integer ix, iy, iz, icell
        integer, dimension(n) :: ordered
        integer, allocatable, dimension(:) :: orders
        integer, dimension(8) :: s, cs

        real(8), dimension(3,n) :: x

        xCell(1:n)=0
        xarray(1:n)=0
        forall(i=1:this%nCells)
            cellx(i)%n=0
            cellx(i)%start=0
        end forall

        x(:,:) = xin(:,:)
        call m_toFractional(this, n, x)

        ! create the ordered array. if placeordered is true the
        ! atom will be placed in 8 cells of equal volume defined
        ! by the box vetices. check placeordered documentation for
        ! more details
        if ( this%placeordered) then
            allocate( orders(n))
            s = 0
            do i = 1, n
                r(1:3) = x(1:3,i)
                j = 1
                if ( r(1) > 0.5d0) j = j + 1
                if ( r(2) > 0.5d0) j = j + 2
                if ( r(3) > 0.5d0) j = j + 4
                orders(i) = j
                s(j) = s(j) + 1
            end do
            if ( sum(s) /= n) &
                stop 'error in m_fillcells (simulation_box module @box_mod.F90)'
            cs(1) = 1
            do i = 2, 8
                cs(i) = cs(i-1) + s(i-1)
            end do
            do i = 1, n
                j = orders(i)
                ordered( cs(j)) = i
                cs(j) = cs(j) + 1
            end do
            deallocate( orders)
        else
            forall(i=1:n) ordered(i) = i
        end if

        do i = 1, n
            r(1:3) = x(1:3,i)
            ! wrap in the cell
            do j=1, 3
            r(j) = r(j) - int(r(j))
            if ( r(j) > 1.d0 ) then
                if ( r(j) - 1.d0 < EPS ) then
                r(j) = 1.d0
                else
                r(j) = r(j) - 1.d0
                end if
            else if ( r(j) < 0.d0 ) then
                if ( r(j) > -EPS ) then
                r(j) = 0.d0
                else
                r(j) = r(j) + 1.d0
                end if
            end if
            end do
            ix = int(r(1)*this%rdx)
            if ( ix == this%nxCells ) ix = ix - 1
            iy = int(r(2)*this%rdy)
            if ( iy == this%nyCells ) iy = iy - 1
            iz = int(r(3)*this%rdz)
            if ( iz == this%nzCells ) iz = iz - 1
            icell = 1 + ix + iy*this%nxCells + iz*this%nxyCells
            xCell(i)=icell
        end do

        ! find the number of x in the cells
        do i = 1, n
            icell = xCell(i)
            cellx(icell)%n = cellx(icell)%n + 1
        end do

        ! find the starting indexes in the x array for cells
        cellx(1)%start=1
        do i = 2, this%nCells
            cellx(i)%start = cellx(i-1)%start+cellx(i-1)%n
        end do

        ! zero the number of x in each cell
        cellx(1:this%nCells)%n=0

        ! place x in the correct order in xarray array
        do ii = 1, n
            i = ordered(ii)
            icell = xCell(i)
            xarray(cellx(icell)%start+cellx(icell)%n) = i
            cellx(icell)%n=cellx(icell)%n+1
        end do

        call m_toCartesian( this, n, x)

    end subroutine m_fillCells

    subroutine m_splitAndfill(this, nx, ny, nz, n, x)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(in) :: x

        this%initialized = .true.

        ! split the box
        call m_split(this,nx,ny,nz)

        ! allocate the arrays for atoms
        if ( this%nAtoms /= n ) then
            this%nAtoms = n
            if ( allocated(this%atomCell) ) deallocate(this%atomCell)
            allocate(this%atomCell(this%nAtoms))
            if ( allocated(this%atoms) ) deallocate(this%atoms)
            allocate(this%atoms(this%nAtoms))
        end if

        ! fill the box with
        call m_fill(this,n,x)

        return
    end subroutine m_splitAndfill

    function m_getAtomCell(this, iat) result(icell)
        class(SimulationBox), target, intent(in) :: this
        integer*4, intent(in) :: iat
        integer*4 :: icell
        icell = this%atomCell(iat)
    end function m_getAtomCell

    subroutine m_getCellAtoms(this, indx, n, a)
        class(SimulationBox), target, intent(in) :: this
        integer, intent(in) :: indx
        integer, intent(out)::  n
        integer, pointer, dimension(:) :: a

        integer start

        if ( indx == -1 .OR. indx > this%nCells ) then
            n=0
            a => NULL()
        else
            n=this%cellAtoms(indx)%n
            start = this%cellAtoms(indx)%start
            a=>this%atoms(start:start+n-1)
        end if

        return
    end subroutine m_getCellAtoms

    subroutine m_getCellPoints(this, indx, n, a)
        class(SimulationBox), target, intent(in) :: this
        integer, intent(in) :: indx
        integer, intent(out)::  n
        integer, pointer, dimension(:) :: a

        integer start

        if ( indx == -1 .OR. indx > this%nCells ) then
            n=0
            a => NULL()
        else
            n=this%cellAtoms(indx)%n
            start = this%cellAtoms(indx)%start
            a=>this%atoms(start:start+n-1)
        end if

        return
    end subroutine m_getCellPoints

    subroutine m_getNeighbourCells(this, indx, n, c)
        class(SimulationBox), target, intent(in) :: this
        integer, intent(in) :: indx
        integer, intent(out)::  n
        integer, dimension(27) :: c

        integer i, j, k
        integer icell, jcell, kcell
        integer ix, iy, iz
        integer tmp
        logical exclude

        ! basic arguments check
        if ( indx == -1 .OR. indx > this%nCells ) then
            n = 0
            forall(i=1:27) c(i)=-1
            return
        end if

        ! find cell indexes i, j, k from index
        tmp = indx-1
        k = tmp/(this%nxyCells)
        tmp = tmp - k * this%nxyCells
        k = k + 1
        j = tmp / this%nxCells
        tmp = tmp - j * this%nxCells
        j = j + 1
        i = tmp + 1

        ! loop over 27 near neighours
        n=0
        do iz = -1, 1
            do iy = -1, 1
                do ix = -1, 1
                    icell = i + ix
                    jcell = j + iy
                    kcell = k + iz
                    exclude=.false.
                    if ( icell == this%nxCells+1 ) then
                        icell=1
                        if ( .NOT. this%periodic ) exclude=.true.
                    else if ( icell == 0 ) then
                        icell=this%nxCells
                        if ( .NOT. this%periodic ) exclude=.true.
                    endif
                    if ( jcell == this%nyCells+1 ) then
                        jcell=1
                        if ( .NOT. this%periodic ) exclude=.true.
                    else if ( jcell == 0 ) then
                        jcell=this%nyCells
                        if ( .NOT. this%periodic ) exclude=.true.
                    endif
                    if ( kcell == this%nzCells+1 ) then
                        kcell=1
                        if ( .NOT. this%periodic ) exclude=.true.
                    else if ( kcell == 0 ) then
                        kcell=this%nzCells
                        if ( .NOT. this%periodic ) exclude=.true.
                    endif
                    if ( exclude ) then
                        tmp = -1
                    else
                        tmp=icell + (jcell-1) * this%nxCells + (kcell-1) * this%nxyCells
                    end if
                    n=n+1
                    c(n) = tmp
                end do
            end do
        end do

      return
    end subroutine m_getNeighbourCells

    subroutine m_toCartesian(this, n, x)
        class(SimulationBox), target, intent(inout) :: this
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(inout) :: x

        integer i

        ! set the transformation matrix
        call m_setMatrix(this)

        do i = 1, n
            x(1,i)=this%matrix(1)*x(1,i)+this%matrix(2)*x(2,i)+this%matrix(3)*x(3,i)+this%r0(1)
            x(2,i)=                      this%matrix(4)*x(2,i)+this%matrix(5)*x(3,i)+this%r0(2)
            x(3,i)=                                            this%matrix(6)*x(3,i)+this%r0(3)
        end do

        return
    end subroutine m_toCartesian

    subroutine m_toFractional(this, n, x)
        class(SimulationBox), intent(inout) :: this
        integer, intent(in) :: n
        real(8), dimension(3,n), intent(inout) :: x

        integer i

        ! set the inverce transformation matrix
        call m_setMatrix(this)
        this%matrix(1)=1.d0/this%matrix(1)
        this%matrix(4)=1.d0/this%matrix(4)
        this%matrix(6)=1.d0/this%matrix(6)

        do i = 1, n
            x(1:3,i) = x(1:3,i)-this%r0(1:3)
            x(3,i)=                              x(3,i)                *this%matrix(6)
            x(2,i)=       (x(2,i)               -x(3,i)*this%matrix(5))*this%matrix(4)
            x(1,i)=(x(1,i)-x(2,i)*this%matrix(2)-x(3,i)*this%matrix(3))*this%matrix(1)
        end do

        return
    end subroutine m_toFractional

    subroutine m_getProjections(this,x,y,z)
        class(SimulationBox), intent(in) :: this
        real(8), intent(out) :: x, y, z

        x = 1.d0/this%px
        y = 1.d0/this%py
        z = 1.d0/this%pz

        return
    end subroutine m_getProjections

    subroutine m_getVectors(this,r0_, a_, b_, c_, nXY_, nYZ_, nZX_)
        class(SimulationBox), intent(in) :: this
            
        real(8), dimension (3), intent(inout) :: r0_
        real(8), dimension (3), intent(inout) :: a_, b_, c_
        real(8), dimension (3), intent(inout) :: nXY_, nYZ_, nZX_

        r0_ = this%r0
        a_ = this%ux
        b_ = this%uy
        c_ = this%uz
        nXY_ = this%nXY
        nYZ_ = this%nYZ
        nZX_ = this%nZX

        return
    end subroutine m_getVectors

    subroutine m_destroy(this)
        class(SimulationBox), intent(inout) :: this
        integer i

        if ( allocated(this%atomCell) ) deallocate(this%atomCell)
        if ( allocated(this%cellAtoms) ) deallocate(this%cellAtoms)
        if ( allocated(this%atoms) ) deallocate(this%atoms)

        return
    end subroutine m_destroy

  end module simulation_box
