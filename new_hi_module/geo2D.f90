module geo2D_mod
    implicit none
    integer,parameter :: rk=8

     type :: point2D
        real(rk) :: x,y
    contains
        procedure :: init=>initP
        procedure :: getArray
        !procedure :: getX
        !procedure :: getY
    end type

    type,extends(point2D) :: point3D
        real(rk),private :: z
    !contains
        !procedure :: init=>initP
        !procedure :: getArray
    end type

    type :: edge2D
        type(point2D) :: head,tail
        real(rk) :: dir
        procedure(getX),pointer,nopass :: getFixed,getUnfixed
        procedure(setX),pointer,nopass :: setFixed,setUnfixed
    contains
        procedure :: init=>initE
        procedure :: associate_ptr
        procedure :: get_direct
    end type

    type :: elem2D
        type(point2D),dimension(:),allocatable :: nodes
        integer,allocatable :: edges(:,:)
    contains
        procedure :: get_std=>std_elem2D
    end type

    public :: setX,setY,getX,getY


contains

    !@ define a standard 2D elem 
    subroutine std_elem2D(this)
        class(elem2D) ::this
        !       4----7----3
        !       |         |
        !       8    9    6
        !       |         |
        !       1----5----2
        allocate(this%nodes(9))
        call this%nodes(1)%init(-1.0d0,-1.0d0)
        call this%nodes(2)%init(1.0d0 ,-1.0d0)
        call this%nodes(3)%init(1.0d0 , 1.0d0)
        call this%nodes(4)%init(-1.0d0, 1.0d0)
        call this%nodes(5)%init( 0.0d0,-1.0d0)
        call this%nodes(6)%init( 1.0d0, 0.0d0)
        call this%nodes(7)%init( 0.0d0, 1.0d0)
        call this%nodes(8)%init(-1.0d0, 0.0d0)
        call this%nodes(9)%init( 0.0d0, 0.0d0)

        allocate(this%edges(4,3))
        this%edges(1,:) = [1,2,5]
        this%edges(2,:) = [2,3,6]
        this%edges(3,:) = [3,4,7]
        this%edges(4,:) = [4,1,8]
    end subroutine


    ! ======================= edge related ==============
    ! @func : initial edge
    subroutine initE(this,p1,p2)
        class(edge2D) :: this
        type(point2D) :: p1,p2
        call this%head%init(getX(p1),getY(p1))
        call this%tail%init(getX(p2),getY(p2))
        call this%associate_ptr()
        call this%get_direct()

    end subroutine

    subroutine updatehead(this,p)
        class(edge2D) :: this
        type(point2D) :: p
        this%head%x = p%x
        this%head%y = p%y
    end subroutine

    subroutine updatetail(this,p)
        class(edge2D) :: this
        type(point2D) :: p
        this%tail%x = p%x
        this%tail%y = p%y
    end subroutine

    subroutine associate_ptr(this)
        implicit none

        class(edge2D) :: this

        if (dabs(this%head%x-this%tail%x)<1e-8) then
            this%getFixed => getX
            this%getUnfixed => getY

            this%setFixed => setX
            this%setUnfixed => setY
        else
            this%getFixed => getY
            this%getUnfixed => getX 

            this%setFixed => setY
            this%setUnfixed => setX 
        end if
    end subroutine

    subroutine get_direct(this)
        class(edge2D) :: this
        this%dir = dsign(1.d0,this%getUnfixed(this%tail)-this%getUnfixed(this%head))
    end subroutine


    ! ============== point2D related ===============

    subroutine initP(this,x,y,z)
        class(point2D) :: this
        real(rk) :: x,y
        real(rk),optional :: z
        this%x = x
        this%y = y
        select type (this)
        type is (point3D)
            if(present(z)) then 
                this%z=z
            else
                print *,"missing z value,exiting prog"
                stop
            end if
        end select
    end subroutine
    
    function getArray(this) result(res)
        class(point2D) :: this
        real(rk),dimension(:),allocatable :: res
        select type(this)
        type is (point2D)
            allocate(res(2))
            res=[this%x,this%y]
        type is (point3D)
            allocate(res(3))
            res=[this%x,this%y,this%z]
        end select
    end function
    ! ============  others related to point2D ==========
    function getX(p) result(res)
        type(point2D) :: p
        real(rk) :: res

        res = p%x
    end function

    function getY(p) result(res)
        type(point2D) :: p
        real(rk) :: res

        res = p%y
    end function

    subroutine setX(p,x) 
        type(point2D) :: p
        real(rk) :: x
        p%x = x
    end subroutine

    subroutine setY(p,y) 
        type(point2D) :: p
        real(rk) :: y
        p%y = y

    end subroutine

    function getDist(p1,p2) result(res)
        type(point2D) :: p1,p2
        real(rk) :: res
        res = norm2([p1%x-p2%x,p1%y-p2%y])
    end function


end module


    
