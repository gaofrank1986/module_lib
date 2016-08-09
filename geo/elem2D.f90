module elem2D_mod
    use pt2D_mod
    use edge2D_mod
    implicit none
    integer,parameter,private :: rk=8

    type :: elem2D
        type(point2D),dimension(:),allocatable :: nodes
        integer,allocatable :: edges(:,:)
    contains
        procedure :: get_std=>std_elem2D
    end type

    !public :: setX,setY,getX,getY
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

end module
