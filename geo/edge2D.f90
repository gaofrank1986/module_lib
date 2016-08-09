module edge2D_mod
    use pt2D_mod
    implicit none
    integer,parameter,private :: rk=8

    type :: edge2D
        type(point2D) :: head,tail
        integer,private :: fixed
        real(rk) :: dir
        procedure(getX2),pointer,nopass :: getFixed,getUnfixed
        procedure(setX2),pointer,nopass :: setFixed,setUnfixed
    contains
        procedure :: init=>initE
        procedure :: associate_ptr
        procedure :: get_direct
    end type
    private :: getX2,getY2,setX2,setY2
contains
    ! ======================= edge related ==============
    ! @func : initial edge
    subroutine initE(this,p1,p2)
        class(edge2D) :: this
        type(point2D) :: p1,p2
        if (getDist(p1,p2)<1e-8) then
            print *, "Error init edge, head and tail concident point"
            stop
        end if
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
            this%fixed=1
            this%getFixed => getX2
            this%getUnfixed => getY2

            this%setFixed => setX2
            this%setUnfixed => setY2
        else if  (dabs(this%head%y-this%tail%y)<1e-8) then
            this%fixed=2
            this%getFixed => getY2
            this%getUnfixed => getX2

            this%setFixed => setY2
            this%setUnfixed => setX2
        else
            this%fixed=3
            print *,"Warning,Edge not horizontal or veritcal"
        end if
    end subroutine

    subroutine get_direct(this)
        class(edge2D) :: this
        if(this%fixed.ne.3) then
            this%dir = dsign(1.d0,this%getUnfixed(this%tail)-this%getUnfixed(this%head))
        else
            this%dir=0
        end if
    end subroutine
    ! =================
    function getX2(p) result(res)
        class(point2D) :: p
        real(rk) :: res

        res = p%x
    end function

    function getY2(p) result(res)
        class(point2D) :: p
        real(rk) :: res

        res = p%y
    end function

    subroutine setX2(p,x) 
        class(point2D) :: p
        real(rk) :: x
        p%x = x
    end subroutine

    subroutine setY2(p,y) 
        class(point2D) :: p
        real(rk) :: y
        p%y = y
    end subroutine


end module
