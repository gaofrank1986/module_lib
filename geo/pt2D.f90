module pt2D_mod
    implicit none
    integer,parameter,private :: rk=8
    !public :: point2D,point3D

    type :: point2D
        real(rk) :: x,y
    contains
        procedure :: init=>initP
        procedure :: getArray
        procedure :: putArray
        procedure :: getX
        procedure :: getY
        procedure :: setX
        procedure :: setY
        procedure :: pprint
    end type

    type,extends(point2D) :: point3D
        real(rk),private :: z
    contains
        procedure :: getZ
        procedure :: setZ
    end type 

    interface getDist
        module procedure getDist2D
        module procedure getDist3D
    end interface

contains
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
    
    subroutine putArray(this,a)
        class(point2D) :: this
        real(rk) :: a(:)
        if (size(a)>3) then
            print *,"Error PutArray,size larger than 3"
            stop
        else
            select type (this)
            type is (point3D)
                if(size(a).eq.3) then
                    this%z=a(3)
                    this%y=a(2)
                    this%x=a(1)
                else
                    print *,"Error PutArray,size not compatble"
                    stop
                end if
                type is (point2D) 
                    if(size(a).eq.2) then
                        this%y=a(2)
                        this%x=a(1)
                    else
                        print *,"Error PutArray,size not compatble"
                        stop
                    end if
                end select
            end if
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

    subroutine pprint(this,s)
        class(point2D) :: this
        character(len=1024),optional :: s
        if (present(s)) then
            write(s,'(3f14.8)') this%getArray()
        else 
            print '(3f14.8)',this%getArray()
        end if
    end subroutine

    ! ============  others related to point2D ==========
    function getX(this) result(res)
        class(point2D) :: this
        real(rk) :: res

        res = this%x
    end function

    function getY(this) result(res)
        class(point2D) :: this
        real(rk) :: res

        res = this%y
    end function

    function getZ(this) result(res)
        class(point3D) :: this
        real(rk) :: res
        select type(this)
        type is(point3d)
            res = this%y
        class default
            print *,"Error getZ,no Z component"
        end select
    end function

    subroutine setX(this,x) 
        class(point2D) :: this
        real(rk) :: x
        this%x = x
    end subroutine

    subroutine setY(this,y) 
        class(point2D) :: this
        real(rk) :: y
        this%y = y
    end subroutine

    subroutine setZ(this,z) 
        class(point3D) :: this
        real(rk) :: z
        select type(this)
        type is(point3D)
            this%z = z
        class default
            print *,"Error setZ,no Z componet"
        end select

    end subroutine
    !===================================================

    function getDist2D(p1,p2) result(res)
        type(point2D) :: p1,p2
        real(rk) :: res
        res = norm2([p1%x-p2%x,p1%y-p2%y])
    end function

    function getDist3D(p1,p2) result(res)
        type(point3D) :: p1,p2
        real(rk) :: res
        res = norm2([p1%x-p2%x,p1%y-p2%y,p1%z-p2%z])
    end function


end module


    

