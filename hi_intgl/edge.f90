module edge_mod
    implicit none
    integer,parameter :: rk=8

     type :: point2D
        real(rk) :: x,y
    contains
        procedure :: init=>initP
        !procedure :: getX
        !procedure :: getY
    end type

    type :: edge2D
        type(point2D) :: head,tail
        real(rk) :: dir
        procedure(getX),pointer,nopass :: getFixed,getUnfixed
        procedure(setX),pointer,nopass :: setFixed,setUnfixed
    contains
        procedure :: init=>initE
        procedure :: associate_ptr
    end type

    private :: setX,setY,getX,getY

             
contains
    subroutine initE(this,p1,p2)
        class(edge2D) :: this
        type(point2D) :: p1,p2
        call this%head%init(getX(p1),getY(p1))
        call this%tail%init(getX(p2),getY(p2))
        call this%associate_ptr()
        call this%get_direct()

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
        class(edge2D),intent(in) :: this
        this%dir = dsign(1.d0,this%getUnfixed(this%tail)-this%getUnfixed(this%head))
    end subroutine



    
    function get_step(side,slice) result(res)

        type(point2D) :: src2D
        
        real(rk) :: tmp
        ! distance from slice's tail to side edge's head
        d0 = this%dir*(this%getUnFixed(slice%tail)-this%getUnfixed(side%head))
        if(d0<1.0e-8) then
            !print *, "if seg_start located out of edge,exit num_converge loop"
            res = -9.0d0
            return
        end if

        !dist_unfixed = edge_direct*(src_lcl(unfixed)-seg_start(unfixed))
        !d1 =   dist_fixed2=(seg_start(fixed)-src_lcl(fixed))**2
        d1 =           side%getFixed(src2D)  - side%getFixed(slice%head)
        d2 = side%dir*(side%getUnfixed(src2D)- side%getUnfixed(slice%head))

        ! compute step size====================

        if (d2 < 1.d-8) then !if src2D outside slice
            !seg_step = FK*norm2(seg_start-src_lcl)  ! fk is a factor?
            seg_step = fk*getDist(src2D,slice%head)
        else 
            tmp=d1**2*+d2**2 - (d1*fk)**2
            if(tmp < 0.d0) then
                seg_step=d2
            else
                seg_step=fk*(fk*d2-dsqrt(tmp))/(fk**2-1.)
            endif
        endif

        bool_expr = side%dir*(side%dir*seg_step+side%getUnfixed(slice%head))+1.d-8 

        if (bool_expr > d0) then
            ! if next step goes beyond the edge, use the end node
            seg_step = d0
        endif
        half_step = 0.5D0*side%dir*seg_step
        res = half_step
    end function


!   subroutine get_direct(this)
        !class(edge2D),intent(in) :: this
           !this%dir = dsign(1.d0,this%getUnfixed(this%tail)-this%getUnfixed(this%head))
   !end subroutine


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

      subroutine initP(this,x,y)
        class(point2D) :: this
        real(rk) :: x,y
        this%x = x
        this%y = y
    end subroutine

end module


    
