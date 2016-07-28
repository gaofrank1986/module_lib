module mesh_mod
    implicit none
    integer,parameter :: rk=8

    type :: node
        private
        integer :: id
        real(rk),dimension(3) :: pos
    contains
        procedure,private :: initNode_v1
        procedure,private :: initNode_v2
        procedure :: pprint=>printNode
        generic :: init=>initNode_v1,initNode_v2
    end type
    
!    interface node
        !procedure:constructor=>const_node
    !end interface

    type,extends(node) :: nrml
        private
        integer :: base
    end type

    type :: elem
        private
        integer :: id,etype
        integer,dimension(8) :: nodlist
        integer,dimension(8) :: nrmlist
    contains
        procedure :: init=>initElem
        procedure :: pprint=>printElem
    end type

    type :: edge
        private
        integer :: id
        integer:: h,t
    contains
        !procedure,private :: swap
    end type

    type :: mesh
        private
        integer :: nnf,nelemf
        integer :: nnode,nelem
        type(node),dimension(:),allocatable :: nodes
        type(nrml),dimension(:),allocatable :: nrmls
        type(elem),dimension(:),allocatable :: elems
        ! nodearray,nrml array,elem array
    contains
        !procedure :: readmesh
    end type 

contains
    subroutine initNode_v1(this,n,x,y,z,b)
        class(node) :: this
        real(rk),intent(in) :: x,y,z
        integer,intent(in) :: n
        integer,intent(in),optional :: b
            this%pos=[x,y,z]
            this%id = n 
        select type(this)
        type is (nrml)
            if(present(b)) then
                this%base=b
            else
                print *,"Mesh Error: Nrml is not initialized with base"
            end if
        end select
    end subroutine

    subroutine initNode_v2(this,n,pos,b)
        class(node) :: this
        real(rk),intent(in),dimension(3) :: pos
        integer,intent(in) :: n
        integer,intent(in),optional :: b
        this%id = n
        this%pos=pos
        select type(this)
        type is (nrml)
            if(present(b)) then
                this%base=b
            else
                print *,"Mesh Error: Nrml is not initialized with base"
            end if
        end select
    end subroutine

    subroutine initElem(this,id,etype,ndlst,nmlst)
        class(elem) :: this
        integer,intent(in) :: etype,id
        integer,dimension(8) :: ndlst,nmlst
        this%id = id
        this%etype=etype
        this%nodlist = ndlst
        this%nrmlist = nmlst
    end subroutine

    subroutine printNode(this,u)
        class(node) :: this
        integer,optional :: u
        integer :: tmp
        if(present(u)) then
            tmp=u
        else
            tmp=6
        end if

        select type(this)
        type is (node)
            write (tmp, '("[mesh]  ",a5,i7,3f14.8)') 'Node',this%id,this%pos
        type is (nrml)
            write (tmp,'("[mesh]  ",a5,2i7,3f14.8)') 'Nrml',this%id,this%base,this%pos
        end select
    end subroutine

    subroutine printElem(this,u)
        class(elem) :: this
        integer,optional :: u
        integer :: tmp
        if(present(u)) then
            tmp=u
        else
            tmp=6
        end if
        
        write (tmp,'("[mesh]  ",a5,2i7)') 'Elem',this%id,this%etype
        write (tmp,'("[mesh]  ",a5,8i7)') ' ',this%nodlist
        write (tmp,'("[mesh]  ",a5,8i7)') ' ',this%nrmlist
    end subroutine
end module 




