module param_mod
    use geo2D_mod
    integer,parameter,private :: rk=8

    interface 
        function f_bar(ndim,nf,cosn,drdx,drdn,shap)
            implicit none
            integer :: ndim,nf
            real(8),intent(in) :: drdx(ndim),drdn,cosn(ndim),shap(nf)
            real(8) :: f_bar(nf)       
        end function
    end interface

    type :: HSCoef
        real(8),allocatable,dimension(:) :: G,H
        real(8),allocatable,dimension(:,:) :: B
    contains
        procedure :: init=>init_coef
    end type

    type :: HSElem! QUAD
        integer :: nsub = 4!,nnode,ndim,nbdm
        real(rk) :: ck(3,8)
        real(rk) :: nk(3,8) 
        type(elem2D) :: mapped
        !type(trinfo),allocatable :: tripole
    contains
        procedure :: map2D
        procedure :: get_SF
        procedure :: get_nrml_at
        procedure :: get_nk
        !procedure :: divide_tri
        !procedure :: swap
        !procedure,private :: get_mapped_CDL
    end type

    type :: HSParams
        integer :: ndim=3
        integer :: nbdm=2
        integer :: node=8
        integer :: ngl=-10
        integer :: ngr=-10
        integer :: npw = 4
        integer  :: nf =8 
        integer :: npowg=4
        integer :: nnode = 8
        integer :: nelem = 1
        logical :: user_nrml=.true.
        real(8) :: xp(3),xip(2),xiq(2),beta
        type(HSCoef) :: mat
        real(8),dimension(10) :: gpl,gwl,gpr,gwr
        procedure(f_bar),pointer,nopass :: f_bar=>null()
    contains 
        procedure :: pprint    
        procedure :: init_mat
    end type

    interface swap_g2t
        module procedure :: swap_g2t_rank1
        module procedure :: swap_g2t_rank2
    end interface
    interface swap_t2g
        module procedure :: swap_t2g_rank1
        module procedure :: swap_t2g_rank2
    end interface

    private :: get_mapped_CDL,SHAPEF,DSHAPE
contains
!    subroutine div_tris(this,node)
        !class(HSElem) :: this
        !select case(node)
        !case(1)
            !this%tris(1) = [1,3,2]
            !this%tris(2) = [1,3,4]
        !case(2)
            !this%tris(1) = [2,4,1]
            !this%tris(2) = [2,4,3]
        !case(3)
            !this%tris(1) = [3,1,2]
            !this%tris(2) = [3,1,4]
        !case(4)
            !this%tris(1) = [4,2,1]
            !this%tris(2) = [4,2,3]
        !case(5)
            !this%tris(1) = [5,1,4]
            !this%tris(2) = [5,4,3]
            !this%tris(3) = [5,2,3]
        !case(6)
            !this%tris(1) = [6,1,2]
            !this%tris(2) = [6,1,4]
            !this%tris(3) = [6,4,3]
        !case(7)
            !this%tris(1) = [7,4,1]
            !this%tris(2) = [7,1,2]
            !this%tris(3) = [7,2,3]
        !case(8)
            !this%tris(1) = [8,4,3]
            !this%tris(2) = [8,1,2]
            !this%tris(3) = [8,2,3]
        !end select
    !end subroutine


    


    function get_mapped_CDL() result(ans)
        implicit none
        real(8),dimension(18) :: ans
        ans =[-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,-1.,0.,0.,0.]
    end function

    !===========================================

    !
    subroutine init_coef(this,npw,npowg,nf)
        class(HSCoef) :: this
        integer :: npw
        allocate(this%G(0:NPOWG),this%H(0:NPW))
        allocate(this%B(0:11,0:nf))
        this%B=0.0d0
        this%G=0.0d0
        this%H=0.0d0
    end subroutine

    ! @func : allocate coeff
    subroutine init_mat(this)
        class(HSParams) ::this
        call this%mat%init(this%npw,this%npowg,this%nf)
    end subroutine

    ! @func : print parameters
    subroutine pprint(this)
        class(HSParams) :: this
        print '(6a6)',"ndim","nbdm","node","nnode","nelem","nf"
        print '(6i6)',this%ndim,this%nbdm,this%node,this%nnode,this%nelem,this%nf
        print *,"npowg",this%npowg
        print *,"src2D",this%xip
        print *,"src3D",this%xp
        print *,"beta",this%beta
    end subroutine


    function map2D(this,x) result(ans)
        implicit none
        class(HSElem) ::this
        real(8),dimension(3) :: ans
        real(8) :: x(2),ri(3),sf(8),pt3D(3)
        pt3D = 0.0d0
        call shapef(3,8,x,this%ck,pt3D,ans,get_mapped_CDL(),sf)
    end function


    function get_SF(this,x) result(res)
        implicit none
        class(HSElem) ::this
        real(8),dimension(8) :: res
        real(8) :: x(2),pt3D(3),ri(3)
        pt3D = 0.0d0
        call shapef(3,8,x,this%ck,pt3D,ri,get_mapped_CDL(),res)
    end function

    !function get_DSF(this,x) result(res)
    !implicit none
    !class(HSElem) ::this
    !real(8),dimension(8) :: res
    !real(8) :: x(2),dir3D(3),ri(3),jacb,gd(3,8)
    !dir3D = 0.0d0
    !call dshapef(3,2,8,x,this%ck,dir3D,FJCB,get_mapped_CDL(),GD)

    !!SUBROUTINE DSHAPE(NDIM,NBDM,NODE,X,CK,COSN,FJCB,C,GD)
    !end function

    ! @func: get nrml at give 2D position[x]
    ! @var :[cosn]  normalized nrml vector
    ! @var :[fjcb]  jacobian determinant
    ! @var :[gd] two tangent vector
    subroutine get_nrml_at(this,x,cosn,fjcb,gd)
        implicit none
        class(HSElem) ::this
        real(8) :: cosn(3),fjcb,gd(3,2)
        real(8) :: x(2)
        call dshape(3,2,8,x,this%ck,cosn,FJCB,get_mapped_CDL(),GD)
    end subroutine

    subroutine get_nk(this)
        class(HSElem) :: this
        real(rk) :: tmp1(3,2),tmp,x(2)
        integer :: i

        if (allocated(this%mapped%nodes)) then

        else
            call this%mapped%get_std()
        end if

        do i=1,8
            x=this%mapped%nodes(i)%getArray()
            
            call this%get_nrml_at(x,this%nk(:,i),tmp,tmp1)
        end do
    end subroutine
    include './elem_ty_func.inc'

    !@ func : swap g2t
    function swap_g2t_rank1(input) result(ans)
        implicit none
        real(rk) :: input(:)
        real(rk),allocatable :: ans(:)
        integer :: tmp(2)
        if(rank(input).ne.1) then
            print *,"Error,not 1D matrix"
            stop
        end if
        if(size(input).ne.8) then
            print *,"Error! not 8 node elem"
            stop
        end if

        allocate(ans,source=input)
        ans=0.0d0

                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                !         7     6     5

                !         8           4

                !         1     2     3
        ans(1) = input(1)
        ans(2) = input(5)
        ans(3) = input(2)
        ans(4) = input(6)
        ans(5) = input(3)
        ans(6) = input(7)
        ans(7) = input(4)
        ans(8) = input(8)

    end function
    function swap_g2t_rank2(input) result(ans)
        implicit none
        real(rk) :: input(:,:)
        real(rk),allocatable :: ans(:,:)
        integer :: tmp(2)
        if(rank(input).ne.2) then
            print *,"Error,not 2D matrix"
            stop
        end if
        tmp = shape(input)
        if(tmp(2).ne.8) then
            print *,"Error! not 8 node elem"
            stop
        end if

        allocate(ans,source=input)
        ans=0.0d0

                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                !         7     6     5

                !         8           4

                !         1     2     3
        ans(:,1) = input(:,1)
        ans(:,2) = input(:,5)
        ans(:,3) = input(:,2)
        ans(:,4) = input(:,6)
        ans(:,5) = input(:,3)
        ans(:,6) = input(:,7)
        ans(:,7) = input(:,4)
        ans(:,8) = input(:,8)
    end function
    function swap_t2g_rank1(input) result(ans)
        implicit none
        real(rk) :: input(:)
        real(rk),allocatable :: ans(:)
        if(rank(input).ne.1) then
            print *,"Error,not 1D matrix"
            stop
        end if
        if(size(input).ne.8) then
            print *,"Error! not 8 node elem"
            stop
        end if       
        allocate(ans,source=input)
        ans=0.0d0

        ans(1) = input(1)
        ans(2) = input(3)
        ans(3) = input(5)
        ans(4) = input(7)
        ans(5) = input(2)
        ans(6) = input(4)
        ans(7) = input(6)
        ans(8) = input(8)
    end function

    function swap_t2g_rank2(input) result(ans)
        implicit none
        real(rk) :: input(:,:)
        real(rk),allocatable :: ans(:,:)
        integer :: tmp(2)
        if(rank(input).ne.2) then
            print *,"Error,not 2D matrix"
            stop
        end if
        tmp = shape(input)
        if(tmp(2).ne.8) then
            print *,"Error! not 8 node elem"
            stop
        end if       
        allocate(ans,source=input)
        ans=0.0d0

        ans(:,1) = input(:,1)
        ans(:,2) = input(:,3)
        ans(:,3) = input(:,5)
        ans(:,4) = input(:,7)
        ans(:,5) = input(:,2)
        ans(:,6) = input(:,4)
        ans(:,7) = input(:,6)
        ans(:,8) = input(:,8)
    end function

    function indx(inp) result(ans)
        implicit none
        integer,intent(in) :: inp
        integer            :: ans

        select case(inp)
        case (1)
            ans=1
        case(2)
            ans=5
        case(3)
            ans=2
        case(4) 
            ans=6
        case(5)
            ans=3
        case(6)
            ans=7
        case(7)
            ans=4
        case(8)
            ans=8
        end select
    end function

end module



