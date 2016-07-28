module body_mod
    !use kinds
    implicit none
    public
    integer,parameter :: rk=8

    type :: rigid_body
        private
        real(rk) :: mass
        ! @var: [displ,velc,accel] disp,velc,acceleration at six direction
        real(rk),dimension(6) :: displ,velc,accel
        ! @var: roation,mass center
        real(rk),dimension(3) :: rot_ctr,mass_ctr
        ! @var: [m_mat] mass matrix
        ! @var: [ks_mat] stiffness matrix
        ! @var: [cdmp_m] damping- matrix
        real(rk),dimension(6,6) :: m_mat,ks_mat,b_mat,cdmp_m
    contains
        procedure,private :: initialize
        !procedure :: calc_mass_matrix
    end type



    type, extends(rigid_body) :: floatBD
        private
        real(rk),dimension(3) :: buoy_ctr
        ! @var: [a_mat] added mass matrix
        ! @var: [kh_mat] hydro-restoring matrix,work like spring, put in k matrix
        real(rk),dimension(6,6) :: a_mat,kh_mat
        
        real(rk) :: volm,area,xf,yf,xk2,yk2,xcf

    contains
        procedure ::  read_body_info
        !procedure,private:: toString_int
        !procedure,private :: toString_real
        !generic :: toString =>toString_int,toString_real

    end type

    interface floatBD 
        procedure :: constructor
        !procedure :: destructor
    end interface

contains

    ! @func : same as C++ constructor
    function constructor()
        type(floatBD) ::constructor
        call constructor%initialize()
    end function


    ! @func : used in constructor to intialized object 
    subroutine initialize(this)
        class (rigid_body) ::  this

        this%displ=0.0d0
        this%velc=0.0d0
        this%accel=0.0d0

        select type(this)

        type is (rigid_body)
            this%m_mat=0.0d0!mass matrix
            this%ks_mat=0.0d0!stiffness
            this%b_mat=0.0d0!damping

        type is (floatBD)
            this%a_mat=0.0d0!added mass
            this%kh_mat=0.0d0!hydro-restoring

        class default
            ! give error for unexpected/unsupported type
            stop 'initialize: unexpected type for sh object!'
        end select

    end subroutine initialize
    

    subroutine read_body_info(this,filename)
        class(floatBD) :: this
        character(len=*),intent(in) :: filename

        character(len=1024) :: tmp
        integer i,j,ifl_mass
        real(rk) :: xia(3,3)

        open(unit=7,file=filename,status='unknown')

        read(7,*)   tmp
        read(7,*)   ifl_mass 

        read(7,*)   tmp
        read(7,*)   (this%rot_ctr(i),i=1,3) 

        read(7,*)   tmp
        read(7,*)   (this%mass_ctr(i),i=1,3) 

        read(7,*)   tmp
        !read(7,*)    this%area,xf,yf,xk2,yk2,xcf

        read(7,*)   tmp
        read(7,*)   this%volm,(this%buoy_ctr(i),i=1,3) 

        if(ifl_mass .eq. 0) then
            read(7,*)   tmp
            read(7,*)   this%mass
            do i=1, 3
                read(7,*)  (xia(i,j),  j=1,3)
            end do
            this%m_mat = calc_mass_matrix(this%mass,this%mass_ctr,this%rot_ctr,xia)

        else if(ifl_mass .eq. 1) then
            read(7,*)   tmp
            do i=1, 6
                read(7,*)  (this%m_mat(i,j),   j=1,6)
            end do
        endif

        ! read hydro-restoring matrix
        read(7,*)   tmp
        do i=1, 6
            read(7,*)  (this%kh_mat(i,j),   j=1,6)
        end do

        ! read stiffness matrix
        read(7,*)   tmp
        do i=1, 6
            read(7,*)  (this%ks_mat(i,j),   j=1,6)
        end do

        ! read damping matrix
        read(7,*)   tmp
        do i=1, 6
            read(7,*)  (this%b_mat(i,j),  j=1,6)
        end do

    end subroutine
    

    !subroutine pprint_matrix2D(mat)

    !end subroutine


    ! @func:
    ! @param: [xmas] body mass
    ! @param: [mc] mass center
    ! @param: [c] rotation center
    ! @output: [m] mass matrix
    function calc_mass_matrix(xms,mc,c,xia) result(m)
        implicit none

        real(rk),dimension(6,6) :: m
        real(rk),intent(in) :: xms,mc(3),c(3),xia(3,3)
        real(rk) :: tmp(3,3),tmp2(3,3)
        integer :: i

        m = 0.0d0

        do i=1,3
            m(i,i) = xms
        end do

        tmp(1,1) = xia(2,2)+xia(3,3)+xms*(mc(2)-c(2))**2+xms*(mc(3)-c(3))**2
        tmp(2,2) = xia(1,1)+xia(3,3)+xms*(mc(1)-c(1))**2+xms*(mc(3)-c(3))**2
        tmp(3,3) = xia(1,1)+xia(2,2)+xms*(mc(2)-c(2))**2+xms*(mc(3)-c(3))**2

        tmp(2,1) = -xia(2,1)-xms*(mc(2)-c(2))*(mc(1)-c(1))
        tmp(1,2) = -xia(1,2)-xms*(mc(1)-c(1))*(mc(2)-c(2))

        tmp(1,3) = -xia(1,3)-xms*(mc(1)-c(1))*(mc(3)-c(3))
        tmp(3,1) = -xia(3,1)-xms*(mc(1)-c(1))*(mc(3)-c(3))

        tmp(2,3) = -xia(2,3)-xms*(mc(2)-c(2))*(mc(3)-c(3))
        tmp(3,2) = -xia(3,2)-xms*(mc(2)-c(2))*(mc(3)-c(3))

        tmp2(1,3) = -xms*(mc(2)-c(2))
        tmp2(1,2) = xms*(mc(3)-c(3))
        tmp2(2,3) = xms*(mc(1)-c(1))

        tmp2=tmp2-transpose(tmp2)
        ! this is negative symmetric

        m(4:6,4:6)=tmp
        m(1:3,4:6)=tmp2
        m(4:6,1:3)=transpose(tmp2)
        ! m right-upper 3x3 matrix is symmetric to left-lower 3x3 matrix

    end function



end module body_mod 

