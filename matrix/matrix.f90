module Matrix_mod
    integer,parameter :: rk = 8
    type :: Matrix2D
        private
        real(rk),dimension(:,:),allocatable :: m
        character(len=10),dimension(:),allocatable :: row_name,col_name
    contains
        procedure :: init
        procedure :: pprint
        procedure :: assign_row_name
    end type

contains
    subroutine init(this,i)
        class(Matrix2D) :: this
        real(rk),dimension(:,:) ::i
        allocate(this%m,source=i)
    end subroutine

    subroutine assign_row_name(this,rn)
        class(Matrix2D) :: this
        character(len=*),dimension(:) :: rn 
        integer :: tmp(2)
        if(not(allocated(this%m))) then
            print *,"Error,matrix not initialized"
            stop
        end if
        tmp = int(shape(this%m))
        if(not(tmp(1).eq.size(rn))) then
            print *,"Error,row_name input size wrong"
            stop
        end if
        allocate(this%row_name(tmp(1)))
        this%row_name=rn
    end subroutine

!    subroutine assign_col_name(this,cn)
        !class(Matrix2D) :: this
        !real(rk) :: tmp(2)
        !if(not(allocated(this%m))) then
            !print *,"Error,matrix not initialized"
            !stop
        !end if
        !tmp = shape(this%m) 
        !if(not(tmp(2).eq.size(cn))) then
            !print *,"Error,col_name input size wrong"
            !stop
        !end if
        !allocate(this%col_name(tmp(2)))
        !this%col_name=cn
    !end subroutine

    subroutine pprint(this)
        class(Matrix2D) :: this
        integer :: i,j,tmp(2)
        character(len=100) :: fm,ts
        tmp = shape(this%m)

        print *,"=========== Printing Matrix ==========="
        write(6,*) ' '

        do i=1,tmp(1)
            if(i.eq.1) then
                write (6,fmt='(a15)',advance='no') ' '
                do j=1,tmp(2)
                    if(not(allocated(this%col_name))) then
                        write (6,fmt='(a4,i4,a1,5x)',advance='no') '[Col',j,']'
                    else
                        write (6,fmt='(a1,a7,a1,5x)',advance='no') '[',this%col_name(j),']'
                    end if
                end do
                write(6,*) ' '          
            end if


            if(not(allocated(this%row_name))) then
                write (6,fmt='(a4,i4,a1,1x)',advance='no') '[Row',i,']'
            else
                write (6,fmt='(a1,a7,a1,1x)',advance='no') '[',this%row_name(i),']'
            end if

            do j=1,tmp(2)
                write (6,fmt='(f14.8)',advance='no') this%m(i,j)
            end do
            write(6,*) ' '          
        end do
        print *,"=================================="
    end subroutine
end module
            




