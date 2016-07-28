module io
    implicit none
    public

    type Ostream
        private
        ! @var : [port] the port number list
        integer,allocatable :: port(:)
        ! @var: [header] prompt at beginning
        character(len=:),allocatable :: header
        ! @var: [namelst] input file names,optional
        character(len=20),dimension(:),allocatable :: namelst
        ! number list
        ! todo name list or a dictionary
        ! @var : [reg_name] file name register to track all files declared
        character(len=30),dimension(100) :: reg_name
    contains
        procedure,private :: initialize
        procedure,private :: dump
        procedure ::  fout
        procedure ::  set_header
        procedure ::  get_header
        procedure,private:: toString_int
        !procedure,private :: toString_real8
        !procedure,private :: toString_real4
        procedure,private :: toString_real
        generic :: toString =>toString_int,toString_real
        !generic :: toString =>toString_int,toString_real8,toString_real4,toString_realk
        !generic :: operator(<) => fout

    end type
    ! @ interface : this interface bound constructor with type Ostream to mimic a constructor behavior
    interface Ostream 
        procedure :: constructor
        procedure :: destructor
    end interface

contains
    ! @func : same as C++ constructor
    function constructor(string,port,namelst)
        type(Ostream) ::constructor
        character(len=*),intent(in)::string
        integer,intent(in) :: port(:)
        character(len=*),dimension(:),intent(in),optional :: namelst

        call constructor%initialize(string,port,namelst)
    end function

    function destructor()
        type(Ostream) ::destructor
        call destructor%dump()
    end function

    ! @func : getter
    function get_header(this)
        character(:),allocatable :: get_header,tmp
        class (Ostream) ::  this
        !write 
        tmp = trim(this%header)
        print *,len(tmp)
        get_header = tmp(2:len(tmp)-1) 

    end function 
    
    ! @func : setter
    subroutine set_header(this,string)
        class (Ostream) ::  this
        character(len=*),intent(in)::string
        this%header='['//string//']'
    end subroutine 

    ! @func : used in constructor to intialized object 
    ! @param: [string] header appear at beginning
    ! @param: [port] which termianl to write output, 6 for stdout, other for file
    subroutine initialize(this,string,port,namelst)
        class (Ostream) ::  this
        character(len=*),intent(in)::string
        integer,intent(in) :: port(:)
        character(len=*),dimension(:),intent(in),optional :: namelst
        
        integer :: tmp,i

        this%header='['//string//']'
        tmp = size(port)
        allocate(this%port(tmp))
        this%port=port
        if(present(namelst)) then
            if(size(namelst).ne.size(port)) then
                print *,"IO Error, name list and port list not same length."
                stop
            else
                allocate(this%namelst(tmp))
                this%namelst = namelst
                do i=1,tmp 
                    print '("[io]  ",a,i3,3x,a5,a)','terminal',this%port(i)&
                        &,'=>   ',this%namelst(i)
                    !print *,'terminal',this%port(i),trim(this%namelst(i))
                    !print *,trim(this%namelst(i))
                    if((this%port(i).eq.6).and.(this%namelst(i).ne.'')) then
                            print '("[io]  ",a)',"IO Error, port 6 &
                                &is stdout cannot be assigned with file name."
                    end if
                end do
            end if
        end if
        this%reg_name=''
    end subroutine initialize

    subroutine dump(this)
        class (Ostream) ::  this

    end subroutine
    ! @func : output given string to designated port
    subroutine fout(this,string)
        class (Ostream),intent(in) ::  this
        character(len=*),intent(in)::string
        !character(len=*),optional::ctl
        integer :: i 
        if(allocated(this%namelst)) then
            !do i=1,size(this%port) 
                !write(this%port(i),fmt='(a)') this%header//' '//string
            !end do
        else
            do i=1,size(this%port) 
                write(this%port(i),fmt='(a)') this%header//' '//string
            end do

        end if
    end subroutine
   
    ! @func : [toString_int] convert a integer into string 
    ! @ param: ctl optional, format of integer output
    ! todo this need to be polymorphism for integer/logical/double/complex etc.
    function toString_int(this,arg,ctl) result(string)
        class (Ostream),intent(in) ::  this
        character(*),optional :: ctl
        integer,intent(in)::arg
        character(:),allocatable :: string
        character(len=1024) :: tmp
        if(present(ctl)) then
            write(tmp,fmt=ctl) arg
        else
            write(tmp,fmt='(i7)') arg ! output to string
        end if
        string=trim(tmp)
    end function

    !function toString_real8(this,arg,ctl) result(string)
        !class (Ostream) ::  this
        !character(:),allocatable :: string
        !character(*),optional :: ctl
        !character(len=1024) :: tmp
        !real(8)::arg
        
        !if(present(ctl)) then
            !write(tmp,fmt=ctl) arg
        !else
            !write(tmp,fmt='(f10.5)') arg 
        !end if
        !string=trim(tmp)
    !end function

    !function toString_real4(this,arg,ctl) result(string)
        !class (Ostream) ::  this
        !character(:),allocatable :: string
        !character(*),optional :: ctl
        !character(len=1024) :: tmp
        !real(4)::arg
        
        !if(present(ctl)) then
            !write(tmp,fmt=ctl) arg
        !else
            !write(tmp,fmt='(f8.4)') arg 
        !end if
        !string=trim(tmp)
    !end function

    ! @func : [toString_real] convert a given real number to a string
    ! @param: [ctl] optional, control the format
    function toString_real(this,arg,ctl) result(string)
        class (Ostream) ::  this
        character(:),allocatable :: string
        character(*),optional :: ctl
        character(len=1024) :: tmp
        real::arg
        
        if(present(ctl)) then
            write(tmp,fmt=ctl) arg
        else
            write(tmp,fmt='(f10.5)') arg 
        end if
        string=trim(tmp)
    end function
end module io

