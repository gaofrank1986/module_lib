program tes
    use io
    type(Ostream) :: test,test2
    !real(4) a
    integer :: a
    character(1024) :: tmp
    character(:),allocatable :: fd
    !call test%initialize("test",6)

    a=2
    !print *,test%title
    test=Ostream("test",[6,13])
    call test%fout(test%toString(a,'(i9)'))
    call test%set_header('new')
    call test%fout(test%toString(a,'(i9)'))
    print *,test%get_header()
    test2=Ostream("test",[6,13],['a.txt','b.txt'])

    !call test<<'2'
    !call test%fout(test%toString(a,'(f6.2)'),'2')
    print *,test%toString(a)
    
    !call tos('string')
    !character(len=:),allocatable :: string
    !string ="A"
    !print *,string
    !print *,getfilename("name",1)
    !tmp=getfilename("name",1)
    !print *,trim(tmp)
    open(101,file=getfilename("name",1))
    write(101,*) "hello world"
    close(101)
    fd=create_folder()
    print *,fd
    print *,timestamp()
    
contains
    subroutine tos(this)
        character(len=*) :: this
        print *,this
    end subroutine
end program
