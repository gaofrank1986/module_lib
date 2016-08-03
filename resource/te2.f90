
          program main
            use iso_c_binding
            implicit none
              type :: mem
                  real(8),allocatable :: b(:,:,:)
              end type

            interface
              subroutine my_routine(p) bind(c,name='getCurrentRSS')
                  import :: c_size_t
                  integer(c_size_t) :: p
              end subroutine
            end interface
            
            real(8) :: a(8,9)
            integer(c_size_t) :: i,res
            type(mem) :: c(4)
            a=10
            do i = 1,4
            allocate(c(i)%b(1024,1024,20))
            c(i)%b = 1.0
            call my_routine(res)
            print *,res/1024/1024,"MB"
        end do
            pause
        do i = 1,4
            deallocate(c(i)%b)
        end do
            call my_routine(res)
            print *,res/1024/1024,"MB"
        
            pause
            allocate(c(i)%b(1024,1024,20))
            c(i)%b = 1.0
            call my_routine(res)
            print *,res/1024/1024,"MB"
            pause


          end program main
