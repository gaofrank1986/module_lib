
          program main
            use iso_c_binding
            implicit none
              type :: mem
                  real(8),allocatable :: b(:,:,:)
              end type

            interface
              subroutine my_routine() bind(c,name='pmem')
              end subroutine
            end interface
            
            real(8) :: a(8,9)
            integer :: i
            type(mem) :: c(4)
            a=10
            do i = 1,4
            allocate(c(i)%b(1024,1024,20))
            c(i)%b = 1.0
            call my_routine()
        end do
            pause
        do i = 1,4
            deallocate(c(i)%b)
        end do
            call my_routine()
        
            pause
            allocate(c(i)%b(1024,1024,20))
            c(i)%b = 1.0
            call my_routine()
            pause


          end program main
