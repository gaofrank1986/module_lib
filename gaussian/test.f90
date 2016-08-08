program test
    use gaussm3
    implicit none
    integer :: ngp 
    real(8) :: a(10),b(10)

    ngp=8
    call gauleg(ngp,a,b)
    print *,a
    print *,b
    print *, qgauss(f,-2.d0,2.d0,8)
    contains
      real(8)  function f(x)
            real(8) :: x
            f=x**2
            end function
    end program
