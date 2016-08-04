module mesh
    integer :: nnode,nelem,nnoded
end module

program test
    use mesh
    real(8) :: m1(3,8)
    nnode = 8
    nelem = 1
    nnoded = 8
    

    call init_hi_var()

    
