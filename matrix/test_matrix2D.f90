program test
    use matrix_mod

    type(matrix2D) :: mat
    real(8) :: a(0:1,2)
    real(8),allocatable :: b(:,:),c(:)
    integer :: tmp(2)
    a = reshape([1,2,3,4],[2,2])
    !print *,a
    !print *,rank(shape(a))
    !print *,shape(a)
    !print *,rank(a)
    tmp=shape(a)
    print *,lbound(a,1),ubound(a,1)
    allocate(c(tmp(1)))
    print *,len('hello')
    call mat%init(a)
    !call mat%assign_row_name(['x','y'])
    call mat%pprint()
end program
