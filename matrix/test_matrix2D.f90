program test
    use matrix_mod

    type(matrix2D) :: mat
    real(8) :: a(2,2)
    real(8),allocatable :: b(:,:),c(:)
    integer :: tmp(2)
    a = reshape([1,2,3,4],[2,2])
    !print *,a
    !print *,rank(shape(a))
    !print *,shape(a)
    !print *,rank(a)
    tmp=shape(a)
    allocate(c(tmp(1)))
    call mat%init(a)
    call mat%assign_row_name(['x','y'])
    call mat%pprint()
end program
