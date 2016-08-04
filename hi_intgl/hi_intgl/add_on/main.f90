program main
    use hi_intg
    !use matrix_mod

    !type(Matrix2D) :: m

    real(8) res(8),r2(8),r3(8)
    call read_model_from_DAT()
    
!    call m%init(full_mesh_matrix(:,:,1))
    !call m%pprint()
    call eval_singular_elem(full_mesh_matrix(:,:,1),res,r2,r3)

    print *,res
    print *,sum(res)
end program
