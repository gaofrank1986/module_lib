program main
    use hi_intg
    use matrix_mod

    type(Matrix2D) :: m

    real(8) res(8),glb(3),glb_ctr(3),r1(8),r2(8)
    call read_model_from_DAT()
    
    call m%init(full_mesh_matrix(:,:,1))
    call m%pprint()
    glb=full_mesh_matrix(:,6,1)
    glb(2) = 0.66d0

    !call preset_src(0.66d0,0.0d0,[0.0d0,0.5d0,0.0d0],[0.0d0,0.0d0,0.0d0])
    call preset_src(0.66d0,0.0d0,glb,glb_ctr)
    call eval_singular_elem(full_mesh_matrix(:,:,1),res,r1,r2)

    print *,res
    print *,sum(res)
end program
