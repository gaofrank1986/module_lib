program test
    use mesh_mod
    type(node) :: n1
    type(nrml) :: n2,n3
    type(elem) :: e1,e2
    real(8) :: a,b,c
    integer d(8)
    d=1
    a=1.0
    b=2.0
    c=3.0

    call n1%init(1,a,b,c)
    call n1%pprint()
    call n1%init(2,1.0d0,1.0d0,1.0d0)
    call n1%pprint()
    call n2%init(1,a,b,c)
    call n2%pprint()
    call n3%init(1,a,b,c,2)
    call n3%pprint()
    call e1%init(1,8,d,d)
    call e1%pprint()

end program
