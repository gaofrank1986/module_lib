program test
    use geo2D_mod
    !procedure(getX),pointer :: p1,p2,p3
    type(point2D) :: p1,p2
    type(point3D) :: p3
    type(edge2D) :: e
    type(elem2D) :: ce

   call p1%init(2.0d0,3.0d0) 
   call p2%init(5.0d0,3.0d0) 
   call e%init(p1,p2)

    !call associate_ptr(p1,p2,p3)
    print *,e%getUnfixed(e%head)
    call e%setUnfixed(e%head,7.0d0)
    print *,e%getUnfixed(e%head)
    call setX(p1,4.5d0)
    print *,p1%getArray()
    !print *,p1%x
    call p3%init(1.0d0,2.0d0,3.0d0)
    print *,p3%getArray()
    !call p3%init(1.0d0,2.0d0)

    call ce%get_const()
    print *, " cnst elem"
    print *,ce%edges(1,:)
    print *,ce%nodes(1)
    print *,ce%nodes(2)

    !print *,p1(p)
    !print *,p%p1
    !print *,p1(p)
    !print *,getX(p)

end program
