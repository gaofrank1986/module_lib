program test_edge
    use edge2D_mod
    type(point2D) :: pt1,pt2
    type(edge2D) :: s1

    call pt1%putArray([1.0d0,3.d0])
    call pt2%putArray([3.0d0,3.d0])
    call pt1%pprint()
    call pt2%pprint()
    call s1%init(pt1,pt2)
    print *,s1%getFixed(s1%head)
    print *,s1%getUnFixed(s1%head)
end program
