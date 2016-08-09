program test_pt2d
    use pt2D_mod
    type(point3D) :: pt1,pt2
    type(point2D) :: pt3,pt4


    call pt3%putArray([1.0d0,2.0d0])
    !call pt3%setX1(2.0D0)
    !call pt4%putArray([4.0d0,6.0d0])
    call pt3%pprint()
    !call pt4%pprint()
    !print *,"dist=",getDist(pt3,pt4)

    !call pt1%putArray([1.0d0,2.0d0,3.0d0])
    !call pt2%putArray([4.0d0,6.0d0,8.0d0])
    !call pt1%pprint()
    !call pt2%pprint()
    !print *,"dist=",getDist(pt1,pt2)

    !call pt1%init(1.0d0,2.0d0,4.0d0)
    !!!print *,pt1
    !!call pt1%pprint()
    !call pt1%setX(2.0d0)
    !print *,pt1
    !!call pt1%setY(3.0d0)
    !!!print *,pt1
    !call pt1%setZ(3.0d0)
    !call pt1%pprint()
    !!!print *,pt1%getArray()
    !!call pt1%putArray([0.d0,0.d0])
    !!print *,pt1
    !call pt1%putArray([0.d0,0.d0,1.d0])
    !!print *,pt1
    !call pt1%pprint()
    !!call pt1%putArray([1.d0])
end program
