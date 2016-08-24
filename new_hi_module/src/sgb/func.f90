module util_funcs
    integer,parameter,private:: rk=8
contains

    !!!
    !! @funct : given a angle, return the positon of intersection at unit square
    !! @param : [theta] input angle range in [0,2*pi]
    !! @param : [tol] optional, tolerance for real comparation
    !!!
    function intersect_unit_square(theta,tol) result(x)
        implicit none
        integer,parameter :: rk=8
        real(rk) theta
        real(rk) :: eps,x(2)
        real(rk) :: pi = 4*atan(1.0d0)
        real(rk),optional :: tol 
        x=0.0d0

        !=============
        if(present(tol)) then
            eps=tol
        else
            eps=10e-6
        end if

        if((theta<0.0d0).or.((theta-2*pi)>10e-8)) then
            print *,"Error! theta fall out defined interval"
        end if

        if(dabs(theta-pi*1/4)<eps) x=[1,1]
        if(dabs(theta-pi*3/4)<eps) x=[-1,1]
        if(dabs(theta-pi*5/4)<eps) x=[-1,-1]
        if(dabs(theta-pi*7/4)<eps) x=[1,-1]
        if((theta-pi*0/4)<eps) x=[0,1]
        if((pi*8/4-theta)<eps) x=[0,1]

        ! 1/4 pi <= theta <= 3/4 pi

        if(((theta-pi*1/4).ge.eps).and.((pi*3/4-theta).ge.eps)) then
            !x(1) = dsign(1.0d0,dcos(theta))/dabs(dtan(theta))
            x(1) = 1/(dtan(theta))
            x(2) = 1

        else if(((theta-pi*3/4).ge.eps).and.((pi*5/4-theta).ge.eps)) then
            x(1) = -1
            !x(2) = dsin(theta)
            !x(2) = -dsign(1.0d0,dsin(theta))*dtan(theta)
            x(2) = -1*dtan(theta)

            ! 3/4 pi <= theta <= 5/4 pi
        else if(((theta-pi*5/4).ge.eps).and.((pi*7/4-theta).ge.eps)) then
            !print *,"zone 3"
            !x(1) = -dsign(1.0d0,dcos(theta))/dtan(theta)
            x(1) = -1/dtan(theta)
            x(2) = -1

        else if(((theta-pi*7/4).ge.eps).or.((theta-pi*0/4).ge.eps)) then
            !print *,"zone 4"
            x(1) = 1
            !x(2) = dsign(1.0d0,dcos(theta))*dtan(theta)
            x(2) = 1*dtan(theta)
        end if

    end function

    !!!
    !! @func : reformat a angle to target range [0,2*pi]
    !! @var  : input angle which can be [-pi,2*pi]
    !!!

    function refmt(angle) result(ans)
        implicit none
        real(8) :: ans,angle,eps
        real(8) :: pi = 4*atan(1.0d0)
        eps=1e-8
        ! if <0 and >-pi
        if((angle.lt.-eps).and.(angle.gt.-eps-pi)) then
            ans=angle+2*pi
        else
            ans=angle
        end if
    end function

   !!!
   !! @func : given index of a unit square, return the range of angle it enclosed
   !! @param: [indx] index number in teng's order
   !!!
   function get_range_t(indx) result(ans)
       implicit none
       real(8) :: ans(2)
       integer :: indx
       real(8) :: pi=4*atan(1.0d0)
       select case(indx)
       case(1)
           ans=[0.d0,pi/2]
       case(3)
           ans=[pi/2,pi]
       case(5)
           ans=[pi,1.5*pi]
       case(7)
           ans=[1.5*pi,2*pi]
       case(2)
           ans=[0.0d0,pi]
       case(4)
           ans=[0.5*pi,1.5*pi]
       case(6)
           ans=[1*pi,2*pi]
       case(8)
           ans=[-0.5*pi,0.5*pi]
       case default
           ans=[0*pi,2*pi]
       end select
    end function

    !!!
    !!  @func : divint the unit square to four corner,four edge and other
    !!!
    function pos_in_sq(pos) result(ans)
        implicit none
        real(rk) :: pos(2)
        real(rk) :: eps=1e-8
        integer :: ans

        if ((1+pos(1)<eps).and.(1+pos(2)<eps)) then
            ans=1
        else if ((1-pos(1)<eps).and.(1+pos(2)<eps)) then
            ans=3
        else if ((1-pos(1)<eps).and.(1-pos(2)<eps)) then
            ans=5
        else if ((1+pos(1)<eps).and.(1-pos(2)<eps)) then
            ans=7
        else if(1+pos(2)<eps) then
            ans=2
        else if(1-pos(1)<eps) then
            ans=4
        else if(1-pos(2)<eps) then
            ans=6
        else if(1+pos(1)<eps) then
            ans=8
        else 
            ans=9
        end if
    end function

    !!!
    !! @func  :: get jocobian determiant given two partial directive at a point
    !! @param :: [gd] gives the two partial direvative has shape (3,2)
    !!!
    function get_jcb_det(gd) result(ans)
        implicit none
        real(rk) :: gd(3,2),gr(3),ans
        real(rk) :: v1(3),v2(3)
        v1=gd(:,1)
        v2=gd(:,2)
        !ans=norm2(cross_product(v1,v2))
        GR(1)=GD(2,1)*GD(3,2)-GD(3,1)*GD(2,2)     ! For 3D normals
        GR(2)=GD(3,1)*GD(1,2)-GD(1,1)*GD(3,2)
        GR(3)=GD(1,1)*GD(2,2)-GD(2,1)*GD(1,2)
        ans=DSQRT(GR(1)*GR(1)+GR(2)*GR(2)+GR(3)*GR(3)) ! 3D JACOBIAN
    end function


end module
