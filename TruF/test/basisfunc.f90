       function Zm(m,z)
        use var_mod
        implicit none
        integer,intent(in) :: m
        real(8) :: Zm,z,wm
        if (m .eq. 0) then
           Zm = dcosh(wk*(h+z))/dcosh(wk*h)
        else
           Wm=wvno1(m)
           Zm=dcos(wm*(h+z))/dcos(wm*h)
        end if
       end function

       function dZm(m,z)
        use var_mod
        implicit none
        integer,intent(in) :: m
        real(8) ::dZm,z,Wm
        if (m .eq. 0) then
           dZm = -wk*dsinh(wk*(h+z))/dcosh(wk*h)
        else
           Wm=wvno1(m)
           dZm=-wm*dsin(wm*(h+z))/dcos(wm*h)
        end if
       end function
       
       function Ym(m,z)
        use var_mod
        implicit none
        integer,intent(in) :: m
        real(8) :: Ym,z,lamda
        if (m .eq. 0) then
           Ym = dsqrt(2.0d0)/2. 
        else
           lamda=m*pi/s
           Ym = dcos(lamda*(z+h))
        end if
       end function

      function dYm(m,z)
        use var_mod
        implicit none
        integer,intent(in) :: m
        real(8) :: dYm,z,lamda
        if (m .eq. 0) then
           dYm = 0.0d0
           else
           lamda=m*pi/s
           dYm = lamda*-dsin(lamda*(z+h))
        end if
       end function


