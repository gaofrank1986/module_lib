C    ! basis eigen function and their derivatives/integration along z direction
       function Zm(m,z,k,h)
        !use var_mod
        implicit none
        integer,intent(in) :: m
        real(8) :: Zm,z,k,h
        if (m .eq. 0) then
           Zm = dcosh(k*(h+z))/dcosh(k*h)
        else
           !Wm=wvno1(m)
           Zm=dcos(k*(h+z))/dcos(k*h)
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

C *********************************************************
C  Integration of  Z(L)  from -T to 0
C  
C   L=0,   Z(L)= [Cosh(K0(Z+d))/Cosh(K0*d)]
C   L>0,   Z(L)= [Cos(Kj(Z+d))/Cos(Kj*d)],   j=L
C   
C *********************************************************
C
	FUNCTION  ZWL(J)
	USE VAR_MOD
      IMPLICIT  NONE
	REAL*8 ZWL,WL
	INTEGER,INTENT(IN):: J


	IF (J .EQ. 0)  THEN
	 ZWL=(DSINH(WK*H)-DSINH(WK*S))/WK/DCOSH(WK*H)
	ELSE
	 WL=WVNO1(J)
	 ZWL=(DSIN(WL*H)-DSIN(WL*S))/WL/DCOS(WL*H)
	END IF

c	  WRITE(6,*) ' J=',J,'  GWL=',GWL

	RETURN
	END

C
C *********************************************************
C  Integration of  z*Z(L)  from -T to 0 (S=H-T)
C  
C   L=0,   Z(L)= [Cosh(K0(Z+H))/Cosh(K0*H)]
C   L>0,   Z(L)= [ Cos(Kj(Z+H))/Cos(Kj*H)],   j=L
C   
C *********************************************************
C
	FUNCTION  ZZWL(J)
	USE VAR_MOD
      IMPLICIT  NONE
	REAL*8 ZZWL,WL
	INTEGER,INTENT(IN):: J


	IF (J .EQ. 0)  THEN
	 ZZWL=(WK*BH*DSINH(WK*S)+DCOSH(WK*S)-DCOSH(WK*H))/WK/WK/DCOSH(WK*H)
	ELSE
	 WL=WVNO1(J)
	 ZZWL=(WL*BH*DSIN(WL*S)-DCOS(WL*S)+DCOS(WL*H))/WL/WL/DCOS(WL*H)
	END IF

c	  WRITE(6,*) ' J=',J,'  GWL=',GWL

	RETURN
	END

C
C
C *********************************************************
C  Integration of  Z(L)Z(L)  from -H to 0
C  
C   L=0,   Z(L)= [Cosh(K0(Z+H))/Cosh(K0*H)]
C   L>0,   Z(L)= [Cos(Kj(Z+H))/Cos(Kj*H)],   j=L
C   
C *********************************************************
C
	FUNCTION  GWL(J)
	USE VAR_MOD
      IMPLICIT  NONE
	REAL*8 GWL,WL
	INTEGER,INTENT(IN):: J


	IF (J .EQ. 0)  THEN
	 GWL=H/2.0d0+DSINH(2.0d0*WK*H)/4.0d0/WK
	 GWL=GWL/DCOSH(WK*H)**2
	ELSE
	 WL=WVNO1(J)
	 GWL=H/2.0d0+DSIN(2.0d0*WL*H)/4.0d0/WL
	 GWL=GWL/DCOS(WL*H)**2
	END IF

c	  WRITE(6,*) ' J=',J,'  GWL=',GWL

	RETURN
	END

C
C
C *********************************************************
C  Integration of  Z(j)*Y(m)  from -d to -T	
C   j=0,       Z(L)=Cosh(K0(Z+d))/Cosh(K0*d)
C   j>0,       Z(L)=Cos(Kj(Z+d)) /Cos(Kj*d),  
C   m=0,       Y(M)=SQRT(2)/2,          Lm=0
C   m>0,       Y(M)=Cos(Lm(Z+d)),       Lm=m*pi/S
C
C *********************************************************
C
	FUNCTION  HCCLM(J,M)
	USE VAR_MOD
      IMPLICIT  NONE
	REAL*8 HCCLM,VM,WL
	INTEGER,INTENT(IN):: J,M

C

	 VM=WVZJ1(M) 
C
	IF 	(M .EQ. 0)  THEN
	 IF (J .EQ. 0)  THEN
	  HCCLM=DSINH(WK*S)/WK/DCOSH(WK*H)/SQRT(2.0D0)
	 ELSE
	  WL=WVNO1(J)
	  HCCLM=DSIN(WL*S)/WL/DCOS(WL*H)/SQRT(2.0D0)
	 END IF
C
	ELSE
	 IF (J .EQ. 0)  THEN
	  HCCLM=(-1.0d0)**M*WK*DSINH(WK*S)/(WK*WK+VM*VM)
	  HCCLM=HCCLM/DCOSH(WK*H)
	 ELSE
	  WL=WVNO1(J)
	  HCCLM=(-1.0d0)**M*WL*DSIN(WL*S)/(WL*WL-VM*VM)
	  HCCLM=HCCLM/DCOS(WL*H)
	 END IF
	END IF
C
	RETURN
	END
