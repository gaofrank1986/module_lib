!module hi_integrand
module intgrd_funcs
        real(8),private :: pi=4*atan(1.d0)
contains

    !dG/(dndy3)
    function f1(ndim,nf,cosn,drdx,drdn,shap) 
        implicit none
        integer,intent(in) :: ndim,nf
        real(8),intent(in) :: drdx(ndim),drdn,cosn(ndim),shap(*)
        real(8) :: f1(nf)       

        real(8) :: x(ndim),tmp
        integer :: id

        tmp =(3.*DRDX(3)*DRDN-COSN(3))/(4.*pi)    ! GUIG 4.2
        !todo neg sign removed in tmp to match teng's result
        do id =1,8
            f1(id) = tmp*shap(id)
        end do
    end function

    ! dG/dy3

     function f2(ndim,nf,cosn,drdx,drdn,shap)
        implicit none
        integer,intent(in) :: ndim,nf
        real(8),intent(in) :: drdx(ndim),drdn,cosn(ndim),shap(*)
        real(8) :: f2(nf)       

        real(8) :: x(ndim),tmp
        integer :: id
        !should be partial G over partial y3      
        tmp =-DRDx(3)/(4.*pi)    ! GUIG 4.2
        do id =1,8
            f2(id) = tmp*shap(id)
        end do
    end function   

    ! dG/dn

!   subroutine f_integrand4(ndim,nf,cosn,drdx,drdn,shap,fb)
        !implicit none
        !integer,intent(in) :: ndim,nf
        !real(8),intent(in) :: drdx(ndim),drdn,cosn(ndim),shap(*)
        !real(8),intent(out) :: fb(nf)       

        !real(8) :: x(ndim),tmp
        !integer :: id
        !!should be partial G over partial y3      
        !tmp =DRDn/(4.*hi_PI)    ! GUIG 4.2
        !!attention neg sign removed in tmp to match teng's result
        !do id =1,8
            !FB(id) = tmp*shap(id)
        !end do
    !end subroutine   

    !G
!     function f3(ndim,nf,shap,fb)
        !implicit none
        !integer,intent(in) :: ndim,nf
        !real(8),intent(in) :: shap(*)
        !real(8) :: f3(nf)       

        !real(8) :: x(ndim),tmp
        !integer :: id
        !!should be partial G over partial y3      
        !tmp =-1/(4.*hi_PI)    ! GUIG 4.2
        !!attention neg sign removed in tmp to match teng's result
        !do id =1,8
            !F3(id) = tmp*shap(id)
        !end do
    !end function   

 !     SUBROUTINE F_BAR(p,DRDX,COSN,R,DRDN,XI,SHAP,XQ,NF,FB)
          !use param_mod
      !IMPLICIT REAL*8 (A-H,O-Z)
      !type(HSParams) :: p
      !DIMENSION X(p%NDIM),XI(p%NBDM),DRDX(p%NDIM),COSN(p%NDIM),SHAP(*),&
     !&          FB(p%NF),XQ(p%NDIM)
      !!COMMON/RIM_COEF/PI,DLT(3,3),CNU
      !integer :: i
      !real(8) :: tmp,pi=4*atan(1.0d0)

!!     COMPUTING Gij (The first example (2D))
!!      FB(1)=-(1.-2.*DRDX(1)*DRDX(1))/(2.*PI)
!!      FB(2)=2.*DRDX(1)*DRDX(2)/(2.*PI)
!!      FB(3)=-(1.-2.*DRDX(2)*DRDX(2))/(2.*PI)


      !tmp=-(3.*DRDX(3)*DRDN-COSN(3))/(4.*PI)    ! GUIG 4.2
!!     CYLINDER SURFACE ELEMENT (The second example (3D))
      !do i=1,p%nf
      !FB(i)=tmp*SHAP(i)    ! GUIG 4.2
        !end do
        !! todo match sign
        !FB=-FB
        !!=========

!!        IF(R.GT.1.D-18) FB(1)=R*DLOG(R)
!!        IF(R.LT.1.D-18) FB(1)=0.


!!	  FB(1)=-DRDX(1)      ! Karami 4.3

!!        FB(1)=DRDX(1)    
!!        FB(2)=DRDX(2)    

!!	  FB=DRDX

!!!      FORM Tij OF 2D
!!      CON=-1./(4.*PI*(1.-CNU)); PR2=1.-2.*CNU; N=0
!!      DO I=1,NDIM; DO J=1,NDIM; N=N+1
!!       FB(N)=CON*(PR2*(COSN(I)*DRDX(J)-COSN(J)*DRDX(I))+               &
!!     &           (PR2*DLT(I,J)+2*DRDX(I)*DRDX(J))*DRDN)
!!      ENDDO; ENDDO
      
!!     FORM Tij OF 3D
     !! NBDM=NDIM-1; PR1=1.-CNU; PR2=1.-2.*CNU
     !! CON=-1./(4.*NBDM*PI*PR1)
     !! N=0
     !! DO I=1,NDIM; DO J=1,NDIM
     !!  N=N+1
     !!  FB(N)=CON*(DRDN*(PR2*DLT(I,J)+NDIM*DRDX(I)*DRDX(J))-             &
     !!& PR2*(DRDX(I)*COSN(J)-DRDX(J)*COSN(I)))
     !! ENDDO; ENDDO

!!      FB(1)=DRDX(2)    ! Karami 4.1
!!     FB(1)=1.D0  ! GUIG 4.1


!!      FB(1)=DRDX(1)  ! KARAMI3

!!       FB(1)=(1-2.*DRDX(1)**2)
!!       FB(2)=(1-2.*DRDX(2)**2)

!!     COMPUTING Gij
!!      FB(1)=(1.-2.*DRDX(1)*DRDX(1))/(2.*PI)
!!      FB(2)=-2.*DRDX(1)*DRDX(2)/(2.*PI)
!!      FB(3)=(1.-2.*DRDX(2)*DRDX(2))/(2.*PI)

!!     COMPUTING Uijk
!!      CON=1./(4.*PI*(1.-CNU)); PR2=1.-2.*CNU
!!      PHI2=SHAP(3)        
!!      FB(1)=CON*PHI2*(PR2*DRDX(1)+2.*DRDX(1)**3)
!!      FB(2)=CON*PHI2*(PR2*DRDX(1)+2.*DRDX(1)*DRDX(2)**2)
!!      FB(3)=CON*PHI2*(PR2*DRDX(2)+2.*DRDX(2)**3)

!!     COMPUTING Tijk 
!!      E=1.; G=E/(2.*(1.+CNU)); CON=G/(2.*PI*(1.-CNU))
!!      C=1.-2.*CNU; D=1.-4.*CNU; K=2
!!      M=0; DO 10 I=1,NDIM; DO 10 J=I,NDIM; M=M+1
!!      TEM=2.*DRDN*(C*DLT(I,J)*DRDX(K)+CNU*(DLT(I,K)*DRDX(J)+DLT(J,K)*DRDX(I))-   &
!!     & 4.*DRDX(I)*DRDX(J)*DRDX(K))+2.*CNU*(COSN(I)*DRDX(J)*DRDX(K)+COSN(J)*DRDX(I)*       &
!!     & DRDX(K))+C*(2.*COSN(K)*DRDX(I)*DRDX(J)+COSN(J)*DLT(I,K)+COSN(I)*          &
!!     & DLT(J,K))-D*COSN(K)*DLT(I,J)
!!      FB(M)=CON*TEM*SHAP(M)       
!!  10  CONTINUE

      !END

end module 
