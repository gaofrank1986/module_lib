
      
      
    ! @var [x] local pos
    ! @var : [xp] for calc distance
    ! @Var : [C] std 2D mtx
    ! @var : [RI] distance vector
    ! @var : [CK] glb pos
      SUBROUTINE SHAPEF(NDIM,NODE,X,CK,XP,RI,C,SP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SP(NODE),CK(3,*),XP(*),RI(*),C(*),X(*)
      IF(NODE.GT.3) GOTO 4
  !                  2-noded line element
      SP(1)=0.5*(1.-X(1)); SP(2)=0.5*(1.+X(1))
      IF(NODE.EQ.2) GOTO 50
  !                  3-noded line element
      SP(1)=-X(1)*SP(1); SP(2)=X(1)*SP(2); SP(3)=1.-X(1)*X(1)
      GOTO 50
  !                  4-noded quadrilateral element
  4   DO I=1,4
       SP(I)=0.25*(1.+C(2*I-1)*X(1))*(1.+C(2*I)*X(2))
      ENDDO
      IF(NODE.EQ.8) THEN
  !                  8 noded-element (square element)
       DO 15 I=1,4; L=2*I-1; SP(I)=SP(I)*(C(L)*X(1)+C(L+1)*X(2)-1.D0)
       WL=C(L+8)*X(1)+C(L+9)*X(2)
  15   SP(I+4)=.5D0*(WL+1.D0)*(1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2)
      ELSEIF(NODE.EQ.9) THEN
  !                  9 noded-element (square element)
       DO 20 I=1,4; L=2*I-1; SP(I)=SP(I)*C(L)*X(1)*C(L+1)*X(2)
       WL=C(L+8)*X(1)+C(L+9)*X(2)
  20   SP(I+4)=0.5D0*WL*(WL+1.D0)*(1.D0-(C(L+8)*X(2))**2-               &
     &                  (C(L+9)*X(1))**2)
       SP(9)=(1.D0-X(1)*X(1))*(1.D0-X(2)*X(2))
      ENDIF
  !                  Calculate r and its vector components
  50  RQ2=0.; DO 70 I=1,NDIM; RI(I)=-XP(I)
      DO ID=1,NODE; RI(I)=RI(I)+SP(ID)*CK(I,ID); ENDDO
  70  RQ2=RQ2+RI(I)*RI(I)
      END

      SUBROUTINE DSHAPE(NDIM,NBDM,NODE,X,CK,COSN,FJCB,C,GD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),CK(3,*),DN(2,NODE),GD(3,*),COSN(*),C(*),GR(3)
      IF(NODE.GT.3) GOTO 5
      DN(1,1)=-0.5; DN(1,2)=0.5        ! 2-noded line element
      IF(NODE.EQ.2) GOTO 30
      DN(1,1)=-0.5*(1.-2.*X(1))        ! 3-noded line element
      DN(1,2)=0.5*(1.+2.*X(1)); DN(1,3)=-2.*X(1); GOTO 30
   5  DO 10 I=1,4; I0=2*(I-1)          ! 4-noded quadr. element
      DN(1,I)=0.25*C(I0+1)*(1.+C(I0+2)*X(2))
  10  DN(2,I)=0.25*C(I0+2)*(1.+C(I0+1)*X(1))
      IF(NODE.EQ.8) THEN
!                  8 noded-element (square element)
       DO 15 I=1,4; L=2*I-1; S=C(L)*X(1)+C(L+1)*X(2)-1.D0
       DN(1,I)=DN(1,I)*S+0.25D0*(1.D0+C(L)*X(1))*(1.D0+C(L+1)*X(2))*C(L)
       DN(2,I)=DN(2,I)*S+0.25D0*(1.D0+C(L)*X(1))*(1.D0+C(L+1)*X(2))*    &
     &                 C(L+1)
       S=1.D0+C(L+8)*X(1)+C(L+9)*X(2)      
       T=1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2
       DN(1,I+4)=0.5D0*C(L+8)*T-C(L+9)*C(L+9)*X(1)*S
  15   DN(2,I+4)=0.5D0*C(L+9)*T-C(L+8)*C(L+8)*X(2)*S
      ELSEIF(NODE.EQ.9) THEN
  !                 9 noded-element (square element)
       DO 20 I=1,4; L=2*I-1
       DN(1,I)=DN(1,I)*C(L+1)*X(2)*(1.+2.D0*C(L)*X(1))
       DN(2,I)=DN(2,I)*C(L)*X(1)*(1.+2.D0*C(L+1)*X(2))
       S=C(L+8)*X(1)+C(L+9)*X(2)      
       SS=S*(1.D0+S)
       T=(1.D0-(C(L+8)*X(2))**2-(C(L+9)*X(1))**2)*(1.D0+2.D0*S)
       XI2=C(L+8)*C(L+8); ET2=C(L+9)*C(L+9)
       DN(1,I+4)=0.5D0*(C(L+8)*T-2.D0*ET2*X(1)*SS)
  20   DN(2,I+4)=0.5D0*(C(L+9)*T-2.D0*XI2*X(2)*SS)
       DN(1,9)=-2.D0*X(1)*(1.D0-X(2)*X(2))
       DN(2,9)=-2.D0*X(2)*(1.D0-X(1)*X(1))
      ENDIF
  30  DO 50 I=1,NDIM; DO 50 J=1,NBDM; GD(I,J)=0.; DO 50 ID=1,NODE
50    GD(I,J)=GD(I,J)+DN(J,ID)*CK(I,ID)
      IF(NDIM.EQ.NBDM) THEN
       GOTO (51,52,53),NDIM
  51   FJCB=DABS(GD(1,1))
       RETURN          ! For internal cell integral
  52   FJCB=GD(1,1)*GD(2,2)-GD(1,2)*GD(2,1)
       RETURN
  53   FJCB=GD(1,1)*(GD(2,2)*GD(3,3)-GD(3,2)*GD(2,3))                   &
     &     -GD(1,2)*(GD(2,1)*GD(3,3)-GD(3,1)*GD(2,3))+                  &
     &      GD(1,3)*(GD(2,1)*GD(3,2)-GD(2,2)*GD(3,1))
       RETURN
      ENDIF
      IF(NODE.GT.3) GOTO 60
      GR(1)=GD(2,1); GR(2)=-GD(1,1)             ! For 2D normals
      FJCB=DSQRT(DOT_PRODUCT(GD(1:NDIM,1),GD(1:NDIM,1))) ! LINE JACOBIAN
      GOTO 70
  60  GR(1)=GD(2,1)*GD(3,2)-GD(3,1)*GD(2,2)     ! For 3D normals
      GR(2)=GD(3,1)*GD(1,2)-GD(1,1)*GD(3,2)
      GR(3)=GD(1,1)*GD(2,2)-GD(2,1)*GD(1,2)
      FJCB=DSQRT(GR(1)*GR(1)+GR(2)*GR(2)+GR(3)*GR(3)) ! 3D JACOBIAN
  70  COSN(1:NDIM)=GR(1:NDIM)/FJCB              ! 2D and 3D Normal
      END

