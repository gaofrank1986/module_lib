module param_mod
    use geo2D_mod
    type :: Coef
        real(8),allocatable,dimension(:) :: G,H
        real(8),allocatable,dimension(:,:) :: B
    contains
        procedure :: init=>init_coef
    end type

    type :: Elem
        integer :: nsub,nnode,ndim,nbdm
        real(rk) :: ck(3,8)
        type(elem2D) :: mapped
    contains
        procedure :: map2D
        procedure :: get_SF
        procedure :: get_nrml_at
        !procedure,private :: get_mapped_CDL
    end type

    type :: params
        integer :: ndim,nbdm,node,npowg,nnode,nelem,nf
        integer :: ngl=-10
        integer :: ngr=-10
        integer :: npw = 4
        real(8) :: xp(3),xip(2),xiq(2),beta
        type(Coef) :: mat
        real(8),dimension(10) :: gpl,gwl,gpr,gwr
    contains 
        procedure :: pprint    
        procedure :: init_mat
    end type
    private :: get_mapped_CDL,SHAPEF,DSHAPE
contains
    subroutine init_mat(this)
        class(params) ::this
        call this%mat%init(this%npw,this%npowg,this%nf)
    end subroutine

    subroutine init_coef(this,npw,npowg,nf)
        class(Coef) :: this
        integer :: npw
        allocate(this%G(0:NPOWG),this%H(0:NPW))
        allocate(this%B(0:11,0:nf))
        this%B=0.0d0
        this%G=0.0d0
        this%H=0.0d0
    end subroutine

        
        subroutine pprint(this)
            class(params) :: this
            print '(5a6)',"ndim","nbdm","node","nnode","nelem","nf"
            print '(5i6)',this%ndim,this%nbdm,this%node,this%nnode,this%nelem,this%nf
            print *,"npowg",this%npowg
            print *,"src2D",this%xip
            print *,"src3D",this%xp
            print *,"beta",this%beta
            print *,"nf",this%nf
        end subroutine


        function get_mapped_CDL() result(ans)
            implicit none
            real(8),dimension(18) :: ans
            ans =[-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,-1.,0.,0.,0.]
        end function

        function map2D(this,x) result(ans)
            implicit none
            class(Elem) ::this
            real(8),dimension(3) :: ans
            real(8) :: x(2),ri(3),sf(8),pt3D(3)
            pt3D = 0.0d0
            call shapef(3,8,x,this%ck,pt3D,ans,get_mapped_CDL(),sf)
        end function



        function get_SF(this,x) result(res)
            implicit none
            class(Elem) ::this
            real(8),dimension(8) :: res
            real(8) :: x(2),pt3D(3),ri(3)
            pt3D = 0.0d0
            call shapef(3,8,x,this%ck,pt3D,ri,get_mapped_CDL(),res)
        end function

!        function get_DSF(this,x) result(res)
            !implicit none
            !class(Elem) ::this
            !real(8),dimension(8) :: res
            !real(8) :: x(2),dir3D(3),ri(3),jacb,gd(3,8)
            !dir3D = 0.0d0
            !call dshapef(3,2,8,x,this%ck,dir3D,FJCB,get_mapped_CDL(),GD)

            !!SUBROUTINE DSHAPE(NDIM,NBDM,NODE,X,CK,COSN,FJCB,C,GD)
        !end function

        ! @var :[cosn]  normalized nrml vector
        ! @var :[fjcb]  jacobian determinant
        ! @var :[gd] two tangent vector
        subroutine get_nrml_at(this,x,cosn,fjcb,gd)
            implicit none
            class(Elem) ::this
            real(8) :: cosn(3),fjcb,gd(3,2)
            real(8) :: x(2)
            call dshape(3,2,8,x,this%ck,cosn,FJCB,get_mapped_CDL(),GD)
        end subroutine

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

end module



