      
      PROGRAM SIEPPEM
          use param_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RIM_COEF/PI,DLT(3,3),CNU
      DIMENSION XP(3)
      ALLOCATABLE CD(:,:),LNDB(:,:),NSEL(:),XIS(:,:),VINT(:)
	  DATA NGR,NGL/-10,-10/  !If |NGR| or |NGL| is bigger than 10, must change subroutine GAUSSV
      type(HSParams) :: par
      OPEN(5,FILE='sieppem.dat',STATUS='OLD')
      OPEN(7,FILE='SIEPPEM.OUT',STATUS='UNKNOWN')
      READ(5,*)NDIM,NTP,NBE,NODE,BETA,NF        ! Card set 1
      ALLOCATE(CD(NDIM,NTP),LNDB(NODE,NBE),NSEL(NBE),XIS(2,NBE),        &
     &         VINT(NF))
 !    Assign values to COMMON BLOCK variables
      DLT=RESHAPE((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
      PI=4.D0*DATAN(1.D0); CNU=0.3; TOLGP=1.D-14
      IF(NDIM.EQ.2) NGL=1
 !    Input nodal coordinates and element connectivity
      DO 10 IP=1,NTP
  10  READ(5,*)NP,(CD(I,NP),I=1,NDIM)                 ! Card set 2
      DO 20 IE=1,NBE
  20  READ(5,*)NE,(LNDB(ID,NE),ID=1,NODE),NSEL(NE)    ! Card set 3
      READ(5,*)(XP(I),I=1,NDIM)                       ! Card set 4  
	DO 30 IE=1,NBE
30    IF(NSEL(IE).LT.0) READ(5,*)(XIS(I,IE),I=1,NDIM-1)  ! Card set 5
      !--------------------------------------------------------------------
!      CALL NASTRAN_MESH(32,'BOUNDARY MESH AND INTERNAL POINTS',NDIM,    &
!     &                  NDIM-1,NTP,NBE,CD, LNDB,NODE,1,1)
      !--------------------------------------------------------------------
 !    Evaluate boundary integrals 
      par%ndim = ndim
      par%nbdm = ndim-1
      par%node = node
      par%beta = beta
      par%nf = nf
      par%xp = xp
      par%xip = xis(:,1)
      par%nnode = ntp
      par%nelem = NBE

      !CALL RIM_ELEMS(NDIM,NODE,BETA,NBE,CD,LNDB,NSEL,XP,XIS,NGR,NGL,    &
     !&               TOLGP,NF,VINT)
       CALL RIM_ELEMS(par,CD,LNDB,NSEL,TOLGP,VINT)  
       write(*,*)'              Results :'
      !WRITE(*,'(1P,8(3E14.6/,6X))')(VINT(I),I=1,NF)
      !WRITE(7,'(1P,8(3E14.6/,6X))')(VINT(I),I=1,NF)
      write (*,'(8f14.8)') vint(1),vint(5),vint(2)&
                    ,vint(6),vint(3),vint(7),vint(4),vint(8)
      print *,sum(vint)
      DEALLOCATE (CD,LNDB,NSEL,XIS,VINT)
      STOP
      end program

       SUBROUTINE RIM_ELEMS(p,CD,LNDB,NSEL,TOLGP,VINT)
          use hs_intg
          use param_mod
          use matrix_mod
 IMPLICIT REAL*8 (A-H,O-Z)
      type(HSParams) :: p
      type(HSElem) :: e
      DIMENSION CD(p%NDIM,*),LNDB(p%NODE,p%nelem),NSEL(p%nelem),VINT(p%NF),     &
     & CDL(18),NODEF(12),CK(3,p%NODE),     &
     &          V1E(p%NF)
      EXTERNAL INT_ELEM
            DATA NODEF/1,2,5, 2,3,6, 3,4,7, 4,1,8/,NPW/4/
      type(Matrix2D) :: mat
     
      IF(p%NDIM.EQ.2) p%NPOWG=(p%NODE/3)*2               ! 0,  2
      IF(p%NDIM.EQ.3) p%NPOWG=p%NODE/2+(p%NODE/9)*2        ! 2,  4,  6 
      call p%init_mat()
      NSUB=2*p%NBDM
      e%nsub = nsub
      NDSID=2+p%NODE/8
      VINT=0.
      DO 20 IE=1,p%nelem
      DO 10 ID=1,p%NODE
          e%ck(1:3,id) = CD(1:p%NDIM,LNDB(ID,IE))
10    CK(1:p%NDIM,ID)=CD(1:p%NDIM,LNDB(ID,IE))
      call e%mapped%get_std()
      !call mat%init(e%ck)
      V1E=0.
      IF(NSEL(IE).EQ.0) THEN    ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
     !  CALL ADAPTINT_ELEM(NDIM,NBDM,NODE,BETA,CP0,XP,CK,CDL,NF,V1E,     &
     !&                    NGR,NGL,GPR,GWR,GPL,GWL,TOLGP,BETA,INT_ELEM)
      ELSE     ! EVALUATE INTEGRAL OVER SINGULAR ELEMENT
       CALL SINGULAR_ELEM(e,p,TOLGP,NDSID,NSEL(IE),V1E)
!       print *,NSEL
!       print *,XIS
      ENDIF
      VINT=VINT+V1E
20    CONTINUE
      !call p%pprnnnint()
      !call mat%pprint()
      print *,nsel(1)
      END
