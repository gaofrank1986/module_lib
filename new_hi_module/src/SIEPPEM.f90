module hs_intg
    use param_mod
    integer,parameter,private :: rk=8

contains
      subroutine singular_elem(e,p,tolgp,ndsid,nsp,v1e)
          use matrix_io_mod
          use gaussm3
          implicit none
          type(HSParams) :: p
          type(HSElem)   :: e

          integer :: KSB(4)
          real(8) :: XIQ(p%NBDM),CSUB(2,2),sf0(p%NODE),XI0(p%NBDM),V1E(p%NF),   &
              &          RINT(p%NF)

          real(8),dimension(18) :: CDL = [-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,          &
              &         -1.,0.,0.,0.]

          real(8) :: tol,wfa,fk,tolgp,sv,vlc2,rs,xiel,sum,fjcbl,rhoq,comt,drdnp
          integer :: nsp,isid,ks,id,jd,ixi,iswp,isub,igl,ndsid,npw


          v1e=0.0d0
          tol=1.d-5
          ksb=(/-2,1,2,-1/)
          if(nsp.lt.0) then

              p%xp = e%map2D(p%xip)
              !sf0 = e%get_SF(p%xip)
          ENDIF

          if(nsp.gt.0) then
              p%xp(1:p%ndim)=e%ck(1:p%ndim,nsp) 
              p%xip=e%mapped%nodes(nsp)%getArray()
          end if
          !call p%pprint()




          IF(p%NDIM.EQ.3) then          
              !   Evaluate 2D SINGULAR INTEGRALS
              !   Evaluate 3D SINGULAR INTEGRALS    
              call gauleg(iabs(p%ngl),p%gpl,p%gwl)

              IF(NSP.GT.0) p%XIP=CDL(2*NSP-1:2*NSP)
              wfa=dsqrt(p%beta*2.d0/3.d0+0.4d0)*dlog(dabs(tolgp)/2.d0)
              fk=3.d0/8.d0*(-10.d0*dble(iabs(p%ngl))/wfa-1.d0)**(4.d0/3.d0)
              DO 80 ISID=1,e%NSUB

                  KS=KSB(ISID)
                  IF(DABS(p%XIP(IABS(KS))-DBLE(KS)/DABS(DBLE(KS))).LT.TOL) GOTO 80  

                  !   Set up a sub-element from side ISID and the singular point NSP
                  DO ID=1,2; 
                      !JD=NODEF(3*(ISID-1)+ID)
                      jd = e%mapped%edges(isid,id)
                      CSUB(1:2,ID)=CDL(2*JD-1:2*JD)
                  ENDDO
                  IXI=1; IF(ISID/2*2.EQ.ISID) IXI=2; ISWP=3-IXI
                  SV=DSIGN(1.D0,CSUB(IXI,2)-CSUB(IXI,1))
                  VLc2=(CSUB(ISWP,2)-p%XIP(ISWP))**2
                  XI0=CSUB(:,1)
                  XIQ(ISWP)=XI0(ISWP)
                  DO ISUB=1,500
                      RS=DSQRT(DOT_PRODUCT(XI0-p%XIP,XI0-p%XIP))
                      IF(SV*XI0(IXI)+1.D-8.GE.SV*CSUB(IXI,2)) GOTO 80
                      IF(SV*XI0(IXI).GE.SV*p%XIP(IXI)) THEN
                          XIEL=FK*RS
                      ELSE
                          SUM=VLc2*(1.-FK*FK)+(p%XIP(IXI)-XI0(IXI))**2
                          IF(SUM.LT.0.D0) THEN
                              XIEL=SV*(p%XIP(IXI)-XI0(IXI))
                          ELSE
                              XIEL=FK*(FK*SV*(p%XIP(IXI)-XI0(IXI))-DSQRT(SUM))/(FK*FK-1.)
                          ENDIF
                      ENDIF
                      IF(SV*(XI0(IXI)+SV*XIEL)+1.D-8.GT.SV*CSUB(IXI,2))                 &
                          &   XIEL=SV*(CSUB(IXI,2)-XI0(IXI))
                      FJCBL=0.5D0*SV*XIEL
                      DO IGL=1,iabs(p%NGL)
                          XIQ(IXI)=XI0(IXI)+FJCBL*(1.D0+p%GPL(IGL))
                          RHOQ=DSQRT(DOT_PRODUCT(XIQ-p%XIP,XIQ-p%XIP))
                          DRDNP=DABS(XIQ(ISWP)-p%XIP(ISWP))/RHOQ
                          CALL COEFS_GH(e,p,XIQ)
                          CALL INT_RHO(e,p,XIQ,RINT)
                          COMT=DABS(FJCBL)*p%GWL(IGL)*DRDNP/RHOQ
                          V1E=V1E+COMT*RINT
                      end do
                      XI0(IXI)=XI0(IXI)+SV*XIEL
                  end do
                      80  CONTINUE
        end if
                  END subroutine

      SUBROUTINE INT_RHO(e,p,XIQ,RINT)      
          use param_mod
          implicit none
      type(HSParams) :: p
      type(HSElem) ::e
      real(8) :: xiq(p%nbdm),slop(p%nbdm),xi(p%nbdm),rint(p%nf)
      integer :: npowf,beta1,k,pw,nbeta
      real(8) :: rhoq

      SLOP=XIQ-p%XIP
      RHOQ=DSQRT(DOT_PRODUCT(SLOP,SLOP))
      SLOP=SLOP/RHOQ
      NPOWF=3+2.1214*RHOQ    ! NPOWF IS FROM 3 TO 9
      IF(NPOWF.LT.p%BETA-p%NBDM) NPOWF=p%BETA-p%NBDM
      !   DETERMINE COEFFICIENTS Bn
      XI=p%XIP+RHOQ*SLOP
      CALL COEF_B(e,p,NPOWF,XI,SLOP,RHOQ)
      !   EVALUATE RESULTS USING Eqs.(3-6-44) AND (3-6-63)
      BETA1=p%BETA-p%NBDM+1
      RINT=0.D0  

      DO  K=0,int(p%BETA)-p%NDIM
          PW=p%BETA-K-p%NBDM
          RINT=RINT-(1.D0/RHOQ**PW-p%mat%H(INT(PW)))/PW*p%mat%B(K,:)
      end do


      NBETA=p%BETA-p%NBDM  
      IF(NBETA.GE.0) RINT=RINT+p%mat%B(NBETA,:)*DLOG(RHOQ/p%mat%H(0))

      DO K=int(p%BETA)-p%NBDM+1,NPOWF
       PW=K-p%BETA+p%NBDM
       RINT=RINT+RHOQ**PW/PW*p%mat%B(K,:)
      ENDDO
      END




      SUBROUTINE COEFS_GH(e,p,XIQ)
      use param_mod
      implicit none
      type(HSParams) :: p 
      type(HSElem) :: e 
      real(8) :: XIQ(p%NBDM),COEFC(0:p%NPW)
      integer :: n,i

      ! DETERMINE COEFFICIENTS Gn
      CALL COEF_G(e,p,XIQ)
      ! DETERMINE COEFFICIENTS Cn USING Eq. (3-6-37)
      COEFC(0)=DSQRT(p%mat%G(0))
      DO N=1,p%NPW
       COEFC(N)=0.D0
       DO I=1,N-1
        COEFC(N)=COEFC(N)-COEFC(I)*COEFC(N-I)     
       ENDDO
       IF(N.LE.p%NPOWG)COEFC(N)=p%mat%G(N)+COEFC(N)   ! Eq. (3-6-37)
       COEFC(N)=COEFC(N)/(2.D0*COEFC(0))
      ENDDO
      ! DETERMINE COEFFICIENTS Hi USING Eq.(3-6-47a)
      COEFC(1:p%NPW)=COEFC(1:p%NPW)/COEFC(0)    ! Eq.(3-6-47b)
      
      p%mat%H(0)=1.D0/COEFC(0)
      p%mat%H(1)=COEFC(1)
      p%mat%H(2)=2.*COEFC(2)-COEFC(1)*COEFC(1)
      p%mat%H(3)=3.*COEFC(3)-3.*COEFC(1)*COEFC(2)+COEFC(1)**3
      p%mat%H(4)=4.*COEFC(4)+4.*COEFC(1)*COEFC(1)*COEFC(2)                &
     &        -4.*COEFC(1)*COEFC(3)-2.*COEFC(2)*COEFC(2)                &
     &        -COEFC(1)**4
      END

      SUBROUTINE COEF_G(e,p,XIQ)
      ! THIS ROUTINE DETERMINES COEFFICIENTS Gm USING Eqs.(3-6-58)~(3-6-60)
      use param_mod
      !IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
      type(HSParams) :: p
      type(HSElem) :: e
      real(8) :: XIQ(p%NBDM),COSN(p%NDIM),        &
     &          GCD(p%NDIM,p%NBDM),XI(p%NDIM),SLOP(p%NBDM),RI(p%NDIM),     &
     &          SHAP(p%NODE),RMAT(p%NPOWG,p%NPOWG)
      real(8) :: xbar2,xxbar,rhoq,vk,rho,fjcb,r2
      integer :: ip,jp

      IF(p%NODE.EQ.2) THEN         ! Eq.(24)
       p%mat%G(0)=0.25*((e%CK(1,1)-e%CK(1,2))**2+(e%CK(2,1)-e%CK(2,2))**2)  
      ELSEIF(p%NODE.EQ.3) THEN             ! FORM Eq.(3-6-33)
       RI=e%CK(:,1)+e%CK(:,2)-2.*e%CK(:,3)     ! Xbar
       XI=e%CK(:,1)-e%CK(:,2)                ! X1i-X2i 
       XBAR2=RI(1)*RI(1)+RI(2)*RI(2)     ! Xbar^2
       XXBAR=XI(1)*RI(1)+XI(2)*RI(2)
       p%mat%G(0)=0.25*(XI(1)*XI(1)+XI(2)*XI(2))-XXBAR*p%XIP(1)             &
     &         +XBAR2*p%XIP(1)*p%XIP(1)
       p%mat%G(1)=(-0.5*XXBAR+XBAR2*p%XIP(1))*XIQ(1)
       p%mat%G(2)=0.25*XBAR2
      ELSE      ! QUADRILATERAL ELEMENTS
       SLOP=XIQ-p%XIP
       RHOQ=DSQRT(DOT_PRODUCT(SLOP,SLOP))
       SLOP=SLOP/RHOQ
       VK=RHOQ/DBLE(p%NPOWG)
       DO 20 IP=0,p%NPOWG
       RHO=VK*DBLE(IP)
       XI(1:p%NBDM)=p%XIP+RHO*SLOP
       
       shap=e%get_SF(xi)
       r2 = norm2(e%map2D(xi)-p%xp)
       IF(IP.NE.0) GOTO 10
       call e%get_nrml_at(xi,cosn,fjcb,gcd)

       RI=MATMUL(GCD,SLOP); p%mat%G(0)=DOT_PRODUCT(RI,RI)  ! Eq.(3-6-58)
       GOTO 20
 10    p%mat%G(IP)=((R2/RHO)**2-p%mat%G(0))/RHO     ! Eq.(3-6-60)
       RMAT(IP,1)=1.D0        ! Eq.(3-6-42)
       DO JP=2,p%NPOWG; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20    CONTINUE
       CALL INVSOLVR(p%NPOWG,p%NPOWG,RMAT,p%NPOWG,-1)       ! INVERSE MATRIX [R] IN Eq.(3-6-59)
       p%mat%G(1:p%NPOWG)=MATMUL(RMAT,p%mat%G(1:p%NPOWG))     ! SOLVING Eq.(3-6-59) FOR {G}
      ENDIF
      END

      SUBROUTINE COEF_B(e,p,NPOWF,XIQ,SLOP,RHOQ)
        use param_mod
      !IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
      type(HSParams) ::p
      type(HSElem) ::e
      
      integer,intent(in) :: npowf

      !! drdx,cosn,gcd,xiq,xi,slop,ri,xi,rmat,coefg,shap,fq,a
      real(rk) :: DRDX(p%NDIM),COSN(p%NDIM),          &
          &          gcd(p%ndim,p%nbdm),xiq(p%nbdm),xi(p%nbdm),slop(p%nbdm),    &
          &          RI(p%NDIM),X(p%NDIM),RMAT(NPOWF,NPOWF),COEFG(0:p%NPOWG),      &
          &   shap(p%node),fq(p%nf),a(p%ndim)
      !real(8),dimension(:,:),allocatable :: B
      real(rk) ::vk,rho,rhoq,r,robar,drdn,gm,fjcb
      integer :: j,ip,m,i,jp


      VK=RHOQ/DBLE(NPOWF)

      DO 20 IP=0,NPOWF
      RHO=VK*DBLE(IP)
      XI=p%XIP+RHO*SLOP
      ROBAR=p%mat%G(0); DO M=1,p%NPOWG; ROBAR=ROBAR+p%mat%G(M)*RHO**M; ENDDO
      ROBAR=DSQRT(ROBAR)       ! Eq.(3-6-28)

      ri = e%map2D(xi)-p%xp
      shap = e%get_sf(xi)
      r = norm2(ri)
      X=p%XP+RI

      call e%get_nrml_at(xi,cosn,fjcb,gcd)
      IF(RHO.GT.1.0D-10)THEN    
          DRDX=RI/R
      ELSE
          DO I=1,p%NDIM; A(I)=0.D0
              DO J=1,p%NBDM            
                  A(I)=A(I)+GCD(I,J)*SLOP(J)
              end do
          end do
          GM=DSQRT(DOT_PRODUCT(A,A))
          DRDX=A/GM                 ! Eq.(3-6-74)
      ENDIF 
      ! todo use self-defined normal
      !==============
      if(p%user_nrml) then
          cosn=matmul(e%nk,shap)
      end if 
      !================
      DRDN=DOT_PRODUCT(COSN,DRDX)
      !CALL F_BAR1(p,DRDX,COSN,R,DRDN,XI,SHAP,X,NF,FQ)
      FQ=p%f_bar(p%ndim,p%nf,cosn,drdx,drdn,shap)


      p%mat%B(IP,:)=FQ*FJCB/ROBAR**p%BETA
      IF(IP.EQ.0) GOTO 20
      p%mat%B(IP,:)=(p%mat%B(IP,:)-p%mat%B(0,:))/RHO
      RMAT(IP,1)=1.D0
      DO JP=2,NPOWF; RMAT(IP,JP)=RMAT(IP,JP-1)*RHO; ENDDO
 20   CONTINUE
      CALL INVSOLVR(NPOWF,NPOWF,RMAT,NPOWF,-1)
      p%mat%B(1:NPOWF,:)=MATMUL(RMAT,p%mat%B(1:NPOWF,:))    ! Bk IN Eq.(3-6-62)
      END subroutine


      SUBROUTINE INVSOLVR(NROW,NCOL,A,N,INDIC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NROW,NCOL),W(N),IROW(N)
      EPS=1.D-8             ! The tolerance of the minimum pivot
      DO 40 K=1,N; KM1=K-1; PIVOT=0.D0
      DO 20 I=1,N; IF(K.EQ.1) GOTO 10
      DO ISCAN=1,KM1; IF(I.EQ.IROW(ISCAN)) GOTO 20; ENDDO
  10  IF(DABS(A(I,K)).LE.DABS(PIVOT)) GOTO 20
      PIVOT=A(I,K); IROW(K)=I
  20  CONTINUE
      IF(DABS(PIVOT).GT.EPS) GOTO 30; STOP 9999
  30  IROWK=IROW(K)
      A(IROWK,1:NCOL)=A(IROWK,1:NCOL)/PIVOT; A(IROWK,K)=1.D0/PIVOT
      DO 40 I=1,N; AIK=A(I,K); IF(I.EQ.IROWK) GOTO 40
      A(I,1:NCOL)=A(I,1:NCOL)-AIK*A(IROWK,1:NCOL); A(I,K)=-AIK/PIVOT
  40  CONTINUE
      IF(INDIC.LT.0) GOTO 60
      DO 50 IX=N+1,NCOL; W(1:N)=A(IROW(1:N),IX); DO 50 I=1,N
  50  A(I,IX)=W(I)
      IF(INDIC.GT.0) RETURN
  60  DO 70 J=1,N; W(1:N)=A(IROW(1:N),J); DO 70 I=1,N
  70  A(I,J)=W(I)
      DO 80 I=1,N; W(IROW(1:N))=A(I,1:N); DO 80 J=1,N
  80  A(I,J)=W(J)
      END



      end module
