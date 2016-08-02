C   
C ******************************************************  
C *                                                    *  
C *  Programme for computing wave height               * 
C *  around a single cylinder                          *                          
C *                                                    *  
C ******************************************************  
C  
        PROGRAM    CYLWH 
        IMPLICIT   REAL*8(A-H,O-Z) 
C   
        INTEGER M 
!         
          REAL*8, ALLOCATABLE:: XP(:),YP(:),ZWHR(:),ZWHI(:), 
     1                        CWHR(:),CWHI(:),DWHR(:),DWHI(:) 
! 
!  CWH: wave height from incident waves 
!  DWH: wave height from diffraction waves 
!  ZWHI: wave height from the total waves 
! 
        COMPLEX*16 CI,DHN,HN,DHNR,HNR,COEF 
        COMPLEX*16 POTN,POTT,WHEC,WHED,WHET 
C 
        DATA G,PI /9.807D0, 3.141592653589793D0 /          
        DATA CI   /(0.0D0,  1.0D0) /   
C  
C------------------------------------------- 
C 
       OPEN(4,   FILE='DPOINT.txt',       STATUS='OLD') 
 
       OPEN(11,   FILE='ZWHEIT.DAT',    STATUS='UNKNOWN') 
       OPEN(12,   FILE='IWHEIT.DAT',    STATUS='UNKNOWN') 
       OPEN(13,   FILE='DWHEIT.DAT',    STATUS='UNKNOWN') 
C 
!        WRITE(6,*) 'Input  H=?   WK=?' 
        READ(4,*)    H, WK 
C 
        WRITE(6, *)  'H=',H,'  WK=',WK 
C 
C -------------------------------------- 
C 
           RAD=1.0d0 
           Z=0.0d0 
C 
c         DO 1000 JJ=1, 41 
 
         READ(4, *)  NPIT 
        WRITE(6, *)  ' Number of points=',NPIT 
 
          ALLOCATE(XP(NPIT),YP(NPIT),ZWHR(NPIT),ZWHI(NPIT), 
     1                 CWHR(NPIT),CWHI(NPIT),DWHR(NPIT),DWHI(NPIT)) 
!!    read data points 
       DO 100 IT=1,  NPIT 
         READ(4, *)  M, XP(IT), YP(IT) 
100   CONTINUE               
! 
! 
!  ======================================= 
! 
         DO 1000  IT=1, NPIT 
          
           X= XP(IT) 
           Y= YP(IT) 
           R= DSQRT(X**2 + Y**2) 
 
           SITA=DATAN2(Y,X)  
!          
C 
          IF (H .LE. 0.0D0) THEN ! infinite deep
            W1=DSQRT(G*WK)! dispersion relationship 
          ELSE               
            W1=DSQRT(G*WK*DTANH(WK*H)) !dispersion relationship
        END IF 
C 
          V1=W1*W1/G 
          TPER=2.D0*PI/W1 
C 
         A   =1.0D0 
         AMP =1.0D0 
         BETA=0.0D0      
C          
           COEF=-CI*G/W1 
         ZFUNC=   DCOSH(WK*(Z+H))/DCOSH(WK*H)  
         DZFNC=WK*DSINH(WK*(Z+H))/DCOSH(WK*H)  
C 
           POTN=(0.0D0, 0.0D0) 
C 
        DO 200  M=0,  30 
           EPS=2.0D0 
           IF(M .EQ. 0) EPS=1.0D0 
C 
         DJN=DBJ(WK, M) 
         DYN=DBY(WK, M) 
         DHN=DCMPLX(DJN,DYN) 
C        
         BJN=BJ(WK, M) 
         BYN=BY(WK, M) 
         HN=DCMPLX(BJN,  BYN) 
C 
         DJNR=DBJ(WK*R, M) 
         DYNR=DBY(WK*R, M) 
         DHNR=DCMPLX(DJNR,DYNR) 
C        
         BJNR=BJ(WK*RAD, M) 
         BYNR=BY(WK*RAD, M) 
         HNR=DCMPLX(BJNR,  BYNR) 
C 
         POTN=POTN-EPS*CI**M*HNR *DJN/DHN*DCOS(M*SITA) 
 
C              
200     CONTINUE 
C 
C         ' DIFFRACTION POTENTIAL' 
C 
          POTN=POTN*ZFUNC*COEF 
            WHED=-CI*W1*POTN/G 
            DWHR(IT)= DBLE(WHED) 
            DWHI(IT)=aIMAG(WHED) 

C         ' INCIDENT POTENTIAL' 
C 
          POTT=COEF*ZFUNC*CDEXP(CI*WK*X) 
            WHEC=-CI*W1*POTN/G 
            CWHR(IT)= DBLE(WHEC) 
            CWHI(IT)=DIMAG(WHEC) 
C 
C         ' TOTAL POTENTIAL' 
C 
            POTN=POTN+POTT 
            WHET=-CI*W1*POTN/G 
            ZWHR(IT)= DBLE(WHET) 
            ZWHI(IT)=DIMAG(WHET) 
 
1000     CONTINUE 
! 
!  ======================================= 
! 
          WRITE(11,*) '  Elevation of total waves' 
          WRITE(12,*) '  Elevation of incident waves' 
          WRITE(13,*) '  Elevation of diffraction waves' 
           
          WRITE(11,1010)  WK,W1 
          WRITE(12,1010)  WK,W1 
          WRITE(13,1010)  WK,W1    
 
        DO IT=1, NPIT 
           WRITE(11,1030)  XP(IT),YP(IT),ZWHR(IT),ZWHI(IT) 
           WRITE(12,1030)  XP(IT),YP(IT),CWHR(IT),CWHI(IT) 
           WRITE(13,1030)  XP(IT),YP(IT),DWHR(IT),DWHI(IT) 
        ENDDO 
! 
!                     
201     FORMAT(I3,2X,6F15.5)  
 
! 
1010    FORMAT(' Wave Number:',F8.3,2x,' Wave Freq.:',F8.3)  
1020    FORMAT(2x,F8.4,4x,6E13.5)  
1030    FORMAT(2x,2F8.4,2x,2E13.5)  
C  
C!1111    FORMAT(//,'  WAVE AMPLITUDE=',F6.2, /, 
C     !1    '  WAVE NUMBER=',F9.5,'   V=',F9.5, /,  
C     !3    '  ANGULAR FREQU.=',F9.5,'  BETA=',F7.3,/)       
C  
        print *,"end"
        END       
C  
! =================================================== 
! 
! 
! =================================================== 
! 
        SUBROUTINE WAVECK(SIGMA,H,WK,IER) 
        IMPLICIT   REAL*8(A-H,O-Z) 
C 
C  H: WATER DEPTH    SIGMA: WAVE FREQUENCY 
C  C: WAVE CELERITY  Wk: wave number 
C 
        DATA G/9.807D0/ 
C 
        IF( SIGMA .LE. 0.0D0 ) GOTO 99 
        IER=0 
        IF( H .GT. 0.0D0) THEN 
          B = G * H 
          Y = SIGMA * SIGMA *H/G 
  12      A=1.0D0/(1.0D0+Y*(0.66667D0+Y*(0.35550D0+ 
     1       Y*(0.16084D0+Y*(0.063201D0+Y*(0.02174D0+ 
     2       Y*(0.00654D0+Y*(0.00171D0+ 
     2       Y*(0.00039D0+Y*0.00011D0))))))))) 
  22      C=SQRT(B/(Y+A)) 
          WK=SIGMA/C 
        ELSE IF( H .LE. 0.0D0) THEN 
          WK=SIGMA**2/G 
        END IF 
        RETURN 
C 
  99    IER=-1 
        WK=0.0D0 
        RETURN 
        END 
 
 
 
 
