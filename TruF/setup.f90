! 
! *************************************************************
! *  Setup AMJI,BMIJ,FMJ,EMI matrices                         *
! *                                                           *
! *************************************************************
!
   SUBROUTINE SETUP
   use VAR_mod

   IMPLICIT NONE
   INTEGER    L,M,N,MJ
   REAL*8     BJ,DBJ,BSL,DBL,BY,DBY,BIE,DBIE,BK,DBK,HCCLM,VID,GWL
   COMPLEX*16 PN,HKL,DHKL
! 
! ----------------------------------
!
      EMI=DCMPLX(0.0d0, 0.0D0)
!
       DO 500 N=0, NNOD
       Print *,'N=',N
!
!   LLOD: eigen-mode in outside domain 
!   MMOD: eigen-mode in inside domain  
!   NNOD: Fourier mode number                                             
!
! --- Setup Matrix AMJI and FMJ
!
        BSL =BJ(WK*A,N)
        DBL =DBJ(WK*A,N)
        HKL =DCMPLX(BSL, BY(WK*A, N))
        DHKL=DCMPLX(DBL, DBY(WK*A,N))        
!
      DO 100 M=0, MMOD
	    L=0
		FMJ(M+1,N)=2.0D0*BSL/S*HCCLM(0,M)
        AMJI(M+1,1,N)=2.0D0*HKL/S*HCCLM(0,M)
       DO 80 L=1, LLOD
        AMJI(M+1,L+1,N)=2.0*BK(WVNO1(L)*A,N)/S*HCCLM(L,M)
80     CONTINUE
100    CONTINUE
!
! --- Setup Matrix BMIJ and EMI
!!
          L=0
          M=0
      EMI(1,N)=DBL/DHKL
      BMIJ(1,1,N)=N/A/DHKL/WK*HCCLM(L,M)/GWL(L)
!
           M=0
      DO 120 L=1, LLOD
        BMIJ(L+1,1,N)=N/A/DBK(WVNO1(L)*A,N)/WVNO1(L)*HCCLM(L,M)/GWL(L)
120    CONTINUE
!
          L=0
      DO 150 M=1, MMOD
        VID=DBIE(WVZJ1(M)*A,WVZJ1(M)*A,N)/BIE(WVZJ1(M)*A,WVZJ1(M)*A,N)
        !VID=WVZJ1(M)*DBIE(WVZJ1(M)*A,WVZJ1(M)*A,N)/BIE(WVZJ1(M)*A,WVZJ1(M)*A,N)
       BMIJ(1,M+1,N)=WVZJ1(M)*VID/DHKL/WK*HCCLM(L,M)/GWL(L)
150   CONTINUE
!
      DO 200 L=1, LLOD
       DO 180 M=1, MMOD
         VID=DBIE(WVZJ1(M)*A,WVZJ1(M)*A,N)/BIE(WVZJ1(M)*A,WVZJ1(M)*A,N)
         !VID=WVZJ1(M)*DBIE(WVZJ1(M)*A,WVZJ1(M)*A,N)/BIE(WVZJ1(M)*A,WVZJ1(M)*A,N)
         BMIJ(L+1,M+1,N)=WVZJ1(M)*VID/DBK(WVNO1(L)*A,N)/WVNO1(L)*HCCLM(L,M)/GWL(L)
180    CONTINUE
200   CONTINUE
! 
500    CONTINUE
       Print *,'After 500'
!
! ----------------------------------
!
! --- Setup matrix CCM
!    [CCM] = [I]- [AMJI] [BMIJ]   
!
        CCM=(0.0D0, 0.0D0)
        DO 600 N=0, NNOD
!
         DO 520 M=1, MMOD+1
          DO 520 MJ=1, MMOD+1
          DO 510 L=1, LLOD+1
           CCM(M,MJ,N)=CCM(M,MJ,N)-AMJI(M,L,N)*BMIJ(L,MJ,N)
510       CONTINUE
520      CONTINUE
!
         DO 540 M=1, MMOD+1
           CCM(M,M,N)=1.0D0+CCM(M,M,N)
540      CONTINUE
!
600    CONTINUE
! ----------------------------------
!
!    [GMJ] = [FMJ]- [AMJI] [EMI]                            
!                                                      


        DO 700 N=0, NNOD
!
         DO 620 M=1, MMOD+1
          GMJ(M,N)=FMJ(M,N)
          DO 610 L=1, LLOD+1
           GMJ(M,N)=GMJ(M,N)-AMJI(M,L,N)*EMI(L,N)
610       CONTINUE
620      CONTINUE
!
700     CONTINUE
!
        DO 800 N=0,  NNOD
         Print *,' N=',N
         CALL SOLVE(N,CCM,GMJ,MLAR,1,MMOD+1,1,NAR)
800     CONTINUE
!
	       BM=GMJ
!
! --- {AM} =[BMiJ] {BM} -{EMi}
!
        DO 900 N=0, NNOD
         DO 820 L=1, LLOD+1
          AM(L,N)=-EMI(L,N)
          DO 810 M=1, MMOD+1
           AM(L,N)=AM(L,N)+BMij(L,M,N)*BM(M,N)
810       CONTINUE
820      CONTINUE
900    CONTINUE
!
! -----------------------------
!
	 DO N=0, NNOD
        WRITE(10,*) N
	  DO M=1, MMOD
        WRITE(10,1010) M,AM(M,N), BM(M,N)
	  END DO
	 END DO

1010	 FORMAT(I4,6E15.5)
1020	 FORMAT(2I4,6E15.5)
     
	 RETURN
     END

!

!C
!C --------------------------------------------------
!C
      SUBROUTINE SOLVE(IP,A,B,NA,NB,N,NRHS,NNOD) 
!C 
!CCC   SOLUTION OF EQUATIONS ( MINIMUM OF PAGING FOR VIRTUAL MACHINES ) 
!C 
!CCC   MATRIX FACTORISATION 
!C 
        COMPLEX*16 A(NA,NA,0:NNOD),B(NA,NB,0:NNOD),ALRLR,ALRLRI,BLRNR,BT,ASUM(1000)   
!C
       N1=N-1 
       IF(N1.NE.0)   GOTO 1 
       IF(CDABS(A(N,N,IP)).LT.1.E-6)  WRITE(6,9999) N,N 
       GOTO 2 
    1  DO 1000 LR=1,N1   
       ALRLR=A(LR,LR,IP) 
       IF(CDABS(ALRLR).LT.1.E-6)  WRITE(6,9999) LR,LR 
       IS=LR+1 
       DO 10 ITBP=IS,N 
   10  A(ITBP,LR,IP)=A(ITBP,LR,IP)/ALRLR 
       DO 100 LRI=IS,N 
       ALRLRI=A(LR,LRI,IP) 
       DO 100 ITB=IS,N 
  100  A(ITB,LRI,IP)=A(ITB,LRI,IP)-A(ITB,LR,IP)*ALRLRI 
!C 
!CCC   RHS FACTORISATION 
!C 
       DO 1000 NR=1,NRHS 
       BLRNR=B(LR,NR,IP) 
       DO 1000 ITB=IS,N 
 1000  B(ITB,NR,IP)=B(ITB,NR,IP)-A(ITB,LR,IP)*BLRNR 
!C 
!CCC   BACK SUBSTITUTION 
!C 
    2  CONTINUE 
       DO 2000 NR=1,NRHS 
       B(N,NR,IP)=B(N,NR,IP)/A(N,N,IP) 
       IF(N1.EQ.0)GOTO 2000 
       BT=B(N,NR,IP) 
       DO 20 ITB=1,N 
   20  ASUM(ITB)=A(ITB,N,IP)*BT 
       DO 200 IBT=N1,1,-1 
       B(IBT,NR,IP)=(B(IBT,NR,IP)-ASUM(IBT))/A(IBT,IBT,IP) 
       BT=B(IBT,NR,IP) 
       DO 200 ITB=1,IBT 
  200  ASUM(ITB)=ASUM(ITB)+A(ITB,IBT,IP)*BT 
 2000  CONTINUE 
 9999  FORMAT(5X,' ** WARNING **',2I3,' TH DIAG. COEF. ZERO IN CROFAC') 
       RETURN 
       END 
 
 
 
