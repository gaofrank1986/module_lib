C
C
C  *************************************************
C  *    The subroutine computes wave number        *
C  *  by wave frequency and water depth.           *
C  *    For infinite water depth, give negative H. *
C  *************************************************

        SUBROUTINE WAVECK(SIGMA,H,WK)
        IMPLICIT  NONE

	  REAL*8 SIGMA,H,WK,B,G,A,Y,C
C
C  H: WATER DEPTH    SIGMA: WAVE FREQUENCY
C  C: WAVE CELERITY  Wk: wave number
C
        DATA G/9.807D0/
C
        IF( SIGMA .LE. 0.0D0 )  THEN
	    Print *,' IN WAVECK:  W=',SIGMA
		STOP
	
        ELSE 

          IF( H .GT. 0.0D0) THEN
           B = G * H
           Y = SIGMA * SIGMA *H/G
  12       A=1.0D0/(1.0D0+Y*(0.66667D0+Y*(0.35550D0+
     1       Y*(0.16084D0+Y*(0.063201D0+Y*(0.02174D0+
     2       Y*(0.00654D0+Y*(0.00171D0+
     2       Y*(0.00039D0+Y*0.00011D0)))))))))
  22       C=SQRT(B/(Y+A))
           WK=SIGMA/C
           ELSE IF( H .LE. 0.0D0) THEN
           WK=SIGMA**2/G
        END IF
        RETURN
C

	 END IF


        RETURN
        END

C
C
C  *******************************************************
C  *                                                     *
C  *    The subroutine computes the real and imaginary   *
C  *  roots of dispersion equation by wave frequency     *
C  *  W and water depth H.                               *
C  *                                                     *
C  *******************************************************
C
C
      SUBROUTINE DISPERS(WAVENO,MNO,W,H)
C
CCC   EVALUATION OF THE ROOTS OF THE FOLLOWING EQUATIONS
CCC   BY NEWTON-RAPHSON METHOD,RESULTS ARE GIVEN IN ARRAY WAVENO
CCC   MNO ARE THE NUMBER OF ROOTS REQUIRED, FIRST ROOT IS FROM EQN. (I)
CCC   THE REST ARE THE FIRST +VE (MNO-1) ROOTS OF (II)
CCC   I) W*W/G = K TANH( KH )
CCC   II) -W*W/G = M TAN( MH )
C

      IMPLICIT NONE
	INTEGER,INTENT(IN)::MNO
	REAL*8, INTENT(IN)::W,H
	REAL*8, INTENT(OUT)::WAVENO(800)
C
	INTEGER I,M,MM,IFLAG
      REAL*8 UPLIM,LOLIM,G,PI
      REAL*8 WWH,FUN,DFUN,TRIAL,EXX,EXXOR,CHECK


        DATA G,PI/9.807d0,3.141592653589793d0/
C
	 IF(MNO .GT. 799) THEN
	   Print *,' MNO=',MNO,' To enlarge WAVENO(*)'
	   STOP
	 END IF

C
      DO 1 I=1,9
    1 WAVENO(I)=0.D0
      WWH=W*W*H
C
CCC   CALCULATION OF WAVE NUMBER (ROOT OF EQN. (I))
C
      TRIAL=WWH/G
      M=0
      IF (TRIAL.GT.1.D1) GO TO 20
      EXX=0.D0
      IFLAG=0
   10 FUN=G*TRIAL - WWH/DTANH(TRIAL)
      DFUN=G + WWH/DSINH(TRIAL)/DSINH(TRIAL)
      TRIAL=TRIAL - FUN/DFUN
      EXXOR=DABS(TRIAL - EXX)
      IF (EXXOR.LE.1.0D-10) GO TO 20
      EXX=TRIAL
      GO TO 10
   20 MM=M + 1
      WAVENO(MM)=TRIAL/H
      CHECK=DABS(W*W/G - WAVENO(MM)*DTANH(TRIAL))
      IF (CHECK.GT.1.0D-5) GO TO 999
      IF (MNO.LE.1) RETURN
C
CCC   CALCULATION OF FIRST +VE (MNO-1) ROOTS OF EQN. (II)
C
      M=1
      IFLAG=1
      EXX=0.D0
      IF (WWH.LE.2.25D2) GO TO 120
      GO TO 110
  100 M=MM
      EXX=0.D0
      IF (MM.EQ.MNO) GO TO 9999
      IF (IFLAG.EQ.2) GO TO 120
  110 TRIAL=(DBLE(FLOAT(M)) - 1.D-1)*PI
      GO TO 140
  120 IFLAG=2
      TRIAL=(DBLE(FLOAT(M)) - 5.D-1)*PI + 1.D-1
  140 IF(IFLAG .EQ. 1)GO TO 160
      IF(IFLAG .EQ. 2)GO TO 170
  150 TRIAL=TRIAL - FUN/DFUN
      EXXOR=DABS(TRIAL - EXX)
      IF (EXXOR.LE.1.D-10) GO TO 180
      EXX=TRIAL
      IF(IFLAG .EQ. 2)GO TO 170
  160 FUN=G*TRIAL + WWH/DTAN(TRIAL)
      DFUN=G - WWH/DSIN(TRIAL)/DSIN(TRIAL)
      GO TO 150
  170 FUN=WWH/(G*TRIAL) + DTAN(TRIAL)
      DFUN=-WWH/(G*TRIAL*TRIAL) + 1.D0/DCOS(TRIAL)/DCOS(TRIAL)
      GO TO 150
  180 UPLIM=(DBLE(FLOAT(M)) + 5.D-1)*PI
      LOLIM=UPLIM - PI
      IF ((TRIAL.GT.UPLIM).OR.(TRIAL.LT.LOLIM)) GO TO 190
      MM=M + 1
      WAVENO(MM)=TRIAL/H
      CHECK=DABS(W*W/G + WAVENO(MM)*DTAN(TRIAL))
      IF (CHECK.GT.1.0D-5) GO TO 999
      GO TO 100
  190 IF (IFLAG.EQ.1) GO TO 120
      WRITE(6,200)M
  200 FORMAT('  **ERROR**  OCCURS AT M',I3)
      STOP
  999 WRITE(6,1000)CHECK
      GO TO 190
 1000 FORMAT('  CHECK=',D11.4)
 9999 CONTINUE

      RETURN
      END





