C 
C  BKE =Kn(X)*EXP(X1)         BIE =In(X)*EXP(-X1)
C  DBKE=Kn'(X)*EXP(X1)        DBIE=In'(X)*EXP(-X1)
C
C 
      DOUBLE PRECISION FUNCTION BKE(X,X1,N) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      CALL BESELKM(X,X1,N,BK0,BK1) 
      BKE=BK0 
      RETURN 
      END 
C 
      DOUBLE PRECISION FUNCTION DBKE(X,X1,N) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      N1=N+1 
      CALL BESELKM(X,X1,N1,BK1,BK0) 
      DBKE=N/X*BK0-BK1 
      RETURN 
      END 
C 
      DOUBLE PRECISION FUNCTION BIE(X,X1,N) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION B(100) 
      NB=N+1 
      IZE=1 
      CALL BESLJI(X,NB,IZE,B,NCALC) 
      BIE=B(NB)*DEXP(-X1) 
      RETURN 
      END 
C 
      DOUBLE PRECISION FUNCTION DBIE(X,X1,N) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION B(100) 
      NB=N+1 
      NB1=NB+1 
      IZE=1 
      CALL BESLJI(X,NB1,IZE,B,NCALC) 
      DBI=N/X*B(NB)+B(NB1) 
      DBIE=DBI*DEXP(-X1)
      RETURN 
      END 
C 
C 
C 
      SUBROUTINE BESELKM(X,X1,INN,BKN1,BKN) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION BK(200) 
      IFAIL=0 
      BK(1)=S18ACFM(X,X1,IFAIL) 
      BK(2)=S18ADFM(X,X1,IFAIL) 
   10 IN=IABS(INN) 
      IF(INN.LT.1)IN=IN+1 
      I=2 
      IF(IN.EQ.1)GO TO 2 
      IN1=IN-1 
      DO 1 J=1,IN1 
      I=I+1 
    1 BK(I)=2.D0*(I-2)/X*BK(I-1) + BK(I-2) 
    2 BKN1=BK(I) 
      BKN =BK(I-1) 
      IF(INN.GE.1)RETURN 
      TEMP=BKN 
      BKN =BKN1 
      BKN1=TEMP 
      RETURN 
      END 
C 
CC 
C JZK0 
      DOUBLE PRECISION FUNCTION S18ACFM(X,X1,IFAIL) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DOUBLE PRECISION X 
      INTEGER IFAIL 
      DOUBLE PRECISION SRNAMEM
      DOUBLE PRECISION EGAM,G,T,XBIG,XVSMAL,Y 
      DOUBLE PRECISION DLOG,DEXP,DSQRT 
      INTEGER P01AAF 
      DATA SRNAMEM/8H S18ACFM / 
      DATA XVSMAL,EGAM/3.2D-9,5.772156649015329D-1/ 
      DATA XBIG/10000.1D0/ 
      IF(X.LE.0.0D0)GOTO120 
      IFAIL=0 
      IF(X.LT.XBIG)GOTO20 
      S18ACFM=0.0D0 
      GOTO140 
20    IF(X.GT.4.0D0)GOTO100 
      IF(X.GT.2.0D0)GOTO80 
      IF(X.GT.1.0D0)GOTO60 
      IF(X.GT.XVSMAL)GOTO40 
      S18ACFM=-(DLOG(0.5D0*X)+EGAM) 
      GOTO140 
40    T=2.0D0*X*X-1.0D0 
      G=+1.128960929454128D+0+T*(+1.329769664783382D-1+T*(+4.07157485171 
     *3890D-3+T*(+5.597023382279154D-5+T*(+4.345626715461582D-7+T*(+2.16 
     *3824118247215D-9+T*(+7.491107368941348D-12+T*(+1.906741975145613D- 
     *14))))))) 
      Y=+2.618418792586871D-1+T*(+1.524369217993952D-1+T*(+6.63513979313 
     *9438D-3+T*(+1.095342926324015D-4+T*(+9.578784932659294D-7+T*(+5.19 
     *9068658006656D-9+T*(+1.924052642197067D-11+T*(+5.168678869463322D- 
     *14+T*(+1.054077181913600D-16)))))))) 
      S18ACFM=-DLOG(X)*G+Y 
	S18ACFM=S18ACFM*DEXP(X1)
      GOTO140 
60    T=2.0D0*X-3.0D0 
      Y=+9.582100532948965D-1+T*(-1.424779101288283D-1+T*(+3.23582010649 
     *6530D-2+T*(-8.277803503516927D-3+T*(+2.247097296177705D-3+T*(-6.32 
     *6783574605949D-4+T*(+1.826524600893428D-4+T*(-5.371012088984418D-5 
     *+T*(+1.601859741497206D-5+T*(-4.831342503369222D-6+T*(+1.470557960 
     *782317D-6+T*(-4.510172923752000D-7+T*(+1.392172702246142D-7+T*(-4. 
     *321850898418341D-8+T*(+1.347904673613401D-8+T*(-4.205973292582499D 
     *-9+T*(+1.320693623859689D-9+T*(-4.333266656187809D-10+T*(+1.379992 
     *680744427D-10+T*(-3.192410591988521D-11+T*(+9.744101522706792D-12+ 
     *T*(-7.837386091085693D-12+T*(+2.574662885758206D-12))))))))))))))) 
     *))))))) 
      S18ACFM=DEXP(-X+X1)*Y 
      GOTO140 
80    T=X-3.0D0 
      Y=+6.977615980438518D-1+T*(-1.088018820849351D-1+T*(+2.56253646031 
     *9603D-2+T*(-6.744596079401692D-3+T*(+1.872929397259624D-3+T*(-5.37 
     *1456229719100D-4+T*(+1.574515162358606D-4+T*(-4.689366538148967D-5 
     *+T*(+1.413765093436227D-5+T*(-4.303738717272685D-6+T*(+1.320522610 
     *589324D-6+T*(-4.078512078621890D-7+T*(+1.266726294175674D-7+T*(-3. 
     *954032557135184D-8+T*(+1.239231378983469D-8+T*(-3.883497052505557D 
     *-9+T*(+1.224249827794330D-9+T*(-4.034246078719601D-10+T*(+1.289055 
     *874799801D-10+T*(-2.977875646332351D-11+T*(+9.111094308330013D-12+ 
     *T*(-7.396727839879332D-12+T*(+2.435382422475375D-12))))))))))))))) 
     *))))))) 
      S18ACFM=DEXP(-X+X1)*Y 
      GOTO140 
100   T=10.0D0/(1.0D0+X)-1.0D0 
      Y=+1.236886647694254D+0+T*(-1.726836523853216D-2+T*(-9.25551464765 
     *6371D-4+T*(-9.025533451874046D-5+T*(-6.316923983337465D-6+T*(-7.69 
     *1776225292729D-7+T*(-4.160448111741146D-8+T*(-9.415553211371761D-9 
     *+T*(+1.753593212735806D-10+T*(-2.228295822888333D-10+T*(+3.4956429 
     *32565460D-11+T*(-1.113917585726476D-11+T*(+2.854812351677059D-12+T 
     **(-7.313444826639319D-13+T*(+2.063288925625549D-13+T*(-1.281083108 
     *269916D-13+T*(+4.437419798865510D-14)))))))))))))))) 
      S18ACFM=DEXP(-X+X1)*Y/DSQRT(X) 
      GOTO140 
120   IFAIL=P01AAF(IFAIL,1,SRNAMEM) 
      S18ACFM=0.0D0 
140   RETURN 
      END 
C 
C JZK1 
      DOUBLE PRECISION FUNCTION S18ADFM(X,X1,IFAIL) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DOUBLE PRECISION X 
      INTEGER IFAIL 
      DOUBLE PRECISION SRNAMEM
      DOUBLE PRECISION G,T,XBIG,XSEST,XSMALL,Y 
      DOUBLE PRECISION DLOG,DEXP,DSQRT 
      INTEGER P01AAF 
      DATA SRNAMEM/8H S18ADFM / 
      DATA XSMALL/7.9D-10/ 
      DATA XBIG,XSEST/10000.1D0,1.4D-36/ 
      IF(X.LE.0.0D0)GOTO120 
      IF(X.LE.XSEST)GOTO140 
      IFAIL=0 
      IF(X.LT.XBIG)GOTO20 
      S18ADFM=0.0D0 
      GOTO160 
20    IF(X.GT.4.0D0)GOTO100 
      IF(X.GT.2.0D0)GOTO80 
      IF(X.GT.1.0D0)GOTO60 
      IF(X.GT.XSMALL)GOTO40 
      S18ADFM=1.0D0/X 
      GOTO160 
40    T=2.0D0*X*X-1.0D0 
      G=+5.319078659133528D-1+T*(+3.257259881371105D-2+T*(+6.71642805873 
     *4987D-4+T*(+6.953002745482062D-6+T*(+4.327648236429978D-8+T*(+1.79 
     *7847923801558D-10+T*(+5.338882686656589D-13+T*(+1.189649624399104D 
     *-15))))))) 
     * 
      Y=+3.518258282893255D-1+T*(+4.504904429669437D-2+T*(+1.20333585658 
     *2190D-3+T*(+1.446124325330061D-5+T*(+9.966866892737815D-8+T*(+4.46 
     *8286284356187D-10+T*(+1.409171030245143D-12+T*(+3.298810580198656D 
     *-15))))))) 
      S18ADFM=1.0D0/X+X*(DLOG(X)*G-Y) 
	S18ADFM=S18ADFM*DEXP(X1)
      GOTO160 
60    T=2.0D0*X-3.0D0 
     * 
      Y=+1.243165873552553D+0+T*(-2.719107143886894D-1+T*(+8.20250220860 
     *6939D-2+T*(-2.625458187294274D-2+T*(+8.573880870674101D-3+T*(-2.82 
     *4507878416560D-3+T*(+9.345941543876429D-4+T*(-3.100076810136266D-4 
     *+T*(+1.029827467000607D-4+T*(-3.424249122119421D-5+T*(+1.139301692 
     *025535D-5+T*(-3.792276988211429D-6+T*(+1.262655783319419D-6+T*(-4. 
     *205071523389350D-7+T*(+1.401383519851855D-7+T*(-4.669289121680201D 
     *-8+T*(+1.544566539090127D-8+T*(-5.137835081403322D-9+T*(+1.8280838 
     *13812054D-9+T*(-6.152114168988951D-10+T*(+1.280440239499463D-10+T* 
     *(-4.025910666270238D-11+T*(+4.274043305687672D-11+T*(-1.4663929178 
     *29485D-11))))))))))))))))))))))) 
      S18ADFM=DEXP(-X+X1)*Y 
      GOTO160 
80    T=X-3.0D0 
      Y=+8.065634801287869D-1+T*(-1.600526112913272D-1+T*(+4.58591528414 
     *0231D-2+T*(-1.423631366844236D-2+T*(+4.558657512067247D-3+T*(-1.48 
     *1854720326885D-3+T*(+4.857071747786637D-4+T*(-1.599948736215991D-4 
     *+T*(+5.287129191231318D-5+T*(-1.750895943540799D-5+T*(+5.806923118 
     *422967D-6+T*(-1.927945869964326D-6+T*(+6.405818140373983D-7+T*(-2. 
     *129692293463103D-7+T*(+7.087233666965699D-8+T*(-2.358556184610253D 
     *-8+T*(+7.794216511448327D-9+T*(-2.590393993080091D-9+T*(+9.2078168 
     *59061105D-10+T*(-3.096673923432451D-10+T*(+6.449134235458942D-11+T 
     **(-2.026804015147359D-11+T*(+2.147367510651332D-11+T*(-7.364782970 
     *504217D-12))))))))))))))))))))))) 
      S18ADFM=DEXP(-X+X1)*Y 
      GOTO160 
100   T=10.0D0/(1.0D0+X)-1.0D0 
      Y=+1.303875736042304D+0+T*(+5.448452543189316D-2+T*(+4.31639434283 
     *4454D-3+T*(+4.299739708987668D-4+T*(+4.047206315284950D-5+T*(+4.32 
     *7764097842352D-6+T*(+4.075638569318435D-7+T*(+4.866514200081540D-8 
     *+T*(+3.827176921214383D-9+T*(+6.776889438575889D-10+T*(+6.97075379 
     *1177314D-12+T*(+1.720260972859309D-11+T*(-2.607745020202711D-12+T* 
     *(+8.582115237135606D-13+T*(-2.192871044418028D-13+T*(+1.3932112294 
     *06003D-13+T*(-4.778502381115802D-14)))))))))))))))) 
      S18ADFM=DEXP(-X+X1)*Y/DSQRT(X) 
      GOTO160 
120   IFAIL=P01AAF(IFAIL,1,SRNAMEM) 
      S18ADFM=0.0D0 
      GOTO160 
140   IFAIL=P01AAF(IFAIL,2,SRNAMEM) 
      S18ADFM=1.0D0/XSEST*DEXP(X1) 
160   RETURN 
      END 
C 
C 
C
C
C *******************************************************
C *                                                     *
C *  BSELJ0 : Bessel function and its first derivative  *
C *  BSELJ1 : Bessel and Hankel functions and their     *
C *           first derivatives                         *
C *  BSELJ2 : Bessel and Hankel functions and their     *
C *           first and second derivatives              *
C *                                                     *
C *******************************************************
C
	SUBROUTINE  BSELJ0(R,U,N,BJ,DBJ)
      	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      	DIMENSION B(100) 
      	X=R*U
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF
C	  
   	NB1=N+2
   	NB0=N+1 
   	IZE=0 
  	CALL BESLJI(X,NB1,IZE,B,NCALC) 
C
	BJ=B(NB0)
	BP1=B(NB1)
	
	IF (N .EQ. 0) THEN
	 BN1=-BP1
 	ELSE 
	 BN1=B(NB0-1)
	END IF
C	
	DBJ=(BN1-BP1)/2.0D0
C
   	RETURN 
	END 
C
C  ------------------------
C
	SUBROUTINE  BSELJ1(R,U,N,BJ,DBJ,BH,DBH)
      	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        COMPLEX*16  BH,DBH
      	DIMENSION B(100) 
      	X=R*U
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF
C	  
   	NB1=N+2
   	NB0=N+1 
   	IZE=0 
  	CALL BESLJI(X,NB1,IZE,B,NCALC) 
C
	BJ=B(NB0)
	BP1=B(NB1)
	
	IF (N .EQ. 0) THEN
	 BN1=-BP1
 	ELSE 
	 BN1=B(NB0-1)
	END IF
C	
	DBJ=(BN1-BP1)/2.0D0
C
	CALL BESELY(X,NB0,BY1,BY0) 
	BY=BY0 
	DBY=N/X*BY0-BY1 

	BH  =DCMPLX(BJ,   BY)
	DBH =DCMPLX(DBJ,  DBY)
C
   	RETURN 
	END 
	
C
C  ------------------------
C
	SUBROUTINE  BSELJ2(R,U,N,BJ,DBJ,DDBJ,BH,DBH,DDBH)
      	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        COMPLEX*16  BH,DBH,DDBH
      	DIMENSION B(100) 
      	X=R*U
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF 
C
   	NB2=N+3
   	NB1=N+2
   	NB0=N+1 
   	IZE=0 
  	CALL BESLJI(X,NB2,IZE,B,NCALC) 
C
	BJ=B(NB0)
	BP1=B(NB1)
	BP2=B(NB2)
	
	IF (N .EQ. 0) THEN
	 BN1=-BP1
	 BN2= BP2
 	ELSE IF(N .EQ. 1) THEN
	 BN1=B(NB0-1)
	 BN2=-BJ
 	ELSE 
	 BN1=B(NB0-1)
	 BN2=B(NB0-2)
	END IF
C	
	DBJ=(BN1-BP1)/2.0D0
	DDBJ=(BN2+BP2-2.0D0*BJ)/4.0D0
C
	CALL BESELY(X,NB0,BY1,BY0) 
	BY=BY0 
	DBY=N/X*BY0-BY1 
	DDBY=-BY*(1.0D0+(1.0D0-DFLOAT(N))*DFLOAT(N)/X/X)+BY1/X

	BH  =DCMPLX(BJ,   BY)
	DBH =DCMPLX(DBJ,  DBY)
	DDBH=DCMPLX(DDBJ, DDBY)
C
   	RETURN 
	END 
	
C
C
C *******************************************************
C *                                                     *
C *  BSELK1 : Modified Bessel function K and its        *
C *           first derivatives                         *
C *  BSELK2 : Modified Bessel function K and its        *
C *           first and second derivatives              *
C *                                                     *
C *******************************************************
C
C
	SUBROUTINE  BSELK1(R,U,N,BK,DBK)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF
C
	X=R*U
	N1=N+1 
	CALL BESELK(X,N1,BK1,BK) 
	DBK=N/X*BK-BK1 
	RETURN 
	END 
	
	

C  ------------------------
C
	SUBROUTINE  BSELK2(R,U,N,BK,DBK,DDBK)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF
C
	X=R*U
	N1=N+1 
	CALL BESELK(X,N1,BK1,BK) 
	DBK=N/X*BK-BK1 
	DDBK=BK*(1.0D0-(1.0D0-DFLOAT(N))*DFLOAT(N)/X/X)+BK1/X

	RETURN 
	END 

C  ------------------------
C
	SUBROUTINE  BSELI0(R,U,N,BSI,DBSI)
      	IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      	DIMENSION B(100) 
      	X=R*U
C
	IF(N .LT. 0) THEN
	  PRINT *,'  J(N<0)'
	  STOP
	END IF
C	  
   	NB1=N+2
   	NB0=N+1 
   	IZE=1
  	CALL BESLJI(X,NB1,IZE,B,NCALC) 
C
	BSI=B(NB0)
	BP1=B(NB1)
	
	IF (N .EQ. 0) THEN
	 BN1=BP1
 	ELSE 
	 BN1=B(NB0-1)
	END IF
C	
	DBSI=(BN1+BP1)/2.0D0
C
   	RETURN 
	END 
