C 
C ************************************************************************* 
C *                                                                       *  
C *  Analytic solution for the first order wave diffraction               *  
C *  from a truncated cylinder                                            *  
C *                                                                       * 
C *  A: radius              BH: draft             H: water depth          * 
C *  AMP: wave amplitude    WK: wave number                               * 
C *                                                                       * 
C *  LLOD: eigen-mode in outer domain                                     *                            
C *  MMOD: eigen-mode in inner domain                                     * 
C *  NNOD: Fourier mode number                                            * 
C *                                                                       * 
C ************************************************************************* 
C  
        PROGRAM    TRUNCYL 
          USE VAR_mod  
          use io_utils
        IMPLICIT   NONE 
C 
          INTEGER    I,J,II,MN,NKNUM,KOUT_FLAG,m,reason
          REAL*8     TPER,DWK,Wk1,Wk2 ,dbi
          REAL*8     WVN(800),V ,pot
          real(8) ::x,y,z,r,theta
          complex*16 ::re,r2,coeff
        COMPLEX*16 FORC1,FORC3,FORC5 
        
C  
C ----------------------------------------  
C 
        !OPEN(1,  FILE='./DATTC.txt',      STATUS='OLD')  
        open(1,file='config.txt',status='old')
        !open(2,file='nodedic.txt',status='old')
        open(2,file='./potential.0000000.out',status='old')

        !open(3,file='out_dif.txt',status='unknown')
        !open(4,file='out_ind.txt',status='unknown')
        !open(5,file='out_dpdr_dif.txt',status='unknown')
        !open(8,file='out_dpdr_ind.txt',status='unknown')

        print *,"done"
        OPEN(9,  FILE='OUTPUT1',    STATUS='UNKNOWN')   
        OPEN(10, FILE='OUTPUT' ,    STATUS='UNKNOWN')   
C                         
        !OPEN(11, FILE='DATBDP1',    STATUS='UNKNOWN')  
        !OPEN(13, FILE='DATBDP3',    STATUS='UNKNOWN')  
C  
C------------------------------------------- 
C 
        READ(1,*)  H, A, BH, AMP 
        !print *,h,a,bh,amp
        READ(1,*)  Wk1, WK2, DWK 
        READ(1,*)  LLOD, MMOD, NNOD 
         READ(1,*) KOUT_FLAG 
         print *,"read input done"
 
! ------------------------------------------------------ 
 
         IF (KOUT_FLAG .EQ.0 ) THEN 
 
 
        OPEN(21, FILE='OFEX1.DAT',  STATUS='UNKNOWN')   
        OPEN(22, FILE='OFEX3.DAT',  STATUS='UNKNOWN')   
        OPEN(23, FILE='OFEX5.DAT',  STATUS='UNKNOWN')   
        OPEN(24, FILE='Oexf.DAT',   STATUS='UNKNOWN')   
        !OPEN(25, FILE='vector.DAT',   STATUS='UNKNOWN')   
        
        write(21, *) ' WK    W1     Re(FORC1)   Im(FORC1)    ABS(FORC1)' 
        write(22, *) ' WK    W1     Re(FORC3)   Im(FORC3)    ABS(FORC3)' 
        write(23, *) ' WK    W1     Re(FORC5)   Im(FORC5)    ABS(FORC5)' 
        write(24, *) ' WK    W1       Fx        Fz       My' 
 
         ELSE    ! Continue to write  
 
        OPEN(21, FILE='OFEX1.DAT',ACCESS='APPEND',STATUS='UNKNOWN')   
        OPEN(22, FILE='OFEX3.DAT',ACCESS='APPEND',STATUS='UNKNOWN')   
        OPEN(23, FILE='OFEX5.DAT',ACCESS='APPEND',STATUS='UNKNOWN')   
        OPEN(24, FILE='Oexf.DAT', ACCESS='APPEND',STATUS='UNKNOWN')   
 
         ENDIF 
 
 
! ------------------------------------------------------ 
 
 
          NKNUM=(WK2-WK1+0.00001)/DWK 
C 
 
        IF(LLOD .GT. MLAR-1)  THEN 
            LLOD=199 
            Print *,' *** LLOD+1 = MLAR=', MLAR 
          END IF 
 
        IF(MMOD .GT. MLAR)  THEN 
            MMOD=199 
            Print *,' *** MMOD+1 = MLAR', MLAR 
          END IF 
 
        IF(NNOD .GT. NAR) THEN 
            NNOD=20 
            Print *,' *** NNOD = NAR', NAR 
          END IF 
 
C 
         WRITE(6,*) ' H=',H,' A=',A 
         WRITE(6,*) ' BH=',BH,' AMP=',AMP 
         WRITE(6,*) ' LLOD=',LLOD,' MMOD=',MMOD,' NNOD=',NNOD 
         print *,"nnod",nnod
! 
 
! 
! ----------------------------------- 
         DO 1000 II=0,  NKNUM 
          WK=WK1+II*DWK 
          print *,"wk",ii,"=",wk
 
          IF (H .LE. 0.0D0) THEN 
            W1=DSQRT(G*WK) 
          ELSE               
            W1=DSQRT(G*WK*DTANH(WK*H)) 
        END IF 
         
         V1=W1*W1/G 
         TPER=2.D0*PI/W1  
!        
          MN=201 
        CALL DISPERS(WVN,MN,W1,H)  
        WVNO1(0)=WVN(1)  
        DO 2 I=2, MN  
          WVNO1(I-1)=WVN(I)  
2       CONTINUE

! get k0 and ki
! 
!       Print *,' After M0MJ' 
! 
! ------------------------------------ 
!                     
        TPER=2.*PI/W1  
        V=W1*W1/G  
!  
        WRITE(6,1111)  AMP,WK,V,W1,TPER  
        WRITE(9,1111)  AMP,WK,V,W1,TPER  
C 
C --------------------------------------------------- 
C 
        S=H-BH 
          DO J=0,  MMOD 
            WVZJ1(J)=J*PI/S 
          END DO 
C 
C --------------------------------------------------- 
C 
        Print *, 'Before SETUP' 
        CALL SETUP 
C 
        Print *,'Before MAXTL' 
        Print *,'After MAXTL' 
        WRITE(6,*) 'AFTER Assembling' 
C 
C -------------------------------------------------- 
C 
        WRITE(6, 1005) 
C 
C ---  Force in surge 
C 
D       Print *,' Before FORSG' 
        CALL FORSG(FORC1) 
D       Print *,'After FORSG' 
C 
C ---  Force in heave 
C 
        CALL FORHV(FORC3) 
C 
C ---  Pitch moment 
C 
        CALL FORPM(FORC5) 
 
        write(21,201) WK, W1, DBLE(FORC1)/RHO/G,  
     1                       DIMAG(FORC1)/RHO/G, CDABS(FORC1)/RHO/G 
        write(22,201) WK, W1, DBLE(FORC3)/RHO/G,  
     1                       DIMAG(FORC3)/RHO/G, CDABS(FORC3)/RHO/G 
       write(23,201) WK, W1, DBLE(FORC5)/RHO/G,  
     1                       DIMAG(FORC5)/RHO/G, CDABS(FORC5)/RHO/G 

        write(24,201) WK, W1,CDABS(FORC1)/RHO/G,CDABS(FORC3)/RHO/G, 
     1                CDABS(FORC5)/RHO/G 
C 
        coeff=-ci/g*w1
!-ci*w1/g
        print *,"w1",w1,"ci",ci,"g",g
        print *,"coeff",coeff!-ci/g*w1
        open(1001,file=getfilename('pot_inc_result',0))
        open(1002,file=getfilename('pot_dif_result',0))
        open(1003,file=getfilename('pot_inc_result_comp',0))
        open(1004,file=getfilename('dpdz_dif_result',0))
        do !i=1,200
            !read before file end
            read (2,*,iostat=reason) m,x,y,z,pot
            if (reason<0) exit
            r = dsqrt(x**2+y**2)
            theta=datan2(y,x)

            call compute_potential(r,theta,z,re,r2)
!            write(3,202) m,dble(re),dimag(re)
            !write(4,202) m,dble(r2),dimag(r2)
            write(1002,'(i6,2f14.8)') m,dble(re),dimag(re)
            write(1003,'(i6,3f14.8)') m,dble(r2),pot,pot-dble(r2)
            !
            write(1001,'(i6,2f14.8)') m,dble(r2),dimag(r2)
            !write(102,'(i6,f14.8)') m,dble(re*coeff)
            re=0.
            r2=0.
            call compute_dpdz(r,theta,z,re,r2)
            write(401,202) m,dble(re),dimag(re)
            write(1004,202) m,dble(r2),dimag(r2)
        end do
202    format(1i8,2f14.10)

!       do !i=1,200
       !read (202,*,iostat=reason) m,x,y,z
        !if (reason<0) exit
       !r = dsqrt(x**2+y**2)
       !theta=datan2(y,x)
       
       !call compute_dpdr(r,theta,z,re,r2)
        !write(5,202) m,dble(re),dimag(re)
        !write(8,202) m,dble(r2),dimag(r2)
        !print *,"..............."
        !print *,re+r2
        !!print *,r2
        !print *," "
        !end do

!       do !i=1,200
       !read (202,*,iostat=reason) m,x,y,z
        !if (reason<0) exit
       !r = dsqrt(x**2+y**2)
       !theta=datan2(y,x)
       
       !!call compute_potential(r,theta,z,re,r2)
         !!print *,re+r2
        !!write(5,202) m,dble(re+r2),dimag(re+r2)
       !call compute_dp2z(r,theta,z,re)
         !print *,re
         !print *," "
        !write(8,202) m,dble(re),dimag(re)
        !end do

!       do 
       !!print *,"dbi", dbi(2.0d0,0)
       !print *,"..........."
       !read (204,*,iostat=reason) m,x,y,z
        !if (reason<0) exit
       !r = dsqrt(x**2+y**2)
       !theta=datan2(y,x)
        !print *,"r,theta=",r,theta 
       !call compute_dpdr(r,theta,z,re,r2)
         !print *,re+r2
        !write(5,202) m,dble(re+r2),dimag(re+r2)
       !call compute_dp2r(r,theta,z,re)
         !print *,re
         !print *," "
        !write(8,202) m,dble(re),dimag(re)
        !end do







201      FORMAT(F10.5,1x,F10.5, 4(1x,E15.7)) 
 
1000     CONTINUE 
C 
C ************************************** 
C 
1005    FORMAT(//,8X,' *****   EXCITING  FORCES  *****')   
1010    FORMAT(F7.3,2x,F7.3,2x,6E12.5)  
1300    FORMAT(/,1X, '       FORCE(',I1,')=',E14.6,2X,E14.6)  
1400    FORMAT(7X,' MAGNITUDE',E14.6)  
1410    FORMAT(F7.3,2x,6(2x,E12.5))  
C  
1111    FORMAT(//,'  WAVE AMPLITUDE=',F6.2, /, 
     1    '  WAVE NUMBER=',F9.5,'   V=',F9.5, /,  
     2    '  ANGULAR FREQU.=',F9.5,'  WAVE PERIOD=',F7.3,/)       
C  
C  
        END       
  
