!
        ! @func : evaluate potential for given (r,theta,z)
        ! @output: [re] diffraction potential
        ! @output: [r2] incident potential
        SUBROUTINE compute_potential(r,theta,z,re,r2)
	     use VAR_mod
        IMPLICIT   NONE
	    INTEGER N, L
	    REAL*8 BSJ,BJ,DBJ,BY,DBY,BK,ZWL,Zm ,r,theta,z,eps
        COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re,r2
        vector = 0 
        do n=0,NNOD!
           eps=2.0d0
           if (N .eq. 0) then
              eps = 1.0d0
           end if
   	    BSJ =BJ(WK*r,n)
            HKL =DCMPLX(BSJ,         BY(WK*r,n) )

            vector(n)=(AM(1,n)*HKL)*Zm(0,z,wk,h)
            vector1(n)=(BSJ)*Zm(0,z,wk,h)
            DO 40 L=1, LLOD
               vector(n)=vector(n)+AM(L+1,n)*BK(WVNO1(L)*r,n)*Zm(L,z,wvno1(l),h)
        
40          CONTINUE

            vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
        end do!n=0,nnod
        re = dot_product(vector,vector2)!diffraction
        r2 = dot_product(vector1,vector2)!incident
        RETURN
        END
!
        SUBROUTINE compute_dp2z(r,theta,z,re)
	     use VAR_MOD
        IMPLICIT   NONE
	    INTEGER N, L,m
	    REAL*8 BI,dYm ,r,theta,z,eps,vj,lamda
        COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re
        vector = 0 
        do n=0,NNOD!
           eps=2.0d0
           if (N .eq. 0) then
              eps = 1.0d0
           end if
            M=0
            vj=(r/a)**n
            vector(n)=(BM(1,n)*vj)*dYm(0,z)
            
            DO 40 M=1, MMOD
                lamda=pi*m/s
                vj = BI(lamda*r,n)/BI(lamda*a,n)
               vector(n)=vector(n)+BM(M+1,n)*vj*dYm(M,z)
        
40          CONTINUE
            vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
        end do!n=0,nnod
        re = dot_product(vector,vector2)!diffraction
        RETURN
        END       

        SUBROUTINE compute_dpdz(r,theta,z,re,r2)
            use VAR_mod
            IMPLICIT   NONE
            INTEGER N, L
            REAL*8 BSJ,BJ,BY,BK,dZm,r,theta,z,eps
            COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re,r2
            vector = 0 
            do n=0,NNOD!
            
                eps=2.0d0
                if (N .eq. 0) then
                    eps = 1.0d0
                end if
                
                BSJ =BJ(WK*r,n)
                HKL =DCMPLX(BSJ,         BY(WK*r,n) )

                vector(n)=(AM(1,n)*HKL)*dZm(0,z)
                vector1(n)=(BSJ)*dZm(0,z)
                
                DO 40 L=1, LLOD
                    vector(n)=vector(n)+AM(L+1,n)*BK(WVNO1(L)*r,n)*dZm(L,z)
                40          CONTINUE

                vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
            end do!n=0,nnod
            re = dot_product(vector,vector2)!diffraction
            r2 = dot_product(vector1,vector2)!incident
            RETURN
        END!


        
        SUBROUTINE compute_p2(r,theta,z,re)
	     use VAR_MOD
        IMPLICIT   NONE
	    INTEGER N, L,m
	    REAL*8 BI,Ym ,r,theta,z,eps,vj,lamda
        COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re
        vector = 0 
        do n=0,NNOD!
           eps=2.0d0
           if (N .eq. 0) then
              eps = 1.0d0
           end if
               !BSJ =BJ(WK*r,n)
            !HKL =DCMPLX(BSJ,         BY(WK*r,n) )
            M=0
            vj=(r/a)**n
            vector(n)=(BM(1,n)*vj)*Ym(0,z)
            
            DO 40 M=1, MMOD
                lamda=pi*m/s
                vj = BI(lamda*r,n)/BI(lamda*a,n)
               vector(n)=vector(n)+BM(M+1,n)*vj*Ym(M,z)
        
40          CONTINUE
            vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
        end do!n=0,nnod
        re = dot_product(vector,vector2)!diffraction
        RETURN
        END
        
         SUBROUTINE compute_dp2r(r,theta,z,re)
	     use VAR_MOD
        IMPLICIT   NONE
	    INTEGER N, L,m
	    REAL*8 BI,Ym ,r,theta,z,eps,lamda,dvj,dbi
            COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re
        vector = 0 
        do n=0,NNOD!
           eps=2.0d0
           if (N .eq. 0) then
              eps = 1.0d0
      end if
               !BSJ =BJ(WK*r,n)
            !HKL =DCMPLX(BSJ,         BY(WK*r,n) )
            M=0
            dvj=n/a*(r/a)**(n-1)
            vector(n)=(BM(1,n)*dvj)*Ym(0,z)
            
            DO 40 M=1, MMOD
                lamda=pi*m/s
                dvj = lamda*DBI(lamda*r,n)/BI(lamda*a,n)
                vector(n)=vector(n)+BM(M+1,n)*dvj*Ym(M,z)
        
40          CONTINUE
            vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
        end do!n=0,nnod
        re = dot_product(vector,vector2)!diffraction
        RETURN
        END       

















        
        SUBROUTINE compute_dpdr(r,theta,z,re,r2)
	     use VAR_mod
             IMPLICIT   NONE
             INTEGER N, L
             REAL*8 DBJ,DBY,BK,ZWL,Zm ,r,theta,z,eps,DJNR,DYNR,DBK
             COMPLEX*16 vector1(0:NNOD),vector(0:NNOD),vector2(0:NNOD),re,r2,DHNR

             vector = 0 
             vector1= 0 
             
             do n=0,NNOD!
                 eps=2.0d0
                 if (N .eq. 0) then
                     eps = 1.0d0
                 end if

                 DJNR=DBJ(WK*r, n) 
                 DYNR=DBY(WK*r, n) 
                 DHNR=DCMPLX(DJNR,DYNR) 

                 vector(n)=wk*(AM(1,n)*DHNR)*Zm(0,z,wk,h)
                 vector1(n)=wk*(DJNR)*Zm(0,z,wk,h)
                 
                 DO 40 L=1, LLOD
                     vector(n)=vector(n)+wvno1(L)*AM(L+1,n)*dBK(WVNO1(L)*r,n)*Zm(L,z,wvno1(L),h)
                 40          CONTINUE
                
               vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
             end do!n=0,nnod
             re = dot_product(vector,vector2)!diffraction
             r2 = dot_product(vector1,vector2)!incident
             RETURN
        END
        !SUBROUTINE compute_potential(r,theta,z,re)
	     !use VAR_mod
        !IMPLICIT   NONE
	    !INTEGER N, L
	    !REAL*8 BSJ,BJ,DBJ,BY,DBY,BK,ZWL,Zm ,r,theta,z,eps
        !COMPLEX*16 vector(0:NNOD),vector2(0:NNOD),HKL,DHKL,re
        !vector = 0 
    !!    r=dble(r)
        !!print *,"r=",r
        !!print *,"theta=",theta
        !!print *,"z=",z
        !do n=0,NNOD!
           !eps=2.0d0
           !if (N .eq. 0) then
              !eps = 1.0d0
           !end if
               !BSJ =BJ(WK*r,n)
            !HKL =DCMPLX(BSJ,         BY(WK*r,n) )

             !!vector(n)=(BSJ+AM(1,n)*HKL)*Zm(0,z)
             !vector(n)=(AM(1,n)*HKL)*Zm(0,z)
             !!print *,Zm(0,z)
            !!DO 40 L=1, LLOD
               !!vector(n)=vector(n)+AM(L+1,n)*BK(WVNO1(L)*r,n)*Zm(L,z)
        
!!40          CONTINUE
            !vector2(n) = -CI*G*AMP/W1*eps*CI**n*dcos(n*theta)
        !end do!n=0,nnod
        !!print *,vector
        !!print *,vector2
       !re = dot_product(vector,vector2)
        !!print *,re
        !!vector=!todo-2.0D0*PI*CI*G*A*RHO*FORC1
!!
       !!     Print *
        !!Print *,   ' F1=', FORC1/RHO/G
        !!Print *,'       ', CDABS(FORC1)/RHO/G
        !!write(9,*) ' F1=', FORC1/RHO/G
        !!write(9,*) '    ', CDABS(FORC1)/RHO/G
        

        !RETURN
        !END
! *************************************************************
! *  Compute surging force                                    *
! *                                                           *
! *************************************************************
!
        SUBROUTINE FORSG(FORC1)
	     use VAR_mod
        IMPLICIT   NONE
	    INTEGER N, L
	    REAL*8 BSJ,BJ,DBJ,BY,DBY,BK,ZWL 
        COMPLEX*16 FORC1,HKL,DHKL
!
        N=1
   	    BSJ =BJ(WK*A,1)
        HKL =DCMPLX(BSJ,         BY(WK*A,1) )
        DHKL=DCMPLX(DBJ(WK*A,1), DBY(WK*A,1))
   
        FORC1=(BSJ+AM(1,1)*HKL)*ZWL(0)
        print *,FORC1
        DO 40 L=1, LLOD
         FORC1=FORC1+AM(L+1,1)*BK(WVNO1(L)*A,1)*ZWL(L)
40      CONTINUE
        FORC1=-2.0D0*PI*CI*G*A*RHO*FORC1
!
	    Print *
        Print *,   ' F1=', FORC1/RHO/G
        Print *,'       ', CDABS(FORC1)/RHO/G
        write(9,*) ' F1=', FORC1/RHO/G
        write(9,*) '    ', CDABS(FORC1)/RHO/G
!        write(21,201) WK, CDABS(FORC1)/RHO/G
201	 FORMAT(F10.5, 2E13.5)

        RETURN
        END
!
! *************************************************************
! *  Compute heaving force                                    *
! *                                                           *
! *************************************************************
!
    SUBROUTINE FORHV(FORC3)
	use VAR_mod
    IMPLICIT   NONE
!
	INTEGER N, M
	REAL*8 WL,BSJ,BJ,DBJ,BY,DBY,BIE,BKE,ZWL 
    COMPLEX*16 FORC3
! 
    N=0
    Print *, ' '
    FORC3=BM(1,0)*A*A/SQRT(2.0D0)/2.0
    DO 40 M=1, MMOD
     WL=WVZJ1(M)
     FORC3=FORC3+BM(M+1,0)*A*BIE(WL*A,WL*A,1)/WL/BIE(WL*A,WL*A,0)*DCOS(WL*S)
40  CONTINUE
    FORC3=2.0D0*PI*RHO*G*FORC3
!
	    Print *
        Print *,'F3=',FORC3/RHO/G
        Print *,'     ',CDABS(FORC3)/RHO/G

        write(9,*) 'F3=',FORC3/RHO/G
        write(9,*) '     ',CDABS(FORC3)/RHO/G
!        write(22,201) WK, CDABS(FORC3)/RHO/G
201	    FORMAT(F10.5, 2E13.5)
        RETURN
        END
!
! *************************************************************
! *  Compute Pitch Moment                                    *
! *                                                          *
! *************************************************************
!
       SUBROUTINE FORPM(FORC5)
	   use VAR_mod

       IMPLICIT   NONE
!
	   INTEGER L,M,N
	   REAL*8 BSJ,BJ,DBJ,BY,DBY,ZZWL,BK,BIE,WL
       COMPLEX*16 PN,FORC5,FORC51,FORC52,HKL,DHKL
!
! Integration on the side
!
        N=1
   	    BSJ= BJ(WK*A,1)
        HKL=DCMPLX(BSJ, BY(WK*A,1))
        DHKL=DCMPLX(DBJ(WK*A,1), DBY(WK*A, 1))
        FORC51=(BSJ+AM(1,1)*HKL)*ZZWL(0)

        DO 40 L=1, LLOD
         FORC51=FORC51+AM(L+1,1)*BK(WVNO1(L)*A,1)*ZZWL(L)
40      CONTINUE
	    FORC51=FORC51*A
!
!  Integration on the bottom
!
        Print *, ' '
!
        FORC52=BM(1,1)*A*A*A/SQRT(2.0)/4.0D0
        DO 60 M=1, MMOD
	    WL=WVZJ1(M)
         FORC52=FORC52+BM(M+1,1)*BIE(WL*A,WL*A,2)/BIE(WL*A,WL*A,1)*A*A/WL*DCOS(WL*S)
60      CONTINUE
  	   FORC5=-2.0D0*PI*CI*G*RHO*(FORC51+FORC52)

        Print *,'F51=',FORC51
        Print *,'F52=',FORC52

	    Print *
        Print *,'F5=',FORC5/RHO/G
        Print *,'     ',CDABS(FORC5)/RHO/G

        write(9,*) 'F5=',FORC5/RHO/G
        write(9,*) '     ',CDABS(FORC5)/RHO/G
!        write(23,201) WK, CDABS(FORC5)/RHO/G
!
201	   FORMAT(F10.5, 2E13.5)

       RETURN
       END

