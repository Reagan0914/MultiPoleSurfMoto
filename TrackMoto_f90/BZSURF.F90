!******************************************************************
!
!	SUBROUTINE QB: CALCULATE THE BALL SPHERICAL COORDINATES
!
!******************************************************************
    SUBROUTINE QBFUN(U,V,QB)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: U, V, QB(3)

	QB(1)=COS(U)*COS(V)
	QB(2)=COS(U)*SIN(V)
   	QB(3)=SIN(U)

 

    RETURN

	END SUBROUTINE QBFUN



!*********************************************************************************************************
!
!	SUBROUTINE Q: CALCULATE THE RADIUS VECTOR ON THE SURFACE R=(TMPX,TMPY,TMPZ)=r(cosucosv,cosusinv,sinu)
!
!********************************************************************************************************

	SUBROUTINE QFUN(U,V,Q)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: U, V, Q(3), TMPR, QB(3)

	INTEGER I, J, K, LNUM

	CALL QBFUN(U,V,QB)

	LNUM=0

	Q(1)=0
	Q(2)=0
	Q(3)=0


	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

         TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL

         Q(1)=Q(1)+TMPR*QB(1)
		 Q(2)=Q(2)+TMPR*QB(2)
		 Q(3)=Q(3)+TMPR*QB(3)


		
		 END DO

	   END DO

	END DO

	RETURN

	END SUBROUTINE QFUN

!****************************************************************************************
!
!	SUBROUTINE QU: CALCULATE THE U-DERIVATIVE VECTOR ON THE SURFACE  Ru=(Xu,Yu,Zu)
!
!****************************************************************************************

	SUBROUTINE QUFUN(U,V,QU)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER I,J,K,LNUM

    DOUBLE PRECISION ::  QB(3)

	DOUBLE PRECISION :: U, V, QU(3),TMPX, TMPY, TMPZ,TMPRX,TMPRY,TMPRZ,TMPR,DERX,DERY,DERZ,DERRU

	CALL QBFUN(U,V,QB)

	LNUM=0

	QU(1)=0
	QU(2)=0
	QU(3)=0

    XQU=-1*SIN(U)*COS(V)  !! matlab QU(бн) 
	YQU=-1*SIN(U)*SIN(V)
	ZQU=COS(U)

	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

         DERX=(K-I)*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL   ! First derivative of R respect to X
		 DERY=(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL 
		 DERZ=J*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL
		 DERRU=DERX*XQU+DERY*YQU+DERZ*ZQU  !!First derivative of R respect to U
         
 !!!!!!! TERM 1
         TMPX=DERRU*QB(1)
		 TMPY=DERRU*QB(2)
         TMPZ=DERRU*QB(3)



!!!!!!!! TERM 2
		 TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL

		 TMPRX=TMPR*XQU !!-r*sinU*cosV
         TMPRY=TMPR*YQU
         TMPRZ=TMPR*ZQU

		 QU(1)=QU(1)+TMPX+TMPRX
		 QU(2)=QU(2)+TMPY+TMPRY
		 QU(3)=QU(3)+TMPZ+TMPRZ

		
		 END DO

	   END DO

	END DO  



	RETURN

	END SUBROUTINE QUFUN

!******************************************************************
!
!	SUBROUTINE QV: CALCULATE THE V-DERIVATIVE VECTOR ON THE SURFACE
!
!******************************************************************

    SUBROUTINE QVFUN(U,V,QV)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER I,J,K,LNUM

    DOUBLE PRECISION ::  QB(3)

	DOUBLE PRECISION :: U, V, QV(3), TMPX, TMPY, TMPZ,TMPRX,TMPRY,TMPRZ,TMPR,DERX,DERY,DERZ,DERRV

	CALL QBFUN(U,V,QB)


	XQV=-1*COS(U)*SIN(V)
	YQV=COS(U)*COS(V)
	ZQV=0

	QV(1)=0
	QV(2)=0
    QV(3)=0

	LNUM=0

	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

		 DERX=(K-I)*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL   ! First derivative of R respect to X
		 DERY=(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL 
		 DERZ=J*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL
		 DERRV=DERX*XQV+DERY*YQV+DERZ*ZQV                                !!First derivative of R respect to U

         TMPX=DERRV*QB(1)
		 TMPY=DERRV*QB(2)
         TMPZ=DERRV*QB(3)



         TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL
         TMPRX=TMPR*XQV
         TMPRY=TMPR*YQV
         TMPRZ=TMPR*ZQV

		 QV(1)=QV(1)+TMPX+TMPRX
		 QV(2)=QV(2)+TMPY+TMPRY
		 QV(3)=QV(3)+TMPZ+TMPRZ  !!TMPRZ==0


		 END DO

	   END DO

	END DO  

 

	RETURN

	END SUBROUTINE QVFUN

!******************************************************************
!
!	SUBROUTINE QUU: CALCULATE THE UU-DERIVATIVE VECTOR ON THE SURFACE
!
!******************************************************************

	SUBROUTINE QUUFUN(U,V,QUU)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER I, K, J, LNUM

    DOUBLE PRECISION ::  QB(3)

	DOUBLE PRECISION :: U, V, QUU(3), DERX,DERY,DERZ,DERXX,DERXY,DERXZ,DERYY,DERYZ,DERZZ,DERRU,DERRUU,TMPX2,TMPY2,TMPZ2,TMPRX,TMPRY,TMPRZ,TMPR
	CALL QBFUN(U,V,QB)


	XQU=-1*SIN(U)*COS(V)  !dx/du
	YQU=-1*SIN(U)*SIN(V)
	ZQU=COS(U)
	
	XQV=-1*COS(U)*SIN(V)
	YQV=COS(U)*COS(V)
	ZQV=0

        
    XQUU=-1*COS(U)*COS(V)  !!!d^2x/du^2 
	YQUU=-1*COS(U)*SIN(V)
	ZQUU=-1*SIN(U)



	XQVV=-1*COS(U)*COS(V)
	YQVV=-1*COS(U)*SIN(V)
	ZQVV=0

	XQUV=SIN(U)*SIN(V)
	YQUV=-1*SIN(U)*COS(V)
	ZQUV=0



	QUU(1)=0
	QUU(2)=0
	QUU(3)=0

	LNUM=0

	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

		 DERX=(K-I)*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL  !!dr/dx
		 DERY=(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL  !!dr/dy
		 DERZ=J*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL    !dr/dz

		 DERXX=(K-I)*(K-I-1)*QB(1)**(K-I-2)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL   !!d^2r/dx^2 
		 DERXY=(K-I)*(I-J)*QB(1)**(K-I-1)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL   !!d^2r/dxdy
		 DERXZ=(K-I)*J*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL     !!!!d^2r/dxdz
         DERYY=(I-J)*(I-J-1)*QB(1)**(K-I)*QB(2)**(I-J-2)*QB(3)**J*A(LNUM,1)/SSCAL   !!d^2r/dy^2 
		 DERYZ=(I-J)*J*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**(J-1)*A(LNUM,1)/SSCAL     !!d^2r/dydz 
		 DERZZ=J*(J-1)*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-2)*A(LNUM,1)/SSCAL       !!d^2r/dz^2 
         
         DERRU=DERX*XQU+DERY*YQU+DERZ*ZQU

		 DERRUU=(DERXX*XQU+DERXY*YQU+DERXZ*ZQU)*XQU+DERX*XQUU+(DERXY*XQU+DERYY*YQU+DERYZ*ZQU)*YQU+DERY*YQUU+(DERXZ*XQU+DERYZ*YQU+DERZZ*ZQU)*ZQU+DERZ*ZQUU

	
		 !!!!! ITEM TWO
         TMPX2=DERRU*XQU
         TMPY2=DERRU*YQU
         TMPZ2=DERRU*ZQU      

         !! Item three
         TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL

		 TMPRX=TMPR*XQUU
		 TMPRY=TMPR*YQUU
		 TMPRZ=TMPR*ZQUU


         !!!! ITEM ONE(d^2r/dtheta^2)+TWO+THREE

		 QUU(1)=QUU(1)+DERRUU*QB(1)+2*TMPX2+TMPRX
		 QUU(2)=QUU(2)+DERRUU*QB(2)+2*TMPY2+TMPRY
		 QUU(3)=QUU(3)+DERRUU*QB(3)+2*TMPZ2+TMPRZ


		 END DO

	   END DO

	END DO  

	RETURN

	END SUBROUTINE QUUFUN

!******************************************************************
!
!	SUBROUTINE QVV: CALCULATE THE VV-DERIVATIVE VECTOR ON THE SURFACE
!
!******************************************************************

	SUBROUTINE QVVFUN(U,V,QVV)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER I, K, J, LNUM

    DOUBLE PRECISION ::  QB(3)

	DOUBLE PRECISION :: U, V, QVV(3), DERX, DERY, DERZ,DERXX,DERXY,DERXZ,DERYY,DERYZ,DERZZ,DERRV,DERRVV,TMPX3,TMPY3,TMPZ3,TMPRX,TMPRY,TMPRZ,TMPR

	CALL QBFUN(U,V,QB)


	XQU=-1*SIN(U)*COS(V)
	YQU=-1*SIN(U)*SIN(V)
	ZQU=COS(U)
	
	XQV=-1*COS(U)*SIN(V)
	YQV=COS(U)*COS(V)
	ZQV=0

        
    XQUU=-1*COS(U)*COS(V)
	YQUU=-1*COS(U)*SIN(V)
	ZQUU=-1*SIN(U)



	XQVV=-1*COS(U)*COS(V)
	YQVV=-1*COS(U)*SIN(V)
	ZQVV=0

	XQUV=SIN(U)*SIN(V)
	YQUV=-1*SIN(U)*COS(V)
	ZQUV=0


	QVV(1)=0
	QVV(2)=0
	QVV(3)=0

	LNUM=0

	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

		 DERX=(K-I)*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL          !!dr/dx
		 DERY=(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL          !!dr/dy
		 DERZ=J*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL            !!dr/dz

		 DERXX=(K-I)*(K-I-1)*QB(1)**(K-I-2)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dx^2 
		 DERXY=(K-I)*(I-J)*QB(1)**(K-I-1)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dxdy
		 DERXZ=(K-I)*J*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL   !!d^2r/dxdz
         DERYY=(I-J)*(I-J-1)*QB(1)**(K-I)*QB(2)**(I-J-2)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dy^2 
		 DERYZ=J*(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**(J-1)*A(LNUM,1)/SSCAL   !!d^2r/dydz 
		 DERZZ=J*(J-1)*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-2)*A(LNUM,1)/SSCAL     !!d^2r/dz^2 

         DERRV=DERX*XQV+DERY*YQV+DERZ*ZQV
		 DERRVV=(DERXX*XQV+DERXY*YQV+DERXZ*ZQV)*XQV+DERX*XQVV+(DERXY*XQV+DERYY*YQV+DERYZ*ZQV)*YQV+DERY*YQVV+(DERXZ*XQV+DERYZ*YQV+DERZZ*ZQV)*ZQV+DERZ*ZQVV



      	  
	     !!!! item two
         TMPX3=DERRV*XQV
		 TMPY3=DERRV*YQV
		 TMPZ3=DERRV*ZQV  !!==0
     

                    !!!!!  ITME THREE 
		 TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL

         TMPRX=TMPR*XQVV
		 TMPRY=TMPR*YQVV
         TMPRZ=TMPR*ZQVV  !!==0

         QVV(1)=QVV(1)+DERRVV*QB(1)+2*TMPX3+TMPRX
         QVV(2)=QVV(2)+DERRVV*QB(2)+2*TMPY3+TMPRY
         QVV(3)=QVV(3)+DERRVV*QB(3)

		 END DO

	   END DO

	END DO  

	RETURN

	END SUBROUTINE QVVFUN

!******************************************************************
!
!	SUBROUTINE QUV: CALCULATE THE UV-DERIVATIVE VECTOR ON THE SURFACE
!
!******************************************************************

	SUBROUTINE QUVFUN(U,V,QUV)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER I, K, J, LNUM

    DOUBLE PRECISION ::  QB(3)

	DOUBLE PRECISION :: U, V, QUV(3), DERXX,DERXY,DERXZ,DERYY,DERYZ,DERZZ,DERX,DERY,DERZ,TMPX2,TMPY2,TMPZ2,TMPX3,TMPY3,TMPZ3,TMPRX,TMPRY,TMPRZ,DERRU,DERRV,DERRUV,TMPR

	CALL QBFUN(U,V,QB)


	XQU=-1*SIN(U)*COS(V)
	YQU=-1*SIN(U)*SIN(V)
	ZQU=COS(U)
	
	XQV=-1*COS(U)*SIN(V)
	YQV=COS(U)*COS(V)
	ZQV=0

        
    XQUU=-1*COS(U)*COS(V)
	YQUU=-1*COS(U)*SIN(V)
	ZQUU=-1*SIN(U)



	XQVV=-1*COS(U)*COS(V)
	YQVV=-1*COS(U)*SIN(V)
	ZQVV=0

	XQUV=SIN(U)*SIN(V)
	YQUV=-1*SIN(U)*COS(V)
	ZQUV=0



	QUV(1)=0
	QUV(2)=0
	QUV(3)=0

	LNUM=0

	DO K=0,22
	  
	  DO I=0,K
	   
	     DO J=0,I

		 LNUM=LNUM+1

		 DERX=(K-I)*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL
		 DERY=(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL
		 DERZ=J*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL


		 DERXX=(K-I)*(K-I-1)*QB(1)**(K-I-2)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dx^2 
		 DERXY=(K-I)*(I-J)*QB(1)**(K-I-1)*QB(2)**(I-J-1)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dxdy
		 DERXZ=(K-I)*J*QB(1)**(K-I-1)*QB(2)**(I-J)*QB(3)**(J-1)*A(LNUM,1)/SSCAL   !!d^2r/dxdz
         DERYY=(I-J)*(I-J-1)*QB(1)**(K-I)*QB(2)**(I-J-2)*QB(3)**J*A(LNUM,1)/SSCAL !!d^2r/dy^2 
		 DERYZ=J*(I-J)*QB(1)**(K-I)*QB(2)**(I-J-1)*QB(3)**(J-1)*A(LNUM,1)/SSCAL   !!d^2r/dydz 
		 DERZZ=J*(J-1)*QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**(J-2)*A(LNUM,1)/SSCAL     !!d^2r/dz^2 

		 !!!! ITEM TWO 
		 DERRU=DERX*XQU+DERY*YQU+DERZ*ZQU

         TMPX2=DERRU*XQV
         TMPY2=DERRU*YQV
         TMPZ2=DERRU*ZQV !!==0


         !!!! ITEM THREE
         DERRV=DERX*XQV+DERY*YQV+DERZ*ZQV
        
		 TMPX3=DERRV*XQU
         TMPY3=DERRV*YQU
		 TMPZ3=DERRV*ZQU



		 !!!! ITME FOUR 
		 TMPR=QB(1)**(K-I)*QB(2)**(I-J)*QB(3)**J*A(LNUM,1)/SSCAL

         TMPRX=TMPR*XQUV
         TMPRY=TMPR*YQUV
		 TMPRZ=TMPR*ZQUV !!==0

       
		 DERRUV=(DERXX*XQV+DERXY*YQV+DERXZ*ZQV)*XQU+DERX*XQUV+(DERXY*XQV+DERYY*YQV+DERYZ*ZQV)*YQU+DERY*YQUV+(DERXZ*XQV+DERYZ*YQV+DERZZ*ZQV)*ZQU+DERZ*ZQUV


		 QUV(1)=QUV(1)+DERRUV*QB(1)+TMPX2+TMPX3+TMPRX
		 QUV(2)=QUV(2)+DERRUV*QB(2)+TMPY2+TMPY3+TMPRY
		 QUV(3)=QUV(3)+DERRUV*QB(3)+TMPZ3


		 END DO

	   END DO

	END DO  

	RETURN

	END SUBROUTINE QUVFUN










!******************************************************************
!
!	SUBROUTINE QN: CALCULATE THE UNIT NORMAL VECTOR AT Q
!
!******************************************************************

	SUBROUTINE QNFUN(U,V,QN)

	IMPLICIT NONE

	DOUBLE PRECISION :: U, V, QN(3)
	DOUBLE PRECISION :: QU(3), QV(3)

	CALL QUFUN(U,V,QU)
	CALL QVFUN(U,V,QV)
	CALL CROSS(QU,QV,QN)

	QN=QN/SQRT(DOT_PRODUCT(QN,QN))

	RETURN

	END SUBROUTINE QNFUN



	


