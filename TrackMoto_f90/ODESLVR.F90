!****************************************************************************
!
!  SUBROUTINE ODESLVR: 
!
!****************************************************************************

	SUBROUTINE ODESLVR()

	USE GLOBAL

	IMPLICIT NONE

	EXTERNAL SURFODE, ORBODE

	INTEGER :: I, EVT, TNUM ,STOPNUMBER
	INTEGER :: FLG(2), COND(2)
	DOUBLE PRECISION :: T, SX(4), SX1(4), OX(6), OX1(6), NFC, DTC
	DOUBLE PRECISION :: UV(2,2), LAM, LAM2, DUV(2), VN(3), VT(3)
	DOUBLE PRECISION :: TMP1(3), TMP2(3), TMP3(3), TF, NF

!! JUMP out of same position  
    STOPNUMBER=0
!MODE,ID(2),EVENT,GROUND_TRACE(2),TIME,STATS(6)
	OPEN(UNIT=138, FILE='MOTOLOG.TXT')

!INITIALIZATION		
	SELECT CASE (MODE)
	CASE (0)
	  SX=SURFINI
	  UV(1,:)=SX(1:2)
	  CALL NFORCE(SX,NFC)
	CASE (1)
	  OX=ORBINI
	  CALL ORBLOCO(OX(1:3),UV(1,:))
	  CALL ALTITD(OX(1:3),UV(1,:),DTC)
	CASE DEFAULT
	  WRITE(*,*) "ERROR 1!"
	  STOP
	END SELECT

!INITIAL STATS CHECKING

	SELECT CASE (MODE)
	  CASE (0)
	    IF (NFC.LT.0.0D0) THEN
		  WRITE(*,*) "INITIAL MODE (SLIDE) AND VELOCITY DONOT MATCH!"
		  GOTO 532
		END IF
	  CASE (1)
	    IF (DTC.LT.0.0D0) THEN
		  WRITE(*,*) "INITIAL MODE (ORBIT) AND POISITION DONOT MATCH!"
		  GOTO 532
		END IF
	  END SELECT
	EVT=0

!TIME INITIALIZATION
	T=0.0D0
	TNUM=0
   
!RECORD INITIALIZATION
	SELECT CASE (MODE)
	CASE (0)
	  WRITE(138,173) MODE,ID,EVT,UV(1,:),T,SX,0.0D0,0.0D0
	CASE (1)
	  WRITE(138,173) MODE,ID,EVT,UV(1,:),T,OX
	END SELECT

!INTEGRATION
	DO WHILE (T.LT.TEND)
	TNUM=TNUM+1
	!!STEP ONE STEP
	  SELECT CASE (MODE)
	  CASE(0)
	    CALL RK4(4,T,SX,DT,SURFODE,SX1)
		UV(2,:)=SX1(1:2)
	  CASE(1)
	    CALL RK4(6,T,OX,DT,ORBODE,OX1)
	    CALL ORBLOCO(OX1(1:3),UV(2,:))
	  END SELECT
	  T=T+DT
	  EVT=0



  
  IF (((OX(1)-OX1(1)).LT.10D-10).AND.((OX(2)-OX1(2)).LT.10D-10).AND.((OX(3)-OX1(3)).LT.10D-10)) THEN 
  STOPNUMBER=STOPNUMBER+1
  END IF 
  IF (STOPNUMBER.GT.3) THEN
 SIDX=1 
  GOTO 532
 END IF 

!!!!!!!!!!!!!!
	    SELECT CASE (MODE)
		CASE (0)
		  CALL NFORCE(SX1,NFC)
		  IF (NFC.LT.0.0D0) THEN
		    CALL LIFT(UV,SX1,OX1)			
			MODE=1
			EVT=8		
		  END IF	  
		CASE (1)
		  CALL ALTITD(OX1(1:3),UV(2,:),DTC)
		  
		  IF (DTC.LT.0.0D0) THEN


			
		    CALL COLLS(UV,FLG,OX,OX1,VN,VT,LAM2)
!			LAM2 = 0
!			IF (LOVE1.EQ.0)THEN
!			SIDX=0
!			GOTO 532
!			END IF
			CALL CMCHK(UV,FLG,OX1(1:3),VN,VT,COND,DUV)
			IF ((COND(1).EQ.0).AND.(COND(2).EQ.0)) THEN
			  MODE=0
			  SX1(1:2)=UV(2,:)
			  SX1(3:4)=DUV
			  EVT=5
			ELSE
			  OX1(4:6)=VN+VT
			  EVT=4
			END IF
!			T=T-(1.0D0-LAM2)*DT
		  ELSEIF(DTC.GT.10D+01)THEN
		  GOTO 532
		  END IF
		CASE DEFAULT
		  WRITE(*,*) "ERROR 7!"
		END SELECT

	  CALL ORBLOCO(OX1(1:3),UV(1,:))
	  CALL ALTITD(OX1(1:3),UV(1,:),DTC)
	  IF (DTC.LT.0.0D00)THEN
	  WRITE (*,*)(DTC)
	  ENDIF

	  
	  
	  SELECT CASE (MODE)
	  CASE (0)
	    WRITE(138,173) MODE,ID,EVT,UV(2,:),T,SX1,0.0D0,0.0D0
		CALL QUFUN(SX1(1),SX1(2),TMP1)
		CALL QVFUN(SX1(1),SX1(2),TMP2)
		TMP3=TMP1*SX1(3)+TMP2*SX1(4)
		IF (SQRT(DOT_PRODUCT(TMP3,TMP3)).LT.SFV) THEN
		  CALL SFSCHK(SX1,TMP1,TMP2,TF,NF)
		  IF (TF.LT.(MUS*NF)) THEN
		    !WRITE(*,*) "------STATIC FRICTION STUCK!"
			SIDX=1
		    GOTO 532
		  END IF
		END IF
	  CASE (1)
	    WRITE(138,173) MODE,ID,EVT,UV(2,:),T,OX1
	  END SELECT

	  UV(1,:)=UV(2,:)
	  FLG(1)=FLG(2)
	  SX=SX1
	  OX=OX1

      

	  IF(TNUM.GT.2.5D+8)THEN
	  SIDX=0
	  GOTO 532
	  END IF

	  SIDX=1

	END DO

532	CLOSE(138)

173 FORMAT(I12,' ',I12,' ',I12,' ',I12,' ',E22.12,' ',E22.12,' ',E22.12,' ',&
           E22.12,' ',E22.12,' ',E22.12,' ',E22.12,' ',E22.12,' ',E22.12) 

	RETURN

	END SUBROUTINE ODESLVR

!*****************************************************************************
!
! RK4 TAKES ONE RUNGE-KUTTA STEP.
!
!  DISCUSSION:
!
!    IT IS ASSUMED THAT AN INITIAL VALUE PROBLEM, OF THE FORM
!
!      DU/DT = F ( T, U )
!      U(T0) = U0
!*****************************************************************************

	SUBROUTINE RK4(DIM,T0,U0,DT,F,U1)

	IMPLICIT NONE

	EXTERNAL F

	INTEGER :: DIM
	DOUBLE PRECISION :: T0, U0(DIM), DT, U1(DIM)
	DOUBLE PRECISION :: F1(DIM), F2(DIM), F3(DIM), F4(DIM)

	CALL F(T0,            U0,               F1)
	CALL F(T0+DT/2.0D+00, U0+DT*F1/2.0D+00, F2)
	CALL F(T0+DT/2.0D+00, U0+DT*F2/2.0D+00, F3)
	CALL F(T0+DT,         U0+DT*F3,         F4)

	U1=U0+DT*(F1+2.0D0*F2+2.0D0*F3+F4)/6.0D0

	RETURN

	END SUBROUTINE RK4

!*****************************************************************************
!
! GETDT: ESTIMATE RUNGE-KUTTA TIME STEP.
!
!*****************************************************************************


	SUBROUTINE GETDT()

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: ACC

	ACC=ABS(4.0D0*PI/3.0D0*KAPPA-4.0D0*PI**2.0D0)
	DT=SQRT(2.0D0*RESO_H/ACC)/RESO_N

	RETURN

	END SUBROUTINE GETDT