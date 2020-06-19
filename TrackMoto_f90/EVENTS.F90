!****************************************************************************
!
!  SUBROUTINE ORBLOCO: FIND LOCAL COORDINATES OF ORBIT
!
!****************************************************************************

	SUBROUTINE ORBLOCO(X,Y)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: X(3), Y(2), RXYZ
	DOUBLE PRECISION :: P1(3), P12(3), P13(3)
	DOUBLE PRECISION :: TMP1(3), TMP2(3), TMP3(3)

	RXYZ=SQRT(X(1)**2+X(2)**2)
    Y(1)=ATAN2(X(3),RXYZ)    !!
    Y(2)=ATAN2(X(2),X(1))


    RETURN
	
	END SUBROUTINE ORBLOCO
	
	
!****************************************************************************
!
!	SUBROUTINE : CHECK THE CONTACT FORCE OF THE SURFACE ON PARTICLE
!
!****************************************************************************	

	SUBROUTINE NFORCE(X,N)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: X(4), N
	DOUBLE PRECISION :: Q(3), QU(3), QV(3), QUU(3), QUV(3), QVV(3), QN(3)
	DOUBLE PRECISION :: GF(3), TMP1(3), TMP2(3), TMP3(3), TMP4(3)

	CALL QFUN(X(1),X(2),Q)
	CALL QUFUN(X(1),X(2),QU)
	CALL QVFUN(X(1),X(2),QV)
	CALL QUUFUN(X(1),X(2),QUU)
	CALL QUVFUN(X(1),X(2),QUV)
	CALL QVVFUN(X(1),X(2),QVV)

	CALL CALGF(Q,GF)
	CALL CROSS(OMG,Q,TMP1)
	CALL CROSS(OMG,TMP1,TMP2)

	TMP1=X(3)*QU+X(4)*QV
	CALL CROSS(OMG,TMP1,TMP3)
	TMP4=QUU*X(3)**2.0D0+2.0D0*QUV*X(3)*X(4)+QVV*X(4)**2.0D0
	TMP1=KAPPA*GF-TMP2-2.0D0*TMP3-TMP4

	CALL CROSS(QU,QV,QN)
	QN=QN/SQRT(DOT_PRODUCT(QN,QN))
	
	N=-DOT_PRODUCT(QN,TMP1)

	RETURN 

	END SUBROUTINE NFORCE

!****************************************************************************
!
!	SUBROUTINE : CHECK THE SIGHT ALTITUDE OF ORBIT
!
!****************************************************************************

	SUBROUTINE ALTITD(X,Y,AL)

	IMPLICIT NONE

	DOUBLE PRECISION :: X(3), Y(2), AL
	DOUBLE PRECISION :: Q(3)

	CALL QFUN(Y(1),Y(2),Q)
	AL=SQRT(DOT_PRODUCT(X,X))-SQRT(DOT_PRODUCT(Q,Q))

	RETURN

	END SUBROUTINE ALTITD
	
!****************************************************************************
!
!	SUBROUTINE CRPH: CROSSING PATCH
!
!****************************************************************************

	SUBROUTINE CRPH(UV,FLG,LAM)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: FLG(2)
	INTEGER :: I, TMP(3), TMP2
	DOUBLE PRECISION :: UV(2,2)
	DOUBLE PRECISION :: U1, V1, U2, V2, DU, DV
	DOUBLE PRECISION :: LAM

	U1=UV(1,1)
	V1=UV(1,2)
	U2=UV(2,1)
	V2=UV(2,2)
	DU=U2-U1
	DV=V2-V1

	DO I=1,3
	  TMP(I)=0
	END DO

	IF ((U1*U2.LT.0.0D0).AND.((-DV*U1+DU*V1)*(-DV*U1-DU*(1.0D0-V1)).LT.0.0D0)) THEN
		TMP(1)=1
	END IF
	IF ((V1*V2.LT.0.0D0).AND.((-DV*U1+DU*V1)*(DV*(1.0D0-U1)+DU*V1).LT.0.0D0)) THEN
	    TMP(2)=1
	END IF
	IF (((U1+V1-1.0D0)*(U2+V2-1.0D0).LT.0.0D0).AND.((DV*(1.0D0-U1)+DU*V1)*(-DV*U1-DU*(1.0D0-V1)).LT.0.0D0)) THEN
		TMP(3)=1
	END IF

	TMP2=SUM(TMP)

	SELECT CASE (FLG(1))
	CASE (0)
      SELECT CASE (TMP2)
	  CASE (0)
		WRITE(*,*) "WARNING 1: CROSSING VERTICES!"
		STOP
	  CASE (1)
	    IF (TMP(1).EQ.1) THEN
		  LAM=-U1/DU
	      FLG(2)=1
		ELSE IF (TMP(2).EQ.1) THEN
		  LAM=-V1/DV
		  FLG(2)=2
		ELSE IF (TMP(3).EQ.1) THEN
		  LAM=(1.0D0-U1-V1)/(DU+DV)
   		  FLG(2)=3
		END IF
	  CASE DEFAULT
	    WRITE(*,*) "ERROR 2!"
		STOP
	  END SELECT
	CASE (1)
	  IF (TMP(2).EQ.1) THEN
	    LAM=-V1/DV
		FLG(2)=2
	  ELSE IF (TMP(3).EQ.1) THEN
		LAM=(1.0D0-U1-V1)/(DU+DV)
		FLG(2)=3
	  ELSE
	    LAM=0.0D0
		FLG(2)=FLG(1)
	  END IF
	CASE (2)
	  IF (TMP(1).EQ.1) THEN
		LAM=-U1/DU
	    FLG(2)=1
	  ELSE IF (TMP(3).EQ.1) THEN
	    LAM=(1.0D0-U1-V1)/(DU+DV)
   		FLG(2)=3
	  ELSE
	    LAM=0.0D0
		FLG(2)=FLG(1)
	  END IF
	CASE (3)
	  IF (TMP(1).EQ.1) THEN
		LAM=-U1/DU
		FLG(2)=1
	  ELSE IF (TMP(2).EQ.1) THEN
	    LAM=-V1/DV
   		FLG(2)=2
	  ELSE
	    LAM=0.0D0
		FLG(2)=FLG(1)
	  END IF
	CASE DEFAULT
      WRITE(*,*) "WARNING 2: CROSSING VERTICES!"
	  STOP
	END SELECT

	UV(2,1)=U1+DU*LAM
	UV(2,2)=V1+DV*LAM

	RETURN

	END SUBROUTINE CRPH

		

!****************************************************************************
!
!	SUBROUTINE SX2OX:
!
!****************************************************************************

	SUBROUTINE SX2OX(SX,OX)

	IMPLICIT NONE

	DOUBLE PRECISION :: SX(4), OX(6)
	DOUBLE PRECISION :: TMP1(3), TMP2(3)

	CALL QFUN(SX(1),SX(2),OX(1:3))
	CALL QUFUN(SX(1),SX(2),TMP1)
	CALL QVFUN(SX(1),SX(2),TMP2)
	OX(4:6)=TMP1*SX(3)+TMP2*SX(4)

	RETURN

	END SUBROUTINE SX2OX


!****************************************************************************
!
!	SUBROUTINE LIFT:
!
!****************************************************************************

	SUBROUTINE LIFT(UV,SX1,OX1)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: UV(2,2),SX1(4),OX1(6)
	DOUBLE PRECISION :: TMP

	CALL SX2OX(SX1,OX1)
	
	DO WHILE (1)
	  CALL ORBLOCO(OX1(1:3),UV(2,:))
	  CALL ALTITD(OX1(1:3),UV(2,:),TMP)
	  IF (TMP.LT.0.0D0) THEN
	    OX1(1:3)=OX1(1:3)*LIF
	  ELSE
	    GOTO 762
	  ENDIF
	END DO
		
762	RETURN

	END SUBROUTINE LIFT

!****************************************************************************
!
!	SUBROUTINE COLLS:
!
!****************************************************************************

	SUBROUTINE COLLS(UV,FLG,OX,OX1,VN,VT,LAM)

	USE GLOBAL
	
	IMPLICIT NONE

	INTEGER :: FLG(2), COLLNUM
	DOUBLE PRECISION :: UV(2,2),OX(6),OX1(6), VN(3), VT(3), LAM
	DOUBLE PRECISION :: NVR(3), DTC, KISS


 !   COLLNUM=0
!	LOVE1=1
!	LAM=1.0D0
!	DO WHILE (1)
!	  COLLNUM=COLLNUM+1
!	  LAM=LAM*CUT
!	  OX1=CUT*OX1+(1.0D0-CUT)*OX
!	  KISS=SQRT((OX(1)-OX1(1))**2+(OX(2)-OX1(2))**2+(OX(3)-OX1(3))**2)
!	  IF (KISS.LT.1.0D-03) THEN
	  OX1 = OX
!	  END IF


!	  CALL ORBLOCO(OX1(1:3),UV(2,:))
!	  CALL ALTITD(OX1(1:3),UV(2,:),DTC)
!	  IF (COLLNUM.GT.5D+3)THEN
!	  	  STOP
!	  LOVE1=0
!	  GOTO 988
!	  END IF
!	  IF (DTC.GE.0.0D0) THEN
!	  WRITE(*,*)(OX(1:3),DTC)
!		GOTO 982
!	  END IF
!	END DO
	
982	CALL QNFUN(UV(2,1),UV(2,2),NVR)
	VN=DOT_PRODUCT(NVR,OX1(4:6))*NVR
	VT=OX1(4:6)-VN

	VN=-0.8*VN
	VT=0.8*VT

988	RETURN

	END SUBROUTINE COLLS

!****************************************************************************
!
!	SUBROUTINE CMCHK: COLLISION MODE CHECK
!
!****************************************************************************

	SUBROUTINE CMCHK(UV,FLG,PT,VN,VT,COND,DUV)

	USE GLOBAL
	
	IMPLICIT NONE

	INTEGER :: COND(2), FLG(2)
	DOUBLE PRECISION :: PT(3), VN(3), VT(3), UV(2,2), DUV(2)
	DOUBLE PRECISION :: TMP1(3), TMP2(3), TMP3(3), GF(3), VNN, ACC
	DOUBLE PRECISION :: NFC, TMP4(4)

	CALL CALGF(PT,GF)
	CALL CROSS(OMG,PT,TMP1)
	CALL CROSS(OMG,TMP1,TMP2)

	TMP1=KAPPA*GF-TMP2
	VNN=SQRT(DOT_PRODUCT(VN,VN))
	TMP3=VN/VNN
	ACC=DOT_PRODUCT(TMP3,TMP1)

	IF (ACC.GE.0.0D0) THEN
	  COND(1)=1
	ELSE IF (-VNN/ACC.GT.CRS*DT) THEN
	  COND(1)=1
	ELSE
	  COND(1)=0
	END IF

	CALL QUFUN(UV(2,1),UV(2,2),TMP1)
	CALL QVFUN(UV(2,1),UV(2,2),TMP2)
	CALL LINSV(TMP1,TMP2,VT,DUV(1),DUV(2))
	TMP4(1:2)=UV(2,:)
	TMP4(3:4)=DUV
	CALL NFORCE(TMP4,NFC)

	IF (NFC.GT.0.0D0) THEN
	  COND(2)=0
	ELSE
	  COND(2)=1
	END IF

	RETURN

	END SUBROUTINE CMCHK


!****************************************************************************
!
!	SUBROUTINE SFSCHK: STATIC FRICTION STICKING CHECKING
!
!****************************************************************************

	SUBROUTINE SFSCHK(X,QU,QV,FR,N)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: X(4), QU(3), QV(3), FR, N
	DOUBLE PRECISION :: Q(3), QUU(3), QUV(3), QVV(3), QN(3)
	DOUBLE PRECISION :: GF(3), TMP1(3), TMP2(3), TMP3(3), TMP4(3)

	CALL QFUN(X(1),X(2),Q)
	CALL QUUFUN(X(1),X(2),QUU)
	CALL QUVFUN(X(1),X(2),QUV)
	CALL QVVFUN(X(1),X(2),QVV)

	CALL CALGF(Q,GF)
	CALL CROSS(OMG,Q,TMP1)
	CALL CROSS(OMG,TMP1,TMP2)

	TMP1=X(3)*QU+X(4)*QV
	CALL CROSS(OMG,TMP1,TMP3)
	TMP4=QUU*X(3)**2.0D0+2.0D0*QUV*X(3)*X(4)+QVV*X(4)**2.0D0
	TMP1=KAPPA*GF-TMP2-2.0D0*TMP3-TMP4

	CALL CROSS(QU,QV,QN)
	QN=QN/SQRT(DOT_PRODUCT(QN,QN))
	
	N=-DOT_PRODUCT(QN,TMP1)
	TMP2=TMP1+QN*N

    FR=SQRT(DOT_PRODUCT(TMP2,TMP2))
	N=ABS(N)

	RETURN

	END SUBROUTINE SFSCHK