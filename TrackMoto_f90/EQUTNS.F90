!****************************************************************************
!
!  SUBROUTINE ORBODE: 
!
!****************************************************************************

	SUBROUTINE ORBODE(T, Y, YP)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: T, Y(6), YP(6)
	DOUBLE PRECISION :: GF(3), TMP1(3), TMP2(3)

	CALL CALGF(Y(1:3),GF)
	CALL CROSS(OMG,Y(1:3),TMP1)
	CALL CROSS(OMG,TMP1,TMP2)
	CALL CROSS(OMG,Y(4:6),TMP1)

	YP(1:3)=Y(4:6)
	YP(4:6)=KAPPA*GF-TMP1*2.0D0-TMP2

	RETURN
	
	END SUBROUTINE ORBODE

!****************************************************************************
!
!  SUBROUTINE SURFODE: 
!	U=Y(1)	V=Y(2)	DU=Y(3)	DV=Y(4)
!
!****************************************************************************

	SUBROUTINE SURFODE(T,Y,YP)

	USE GLOBAL

	IMPLICIT NONE

	DOUBLE PRECISION :: T, Y(4), YP(4)
	DOUBLE PRECISION :: Q(3), QU(3), QV(3), QUU(3), QUV(3), QVV(3) 
	DOUBLE PRECISION :: GF(3), TMP1(3), TMP2(3), TMP3(3), TMP4(3), F(3), R1, R2

	CALL QFUN(Y(1),Y(2),Q)
	CALL QUFUN(Y(1),Y(2),QU)
	CALL QVFUN(Y(1),Y(2),QV)
	CALL QUUFUN(Y(1),Y(2),QUU)
	CALL QUVFUN(Y(1),Y(2),QUV)
	CALL QVVFUN(Y(1),Y(2),QVV)

	CALL CALGF(Q,GF)
	CALL CROSS(OMG,Q,TMP1)
	CALL CROSS(OMG,TMP1,TMP2)

	TMP1=Y(3)*QU+Y(4)*QV
	CALL CROSS(OMG,TMP1,TMP3)
	TMP4=QUU*Y(3)**2.0D0+2.0D0*QUV*Y(3)*Y(4)+QVV*Y(4)**2.0D0
	F=KAPPA*GF-TMP2-MU*TMP1-2.0D0*TMP3-TMP4

	CALL LINSV(QU,QV,F,R1,R2)

	YP(1:2)=Y(3:4)		
	YP(3)=R1
	YP(4)=R2

	RETURN
	
	END SUBROUTINE SURFODE


