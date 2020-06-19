!****************************************************************************
!
!  SUBROUTINE RESET: RESET THE COORDINATES TO MASS CENTER PRINCIPLE AXIS FRAME
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE RESET()

	USE GLOBAL

	USE IMSL

	IMPLICIT NONE

	INTEGER :: I

	DOUBLE PRECISION :: JMX1(3,3), RCMX(3,1), RCLEN, TMP1(3,3), TMP2(3,3), TRANMX(3,3), EIGS(3), TMPV(3), TMPE

	DO I=1,3

	  RCMX(I,1)=RC0(I)

	END DO

	RCLEN=NORM(RC0)

	TMP1=(RCLEN**2)*EYE(3)

	TMP2=RCMX.XT.RCMX

	JMX1=JMX0-MASS*(TMP1-TMP2)

	EIGS=EIG(JMX1,V=TRANMX)

	IF ((EIGS(1)>EIGS(2)).AND.(EIGS(1)>EIGS(3))) THEN

	  TMPV=TRANMX(:,3)
	  TRANMX(:,3)=TRANMX(:,1)
	  TRANMX(:,1)=TMPV

	  TMPE=EIGS(3)
	  EIGS(3)=EIGS(1)
	  EIGS(1)=TMPE

	  IF (EIGS(2)<EIGS(1)) THEN

	  TMPV=TRANMX(:,2)
	  TRANMX(:,2)=TRANMX(:,1)
	  TRANMX(:,1)=TMPV

	  TMPE=EIGS(2)
	  EIGS(2)=EIGS(1)
	  EIGS(1)=TMPE

	  END IF

	ELSE IF (EIGS(2)>EIGS(3)) THEN
	  
	  TMPV=TRANMX(:,3)
	  TRANMX(:,3)=TRANMX(:,2)
	  TRANMX(:,2)=TMPV

	  TMPE=EIGS(3)
	  EIGS(3)=EIGS(2)
	  EIGS(2)=TMPE

	  IF (EIGS(2)<EIGS(1)) THEN

	  TMPV=TRANMX(:,2)
	  TRANMX(:,2)=TRANMX(:,1)
	  TRANMX(:,1)=TMPV

      TMPE=EIGS(2)
	  EIGS(2)=EIGS(1)
	  EIGS(1)=TMPE

	  END IF

	ELSE
	  
	  IF(EIGS(2)<EIGS(1)) THEN

	  TMPV=TRANMX(:,2)
	  TRANMX(:,2)=TRANMX(:,1)
	  TRANMX(:,1)=TMPV

	  TMPE=EIGS(2)
	  EIGS(2)=EIGS(1)
	  EIGS(1)=TMPE

	  END IF

	END IF

	TMPE=DET(TRANMX)

	IF (TMPE<0.0) THEN

	  TRANMX(:,1)=-TRANMX(:,1)

	END IF

	JMX=DIAG(EIGS)

	DO I=1,N_V

	  TMPV=V_DATA(I,:)
	  TMPV=TMPV-RC0
	  V_DATA(I,:)=TRANMX.TX.TMPV

	END DO

	RETURN

	END SUBROUTINE RESET
