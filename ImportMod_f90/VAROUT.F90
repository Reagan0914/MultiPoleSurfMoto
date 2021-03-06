
!****************************************************************************
!
!  SUBROUTINE PREVARSOUT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE PREVARSOUT()

	USE GLOBAL 

	IMPLICIT NONE

	INTEGER :: I

	OPEN(UNIT=127, FILE="NUMS.TXT")

	WRITE(127,"(3(I12))") N_V, N_E, N_F 

	CLOSE(127)

	OPEN(UNIT=177, FILE="VS.TXT")

	DO I=1,N_V

	  WRITE(177,"(3(F22.12))") V_DATA(I,:)

	END DO

	CLOSE(177)
	
	OPEN(UNIT=178, FILE="NES.TXT")

	DO I=1,N_E

	  WRITE(178,"(6(F22.12))") NE_DATA(I,:)

	END DO

	CLOSE(178)

	OPEN(UNIT=179, FILE="LES.TXT")

	DO I=1,N_E

	  WRITE(179,"(F22.12)") LE_DATA(I)

	END DO

	CLOSE(179)

	OPEN(UNIT=180, FILE="NFS.TXT")

	DO I=1,N_F

	  WRITE(180,"(3(F22.12))") NF_DATA(I,:)

	END DO

	CLOSE(180)
	
	OPEN(UNIT=181, FILE="FS.TXT")

	DO I=1,N_F

	  WRITE(181,"(3(I12))") F_DATA(I,:)

	END DO

	CLOSE(181)

	OPEN(UNIT=182, FILE="ES.TXT")

	DO I=1,N_E

	  WRITE(182,"(4(I12))") E_DATA(I,:)

	END DO

	CLOSE(182)

	OPEN(UNIT=183, FILE="SECS.TXT")

	DO I=1,N_F

	  WRITE(183,"(6(I12))") SEC_DATA(I,1,1),SEC_DATA(I,1,2),SEC_DATA(I,2,1),&
	                        SEC_DATA(I,2,2),SEC_DATA(I,3,1),SEC_DATA(I,3,2)
	END DO

	CLOSE(183)

	RETURN

	END SUBROUTINE PREVARSOUT

!****************************************************************************
!
!  SUBROUTINE PRELOADLOG: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE PRELOADLOG()

	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: I

	TRULER=PERIOD
	KAPPA=TRULER**2.0D0*G*DENSITY

	OPEN(UNIT=164, FILE="LOG.TXT")

	WRITE(164,*) "NUMBERS OF VERTEX, FACE, EDGE OF THE POLYHEDRON:"
	WRITE(164,"(3(I10))") N_V,N_F,N_E

	WRITE(164,*) "SPATIAL SCALE, TIME SCALE, DYNAMICAL KAPPA:"
	WRITE(164,"(3(3X,E20.12E3))") SRULER,TRULER,KAPPA

	WRITE(164,*) "IN INITIAL FIXED FRAME, ZEROTH MOMENT T0:"
	WRITE(164,"(3X,E20.12E3)") T0
	
	WRITE(164,*) "IN INITIAL FIXED FRAME, FIRST MONENT T1=(TX,TY,TZ):"
	WRITE(164,"(3(3X,E20.12E3))") T1
	
	WRITE(164,*) "IN FIXED FRAME, SECOND MONENT T2=(TXX,TYY,TZZ,TXY,TYZ,XZ):"	
	WRITE(164,"(6(3X,E20.12E3))") T2

	WRITE(164,*) "THE VOLUME OF POLYHEDRON:"
	WRITE(164,"(3X,E20.12E3)") T0

	WRITE(164,*) "THE BULK DENSITY OF POLYHEDRON:"
	WRITE(164,"(3X,E20.12E3)") DENSITY

	WRITE(164,*) "THE MASS OF POLYHEDRON:"
	WRITE(164,"(3X,E20.12E3)") MASS

	WRITE(164,*) "IN INITIAL FIXED FRAME, MASS CENTER COORDINATES RC0=(XC0,YC0,ZC0):"
	WRITE(164,"(3(3X,E20.12E3))") RC0

	WRITE(164,*) "IN INITIAL FIXED FRAME, INERTIAL MATRIX J0:"

	DO I=1,3

	  WRITE(164,"(3(3X,E20.12E3))") JMX0(I,:)

	END DO

	WRITE(164,*) "IN RESETED FIXED FRAME, MASS CENTER COORDINATES RC=(XC,YC,ZC):"
	WRITE(164,"(3(3X,E20.12E3))") 0.0D0, 0.0D0, 0.0D0

	WRITE(164,*) "IN RESETED FIXED FRAME, INERTIAL MATRIX J:"

	DO I=1,3

	  WRITE(164,"(3(3X,E20.12E3))") JMX(I,:)

	END DO

	CLOSE(164)

	RETURN

	END SUBROUTINE PRELOADLOG



