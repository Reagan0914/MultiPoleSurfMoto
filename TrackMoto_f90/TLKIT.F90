!****************************************************************************
!
!  SUBROUTINE CROSS: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE CROSS(X,Y,Z)

	IMPLICIT NONE

	DOUBLE PRECISION :: X(3),Y(3),Z(3)

	Z=(/X(2)*Y(3)-X(3)*Y(2),X(3)*Y(1)-X(1)*Y(3),X(1)*Y(2)-X(2)*Y(1)/)

	RETURN

	END SUBROUTINE CROSS

!****************************************************************************
!
!  SUBROUTINE MULPT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE MULPT(X,Y,Z)

	IMPLICIT NONE

	INTEGER I

	DOUBLE PRECISION :: X(2300),Y(2300,1),TMP(1),TMPZZ,Z

	Z=0

	DO I=1,2300

	TMPZZ=X(I)*Y(I,1)

	Z=Z+TMPZZ

	END DO

	RETURN

	END SUBROUTINE MULPT
!****************************************************************************
!
!  SUBROUTINE LINSV: SOLVE LINEAR SYSTEM A(3)*X+B(3)*Y=C(3)
!                   WHERE A, B, C ARE COPLANAR, A, B ARE INDEPENDENT 
!
!****************************************************************************

	SUBROUTINE LINSV(AV,BV,CV,X,Y)

	IMPLICIT NONE

	DOUBLE PRECISION :: AV(3),BV(3),CV(3),X,Y
	DOUBLE PRECISION :: DET,T1,T2,T3,F1,F2

	T1=DOT_PRODUCT(AV,AV)
	T2=DOT_PRODUCT(AV,BV)
	T3=DOT_PRODUCT(BV,BV)

	F1=DOT_PRODUCT(AV,CV)
	F2=DOT_PRODUCT(BV,CV)

	DET=T1*T3-T2**2
		
	X=(T3*F1-T2*F2)/DET
	Y=(-T2*F1+T1*F2)/DET

	RETURN

	END SUBROUTINE LINSV


