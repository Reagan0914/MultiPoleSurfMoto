!****************************************************************************
!
!  SUBROUTINE CROSS: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE CROSS(A,B,C)

	IMPLICIT NONE

	DOUBLE PRECISION :: A(3),B(3),C(3)

	C=(/A(2)*B(3)-A(3)*B(2),A(3)*B(1)-A(1)*B(3),A(1)*B(2)-A(2)*B(1)/)

	RETURN

	END SUBROUTINE CROSS