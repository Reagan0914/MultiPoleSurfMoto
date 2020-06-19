!****************************************************************************
!
!  SUBROUTINE FILELOAD: LOAD THE VERTICES COORDINATES AND FACES VERTICES DATA
!
!  INPUT: FNAME-THE ASTEROID POLYHEDRON FILE NAME
!
!  OUTPUT: 
!  THERE MUST NOT BE EMPTY ROW IN THE FILE
!
!****************************************************************************

	SUBROUTINE FILEIN()

	USE GLOBAL	
	
	IMPLICIT NONE

	INTEGER :: I,TMP

	OPEN(UNIT=5,FILE=FILENAME,STATUS='OLD')

	READ(5,*) N_V,N_F

	ALLOCATE(V_DATA(N_V,3))

	ALLOCATE(F_DATA(N_F,3))

	DO I=1,N_V

	   READ(5,*) TMP,V_DATA(I,1),V_DATA(I,2),V_DATA(I,3)

	END DO

	DO I=1,N_F

       READ(5,*) TMP,F_DATA(I,1),F_DATA(I,2),F_DATA(I,3)

	END DO

	CLOSE(5)

	RETURN

	END SUBROUTINE FILEIN
