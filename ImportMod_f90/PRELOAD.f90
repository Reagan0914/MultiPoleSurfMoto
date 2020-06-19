!****************************************************************************
!
!  PROGRAM: PRELOAD
!
!  PURPOSE: 
!
!****************************************************************************

	PROGRAM PRELOAD

	USE GLOBAL

	IMPLICIT NONE

	WRITE(*,*) "MODEL DATA FILE INPUT ..."
	CALL FILEIN

	WRITE(*,*) "VOLUME INTEGRATION COMPUTATING ..."
	CALL VOLINT

	WRITE(*,*) "COORDINATES RESET TO MASS CENTER PRINCIPLE AXIS FRAME ..."
	CALL RESET

	WRITE(*,*) "GENERATE VARIANTS IN ADVANCE ..."
	CALL GENVAR
	 
    WRITE(*,*) "PRINT PRELOADLOG ..."
	CALL PRELOADLOG
	
	WRITE(*,*) "PRINT PREGENERATED VARIANTS ..."
	CALL PREVARSOUT

    WRITE(*,*) "PROGRAM PRELOAD FINISHED!"

	END PROGRAM PRELOAD




	
