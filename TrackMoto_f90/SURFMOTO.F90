!****************************************************************************
!
!	PROGRAM: 
!
!	PURPOSE: 
!
!****************************************************************************

	PROGRAM SURFMOTO

	IMPLICIT NONE

	WRITE(*,*) "LOADING NUMBERS ..."
	CALL LOADNUMS

	WRITE(*,*) "LOADING SURFACE VARIANTS ..."
	CALL LOADSVAR

	WRITE(*,*) "LOADING GRAVITY VARIANTS ..."
	CALL LOADGVAR
	
    WRITE(*,*) "ESTIMATING VALUE OF DELTA TIME ..."
    CALL GETDT

	WRITE(*,*) "CALCULATING THE LOCAL RESPOSE OF SEISMIC ..."
	CALL LOCSEISMIC()

	WRITE(*,*) "FINISHED!"

	END PROGRAM SURFMOTO

	
	
	
