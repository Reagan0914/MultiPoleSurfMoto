!****************************************************************************
!
!  SUBROUTINE CALGF: CALCULATE THE GRAVITY, UNDIMENSIONAL
!
!****************************************************************************

	SUBROUTINE CALGF(FPT,GF)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: I, J
	INTEGER :: IDCR(N_F)

	DOUBLE PRECISION :: FPT(3), GF(3)	
	DOUBLE PRECISION :: LPT(N_V,4), ERLE(N_E,3), FROF(N_F,3)	
	DOUBLE PRECISION :: TMP, LE, RE(3), TMPF(3,4), FR(3,3), OMEGAF, TMPV(3), F1(3), F2(3), D

	DO I=1,N_V

	  LPT(I,1:3)=VS(I,:)-FPT
	  LPT(I,4)=SQRT(DOT_PRODUCT(LPT(I,1:3),LPT(I,1:3)))

	END DO

	DO I=1,N_F
	  
	  IDCR(I)=0
	
	END DO

	DO I=1,N_E

	  CALL SINGCHK(I,FPT,D)	 
	   
	  IF (D.LT.SERR) THEN

	    DO J=1,3		  
		  ERLE(I,J)=0.0D0
		END DO
		DO J=3,4
		  IDCR(ES(I,J))=1
		END DO	
				  
	  ELSE

	    TMP=SUM(LPT(ES(I,1:2),4))
		LE=LOG((TMP+LES(I))/(TMP-LES(I)))
		RE=LPT(ES(I,1),1:3)
	    ERLE(I,:)=NFS(ES(I,3),:)*DOT_PRODUCT(RE,NES(I,1:3))&
	              +NFS(ES(I,4),:)*DOT_PRODUCT(RE,NES(I,4:6))
        ERLE(I,:)=ERLE(I,:)*LE

	  END IF	  

	END DO

	DO I=1,3

	  F1(I)=SUM(ERLE(:,I))

	END DO

	DO I=1,N_F

	  IF (IDCR(I).EQ.0) THEN

	    TMPF=LPT(FS(I,:),:)

	    DO J=1,3
	      FR(J,:)=TMPF(J,1:3)/TMPF(J,4)
	    END DO

	    CALL CROSS(FR(2,:),FR(3,:),TMPV)

	    TMP=1.0D0+DOT_PRODUCT(FR(1,:),FR(2,:))+DOT_PRODUCT(FR(2,:),FR(3,:))+DOT_PRODUCT(FR(3,:),FR(1,:))
	    OMEGAF=2.0D0*ATAN2(DOT_PRODUCT(TMPV,FR(1,:)),TMP)
	    FROF(I,:)=DOT_PRODUCT(TMPF(1,1:3),NFS(I,:))*NFS(I,:)
	    FROF(I,:)=FROF(I,:)*OMEGAF

	  ELSE

	    DO J=1,3
		  FROF(I,J)=0.0D0
		END DO

	  END IF

	END DO

	DO I=1,3

	  F2(I)=SUM(FROF(:,I))

	END DO

	GF=-F1+F2

	RETURN

	END SUBROUTINE CALGF

!****************************************************************************
!
!  SUBROUTINE CALPE: CALCULATE THE SURFACE GRAVITAIONAL POTENTIAL, UNDIMENSIONAL
!
!****************************************************************************

	SUBROUTINE CALPE(FPT,PE)
	
	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: I, J
	INTEGER :: IDCR(N_F)

	DOUBLE PRECISION :: FPT(3), PE	
	DOUBLE PRECISION :: LPT(N_V,4), RERE(N_E), LE(N_E), RFRF(N_F), OMEGAF(N_F)	
	DOUBLE PRECISION :: U1, U2, TMP, RE(3), TMPF(3,4), FR(3,3), TMPV(3), D

	DO I=1,N_V

	  LPT(I,1:3)=VS(I,:)-FPT
	  LPT(I,4)=SQRT(DOT_PRODUCT(LPT(I,1:3),LPT(I,1:3)))

	END DO

	DO I=1,N_F
	  
	  IDCR(I)=0
	
	END DO

	DO I=1,N_E

	  CALL SINGCHK(I,FPT,D)	 
	   
	  IF (D.LT.SERR) THEN
		  
		RERE(I)=0.0D0
		LE(I)=0.0D0

		DO J=3,4
		  IDCR(ES(I,J))=1
		END DO	
				  
	  ELSE

	    TMP=SUM(LPT(ES(I,1:2),4))
	    LE(I)=LOG((TMP+LES(I))/(TMP-LES(I)))
	    RE=LPT(ES(I,1),1:3)
	    RERE(I)=DOT_PRODUCT(NES(I,1:3),RE)*DOT_PRODUCT(NFS(ES(I,3),:),RE)&
			   +DOT_PRODUCT(NES(I,4:6),RE)*DOT_PRODUCT(NFS(ES(I,4),:),RE)

	  END IF	  

	END DO

	U1=DOT_PRODUCT(RERE,LE)

	DO I=1,N_F

	  IF (IDCR(I).EQ.0) THEN

	    TMPF=LPT(FS(I,:),:)

	    DO J=1,3
	      FR(J,:)=TMPF(J,1:3)/TMPF(J,4)
	    END DO

	    CALL CROSS(FR(2,:),FR(3,:),TMPV)

	    TMP=1.0D0+DOT_PRODUCT(FR(1,:),FR(2,:))+DOT_PRODUCT(FR(2,:),FR(3,:))+DOT_PRODUCT(FR(3,:),FR(1,:))
	    OMEGAF(I)=2.0D0*ATAN2(DOT_PRODUCT(TMPV,FR(1,:)),TMP)
	    RFRF(I)=DOT_PRODUCT(TMPF(1,1:3),NFS(I,:))**2

	  ELSE

	    OMEGAF(I)=0.0D0
	    RFRF(I)=0.0D0

	  END IF

	END DO

	U2=DOT_PRODUCT(RFRF,OMEGAF)

	PE=-0.5D0*(U1-U2)

	RETURN

	END SUBROUTINE CALPE

!****************************************************************************
!
!  SUBROUTINE SINGCHK: SINGULARITY CHECK 
!
!****************************************************************************

	SUBROUTINE SINGCHK(I,FPT,D)

	USE GLOBAL
	
	IMPLICIT NONE

	INTEGER :: I

	DOUBLE PRECISION :: LAM, FPT(3), D, TMP1, TMP2(3)

	LAM=DOT_PRODUCT(FPT-VS(ES(I,1),:),VS(ES(I,2),:)-VS(ES(I,1),:))	  

	LAM=LAM/LES(I)**2.0D0	  

	IF (LAM.LT.0.0D0) THEN

       TMP1=-LAM*LES(I)	
	
	ELSE IF (LAM.GT.1.0D0) THEN
	     
	   
	   TMP1=(LAM-1.0D0)*LES(I)
	
	ELSE
	   
	   TMP1=0.0D0
	
	END IF

	TMP2=FPT-LAM*VS(ES(I,2),:)+(LAM-1.0D0)*VS(ES(I,1),:)

	D=SQRT(TMP1**2.0D0+DOT_PRODUCT(TMP2,TMP2))
	
	RETURN

	END SUBROUTINE SINGCHK