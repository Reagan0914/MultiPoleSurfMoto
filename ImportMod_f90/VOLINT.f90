!****************************************************************************
!
!  SUBROUTINE VOLINT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE VOLINT()

	USE GLOBAL

	IMPLICIT NONE

	CALL COMPVOLUMEINT

	MASS=DENSITY*T0

	RC0=T1/T0

	JMX0(1,1)=DENSITY*(T2(2)+T2(3))
	JMX0(2,2)=DENSITY*(T2(3)+T2(1))
	JMX0(3,3)=DENSITY*(T2(1)+T2(2))
	JMX0(1,2)=-DENSITY*T2(4)
	JMX0(2,1)=-DENSITY*T2(4)
	JMX0(2,3)=-DENSITY*T2(5)
	JMX0(3,2)=-DENSITY*T2(5)
	JMX0(1,3)=-DENSITY*T2(6)
	JMX0(3,1)=-DENSITY*T2(6)

	RETURN

	END SUBROUTINE VOLINT

!****************************************************************************
!
!  SUBROUTINE COMPVOLUMEINT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE COMPVOLUMEINT()

	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: I,A,B,C

	DOUBLE PRECISION :: NX,NY,NZ,TMP

	DOUBLE PRECISION :: WFS(N_F), NFS(N_F,3), F(12)!F=(FA,FB,FC,FAA,FBB,FCC,FAAA,FBBB,FCCC,FAAB,FBBC,FCCA)

	T0=0.0D0
	T1=(/0.0D0, 0.0D0, 0.0D0/)
	T2=(/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)

	CALL COMPWFNF(WFS,NFS)

	DO I=1,N_F

	  NX=ABS(NFS(I,1))
	  NY=ABS(NFS(I,2))
	  NZ=ABS(NFS(I,3))
	  IF ((NX>NY).AND.(NX>NZ)) THEN 
	    C=1
	  ELSE IF (NY>NZ) THEN 
	    C=2
	  ELSE 
	    C=3
	  END IF
	  A=MOD(C,3)+1
	  B=MOD(A,3)+1

	  CALL COMPFACEINT(WFS(I),NFS(I,:),I,A,B,C,F)

	  IF (A==1) THEN
	    TMP=F(1)
	  ELSE IF (B==1) THEN
	    TMP=F(2)
	  ELSE
	    TMP=F(3)
	  END IF

	  T0=T0+NFS(I,1)*TMP

	  T1(A)=T1(A)+NFS(I,A)*F(4)
	  T1(B)=T1(B)+NFS(I,B)*F(5)
	  T1(C)=T1(C)+NFS(I,C)*F(6)
	  T2(A)=T2(A)+NFS(I,A)*F(7)
	  T2(B)=T2(B)+NFS(I,B)*F(8)
	  T2(C)=T2(C)+NFS(I,C)*F(9)
	  T2(3+A)=T2(3+A)+NFS(I,A)*F(10)
	  T2(3+B)=T2(3+B)+NFS(I,B)*F(11)
	  T2(3+C)=T2(3+C)+NFS(I,C)*F(12)

	END DO

	T1=T1/2.0D0
	T2(1:3)=T2(1:3)/3.0D0
	T2(4:6)=T2(4:6)/2.0D0

	RETURN

	END SUBROUTINE COMPVOLUMEINT

!****************************************************************************
!
!  SUBROUTINE COMPWFNF: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE COMPWFNF(WFS,NFS)

	USE GLOBAL

	USE IMSL

	IMPLICIT NONE

	INTEGER :: I

	DOUBLE PRECISION :: WFS(N_F),NFS(N_F,3),RA(3),RB(3),RC(3),TMP(3)

	DO I=1,N_F
	  
	  RA=V_DATA(F_DATA(I,1),:)
	  RB=V_DATA(F_DATA(I,2),:)
	  RC=V_DATA(F_DATA(I,3),:)
	  CALL CROSS(RB-RA,RC-RB,TMP)
	  TMP=TMP/SQRT(DOT_PRODUCT(TMP,TMP))
	  WFS(I)=-DOT_PRODUCT(RA,TMP)
	  NFS(I,:)=TMP
	  
	END DO

	RETURN

	END SUBROUTINE COMPWFNF

!****************************************************************************
!
!  SUBROUTINE COMPFACEINT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE COMPFACEINT(W,NF,ID,A,B,C,F)

	IMPLICIT NONE

	INTEGER :: A,B,C,ID

	DOUBLE PRECISION :: W,NF(3),F(12),P(10)

	DOUBLE PRECISION :: K(4)

	CALL COMPPROJECTIONINT(A,B,ID,P)

	K(1)=1.0D0/NF(C)
	K(2)=K(1)*K(1)
	K(3)=K(2)*K(1)
	K(4)=K(3)*K(1)

	F(1)=K(1)*P(2)
	F(2)=K(1)*P(3)
	F(3)=-K(2)*(NF(A)*P(2)+NF(B)*P(3)+W*P(1))

	F(4)=K(1)*P(4)
	F(5)=K(1)*P(6)
	F(6)=K(3)*((NF(A)**2)*P(4)+2.0D0*NF(A)*NF(B)*P(5)&
		 +(NF(B)**2)*P(6)+W*(2.0D0*(NF(A)*P(2)+NF(B)*P(3))+W*P(1)))

	F(7)=K(1)*P(7)
	F(8)=K(1)*P(10)
	F(9)=-K(4)*((NF(A)**3)*P(7)+3.0D0*(NF(A)**2)*NF(B)*P(8)& 
	     +3.0D0*NF(A)*(NF(B)**2)*P(9)+(NF(B)**3)*P(10)&
	     +3.0D0*W*((NF(A)**2)*P(4)+2.0D0*NF(A)*NF(B)*P(5)+(NF(B)**2)*P(6))&
	     +W*W*(3.0D0*(NF(A)*P(2)+NF(B)*P(3))+W*P(1)))

	F(10)=K(1)*P(8)
	F(11)=-K(2)*(NF(A)*P(9)+NF(B)*P(10)+W*P(6))
	F(12)=K(3)*((NF(A)**2)*P(7)+2.0D0*NF(A)*NF(B)*P(8)+(NF(B)**2)*P(9)&
	      +W*(2.0D0*(NF(A)*P(4)+NF(B)*P(5))+W*P(2)))

	RETURN

	END SUBROUTINE COMPFACEINT

!****************************************************************************
!
!  SUBROUTINE COMPPROJECTIONINT: 
!
!  INPUT: 
!
!  OUTPUT: 
!
!****************************************************************************

	SUBROUTINE COMPPROJECTIONINT(A,B,ID,P)

	USE GLOBAL

	IMPLICIT NONE

	INTEGER :: A,B,ID,I
	DOUBLE PRECISION :: P(10)
	DOUBLE PRECISION :: A0,A1,DA
	DOUBLE PRECISION :: B0,B1,DB
	DOUBLE PRECISION :: A0_2,A0_3,A0_4,B0_2,B0_3,B0_4
	DOUBLE PRECISION :: A1_2,A1_3,B1_2,B1_3
	DOUBLE PRECISION :: C1,CA,CAA,CAAA,CB,CBB,CBBB
	DOUBLE PRECISION :: CAB,KAB,CAAB,KAAB,CABB,KABB

	P=(/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)

	DO I=1,3

	  A0=V_DATA(F_DATA(ID,I),A)
      B0=V_DATA(F_DATA(ID,I),B)
      A1=V_DATA(F_DATA(ID,MOD(I,3)+1),A)
      B1=V_DATA(F_DATA(ID,MOD(I,3)+1),B)
      DA=A1-A0
      DB=B1-B0

      A0_2=A0*A0
	  A0_3=A0_2*A0
	  A0_4=A0_3*A0

      B0_2=B0*B0
	  B0_3=B0_2*B0
	  B0_4=B0_3*B0

      A1_2=A1*A1
	  A1_3=A1_2*A1 
      B1_2=B1*B1
	  B1_3=B1_2*B1

      C1=A1+A0
      CA=A1*C1+A0_2
	  CAA=A1*CA+A0_3
	  CAAA=A1*CAA+A0_4
      CB=B1*(B1+B0)+B0_2
	  CBB=B1*CB+B0_3
	  CBBB=B1*CBB+B0_4
      CAB=3.0D0*A1_2+2.0D0*A1*A0+A0_2
	  KAB=A1_2+2.0D0*A1*A0+3.0D0*A0_2
      CAAB=A0*CAB+4.0D0*A1_3
	  KAAB=A1*KAB+4.0D0*A0_3
      CABB=4.0D0*B1_3+3.0D0*B1_2*B0+2.0D0*B1*B0_2+B0_3
      KABB=B1_3+2.0D0*B1_2*B0+3.0D0*B1*B0_2+4.0D0*B0_3

!P=(P1,PA,PB,PAA,PAB,PBB,PAAA,PAAB,PABB,PBBB)

      P(1)=P(1)+DB*C1
      P(2)=P(2)+DB*CA
      P(4)=P(4)+DB*CAA
      P(7)=P(7)+DB*CAAA
      P(3)=P(3)+DA*CB
      P(6)=P(6)+DA*CBB
      P(10)=P(10)+DA*CBBB
      P(5)=P(5)+DB*(B1*CAB+B0*KAB)
      P(8)=P(8)+DB*(B1*CAAB+B0*KAAB)
      P(9)=P(9)+DA*(A1*CABB+A0*KABB)

	END DO

	P(1)=P(1)/2.0D0
    P(2)=P(2)/6.0D0
    P(4)=P(4)/12.0D0
    P(7)=P(7)/20.0D0
    P(3)=-P(3)/6.0D0
    P(6)=-P(6)/12.0D0
    P(10)=-P(10)/20.0D0
    P(5)=P(5)/24.0D0
    P(8)=P(8)/60.0D0
    P(9)=-P(9)/60.0D0

	RETURN

	END SUBROUTINE COMPPROJECTIONINT
