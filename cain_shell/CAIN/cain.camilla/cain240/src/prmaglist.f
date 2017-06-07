      SUBROUTINE PRMAGLIST(P00,FILE)
C  Print list of all magnets
	USE BEAMLN
	IMPLICIT NONE
	INTEGER FILE
	REAL(8) P00
	CHARACTER(1000) TEXT
	INTEGER NC,NC1,I
	REAL(8) X
	REAL(8) PI/3.141592653589793238D0/,CVEL/2.99792458D8/

	IF(NMAG.EQ.0) RETURN
	TEXT=' name     L(m)   ang(rad)  K1(1/m) rot(deg)'
     %  //' tan(edge(1,2))   apert(mm) '
	NC=71
	IF(P00.NE.0) THEN
	  TEXT(NC+1:NC+16)="  B(T)   B'(T/m)"
	  NC=NC+16
	ENDIF
	WRITE(FILE,100) TEXT(1:NC)
100   FORMAT(' *** Defined Magnets ***',/,A)

	DO I=1,NMAG
	  TEXT=MAG(I)%NAME
	  NC=MAXMAGNAME
	  IF(MAG(I)%LENGTH%X.NE.0) 
     %    WRITE(TEXT(NC+1:NC+8),'(F8.4)') MAG(I)%LENGTH%X
	  NC=NC+8
	  IF(MAG(I)%ANGLE.NE.0)
     %    WRITE(TEXT(NC+1:NC+9),'(F9.5)') MAG(I)%ANGLE
	  NC=NC+9
	  IF(MAG(I)%K1%X.NE.0)
     %    WRITE(TEXT(NC+1:NC+10),'(F10.5)') MAG(I)%K1%X
	  NC=NC+10
	  IF(MAG(I)%ROTATE.NE.0)
     %    WRITE(TEXT(NC+1:NC+6),'(F6.1)') MAG(I)%ROTATE*180/PI
	  NC=NC+6
	  IF(MAG(I)%TEDGE(1).NE.0)
     %    WRITE(TEXT(NC+1:NC+9),'(F9.5)') MAG(I)%TEDGE(1)
	  NC=NC+9
	  IF(MAG(I)%TEDGE(2).NE.0)
     %    WRITE(TEXT(NC+1:NC+9),'(F9.5)') MAG(I)%TEDGE(2)
	  NC=NC+9
	  IF(MAG(I)%APERT(1).NE.0)
     %    WRITE(TEXT(NC+1:NC+6),'(F6.1)') MAG(I)%APERT(1)*1000
	  NC=NC+6
	  IF(MAG(I)%APERT(2).NE.0)
     %    WRITE(TEXT(NC+1:NC+6),'(F6.1)') MAG(I)%APERT(2)*1000
	  NC=NC+6
	  IF(P00.NE.0) THEN
	    IF(MAG(I)%ANGLE.NE.0.AND.MAG(I)%LENGTH%X.NE.0) THEN
	      X=P00*MAG(I)%ANGLE/(CVEL*MAG(I)%LENGTH%X)
	      WRITE(TEXT(NC+1:NC+8),'(F8.5)') X
	    ENDIF
	    NC=NC+8
	    IF(MAG(I)%K1%X.NE.0.AND.MAG(I)%LENGTH%X.NE.0) THEN
	      X=P00*MAG(I)%K1%X/(CVEL*MAG(I)%LENGTH%X)
	      WRITE(TEXT(NC+1:NC+8),'(F8.3)') X
	    ENDIF
	    NC=NC+8
	  ENDIF
	  IF(MAG(I)%LENGTH%L.EQ.2) THEN
	    NC1=10+MAG(I)%LENGTH%C(2)-MAG(I)%LENGTH%C(1)+1+1
	    TEXT(NC+1:NC+NC1)='  length="'//
     %            GSTRMG(MAG(I)%LENGTH%C(1):MAG(I)%LENGTH%C(2))//'"'
	    NC=NC+NC1
	  ENDIF
	  IF(MAG(I)%K1%L.EQ.2) THEN
	    NC1=6+MAG(I)%K1%C(2)-MAG(I)%K1%C(1)+1+1
	    TEXT(NC+1:NC+NC1)='  K1="'//
     %            GSTRMG(MAG(I)%K1%C(1):MAG(I)%K1%C(2))//'"'
	    NC=NC+NC1
	  ENDIF
	  WRITE(FILE,'(A)') TEXT(1:NC)
	ENDDO
	RETURN
	END

