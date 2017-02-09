	SUBROUTINE IFARRAY(NAM,ID)
C        ID<0 if character array
	USE ARRAYMOD
	IMPLICIT NONE
	CHARACTER(*) NAM
	INTEGER ID
	INTEGER I,NC,I0,I1

	ID=0
	NC=LEN(NAM)
	IF(NC.LE.0) RETURN
	I0=0
	DO I=1,NC
	  IF(NAM(I:I).NE.' ') THEN
	    IF(I0.EQ.0) I0=I
	    I1=I
	  ENDIF
	ENDDO
	IF(I0.EQ.0) RETURN
	IF(NARRAY.GT.0) THEN
	  DO I=1,NARRAY
	    IF(ARR(I)%TYPE.NE.0) THEN
	      IF(NAM(I0:I1).EQ.ARR(I)%NAME) THEN
	        ID=I
	        IF(ARR(I)%TYPE.EQ.2) ID=-I
	        RETURN
	      ENDIF
	    ENDIF
	  ENDDO
	ENDIF
	RETURN
	END
