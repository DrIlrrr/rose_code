	SUBROUTINE OPENLOGO
	INCLUDE 'include/headcm.h'
	INCLUDE 'include/ctrlcm.h'
	CHARACTER(512) FN
	CHARACTER(80) ERR
	INTEGER I,NC,NC1,NC2,ID1,ID2,IDUMMY

	WRITE(MSGFL,100) VERSION,VERSIONDATE
100   FORMAT(' ***** ',A,5X,A,' *****')
	INQUIRE(RDFL,NAME=FN)
	IF(OSID.EQ.1) CALL WINFN2UNIX(FN)
C           Change backslash to forwardslash on Windows.
	NC=0
	DO I=LEN(FN),1,-1
	  IF(FN(I:I).NE.' ') THEN
	    NC=I
	    EXIT
	  ENDIF
	ENDDO
	IF(NC.GT.0) THEN
	  NC1=0
	  NC2=NC
	  DO I=NC,1,-1
	    IF(FN(I:I).EQ.'/') THEN
	      NC2=NC-I
	      NC1=I
	      EXIT
	    ENDIF
	  ENDDO
	  CALL IFARRAY('$InFilePath',ID1)
	  CALL IFARRAY('$InFileName',ID2)
	  ID1=-ID1
	  ID2=-ID2
	  CALL EVARRCSET(ID2,0,IDUMMY,FN(NC-NC2+1:NC),ERR)
	  IF(NC1.NE.0) CALL EVARRCSET(ID1,0,IDUMMY,FN(1:NC1),ERR)
	  WRITE(MSGFL,200) FN(1:NC)
200     FORMAT(15X,'Input File = ',A)
	ENDIF
      RETURN
	END
