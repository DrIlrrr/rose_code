	SUBROUTINE MATCHSET(IRTN)
	USE FLCHTYP
	USE BEAMLN
	USE ARRAYMOD
	USE MATCHMOD
	IMPLICIT NONE
	INTEGER IRTN
	INCLUDE 'include/ctrlcm.h'
	REAL(8) BAI(2,2),BAF(2,2),EPI(4),EPF(4)
	TYPE(FLCHTYPE) FC
	CHARACTER(80) ERR
	INTEGER I,J,K

	IF(LPER.EQ.0) THEN
	  DO I=1,8
	    IF(TWISSIN(I)%L.EQ.2) THEN
	      CALL EVAL0(GSTRMATCH(TWISSIN(I)%C(1):TWISSIN(I)%C(2)),
     %        FC,ERR)
	      IF(ERR.NE.' ') GOTO 990
	      TWISSIN(I)%X=FC%X
	    ENDIF
	    IF(I.LE.2.AND.TWISSIN(I)%X.LE.0) GOTO 940
	  ENDDO
	ENDIF
	DO I=1,2
	  BAI(1,I)=TWISSIN(I)%X
	  BAI(2,I)=TWISSIN(I+2)%X
	  EPI(2*I-1)=TWISSIN(I+4)%X
	  EPI(2*I)=TWISSIN(I+6)%X
	ENDDO
C--- Set magnet parameters
	CALL MAGSET(IDBLMATCH,ERR)
	IF(ERR.NE.' ') GOTO 990
C--- Compute optics
	CALL BLOPTICS(LPER,BLNAM,BAI,BAF,EPI,EPF,IRTN)
	IF(IRTN.NE.0) GOTO 1000
	IRTN=0
	GOTO 1000

940   IRTN=1040
	IF(MSGLVL.GE.0) WRITE(MSGFL,941)
941   FORMAT(' (SUBR.MATCHSET) Invalid beta function at entrance.')
	GOTO 1000
990   IRTN=1090
	IF(MSGLVL.GE.0) WRITE(MSGFL,992) ERR
992   FORMAT(' (SUBR.MATCHSET) ',A)
	GOTO 1000
1000  RETURN
	END