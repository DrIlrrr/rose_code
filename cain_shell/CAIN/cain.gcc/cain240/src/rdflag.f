      SUBROUTINE RDFLAG(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INCLUDE 'include/ctrlcm.h'
      INCLUDE 'include/ctrlcm2.h'
C      INCLUDE 'include/readcm.h'
      INTEGER MOP
      PARAMETER (MOP=MFLAG+2)
      INTEGER NFF(MOP)
      DATA NFF/MOP*-1/
      INTEGER ONOFF,I,J
C
      IF(LN(1,2).EQ.0) GOTO 400
      CALL CMDBLK('FLAG',MOP,FLAGNM,100,NFF,MBL,NBL,KOP,KBL,REL,
     %     LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.2) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) GOTO 400
      ONOFF=1
      DO 300 J=1,NBL
        I=KBL(J)
        IF(FLAGNM(I).EQ.'ON') THEN
          ONOFF=1
        ELSEIF(FLAGNM(I).EQ.'OFF') THEN
          ONOFF=0
        ELSE
          IFLAG(I-2)=ONOFF
        ENDIF
 300  CONTINUE
 400  IRTN=0
      RETURN
 990  RETURN
      END
