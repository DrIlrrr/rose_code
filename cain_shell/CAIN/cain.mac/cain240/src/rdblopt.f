      SUBROUTINE RDBLOPT(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE BEAMCM
      IMPLICIT NONE
      INTEGER IC,LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INTEGER MOP
      PARAMETER (MOP=6)
      CHARACTER*18 OP(MOP)/'BEAMLINE',
     %   'PERIODIC','BETA','ALPHA','ETA','ETAPRIME'/
      INTEGER NFF(MOP)/1,0,2,2,2,2/
	INTEGER IDBETA,IDALPHA,IDETA,IDETAP,IPPER
	INTEGER IOP(MOP)
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
      INCLUDE 'include/nestcm.h'
      INTEGER J,I,NF,NC,K,LPER,IRTN1
      INTEGER NCBLNAM
	CHARACTER(32) BLNAM
	REAL(8) BAI(2,2),BAF(2,2),EPI(4),EPF(4)
C
      IRTN=0
      IF(LN(1,2).EQ.0) GOTO 900
      CALL CMDBLK('BLOPTICS',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) GOTO 900
	
	IOP=0
	NCBLNAM=0
      J=1
      DO I=1,MOP
        ID(I)=J
	  IF(OP(I).EQ.'BETA') IDBETA=ID(I)
	  IF(OP(I).EQ.'ALPHA') IDALPHA=ID(I)
	  IF(OP(I).EQ.'ETA') IDETA=ID(I)
	  IF(OP(I).EQ.'ETAPRIME') IDETAP=ID(I)
	  IF(OP(I).EQ.'PERIODIC') IPPER=I
        J=J+MAX(0,NFF(I))
      ENDDO
      DO I=1,J-1
        PAR(I)=UNDEF
      ENDDO
	NGSTRRD=0
C
      DO 300 J=1,NBL
        I=KBL(J)
	  IF(NFF(I).EQ.0) THEN
	    IOP(I)=1
        ELSEIF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
	    IF(OP(I).EQ.'BEAMLINE') THEN
	      CALL APSOFF(TEXT(1:NC),NC)
	      NCBLNAM=NC
            BLNAM=TEXT(1:NC)
	    ELSE
            DO K=1,NFF(I)
	        FFF(K)=UNDEF
            ENDDO
            CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
            IF(IRTN.NE.0) GOTO 990
            IF(NF.GT.NFF(I)) GOTO 910
            IF(NF.GE.1) THEN
              DO K=1,NF
                IF(FFF(K).NE.UNDEF) THEN
                  PAR(ID(I)+K-1)=FFF(K)
                ENDIF
              ENDDO
            ENDIF
	    ENDIF
        ENDIF
 300  CONTINUE
	IF(NCBLNAM.LE.0) GOTO 920
	LPER=IOP(IPPER)
	IF(LPER.EQ.0) THEN
	  DO I=1,2
	    BAI(1,I)=0
	    EPI(2*I-1)=0
	    EPI(2*I)=0
	    IF(PAR(IDBETA+I-1).NE.UNDEF) BAI(1,I)=PAR(IDBETA+I-1)%X
	    IF(PAR(IDALPHA+I-1).NE.UNDEF) BAI(2,I)=PAR(IDALPHA+I-1)%X
	    IF(PAR(IDETA+I-1).NE.UNDEF) EPI(2*I-1)=PAR(IDETA+I-1)%X
	    IF(PAR(IDETAP+I-1).NE.UNDEF) EPI(2*I)=PAR(IDETAP+I-1)%X
	  ENDDO
	  IF(BAI(1,1).LE.0.OR.BAI(1,2).LE.0) GOTO 950
	ENDIF
	CALL BLOPTICS(LPER,BLNAM(1:NCBLNAM),BAI,BAF,EPI,EPF,IRTN1)
	IF(IRTN1.NE.0) GOTO 990

      RETURN
C
 900  IRTN=1000
      WRITE(MSGFL,905)
 905  FORMAT(' (SUBR.RDBLOPT) No parameter specified.')
      RETURN
 910  IRTN=1001
      WRITE(MSGFL,915) OP(I)
 915  FORMAT(' (SUBR.RDBLOPT) Too many numbers for ',
     %  'operand "',A,'".')
      RETURN
 920  IRTN=1002
      WRITE(MSGFL,925)
 925  FORMAT(' (SUBR.RDBLOPT) Beamline name not specified.')
      RETURN
 950  IRTN=1005
      WRITE(MSGFL,955)
 955  FORMAT(' (SUBR.RDBLOPT)  Invalid Twiss parameter.')
      RETURN
 990  IRTN=1009
      RETURN
      END
