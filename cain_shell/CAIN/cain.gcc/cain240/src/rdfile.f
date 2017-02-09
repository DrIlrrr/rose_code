      SUBROUTINE RDFILE(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE ARRAYMOD
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INTEGER MOP
      PARAMETER (MOP=7)
      CHARACTER*12 OP(MOP)/'OPEN','CLOSE','REWIND','ENDFILE',
     %   'UNIT','STATUS','NAME'/
      INTEGER NFF(MOP)/-1,-1,-1,-1,1,1,1/
      INTEGER IDUNIT,IDSTAT,IDNAME
      INTEGER UNIT,NCFN,IOSTAT
	TYPE(FLCHTYPE) FC
      CHARACTER*7 STATUS
      CHARACTER*512 FNAME
	CHARACTER(80) ERR
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
      INTEGER J,NF,I,NC,K,II
C
      IRTN=0
      IF(LN(1,2).EQ.0) RETURN
      CALL CMDBLK('FILE',MOP,OP,1,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) RETURN
      J=1
      DO 180 I=1,MOP
        ID(I)=J
        J=J+MAX(0,NFF(I))
        IF(OP(I).EQ.'UNIT') IDUNIT=ID(I)
        IF(OP(I).EQ.'STATUS') IDSTAT=ID(I)
        IF(OP(I).EQ.'NAME') IDNAME=ID(I)
 180  CONTINUE
      DO 190 I=1,J-1
        PAR(I)=UNDEF
 190  CONTINUE
	NGSTRRD=0
      STATUS=' '
      FNAME=' '
      UNIT=0
C
      II=KBL(1)
      DO 300 J=2,NBL
        I=KBL(J)
        IF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          IF(OP(I).EQ.'STATUS') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' '.OR.FC%L.NE.2) GOTO 930
            STATUS=GSTR2(EVALLAST)(FC%C(1):FC%C(2))
          ELSEIF(OP(I).EQ.'NAME') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' ') GOTO 970
	      IF(FC%L.NE.2) GOTO 980
	      FNAME=GSTR2(EVALLAST)(FC%C(1):FC%C(2))
	      NCFN=FC%C(2)-FC%C(1)+1
	      GOTO 300
          ELSE
            DO 220 K=1,NFF(I)
	        FFF(K)=UNDEF
 220        CONTINUE
            CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
            IF(IRTN.NE.0) GOTO 990
            IF(NF.GT.NFF(I)) GOTO 900
            IF(NF.GE.1) THEN
              DO 240 K=1,NF
                IF(FFF(K).NE.UNDEF) PAR(ID(I)+K-1)=FFF(K)
 240          CONTINUE
            ENDIF
          ENDIF
        ENDIF
 300  CONTINUE
      IF(PAR(IDUNIT).NE.UNDEF) UNIT=NINT(PAR(IDUNIT)%X)
      IF(UNIT.EQ.0) GOTO 910
      GOTO (320,340,360,380), II
C-- OPEN
 320  IF(STATUS.EQ.' ') STATUS='UNKNOWN'
      CALL OPENFL(UNIT,FNAME(1:NCFN),STATUS,0,NCFN,IRTN)
      IF(IRTN.NE.0) GOTO 920
      IF(MSGLVL.GE.1) THEN
        WRITE(MSGFL,330) UNIT,STATUS
 330    FORMAT(' +++ OPEN ',I3,' STATUS=',A)
        IF(FNAME.NE.' ') WRITE(MSGFL,335) FNAME(1:NCFN)
 335    FORMAT('       Filename=',A)
      ENDIF
      GOTO 800
C-- CLOSE
 340  IF(STATUS.EQ.' ') STATUS='KEEP'
      CLOSE(UNIT=UNIT,STATUS=STATUS,IOSTAT=IOSTAT)
      IF(IOSTAT.NE.0) GOTO 940
      IF(MSGLVL.GE.1) WRITE(MSGFL,350) UNIT,STATUS
 350  FORMAT(' +++ CLOSE ',I3,' STATUS=',A)
      GOTO 800
C-- REWIND
 360  REWIND (UNIT=UNIT,IOSTAT=IOSTAT)
      IF(IOSTAT.NE.0) GOTO 950
      IF(MSGLVL.GE.1) WRITE(MSGFL,370) UNIT
 370  FORMAT(' +++ REWIND ',I3,' +++')
      GOTO 800
C-- ENDFILE
 380  ENDFILE (UNIT=UNIT,IOSTAT=IOSTAT)
      IF(IOSTAT.NE.0) GOTO 960
      IF(MSGLVL.GE.1) WRITE(MSGFL,390) UNIT
 390  FORMAT(' +++ ENDFILE ',I3,' +++')
      GOTO 800
C---
 800  IRTN=0
      RETURN
C
 900  IRTN=1000
      WRITE(MSGFL,905) OP(I)
 905  FORMAT(' (SUBR.RDFILE) Too many numbers for ',
     %  'operand "',A,'".')
      RETURN
 910  IRTN=1001
      WRITE(MSGFL,915)
 915  FORMAT(' (SUBR.RDFILE) UNIT operand not specified.')
      RETURN
 920  IRTN=1002
      RETURN
 930  IRTN=1003
      WRITE(MSGFL,935) TEXT(1:NC)
 935  FORMAT(' (SUBR.RDFILE) Invalid character expression "',A,
     %   '" for STATUS.')
      RETURN
 940  IRTN=1004
      WRITE(MSGFL,945) IOSTAT
 945  FORMAT(' (SUBR.RDFILE) CLOSE failed.  IOSTAT=',I5)
      WRITE(MSGFL,350) UNIT,STATUS
      RETURN
 950  IRTN=1005
      WRITE(MSGFL,955) UNIT,IOSTAT
 955  FORMAT(' (SUBR.RDFILE) REWIND unit=',I3,' failed.  IOSTAT=',I5)
      RETURN
 960  IRTN=1006
      WRITE(MSGFL,965) UNIT,IOSTAT
 965  FORMAT(' (SUBR.RDFILE) ENDFILE unit=',I3,' failed.  IOSTAT=',I5)
      RETURN
 970  IRTN=1007
      WRITE(MSGFL,975) ERR
 975  FORMAT(' (SUBR.RDFILE) Invalid character expression for file ',
     %  'name.',/,3X,A)
      RETURN
 980  IRTN=1008
      WRITE(MSGFL,985)
 985  FORMAT(' (SUBR.RDFILE) Invalid character expression for file ',
     %  'name.')
      RETURN
 990  IRTN=1009
      RETURN
      END