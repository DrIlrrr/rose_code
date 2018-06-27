      SUBROUTINE RDTRANS(IC,LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE BEAMCM
      IMPLICIT NONE
      INTEGER IC,LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INTEGER MOP
      PARAMETER (MOP=10)
      CHARACTER*18 OP(MOP)/'BEAMLINE','MOMENTUM','CHARGE',
     %   'TXYS','E3','E1','RIGHT','LEFT','KIND','LOSSMONITOR'/
      INTEGER NFF(MOP)/1,1,1,4,3,3,0,0,3,0/
	INTEGER IDMOM,IDCHARGE,IDTXYS,IDE3,IDE1,IDKIND,IPRIGHT,IPLEFT
	INTEGER IOP(MOP)
      INCLUDE 'include/ctrlcm.h'
	INCLUDE 'include/cnstcm.h'
C      INCLUDE 'include/readcm.h'
      INCLUDE 'include/nestcm.h'
      INTEGER J,I,N,L,NF,NC,K,IFIRST,LSPIN,LOSSMON
      INTEGER LR(2),KIN(3),NCBLNAM
	CHARACTER(32) BLNAM
	INTEGER IQ000
	REAL(8) E3(3),E1(3),TXYS0(0:3),P000
C
      IRTN=0
      IF(NESTRT.EQ.0) THEN
        IF(NESTLV.NE.0) THEN
          DO 120 I=1,NESTLV-1
            IF(NEST(I).EQ.NEST_TRANS.OR.NEST(I).EQ.NEST_PUSH) GOTO 940
 120      CONTINUE
        ENDIF
        IF(NESTLV.GE.MNEST) GOTO 950
        NESTLV=NESTLV+1
        NEST(NESTLV)=NEST_TRANS
        ICNEST(NESTLV)=IC
      ELSE
        NESTRT=0
        IFIRST=0
        GOTO 500
      ENDIF
      IRTN=0
      IF(LN(1,2).EQ.0) GOTO 900
      CALL CMDBLK('TRANSPORT',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) GOTO 900
	
	IOP=0
	KIN=0
	NCBLNAM=0
	LOSSMON=0
      J=1
      DO I=1,MOP
        ID(I)=J
	  IF(OP(I).EQ.'TXYS') IDTXYS=ID(I)
	  IF(OP(I).EQ.'MOMENTUM') IDMOM=ID(I)
	  IF(OP(I).EQ.'CHARGE') IDCHARGE=ID(I)
	  IF(OP(I).EQ.'E3') IDE3=ID(I)
	  IF(OP(I).EQ.'E1') IDE1=ID(I)
	  IF(OP(I).EQ.'KIND') IDKIND=ID(I)
	  IF(OP(I).EQ.'RIGHT') IPRIGHT=I
	  IF(OP(I).EQ.'LEFT') IPLEFT=I
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
	    IF(OP(I).EQ.'LOSSMONITOR') LOSSMON=1
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
	LR(1)=IOP(IPRIGHT)
      LR(2)=IOP(IPLEFT)
      CALL ALTRUE(LR,2)
      DO I=1,3
        IF(PAR(IDKIND+I-1).NE.UNDEF) THEN
          K=NINT(PAR(IDKIND+I-1)%X)
          IF(K.GE.1.AND.K.LE.3) KIN(K)=1
        ENDIF
      ENDDO
      CALL ALTRUE(KIN,3)
	KIN(1)=0
C        exclude photons
	P000=0
	IQ000=0
	IF(PAR(IDMOM).NE.UNDEF) P000=PAR(IDMOM)%X
	IF(PAR(IDCHARGE).NE.UNDEF) IQ000=NINT(PAR(IDCHARGE)%X)
	IF(P000.EQ.0.OR.IQ000.EQ.0) GOTO 930
	DO I=0,3
	  TXYS0(I)=0
	  IF(PAR(IDTXYS+I).NE.UNDEF) TXYS0(I)=PAR(IDTXYS+I)%X
	ENDDO
	E3=0
	E1=0
	E3(3)=1
	E1(1)=1
	DO I=1,3
	  IF(PAR(IDE3+I-1).NE.UNDEF) E3(I)=PAR(IDE3+I-1)%X
	  IF(PAR(IDE1+I-1).NE.UNDEF) E1(I)=PAR(IDE1+I-1)%X
	ENDDO
C  Asign particles to be excluded. ISBIN is used as a work area
C  untill the TRANSPORT loop finishes.
	DO N=1,NP
	  ISBIN(N)=1
	  IF(KIN(KIND(N)).EQ.0) THEN
	    ISBIN(N)=0
	  ELSE
	    L=1
	    IF(EP(3,N).LE.0) L=2
	    IF(LR(L).EQ.0) THEN
	      ISBIN(N)=0
	    ENDIF
	  ENDIF
	ENDDO
	LSPIN=ISPIN
      IFIRST=1
      INPUSH=1
 500  CALL TRANSPORT(IFIRST,LSPIN,BLNAM(1:NCBLNAM),LOSSMON,
     %   IQ000,P000,EMASS,TXYS0,E1,E3,
     %   NP,ISBIN,KIND,LOST,PNAME,TXYS,EP,SPIN,IRTN)
      IF(IRTN.NE.0) GOTO 990
      RETURN
C
 900  IRTN=1000
      WRITE(MSGFL,905)
 905  FORMAT(' (SUBR.RDTRANS) No parameter specified.')
      RETURN
 910  IRTN=1001
      WRITE(MSGFL,915) OP(I)
 915  FORMAT(' (SUBR.RDTRANS) Too many numbers for ',
     %  'operand "',A,'".')
      RETURN
 920  IRTN=1002
      WRITE(MSGFL,925)
 925  FORMAT(' (SUBR.RDTRANS) Beamline name not specified.')
      RETURN
 930  IRTN=1003
      WRITE(MSGFL,935)
 935  FORMAT(' (SUBR.RDTRANS) MOMENTUM and/or CHARGE not specified.')
      RETURN
 940  IRTN=1004
      WRITE(MSGFL,945)
 945  FORMAT(' (SUBR.RDTRANS) PUSH/TRANSPORT loop cannot nest.')
      RETURN
 950  IRTN=1008
      WRITE(MSGFL,955) MNEST
 955  FORMAT(' (SUBR.RDTRANS) Nest level too deep. ',/,
     %  ' Sum of PUSH, TRANSPORT, IF, DO nest must be <=',I2)
      RETURN
 990  IRTN=1009
      RETURN
      END
