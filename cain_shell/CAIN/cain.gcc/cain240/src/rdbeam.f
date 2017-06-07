      SUBROUTINE RDBEAM(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE BEAMCM
	USE ARRAYMOD
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
C      INCLUDE 'include/beamcm.h'
      INCLUDE 'include/nestcm.h'
      INCLUDE 'include/pushcm.h'
      INCLUDE 'include/tstpcm.h'
      INCLUDE 'include/cnstcm.h'
      INTEGER MOP
      PARAMETER (MOP=42)
      CHARACTER*16 OP(MOP)/
     %  'RIGHT','LEFT','TESTPARTICLE','SINGLEPARTICLE',
     %  'ELLIPTIC','TUNIFORM','EUNIFORM','NAME',
     %  'FILE','KIND','NP','AN','TXYS',
     %  'BETA','ALPHA','EMIT','E0','WEIGHT',
     %  'SIGE','SIGT','GCUT','GCUTT','GCUTE','GAUSSWEIGHT',
     %  'SLOPE','CRAB','ETA','ETAPRIME','ESLOPE','XYROLL',
     %  'DALPHADE','DALPHADT','SPIN','P','NAMELIST','USERDEFINED',
     %  'BEGIN','END','TERMINATE',
C           The order of these 3 keyword must be kept
     %  'COMMENT','KEYWORD','CONVERSION'/
      INTEGER NFF(MOP)/0,0,0,0,
     %   0,0,0,1, 
     %   1,1,1,1,4,
     %   2,2,2,1,1,
     %   1,1,2,1,1,1,
     %   2,2,2,2,1,1,
     %   2,2,3,3,0,0,
     %   1,1,1,
     %   1,1,1/
      INTEGER IDELLI,IDKIND,
     %   IDNP,IDAN,IDTXYS,IDBETA,IDALPH,
     %   IDEMIT,IDE0,IDWGT,IDSIGE,IDSIGT,IDGCUT,IDGCTT,IDGCTE,IDGWGT,
     %   IDSLOP,IDCRAB,IDETA,IDETAP,IDESLP,IDXYRL,IDDALPHA,IDDALPHT,
     %   IDSPIN,IDP,
     %   IDFLG,IDCMNT,IPFLG
      INTEGER ISINGLE,ITEST,LR,IGWGT,LR1(2),LDIST(3),KIND1,FILE,NP1,
     %   NAMLST,LUSER,NCTXT2(2),ISTAT,NP2
      REAL*8 AN,TXYS0(0:3),BETA(2),ALPHA(2),EMIT(2),E0,SIGE,SIGT,
     %   GCUT(4),SLOPE(2),CRAB(2),ETA(2),ETAP(2),ESLOPE,
     %   XYROLL,DALPHADE(2),DALPHADT(2),
     %   SPIN1(3),EPTST(0:3),DUMMY32(3,2),WGT1
      INTEGER LR2(2)/1,1/,KIN2(3)/1,1,1/
      REAL*8 BFL1(3,2)/0,0,0,0,0,0/,GCUT00(4)/3.5D0,3.5D0,3.5D0,3.5D0/
C      
      INTEGER N,J,NF,I,NC,K,NCFN
	TYPE(FLCHTYPE) FC
      CHARACTER*4 NAM1
      CHARACTER*512 FILENAME
	CHARACTER(MCTEXT) TEXT2(2)
      CHARACTER*256 ERR
	CHARACTER(1) CMNT
	INTEGER, PARAMETER :: MFLG=3
	INTEGER NCFLG(MFLG)
	CHARACTER(16) FLG(MFLG)
      CHARACTER*8 KIN(3)/'photon','electron','positron'/
C
      IRTN=0
      IF(LN(1,2).EQ.0) GOTO 1000
      CALL CMDBLK('BEAM',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %   LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) GOTO 1000
C set default
      ITEST=0
	ISINGLE=0
      NAM1='    '
      FILENAME=' '
	NCFN=1
      NAMLST=0
	LUSER=0
      LR1(1)=0
      LR1(2)=0
      LDIST(1)=0
      LDIST(2)=0
      LDIST(3)=0
	NCTXT2=0
	NCFLG=0
	FILE=0
      J=1
      DO I=1,MOP
        ID(I)=J
        IF(OP(I).EQ.'ELLIPTIC') IDELLI=ID(I)
        IF(OP(I).EQ.'KIND') IDKIND=ID(I)
        IF(OP(I).EQ.'NP') IDNP=ID(I)
        IF(OP(I).EQ.'AN') IDAN=ID(I)
        IF(OP(I).EQ.'TXYS') IDTXYS=ID(I)
        IF(OP(I).EQ.'BETA') IDBETA=ID(I)
        IF(OP(I).EQ.'ALPHA') IDALPH=ID(I)
        IF(OP(I).EQ.'EMIT') IDEMIT=ID(I)
        IF(OP(I).EQ.'E0') IDE0=ID(I)
	  IF(OP(I).EQ.'WEIGHT') IDWGT=ID(I)
        IF(OP(I).EQ.'SIGE') IDSIGE=ID(I)
        IF(OP(I).EQ.'SIGT') IDSIGT=ID(I)
        IF(OP(I).EQ.'GCUT') IDGCUT=ID(I)
        IF(OP(I).EQ.'GCUTT') IDGCTT=ID(I)
        IF(OP(I).EQ.'GCUTE') IDGCTE=ID(I)
        IF(OP(I).EQ.'GAUSSWEIGHT') IDGWGT=ID(I)
        IF(OP(I).EQ.'SLOPE') IDSLOP=ID(I)
        IF(OP(I).EQ.'CRAB') IDCRAB=ID(I)
        IF(OP(I).EQ.'ETA') IDETA=ID(I)
        IF(OP(I).EQ.'ETAPRIME') IDETAP=ID(I)
        IF(OP(I).EQ.'ESLOPE') IDESLP=ID(I)
        IF(OP(I).EQ.'XYROLL') IDXYRL=ID(I)
	  IF(OP(I).EQ.'DALPHADE') IDDALPHA=ID(I)
	  IF(OP(I).EQ.'DALPHADT') IDDALPHT=ID(I)
        IF(OP(I).EQ.'SPIN') IDSPIN=ID(I)
        IF(OP(I).EQ.'P') IDP=ID(I)
	  IF(OP(I).EQ.'BEGIN') IDFLG=ID(I)
	  IF(OP(I).EQ.'BEGIN') IPFLG=I
	  IF(OP(I).EQ.'COMMENT') IDCMNT=ID(I)
        J=J+MAX(0,NFF(I))
      ENDDO
      DO I=1,J-1
        PAR(I)=ZERO
      ENDDO
	NGSTRRD=0
C
      DO 300 J=1,NBL
        I=KBL(J)
        IF(NFF(I).EQ.0) THEN
          IF(I.LE.2) THEN
            LR1(I)=1
	    ELSEIF(OP(I).EQ.'SINGLEPARTICLE') THEN
	      ISINGLE=1
          ELSEIF(OP(I).EQ.'TESTPARTICLE') THEN
            ITEST=1
          ELSEIF(OP(I).EQ.'NAMELIST') THEN
            NAMLST=1
	    ELSEIF(OP(I).EQ.'USERDEFINED') THEN
            LUSER=1
          ELSEIF(I.GE.IDELLI.AND.I.LE.IDELLI+2) THEN
            LDIST(I-IDELLI+1)=1
          ENDIF
        ELSEIF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          IF(OP(I).EQ.'FILE') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      IF(FC%L.EQ.1) THEN
	        FILE=NINT(FC%X)
	      ELSE
	        FILENAME=GSTR2(EVALLAST)(FC%C(1):FC%C(2))
	        NCFN=FC%C(2)-FC%C(1)+1
            ENDIF
	      GOTO 300
	    ELSEIF(OP(I).EQ.'NAME') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      IF(FC%L.EQ.1) THEN
	        WRITE(NAM1,'(A,I3)') 'T',NINT(FC%X)
	      ELSE
	        NAM1='T'//GSTR2(EVALLAST)(FC%C(1):FC%C(2))
            ENDIF
	      GOTO 300
	    ELSEIF(OP(I).EQ.'KEYWORD') THEN
	      NCTXT2(1)=NC
	      TEXT2(1)=TEXT(1:NC)
	      GOTO 300
	    ELSEIF(OP(I).EQ.'CONVERSION') THEN
	      NCTXT2(2)=NC
	      TEXT2(2)=TEXT(1:NC)
	      GOTO 300
          ENDIF
          DO 220 K=1,NFF(I)
            FFF(K)=UNDEF
 220      CONTINUE
          CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
          IF(IRTN.NE.0) GOTO 990
          IF(NF.GT.NFF(I)) GOTO 900
          IF(NF.GE.1) THEN
            DO 240 K=1,NF
              IF(FFF(K).NE.UNDEF) PAR(ID(I)+K-1)=FFF(K)
 240        CONTINUE
          ENDIF
        ENDIF
 300  CONTINUE
      IF(ITEST.NE.0.OR.ISINGLE.NE.0) GOTO 700
      IF(FILE.NE.0.OR.FILENAME.NE.' ') GOTO 600
      IF(LR1(1).EQ.LR1(2)) GOTO 920
      IF(LR1(1).EQ.1) LR=1
      IF(LR1(2).EQ.1) LR=2
      KIND1=NINT(PAR(IDKIND)%X)
      NP1=NINT(PAR(IDNP)%X)
      AN=PAR(IDAN)%X
      DO 320 J=0,3
        TXYS0(J)=PAR(IDTXYS+J)%X
 320  CONTINUE
      DO 340 J=1,2
        BETA(J)=PAR(IDBETA+J-1)%X
        ALPHA(J)=PAR(IDALPH+J-1)%X
        EMIT(J)=PAR(IDEMIT+J-1)%X
        GCUT(J)=PAR(IDGCUT+J-1)%X
        SLOPE(J)=PAR(IDSLOP+J-1)%X
        CRAB(J)=PAR(IDCRAB+J-1)%X
        ETA(J)=PAR(IDETA+J-1)%X
        ETAP(J)=PAR(IDETAP+J-1)%X
	  DALPHADE(J)=PAR(IDDALPHA+J-1)%X
	  DALPHADT(J)=PAR(IDDALPHT+J-1)%X
 340  CONTINUE
      E0=PAR(IDE0)%X
      SIGE=PAR(IDSIGE)%X
      SIGT=PAR(IDSIGT)%X
      GCUT(3)=PAR(IDGCTT)%X
      GCUT(4)=PAR(IDGCTE)%X
      IGWGT=NINT(PAR(IDGWGT)%X)
      XYROLL=PAR(IDXYRL)%X
      ESLOPE=PAR(IDESLP)%X
      DO 360 J=1,3
        SPIN1(J)=PAR(IDSPIN+J-1)%X
 360  CONTINUE
      ERR=' '
      N=0
      IF(KIND1.EQ.0) CALL ERMSG1(ERR,N,'KIND')
      IF(NP1.EQ.0) CALL ERMSG1(ERR,N,'NP')
      IF(AN.EQ.0) CALL ERMSG1(ERR,N,'AN')
      IF(E0.EQ.0) CALL ERMSG1(ERR,N,'E0')
      IF(BETA(1).LE.0.OR.BETA(2).LE.0) CALL ERMSG1(ERR,N,'BETA')
C      IF(EMIT(1).LE.0.OR.EMIT(2).LE.0) CALL ERMSG1(ERR,N,'EMIT')
      IF(ERR.NE.' ') GOTO 910
      DO 440 I=1,4
        IF(GCUT(I).LE.0) GCUT(I)=GCUT00(I)
 440  CONTINUE
      CALL BMINI(LR,KIND1,NP1,AN,TXYS0,LDIST,BETA,ALPHA,
     %  EMIT,E0,SIGE,SIGT,SLOPE,CRAB,GCUT,IGWGT,ETA,ETAP,ESLOPE,
     %  XYROLL,DALPHADE,DALPHADT,SPIN1,IRTN)
      IF(IRTN.NE.0) GOTO 990
      IRTN=0
      GOTO 1000
C--- Read from File ---
 600  NP1=NINT(PAR(IDNP)%X)
	IF(FILENAME(1:NCFN).NE.' ') FILE=98
	CALL OPENFL(FILE,FILENAME(1:NCFN),'OLD',0,NC,IRTN)
	IF(IRTN.NE.0) GOTO 940
	IF(LUSER.EQ.0) THEN
        CALL BMFILE(NP1,FILE,NAMLST,IRTN)
	ELSE
	  IF(PAR(IDCMNT).NE.ZERO) THEN
	    IF(PAR(IDCMNT)%L.NE.2) GOTO 932
	    NC=PAR(IDCMNT)%C(2)-PAR(IDCMNT)%C(1)+1
	    IF(NC.GE.2) GOTO 932
	    IF(NC.NE.0) CMNT=GSTRRD(PAR(IDCMNT)%C(1):PAR(IDCMNT)%C(2))
	    IF(CMNT.EQ.' '.OR.CMNT.EQ.',') GOTO 932
	  ELSE
	    NC=0
	  ENDIF
	  DO I=1,MFLG
	    IF(PAR(IDFLG+I-1).NE.ZERO) THEN
	      IF(PAR(IDFLG+I-1)%L.NE.2) GOTO 934
	      NCFLG(I)=PAR(IDFLG+I-1)%C(2)-PAR(IDFLG+I-1)%C(1)+1
	      FLG(I)=GSTRRD(PAR(IDFLG+I-1)%C(1):PAR(IDFLG+I-1)%C(2))
	      DO J=1,NCFLG(I)
	        IF(FLG(I)(J:J).EQ.' '.OR.FLG(I)(J:J).EQ.',') GOTO 936
	      ENDDO
	    ENDIF
	  ENDDO
C	  CALL BMFILE2(NP1,FILE,TEXT2,NCTXT2,CMNT(1:NC),MFLG,FLG,NCFLG,
C     %    GSTRRD,NGSTRRD,TEXT,NP2,IRTN)
	  CALL BMFILE2C(NP1,FILE,TEXT2,NCTXT2,CMNT(1:NC),MFLG,FLG,NCFLG,
     %    TEXT,NP2,IRTN)
	  IF(IRTN.EQ.0.AND.MSGLVL.GE.1) THEN
	    IF(FILENAME(1:NCFN).NE.' ') THEN
		    WRITE(MSGFL,610) NP2,FILENAME(1:NCFN)
610         FORMAT(' +++ BEAM',I6,' particles read from file "',A,'"')
          ELSE
		    WRITE(MSGFL,620) NP2,FILE
620         FORMAT(' +++ BEAM',I6,' particles read from file #',I2)
          ENDIF
	  ENDIF
	ENDIF
	IF(FILENAME(1:NCFN).NE.' ') CLOSE(FILE,IOSTAT=ISTAT)
      IF(IRTN.NE.0) GOTO 990
      GOTO 1000
C--- Test Particle or Single Particle ---
 700  K=NINT(PAR(IDKIND)%X)
	IF(K.LE.0.OR.K.GE.4) GOTO 950
	DO J=0,3
        TXYS0(J)=PAR(IDTXYS+J)%X
      ENDDO
	EPTST(0)=MASS(K)**2
	DO J=1,3
	  EPTST(J)=PAR(IDP+J-1)%X
	  EPTST(0)=EPTST(0)+EPTST(J)**2
	  SPIN1(J)=PAR(IDSPIN+J-1)%X
	ENDDO
	EPTST(0)=SQRT(EPTST(0))
	IF(ITEST.NE.0) THEN
        IF(NAM1.EQ.'    ') GOTO 930
        CALL ADDTSTP(K,NAM1,TXYS0,EPTST(0),SPIN1,IRTN)
	ELSE
	  WGT1=PAR(IDWGT)%X
	  IF(WGT1.LT.0) GOTO 954
	  CALL ADDONE(0,K,1,'    ',0,WGT1,TXYS0,EPTST,SPIN1,
     %       0,DUMMY32,IRTN)
	ENDIF
      IF(IRTN.NE.0) GOTO 990
      GOTO 1000
C-----------------------
 900  IRTN=1000
      WRITE(MSGFL,905) OP(I)
 905  FORMAT(' (SUBR.RDBEAM) Too many numbers for ',
     %  'operand "',A,'".')
      GOTO 1000
 910  IRTN=1010
      WRITE(MSGFL,915) ERR(1:N)
 915  FORMAT(' (SUBR.RDBEAM) Following beam parameters ',
     %  'not specified.',/,5X,A)
      GOTO 1000
 920  IRTN=1020
      WRITE(MSGFL,925)
 925  FORMAT(' (SUBR.RDBEAM) Either one of RIGHT,LEFT,FILE ',
     %  'must specified.',/,5X,A)
      GOTO 1000
 930  IRTN=1030
      WRITE(MSGFL,931)
 931  FORMAT(' (SUBR.RDBEAM) Particle name not given for BEAM ',
     %  'TESTPARTICLE')
      GOTO 1000
 932  IRTN=1032
      WRITE(MSGFL,933)
 933  FORMAT(' (SUBR.RDBEAM) Invalid COMMENT character.')
      GOTO 1000
 934  IRTN=1034
      WRITE(MSGFL,935) OP(IPFLG+I-1)
 935  FORMAT(' (SUBR.RDBEAM) Invalid string for ',A)
      GOTO 1000
 936  IRTN=1036
      WRITE(MSGFL,937) FLG(I)(1:NCFLG(I)),OP(IPFLG+I-1)
 937  FORMAT(' (SUBR.RDBEAM) String "',A,'" must not contain a comma ',
     %  'and a blanck space',/,'  for operand ',A)
      GOTO 1000
 940  IRTN=1040
      IF(FILENAME(1:NCFN).EQ.' ') THEN
	  WRITE(MSGFL,942) FILE
 942    FORMAT(' (SUBR.RDBEAM) File#',I3,' open error.')
	ELSE
	  WRITE(MSGFL,944) FILENAME(1:NCFN)
 944    FORMAT(' (SUBR.RDBEAM) File "',A,'" open error.')
	ENDIF
	GOTO 1000
 950  IRTN=1050
      WRITE(MSGFL,952) K
 952  FORMAT(' (SUBR.RDBEAM) Invalid particle KIND=',I3)
      GOTO 1000
 954  IRTN=1054
      WRITE(MSGFL,956) WGT1
 956  FORMAT(' (SUBR.RDBEAM) Invalid particle WEIGHT=',1PD12.4)
      GOTO 1000
 960  IRTN=1060
      WRITE(MSGFL,965) ERR(1:80)
 965  FORMAT(' (SUBR.RDBEAM) Invalid expression for file name.',/,
     %   3X,A)
      GOTO 1000
 990  IRTN=1090
1000  RETURN
      END
