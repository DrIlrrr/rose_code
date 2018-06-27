      SUBROUTINE RDREST(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE ARRAYMOD
	USE LUMCOM
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INTEGER MOP
      PARAMETER (MOP=2)
      CHARACTER*12 OP(MOP)/'FILE','LUMINOSITY'/
      INTEGER NFF(MOP)/1,0/
      INTEGER UNIT,NCFN,IOSTAT
      CHARACTER*7 STATUS
      CHARACTER*512 FNAME
	CHARACTER(80) ERR
      CHARACTER*20 HEADER
      INCLUDE 'include/ctrlcm.h'
C      INCLUDE 'include/readcm.h'
      INCLUDE 'include/lumcom2.h'
      INCLUDE 'include/stdstor.h'
      INTEGER J,NF,I,NC,K,LLUM,IRTN1
	TYPE(FLCHTYPE) FC
C
      IRTN=0
      CALL CMDBLK('STORE',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      J=1
      DO 180 I=1,MOP
        ID(I)=J
        J=J+MAX(0,NFF(I))
 180  CONTINUE
      DO 190 I=1,J-1
        PAR(I)=UNDEF
 190  CONTINUE
	NGSTRRD=0
      STATUS=' '
      FNAME=' '
	NCFN=1
	UNIT=0
      LLUM=0
C
      IF(NBL.EQ.0) GOTO 320
      DO 300 J=1,NBL
        I=KBL(J)
        IF(NFF(I).EQ.0) THEN
          IF(OP(I).EQ.'LUMINOSITY') THEN
            LLUM=1
          ENDIF
        ELSEIF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          IF(OP(I).EQ.'FILE') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      IF(FC%L.EQ.1) THEN
	        UNIT=NINT(FC%X)
	      ELSE
	        FNAME=GSTR2(EVALLAST)(FC%C(1):FC%C(2))
	        NCFN=FC%C(2)-FC%C(1)+1
            ENDIF
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
C-- Open file
 320  IF(UNIT.LE.0) THEN
        UNIT=98
        IF(FNAME.EQ.' ') THEN
          IF(LLUM.EQ.0) THEN
            FNAME=STDSTFL
          ELSE
            FNAME=STDSTFLL
          ENDIF
	    CALL SPCOFF(FNAME,3,NCFN)
	  ENDIF
		CALL OPENFL(UNIT,FNAME(1:NCFN),'UNKNOWN',0,NCFN,IRTN)
	  IF(IRTN.NE.0) GOTO 920
        IF(MSGLVL.GE.1) THEN
          IF(LLUM.EQ.0) THEN
            WRITE(MSGFL,410) ':',FNAME(1:NCFN)
 410        FORMAT(' RESTORE',A,' FILE=',1H',A,1H')
          ELSE
            WRITE(MSGFL,410) ' LUMINOSITY:',FNAME(1:NCFN)
          ENDIF
        ENDIF
      ENDIF
C-- Read
      IF(LLUM.EQ.0) THEN
	  CALL RSTORVAR(UNIT,IRTN)
      ELSE
        IF(NLUM.GE.1.AND.MSGLVL.GE.1) WRITE(MSGFL,490)
 490    FORMAT(' (SUBR.RDREST) RESTORE LUMINOSITY: Luminosity ',
     %   'data ',/,'       in this run will be deleted.')
        CALL RESTOLUM(UNIT,HEADER,IRTN1)
        IF(MSGLVL.GE.1) WRITE(MSGFL,500) HEADER
 500    FORMAT(' +++ Luminosity File Header:',A)
        IF(IRTN1.NE.0) THEN
          IF(MSGFL.GE.0) WRITE(MSGFL,510)
 510      FORMAT(' (SUBR.RDREST) Luminosity file maybe broken.')
          GOTO 990
        ENDIF
      ENDIF
C-- Close file
 600  IF(FNAME.NE.' ') CLOSE(UNIT=UNIT,IOSTAT=IOSTAT)
      IF(IOSTAT.NE.0) GOTO 940
      GOTO 800
C---
 800  IRTN=0
      RETURN
C
 900  IRTN=1000
      WRITE(MSGFL,905) OP(I)
 905  FORMAT(' (SUBR.RDREST) Too many numbers for ',
     %  'operand "',A,'".')
      RETURN
 910  IRTN=1001
      WRITE(MSGFL,915)
 915  FORMAT(' (SUBR.RDREST) Invalid file name.')
      RETURN
 920  IRTN=1002
      RETURN
 940  IRTN=1004
      WRITE(MSGFL,945) UNIT,IOSTAT,FNAME(1:NCFN)
 945  FORMAT(' (SUBR.RDREST) Close file failed. Unit=',I2,/,
     %   '   IOSTAT=',I5,'  File name=',A)
      RETURN
 960  IRTN=1006
      IF(MSGLVL.GE.0) WRITE(MSGFL,965) ERR
 965  FORMAT(' (SUBR.RDREST) Invalid expression for file name.',/,
     %  3X,A)
      RETURN
 990  IRTN=1009
      RETURN
      END

	SUBROUTINE RSTORVAR(UNIT,IRTN)
	USE ARRAYMOD
	IMPLICIT NONE
	INCLUDE 'include/ctrlcm.h'
	INTEGER UNIT,IRTN
	INTEGER, PARAMETER:: MRANK=10
	INTEGER NN,NA,NV,I,NC,IRTN1,TYP,RANK,NC1,N,J,L,ID,ISTAT
	INTEGER DIM(3,MRANK),IND(MRANK)
	REAL*8 VAL
	CHARACTER*16 NAM
	CHARACTER(4) HEAD
	CHARACTER(80) ERR
	CHARACTER*9 MSG(0:1)/'(new)','(revised)'/
	CHARACTER(MGSTR2), POINTER::  CHBUF

	IRTN=0
	ALLOCATE(CHBUF,STAT=ISTAT)
	IF(ISTAT.NE.0) GOTO 950
	IF(MSGLVL.GE.1) WRITE(MSGFL,90)
  90  FORMAT(' --- RESTORE ---')
	NA=0
	NV=0
 100  READ(UNIT,120,END=500) HEAD,NN
 120  FORMAT(A4,I10)
	IF(NN.LE.0) GOTO 900
	IF(HEAD.EQ.'++++') THEN
	  DO I=1,NN
          READ(UNIT,200,END=900) NAM,VAL
 200      FORMAT(A16,1PD25.17)
          CALL SPCOFF(NAM,2,NC)
          IF(NC.LE.0) CYCLE
          CALL SET00(NAM(1:NC),VAL,'',-1,IRTN1)
          IF(IRTN1.GE.2.AND.IRTN1.LE.1000) CYCLE
          IF(IRTN1.GE.1001) THEN
            CALL SET00(NAM(1:NC),VAL,'',MSGLVL,IRTN1)
            GOTO 920
          ENDIF
	    NV=NV+1
          IF(MSGLVL.GE.1) THEN
            WRITE(MSGFL,220) NAM,VAL,MSG(IRTN1)
 220        FORMAT(1X,A,1X,1PD25.17,1X,A)
          ENDIF
	  ENDDO
	ELSEIF(HEAD.EQ.'----') THEN
	  DO I=1,NN
	    READ(UNIT,300,END=900) NAM,NC,TYP,RANK
 300      FORMAT(A16,3I10)
          IF(TYP.LE.0.OR.TYP.GE.3.OR.RANK.LT.0.OR.NC.GT.16) GOTO 900
          IF(RANK.GT.MRANK) GOTO 930
	    N=1
	    IF(RANK.GE.1) THEN
	      DO J=1,RANK
	        READ(UNIT,320,END=900) (DIM(L,J),L=1,3)
 320          FORMAT(I10)
              N=N*DIM(3,J)
	        IF(DIM(1,J).GT.DIM(2,J).OR.
     %             DIM(3,J).NE.(DIM(2,J)-DIM(1,J)+1)) GOTO 900
            ENDDO
	    ENDIF
	    CALL EVDEFARR(NAM(1:NC),0,RANK,3,DIM,0D0,ID,IRTN1,ERR)
	    IF(ERR.NE.' ') GOTO 940
	    IF(MSGLVL.GE.1) THEN
	      CALL ARRAYSTR(NAM(1:NC),RANK,2,3,DIM,CHBUF,NC1)
	      WRITE(MSGFL,330) CHBUF(1:NC1),MSG(IRTN1)
 330        FORMAT(1X,A,T50,A)
          ENDIF
	    DO J=1,N
	      CALL ARRN2IND(IND,RANK,DIM,J,IRTN1)
	      IF(TYP.EQ.1) THEN
		      READ(UNIT,340,END=900) VAL
 340          FORMAT(1PD25.17)
	        CALL EVARRSET(ID,RANK,IND,VAL,ERR)
	      ELSE
	        READ(UNIT,360,END=900) NC1,CHBUF
 360          FORMAT(I10,A)
              CALL EVARRCSET(ID,RANK,IND,CHBUF(1:NC1),ERR)
	      ENDIF
	      IF(ERR.NE.' ') GOTO 960
	    ENDDO
	    NA=NA+1
	  ENDDO
	ELSE
	  GOTO 900
	ENDIF
	GOTO 100

 500  IF(MSGLVL.GE.1) THEN
        IF(NV.GT.0) WRITE(MSGFL,520) NV,'variables'
        IF(NA.GT.0) WRITE(MSGFL,520) NA,'arrays'
 520    FORMAT('   RESTORE:',I4,' ',A,' are added/revised.')
	ENDIF
	IRTN=0
      GOTO 1000

 900  IRTN=1000
	WRITE(MSGFL,905)
 905  FORMAT(' (SUBR.RDREST) RESTORE failed. File is broken.')
	GOTO 1000
 920  IRTN=1002
	GOTO 1000
 930  IRTN=1003
	WRITE(MSGFL,935)
 935  FORMAT(' (SUBR.RDREST) RESTORE failed. Max array rank exceeded.')
	GOTO 1000
 940  IRTN=1004
	WRITE(MSGFL,945) ERR
 945  FORMAT(' (SUBR.RDREST) RESTORE failed in array def.',/,3X,A)
	GOTO 1000
 950  IRTN=1005
	WRITE(MSGFL,955)
 955  FORMAT(' (SUBR.RDREST) RESTORE failed in memory allocation.')
	GOTO 1000
 960  IRTN=1006
	WRITE(MSGFL,965) ERR
 965  FORMAT(' (SUBR.RDREST) RESTORE failed in array fill.',/,3X,A)
	  write(MSGFL,999) rank,(dim(1,i),dim(2,i),i=1,rank)
999     format(' RANK=',I1,' dim= ',5('(',i2,':',i2,')'))
        write(MSGFL,998) (ind(i),i=1,rank)
998     format('      ',1x,' ind= ',5i4)
        write(MSGFL,'(5x,A)') CHBUF(1:NC1)
	GOTO 1000
1000  DEALLOCATE(CHBUF,STAT=ISTAT)
	RETURN
	END