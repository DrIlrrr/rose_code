	SUBROUTINE DEFFILES(IRTN)
C  Define files for Windows envioronment
C  Assume RDFL, OUTFL and TDFL will not be changed after calling subr.INITLZ.
C  Also assume OUTFL2=OUTFL and MSGFL=6.
c	USE MSFLIB
	USE DFLIB
	IMPLICIT NONE
	INTEGER IRTN
	INCLUDE '../include/ctrlcm.h'
	INCLUDE '../include/headcm.h'
	INTEGER status
      TYPE (QWINFO) winfo
	TYPE (windowconfig) wc
	LOGICAL lstatus
	INTEGER NARG,IOSTAT
	CHARACTER*260 NAM,NAM2,INFILE,OUTFILE,TDFILE,PATH,EXT
	INTEGER NCNAME,NCPATH,NCEXT,NCINFL,NCOUTFL,NCTDFL
	INTEGER*2 II
	INTEGER NBUFLINE
C
      OSID=1
C  Readin inifile
	CALL READINI(NBUFLINE,MSGDEST)
C  Maximize the main frame window
C	status = GETWSIZEQQ(QWIN$FRAMEWINDOW, QWIN$SIZEMAX, winfo)
C	winfo.TYPE = QWIN$SET
C	status = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	status = ABOUTBOXQQ(VERSION//' for Windows '//VERSIONDATE)
	
	MSGFL=199
	IF(MSGDEST.EQ.1) THEN
	  MSGFL0=199
	ELSEIF(MSGDEST.GE.2) THEN
	  MSGFL0=6
	ENDIF
      OPEN(MSGFL0,FILE='CON',TITLE='Cain Console',
     %    CARRIAGECONTROL='LIST')
c	II=INITIALIZEFONTS( )
c	ii=SETFONT ('b''FixedSys''h18w8')
      lstatus=GETWINDOWCONFIG(wc)
	wc.numxpixels  = 1024
	wc.numypixels  = -1
	wc.numtextcols = 128
      wc.numtextrows = NBUFLINE
	wc.numcolors   = 16
	wc.title       = "CainW"C
      wc.fontsize    = #00080010
c	wc.fontsize = #0006000C
	lstatus = SETWINDOWCONFIG(wc)
	if (.NOT.lstatus) lstatus = SETWINDOWCONFIG(wc)
C	print *,' numcolors=',wc.numcolors
C  Maximize the console window
	winfo.TYPE = QWIN$MAX
	status = SETWSIZEQQ(MSGFL, winfo)

	NARG=NARGS()
C	print *,' NARG=',narg
	IF(NARG.LE.1) THEN
	  OPEN(RDFL,FILE=' ',ACTION='READ',
     %	   IOSTAT=IOSTAT,STATUS='OLD')
	  IF(IOSTAT.NE.0) THEN
	    IRTN=100
	    WRITE(6,100)
100       FORMAT(' Error: No input file specified.')
          RETURN
	  ENDIF
	  INQUIRE(RDFL,NAME=NAM)
	  NCINFL=FULLPATHQQ(NAM,INFILE)
	ELSE
	  II=1
	  CALL GETARG(II,NAM)
	  NCINFL=FULLPATHQQ(NAM,INFILE)
	  IF(NCINFL.EQ.0) THEN
	    IRTN=101
	    WRITE(6,110)
110       FORMAT(' Error: Invalid file name.')
	    RETURN
	  ENDIF
	  OPEN(RDFL,FILE=INFILE(1:NCINFL),ACTION='READ',
     %	   IOSTAT=IOSTAT,STATUS='OLD')
      ENDIF
C	WRITE(6,300) RDFL,INFILE(1:NCINFL)
	IF(IOSTAT.NE.0) THEN
	  IRTN=102
 	  WRITE(6,320) INFILE(1:NCINFL)
320     FORMAT(' Error: Input file ',A,' open error.')
	  RETURN
	ENDIF
	CALL SPLITPATH(INFILE,PATH,NCPATH,NAM2,NCNAME,EXT,NCEXT)
	status = CHANGEDIRQQ(PATH(1:NCPATH))
C	INFILE=PATH(1:NCPATH)//NAM2(1:NCNAME)//EXT(1:NCEXT)
C	NCINFL=NCPATH+NCNAME+NCEXT
      OUTFILE=PATH(1:NCPATH)//NAM2(1:NCNAME)//'.dat'
	NCOUTFL=NCPATH+NCNAME+4
	TDFILE=PATH(1:NCPATH)//NAM2(1:NCNAME)//'.tdr'
	NCTDFL=NCPATH+NCNAME+4

300	FORMAT(' (SUBR.CMDLIN) Open UNIT=',I2,'  FILE=',A)

C  Open output files
	IF(MSGDEST.GE.2) THEN
	  OPEN(MSGFL,FILE=PATH(1:NCPATH)//NAM2(1:NCNAME)//'.txt',
     %			IOSTAT=IOSTAT,STATUS='REPLACE')
C        In the case MSGDEST=2 (write(6) goto file), the error
C        messages before the above OPEN goes to console
        CALL FILEECHO(0)
      ENDIF
	OPEN(OUTFL,FILE=OUTFILE(1:NCOUTFL),ACTION='WRITE',
     %			IOSTAT=IOSTAT,STATUS='REPLACE')
C	WRITE(6,300) OUTFL,OUTFILE(1:NCOUTFL)
	IF(IOSTAT.NE.0) THEN
	  IRTN=103
	  WRITE(MSGFL,330) OUTFILE(1:NCOUTFL)
330     FORMAT(' Error: Output file ',A,' open error.')
	  RETURN
	ENDIF
	OPEN(TDFL,FILE=TDFILE(1:NCTDFL),ACTION='WRITE',
     %			IOSTAT=IOSTAT,STATUS='REPLACE')
C	WRITE(MSGFL,300) TDFL,TDFILE(1:NCTDFL)
	IF(IOSTAT.NE.0) THEN
	  IRTN=104
	  WRITE(MSGFL,340) TDFILE(1:NCTDFL)
340     FORMAT(' Error: TopDrawer file ',A,' open error.')
	  RETURN
	ENDIF
	IRTN=0
	RETURN
	END

	SUBROUTINE FILEECHO(K)
C       K=0: initialize
C         1: copy the appended part of MSGFL into MSGFL0
	USE DFPORT
	IMPLICIT NONE
	INTEGER K
	INCLUDE '../include/ctrlcm.h'
	CHARACTER(1024) FN
	INTEGER L,RESUL,NCH
	CHARACTER(2048) LINE
	
	IF(K.EQ.0) THEN
	  MSGCOUNT=FTELL(MSGFL)
	  RETURN
	ENDIF
	L=FTELL(MSGFL)
	IF(L.EQ.MSGCOUNT) RETURN
	INQUIRE(MSGFL,NAME=FN)
	CLOSE(MSGFL)
	OPEN(MSGFL,FILE=FN)
	RESUL=FSEEK(MSGFL,MSGCOUNT,0)
200	READ(MSGFL,220,END=300) NCH,LINE(1:NCH)
220   FORMAT(Q,A)
	WRITE(MSGFL0,240) LINE(1:NCH)
240   FORMAT(A)
	GOTO 200
300	MSGCOUNT=L
	RETURN
	END

	SUBROUTINE WINFN2UNIX(FN)
	IMPLICIT NONE
	CHARACTER(*) FN
	INTEGER N,I
	N=LEN(FN)
	DO I=1,N
	  IF(FN(I:I).EQ.'\') FN(I:I)='/'
	ENDDO
	RETURN
	END