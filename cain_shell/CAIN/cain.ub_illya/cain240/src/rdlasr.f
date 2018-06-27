      SUBROUTINE RDLASR(LN,LINE,NCHAR,IRTN)
	USE FLCHTYP
	USE READMOD
	USE ARRAYMOD
	USE LASRDATA
      IMPLICIT NONE
      INTEGER LN(2,3),NCHAR(*),IRTN
      CHARACTER*(*) LINE(*)
      INCLUDE 'include/ctrlcm.h'
	INCLUDE 'include/cnstcm.h'
C      INCLUDE 'include/readcm.h'
      INCLUDE 'include/lasrcm.h'
C
      INTEGER MOP
      PARAMETER (MOP=28)
      CHARACTER*12 OP(MOP)/'RIGHT','LEFT','WAVELENGTH',
     %  'POWERDENSITY','TXYS','E1','E3',
     %  'RAYLEIGH','SIGT','GCUTT','TTOT','TEDGE','GCUT',
     %  'STOKES','TDL','ZMAX',
     %  'FILE','SHIFTT','SHIFTZ','PLOTPROFILE','TPROFILE','NPROFILE',
     %  'DONUT','AAXICON','BAXICON','FOCALLENGTH','SIGMA0','RMAX'
C            From AAXICON to RMAX must be in this order
     %   /
      INTEGER NFF(MOP)/0,0,1,
     %    1,4,3,3, 
     %    2,1,1,1,1,1, 
     %    3,2,1, 
     %    1,1,1,0,2,1,
     %    0,1,1,1,1,1/
      INTEGER IDRL,IDWL,IDPOW,IDTXYS,IDE1,IDE3,IDSTK,IDRAY,IDGAUT,
     %    IDTRPZ,IDGCUT,IDTDL,IDZMAX,IDFILE,IDSHIFTT,IDTPROF,IDNPROF,
     %    IDAAX
C
      INTEGER LR(2),NC,K,J,NF,I,N,IFL,NCFN,NPROF
	TYPE(FLCHTYPE) FC
	REAL*8 POWMAX,POWMAX1,BETA1,TPROF(2),SHIFTTZ(2)
      CHARACTER*256 ERR
	CHARACTER*512 FILENAME
	LOGICAL USEFILE,LTDEF(2),LPLOT,DONUT
      REAL*8 PI/3.14159 26535 89793 238D0/
C
      IRTN=0
      IF(LN(1,2).EQ.0) RETURN
      IF(NLSR.GE.MLSR) GOTO 980
      CALL CMDBLK('LASER',MOP,OP,0,NFF,MBL,NBL,KOP,KBL,REL,
     %    LNKW,LNBL,LN,LINE,NCHAR,IRTN)
      IF(IRTN.GT.1) GOTO 990
      IF(IRTN.EQ.1.OR.NBL.LE.0) RETURN
      NLSR=NLSR+1
	NLSRFL(NLSR)=0
      J=1
      DO 120 I=1,MOP
        ID(I)=J
	  IF(OP(I).EQ.'RIGHT') IDRL=ID(I)
        IF(OP(I).EQ.'WAVELENGTH') IDWL=ID(I)
        IF(OP(I).EQ.'POWERDENSITY') IDPOW=ID(I)
        IF(OP(I).EQ.'TXYS') IDTXYS=ID(I)
        IF(OP(I).EQ.'E1') IDE1=ID(I)
        IF(OP(I).EQ.'E3') IDE3=ID(I)
        IF(OP(I).EQ.'STOKES') IDSTK=ID(I)
        IF(OP(I).EQ.'RAYLEIGH') IDRAY=ID(I)
        IF(OP(I).EQ.'SIGT') IDGAUT=ID(I)
        IF(OP(I).EQ.'TTOT') IDTRPZ=ID(I)
        IF(OP(I).EQ.'GCUT') IDGCUT=ID(I)
	  IF(OP(I).EQ.'TDL') IDTDL=ID(I)
	  IF(OP(I).EQ.'ZMAX') IDZMAX=ID(I)
	  IF(OP(I).EQ.'FILE') IDFILE=ID(I)
	  IF(OP(I).EQ.'SHIFTT') IDSHIFTT=ID(I)
	  IF(OP(I).EQ.'TPROFILE') IDTPROF=ID(I)
	  IF(OP(I).EQ.'NPROFILE') IDNPROF=ID(I)
	  IF(OP(I).EQ.'AAXICON') IDAAX=ID(I)
        J=J+MAX(0,NFF(I))
 120  CONTINUE
      DO 180 I=1,J-1
	  PAR(I)%L=1
        PAR(I)%X=0
 180  CONTINUE
	PAR(IDPOW)=UNDEF
	PAR(IDTPROF)=UNDEF
	NGSTRRD=0

      LR(1)=0
      LR(2)=0
	FILENAME=' '
	NCFN=1
	IFL=0
	LPLOT=.FALSE.
	DONUT=.FALSE.
      DO 300 J=1,NBL
        I=KBL(J)
        IF(NFF(I).GE.1) THEN
          CALL BLKREC(LNBL(1,1,J),LINE,NCHAR,TEXT,NC,IRTN,MSGFL)
          IF(IRTN.NE.0) GOTO 990
          IF(OP(I).EQ.'FILE') THEN
	      CALL EVAL0(TEXT(1:NC),FC,ERR)
	      IF(ERR.NE.' ') GOTO 960
	      IF(FC%L.EQ.1) THEN
	        IFL=NINT(FC%X)
	      ELSE
	        FILENAME=GSTR2(EVALLAST)(FC%C(1):FC%C(2))
	        NCFN=FC%C(2)-FC%C(1)+1
            ENDIF
	      GOTO 300
          ENDIF
          DO 220 K=1,NFF(I)
	      FFF(K)=UNDEF
 220      CONTINUE
          CALL BLKRHS(TEXT(1:NC),FFF,MFFF,NF,GSTRRD,NGSTRRD,IRTN)
          IF(IRTN.NE.0) GOTO 990
          IF(NF.GT.NFF(I)) GOTO 940
          IF(NF.GE.1) THEN
            DO 240 K=1,NF
              IF(FFF(K).NE.UNDEF) PAR(ID(I)+K-1)=FFF(K)
 240        CONTINUE
          ENDIF
        ELSE
	    IF(I.GE.IDRL.AND.I.LE.IDRL+1) THEN
            LR(I-IDRL+1)=1
	    ELSEIF(OP(I).EQ.'PLOTPROFILE') THEN
	      LPLOT=.TRUE.
	    ELSEIF(OP(I).EQ.'DONUT') THEN
	      DONUT=.TRUE.
	    ENDIF
        ENDIF
 300  CONTINUE

C  General parameters
	USEFILE=IFL.NE.0.OR.FILENAME.NE.' '
	IF(LR(1).EQ.1.AND.LR(2).EQ.1) GOTO 950
	LRLSR(NLSR)=0
      IF(LR(1).EQ.1) LRLSR(NLSR)=1
      IF(LR(2).EQ.1) LRLSR(NLSR)=2
      WLLSR(NLSR)=PAR(IDWL)%X
	PLSR(NLSR)=PAR(IDPOW)%X
	IF(PAR(IDPOW).EQ.UNDEF.AND.USEFILE) PLSR(NLSR)=1
      DO 420 I=0,3
        TXYSLS(I,NLSR)=PAR(IDTXYS+I)%X
 420  CONTINUE
      DO 440 I=1,3
        EVLSR(I,1,NLSR)=PAR(IDE1+I-1)%X
        EVLSR(I,3,NLSR)=PAR(IDE3+I-1)%X
        STKSLS(I,NLSR)=PAR(IDSTK+I-1)%X
 440  CONTINUE

C  Spatial structure parameters
	ZETAMAX(NLSR)=PAR(IDZMAX)%X
      DO I=1,2
	  TDL(I,NLSR)=MAX(1D0,PAR(IDTDL+I-1)%X)
	ENDDO
	IF(DONUT) THEN
	  LSPAR(NLSR)=DONUT_SHAPE
	  DO 450 I=1,5
	    SPAR(I,NLSR)=PAR(IDAAX+I-1)%X
 450    CONTINUE
	ELSE
	  LSPAR(NLSR)=GAUSSIAN
        DO 460 I=1,2
          SPAR(I,NLSR)=PAR(IDRAY+I-1)%X
 460    CONTINUE
        SPAR(5,NLSR)=PAR(IDGCUT)%X
	  IF(SPAR(5,NLSR).LE.0) SPAR(5,NLSR)=3.5D0
	ENDIF
C  Time structure parameters
      IF(PAR(IDGAUT)%X.NE.0) THEN
        LTPAR(NLSR)=GAUSSIAN
        TPAR(1,NLSR)=PAR(IDGAUT)%X
        TPAR(2,NLSR)=PAR(IDGAUT+1)%X
      ELSEIF(PAR(IDTRPZ)%X.NE.0) THEN
        LTPAR(NLSR)=TRAPEZOIDAL
        TPAR(1,NLSR)=PAR(IDTRPZ)%X
        TPAR(2,NLSR)=PAR(IDTRPZ+1)%X
      ELSE
        LTPAR(NLSR)=0
      ENDIF
C       initialize Lorentz transformation
      LTRLSR(NLSR)=0
      CALL UNIT8(TRLSR(0,0,NLSR),4,4)
      CALL UNIT8(TRILSR(0,0,NLSR),4,4)
      CALL CLEAR8(TR0LSR(0,NLSR),4)

	IF(USEFILE) THEN
	  IF(FILENAME.NE.' ') THEN
	    IFL=98
	    CALL OPENFL(IFL,FILENAME(1:NCFN),'OLD',0,NCFN,IRTN)
	    IF(IRTN.NE.0) GOTO 920
	  ELSE
	    NCFN=1
	  ENDIF
	  SHIFTTZ(1)=PAR(IDSHIFTT)%X
	  SHIFTTZ(2)=PAR(IDSHIFTT+1)%X
	  CALL LSRRDFL(NLSR,IFL,LTDEF,SHIFTTZ,POWMAX1,
     %     MSGLVL,MSGFL,MSGDEST,IRTN)
C         ZETAMAX needed only for file with ORDER=XY
	  IF(FILENAME.NE.' ') CLOSE(IFL)
	  IF(IRTN.NE.0) GOTO 990
	  IF(LTDEF(1)) LTPAR(NLSR)=DEFINED_BY_FILE
	  IF(LTDEF(2)) LSPAR(NLSR)=DEFINED_BY_FILE
	  POWMAX=POWMAX1*PLSR(NLSR)
	ELSE
	  LTDEF(1)=.FALSE.
	  LTDEF(2)=.FALSE.
	  NCFN=1
	  POWMAX=PLSR(NLSR)
	ENDIF
C  Check general parameters
      ERR=' '
      N=0
      IF(WLLSR(NLSR).LE.0) CALL ERMSG1(ERR,N,'WAVEL')
      IF(PLSR(NLSR).LE.0) CALL ERMSG1(ERR,N,'POWERD')
      IF(EVLSR(1,3,NLSR).EQ.0.AND.EVLSR(2,3,NLSR).EQ.0.AND.
     %  EVLSR(3,3,NLSR).EQ.0) CALL ERMSG1(ERR,N,'E3')
      IF(EVLSR(1,1,NLSR).EQ.0.AND.EVLSR(2,1,NLSR).EQ.0.AND.
     %  EVLSR(3,1,NLSR).EQ.0) CALL ERMSG1(ERR,N,'E1')
C  Check spatioal params
	IF(LSPAR(NLSR).EQ.GAUSSIAN) THEN
	  IF(SPAR(1,NLSR).LE.0) THEN
          CALL ERMSG1(ERR,N,'RAYLEIGH')
	  ELSEIF(SPAR(2,NLSR).LE.0) THEN
	    SPAR(2,NLSR)=SPAR(1,NLSR)
	  ENDIF
	ELSEIF(LSPAR(NLSR).EQ.DONUT_SHAPE) THEN
	  IF(SPAR(2,NLSR).LE.0) THEN
          CALL ERMSG1(ERR,N,'BAXICON')
	  ELSEIF(SPAR(1,NLSR).LE.SPAR(2,NLSR)) THEN
	    CALL ERMSG1(ERR,N,'AAXICON')
	  ELSEIF(SPAR(3,NLSR).LE.0) THEN
	    CALL ERMSG1(ERR,N,'FOCALLENGTH')
	  ELSEIF(SPAR(4,NLSR).LE.0) THEN
	    CALL ERMSG1(ERR,N,'SIGMA0')
	  ELSEIF(SPAR(5,NLSR).LE.0) THEN
	    CALL ERMSG1(ERR,N,'RMAX')
	  ELSEIF(ZETAMAX(NLSR).LE.0) THEN
	    CALL ERMSG1(ERR,N,'ZMAX')
	  ENDIF
	ENDIF
C  Check time params
	IF(LTPAR(NLSR).EQ.0) THEN
        CALL ERMSG1(ERR,N,'SIGT/TTOT')
      ELSEIF(LTPAR(NLSR).EQ.GAUSSIAN) THEN
        IF(TPAR(2,NLSR).EQ.0) TPAR(2,NLSR)=3.5
      ELSEIF(LTPAR(NLSR).EQ.TRAPEZOIDAL) THEN
        IF(2*TPAR(2,NLSR).GT.TPAR(1,NLSR))
     %    CALL ERMSG1(ERR,N,'TTOT<2*TEDGE')
      ENDIF
	IF(ERR.NE.' ') GOTO 960

      CALL ORTHVC(EVLSR(1,3,NLSR),EVLSR(1,1,NLSR),
     %     EVLSR(1,2,NLSR),IRTN)
      IF(IRTN.NE.0) GOTO 970
      WLLSR(NLSR)=WLLSR(NLSR)/(2*PI)
C        note: WLLSR= wavelength/(2*pi)
	OMGLSR(NLSR)=HBARC/WLLSR(NLSR)
C  Define laser range
C       LSRRNG(*,k,*) is the range used for plotting the profile and
C       evaluating the flush energy. Not used for tracking.
C          k=0:  t-z,  k=1: x,  k=2: y,  k=3: z.
	IF(LTPAR(NLSR).EQ.GAUSSIAN) THEN
	  LSRRNG(1,0,NLSR)=-TPAR(1,NLSR)*TPAR(2,NLSR)
	  LSRRNG(2,0,NLSR)=+TPAR(1,NLSR)*TPAR(2,NLSR)
      ELSEIF(LTPAR(NLSR).EQ.TRAPEZOIDAL) THEN
	  LSRRNG(1,0,NLSR)=-TPAR(1,NLSR)/2
	  LSRRNG(2,0,NLSR)=+TPAR(1,NLSR)/2
      ENDIF
	IF(.NOT.LTDEF(2)) THEN
	  LSRRNG(1,3,NLSR)=LSRRNG(1,0,NLSR)
	  LSRRNG(2,3,NLSR)=LSRRNG(2,0,NLSR)
	ENDIF
	IF(LSPAR(NLSR).EQ.GAUSSIAN) THEN
	  DO 520 I=1,2
	    BETA1=SPAR(I,NLSR)
     %       +MAX(LSRRNG(1,3,NLSR)**2,LSRRNG(2,3,NLSR)**2)/SPAR(I,NLSR)
	    LSRRNG(2,I,NLSR)=SPAR(5,NLSR)
     %          *SQRT(SPAR(I+2,NLSR)*0.5D0*WLLSR(NLSR)*BETA1)
	    LSRRNG(1,I,NLSR)=-LSRRNG(2,I,NLSR)
520     CONTINUE
	ELSEIF(LSPAR(NLSR).EQ.DONUT_SHAPE) THEN
        LSRRNG(1,1,NLSR)=-SPAR(5,NLSR)
	  LSRRNG(2,1,NLSR)=+SPAR(5,NLSR)
	  LSRRNG(1,2,NLSR)=-SPAR(5,NLSR)
	  LSRRNG(2,2,NLSR)=+SPAR(5,NLSR)
	  LSRRNG(1,3,NLSR)=-ZETAMAX(NLSR)
	  LSRRNG(2,3,NLSR)=+ZETAMAX(NLSR)
	ENDIF
      
	IF(LSPAR(NLSR).EQ.DONUT_SHAPE) THEN
        CALL DONUTTABLE(NLSR,POWMAX1,MSGLVL,MSGFL,IRTN)
	  IF(IRTN.NE.0) GOTO 990
	  POWMAX=POWMAX1*POWMAX
	ENDIF

      IF(MSGLVL.GE.1)
     %  CALL RDLASR1(NLSR,LRLSR(NLSR),WLLSR(NLSR),
     %   OMGLSR(NLSR),PLSR(NLSR),
     %   TXYSLS(0,NLSR),EVLSR(1,1,NLSR),STKSLS(1,NLSR),
     %   LTPAR(NLSR),TPAR(1,NLSR),LSPAR(NLSR),SPAR(1,NLSR),
     %   TDL(1,NLSR),ZETAMAX(NLSR),
     %   POWMAX,LSRRNG(1,0,NLSR),IFL,FILENAME(1:NCFN),MSGFL,
     %   GAUSSIAN,TRAPEZOIDAL,DONUT_SHAPE,DEFINED_BY_FILE)
	IF(LPLOT) THEN
	  TPROF(1)=PAR(IDTPROF)%X
	  TPROF(2)=PAR(IDTPROF)%X
	  IF(PAR(IDNPROF)%X.EQ.UNDEF%X) THEN
	    IF(TPROF(1).EQ.UNDEF%X) THEN
		    TPROF(1)=0.5D0*(LSRRNG(1,0,NLSR)+LSRRNG(2,0,NLSR)
     %                  +LSRRNG(1,3,NLSR)+LSRRNG(2,3,NLSR))
	      TPROF(2)=TPROF(1)
	      NPROF=1
	    ELSEIF(TPROF(2).EQ.UNDEF%X) THEN
	      TPROF(2)=TPROF(1)
	      NPROF=1
	    ELSE
	      NPROF=2
	    ENDIF
	  ELSE
	    NPROF=NINT(PAR(IDNPROF)%X)
	    IF(TPROF(1).EQ.UNDEF%X) THEN
	      IF(NPROF.EQ.1) THEN
		    TPROF(1)=0.5D0*(LSRRNG(1,0,NLSR)+LSRRNG(2,0,NLSR)
     %                  +LSRRNG(1,3,NLSR)+LSRRNG(2,3,NLSR))
	        TPROF(2)=TPROF(1)
	      ELSE
	        TPROF(1)=LSRRNG(1,3,NLSR)
     %            +0.5D0*(LSRRNG(1,0,NLSR)+LSRRNG(2,0,NLSR))
	        TPROF(2)=LSRRNG(2,3,NLSR)
     %            +0.5D0*(LSRRNG(1,0,NLSR)+LSRRNG(2,0,NLSR))
	      ENDIF
	    ELSEIF(TPROF(2).EQ.UNDEF%X) THEN
	      TPROF(2)=TPROF(1)
	      NPROF=1
	    ELSEIF(TPROF(1).EQ.TPROF(2)) THEN
	      NPROF=1
	    ENDIF
	  ENDIF
	  CALL PLTLASER1(NLSR,NPROF,TPROF,LSRRNG(1,0,NLSR),POWMAX)
	  CALL PLTLASER2(NLSR,NPROF,TPROF,LSRRNG(1,0,NLSR),POWMAX)
	  CALL PLTLASER3(NLSR,NPROF,TPROF,LSRRNG(1,0,NLSR))
	ENDIF
      IRTN=0
      RETURN
 910  IRTN=1001
      WRITE(MSGFL,915) ERR(1:80)
 915  FORMAT(' (SUBR.RDLASR) Invalid expression for file name.',/,
     %   3X,A)
	GOTO 995
 920  IRTN=1002
      WRITE(MSGFL,925) FILENAME(1:NCFN)
 925  FORMAT(' (SUBR.RDLASR) Openfile "',A,'" failed.')
      GOTO 995
 940  IRTN=1004
      WRITE(MSGFL,945) OP(I)
C              corrected Mar.18.2002. used to be OP(J)
 945  FORMAT(' (SUBR.RDLASR) Too many numbers for ',
     %  'operand "',A,'".')
      GOTO 995
 950  IRTN=1005
      WRITE(MSGFL,955)
 955  FORMAT(' (SUBR.RDLASR) Both RIGHT and LEFT ',
     %  'specified.')
      GOTO 995
 960  IRTN=1006
      WRITE(MSGFL,965) ERR(1:N)
 965  FORMAT(' (SUBR.RDLASR) Following laser parameters ',
     %  'not specified.',/,5X,A)
      GOTO 995
 970  IRTN=1007
      WRITE(MSGFL,975) ERR(1:N)
 975  FORMAT(' (SUBR.RDLASR) E3 and E1 vectors are ',
     %  'parallel to each other.')
      GOTO 995
 980  IRTN=1008
      WRITE(MSGFL,985)
 985  FORMAT(' (SUBR.RDLASR) Too many lasers.')
      GOTO 995
 990  IRTN=1009
      GOTO 995
 995  IF(NLSRFL(NLSR).NE.0) CALL FREELSRDT(NLSR)
	NLSR=NLSR-1
	RETURN
      END

      SUBROUTINE RDLASR1(LSR,LR,WLBAR,OMEGA,P,TXYS,EV,STKS,
     %     LT,TPAR,LS,SPAR,TDL,ZETAMAX,POWMAX,REGION,IFL,FILENAME,MSGFL,
     %     GAUSSIAN,TRAPEZOIDAL,DONUT_SHAPE,DEFINED_BY_FILE)
      IMPLICIT NONE
      INTEGER LSR,LR,LT,LS,IFL,MSGFL,
     %   GAUSSIAN,TRAPEZOIDAL,DONUT_SHAPE,DEFINED_BY_FILE
      REAL*8 WLBAR,OMEGA,P,TXYS(0:3),EV(3,3),
     %     TPAR(3),SPAR(5),TDL(2),ZETAMAX,STKS(3),POWMAX,REGION(2,0:3)
	CHARACTER*(*) FILENAME
      INTEGER I,J
	LOGICAL LFILE
	REAL*8 XIMAX,SIGXY0(2),EFLUSH,EFLUSH1(3),ZFLUSH(3)
	REAL*8 DERFC
      CHARACTER*5 LLRR(2)/'Right','Left'/
	CHARACTER*50 BLANCK/' '/
      REAL*8 PI/3.14159 26535 89793 238D0/
      INCLUDE 'include/cnstcm.h'
C
      XIMAX=WLBAR/MASS(2)*SQRT(4*PI*1D-7*CVEL*POWMAX)
	LFILE=LT.EQ.DEFINED_BY_FILE.OR.LS.EQ.DEFINED_BY_FILE
	IF(LS.EQ.GAUSSIAN) THEN
	  DO 80 I=1,2
	    SIGXY0(I)=SQRT(TDL(I)*WLBAR/2*SPAR(I))
 80	  CONTINUE
	  EFLUSH=P*2*PI*SIGXY0(1)*SIGXY0(2)
	  IF(SPAR(5).LE.10) EFLUSH=EFLUSH*(1-EXP(-SPAR(5)**2/2))
	ELSEIF(LS.EQ.DONUT_SHAPE) THEN
	  EFLUSH=(SPAR(1)**2-SPAR(2)**2)/(2*SPAR(4)**2)
	  IF(EFLUSH.GE.30) THEN
	    EFLUSH=1
	  ELSE
	    EFLUSH=1-EXP(-EFLUSH)
	  ENDIF
	  EFLUSH=2*PI*SPAR(4)**2*P*EFLUSH
	ELSEIF(LS.EQ.DEFINED_BY_FILE) THEN
	  EFLUSH=1
	ENDIF
	IF(LT.EQ.GAUSSIAN) THEN
	  EFLUSH=EFLUSH*SQRT(2*PI)*TPAR(1)/CVEL
	  IF(TPAR(2).LE.10) EFLUSH=EFLUSH*(1-DERFC(TPAR(2)/SQRT(2D0)))
	ELSEIF(LT.EQ.TRAPEZOIDAL) THEN
	  EFLUSH=EFLUSH*(TPAR(1)-TPAR(2))/CVEL
	ELSEIF(LT.EQ.DEFINED_BY_FILE) THEN

	ENDIF
	CALL LSRFLUSH(LSR,REGION,3,ZFLUSH,EFLUSH1)
      WRITE(MSGFL,100)
 100  FORMAT(' +++ Laser parameter defined')
      IF(LR.NE.0) WRITE(MSGFL,110) LLRR(3-LR)
 110  FORMAT('    Interact with ',A,' going beam only')
      WRITE(MSGFL,120) WLBAR*2*PI,OMEGA,POWMAX,XIMAX,
     %  (TXYS(I),I=0,3),((EV(I,J),I=1,3),J=1,3)
 120  FORMAT(
     % T5,'Wavelength',T40,6PF10.5,' micron',/,
     % T5,'Photon energy',T40,0PF10.6,' eV',/,
     % T5,'Peak power density',T40,1PD10.3,' Watt/m**2',/,
     % T5,'Maximum Xi parameter',T40,0PF10.5,/,
     % T5,'Focus (t,x,y,s)',T30,1P4D10.3,' meter',/,
     % T5,'Direction of e1 vector',T40,0P3F10.6,/,
     % T5,'Direction of e2 vector',T40,0P3F10.6,/,
     % T5,'Propagation direction (e3 vector)',T40,0P3F10.6)
	WRITE(MSGFL,130) (ZFLUSH(I),I=1,3),(EFLUSH1(I),I=1,3)
 130  FORMAT(T5,'Flush energy (num.integ. at z=)',
     %   T40,1P3D10.3,' m',/, T40,1P3D10.3,' Joule')
	IF(.NOT.LFILE) WRITE(MSGFL,140) EFLUSH
 140  FORMAT(T5,'Flush energy (analytic)',T40,1PD10.3,' Joule')
      WRITE(MSGFL,150) (STKS(I),I=1,3)
 150  FORMAT(T5,'Stokes parameter',T40,'(',0P3F9.5,')')
	
      IF(LT.EQ.GAUSSIAN) THEN
        WRITE(MSGFL,160) (TPAR(I),I=1,2)
 160    FORMAT(T5,'Time profile of pulse',T40,'Gaussian',/,
     % T5,'R.m.s. pulse length',T40,3PF10.4,' mm',/,
     % T5,'Longitudinal Gaussian tail cutoff',T40,0PF10.3,' sigmas')
      ELSEIF(LT.EQ.TRAPEZOIDAL) THEN
        WRITE(MSGFL,170) (TPAR(I),I=1,2)
 170    FORMAT(T5,'Time profile of pulse',T40,'Trapezoidal',/,
     % T5,'Total pulse length',T40,6PF10.3,' micron',/,
     % T5,'Edge length',T40,6PF10.3,' micron')
	ELSEIF(LT.EQ.DEFINED_BY_FILE) THEN
	  WRITE(MSGFL,180)
 180    FORMAT(T5,'Time profile of pulse',T40,'from file')
	ENDIF
	IF(LS.EQ.GAUSSIAN) THEN
        WRITE(MSGFL,200) (SPAR(I),I=1,2),SPAR(5),(TDL(I),I=1,2),
     %    (SIGXY0(I),I=1,2)
 200    FORMAT(T5,'Spatial profile of pulse',T40,'Gaussian',/,
     %   T5,'Rayleigh length in e1,e2 direction',
     %      T40,3P2F10.4,' mm',/,
     %   T5,'Transverse Gaussian tail cutoff',T40,0PF10.3,' sigmas',/,
     %   T5,'Emittance dilution (x,y)',T40,0P2F10.5,' TDL',/,
     %   T5,'Rms beam size (x,y) at focus',T40,6P2F10.2,' micron')
	ELSEIF(LS.EQ.DONUT_SHAPE) THEN
	  WRITE(MSGFL,210) P,P*2*PI*SPAR(4)**2,(SPAR(I),I=1,5),
     %   ZETAMAX
 210    FORMAT(T5,'Spatial profile of pulse',T40,'Donut-shape',/,
     %   T5,'Input peak power density',T40,1PD10.3,' Watt/m**2',/,
     %   T5,'Input peak power',T40,1PD10.3,' Watt',/,
     %   T5,'Outer radius of axicon mirror',T40,3PF10.4,' mm',/,
     %   T5,'Inner radius of axicon mirror',T40,3PF10.4,' mm',/,
     %   T5,'Focal length',T40,0PF10.5,' m',/,
     %   T5,'Rms input laser radius',T40,3PF10.4,' mm',/,
     %   T5,'Field region: max(r)',T40,6PF10.2,' micron',/,
     %   T5,'              max(abs(z))',T40,3PF10.4,' mm')
      ELSEIF(LS.EQ.DEFINED_BY_FILE) THEN
	  WRITE(MSGFL,220)
 220    FORMAT(T5,'Space profile of pulse',T40,'from file')
      ENDIF
	IF(LT.EQ.DEFINED_BY_FILE.OR.LS.EQ.DEFINED_BY_FILE) THEN
	  IF(FILENAME.EQ.' ') THEN
	    WRITE(MSGFL,240) IFL
 240      FORMAT(T5,'File unit number',T40,I10)
	  ELSE
	    WRITE(MSGFL,250) BLANCK(1:MAX(20,50-LEN(FILENAME)))//FILENAME
 250      FORMAT(T5,'File name',A)
        ENDIF
	ENDIF

      RETURN
      END
