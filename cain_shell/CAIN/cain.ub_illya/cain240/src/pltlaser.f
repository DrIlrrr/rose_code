	SUBROUTINE PLTLASER1(LSR,NPROF,TPROF,REGION,PDMAX)
C  Contour plot
	USE LASRDATA
	IMPLICIT NONE
	INTEGER LSR,NPROF
	REAL*8 TPROF(2),REGION(2,0:3),PDMAX
	INCLUDE 'include/ctrlcm.h'
      INCLUDE 'include/lasrcm.h'
      INCLUDE 'include/topdraw.h'
	INTEGER IPROF,NZ,NX,IZ,IX,ISTAT,I,J
	REAL*8, ALLOCATABLE :: PP(:,:)
	INTEGER IEX(2)
	REAL*8 UVW(0:3),TXYS(0:3),T0,Y0,PD,PD0,EV(3,3),OMG,
     %       ZXMM(2,2),UNIT(2)
	REAL*8 XL/4.0/,XH/12.0/,YL/1.0/,YH/9.0/,X1,Y1
	INTEGER MFV
	PARAMETER (MFV=13)
	INTEGER LFV(MFV)
	REAL*8  FV(MFV)/0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,
     %   0.6,0.7,0.8,0.9,0.95/
	INTEGER LM(MFV)/9,3,9,1,2,3,2,1,2,3,2,3,9/
	INTEGER IC(MFV)/4,6,4,2,4,6,4,2,4,6,4,6,4/
	
	IF(NPROF.LE.0) RETURN
	CALL CPUTIM('PLTLASER',1)
	NZ=100
	NX=100
	ALLOCATE(PP(0:NZ,0:NX),STAT=ISTAT)
	IF(ISTAT.NE.0) THEN
	  WRITE(MSGFL,160)
160     FORMAT(' *** Memory allocation for laser profile plot failed.')
        GOTO 600
	ENDIF
	DO 500 IPROF=1,NPROF
	IF(NPROF.EQ.1) THEN
	  T0=TPROF(1)
	ELSE
	  T0=TPROF(1)+(TPROF(2)-TPROF(1))*(IPROF-1)/(NPROF-1)
	ENDIF
	Y0=0.5D0*(REGION(1,2)+REGION(2,2))
	ZXMM(1,1)=MAX(REGION(1,3),T0-REGION(2,0))
	ZXMM(2,1)=MIN(REGION(2,3),T0-REGION(1,0))
	IF(ZXMM(1,1).GE.ZXMM(2,1)) THEN
	  WRITE(MSGFL,120) T0
120     FORMAT(' (PLTLASER) Laser profile out of laser range for t=',
     %     1PD11.4,' m')
        GOTO 500
	ENDIF
	ZXMM(1,2)=REGION(1,1)
	ZXMM(2,2)=REGION(2,1)

	DO 180 I=1,2
	  IEX(I)=3*(INT(LOG10(MAX(ABS(ZXMM(1,I)),ABS(ZXMM(2,I))))/3D0
     %          +100D0)-100)
	  UNIT(I)=10D0**IEX(I)
180   CONTINUE
	UVW(0)=T0
	UVW(2)=Y0
	DO 300 IX=0,NX
	  UVW(1)=ZXMM(1,2)+(ZXMM(2,2)-ZXMM(1,2))/NX*IX
	  DO 260 IZ=0,NZ
	    UVW(3)=ZXMM(1,1)+(ZXMM(2,1)-ZXMM(1,1))/NZ*IZ
	    DO 200 I=0,3
	      TXYS(I)=UVW(I)
200       CONTINUE
	    CALL LSRCOORD(LSR,TXYS,1)
          CALL LSRGEO(LSR,TXYS,PD,PD0,EV,OMG,0)
	    PP(IZ,IX)=PD/PDMAX
260     CONTINUE
300   CONTINUE

	WRITE(TDFL,320)
320   FORMAT(' NEWFRAME; SET FONT DUPLEX')
      CALL TDHEAD(TDFL)
	DO 360 I=1,MFV
	  LFV(I)=10*IC(I)+LM(I)
360   CONTINUE
	WRITE(TDFL,400) XL,XH,YL,YH,((ZXMM(J,I)/UNIT(I),J=1,2),I=1,2),
     %    (XL+XH)/2,YH+0.5,T0,Y0,(XL+XH)/2,0.4,IEX(1),
     %    XL-1.0,(YL+YH)/2,IEX(2)
400   FORMAT(' SET WINDOW X',0P2F7.3,' Y',2F7.3,/,
     %  ' SET LIMIT X',1P2D12.4,' Y',2D12.4,/,
     %  ' TITLE',0P2F7.3,' SIZE 1.8 CENTER ',2H'',/,
     %  ' MORE ',1H','Laser Profile at T=',1PD9.2,
     %                                'm  H=',1PD9.2,'m',1H',/,
     %  ' CASE ',1H','                 G ',9X,
     %                                '   G ',9X,    ' ',1H',/,
     %  ' TITLE',0P2F7.3,' SIZE 2.4 ',2H'',/,
     %  ' MORE ',1H','Z (102',I3,'3m)',1H',/,
     %  ' CASE ',1H','G    X',3X,'X  ',1H',/,
     %  ' TITLE',0P2F7.3,' SIZE 2.4 ANGLE 90 ',2H'',/,
     %  ' MORE ',1H','X (102',I3,'3m)',1H',/,
     %  ' CASE ',1H','G    X',3X,'X  ',1H',/,
     %  ' PLOT AXIS')
	CALL TDCNT1(PP,NZ,NZ,NX,ZXMM(1,1)/UNIT(1),ZXMM(2,1)/UNIT(1),
     %  ZXMM(1,2)/UNIT(2),ZXMM(2,2)/UNIT(2),
     %  MFV,FV,LFV,TDFL,MSGFL)
	Y1=YH-0.3
	X1=0.3
	WRITE(TDFL,410) X1+1.0,Y1
410   FORMAT(' TITLE',2F7.3,' SIZE 1.5 ',1H','P/P0max1',1H',/,
     %     '         CASE ',             1H','   X   X',1H')
	Y1=Y1-0.4
	WRITE(TDFL,420) X1+0.3,Y1,PDMAX
420   FORMAT(' TITLE',2F7.3,' SIZE 1.5 ',2H'',/,
     %  ' MORE ',1H','P0max1=',1PD8.2,'W/m223',1H',/,
     %  ' CASE ',1H',' X   X ',8X,    '   X X',1H',';PLOT AXIS')
	DO 460 I=1,MFV
	  Y1=Y1-0.4
	  WRITE(TDFL,430) X1,Y1,X1+0.8,Y1,COLOR(IC(I)),PATTRN(LM(I))
430     FORMAT(2(F6.3,F7.3,';'),'SET COLOR ',A,';SET PATTERN ',A,/,
     %    ' JOIN 1 TEXT PATTERNED; PLOT AXIS')
	  WRITE(TDFL,440) X1+1.0,Y1,FV(I)
440     FORMAT(' SET COLOR WHITE',/,
     %   ' TITLE',2F7.3,' SIZE 1.5 ',1H',F5.3,1H',';PLOT AXIS')
460   CONTINUE 
500   CONTINUE	
590	DEALLOCATE(PP,STAT=ISTAT)
600	CALL CPUTIM('PLTLASER',2)
	RETURN
	END
	SUBROUTINE PLTLASER2(LSR,NPROF,TPROF,REGION,PDMAX)
C  Radial profile
	USE LASRDATA
	IMPLICIT NONE
	INTEGER LSR,NPROF
	REAL*8 TPROF(2),REGION(2,0:3),PDMAX
	INCLUDE 'include/ctrlcm.h'
      INCLUDE 'include/lasrcm.h'
      INCLUDE 'include/topdraw.h'
	INTEGER IPROF,NZ,NX,IZ,IX,ISTAT,I,J,LM,IC
	REAL*8, ALLOCATABLE :: PP(:)
	INTEGER IEXX,IEXP
	REAL*8 UVW(0:3),TXYS(0:3),T0,Y0,EV(3,3),OMG,XMIN,XMAX,ZMIN,ZMAX,
     %   UNITX,UNITP,PD0,PDMIN
	REAL*8 XL/4.0/,XH/12.0/,YL/1.0/,YH/9.0/,X1,Y1
	
	IF(NPROF.LE.0) RETURN
	CALL CPUTIM('PLTLASER',1)
	NZ=10
	NX=200
	ALLOCATE(PP(0:NX),STAT=ISTAT)
	IF(ISTAT.NE.0) THEN
	  WRITE(MSGFL,160)
160     FORMAT(' *** Memory allocation for laser profile plot failed.')
        GOTO 600
	ENDIF

	DO 500 IPROF=1,NPROF
	IF(NPROF.EQ.1) THEN
	  T0=TPROF(1)
	ELSE
	  T0=TPROF(1)+(TPROF(2)-TPROF(1))*(IPROF-1)/(NPROF-1)
	ENDIF
	Y0=0.5D0*(REGION(1,2)+REGION(2,2))
	ZMIN=MAX(REGION(1,3),T0-REGION(2,0))
	ZMAX=MIN(REGION(2,3),T0-REGION(1,0))
	IF(ZMIN.GE.ZMAX) THEN
	  WRITE(MSGFL,120) T0
120     FORMAT(' (PLTLASER2) Laser profile out of laser range for t=',
     %     1PD11.4,' m')
        GOTO 500
	ENDIF
	XMIN=REGION(1,1)
	XMAX=REGION(2,1)

	IEXX=3*(INT(LOG10(MAX(ABS(XMIN),ABS(XMAX)))/3D0+100D0)-100)
	UNITX=10D0**IEXX
	IEXP=3*(INT(LOG10(PDMAX)/3D0+100D0)-100)
	UNITP=10D0**IEXP
	PDMIN=PDMAX/1D4

	WRITE(TDFL,200)
200   FORMAT(' NEWFRAME; SET FONT DUPLEX')
      CALL TDHEAD(TDFL)

	WRITE(TDFL,220) XL,XH,YL,YH,XMIN/UNITX,XMAX/UNITX,
     %    PDMIN/UNITP,PDMAX/UNITP,
     %    (XL+XH)/2,YH+0.5,T0,Y0,(XL+XH)/2,0.4,IEXX,
     %    XL-1.0,(YL+YH)/2,IEXP
220   FORMAT(' SET WINDOW X',0P2F7.3,' Y',2F7.3,/,
     %  ' SET LIMIT X',1P2D12.4,' Y',2D12.4,/,
     %  ' SET SCALE Y LOG',/,
     %  ' TITLE',0P2F7.3,' SIZE 1.8 CENTER ',2H'',/,
     %  ' MORE ',1H','Laser Radial Profile at T=',1PD9.2,
     %                                'm  H=',1PD9.2,'m',1H',/,
     %  ' CASE ',1H','                        G ',9X,
     %                                '   G ',9X,    ' ',1H',/,
     %  ' TITLE',0P2F7.3,' SIZE 2.4 CENTER ',2H'',/,
     %  ' MORE ',1H','X (102',I3,'3m)',1H',/,
     %  ' CASE ',1H','G    X',3X,'X  ',1H',/,
     %  ' TITLE',0P2F7.3,' SIZE 1.8 ANGLE 90 CENTER ',2H'',/,
     %  ' MORE ',1H','Power Density (102',I2,'3 W/m223)',1H',/,
     %  ' CASE ',1H','                 X',2X,'X    X X ',1H',/,
     %  ' PLOT AXIS')
	Y1=YH-0.3
	X1=0.3
	WRITE(TDFL,240) X1+1.0,Y1
240   FORMAT(' SET COLOR WHITE',/,
     %   ' TITLE',2F7.3,' SIZE 1.5 ',1H','z (mm)',1H',';PLOT AXIS')
	DO 480 IZ=0,NZ
	  UVW(0)=T0
	  UVW(2)=Y0
	  UVW(3)=ZMIN+(ZMAX-ZMIN)/NZ*IZ
	  DO 300 IX=0,NX
	    UVW(1)=XMIN+(XMAX-XMIN)/NX*IX
	    DO 260 I=0,3
	      TXYS(I)=UVW(I)
260       CONTINUE
	    CALL LSRCOORD(LSR,TXYS,1)
          CALL LSRGEO(LSR,TXYS,PP(IX),PD0,EV,OMG,0)
300     CONTINUE
	  LM=MOD(IZ,MLMODE)+1
	  IC=MOD(IZ,MCOLR)+1
	  WRITE(TDFL,400) ((XMIN+(XMAX-XMIN)/NX*IX)/UNITX,
     %       PP(IX)/UNITP,IX=0,NX)
400     FORMAT(3(1PD11.4,D12.4,';'))
        WRITE(TDFL,420) COLOR(IC),PATTRN(LM)
420     FORMAT(' SET COLOR ',A,';SET PATTERN ',A,';JOIN 1 PATTERNED')
	  Y1=Y1-0.4
	  WRITE(TDFL,440) X1,Y1,X1+0.8,Y1,COLOR(IC),PATTRN(LM)
440     FORMAT(2(F6.3,F7.3,';'),'SET COLOR ',A,';SET PATTERN ',A,/,
     %    ' JOIN 1 TEXT PATTERNED; PLOT AXIS')
	  WRITE(TDFL,460) X1+1.0,Y1,
     %    (ZMIN+(ZMAX-ZMIN)/NZ*IZ)*1D3
460     FORMAT(' SET COLOR WHITE',/,
     %   ' TITLE',2F7.3,' SIZE 1.5 ',1H',F7.3,1H',';PLOT AXIS')
480   CONTINUE 
500   CONTINUE	
	DEALLOCATE(PP,STAT=ISTAT)
600	CALL CPUTIM('PLTLASER',2)
	RETURN
	END

	SUBROUTINE PLTLASER3(LSR,NPROF,TPROF,REGION)
C  Longitudinal profile
	USE LASRDATA
	IMPLICIT NONE
	INTEGER LSR,NPROF
	REAL*8 TPROF(2),REGION(2,0:3)
	INCLUDE 'include/ctrlcm.h'
      INCLUDE 'include/lasrcm.h'
      INCLUDE 'include/topdraw.h'
	INTEGER IT,NZ,IZ,ISTAT,K,LM,IC
	REAL*8, ALLOCATABLE :: FF(:,:),ZTMM(:,:),TAU(:)
	INTEGER IEXZ,IEXP
	REAL*8 TXYS(0:3),ZETA,ZT,ZTMIN,ZTMAX,ZTMM1,
     %   UNITZ,UNITP,FMIN,FMAX,FMAX1,V(3,3),OMG,PD0
	REAL*8 XL/4.0/,XH/12.0/,YL/1.0/,YH/9.0/,X1,Y1
	CHARACTER*23 TTL(2) /'Power Density at X=H=0 ',
     %   'Instantaneous Power '/
      CHARACTER*23 TTLC(2)/'                 G G   ',
     %   '                    '/
      CHARACTER*6 LUNIT(2) /'W/m223','W'/
	CHARACTER*6 LUNITC(2)/'   X X',' '/
	
	IF(NPROF.LE.0) RETURN
	CALL CPUTIM('PLTLASER',1)

	DO 500 K=1,2
	IF(K.EQ.1) THEN
	  NZ=500
	ELSE
	  NZ=50
	ENDIF
	ALLOCATE(FF(0:NZ,NPROF),TAU(NPROF),ZTMM(2,NPROF),STAT=ISTAT)
	IF(ISTAT.NE.0) THEN
	  WRITE(MSGFL,160) NZ,NPROF
160     FORMAT(' (SUBR.PLTLASER3) Memory allocation for laser profile ',
     %     'plot failed.',/,'    NZ=',I4,'  NPROF=',I4)
        GOTO 600
	ENDIF

	FMAX1=0
	ZTMIN=1D60
	ZTMAX=-1D60
	DO IT=1,NPROF
	  IF(NPROF.EQ.1) THEN
	    TAU(IT)=0.5D0*(TPROF(2)+TPROF(1))
	  ELSE
	    TAU(IT)=TPROF(1)+(TPROF(2)-TPROF(1))*(IT-1)/(NPROF-1)
	  ENDIF
	  ZTMM(1,IT)=MAX(REGION(1,3)-TAU(IT),-REGION(2,0))
	  ZTMM(2,IT)=MIN(REGION(2,3)-TAU(IT),-REGION(1,0))
	  ZTMM1=ZTMM(2,IT)-ZTMM(1,IT)
	  ZTMM(1,IT)=ZTMM(1,IT)-0.04*ZTMM1
	  ZTMM(2,IT)=ZTMM(2,IT)+0.04*ZTMM1
	  ZTMIN=MIN(ZTMIN,ZTMM(1,IT))
	  ZTMAX=MAX(ZTMAX,ZTMM(2,IT))
	  DO IZ=0,NZ
	    ZT=ZTMM(1,IT)+(ZTMM(2,IT)-ZTMM(1,IT))*IZ/NZ
	    ZETA=ZT+TAU(IT)
	    IF(K.EQ.1) THEN
	      TXYS(0)=TAU(IT)
	      TXYS(1)=0.5D0*(REGION(2,1)+REGION(1,1))
	      TXYS(2)=0.5D0*(REGION(2,2)+REGION(1,2))
	      TXYS(3)=ZETA
	      CALL LSRCOORD(LSR,TXYS,1)
	      CALL LSRGEO(LSR,TXYS,FF(IZ,IT),PD0,V,OMG,0)
	    ELSE
	      CALL LSRINTEG2(LSR,REGION,TAU(IT),ZETA,FF(IZ,IT))
	    ENDIF
CC           Visual Fortran 5.1 causes errors FF(IZ,IT) is written directly
CC           in the call of LSRINTEG2.  Compiler error? of ALLOCATE?
	    FMAX1=MAX(FMAX1,FF(IZ,IT))
	  ENDDO
	ENDDO
	IF(FMAX1.LE.0) THEN
	  WRITE(MSGFL,180)
180     FORMAT(' (SUBR.PLTLASER3) Laser field 0. Plot ignored.')
        GOTO 500
	ENDIF
	FMAX1=1D-30*FMAX
	FMAX=0
	FMIN=1D60
	DO IT=1,NPROF
	  DO IZ=0,NZ
	    IF(FF(IZ,IT).LT.FMAX1) THEN
			  FF(IZ,IT)=0
	    ELSE
		    FMAX=MAX(FMAX,FF(IZ,IT))
	      FMIN=MIN(FMIN,FF(IZ,IT))
	    ENDIF
	  ENDDO
	ENDDO
	FMAX=1.05*FMAX
	IEXP=3*(INT(LOG10(FMAX)/3D0+100D0)-100)
	FMIN=MAX(FMIN,FMAX/1D4)
	UNITP=10D0**IEXP
	IEXZ=3*(INT(LOG10(MAX(ABS(ZTMIN),ABS(ZTMAX)))/3D0+100D0)-100)
	UNITZ=10D0**IEXZ

	WRITE(TDFL,200)
200   FORMAT(' NEWFRAME; SET FONT DUPLEX')
      CALL TDHEAD(TDFL)

	WRITE(TDFL,220) XL,XH,YL,YH,ZTMIN/UNITZ,ZTMAX/UNITZ,
     %   FMIN/UNITP,FMAX/UNITP,(XL+XH)/2,YH+0.5,(XL+XH)/2,0.4,IEXZ
220   FORMAT(' SET WINDOW X',0P2F7.3,' Y',0P2F7.3,/,
     %  ' SET LIMIT X',1P2D12.4,' Y',1P2D12.4,/,
     %  ' SET SCALE Y LOG',/,
     %  ' TITLE',0P2F7.3,' SIZE 1.8 CENTER ',2H'',/,
     %  ' MORE ',1H','Laser Longitudinal Profile',1H',/,
     %  ' CASE ',1H','                          ',1H',/,
     %  ' TITLE',0P2F7.3,' SIZE 2.4 CENTER ',2H'',/,
     %  ' MORE ',1H','Z-T (102',I3,'3m)',1H',/,
     %  ' CASE ',1H','G G    X',3X,'X  ',1H')
	WRITE(TDFL,230) 
     %    XL-1.0,(YL+YH)/2,TTL(K),IEXP,LUNIT(K),TTLC(K),LUNITC(K)
230   FORMAT(
     %  ' TITLE',0P2F7.3,' SIZE 1.8 ANGLE 90 CENTER ',2H'',/,
     %  ' MORE ',1H',A,'102',I2,'3 ',A,1H',/,
     %  ' CASE ',1H',A,'  X',2X,'X ',A,1H',/,
     %  ' PLOT AXIS')
	  Y1=YH-0.3
	  X1=0.3
	  WRITE(TDFL,240) X1+1.0,Y1
240     FORMAT(' SET COLOR WHITE',/,
     %     ' TITLE',2F7.3,' SIZE 1.5 ',1H','T (mm)',1H',/,
     %     '  CASE ',                  1H','G     ',1H',';PLOT AXIS')
	DO IT=1,NPROF
	  LM=MOD(IT-1,MLMODE)+1
	  IC=MOD(IT-1,MCOLR)+1
	  WRITE(TDFL,400) 
     %     ((ZTMM(1,IT)+(ZTMM(2,IT)-ZTMM(1,IT))*IZ/NZ)/UNITZ,
     %       FF(IZ,IT)/UNITP,IZ=0,NZ)
400     FORMAT(3(1PD11.4,D12.4,';'))
        WRITE(TDFL,420) COLOR(IC),PATTRN(LM)
420     FORMAT(' SET COLOR ',A,';SET PATTERN ',A,';JOIN 1 PATTERNED')
	  Y1=Y1-0.4
	  WRITE(TDFL,440) X1,Y1,X1+0.8,Y1,COLOR(IC),PATTRN(LM)
440     FORMAT(2(F6.3,F7.3,';'),'SET COLOR ',A,';SET PATTERN ',A,/,
     %    ' JOIN 1 TEXT PATTERNED; PLOT AXIS')
	  WRITE(TDFL,460) X1+1.0,Y1,TAU(IT)*1D3
460     FORMAT(' SET COLOR WHITE',/,
     %   ' TITLE',2F7.3,' SIZE 1.5 ',1H',F7.3,1H',';PLOT AXIS')
      ENDDO

500   DEALLOCATE(FF,TAU,ZTMM,STAT=ISTAT)
	
600	CALL CPUTIM('PLTLASER',2)
	RETURN
	END

	SUBROUTINE LSRINTEG2(LSR,REGION,TAU,ZETA,PTOT)
C  2D numerical integration of laser energy on a given plane zeta=const, tau=const.
C  Used only for profile plot. Not used for tracking.
C  Must be called before any Lorentz transformation after laser definition.
C  Output:  PTOT= integrated power  in Watt
	IMPLICIT NONE
	INTEGER LSR
	REAL*8 REGION(2,0:3),TAU,ZETA,PTOT
	INTEGER IX,IY,IX1,IY1,IT1,ISTAT,ND
	REAL*8 TXYS(0:3),PD,V(3,3),OMG,PD0,SUM,F1,DX,DY,FACX,FACY
	INTEGER MM,MAXDIV
	PARAMETER (MM=50,MAXDIV=10)
	REAL*8, ALLOCATABLE :: FF(:,:)
      INCLUDE 'include/ctrlcm.h'
	
	CALL CPUTIM('LSRINTEG2',1)
	PTOT=0
	ALLOCATE(FF(0:MM,0:MM), STAT=ISTAT)
	IF(ISTAT.NE.0) THEN
	  IF(MSGLVL.GE.0) THEN
	    WRITE(MSGFL,100)
100       FORMAT(' (SUBR.LSRINTEG2) Allocation error.')
        ENDIF
	  GOTO 600
	ENDIF
	DX=(REGION(2,1)-REGION(1,1))/MM
	DY=(REGION(2,2)-REGION(1,2))/MM
C   start from a coarse mesh
	SUM=0
	DO 280 IX=0,MM
	  DO 260 IY=0,MM
	    TXYS(0)=TAU
	    TXYS(1)=REGION(1,1)+DX*IX
	    TXYS(2)=REGION(1,2)+DY*IY
	    TXYS(3)=ZETA
	    CALL LSRCOORD(LSR,TXYS,1)
	    CALL LSRGEO(LSR,TXYS,FF(IX,IY),PD0,V,OMG,0)
	    SUM=SUM+FF(IX,IY)
260     CONTINUE
280   CONTINUE
	IF(SUM.LE.0) GOTO 580
	DO 480 IX=0,MM-1
	  DO 460 IY=0,MM-1
	    F1=(FF(IX+1,IY+1)+FF(IX,IY+1)+FF(IX+1,IY)+FF(IX,IY))/4D0
	    ND=MAX(1,MIN(MAXDIV,NINT(SQRT(F1/SUM*MM**2))))
	    IF(ND.GE.2) THEN
	      F1=0
	      DO 420 IX1=0,ND
	        FACX=1D0
	        IF(IX1.EQ.0.OR.IX1.EQ.ND) FACX=0.5D0
	        DO 400 IY1=0,ND
	          TXYS(0)=TAU
	          TXYS(1)=REGION(1,1)+DX*(IX*ND+IX1)/ND
	          TXYS(2)=REGION(1,2)+DY*(IY*ND+IY1)/ND
	          TXYS(3)=ZETA
	          FACY=1D0
	          IF(IY1.EQ.0.OR.IY1.EQ.ND) FACY=0.5D0
	          CALL LSRCOORD(LSR,TXYS,1)
	          CALL LSRGEO(LSR,TXYS,PD,PD0,V,OMG,0)
	          F1=F1+PD*FACX*FACY
400           CONTINUE
420         CONTINUE
	      F1=F1/ND**2
	    ENDIF
	    PTOT=PTOT+F1
460     CONTINUE
480   CONTINUE
	PTOT=PTOT*DX*DY

580   DEALLOCATE(FF,STAT=ISTAT)
600	CALL CPUTIM('LSRINTEG2',2)
	RETURN
	END