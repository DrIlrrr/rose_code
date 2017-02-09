      SUBROUTINE PLBLOPT(FILE,NAM,DSINTP,IRTN)
C  Plot optics of the beamline NAM in file IFL
C  Input
C    NAM         beamline name
C    FILE        output file unit number
C    DSINTP      If>0, interpolation to length DSINTP meter
C  Output
C    IRTN     0   normal
C          1001   Beamline nonexisting
C          1002   optics not ready
C
	USE BEAMLN
	IMPLICIT NONE
	INTEGER FILE,IRTN
	CHARACTER(*) NAM
	REAL(8) DSINTP
      INCLUDE 'include/ctrlcm.h'
	INTEGER II,JJ,I,J,N,NN,IBE,N1,K,MDIV,NDIV,IRTN1,ISTAT
	INTEGER LET(2)
	REAL*8 APERT1(2),BTMX(2),ETMM(2,2),BTMAX,ETMIN,ETMAX,STOT,DS
	LOGICAL INTP
	REAL(8) XWL/1.5/,XWH/11.5/,YWL/3.0/,YWH/9.0/,X1,Y1
	CHARACTER(1) LXY(2)/'x','y'/
	CHARACTER(5) LBE(2)/'SqrtB','  Eta'/
	CHARACTER(5) COLOR(2)/'RED','BLUE'/
	CHARACTER(7) LMOD(2)/'SOLID','DASHES'/
	REAL(8), ALLOCATABLE:: TWISS(:,:)

	CALL STARTBL(NAM,II,APERT1)
	IF(II.LE.0) GOTO 910
	NN=BL(II)%NEXP
	IF(NN.LE.0) GOTO 920
	IF(.NOT.BL(II)%LTWISS) GOTO 920
	DO J=1,2
	  BTMX(J)=BL(II)%TWISS(1+5*(J-1),0)
	  DO I=1,2
	    ETMM(I,J)=BL(II)%TWISS(3+5*(J-1),0)
	  ENDDO
	ENDDO
	DO N=1,NN
	  DO J=1,2
	    BTMX(J)=MAX(BTMX(J),BL(II)%TWISS(1+5*(J-1),N))
	    ETMM(1,J)=MIN(ETMM(1,J),BL(II)%TWISS(3+5*(J-1),N))
	    ETMM(2,J)=MAX(ETMM(2,J),BL(II)%TWISS(3+5*(J-1),N))
	  ENDDO
	ENDDO
	BTMAX=SQRT(MAX(BTMX(1),BTMX(2)))*1.05
	ETMIN=0
	ETMAX=0
	DO J=1,2
	  IF(ETMM(1,J).GE.-1D-4.AND.ETMM(2,J).LE.1D-4) THEN
	    LET(J)=0
	  ELSE
	    LET(J)=1
	    IF(ETMM(1,J).LE.-1D-4) ETMIN=MIN(ETMIN,ETMM(1,J))
	    IF(ETMM(2,J).GE.1D-4) ETMAX=MAX(ETMAX,ETMM(2,J))
	  ENDIF
	ENDDO
      IF(ETMIN.EQ.0) THEN
	  ETMAX=1.05*ETMAX
	ELSEIF(ETMAX.EQ.0) THEN
	  ETMIN=1.05*ETMIN
	ELSE
	  ETMIN=ETMIN-0.05*(ETMAX-ETMIN)
	  ETMAX=ETMAX+0.05*(ETMAX-ETMIN)
	ENDIF
	STOT=BL(II)%SBL(NN)

	INTP=.FALSE.
	IF(DSINTP.GT.0) THEN
	  MDIV=1
		DO N=1,NN
	    JJ=BL(II)%MAGID(N)
	    DS=MAG(JJ)%LENGTH%X
	    IF(DS.GT.DSINTP) MDIV=MAX(MDIV,INT(DS/DSINTP)+1)
	  ENDDO
	  IF(MDIV.GE.2) THEN
	    ALLOCATE(TWISS(MTWISS,0:MDIV),STAT=ISTAT)
	    IF(ISTAT.NE.0) GOTO 930
	    INTP=.TRUE.
	  ENDIF
	ENDIF
	WRITE(FILE,200)
200   FORMAT(' NEWFRAME; SET FONT DUPLEX')
      CALL TDHEAD(FILE)

	DO IBE=1,2
	  IF(IBE.EQ.1) THEN
	    WRITE(FILE,220) XWL,XWH,YWL,YWH,0D0,STOT,0D0,BTMAX,
     %      (XWL+XWH)/2,YWH+0.5,NAM
220       FORMAT(
     %    ' SET WINDOW X',0P2F7.3,' Y',2F7.3,/,
     %    ' SET LIMIT X',1P2D12.4,' Y',2D12.4,/,
     %    ' TITLE',0P2F7.3,' CENTER SIZE 2.0 ',2H'',/,
     %    ' MORE ',1H','Optics of Beamline ',A,1H')
          IF(LET(1).NE.0.OR.LET(2).NE.0) WRITE(FILE,230)
230       FORMAT(' SET AXIS RIGHT OFF')
          Y1=YWL+0.5
	    X1=XWL-1.1
        ELSE
	    IF(LET(1).EQ.0.AND.LET(2).EQ.0) CYCLE
	    WRITE(FILE,400) 0D0,STOT,ETMIN,ETMAX,XWH+1.1,(YWL+YWH)/2
400       FORMAT(' SET COLOR WHITE',/,
     %    ' SET LIMIT X',1P2D12.4,' Y',2D12.4,/,
     %    ' SET AXIS ALL OFF RIGHT ON',/,
     %    ' SET TICKS ALL OFF RIGHT ON',/,
     %    ' SET LABELS ALL OFF RIGHT ON',/,
     %    ' TITLE',0P2F7.3,' ANGLE 90 CENTER SIZE 2.0 ',2H'')
          X1=XWH+1.1
	    Y1=YWL+0.5
          IF(LET(1).EQ.0.OR.LET(2).EQ.0) Y1=YWL+2.5
	  ENDIF
        WRITE(FILE,250)
250     FORMAT(' PLOT AXIS')
      
        DO J=1,2
	    IF(IBE.EQ.2.AND.LET(J).EQ.0) CYCLE
	    WRITE(FILE,260) LBE(IBE),LXY(J)
260       FORMAT('(  ',A,A)
          WRITE(FILE,280) COLOR(IBE)
280       FORMAT(' SET COLOR ',A)
          IF(INTP) THEN
	      DO N=1,NN
	        JJ=BL(II)%MAGID(N)
	        DS=MAG(JJ)%LENGTH%X
	        TWISS(1:MTWISS,0)=BL(II)%TWISS(1:MTWISS,N-1)
	        NDIV=MIN(MDIV,INT(DS/DSINTP)+1)
	        DS=DS/NDIV
	        N1=NDIV
	        IF(N.NE.NN) N1=MAX(0,NDIV-1)
					IF(N1.GT.0) CALL OPTICSINT(JJ,TWISS,NDIV,IRTN1)
	        IF(IBE.EQ.1) THEN
                WRITE(FILE,300) (BL(II)%SBL(N-1)+K*DS,
     %            SQRT(TWISS(1+5*(J-1),K)),K=0,N1)
	        ELSE
	          WRITE(FILE,300) (BL(II)%SBL(N-1)+K*DS,
     %            TWISS(3+5*(J-1),K),K=0,N1)
	        ENDIF
	      ENDDO
	    ELSE
            IF(IBE.EQ.1) THEN
              WRITE(FILE,300) (BL(II)%SBL(N),
     %          SQRT(BL(II)%TWISS(1+5*(J-1),N)),N=0,NN)
	      ELSE
	        WRITE(FILE,300) (BL(II)%SBL(N),
     %          BL(II)%TWISS(3+5*(J-1),N),N=0,NN)
	      ENDIF
	    ENDIF
300       FORMAT(3(1PD11.4,D12.4,';'))
          WRITE(FILE,320) LMOD(J)
320       FORMAT(' JOIN 1 ',A)
          WRITE(FILE,340) X1,Y1,X1,Y1+1.0,LMOD(J)
340       FORMAT(2(2F8.3,';'),'JOIN 1 TEXT ',A)
          IF(IBE.EQ.1) THEN
            WRITE(FILE,360) X1,Y1+1.2,LXY(J)
360         FORMAT(
     %      ' TITLE',0P2F7.3,' ANGLE 90 SIZE 2.0 ',2H'',/,
     %      ' MORE ',1H','2B0',A1,'1 (2m)',1H',/,
     %      ' CASE ',1H','MGX',1X,'X  M  ',1H',/,' PLOT AXIS')
	    ELSE
	      WRITE(FILE,380) X1,Y1+1.2,LXY(J)
380         FORMAT(
     %      ' TITLE',0P2F7.3,' ANGLE 90 SIZE 2.5 ',2H'',/,
     %      ' MORE ',1H','H0',A1,'1 (m)',1H',/,
     %      ' CASE ',1H','GX',1X,'X    ',1H',/,' PLOT AXIS')
	    ENDIF
          Y1=Y1+3.0
        ENDDO
	ENDDO
      CALL BLMAGPIC(FILE,II,XWL,XWH,0.1D0,YWL-0.5)
      IRTN=0
	GOTO 1000
910   IRTN=1001
	IF(MSGLVL.GE.0) WRITE(MSGFL,915) NAM
915   FORMAT(' (SUBR.PRBLOPTICS) Beamline "',A,'" does not exist.')
      GOTO 1000
920   IRTN=1002
	IF(MSGLVL.GE.0) WRITE(MSGFL,925) NAM
925   FORMAT(' (SUBR.PRBLOPTICS) Optics for beamline "',A,'" not ',
     %   'ready. Use BLOPTICS command.')
      GOTO 1000
930   IRTN=1003
	IF(MSGLVL.GE.0) WRITE(MSGFL,935)
935   FORMAT(' (SUBR.PRBLOPTICS) Work are allocation failed.',/,
     %   '    INTERPOLATE parameter too small?')
	GOTO 1000
1000  IF(INTP) DEALLOCATE(TWISS,STAT=ISTAT)
	RETURN
      END

      SUBROUTINE BLMAGPIC(FILE,II,XWL,XWH,YWL,YWH)
	USE BEAMLN
	IMPLICIT NONE
	INTEGER FILE,II
	REAL(8) XWL,XWH,YWL,YWH
	INTEGER JJ,NN,N,N0
	REAL(8) YWM,YC,YN,DYB,DYQ,STOT,DXDS,APERT1(2),X1,X2,Y1,DX,DY
C       YC:  beam center line,  YN: first char of magnet name
      REAL(8)  RATIO/0.3D0/   !  ratio of the vertical space for magnet pictures
	REAL(8) CSIZE/0.15D0/   !  character size

	YWM=RATIO*(YWH-YWL)
      YC=YWH-0.5*YWM
	DYQ=0.5*YWM*0.8
	DYB=DYQ*0.7
      YN=YWH-YWM-0.5*CSIZE

      NN=BL(II)%NEXP
	STOT=BL(II)%STOT
	DXDS=(XWH-XWL)/STOT
	X1=XWL
	Y1=YC
	WRITE(FILE,100)
100   FORMAT('(  Magnet picture',/,' SET COLOR WHITE ')
      WRITE(FILE,120) X1,Y1
120   FORMAT(2F8.3)
      DO N=1,NN
	  JJ=BL(II)%MAGID(N)
	  X2=X1+DXDS*MAG(JJ)%LENGTH%X
	  IF(MAG(JJ)%ANGLE.NE.0.OR.MAG(JJ)%K1%X.NE.0) THEN
	    IF(MAG(JJ)%ANGLE.NE.0) THEN
	      DY=DYB
	    ELSE
	      DY=DYQ
	    ENDIF
	    WRITE(FILE,220) X1,Y1+DY,X2,Y1+DY,X2,Y1-DY,X1,Y1-DY,X1,Y1
220       FORMAT(5(2F8.3,';'),' JOIN 1 TEXT')
          WRITE(FILE,120) X2,Y1
        ELSEIF(MAG(JJ)%LENGTH%X.NE.0) THEN
	    WRITE(FILE,120) X2,Y1
        ENDIF
	  X1=X2
	ENDDO
	WRITE(FILE,240)
240   FORMAT(' JOIN 1 TEXT')
C   magnet name  
C     *  exclude drift-space names
C     *  draw only once when same name appears more than once successively
      N0=0
	DO N=1,NN
	  JJ=BL(II)%MAGID(N)
	  IF(MAG(JJ)%K1%X.EQ.0.AND.MAG(JJ)%ANGLE.EQ.0) THEN
	    IF(MAG(JJ)%LENGTH%X.NE.0) THEN
	      IF(N0.NE.0) THEN
              CALL BLMAGPIC1(FILE,II,N0,N,CSIZE,XWL,DXDS,YN)
	        N0=0
	      ENDIF
	    ENDIF
	  ELSE
	    IF(N0.EQ.0) THEN
			  N0=N
	    ELSE
	      IF(JJ.NE.BL(II)%MAGID(N0)) THEN
	        CALL BLMAGPIC1(FILE,II,N0,N,CSIZE,XWL,DXDS,YN)
	        N0=N
	      ENDIF
	    ENDIF
	  ENDIF
	ENDDO
	IF(N0.NE.0) CALL BLMAGPIC1(FILE,II,N0,NN,CSIZE,XWL,DXDS,YN)
	RETURN
      END

	SUBROUTINE BLMAGPIC1(FILE,II,N0,N,CSIZE,XWL,DXDS,YN)
	USE BEAMLN
	IMPLICIT NONE
	INTEGER FILE,II,N0,N
	REAL(8) CSIZE,XWL,DXDS,YN
	INTEGER J,NC,I
	REAL(8) X1,Y1

	J=BL(II)%MAGID(N0)
	NC=MAG(J)%NC
	IF(N.EQ.N0) THEN
	  X1=XWL+DXDS*0.5*(BL(II)%SBL(N0-1)+BL(II)%SBL(N))
	ELSE
	  X1=XWL+DXDS*0.5*(BL(II)%SBL(N0-1)+BL(II)%SBL(N-1))
	ENDIF
	DO I=1,NC
	  Y1=YN-(I-1)*(CSIZE*1.25D0)
	  WRITE(FILE,100) X1,Y1,10*CSIZE,MAG(J)%NAME(I:I)
100     FORMAT(' TITLE',2F7.3,' SIZE',F5.2,' ',1H',A,1H')
      ENDDO
	RETURN
	END

	  

