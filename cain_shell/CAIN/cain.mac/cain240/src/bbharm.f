      SUBROUTINE BBHARM(NMOM,R00,QMOM,NP,IFLG,FXY)
C Harmonic expansion
      IMPLICIT NONE
      INTEGER NMOM,NP,IFLG(NP)
      REAL*8 R00,FXY(2,NP)
      COMPLEX*16 QMOM(0:NMOM)
      INTEGER N,NM
      COMPLEX*16 CZ,CF
C
      CALL CPUTIM('BBHARM',1)
      DO 300 N=1,NP
        IF(IFLG(N).NE.2) GOTO 300
        CZ=1/DCMPLX(FXY(1,N)/R00,FXY(2,N)/R00)
        CF=QMOM(NMOM)
        IF(NMOM.GE.1) THEN
          DO 220 NM=NMOM-1,0,-1
            CF=CF*CZ+QMOM(NM)
 220      CONTINUE
        ENDIF
        CF=CF*CZ
        FXY(1,N)=-DREAL(CF)
        FXY(2,N)=+DIMAG(CF)
 300  CONTINUE
      CALL CPUTIM('BBHARM',2)
      RETURN
      END
