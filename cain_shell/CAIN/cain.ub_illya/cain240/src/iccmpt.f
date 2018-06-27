      SUBROUTINE ICCMPT(P,K,KIN,TXYSE,TXYSG,IS,NAM,GN,WG,ENH,EMIN,
     %    VOL,DT,MSGFL,IRTN)
      IMPLICIT NONE
      INTEGER KIN,IS,GN,MSGFL,IRTN
      CHARACTER*(*) NAM
      REAL*8 P(0:3),K(0:3),TXYSE(0:3),TXYSG(0:3),WG,ENH,EMIN,
     %    VOL,DT
      INTEGER I,IDIV,NDIV,NEV,L,IRTN1
      REAL*8 N(3),PF(0:3),KF(0:3),SPF(3),STKF(3),PD,PROB,DT1,TXYS1(0:3)
      REAL*8 SP(3)/0,0,0/,UV1(3)/0,0,0/,UV2(3)/0,0,0/,STK(3)/0,0,0/,
     %   BFL1(6)/0,0,0,0,0,0/
      REAL*8 PMAX/0.1D0/
      INCLUDE 'include/cnstcm.h'
C
      DO 200 I=1,3
        N(I)=K(I)/K(0)
 200  CONTINUE
      PD=ECHARG*K(0)*CVEL/VOL
      NEV=0
      NDIV=1
      DT1=DT*ENH
      IDIV=0
 300  IDIV=IDIV+1
 320  CALL LNCPGN(P(1),K(0),N,SP,UV1,UV2,STK,PD,DT1,PMAX,0,0,0,
     %   L,PF(1),KF(1),SPF,STKF,0,PROB,IRTN1)
      IF(IDIV.EQ.1.AND.IRTN1.NE.0) THEN
        NDIV=INT(1.5D0*PROB/PMAX)+1
        DT1=DT1/NDIV
        GOTO 320
      ENDIF
      IF(L.NE.0) THEN
        NEV=NEV+1
        IF(NEV.EQ.1) THEN
          DO 340 I=0,3
            TXYS1(I)=TXYSE(I)+TXYSG(I)
 340      CONTINUE
        ENDIF
        PF(0)=SQRT(EMASS**2+PF(1)**2+PF(2)**2+PF(3)**2)
        KF(0)=SQRT(KF(1)**2+KF(2)**2+KF(3)**2)
        IF(PF(0).GT.EMIN) THEN
          CALL ADDONE(0,KIN,GN,NAM,IS,WG,TXYS1,PF,SPF,0,BFL1,IRTN1)
          IF(IRTN1.NE.0) GOTO 900
        ENDIF
        CALL ADDONE(0,1,GN,NAM,IS,WG,TXYS1,KF,STKF,0,BFL1,IRTN1)
        IF(IRTN1.NE.0) GOTO 900
      ENDIF
      IF(IDIV.LT.NDIV) GOTO 300
      IRTN=0
      RETURN
 900  IRTN=1000
      WRITE(MSGFL,910)
 910  FORMAT(' (SUBR.ICCMPT) Too many Bremsstrahlung in a step.')
      RETURN
      END

