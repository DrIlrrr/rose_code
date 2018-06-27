      SUBROUTINE EXFLRZ(NTRANS,ITRANS,TXYS0,BG,EV,ANGLE,AXIS)
C Lorentz transformation of the external field      
      IMPLICIT NONE
      INTEGER NTRANS,ITRANS(NTRANS)
      REAL*8 TXYS0(0:3),BG,EV(3),ANGLE,AXIS(3)
      INCLUDE 'include/extfcm.h'
      INTEGER I,J,II
      REAL*8 R(3,3),GAM,BET,T(0:3,0:3)
      REAL*8 C,SEB(2),VEB(3,2),EBT(3,2),EBL(3,2)
      INTEGER I2(3)/2,3,1/,I3(3)/3,1,2/
C
      IF(LEXTF.LE.0) RETURN
      DO 500 II=1,NTRANS
      GOTO (200,300,400), ITRANS(II)
C  Origin shift
 200  IF(LEXTFB(1).GE.1.OR.LEXTFB(2).GE.1) THEN
        C=EXTFBV(0)*TXYS0(0)
        DO 220 I=1,3
          C=C-EXTFBV(I)*TXYS0(I)
 220    CONTINUE
        DO 240 I=1,2
          IF(LEXTFB(I).GE.1) EXTFSS(I)=EXTFSS(I)+C
 240    CONTINUE
      ENDIF
      GOTO 500
C  Rotation
 300  CALL ROTMAT(ANGLE,AXIS,1,3,1,3,R)
      IF(LEXTFB(1).GE.1.OR.LEXTFB(2).GE.1) 
     %     CALL MATVEC(R,3,3,EXTFBV(1))
      CALL MATVEC(R,3,3,EXTFEB(1,1))
      CALL MATVEC(R,3,3,EXTFEB(1,2))
      GOTO 500
C  Lorentz boost
 400  GAM=SQRT(1+BG**2)
      BET=BG/GAM
      IF(LEXTFB(1).GE.1.OR.LEXTFB(2).GE.1) THEN
        CALL LBOOST(BG,EV,T)
        CALL MATVEC(T,4,4,EXTFBV)
      ENDIF
      DO 460 J=1,2
        SEB(J)=0
        DO 430 I=1,3
          VEB(I,J)=EV(I2(I))*EXTFEB(I3(I),J)-EV(I3(I))*EXTFEB(I2(I),J)
          SEB(J)=SEB(J)+EV(I)*EXTFEB(I,J)
 430    CONTINUE
        DO 440 I=1,3
          EBL(I,J)=SEB(J)*EV(I)
          EBT(I,J)=EXTFEB(I,J)-EBL(I,J)
 440    CONTINUE
 460  CONTINUE
      DO 480 I=1,3
        EXTFEB(I,1)=GAM*(EBT(I,1)+BET*VEB(I,2))+EBL(I,1)
        EXTFEB(I,2)=GAM*(EBT(I,2)-BET*VEB(I,1))+EBL(I,2)
 480  CONTINUE
 500  CONTINUE
C
      RETURN
      END

