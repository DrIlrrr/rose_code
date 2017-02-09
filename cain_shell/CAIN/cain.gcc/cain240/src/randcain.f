C-------------------- RANDCAIN ------------------------------------
      FUNCTION RANDCAIN()
C   Generates platform-independent random number
C   Must be initialized by RANINI
      IMPLICIT NONE
      REAL*8 RANDCAIN
C
      COMMON /RAN11C/MCGN,SRGN
      INTEGER*4 MCGN,SRGN
C
      INTEGER*4 TEMP,MCGN_A,MCGN_B,MULT,M2
      PARAMETER (MULT=69069)
      PARAMETER (M2=MULT-65536)
      REAL*8 RANGEX
      PARAMETER (RANGEX=4294967296.0D0)
C
      TEMP= IEOR(SRGN,ISHFT(SRGN,-15))
      SRGN= IEOR(TEMP,ISHFT(TEMP,17))
C Following 6 lines are equivalent to MCGN= (MCGN*MULT) MOD 2**32
      MCGN_A= ISHFT(MCGN,-16)
      MCGN_B= IAND(MCGN,65535)
      TEMP= M2*MCGN_B
      MCGN_A= MCGN_B+M2*MCGN_A+ISHFT(TEMP,-16)
      MCGN_B= IAND(TEMP,65535)
      MCGN= ISHFT(MCGN_A,16)+MCGN_B
      TEMP= IEOR(MCGN,SRGN)
      RANDCAIN= DBLE(TEMP)/RANGEX+0.5
      RETURN
      END
C-------------------- RANDS ----------------------------------------
      FUNCTION RANDS(IR)
      IMPLICIT NONE
      INTEGER*4 IR
      REAL*4 RANDS
C
      COMMON /RAN11C/MCGN,SRGN
      INTEGER*4 MCGN,SRGN
C
      INTEGER*4 TEMP,MCGN_A,MCGN_B,MULT,M2
      PARAMETER (MULT=69069)
      PARAMETER (M2=MULT-65536)
      REAL*4 RANGEX
      PARAMETER (RANGEX=4294967296.0E0)
C
      TEMP= IEOR(SRGN,ISHFT(SRGN,-15))
      SRGN= IEOR(TEMP,ISHFT(TEMP,17))
C Following 6 lines are equivalent to MCGN= (MCGN*MULT) MOD 2**32
      MCGN_A= ISHFT(MCGN,-16)
      MCGN_B= IAND(MCGN,65535)
      TEMP= M2*MCGN_B
      MCGN_A= MCGN_B+M2*MCGN_A+ISHFT(TEMP,-16)
      MCGN_B= IAND(TEMP,65535)
      MCGN= ISHFT(MCGN_A,16)+MCGN_B
      TEMP= IEOR(MCGN,SRGN)
      RANDS= TEMP/RANGEX+0.5
      RETURN
      END
C------------------- RANDN -----------------------------------------
      SUBROUTINE RANDN(X,N)
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 X(N)
      COMMON /RAN11C/MCGN,SRGN
      INTEGER*4 MCGN,SRGN
C
      INTEGER*4 TEMP,MCGN_A,MCGN_B,MULT,M2,I
      PARAMETER (MULT=69069)
      PARAMETER (M2=MULT-65536)
      REAL*8 RANGEX
      PARAMETER (RANGEX=4294967296.0D0)
C
      DO 200 I=1,N
        TEMP= IEOR(SRGN,ISHFT(SRGN,-15))
        SRGN= IEOR(TEMP,ISHFT(TEMP,17))
C Following 6 lines are equivalent to MCGN= (MCGN*MULT) MOD 2**32
        MCGN_A= ISHFT(MCGN,-16)
        MCGN_B= IAND(MCGN,65535)
        TEMP= M2*MCGN_B
        MCGN_A= MCGN_B+M2*MCGN_A+ISHFT(TEMP,-16)
        MCGN_B= IAND(TEMP,65535)
        MCGN= ISHFT(MCGN_A,16)+MCGN_B
        TEMP= IEOR(MCGN,SRGN)
        X(I)= DBLE(TEMP)/RANGEX+0.5D0
 200  CONTINUE
      RETURN
      END
C

