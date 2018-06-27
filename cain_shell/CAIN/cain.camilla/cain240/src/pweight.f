      FUNCTION PWEIGHT(N)
	USE BEAMCM
      IMPLICIT NONE
      INTEGER N
      REAL*8 PWEIGHT
C      INCLUDE 'include/beamcm.h'
      PWEIGHT=0
      IF(N.LE.0.OR.N.GT.NP) RETURN
      IF(LOST(N).NE.0.AND.LOST(N).NE.2) RETURN
      PWEIGHT=WGT(N)
      RETURN
      END
