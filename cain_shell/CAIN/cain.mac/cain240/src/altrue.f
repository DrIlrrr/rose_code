      SUBROUTINE ALTRUE(L,N)
      IMPLICIT NONE
      INTEGER N,L(N)
      INTEGER I
      DO 200 I=1,N
        IF(L(I).NE.0) RETURN
 200  CONTINUE
      DO 220 I=1,N
        L(I)=1
 220  CONTINUE
      RETURN
      END
