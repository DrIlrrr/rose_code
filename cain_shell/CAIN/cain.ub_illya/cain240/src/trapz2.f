      SUBROUTINE TRAPZ2(NDIV,XYS0,EP0,DT,FXYS0,BBFLD0,EXTFXY,
     %   BBFLAV,L,DS,KIN,NXY,BBQ)
      IMPLICIT NONE
      INTEGER NDIV,L,KIN,NXY(2)
      REAL*8 XYS0(3),EP0(0:3),DT,FXYS0(3),BBFLD0(3,2),EXTFXY(3),
     %       BBFLAV(3,2),DS,BBQ(*)
      INTEGER IDIV,I,IOUT
      REAL*8 DT1,WG0,WG,XYS1(3),EP1(0:3),FXYS1(3),FF(3),BBFLD(3,2)
      REAL*8 EMASS2/0.2611200393D12/
C
      DT1=DT/NDIV
      WG0=1/DFLOAT(NDIV)
      DO 200 I=1,3
        BBFLAV(I,1)=WG0*0.5D0*BBFLD0(I,1)
        BBFLAV(I,2)=WG0*0.5D0*BBFLD0(I,2)
 200  CONTINUE
      DO 400 IDIV=1,NDIV
        WG=WG0
        IF(IDIV.EQ.NDIV) WG=0.5D0*WG0
        DO 340 I=1,2
          XYS1(I)=XYS0(I)+EP0(I)/EP0(0)*DT1+FXYS0(I)/EP0(0)*DT1**2/2
 340    CONTINUE
        CALL BBKICK0(L,DS,KIN,XYS1,EP0,BBFLD,FF,IOUT,NXY,BBQ)
        DO 350 I=1,3
          FXYS1(I)=EXTFXY(I)+FF(I)
          BBFLAV(I,1)=BBFLAV(I,1)+WG*BBFLD(I,1)
          BBFLAV(I,2)=BBFLAV(I,2)+WG*BBFLD(I,2)
 350    CONTINUE
        DO 360 I=1,3
          EP1(I)=EP0(I)+(FXYS0(I)+FXYS1(I))/2*DT1
 360    CONTINUE
        EP1(0)=SQRT(EMASS2+EP1(1)**2+EP1(2)**2+EP1(3)**2)
        DO 370 I=1,3
          XYS1(I)=XYS0(I)+(EP0(I)/EP0(0)+EP1(I)/EP1(0))/2*DT1
          XYS0(I)=XYS1(I)
          FXYS0(I)=FXYS1(I)
 370    CONTINUE
        DO 380 I=0,3
          EP0(I)=EP1(I)
 380    CONTINUE
 400  CONTINUE
      RETURN
      END
