      MODULE BBCOM
	  IMPLICIT NONE
	  INTEGER MXY
	  INTEGER BBON,BBFLG(2),NXY(2),NMOM,LBBEL,NBBPL,IFLBPL
	  REAL(8), ALLOCATABLE:: BBQ(:),BBWORK(:)
C            REAL(8) BBQ(MXY**2*2),BBWORK((2*MXY)**2)
	  REAL(8) XYMIN(2,2),XYMAX(2,2),XYCENT(2,2),
     %   BBDXY(2,2),BBXYM(2,2,2),PSIZE,
     %   SBBPL(5),EMAX(2),WGTOUT(2),
     %   BBR00(2),BBEL(2),BBU00
      END MODULE BBCOM