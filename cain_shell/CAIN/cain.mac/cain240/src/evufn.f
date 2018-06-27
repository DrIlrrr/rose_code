      RECURSIVE SUBROUTINE EVUFN(FNAME,X,NV,GSTRX,Z,ERRMSG)
	USE FLCHTYP
      IMPLICIT NONE
      INTEGER NV
      TYPE(FLCHTYPE) Z,X(NV)
      CHARACTER*(*) FNAME,GSTRX,ERRMSG
      INCLUDE 'include/ctrlcm3.h'
      INCLUDE 'include/ctrlcm.h'
	TYPE(FLCHTYPE) Y
	INTEGER IERARG,IRTN,I,K

	Y%L=1
      Y%X=0
	Y%C(1)=1
	Y%C(2)=0
      DO I=1,MUFN
        K=I
        IF(FNAME.EQ.UFNAME(I)) GOTO 200
      ENDDO
      GOTO 900
 200	IF(NARGUFN(K).LT.1000) THEN
	  IF(NV.NE.NARGUFN(K)) GOTO 910
	ENDIF
	IF(K.LT.IDUFLM) THEN
C--Beam statistics functions
	  CALL EVUFNBS(FNAME,K-IDUFBSF+1,X,NV,GSTRX,Y,IRTN,IERARG,ERRMSG)
C--Luminosity functions
      ELSEIF(K.LT.IDUFLSR) THEN
	  CALL EVUFNLUM(FNAME,K-IDUFLM+1,X,NV,Y,IRTN,IERARG,ERRMSG)
C  Laser-related functions
 	ELSEIF(K.LT.IDUFSPF) THEN
	  CALL EVUFNLASER(FNAME,K-IDUFLSR+1,X,NV,Y,IRTN,IERARG,ERRMSG)
C  Special functions
      ELSEIF(K.LT.IDUFCH) THEN
	  CALL EVUFNMATH(FNAME,K-IDUFSPF+1,X,NV,Y,IRTN,IERARG,ERRMSG)
C  Character-related functons
      ELSEIF(K.LT.IDUFTSTP) THEN
	  CALL EVUFNCHAR(FNAME,X,NV,GSTRX,Y,IRTN,IERARG,ERRMSG)
C  Test particle functions
      ELSEIF(K.LT.IDUFTWS) THEN
	  CALL EVUFNTST(K-IDUFTSTP+1,X,NV,Y,GSTRX,IRTN,IERARG,ERRMSG)
C  Beamline Functions
	ELSEIF(K.LT.IDUFTMAT) THEN
        CALL EVUFNBL(FNAME,1,K-IDUFTWS+1,X,NV,GSTRX,Y,IRTN,IERARG,
     %    ERRMSG)
	ELSE
	  CALL EVUFNBL(FNAME,2,K-IDUFTMAT+1,X,NV,GSTRX,Y,IRTN,IERARG,
     %    ERRMSG)
	ENDIF
      IF(IRTN.GE.200) THEN
        GOTO 990
	ELSEIF(IRTN.GE.101) THEN
	  GOTO (910,920,930), IRTN-100
	ENDIF
	IRTN=0
	Z=Y
	ERRMSG=' '
      RETURN
 900  ERRMSG='Function '//FNAME//' undefined.'
      GOTO 990
 910	ERRMSG='Wrong number of arguments for function '//FNAME
      GOTO 990
 920  WRITE(ERRMSG,925) IERARG,FNAME
 925  FORMAT('Range invalid:',I2,'th argument for function ',A)
      GOTO 990
 930  WRITE(ERRMSG,935) IERARG,FNAME
 935  FORMAT('Range invalid:',I2,'th argument for function ',A)
      GOTO 990

 990  WRITE(MSGFL,995) ERRMSG
 995  FORMAT(1X,A)
      RETURN
      END
