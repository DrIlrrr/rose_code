      INTEGER MSGFL,OUTFL,OUTFL2,TDFL,
     %  MSGLVL,JRAND,IDEBUG
      REAL*8 SMESH,EQMERR
      COMMON/CTRLCM/SMESH,EQMERR,MSGFL,OUTFL,OUTFL2,TDFL,
     %  MSGLVL,JRAND,IDEBUG
      INTEGER ECHO,ISPIN
      COMMON/FLAGCM/ECHO,ISPIN
      INTEGER RDFL
      COMMON/RDFLCM/RDFL
      INTEGER LRAN11
      COMMON/SYSFLG/LRAN11
      INTEGER OSID,MSGDEST,MSGFL0,MSGCOUNT
C       OSID     0: unix, 1: windows
C       MSGDEST  1: console, 2: file (0 for unix)
C         When MSGDEST=2, MSGFL goes to a file and its copy goes to
C         MSGFL0 (console)
      COMMON/SYSNAME/OSID,MSGDEST,MSGFL0,MSGCOUNT
