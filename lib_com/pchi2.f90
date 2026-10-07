SUBROUTINE PCHI2(N,Q,L,XX)
!!
!! PURPOSE    :  COMPUTE NORMAL FUNCTION DEVIATESS
!!
!! PARAMETERS :
!!         IN :  N : DERGEE OF FREEDOM                                   I*4
!!               Q : UP STAT. (AF)                                       R*4
!!               L : L=O :                                               I*4
!!                   L=1 : USE NUTON INTERATION MOTHED
!!        OUT :  XX: VALUE OF DENCITY DEVIATESS                          R*8
!!
!*
USE ISO_FORTRAN_ENV
IMPLICIT REAL*8 (A-H,O-Z)

  IF(Q .LE. 0.D0) THEN
    WRITE(ERROR_UNIT,'(a,f5.1,a)') ' ****ERROR(pchi2): INPUT AF = ',Q,' <= 0'
    CALL exit(1)
  END IF

  IF (N.EQ.1) GOTO 200
  IF (N.EQ.2) GOTO 300

  P=1.D0-Q
  CALL PNORMAL(Q,X)
  W=2.D0/(9.D0*DBLE(N))
  IF (W.LT.0.D0) W=0.D0
  X0=DBLE(N)*(1.D0-W+X*DSQRT(W))**3

  IF(L.EQ.0) THEN
    XX=X0
    GOTO 400
  END IF

  K=0

100 CALL CHI2(N,X0,PP,D)

  IF (ABS((D-0.D0)).LT.1.D-20) THEN
    XX=X0
    GOTO 400
  END IF

  XX=X0-(PP-P)/D
  IF (ABS(X0-XX).LE.((1.D-4)*ABS(XX))) GOTO 400

  K=K+1
  IF (K.GE.30) GOTO 400

  X0=XX
  GOTO 100

  200 CALL PNORMAL(Q/2.D0,X)
  XX=X*X
  GOTO 400

300 XX=-2.D0*DLOG(Q)

400 CONTINUE

  RETURN

END SUBROUTINE
