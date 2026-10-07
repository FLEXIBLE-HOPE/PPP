!*
SUBROUTINE oi_accelerator(lpart,flnacc,mjd,sod,npar,pname,xics,rot,cmat,acc)
!!
!*
USE par
IMPLICIT NONE

!*
! Declare
!!-------------------
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar,mjd
CHARACTER(LEN=*) :: pname(1:*), flnacc
REAL(RL) :: xics(1:*),sod,cmat(1:*),acc(1:*),rot(3,3)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: rmat(3,3),f(3),qcof(4),ac(3)

  LOGICAL(LG) :: lquant

  INTEGER(IT), PARAMETER :: MAXPARLOC=6
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(LEN_ORBPAR) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname &
      /'ACC_Bx    ','ACC_By    ','ACC_Bz    ', &
       'ACC_Kx    ','ACC_Ky    ','ACC_Kz    '/

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the executable code
  !!-----------------------------

  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO

  DO i=1, MAXPARLOC
    param(i)=0.d0
    IF (i .GT. 3) param(i)=1.d0
    IF (ltog(i) .NE. 0) param(i)=param(i)+xics(ltog(i))
  END DO

  ! Interpolate acc and q.
  ! CALL readLeoAcc(flnacc,mjd,sod,ac,qcof,lquant)

  ! Fotation matrix (if quanternion in acclerat file use it, else use input one)
  IF (lquant .EQ. .TRUE.) THEN
    ! CALL quanternion(qcof,rmat)
  ELSE
    DO i=1, 3
      DO j=1, 3
        rmat(i,j)=rot(i,j)
      END DO
    END DO
  END IF

  ! Force model parameter correction and convert unit from m/s/s to km/s/s
  DO i=1, 3
    f(i)=(param(i)+param(i+3)*ac(i))*1.d-3
  END DO

  CALL matmpy(rmat,f,f,3,3,1)

  DO i=1, 3
    acc(i)=acc(i)+f(i)
  END DO

  IF (.NOT. lpart) RETURN

  DO i=1, MAXPARLOC
    IF (ltog(i) .EQ. 0) CYCLE
    j=(ltog(i)-6 -1)*3
    DO k=1, 3
      IF (i .LE. 3) THEN
        cmat(j+k)=rmat(k,i)*1.d-3
      ELSE
        cmat(j+k)=rmat(k,i-3)*ac(i-3)*1.d-3
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE
