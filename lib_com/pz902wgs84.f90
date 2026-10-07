!*
SUBROUTINE pz902wgs84(mjd,sod,pos,xsat,trans)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: mjd
CHARACTER(LEN=*) :: trans
REAL(RL) :: sod,pos(1:*),xsat(1:*)

  !*
  ! Local variables
  !!-------------------------------
  INTEGER(IT) :: i, j
  REAL(RL) :: biase(3),mat(3,3),factor,vec(3)

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!------------------------------

  IF (timdif(mjd,sod,54362,86367.d0) .GE. 0.d0 ) THEN
    factor=1.0d0
    biase(1)=-0.36d-3
    biase(2)=0.08d-3
    biase(3)=0.18d-3
    DO i=1, 3
      DO j=1, 3
        IF (i .EQ. j) THEN
          mat(i,j)=1.d0
        ELSE
          mat(i,j)=0.d0
        END IF
      END DO
    END DO
  ELSE
    SELECT CASE(TRIM(trans))
      CASE ("MCC")
        factor=1.0d0+22.0d-9
        biase(1)=-0.47d-3
        biase(2)=-0.51d-3
        biase(3)=-1.56d-3
        mat(1,1)=1.0d0
        mat(1,2)=-1.728d-6
        mat(1,3)=-1.7d-8
        mat(2,1)=1.728d-6
        mat(2,2)=1.0d0
        mat(2,3)=7.6d-8
        mat(3,1)=1.7d-8
        mat(3,2)=-7.6d-8
        mat(3,3)=1.0d0
      CASE ("RUS")
        factor=1.0d0-1.2d-7
        biase(1)=-1.1d-3
        biase(2)=-0.3d-3
        biase(3)=-0.9d-3
        mat(1,1)=1.0d0
        mat(1,2)=-8.2d-7
        mat(1,3)=0.0d0
        mat(2,1)=8.2d-7
        mat(2,2)=1.0d0
        mat(2,3)=0.0d0
        mat(3,1)=0.0d0
        mat(3,2)=0.0d0
        mat(3,3)=1.0d0
      CASE ("LMU")
        factor=1.0d0
        biase(1:3) = 0.0d0
        mat(1,1)=1.0d0
        mat(1,2)=-1.6d-6
        mat(1,3)=0.0d0
        mat(2,1)=1.6d-6
        mat(2,2)=1.0d0
        mat(2,3)=0.0d0
        mat(3,1)=0.0d0
        mat(3,2)=0.0d0
        mat(3,3)=1.0d0
       CASE ("MIT")
        factor=1.0d0
        biase(1)=0.0d0
        biase(2)=2.5d-3
        biase(3)=0.0d0
        mat(1,1)=1.0d0
        mat(1,2)=-1.9d-6
        mat(1,3)=0.0d0
        mat(2,1)=1.9d-6
        mat(2,2)=1.0d0
        mat(2,3)=0.0d0
        mat(3,1)=0.0d0
        mat(3,2)=0.0d0
        mat(3,3)=1.0d0
      CASE DEFAULT
        WRITE(ERROR_UNIT,'(A)') '***ERROR(pz902wgs84): unknown transformation type '//TRIM(trans)
        CALL exit(1)
    END SELECT
  ENDIF

  CALL matmpy(mat,pos,vec,3,3,1)

  DO i=1, 3
    xsat(i)=biase(i)+factor*vec(i)
  END DO

  !*
  ! RETURN
  !---------------
  RETURN

END SUBROUTINE
