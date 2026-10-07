!*
SUBROUTINE normalize_coeff(dir,ndegree,morder,ndim,mdim,norm)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
CHARACTER(LEN=*) :: dir
INTEGER(IT) :: ndegree,morder,ndim,mdim
REAL(RL) :: norm(ndim,0:mdim)

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: n,m

  !*
  ! Start the exectuable code
  !!--------------------------

  DO n=1, ndegree
    norm(n,0)=DSQRT(2.d0*n+1.d0)
    norm(n,1)=DSQRT((2.d0*n+1.d0)/(n*(n+1.d0)/2.d0))
  END DO

  DO m=2, morder
    DO n=m, ndegree
      norm(n,m)=norm(n,m-1)/DSQRT(DBLE((n+m)*(n-m+1)))
    END DO
  END DO

  IF (dir(1:2) .EQ. 'de') THEN
    DO n=1, ndegree
      DO m=0, morder
        IF (m .LE. n) norm(n,m)=1.d0/norm(n,m)
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE
