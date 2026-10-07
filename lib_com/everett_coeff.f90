!*
SUBROUTINE everett_coeff(ndim,n,ec)
!!
!! PANDA subroutine
!!
!! purpose  : compute coefficients of the everett interpolation
!!
!! parameter: n  -- degree of the interpolation, with 2*n+1 points
!!            ec -- everett coefficients
!!
!! author   : Maorong Ge, Feb. 1993
!!
!! last mod.:
!*
USE PAR
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: n,ndim
REAL(RL) :: ec(0:ndim,0:ndim)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: a(0:n,0:n),factor(0:2*n+1),sign

  !*
  ! Start the exectuable code
  !!--------------------------

  factor(0)=1.d0
  DO i=1, 2*n+1
    factor(i)=factor(i-1)*i
  END DO

  a(0,0)=1.d0
  DO i=1,n
    a(i,0)=-a(i-1,0)*i*i
    DO j=1,i-1
      a(i,j)=a(i-1,j-1)-i*i*a(i-1,j)
    END DO
    a(i,i)=1.d0
  END DO

  !! everett coefficient
  DO k=0,n
    DO j=0,n
      ec(k,j)=0.d0
      DO i=MAX(k,j),n
        IF (mod(i+j,2) .EQ. 1) THEN
          sign=-1.d0
        ELSE
          sign=1.d0
        END IF
        ec(k,j)=ec(k,j)+sign*a(i,k)/(2.d0*i+1.d0)/factor(i+j)/factor(i-j)
      END DO
      IF (j.EQ.0) ec(k,j)=ec(k,j)/2.d0
    END DO
  END DO

  RETURN

END SUBROUTINE
