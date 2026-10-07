!*
SUBROUTINE oi_adams_coef(n,pbeta,cbeta)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!---------------------
INTEGER(IT) :: n
REAL(RL) :: pbeta(0:*), cbeta(0:*)

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i,j
  REAL(RL):: pgama(0:n+1),cgama(0:n+1),bino(0:n+1,0:n+1)
  REAL(RL):: s,factor(0:n+1)

  !*
  ! Start the exectuable code
  !!----------------------------

  ! Integer factor for bino. computation
  factor(0)=1.d0
  DO i=1, n
    factor(i)=factor(i-1)*i
  END DO

  !! Bion.
  DO i=0, n
    DO j=0, i
      bino(i,j)=factor(i)/factor(j)/factor(i-j)
    END DO
  END DO

  !! Gama
  DO i=0, n
    pgama(i) =1.d0
    cgama(i)=0.d0
    DO j=0, i-1
      pgama(i)=pgama(i)-pgama(j)/(i-j+1.d0)
      cgama(i)=cgama(i)-cgama(j)/(i-j+1.d0)
    END DO
    IF (i .EQ. 0) cgama(i)=1.d0
  END DO

  ! Beta for prediction (beta) and correction (beta1)
  s=1.d0
  DO i=0, n-1
    pbeta(i)=0.d0
    cbeta(i)=0.d0
    DO j=i, n-1
      pbeta(i)=pbeta(i)+bino(j,i)*pgama(j)
      cbeta(i)=cbeta(i)+bino(j,i)*cgama(j)
    END DO
    pbeta(i)=s*pbeta(i)
    cbeta(i)=s*cbeta(i)
    s=-s
  END DO

  RETURN

END SUBROUTINE
