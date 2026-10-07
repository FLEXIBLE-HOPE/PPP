!*
SUBROUTINE oi_rkf(mjd,t,x,h,xj,nequ,acc,amat,bmat,cmat)
!!
!*
USE orbit
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: mjd,nequ
REAL(RL) :: t,h,x(1:*),xj(1:*)
REAL(RL) :: acc(1:*),amat(3,3),bmat(3,3),cmat(1:*)

  !*
  ! The local varible
  !!---------------
  TYPE(RKFCOEF) :: rkfc

  INTEGER(IT) :: i,j,k,kk,iequ
  REAL(RL) :: ti,f(0:10,MAXEQUS),xi(MAXEQUS)

  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./
  SAVE rkfc,lfirst

  !*
  ! Start the executable code
  !!------------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    CALL oi_rkf_coef(rkfc)
    lfirst=.FALSE.
  END IF

  f=0.d0

  ! Compute right function value
  DO i=0, rkfc.m
    ti=t+rkfc.alpha(i)*h
    DO iequ=1, nequ
      xi(iequ)=x(iequ)
      DO j=0, i-1
        xi(iequ)=xi(iequ)+h*rkfc.beta(i,j)*f(j,iequ)
      END DO
    END DO
    CALL oi_fright_acc(mjd,ti,xi,acc,amat,bmat,cmat)
    DO j=1, nequ/6
      kk=(j-1)*6
      DO k=1, 3
        f(i,kk+k)=xi(kk+k+3)
        f(i,kk+3+k)=acc((j-1)*3+k)
      END DO
    END DO
  END DO

  DO iequ=1, nequ
    xi(iequ)=x(iequ)
    xj(iequ)=x(iequ)
    DO i=0, rkfc.m
      xi(iequ)=xi(iequ)+h*rkfc.c(i)*f(i,iequ)
      xj(iequ)=xj(iequ)+h*rkfc.d(i)*f(i,iequ)
    END DO
  END DO

  RETURN

END SUBROUTINE
