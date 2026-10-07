!*
SUBROUTINE shcoef_order_wise(lmax,lmin,ltog)
!
! Procedure for computing the sphere harmonic coefficients by order wise.
! And the S-part follwing the C-part with same degree and order.
!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: lmin, lmax
INTEGER(IT) :: ltog(0:MAXGRADEG,0:MAXGRADEG,2)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: lmin_cur,k,m,l

  !*
  ! Start the exectuable code
  !!-------------------------

  ltog=0
  k = 0
  DO m = 0, lmax
    lmin_cur = MAX(lmin,m)
    DO l = lmin_cur, lmax
      k = k+1
      ltog(l,m,1) = k
      IF (m .GT. 0) THEN
        k = k+1
        ltog(l,m,2) = k
      END IF
    END DO
  END DO

  !*
  ! Return
  !!----------
  RETURN
 
END SUBROUTINE shcoef_order_wise
