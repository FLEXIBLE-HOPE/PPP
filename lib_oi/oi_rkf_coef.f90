!*
SUBROUTINE oi_rkf_coef(rkfc)
!!
!*
USE orbit
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(RKFCOEF) :: rkfc

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT)  :: i, j

  !*
  ! Start the exectable code
  !!----------------------------

  rkfc.m=9
  DO i=0, rkfc.m
    rkfc.alpha(i) = 0.d0
    rkfc.c(i)     = 0.d0
    rkfc.d(i)     = 0.d0
    DO j=0, rkfc.m
      rkfc.beta(i,j) = 0.d0
    END DO
  END DO

  rkfc.alpha(1) = 2.d0/33.d0
  rkfc.alpha(2) = 4.d0/33.d0
  rkfc.alpha(3) = 2.d0/11.d0
  rkfc.alpha(4) = 0.5d0
  rkfc.alpha(5) = 2.d0/3.d0
  rkfc.alpha(6) = 6.d0/7.d0
  rkfc.alpha(7) = 1.d0
  rkfc.alpha(9) = 1.d0

  rkfc.beta(1,0) = 2.d0/33.d0
  rkfc.beta(3,0) = 1.d0/22.d0
  rkfc.beta(4,0) = 43.d0/64.d0
  rkfc.beta(5,0) = -2383.d0/486.d0
  rkfc.beta(6,0) = 10077.d0/4802.d0
  rkfc.beta(7,0) = -733.d0/176.d0
  rkfc.beta(8,0) = 15.d0/352.d0
  rkfc.beta(9,0) = -1833.d0/352.d0
  rkfc.beta(2,1) = 4.d0/33.d0
  rkfc.beta(3,2) = 3.d0/22.d0
  rkfc.beta(4,2) = -165.d0/64.d0
  rkfc.beta(5,2) = 1067.d0/54.d0
  rkfc.beta(6,2) = -5643.d0/686.d0
  rkfc.beta(7,2) = 141.d0/8.d0
  rkfc.beta(9,2) = 141.d0/8.d0
  rkfc.beta(4,3) = 77.d0/32.d0
  rkfc.beta(5,3) = -26312.d0/1701.d0
  rkfc.beta(6,3) = 116259.d0/16807.d0
  rkfc.beta(7,3) = -335763.d0/23296.d0
  rkfc.beta(8,3) = -5445.d0/46592.d0
  rkfc.beta(9,3) = -51237.d0/3584.d0
  rkfc.beta(5,4) = 2176.d0/1701.d0
  rkfc.beta(6,4) = -6240.d0/16807.d0
  rkfc.beta(7,4) = 216.d0/77.d0
  rkfc.beta(8,4) = 18.d0/77.d0
  rkfc.beta(9,4) = 18.d0/7.d0
  rkfc.beta(6,5) = 1053.d0/2401.d0
  rkfc.beta(7,5) = -4617.d0/2816.d0
  rkfc.beta(8,5) = -1215.d0/5632.d0
  rkfc.beta(9,5) = -729.d0/512.d0
  rkfc.beta(7,6) = 7203.d0/9152.d0
  rkfc.beta(8,6) = 1029.d0/18304.d0
  rkfc.beta(9,6) = 1029.d0/1408.d0
  rkfc.beta(9,8) = 1.d0

  rkfc.c(0) = 77.d0/1440.d0
  rkfc.c(3) = 1771561.d0/6289920.d0
  rkfc.c(4) = 32.d0/105.d0
  rkfc.c(5) = 243.d0/2560.d0
  rkfc.c(6) = 16807.d0/74880.d0
  rkfc.c(7)  = 1.d0/270.d0

  rkfc.d(0) = 11.d0/864.d0
  rkfc.d(3) = rkfc%C(3)
  rkfc.d(4) = rkfc%C(4)
  rkfc.d(5) = rkfc%C(5)
  rkfc.d(6) = rkfc%C(6)
  rkfc.d(8) = 11.d0/270.d0
  rkfc.d(9) = 11.d0/270.d0

  RETURN

END SUBROUTINE
