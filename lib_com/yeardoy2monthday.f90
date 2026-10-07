!*
SUBROUTINE yeardoy2monthday(iyear,idoy,imonth,iday)
!!
!! purpose   : get month and day given year and day of year 
!! parameters: iyear,idoy -- year and day of year
!!             imonth,iday -- month and day in the yar
!!
!! created by: Maorong Ge
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: iyear,idoy,imonth,iday

  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: days_in_month(12),id
  DATA days_in_month/31,28,31,30,31,30,31,31,30,31,30,31/

  !*
  ! Start the exectuable code
  !!--------------------------

  days_in_month(2)=28
  IF (MOD(iyear,4).EQ.0 .AND. (MOD(iyear,100).NE.0 .OR. MOD(iyear,400).EQ.0)) days_in_month(2)=29
  id  =idoy
  DO imonth=1,12
    id=id-days_in_month(imonth)
    IF (id .GT. 0) CYCLE
    iday=id+days_in_month(imonth)
    EXIT
  END DO

  RETURN

END SUBROUTINE
