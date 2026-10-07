!*
INTEGER(IT) FUNCTION day_of_year(iday,imonth,iyear)
!!
!! purpose   : given year,month and day, return day of the year
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) iday,imonth,iyear

  !*
  ! The function called
  !!-------------------------
  INTEGER(IT) modified_julday

  !*
  ! Start the exectuable code
  !!--------------------------

  day_of_year=modified_julday(iday,imonth,iyear)-modified_julday(1,1,iyear)+1

  RETURN

END FUNCTION
