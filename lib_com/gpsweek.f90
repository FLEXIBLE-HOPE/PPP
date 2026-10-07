!*
SUBROUTINE gpsweek(day,month,year,week,wd)
!!
!! TO GET GPS WEEK NUMBER AND DAY OF WEEK GIVEN PARTICULAR
!! DAY MONTH AND YEAR OR DOY AND YEAR
!!
!! MAORONG GE: CREATED
!!
!*
USE par
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
INTEGER(IT) :: year,month,day,week,wd
REAL(RL) :: sod,sow

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  INTEGER(IT) :: mjd

  !*
  ! THE FUNCTIONS CALLED
  !!------------------------------
  INTEGER(IT) :: modified_julday

  !*
  ! THE EXECTUABLE CODES
  !!------------------------------

  IF (year .NE. 0) THEN
    IF (month .NE. 0) THEN
      mjd=modified_julday(day,month,year)
    ELSE
      mjd=modified_julday(1,1,year)+day-1
    END IF
  ELSE
    mjd=day
  END IF

  week=(mjd-44244)/7
  wd  =mjd-44244-week*7

  RETURN

END SUBROUTINE
