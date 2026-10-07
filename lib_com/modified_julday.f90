!*
INTEGER(IT) FUNCTION modified_julday(iday,imonth,iyear)
!!
!! TO COMPUTE THE MODIFIED JULIAN DAY
!!
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
INTEGER(IT) :: imonth,iday,iyear

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  INTEGER(IT) :: iyr, doy_of_month(12)
  DATA doy_of_month /0,31,59,90,120,151,181,212,243,273,304,334/

  !*
  ! THE EXECTUABLE CODE
  !!------------------------------

  ! CHECK THE INPUT DATA
  IF ((iyear.LT.0 .OR. imonth.LT.0 .OR. iday.LT.0 .OR. imonth.GT.12 .OR. iday.GT.366) .OR. (imonth.NE.0 .AND. iday.GT.31)) THEN
    WRITE(ERROR_UNIT,'(A,I4,2(1X,I2))') '***ERROR(modified_julday): incorrect date(year,month,day): ',iyear,imonth,iday
    CALL exit(1)
  END IF

  iyr=iyear
  IF (imonth .LE. 2) iyr=iyr-1
  modified_julday=365*iyear-678941+iyr/4-iyr/100+iyr/400+iday
  IF (imonth .NE. 0) modified_julday=modified_julday+doy_of_month(imonth)

  RETURN

END FUNCTION
