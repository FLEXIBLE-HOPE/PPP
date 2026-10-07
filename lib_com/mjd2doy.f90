!*
SUBROUTINE mjd2doy(jd,iyear,idoy)
!!
!! TO GET THE YEAR AND DAY OF YEAR (DOY) GIVEN MODIFIED JULIAN DAY
!! THIS ALGORITHM IS A LITTLE BIT SOLOWER BUT VALIDATES FOR LONG TIME
!!
!*
USE par
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
INTEGER(IT) :: jd,iyear,idoy

  !*
  ! THE FUNCTIONS CALLED
  !!-----------------------------
  INTEGER(IT) :: modified_julday

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  iyear=(jd+678940)/365
  idoy=jd-modified_julday(1,1,iyear)
  DO WHILE(idoy .LE. 0)
    iyear=iyear-1
    idoy=jd-modified_julday(1,1,iyear)+1
  END DO

  RETURN

END SUBROUTINE
