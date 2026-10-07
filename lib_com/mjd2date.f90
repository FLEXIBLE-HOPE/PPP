!*
SUBROUTINE mjd2date(jd,sod,iy,imon,id,ih,im,is)
!!
!! purpose  : transform modified julday to date
!! parameter:
!!    input : jd,sod -- modified julday
!!    output: iy,imon,id,ih,im,is -- date
!! author   : Geng J
!! created  : Nov. 6, 2008
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: jd,iy,imon,id,ih,im
REAL(RL) :: sod,is

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: mjd,doy
  REAL(RL) :: msod

  !! check type of modified julday
  IF (jd .NE. 0) THEN
    mjd=jd
    msod=sod
  ELSE
    mjd=int(sod)
    msod=(sod-mjd)*86400.d0
  END IF

  !! transformation
  CALL mjd2doy(mjd,iy,doy)
  CALL yeardoy2monthday(iy,doy,imon,id)
  CALL sod2hms(msod,ih,im,is)

  RETURN

END SUBROUTINE
