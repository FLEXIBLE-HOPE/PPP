!*
SUBROUTINE tim_tag(mjd,sod,cweek,cwkd,chour)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!--------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod
CHARACTER(LEN=*) cweek,cwkd,chour

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: tmjd,iweek,iwkd,ihour
  REAL(RL) :: tsod,sow

  !*
  ! Start the exectuable code
  !!--------------------------------

  !! 10 min decrement in case of orbit latency
  CALL timinc(mjd,sod,-600.d0,tmjd,tsod)

  !! GPS week and day of the week
  CALL mjd2wksow(tmjd,tsod,iweek,sow)
  iwkd = INT(sow/86400.d0)

  !! seconds of the day
  sow = dmod(sow,86400.d0)
  ihour = INT(sow/3600.d0)

  !! decide which day
  IF (ihour .LT. 3) THEN
    ihour=18
    IF (iwkd .EQ. 0) THEN
      iwkd =6
      iweek=iweek-1
    ELSE
      iwkd =iwkd-1
    END IF
  ELSE IF (ihour.GE.3 .AND. ihour.LT.9) THEN
    ihour=0
  ELSE IF (ihour.GE.9 .AND. ihour.LT.15) THEN
    ihour=6
  ELSE IF (ihour.GE.15 .AND. ihour.LT.21) THEN
    ihour=12
  ELSE
    ihour=18
  END IF

  !! output to strings
  WRITE(cweek,'(I4.4)') iweek
  WRITE(cwkd, '(I1.1)') iwkd
  WRITE(chour,'(I2.2)') ihour

  RETURN

END SUBROUTINE
