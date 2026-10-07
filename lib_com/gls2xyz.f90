!*
SUBROUTINE gls2xyz(cprn,neph,eph,mjd,sod,xsat,clk)
!!
!*
USE brdeph
IMPLICIT NONE

!*
! The arguments
!!--------------------
CHARACTER(LEN_PRN) :: cprn
INTEGER(IT) :: neph,mjd
TYPE(GLONASS_BRDEPH) :: eph(1:*)
REAL(RL) :: sod,xsat(1:*),clk


  !*
  ! The local variables
  !!------------------------------
  INTEGER(IT) :: i,ieph,nsign
  REAL(RL) :: dt,dtmin,dt1,dt2,x(6)

  !*
  ! Start of executable code
  !!------------------------------

  ieph=0
  dtmin=12.d0    ! hours
  dt=0.d0
  DO i=1, neph
    IF (eph(i).cprn .EQ. cprn) THEN
      dt=(mjd-eph(i).mjd)*24.d0+(sod-eph(i).sod)/3600.d0
      IF(ABS(dt) .LE. dtmin) THEN
        dtmin=ABS(dt)
        ieph=i
      END IF
    END IF
  END DO
  IF(ieph .EQ. 0) RETURN

  dt=(mjd-eph(ieph).mjd)*86400.d0+sod-eph(ieph).sod
  clk=eph(ieph).tau

  dt1=dt-INT(dt/60.d0)*60.d0
  dt2=dt1
  nsign=1
  IF (dt .LT. 0.d0) nsign=-1
  CALL glsinit(x,eph(ieph))

  IF (dt1 .NE. 0.d0 ) THEN
    CALL glsrkf4(dt1,x,eph(ieph))
  END IF

  DO WHILE(DABS(dt-dt1) .GE. DABS(60.d0) )
    IF ((dt-dt1) .NE. 0.d0 ) THEN
      CALL glsrkf4(nsign*60.d0,x,eph(ieph))
      dt1=dt1+nsign*60.d0
    END IF
  END DO

  CALL pz902wgs84(mjd,sod,x,xsat,'MCC')
  xsat(1:6)=xsat(1:6)*1.d3

  RETURN

END SUBROUTINE
