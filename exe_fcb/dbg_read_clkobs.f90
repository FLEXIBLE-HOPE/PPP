!*
SUBROUTINE dbg_read_clkobs(mjd,sod,i,snam,cprn,ifreq,clk,obs)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: mjd, i, ifreq(MAXSAT)
REAL(RL) :: sod, clk(MAXSAT), obs(MAXSAT,2*MAXFREQ,MAXSIT)
CHARACTER(LEN_SITENAME) :: snam(MAXSIT)
CHARACTER(LEN_PRN) :: cprn(MAXSAT)

  !*
  ! The local variables
  !!--------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: lfn, ierr
  DATA lfirst /.TRUE./
  SAVE lfirst, lfn

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst = .FALSE.
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE='clk_obs_dbg',FORM='UNFORMATTED')
  END IF

  READ(lfn) mjd,sod,i,snam,cprn,ifreq,clk,obs

  RETURN

END SUBROUTINE

