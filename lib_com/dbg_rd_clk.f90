!*
SUBROUTINE dbg_rd_clk(mjd,sod,CKF,SAT)
!*
USE ckdctrl
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)

  !*
  ! The local variables
  !!-----------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,j,mjdf,mjdx,ierr
  INTEGER(IT) :: lfnclk,iy,imon,id,ih,im,isat,nsat
  CHARACTER(LEN_PRN) :: cprn
  REAL(RL) :: sodf,sodx,sec,dt,dui(2)
  CHARACTER(LEN_STRING) :: line
  LOGICAL(LG) :: lexist

  DATA lfirst /.TRUE./
  SAVE lfirst,lfnclk,mjdf,sodf

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! open file
  IF (lfirst) THEN

    lfirst=.FALSE.

    !! open clock file
    lfnclk=0
    INQUIRE(FILE=CKF.flncck,EXIST=lexist)
    IF (lexist .EQ. .TRUE.) THEN
      lfnclk=get_valid_unit(10)
      OPEN(UNIT=lfnclk,FILE=CKF.flncck,STATUS='OLD')
      DO WHILE(.TRUE.)
        READ(lfnclk,'(A)') line
        IF (INDEX(line,'END OF HEADER') .NE. 0) EXIT
      END DO
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***ERROR(dbg_rd_clkfcb): the clock file is not exist'
      CALL exit(1)
    END IF

  END IF

  DO isat=1, CKF.nprn
    SAT(isat).sclock=0.d0
  END DO

  !! read a epoch of clocks
  DO WHILE(.TRUE. .AND. lfnclk.NE.0)
    READ(lfnclk,'(A)',END=100) line
    READ(line(9:),*) iy,imon,id,ih,im,sec
    IF (mjd .EQ. 0) THEN
      mjd=modified_julday(id,imon,iy)
      sod=ih*3600.d0+im*60.d0+sec
    ELSE
      mjdx=modified_julday(id,imon,iy)
      sodx=ih*3600.d0+im*60.d0+sec
      ! IF (timdif(mjd,sod,mjdx,sodx) .NE. 0.d0) THEN
      IF (timdif(mjd,sod,mjdx,sodx) .LT. 0.d0) THEN
        BACKSPACE(lfnclk)
        EXIT
      END IF
    END IF
    READ(line(4:),*) cprn,iy,imon,id,ih,im,sec,i,dt
    isat=pointer_string(CKF.nprn,CKF.cprn,cprn)
    IF (isat .NE. 0) SAT(isat).sclock=dt
  END DO

100 CONTINUE

  RETURN

END SUBROUTINE
