! read widelane FCBs
! Jianghui Geng
! March 2 2012

SUBROUTINE fcb_rd_wide(CKF,UPD,bmjd,bsod)
!!
!*
USE ckdctrl
USE ambiguity
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(CKDCFG) :: CKF
TYPE(FCB) :: UPD
INTEGER(IT) :: bmjd
REAL(RL) :: bsod

  !*
  ! The local variables
  !!---------------------
  LOGICAL(LG) :: lexist
  INTEGER(IT) :: lfn,mjd,iy,imon,id,istat,ibk,ierr,ip
  REAL(RL) :: sod,dt
  CHARACTER(LEN_PRN) :: prn
  CHARACTER(LEN_STRING) :: line,ccmd

  DATA ibk/0/
  SAVE ibk

  !*
  ! The function called
  !!------------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: system
  INTEGER(IT) :: modified_julday
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!----------------------------

  !! estimate new wide-lane FCBs
  INQUIRE(FILE=CKF.flnmw,EXIST=lexist)
  IF (lexist .EQ. .TRUE.) THEN
    IF (CKF.imw .NE. -1) THEN
      CLOSE(CKF.imw)
      CKF.imw=-1
    END IF
    istat=system('fcbwl '//CKF.flnmw)
    IF (istat .EQ. 0) THEN
      INQUIRE(FILE='fcb_'//CKF.flnmw(4:LEN_TRIM(CKF.flnmw))//'_wl',EXIST=lexist)
      IF (lexist) THEN
        INQUIRE(FILE='fcb_wide',EXIST=lexist)
        IF (lexist .EQ. .TRUE.) THEN
          ibk=ibk+1
          WRITE(ccmd,'(a,i3.3)') 'mv -f fcb_wide fcb_wide_',ibk
          istat=system(ccmd)
        END IF
        istat=system('mv -f fcb_'//CKF.flnmw(4:LEN_TRIM(CKF.flnmw))//'_wl fcb_wide')
      END IF
      istat=system('rm -f '//CKF.flnmw)
    END IF
  END IF

  !! read widelane FCBs
  INQUIRE(FILE='fcb_wide',EXIST=lexist)
  IF (lexist) THEN
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE='fcb_wide',STATUS='OLD')
    line=' '
    DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
      READ(lfn,'(A)') line
      IF (INDEX(line,'TIME END') .NE. 0) THEN
        READ(line,*) mjd,sod
      END IF
    END DO
    dt=timdif(CKF.mjd,CKF.sod,mjd,sod)
    IF (dt .LT. 86400.d0*5.d0) THEN
      UPD.wfcb=10.d0
      UPD.wsl =10.d0
      DO WHILE(.TRUE.)
        READ(lfn,'(A)',IOSTAT=ierr) line
        IF (ierr .NE. 0) EXIT
        IF (line(1:1) .NE. '*') CYCLE
        READ(line(2:),*) prn
        ip=pointer_string(CKF.nprn,CKF.cprn,prn)
        IF (ip .EQ. 0) CYCLE
        READ(line(5:),*) UPD.wfcb(ip),UPD.wsl(ip)
      END DO
      bmjd=CKF.mjd
      bsod=CKF.sod
      WRITE(OUTPUT_UNIT,'(A,I5,F9.1)') '%%%MESSAGE: Wide-lane FCBs updated to ',mjd,sod
      WRITE(OUTPUT_UNIT,'(A,I5,F9.1)') '             processing is started at ',bmjd,bsod
    END IF
    CLOSE(lfn)
  ELSE
    bmjd=CKF.mjd
    bsod=CKF.sod
  END IF

  RETURN

END SUBROUTINE
