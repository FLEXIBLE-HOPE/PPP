!*
SUBROUTINE read_ics(flnics,SAT,ICS)
!!
!*
USE orbit
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
CHARACTER(LEN=*) :: flnics
TYPE(SATE) :: SAT
TYPE(SATEPAR) :: ICS

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: lfn
  INTEGER(IT) :: i,j,ntot,npar
  DATA lfn /0/

  CHARACTER(LEN_PRN) :: cprn
  CHARACTER(LEN_PARNAME) :: pname
  CHARACTER(LEN_STRING) :: line

  LOGICAL(LG) :: lfound

  SAVE lfn

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Satrt the executable code
  !!---------------------------

  IF (lfn .EQ. 0) THEN
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnics,STATUS='OLD')
  END IF

  REWIND(lfn)
  READ(lfn,'(A/A)') line,line
  DO WHILE(.TRUE.)
    READ(lfn,'(A)',ERR=100,END=200) line
    IF (INDEX(line,'END of FILE') .NE. 0) GOTO 200
    READ(line(11:),*) cprn
    IF (SAT.cprn .EQ. cprn) EXIT
    DO WHILE (INDEX(line,'END of SAT') .EQ. 0)
      READ(lfn,'(A)') line
    END DO
  END DO

  ! Get initial values and time limit
  npar=SAT.npar
  ntot=0
  ICS.npwc=0
  READ(lfn,'(A)',END=200) line
  DO WHILE(INDEX(line,'END of SAT' ) .EQ. 0)
    i=INDEX(line,' ')
    READ(line(1:i-1),'(A)') pname
    lfound=.FALSE.
    DO j=1, SAT.npar
      IF (INDEX(SAT.pname(j),TRIM(pname)) .NE. 0) THEN
        lfound=.TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT.lfound .AND. npar.EQ.0) THEN
      SAT.npar=SAT.npar+1
      SAT.pname(SAT.npar)=TRIM(pname)
      j=SAT.npar
      lfound=.TRUE.
    END IF

    IF (lfound .EQ. .TRUE.) THEN
      ntot=ntot+1
      ICS.npwc(j)=ICS.npwc(j)+1
      READ(line(i+1:),*,ERR=100) ICS.val(ntot),ICS.ptime(1,ntot),ICS.ptime(2,ntot)
    ELSE
      WRITE(OUTPUT_UNIT,'(A)') '###WARNING(read_ics): this ics is not in requested list '//TRIM(pname)
    END IF
    READ(lfn,'(A)',END=200) line
  END DO

  DO i=1, SAT.npar
    IF (ICS.npwc(i) .EQ. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ics): requested ics not in ics-file '//TRIM(SAT.pname(i))
      CALL exit(1)
    END IF
  END DO

  ICS.npar=ntot

  RETURN

200 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ics): there is no ics for satellite '//SAT.cprn
  CALL exit(1)

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ics): read ics values'
  CALL exit(1)

END SUBROUTINE
