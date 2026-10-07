!*
SUBROUTINE read_sat_asc(prn,flnasc,mjd,qm,lfind)
!!
!! READ THE ADVANCED STERN CAMMERA DATA OF SATELLITES IN PANDA FORMAT
!! IT IS QUERTIRION FROM SPACECFAFT SYSTEM TO INTERTIAL SYSTEM
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
CHARACTER(LEN_PRN) :: prn
CHARACTER(LEN=*) :: flnasc
REAL(RL) :: mjd, qm(1:*)
LOGICAL(LG) :: lfind

TYPE ASCTAB
  LOGICAL(LG) :: lfirst=.TRUE.
  LOGICAL(LG) :: lexist=.FALSE.
  INTEGER(IT) :: lfn,id0,id1
  REAL(RL) :: mjd0,mjd1
  REAL(RL) :: mjds,mjde
  REAL(RL) :: q0(4),q1(4)
END TYPE

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,j,ierr
  INTEGER(IT) :: isat,nprn

  REAL(RL) :: sod,alpha

  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  DATA nprn /0/

  TYPE(ASCTAB) :: ASC(MAXSAT)

  SAVE nprn,cprn,ASC

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  isat=pointer_string(nprn,cprn,prn)
  IF (isat .EQ. 0) THEN
    nprn=nprn+1
    cprn(nprn)=prn
    isat=nprn
  END IF

  IF (ASC(isat).lfirst .EQ. .TRUE.) THEN
    ASC(isat).lfirst=.FALSE.

    INQUIRE(FILE=flnasc,EXIST=ASC(isat).lexist)
    IF (ASC(isat).lexist .EQ. .TRUE.) THEN
      ASC(isat).lfn=get_valid_unit(10)
      OPEN(UNIT=ASC(isat).lfn,FILE=flnasc,STATUS='OLD',IOSTAT=ierr)
      IF (ierr .NE. 0) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): open file '//TRIM(flnasc)
        CALL exit(1)
      END IF

      READ(ASC(isat).lfn,'(A/A/A)',END=100) line,line,line
      IF (INDEX(line,'%% Start and stop time :') .EQ. 0) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): no information for begin and end time'
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): '//TRIM(line)
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): '//TRIM(flnasc)
        CALL exit(1)
      END IF

      READ(line(25:),*,ERR=200) i,ASC(isat).mjds,j,ASC(isat).mjde
      ASC(isat).mjds=i+ASC(isat).mjds/86400.d0
      ASC(isat).mjde=j+ASC(isat).mjde/86400.d0

      READ(ASC(isat).lfn,*) j,sod,(ASC(isat).q0(i),i=1,4),ASC(isat).id0
      ASC(isat).mjd0=j+sod/86400.d0
      READ(ASC(isat).lfn,*) j,sod,(ASC(isat).q1(i),i=1,4),ASC(isat).id1
      ASC(isat).mjd1=j+sod/86400.d0
    END IF
  END IF

  ! If there is no attitude file, the nominal will be used
  IF (ASC(isat).lexist .EQ. .FALSE.) THEN
    lfind=.FALSE.
    RETURN
  END IF

  ! check data span
  IF (mjd.LT.ASC(isat).mjds .OR. mjd.GT.ASC(isat).mjde) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): the requested epoch is out of the input file'
    CALL exit(1)
  END IF

  ! check the current record
  lfind=.FALSE.
  DO WHILE(lfind .EQ. .FALSE.)
    IF (mjd.GE.ASC(isat).mjd0-MAXWND/86400.d0 .AND. mjd.LE.ASC(isat).mjd1+MAXWND/86400.d0) then
      lfind=.TRUE.
    ELSE IF (mjd .GT. ASC(isat).mjd0-MAXWND/86400.d0) THEN
      ASC(isat).mjd0=ASC(isat).mjd1
      ASC(isat).id0=ASC(isat).id1
      DO i=1,4
        ASC(isat).q0(i)=ASC(isat).q1(i)
      END DO
      READ(ASC(isat).lfn,*,ERR=200,END=100) j,sod,(ASC(isat).q1(i),i=1,4),ASC(isat).id1
      ASC(isat).mjd1=j+sod/86400.d0
    ELSE
      ASC(isat).mjd1=ASC(isat).mjd0
      ASC(isat).id1=ASC(isat).id0
      DO i=1,4
        ASC(isat).q1(i)=ASC(isat).q0(i)
      END DO
      BACKSPACE(ASC(isat).lfn)
      BACKSPACE(ASC(isat).lfn)
      READ(ASC(isat).lfn,*,ERR=200,END=100) j,sod,(ASC(isat).q0(i),i=1,4),ASC(isat).id0
      ASC(isat).mjd0=j+sod/86400.d0
    END IF
  END DO

  IF (lfind .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): the requested epoch is out of the input file'
    CALL exit(1)
  END IF

  ! linear interpolation
  alpha=(mjd-ASC(isat).mjd0)/(ASC(isat).mjd1-ASC(isat).mjd0)
  DO i=1, 4
    qm(i)=ASC(isat).q0(i)+alpha*(ASC(isat).q1(i)-ASC(isat).q0(i))
  END DO

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): end of file'
  CALL exit(1)

200 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_asc): read file '//TRIM(flnasc)
  CALL exit(1)

END SUBROUTINE
