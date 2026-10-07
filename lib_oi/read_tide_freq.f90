!*
SUBROUTINE read_tide_freq(conv,estf)
!!
!*
USE orbit
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------------
CHARACTER(LEN=*) :: conv
TYPE(SOLID_TIDE_FREQ) :: estf(0:2)

  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: i,k,lfn,ierr

  CHARACTER(LEN_STRING) :: line
  LOGICAL(LG) :: lexist

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  SELECT CASE(TRIM(conv))
    CASE('IERS1996')
      line=f_tablefilename('estf96')
    CASE('IERS2003','IERS2010')
      line=f_tablefilename('estf10')
    CASE DEFAULT
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): unknon conventions '//TRIM(conv)
      CALL exit(1)
  END SELECT

  INQUIRE(FILE=line,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): '//TRIM(line)//' is not exist'
    CALL exit(1)
  END IF

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=line,STATUS='OLD')

  line=''
  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=100,ERR=200) line
    IF (line(1:10) .EQ. '----------') THEN
      BACKSPACE(lfn)
      BACKSPACE(lfn)
      READ(lfn,'(A)',END=100,ERR=200) line
      SELECT CASE(line(1:3))
        CASE('K21')
          k=1
        CASE('K20')
          k=0
        CASE('K22')
          k=2
        CASE DEFAULT
          WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): unknown frequention correction for love number '//TRIM(line(1:3))
          CALL exit(1)
      END SELECT

      i=INDEX(line,'NR=')
      IF (i .EQ. 0) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): there is no number for tides'
        CALL exit(1)
      END IF
      READ(line(i+3:),*) estf(k).n
      ALLOCATE(estf(k).etf(6,estf(k).n),STAT=ierr)
      IF (ierr .NE. 0) GOTO 300
      ALLOCATE(estf(k).kf(2,estf(k).n),STAT=ierr)
      IF (ierr .NE. 0) GOTO 300
      ALLOCATE(estf(k).amp(2,estf(k).n),STAT=ierr)
      IF (ierr .NE. 0) GOTO 300
      READ(lfn,'(A)',END=100,ERR=200) line
      DO i=1, estf(k).n
        READ(lfn,'(A)',END=100,ERR=200) line
        READ(line,'(21X,6I3,15X,F6.0,F7.0,F6.1,F6.1)',ERR=200) estf(k).etf(1:6,i),estf(k).kf(1:2,i),estf(k).amp(1:2,i)
        estf(k).kf(1:2,i)=estf(k).kf(1:2,i)*1.d-5
      END DO
      READ(lfn,'(A)',END=100,ERR=200) line
    END IF
  END DO

100 CONTINUE

  CLOSE(lfn)

  RETURN

200 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): read '//TRIM(line)
  CALL exit(1)

300 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_tide_freq): memory allocation'
  CALL exit(1)

END SUBROUTINE
