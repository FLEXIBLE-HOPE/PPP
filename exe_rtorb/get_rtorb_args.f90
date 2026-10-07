!*
SUBROUTINE get_rtorb_args(lpost,erpfile,CKF)
!!
!*
USE const
USE orbit
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
CHARACTER(LEN=*) erpfile
TYPE(ORBHDR) :: CKF
LOGICAL(LG) :: lpost

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: nargs,lfn,i,ierr
  REAL(RL) :: mjderp0,mjderp1
  CHARACTER(LEN_STRING) :: line = ''
  CHARACTER(LEN_FILENAME) :: sp3file = ''

  TYPE(T_FILETABLE) :: FT

  !*
  ! The functions called
  !!-----------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  ! Get the GNSS satellite information as well as start and end time accoring the erpfile

  !! read command arguments
  nargs=iargc()
  lpost=.FALSE.
  IF (nargs .NE. 3) THEN
    WRITE(OUTPUT_UNIT,'(A)') ' Position And Navigation Data Analyst (@PANDA) software '
    WRITE(OUTPUT_UNIT,'(A)') TRIM(VERSION)
    WRITE(OUTPUT_UNIT,'(A)') ' TOOLS: Orbit format transformation for Multi-GNSS systems'
    WRITE(OUTPUT_UNIT,'(A)') ' (including GPS, GLONASS, GALILEO, BEIDOU, QZSS, and so on)'
    WRITE(OUTPUT_UNIT,'(A)') ' '
    WRITE(OUTPUT_UNIT,'(A)') ' USAGE: rtorb file_table sp3file [erpfile] [-post]'
    CALL exit(1)
  END IF
  CALL getarg(1,sp3file)
  CALL read_filetable(sp3file,FT)

  CALL getarg(2,sp3file)
  CALL getarg(3,erpfile)
  IF (erpfile(1:5) .EQ. '-post') THEN
    lpost=.TRUE.
  END IF

  !! time span for the erp
  IF (lpost .EQ. .FALSE.) THEN
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=erpfile,STATUS='OLD',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(2A)') '***ERROR(get_rtorb_args): open ', TRIM(erpfile)
      CALL exit(1)
    END IF

    line=' '
    DO WHILE(INDEX(line,'MJD').EQ.0 .OR. INDEX(line,'UT1').EQ.0 .OR. &
             INDEX(line,'UTC').EQ.0 .OR. INDEX(line,'LOD').EQ.0)
      READ(lfn,'(A)',END=50) line
    END DO
    READ(lfn,'(A)') line
    mjderp0=0.d0
    mjderp1=0.d0
    DO WHILE(.TRUE.)
      READ(lfn,'(A)',END=50) line
      IF (mjderp0 .EQ. 0.d0) THEN
        READ(line,*) mjderp0
      END IF
      READ(line,*) mjderp1
    END DO
50  CLOSE(lfn)
  END IF

  !! read sp3 header
  CALL rdsp3h(sp3file,CKF.mjd0,CKF.sod0,CKF.mjd1,CKF.sod1,CKF.dintv,CKF.nprn,CKF.cprn)
  IF (lpost .EQ. .TRUE.) THEN
    mjderp0=CKF.mjd0+CKF.sod0/86400.d0
    mjderp1=CKF.mjd1+CKF.sod1/86400.d0
  END IF
  IF (CKF.mjd0+CKF.sod0/86400.d0 .LE. mjderp0) THEN
    CKF.mjd0 =INT(mjderp0)
    CKF.sod0=(mjderp0-CKF.mjd0)*86400.d0
  END IF
  IF (CKF.mjd1+CKF.sod1/86400.d0 .GT. mjderp1) THEN
    CKF.mjd1 =INT(mjderp1)
    CKF.sod1=(mjderp1-CKF.mjd1)*86400.d0
  END IF
  CKF.nequ=3
  CKF.rmjd=CKF.mjd0
  CKF.rsod=CKF.sod0

  RETURN

END SUBROUTINE
