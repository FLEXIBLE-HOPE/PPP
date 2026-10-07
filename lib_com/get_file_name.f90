!*
SUBROUTINE get_file_name(ldefined,keyword,param_list,iyear,imonth,iday,ihour,name)
!!
!! TO GET THE FILE NAME ACCORDING PREDEFINED IN THE PANDA_FILE_NAME FILE
!!
!*
USE par
USE tables
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
LOGICAL(LG) :: ldefined
INTEGER(IT) :: iyear,imonth,iday,ihour
CHARACTER(LEN=*) :: keyword,param_list,name

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  TYPE(T_FileTable) :: FT
  LOGICAL(LG) :: lfirst,lexist,lfound
  CHARACTER(LEN_STRING) :: form,cfilename
  INTEGER(IT) :: i,j,k,l,mjd,week,wd,idoy,ierr,lun

  INTEGER(IT) :: nfile,MAXFILE,npar,nword
  PARAMETER(MAXFILE=100)
  CHARACTER(LEN=10) :: fname(MAXFILE)
  CHARACTER(LEN=30) :: param(40)
  CHARACTER(LEN=256) :: word(40)
  CHARACTER(LEN=256) :: fform(MAXFILE)
  DATA lfirst /.TRUE./

  SAVE lfirst,nfile,fname,fform

  !*
  ! THE FUNCTIONS CALLED
  !!-----------------------------
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: get_valid_unit

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  ! TO GET PANDA_FILE_NAME IN FILE_TABLE
  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    cfilename=f_tableFileName('pndfln')
    INQUIRE(FILE=cfilename,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(get_file_name): '//TRIM(cfilename)//' is not exist'
      CALL exit(1)
    END IF

    lun=get_valid_unit(10)
    OPEN(UNIT=lun,FILE=cfilename)
    nfile=0
    DO WHILE (1 .EQ. 1)
      READ(lun,'(A10,1X,A)',ERR=100,END=100) fname(nfile+1),fform(nfile+1)
      nfile=nfile+1
      CALL left_justify_string(fform(nfile+1))
      IF (nfile .GT. MAXFILE) THEN
        WRITE(ERROR_UNIT,'(A,2I5)') '***ERROR(get_file_name): too manay files defined in panda_file_name', nfile, MAXFILE
        CALL exit(1)
      END IF
    END DO
100 CONTINUE
    CLOSE(lun)
  END IF

  ! FIRSTLY FIND OUT THE FORMAT
  IF (ldefined .EQ. .TRUE.) THEN
    form=name
  ELSE
    form=' '
    DO i=1, nfile
      IF(fname(i)(1:LEN_TRIM(fname(i))) .EQ. keyword(1:LEN_TRIM(keyword))) THEN
        form=fform(i)
        EXIT
      END IF
    END DO
    IF (form(1:1) .EQ. ' ') THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(get_file_name): not defined in panda_file_name '//TRIM(keyword)
      CALL exit(1)
    END IF
  END IF

  ! TO FIND OUT THE VARIABLES IN THE FORMAT
  CALL split_string(.FALSE.,form,' ',' ','-',nword,word)
  CALL split_string(.TRUE.,param_list,' ',' ',':',npar,param)

  ! NAME DIFINATION
  mjd=modified_julday(iday,imonth,iyear)
  CALL gpsweek(iday,imonth,iyear,week,wd)
  CALL mjd2doy(mjd,iyear,idoy)
  WRITE(param(npar+1),'(A,I4.4)') 'YYYY=',iyear
  WRITE(param(npar+2),'(A,I3.3)') 'DDD=',idoy
  WRITE(param(npar+3),'(A,I2.2)') 'HH=',ihour
  WRITE(param(npar+4),'(A,I4.4)') 'WWWW=',week
  WRITE(param(npar+4),'(A,I4.4)') 'GPSW=',week
  WRITE(param(npar+5),'(A,I5.5)') 'WWWWD=',week*10+wd
  WRITE(param(npar+6),'(A,I2.2)') 'YY=',iyear-INT(iyear/100)*100
  WRITE(param(npar+7),'(A,I1.1)') 'Y=',iyear-INT(iyear/10)*10
  WRITE(param(npar+8),'(A,I5.5)') 'MJD=',mjd
  npar=npar+8

  name=' '
  DO i=1, nword
    IF(LEN_TRIM(word(i)) .EQ. 0) CYCLE
    IF(INT(i/2)*2 .EQ. i) THEN
      l=LEN_TRIM(word(i))
      lfound=.FALSE.
      DO j=1, npar
        IF(word(i)(1:l)//'=' .EQ. param(j)(1:l+1)) THEN
          k=LEN_TRIM(name)
          name(k+1:)=param(j)(l+2:)
          lfound=.TRUE.
          EXIT
        END IF
      END DO
      IF(lfound .EQ. .FALSE.) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(get_file_name): variable not defined '//word(i)(1:l)
        CALL exit(1)
      END IF
    ELSE
      l=LEN_TRIM(name)
      name(l+1:)=word(i)
    END IF
  END DO

  RETURN

END SUBROUTINE
