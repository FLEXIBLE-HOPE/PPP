!*
SUBROUTINE bds_s2f_epoch(mjd,sod,cprn,lfix)
!!
!*
USE tables
IMPLICIT NONE

!*
! The arguments
!!------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod,beta,u
CHARACTER(LEN=*) :: cprn
LOGICAL(LG) :: lfix,ls2f

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: lfn,lfnt
  INTEGER(IT) :: iy(2),im(2),id(2)
  INTEGER(IT) :: ih(2),imin(2)
  REAL(RL) :: sec(2),mjdt(2),b(2),mu(2)

  CHARACTER(LEN_STRING) :: line

  LOGICAL(LG) :: lfirst,lexist,lfind
  DATA lfirst,lexist /.TRUE.,.TRUE./

  SAVE lfirst,lexist,lfn

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!-----------------------------

  lfix=.FALSE.

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    line=f_tablefilename('bdss2f')
    INQUIRE(FILE=line,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) RETURN

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=line)

  END IF

  IF (lexist .EQ. .FALSE.) RETURN

  REWIND(lfn)
  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=100) line
    IF (line(1:1) .NE. ' ') CYCLE
    IF (line(2:4) .NE. cprn(1:3)) CYCLE
    READ(line(5:),*) iy(1),im(1),id(1),ih(1),imin(1),sec(1),b(1),mu(1),iy(2),im(2),id(2),ih(2),imin(2),sec(2),b(2),mu(2)
    mjdt(1)=modified_julday(id(1),im(1),iy(1))+ih(1)/24.d0+imin(1)/1440.d0+sec(1)/86400.d0
    mjdt(2)=modified_julday(id(2),im(2),iy(2))+ih(2)/24.d0+imin(2)/1440.d0+sec(2)/86400.d0
    IF (mjd+sod/86400.d0.GE.mjdt(1)-MAXWND/86400.d0 .AND. mjd+sod/86400.d0.LT.mjdt(2)+MAXWND/86400.d0) THEN
      lfix=.TRUE.
      EXIT
    END IF
  END DO

100 CONTINUE

  RETURN

ENTRY bds_s2f_out(ls2f,mjd,sod,cprn,beta,u)

  IF (lexist .EQ. .FALSE.) THEN
    lexist=.TRUE.
    line=f_tablefilename('bdss2f')
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=line)
    WRITE(lfn,'(A)') '#CPRN YYYY MM DD HH MI SECD  BETA   MU    YYYY MM DD HH MI SECD  BETA   MU'
  END IF

  lfnt=get_valid_unit(10)
  OPEN(UNIT=lfnt,STATUS='SCRATCH')

  lfind=.FALSE.
  REWIND(lfn)
  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=200) line
    IF (line(1:1) .EQ. ' ') THEN
      IF (line(2:4) .EQ. cprn(1:3)) THEN
        READ(line(5:),*) iy(1),im(1),id(1),ih(1),imin(1),sec(1),b(1),mu(1),iy(2),im(2),id(2),ih(2),imin(2),sec(2),b(2),mu(2)
        IF (iy(1).EQ.9999 .AND. ls2f.EQ..TRUE.) THEN
          lfind=.TRUE.
          CALL mjd2date(mjd,sod,iy(1),im(1),id(1),ih(1),imin(1),sec(1))
          WRITE(line,'(1X,A3,2(1X,I4,4I3,F5.1,2F7.1))') cprn,iy(1),im(1),id(1),ih(1),imin(1),sec(1),beta,u, &
                                                             iy(2),im(2),id(2),ih(2),imin(2),sec(2),b(2),mu(2)
        END IF
        IF (iy(2).EQ.9999 .AND. ls2f.EQ..FALSE.) THEN
          lfind=.TRUE.
          CALL mjd2date(mjd,sod,iy(2),im(2),id(2),ih(2),imin(2),sec(2))
          WRITE(line,'(1X,A3,2(1X,I4,4I3,F5.1,2F7.1))') cprn,iy(1),im(1),id(1),ih(1),imin(1),sec(1),b(1),mu(1), &
                                                             iy(2),im(2),id(2),ih(2),imin(2),sec(2),beta,u
        END IF
      END IF
    END IF
    WRITE(lfnt,'(A)') TRIM(line)
  END DO

200 CONTINUE

  IF (lfind .EQ. .FALSE.) THEN
    IF (ls2f .EQ. .TRUE.) THEN
      CALL mjd2date(mjd,sod,iy(1),im(1),id(1),ih(1),imin(1),sec(1))
      WRITE(lfnt,'(1X,A3,2(1X,I4,4I3,F5.1,2F7.1))') cprn,iy(1),im(1),id(1),ih(1),imin(1),sec(1),beta,u, &
                                                         9999,1,1,1,1,0.0,0.0,0.0
    ELSE
      CALL mjd2date(mjd,sod,iy(2),im(2),id(2),ih(2),imin(2),sec(2))
      WRITE(lfnt,'(1X,A3,2(1X,I4,4I3,F5.1,2F7.1))') cprn,9999,1,1,1,1,0.0,0.0,0.0, &
                                                         iy(2),im(2),id(2),ih(2),imin(2),sec(2),beta,u
    END IF
  END IF

  REWIND(lfn)
  REWIND(lfnt)
  DO WHILE(.TRUE.)
    READ(lfnt,'(A)',END=300) line
    WRITE(lfn,'(A)') TRIM(line)
  END DO

300 CLOSE(lfnt)

  RETURN

END SUBROUTINE
