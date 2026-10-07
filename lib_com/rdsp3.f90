!*
SUBROUTINE rdsp3h(sp3file,jd0,sod0,jd1,sod1,dintv,nprn,cprn)
!!
!! purpose   : read IGS sp3 orbit file
!!
!! parameters: fln -- sp3 file name
!!             jd0,sod0,jd1,sod1,dintv -- start and stop time and interval
!!             nprn,prn -- number of satellites and satellite list
!!             jd,sod -- request epoch time
!!             xsp3   -- satellite position and velocity
!!
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: jd, jd1, nprn
CHARACTER(LEN_PRN) :: cprn(1:*)
REAL(RL) :: sod,sod1,dintv,xsp3(6,MAXSAT)
CHARACTER(LEN=*) :: sp3file

  !*
  ! The local variables
  !!-----------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: iy,im,id,ih,imin,nwk,mjd,i,j,k,l,lfn,ierr,iflag
  INTEGER(IT) :: nepo,jd0,nprn0,jdx
  CHARACTER(LEN_PRN) :: cprn0(MAXSAT*2), prn
  REAL(RL) :: fmjd,sec,sod0,x,y,z,t,dt,sodx
  CHARACTER(LEN_STRING) :: line

  DATA lfirst/.TRUE./
  SAVE lfirst,nprn0,lfn

  !*
  ! The function called
  !!-----------------------
  INTEGER(IT) :: modified_julday,pointer_string,get_valid_unit
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!--------------------------


  lfn=0
  IF (lfirst .EQ. .TRUE.) THEN
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=sp3file,STATUS='OLD',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(rdsp3h): open file '//TRIM(sp3file)
      CALL exit(1)
    END IF
    jd0=0
    sod0=0.d0

    !! scan sp3 file to load all gps satellites
    DO WHILE(.TRUE.)

      !! second line with interval
      READ(lfn,'(A)',END=33) line
      IF (line(1:2).NE.'#c' .AND. line(1:2).NE.'#d') CYCLE
      READ(lfn,'(3X,I4,17X,F14.8,1X,I5,1X,F15.13)') nwk,dintv,mjd,fmjd

      !! third to 12th lines with satellite prns
      READ(lfn,'(3X,I3,3X,17A3,9(/,9X,17A3))') nprn0, (cprn0(i),i=1,MAXSAT*2)
      IF(nprn0 .GT. MAXSAT) THEN
        WRITE(ERROR_UNIT,'(A,I3)') '***ERROR(rdsp3h): too many satellites',nprn0
        CALL exit(1)
      END IF

      nprn = nprn0
      cprn(1:nprn0) = cprn0(1:nprn0)
      DO i=1, nprn
        IF (cprn(i)(1:1).EQ.' ') cprn(i)(1:1)='G'
      END DO
    END DO

33  CONTINUE

    REWIND(lfn)
    DO WHILE(.TRUE.)
      READ(lfn,'(A)',END=5) line
      IF (line(1:1) .EQ.'*') THEN
        READ(line(2:),*) iy,im,id,ih,imin,sec
        IF (jd0 .EQ. 0) THEN
          jd0=modified_julday(id,im,iy)
          sod0=ih*3600.d0+imin*60.d0+sec
        END IF
      END IF
    END DO
5   REWIND(lfn)
    jd1=modified_julday(id,im,iy)
    sod1=ih*3600.d0+imin*60.d0+sec
    lfirst=.FALSE.

    RETURN

  END IF

  !! read one record in SP3 file
ENTRY rdsp3i(jd,sod,nprn,cprn,xsp3,iflag)

  IF (lfn .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rdsp3i): sp3file not open'
    CALL exit(1)
  END IF

  !! check if all requested satellites in sp3-file
  iflag=1
  DO j=1,nprn0
    DO i=1,3
      xsp3(i,j)=1.d15
    END DO
  END DO

  !! read the required record
  k=0
  line=' '
  DO WHILE(line(1:1) .NE. '*')
    READ(lfn,'(a)',end=100) line
  END DO
10 CONTINUE

  READ(line(2:),*) iy,im,id,ih,imin,sec
  CALL yr2year(iy)
  jdx=modified_julday(id,im,iy)
  sodx=ih*3600.d0+imin*60.d0+sec

  dt=timdif(jdx,sodx,jd,sod)
  IF (DABS(dt) .LT. MAXWND) THEN
    !jd=jdx
    !sod=sodx
    iflag=0
    line=' '
    DO WHILE(line(1:1) .NE. '*')
      READ(lfn,'(A)',END=100) line
      IF (line(1:1) .NE. 'P') CYCLE
      READ(line(2:),'(A3)') prn
      IF (prn(1:1) .EQ. ' ') prn(1:1)='G'
      READ(line(5:),*) x,y,z

      !! nprn <= 0, get all the satellite in the file, otherwise,
      !! get data of the requested satellites
      IF (nprn .GT. 0) THEN
        k=pointer_string(nprn,cprn,prn)
      ELSE
        k=k+1
        cprn(k)=prn
      END IF

      !! k == 0, this satellite is not in the requested list
      IF (k .NE. 0) THEN
        xsp3(1,k)=x
        xsp3(2,k)=y
        xsp3(3,k)=z
      END IF
    END DO
  ELSE IF (dt .LE. -MAXWND) THEN
    line=' '
    DO WHILE(line(1:1) .NE. '*')
      READ(lfn,'(A)',END=100) line
    END DO
    GOTO 10
  END IF
  BACKSPACE(lfn)

100 CONTINUE
  IF (nprn .LE. 0) nprn=k

  !! check if all prns have data
  DO i=1, nprn
    IF (dabs(xsp3(1,i))+dabs(xsp3(2,i))+dabs(xsp3(3,i)) .EQ. 0.d0) THEN
      xsp3(1,i)=1.d15
      xsp3(2,i)=1.d15
      xsp3(3,i)=1.d15
    END IF
  END DO

  RETURN

!! close file
ENTRY rdsp3c()

  CLOSE(lfn)
  lfn=0
  lfirst=.TRUE.

  RETURN

!! read one record in SP3 file
ENTRY rdsp3i_kin(jd,sod,nprn,cprn,xsp3,iflag)

  IF (lfn .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rdsp3i): sp3file not open'
    CALL exit(1)
  END IF

  !! check if all requested satellites in sp3-file
  iflag=1
  DO j=1,nprn0
    DO i=1,3
      xsp3(i,j)=1.d15
    END DO
  END DO

  !! read the required record
  k=0
  line=' '
  DO WHILE(line(1:1) .NE. '*')
    READ(lfn,'(a)',end=200) line
  END DO
20 CONTINUE

  READ(line(2:),*) iy,im,id,ih,imin,sec
  CALL yr2year(iy)
  jdx=modified_julday(id,im,iy)
  sodx=ih*3600.d0+imin*60.d0+sec

  dt=timdif(jdx,sodx,jd,sod)
  IF (DABS(dt) .LT. MAXWND) THEN
    jd=jdx
    sod=sodx
    iflag=0
    line=' '
    DO WHILE(line(1:1) .NE. '*')
      READ(lfn,'(A)',END=200) line
      IF (line(1:1) .NE. 'P') CYCLE
      READ(line(2:),'(A3)') prn
      IF (prn(1:1) .EQ. ' ') prn(1:1)='G'
      READ(line(5:),*) x,y,z

      !! nprn <= 0, get all the satellite in the file, otherwise,
      !! get data of the requested satellites
      IF (nprn .GT. 0) THEN
        k=pointer_string(nprn,cprn,prn)
      ELSE
        k=k+1
        cprn(k)=prn
      END IF

      !! k == 0, this satellite is not in the requested list
      IF (k .NE. 0) THEN
        xsp3(1,k)=x
        xsp3(2,k)=y
        xsp3(3,k)=z
      END IF
    END DO
  ELSE IF (dt .LE. -MAXWND) THEN
    line=' '
    DO WHILE(line(1:1) .NE. '*')
      READ(lfn,'(A)',END=200) line
    END DO
    GOTO 20
  END IF
  BACKSPACE(lfn)

200 CONTINUE
  IF (nprn .LE. 0) nprn=k

  !! check if all prns have data
  DO i=1, nprn
    IF (dabs(xsp3(1,i))+dabs(xsp3(2,i))+dabs(xsp3(3,i)) .EQ. 0.d0) THEN
      xsp3(1,i)=1.d15
      xsp3(2,i)=1.d15
      xsp3(3,i)=1.d15
    END IF
  END DO

  RETURN

END SUBROUTINE
