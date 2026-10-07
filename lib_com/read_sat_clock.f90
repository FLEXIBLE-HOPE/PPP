!*
SUBROUTINE read_sat_clock(mjd,sod,CKF,SAT)
!!
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
  !!--------------------------
  LOGICAL(LG) :: lfirst,lexist
  REAL(RL) :: sec
  REAL(RL) :: coef(6),a0(2,MAXSAT)
  REAL(RL) :: a1(2,MAXSAT),a2(2,MAXSAT)
  REAL(RL) :: dt1,dt2,alpha

  INTEGER(IT) :: lfnclk,ierr
  INTEGER(IT) :: i,j,k,iy,imon,id,ih,im
  INTEGER(IT) :: mjdf(2), mjdx
  REAL(RL) :: sodf(2), sodx
  CHARACTER(LEN_STRING) :: line

  INTEGER(IT) :: nprn(2)
  CHARACTER(LEN_PRN) :: cprn(MAXSAT,2)
  DATA lfirst /.TRUE./
  SAVE lfirst,lfnclk,nprn,cprn,a0,a1,mjdf,sodf

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  DO k=1, CKF.nprn
    SAT(k).sclock=0.d0
  END DO

  IF (lfirst .EQ. .TRUE.) THEN

    lfirst=.FALSE.

    INQUIRE(FILE=CKF.flncck,EXIST=lexist)
    IF (lexist .EQ. .TRUE.) THEN
      lfnclk=get_valid_unit(10)
      OPEN(UNIT=lfnclk,FILE=CKF.flncck,STATUS='OLD')

      line=' '
      DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
        READ(lfnclk,'(a)',END=100) line
      END DO

    ELSE
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_clock): the clock file is not exist '//TRIM(CKF.flncck)
      CALL exit(1)
    END IF

    !! read two epoch clocks
    k=1
    mjdf=0
    nprn=0
    DO WHILE(k .LE. 2)
      line=' '
      DO WHILE(line(1:3) .NE. 'AS ')
        READ(lfnclk,'(A)',END=200,ERR=100) line
      END DO
      DO WHILE(.TRUE.)

        IF (line(1:3) .NE. 'AS ') THEN
          READ(lfnclk,'(A)',END=200,ERR=100) line
          CYCLE
        END IF

        coef=0.d0
        READ(line(7:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,coef(1),coef(2)

        !! check time tag
        CALL yr2year(iy)
        mjdx=modified_julday(id,imon,iy)
        sodx=ih*3600.d0+im*60.d0+sec

        IF (mjdf(k) .EQ. 0) THEN
          mjdf(k)=mjdx
          sodf(k)=sodx
        ELSE IF (timdif(mjdx,sodx,mjdf(k),sodf(k)) .GT. 0.d0) THEN
          BACKSPACE(lfnclk)
          GOTO 202
        END IF

        !! read continuing lines
        IF (j .EQ. 4) THEN
          READ(lfnclk,'(A)',END=200,ERR=100) line
          READ(line,*,IOSTAT=ierr) coef(3),coef(4)
        ELSE IF (j .EQ.6) THEN
          READ(lfnclk,'(A)',END=200,ERR=100) line
          READ(line,*,IOSTAT=ierr) coef(3),coef(4),coef(5),coef(6)
        END IF

        !! set prn and clocks
        nprn(k)=nprn(k)+1
        cprn(nprn(k),k)=line(4:6)
        a0(k,nprn(k))=coef(1)
        a1(k,nprn(k))=coef(3)
        a2(k,nprn(k))=coef(5)
        READ(lfnclk,'(A)',END=200,ERR=100) line
      END DO
202   k=k+1
    END DO

  END IF

  !! check time tag
10 dt1=timdif(mjd,sod,mjdf(1),sodf(1))
  dt2=timdif(mjd,sod,mjdf(2),sodf(2))
  IF (dt1 .LT. 0.d0)  THEN
    RETURN

  ELSE IF (dt2 .GT. 0.d0) THEN

    !! transfer clocks
    nprn(1)=nprn(2)
    DO i=1,nprn(2)
      cprn(i,1)=cprn(i,2)
      a0(1,i)=a0(2,i)
      a1(1,i)=a1(2,i)
      a2(1,i)=a2(2,i)
    END DO

    mjdf(1)=mjdf(2)
    sodf(1)=sodf(2)
    mjdf(2)=0
    nprn(2)=0

    !! read next epoch clocks
    line=' '
    DO WHILE(line(1:3) .NE. 'AS ')
      READ(lfnclk,'(A)',END=200,ERR=100) line
    END DO
    DO WHILE(.TRUE.)

      IF (line(1:3) .NE. 'AS ') THEN
        READ(lfnclk,'(A)',END=200,ERR=100) line
        CYCLE
      END IF

      coef=0.d0
      READ(line(9:),*,IOSTAT=ierr) iy,imon,id,ih,im,sec,j,coef(1),coef(2)

      !! check time tag
      CALL yr2year(iy)

      mjdx=modified_julday(id,imon,iy)
      sodx=ih*3600.d0+im*60.d0+sec
      IF (mjdf(2) .EQ. 0) THEN
        mjdf(2)=mjdx
        sodf(2)=sodx
      ELSE IF (timdif(mjdx,sodx,mjdf(2),sodf(2)) .GT. 0.d0) THEN
        BACKSPACE(lfnclk)
        EXIT
      END IF

      !! read continuing lines
      IF (j .EQ. 4) THEN
        READ(lfnclk,'(A)',END=200,ERR=100) line
        READ(line,*,IOSTAT=ierr) coef(3),coef(4)
      ELSE IF (j .EQ.6) THEN
        READ(lfnclk,'(A)',END=200,ERR=100) line
        READ(line,*,IOSTAT=ierr) coef(3),coef(4),coef(5),coef(6)
      END IF

      !! set prn and clocks
      nprn(2)=nprn(2)+1
      cprn(nprn(2),2)=line(4:6)
      a0(2,nprn(2))=coef(1)
      a1(2,nprn(2))=coef(3)
      a2(2,nprn(2))=coef(5)
      READ(lfnclk,'(A)',END=200,ERR=100) line
    END DO
    GOTO 10
  ELSE

    DO i=1, CKF.nprn
      j=pointer_string(nprn(1),cprn(1,1),CKF.cprn(i))
      k=pointer_string(nprn(2),cprn(1,2),CKF.cprn(i))

      !! interpolation
      IF (j.NE.0 .AND. k.NE.0) THEN
        alpha=(a0(2,k)-a0(1,j))/(timdif(mjdf(2),sodf(2),mjdf(1),sodf(1)))
        SAT(i).sclock=a0(1,j)+alpha*dt1
      ELSE IF (j.EQ.0 .AND. k.NE.0) THEN
        SAT(i).sclock=a0(2,k)
      END IF
    END DO

  END IF

  RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_sat_clock): read the clock file error '//TRIM(line)
  CALL exit(1)

200 WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(read_sat_clock): end of the clock file'

  DO i=1, CKF.nprn
    j=pointer_string(nprn(1),cprn(1,1),CKF.cprn(i))
    k=pointer_string(nprn(2),cprn(1,2),CKF.cprn(i))

    !! interpolation
    IF (j.NE.0 .AND. k.NE.0) THEN
      alpha=(a0(2,k)-a0(1,j))/(timdif(mjdf(2),sodf(2),mjdf(1),sodf(1)))
      SAT(i).sclock=a0(1,j)+alpha*dt1
    END IF
  END DO

  RETURN

ENTRY sat_clock_reset()

  lfirst=.TRUE.
  lexist=.FALSE.
  nprn=0
  cprn=''
  mjdf=0
  sodf=0.d0
  CLOSE(lfnclk)

  RETURN

END SUBROUTINE
