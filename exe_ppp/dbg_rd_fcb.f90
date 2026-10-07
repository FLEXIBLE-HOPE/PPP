!*
SUBROUTINE dbg_rd_fcb(fcbt,CKF,UPD)
!*
USE ckdctrl
USE ambiguity
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(CKDCFG) :: CKF
TYPE(FCB) :: UPD
CHARACTER(LEN=*) :: fcbt

  !*
  ! The local variables
  !!-----------------------------
  LOGICAL(LG) :: lfirst(5)

  INTEGER(IT) :: mjdf(2,5),mjdx,i,j,k
  INTEGER(IT) :: ifcb,lfn(5),isat,nsat,ierr
  REAL(RL) :: xfcb(MAXSAT,2,5),sfcb(MAXSAT,2,5)
  REAL(RL) :: sodf(2,5),sodx,dt1,dt2,dui(2),alpha
  CHARACTER(LEN_PRN) :: cprn
  CHARACTER(LEN_STRING) :: line
  LOGICAL(LG) :: lexist

  DATA lfirst /.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./
  SAVE lfirst,lfn,mjdf,sodf,xfcb,sfcb

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

  IF (fcbt(1:2) .EQ. 'WL') THEN
    ifcb=1
  ELSE IF (fcbt(1:2) .EQ. 'NL') THEN
    ifcb=2
  ELSE IF (fcbt(1:3) .EQ. 'EWL') THEN
    ifcb=3
  ELSE IF (fcbt(1:4) .EQ. 'EEWL') THEN
    ifcb=4
  ELSE IF (fcbt(1:4) .EQ. 'HEWL') THEN
    ifcb=5
  ELSE
    WRITE(ERROR_UNIT,'(A)') '####ERROR(dbg_rd_fcb): unknown type '//TRIM(fcbt)
    CALL exit(1)
  END IF

  !! open file
  IF (lfirst(ifcb)) THEN

    xfcb(:,:,ifcb)=10.d0
    sfcb(:,:,ifcb)=10.d0

    lfirst(ifcb)=.FALSE.

    !! open wide-lane FCB file
    lfn(ifcb)=0
    IF (fcbt(1:2) .EQ. 'WL') THEN
      INQUIRE(FILE=CKF.flnfwl,EXIST=lexist)
    ELSE IF (fcbt(1:2) .EQ. 'NL') THEN
      INQUIRE(FILE=CKF.flnfnl,EXIST=lexist)
    ELSE IF (fcbt(1:3) .EQ. 'EWL') THEN
      INQUIRE(FILE=CKF.flnfewl,EXIST=lexist)
    ELSE IF (fcbt(1:4) .EQ. 'EEWL') THEN
      INQUIRE(FILE=CKF.flnfeewl,EXIST=lexist)
    ELSE IF (fcbt(1:4) .EQ. 'HEWL') THEN
      INQUIRE(FILE=CKF.flnfhewl,EXIST=lexist)
    END IF
    IF (lexist .EQ. .TRUE.) THEN
      lfn(ifcb)=get_valid_unit(10)
      IF (fcbt(1:2) .EQ. 'WL') THEN
        OPEN(UNIT=lfn(ifcb),FILE=CKF.flnfwl,STATUS='OLD')
      ELSE IF (fcbt(1:2) .EQ. 'NL') THEN
        OPEN(UNIT=lfn(ifcb),FILE=CKF.flnfnl,STATUS='OLD')
      ELSE IF (fcbt(1:3) .EQ. 'EWL') THEN
        OPEN(UNIT=lfn(ifcb),FILE=CKF.flnfewl,STATUS='OLD')
      ELSE IF (fcbt(1:4) .EQ. 'EEWL') THEN
        OPEN(UNIT=lfn(ifcb),FILE=CKF.flnfeewl,STATUS='OLD')
      ELSE IF (fcbt(1:4) .EQ. 'HEWL') THEN
        OPEN(UNIT=lfn(ifcb),FILE=CKF.flnfhewl,STATUS='OLD')
      END IF
    ELSE
      WRITE(OUTPUT_UNIT,'(A)') '***WARNING(dbg_rd_fcb): the '//TRIM(fcbt)//'  FCB file is not exist'
      RETURN
    END IF

    !! read two epoch clocks
    k=1
    mjdf(1:2,ifcb)=0
    sodf(1:2,ifcb)=0.d0
    DO WHILE(k .LE. 2)
      line=' '
      DO WHILE(line(1:3) .NE. 'TIM')
        READ(lfn(ifcb),'(A)',END=200,ERR=100) line
      END DO

      READ(line,'(30X,I7,F10.2,I4)') mjdx,sodx,i
      IF (mjdf(k,ifcb) .EQ. 0) THEN
        mjdf(k,ifcb)=mjdx
        sodf(k,ifcb)=sodx
      ELSE IF (timdif(mjdx,sodx,mjdf(k,ifcb),sodf(k,ifcb)) .GT. 0.d0) THEN
        BACKSPACE(lfn(ifcb))
        GOTO 202
      END IF

      !! read continuing lines
      DO WHILE(.TRUE.)
        READ(lfn(ifcb),'(a)',IOSTAT=ierr) line
        IF (ierr .NE. 0) EXIT      ! end of file
        IF (line(1:3) .EQ. 'TIM') THEN
          BACKSPACE(lfn(ifcb))
          EXIT
        END IF
        line='   '//TRIM(line)
        nsat=INT(LEN_TRIM(line)/22)
        DO i=1, nsat
          READ(line((i-1)*22+1:),'(3X,A3,2F8.3)') cprn,dui(1:2)
          isat=pointer_string(CKF.nprn,CKF.cprn,cprn)
          IF (isat .NE. 0) THEN
            xfcb(isat,k,ifcb)=dui(1)
            sfcb(isat,k,ifcb)=dui(2)
          END IF
        END DO
      END DO
202   k=k+1
    END DO

  ELSE
    
    IF (lfn(ifcb) .EQ. 0) RETURN

  END IF

  !! check time tag
10 dt1=timdif(CKF.mjd,CKF.sod,mjdf(1,ifcb),sodf(1,ifcb))
  dt2=timdif(CKF.mjd,CKF.sod,mjdf(2,ifcb),sodf(2,ifcb))
  IF (dt1 .LT. 0.d0)  THEN
    RETURN

  ELSE IF (dt2 .GT. 0.d0) THEN

    !! transfer clocks
    DO i=1,CKF.nprn
      xfcb(i,1,ifcb)=xfcb(i,2,ifcb)
      sfcb(i,1,ifcb)=sfcb(i,2,ifcb)
    END DO

    mjdf(1,ifcb)=mjdf(2,ifcb)
    sodf(1,ifcb)=sodf(2,ifcb)
    mjdf(2,ifcb)=0

    !! read next epoch clocks
    line=' '
    DO WHILE(line(1:3) .NE. 'TIM')
      READ(lfn(ifcb),'(A)',END=200,ERR=100) line
    END DO

    READ(line,'(30X,I7,F10.2,I4)') mjdx,sodx,i
    IF (mjdf(2,ifcb) .EQ. 0) THEN
      mjdf(2,ifcb)=mjdx
      sodf(2,ifcb)=sodx
    ELSE IF (timdif(mjdx,sodx,mjdf(2,ifcb),sodf(2,ifcb)) .GT. 0.d0) THEN
      BACKSPACE(lfn(ifcb))
      GOTO 10
    END IF

    !! read continuing lines
    DO WHILE(.TRUE.)
      READ(lfn(ifcb),'(a)',IOSTAT=ierr) line
      IF (ierr .NE. 0) EXIT      ! end of file
      IF (line(1:3) .EQ. 'TIM') THEN
        BACKSPACE(lfn(ifcb))
        EXIT
      END IF
      line='   '//TRIM(line)
      nsat=INT(LEN_TRIM(line)/22)
      DO i=1, nsat
        READ(line((i-1)*22+1:),'(3X,A3,2F8.3)') cprn,dui(1:2)
        isat=pointer_string(CKF.nprn,CKF.cprn,cprn)
        IF (isat .NE. 0) THEN
          xfcb(isat,2,ifcb)=dui(1)
          sfcb(isat,2,ifcb)=dui(2)
        END IF
      END DO
    END DO
    GOTO 10
  ELSE

    DO i=1, CKF.nprn
      !! interpolation
      IF (xfcb(i,1,ifcb).NE.10.d0 .AND. xfcb(i,2,ifcb).NE.10.d0) THEN
        alpha=(xfcb(i,2,ifcb)-xfcb(i,1,ifcb))/(timdif(mjdf(2,ifcb),sodf(2,ifcb),mjdf(1,ifcb),sodf(1,ifcb)))
        IF (fcbt(1:2) .EQ. 'WL') THEN
          UPD.wfcb(i)=xfcb(i,1,ifcb)+alpha*dt1
          UPD.wsl(i)=sfcb(i,1,ifcb)
        ELSE IF (fcbt(1:2) .EQ. 'NL') THEN
          UPD.nfcb(i)=xfcb(i,1,ifcb)+alpha*dt1
          UPD.nsl(i)=sfcb(i,1,ifcb)
        ELSE IF (fcbt(1:3) .EQ. 'EWL') THEN
          UPD.ewfcb(i)=xfcb(i,1,ifcb)+alpha*dt1
          UPD.ewsl(i)=sfcb(i,1,ifcb)
        ELSE IF (fcbt(1:4) .EQ. 'EEWL') THEN
          UPD.eewfcb(i)=xfcb(i,1,ifcb)+alpha*dt1
          UPD.eewsl(i)=sfcb(i,1,ifcb)
        ELSE IF (fcbt(1:4) .EQ. 'HEWL') THEN
          UPD.hewfcb(i)=xfcb(i,1,ifcb)+alpha*dt1
          UPD.hewsl(i)=sfcb(i,1,ifcb)
        END IF
      END IF
    END DO

  END IF

  RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(dbg_rd_fcb): read the fcb file error '//TRIM(line)
  CALL exit(1)

200 WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(dbg_rd_fcb): end of the fcb file'

END SUBROUTINE
