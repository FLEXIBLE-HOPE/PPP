!*
SUBROUTINE oi_sat_panel(name,PAN)
!!
!*
USE tables
USE satellite
IMPLICIT NONE

!*
! The arguments
!!--------------------
CHARACTER(LEN=*) :: name
TYPE(SATEPAN) :: PAN

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j,lfn

  LOGICAL(LG) :: lfirst,lexist

  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_STRING) :: msg

  DATA lfirst /.TRUE./
  SAVE lfirst,lfn

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit
  CHARACTER(LEN_STRING) :: findkey

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    line=f_tableFileName('paninf')

    lexist=.TRUE.
    INQUIRE(FILE=line,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_sat_panel): '//TRIM(line)//' is not exist.'
      CALL exit(1)
    END IF
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=line)
  END IF

  !! PANELS for GNSS satellites
  msg='+PANELS OF '//TRIM(name)
  line=findkey(lfn,msg,'')
  PAN.npan=0
  DO WHILE(INDEX(line,'-PANELS OF '//TRIM(name)) .EQ. 0)
    READ(lfn,'(A)',END=100,ERR=100) line
    IF (line(1:1) .EQ. ' ') THEN
      PAN.npan=PAN.npan+1
      READ(line,*,ERR=100) i,PAN.area(PAN.npan),PAN.norm(1:3,PAN.npan),PAN.rho(PAN.npan),PAN.delta(PAN.npan),PAN.drag(PAN.npan),PAN.rhoi(PAN.npan),PAN.deltai(PAN.npan)
      PAN.alpha(PAN.npan)=1.d0-PAN.rho(PAN.npan)-PAN.delta(PAN.npan)
      j=INDEX(line,'#')
      IF (j .GT. 0 .AND. j .LT. LEN_TRIM(line)) PAN.name(PAN.npan)=line(j+1:LEN_TRIM(line))
    END IF
  END DO

  RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_sat_panel): read '//TRIM(line)
  CALL exit(1)

END SUBROUTINE
