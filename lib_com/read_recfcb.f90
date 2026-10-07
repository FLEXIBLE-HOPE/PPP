!*
SUBROUTINE read_recfcb(CKF,UPD)
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

  !*
  ! The local variables
  !!-----------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: lfnrecbias,nrectype,isat,nsat,ierr,imode,i
  CHARACTER(LEN_PRN) :: cprn
  CHARACTER(LEN_STRING) :: line,rectype
  LOGICAL(LG) :: lexist

  DATA lfirst /.TRUE./
  SAVE lfirst,lfnrecbias

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! open file
  IF (lfirst) THEN
    lfirst =.FALSE.
    INQUIRE(FILE=CKF%flnrbias,EXIST=lexist)
    IF (lexist .EQ. .TRUE.) THEN
        lfnrecbias = get_valid_unit(10)
        OPEN(UNIT=lfnrecbias,FILE=CKF%flnrbias,STATUS='OLD')
        line = ''
        nrectype = 0
        DO WHILE(INDEX(line,'END OF COMMMENT') .EQ. 0)
            READ(lfnrecbias,'(a)',END=100) line
            IF (INDEX(line,'Receiver Type') .NE. 0) THEN
                nrectype = nrectype + 1
                READ(line(32:),'(A40)',IOSTAT=ierr) rectype
                UPD%rectype(nrectype) = TRIM(rectype)
            END IF
        END DO
        ! DO i=1,nrectype
        !     WRITE(*,*)i,UPD%rectype(i)
        ! END DO
        DO WHILE(INDEX(line,'END OF BIAS') .EQ. 0)
            READ(lfnrecbias,'(a)',END=100) line
            IF (line(1:4) .EQ. '  NL') THEN
                imode = 1
            ELSE IF (line(1:4) .EQ. '  WL') THEN
                imode = 2
            ELSE IF (line(1:4) .EQ. ' EWL') THEN
                imode = 3
            ELSE IF (line(1:4) .EQ. 'EEWL') THEN
                imode = 4
            ELSE IF (line(1:4) .EQ. 'HEWL') THEN
                imode = 5
            ELSE
                CYCLE
            END IF

            cprn = line(6:8)
            isat = pointer_string(CKF.nprn,CKF.cprn,cprn)
            IF (isat .EQ. 0) CYCLE
            DO i=1,nrectype
                READ(line( 9+(i-1)*51 :  9+(i-1)*51+39),'(F40.3)') UPD%recfcb(isat,imode,i)
                READ(line(49+(i-1)*51 : 49+(i-1)*51+9 ),'(F10.3)') UPD%recfcbstd(isat,imode,i)
                ! WRITE(*,'(A3,2X,I2,2X,F5.2,2X,F5.2)')cprn,i,UPD%recfcb(isat,imode,i),UPD%recfcbstd(isat,imode,i)
                
                ! xsy: the therold is 0.05
                IF (UPD%recfcbstd(isat,imode,i) .GE. 0.15) THEN
                    UPD%recfcbstd(isat,imode,i) = 0.D0
                END IF
            END DO
        END DO
    ELSE
    !   WRITE(ERROR_UNIT,'(A)') '###WARNING(read_recfcb): the recfcb file is not exist '//TRIM(CKF%flnrbias)
      lfnrecbias=0
      RETURN        
    END IF
  END IF

  RETURN

100 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_recfcb): read the fcb file error '//TRIM(line)
  CALL exit(1)

200 WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(read_recfcb): end of the fcb file'

END SUBROUTINE
