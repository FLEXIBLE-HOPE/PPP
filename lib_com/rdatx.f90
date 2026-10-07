!*
SUBROUTINE rdatx(nfreq,freq,fjd_beg,fjd_end,AT)
!!
!! purpose   : read antenna atx file
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!    output : ATX -- antenna correction
!! author    : Geng J
!! created   : Oct. 2, 2007
!!
!*
USE atx
USE const
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: nfreq(MAXSYS)
CHARACTER(LEN_FREQ) :: freq(MAXFREQ,MAXSYS)
REAL(RL) :: fjd_beg,fjd_end
TYPE(antatx) :: AT

  !*
  ! The local variables
  !!---------------------------
  LOGICAL(LG) :: lfirst,lfound
  INTEGER(IT) :: i,j,k,lfn,iy,imon,id,ih,im,ncol,nrow,icol,irow,ierr,isys
  REAL(RL) :: sod,fjd0,fjd1
  CHARACTER :: atxtyp
  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_FILENAME) :: atxfile

  DATA lfirst/.TRUE./
  SAVE lfirst,lfn,atxtyp

  TYPE(T_FILETABLE) :: FT

  !*
  ! The function called
  !!------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst) THEN
    lfirst=.FALSE.
    atxfile = f_tableFileName("antpcv")
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=atxfile,STATUS='old',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A,I4)') '***ERROR(rdatx): open antenna file abs_igs.atx',ierr
      CALL exit(1)
    END IF

    !! read header
    line=' '
    DO WHILE(INDEX(line,'END OF HEADER') .NE. 61)
      READ(lfn,'(A)',END=100) line
      IF (INDEX(line,'PCV TYPE / REFANT') .EQ. 61) THEN
        READ(line,'(A1)') atxtyp
      !! the SINEX CODE
      ELSE IF (line(1:8) .EQ. 'Changes:') THEN
        READ(lfn,'(A)',END=100) line
        IF (line(3:6).EQ.'week' .AND. line(61:67).EQ.'COMMENT') THEN
          ATXWEEK(7:10)=line(8:11)
        END IF
      END IF
    END DO
  END IF

  IF (atxtyp .NE. 'A') THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rdatx): absolute antenna mode in abs_igs.atx'
    CALL exit(1)
  END IF

  !! find the right antenna
10 CONTINUE
  REWIND(lfn)

  lfound=.FALSE.
  DO WHILE(.TRUE.)
    line=' '
    DO WHILE(INDEX(line,'TYPE / SERIAL NO') .NE. 61)
      READ(lfn,'(A)',END=100) line
    END DO

    IF (AT.antnam .NE. line(1:20)) CYCLE

    j = LEN_TRIM(AT.antnum)
    IF (j.NE.0) THEN
      IF (AT.antnum(1:j) .NE. line(21:21+j-1)) CYCLE
    END IF

    AT.dazi =0.d0
    AT.dzen =0.d0
    AT.neu  =0.d0
    AT.nfreq = 0
    AT.pcv=0.d0
    fjd0=0.d0
    fjd1=1.d10
    DO WHILE(INDEX(line,'START OF FREQUENCY') .NE. 61)
      READ(lfn,'(A)',END=100) line
      IF (INDEX(line,'DAZI') .EQ. 61) THEN
        READ(line,*,ERR=200) AT.dazi
      ELSE IF (INDEX(line,'ZEN1 / ZEN2 / DZEN') .EQ. 61) THEN
        READ(line,*,ERR=200) AT.zen1,AT.zen2,AT.dzen
      ELSE IF (INDEX(line,'VALID FROM') .EQ. 61) THEN
        READ(line,*,ERR=200) iy,imon,id,ih,im,sod
        fjd0=modified_julday(id,imon,iy)+ih/24.d0+im/1440.d0+sod/86400.d0
      ELSE IF (INDEX(line,'VALID UNTIL') .EQ. 61) THEN
        READ(line,*,ERR=200) iy,imon,id,ih,im,sod
        fjd1=modified_julday(id,imon,iy)+ih/24.d0+im/1440.d0+sod/86400.d0
      ELSE IF (INDEX(line,'# OF FREQUENCIES') .EQ. 61) THEN
        READ(line,*,ERR=200) AT.nfreq
      END IF
    END DO
    IF ((fjd0-fjd_beg)*86400.d0.GT.MAXWND .OR. (fjd_end-fjd1)*86400.d0.GT.MAXWND) CYCLE
    lfound=.TRUE.
    EXIT
  END DO

100 CONTINUE
  BACKSPACE(lfn)
  IF (.NOT.lfound) THEN
    WRITE(ERROR_UNIT,'(3A)') '***ERROR(rdatx): Antenna not found ',AT.antnam,AT.antnum
    CALL exit(1)
  END IF

  !! The start and end epoch for SINEX output
  AT.mjd1=fjd0
  IF (fjd1 .NE. 1.d10) AT.mjd2=fjd1

  !! read offset and phase center variation
  ncol=0
  nrow=0
  IF (AT.dazi .NE. 0.d0) nrow=INT(360.d0/AT.dazi)+1
  IF (AT.dzen .NE. 0.d0) ncol=INT((AT.zen2-AT.zen1)/AT.dzen)+1

  DO i=1, AT.nfreq

    !* find frequency
101 line=' '
    DO WHILE(line(61:78) .NE. 'START OF FREQUENCY')
      READ(lfn,'(A)') line
      IF (INDEX(line,'END OF ANTENNA') .NE. 0) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(rdatx): frquency not found '
        CALL exit(1)
      END IF
    END DO

    isys=0
    k=0
    SELECT CASE(line(4:6))
      CASE('G01')
        isys=INDEX(SYS,'G')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('G02')
        isys=INDEX(SYS,'G')
        k = pointer_string(nfreq(isys),freq(:,isys),'L2')
      CASE('G05')
        isys=INDEX(SYS,'G')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
      CASE('R01')
        isys=INDEX(SYS,'R')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('R02')
        isys=INDEX(SYS,'R')
        k = pointer_string(nfreq(isys),freq(:,isys),'L2')
      CASE('E01')
        isys=INDEX(SYS,'E')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('E5A','E5a','E05')
        isys=INDEX(SYS,'E')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
      CASE('E5B','E5b','E07')
        isys=INDEX(SYS,'E')
        k = pointer_string(nfreq(isys),freq(:,isys),'L7')
      CASE('EAB','Eab','E5X','E5x','E08')
        isys=INDEX(SYS,'E')
        k = pointer_string(nfreq(isys),freq(:,isys),'L8')
      CASE('E06')
        isys=INDEX(SYS,'E')
        k = pointer_string(nfreq(isys),freq(:,isys),'L6')
      CASE('B1c','B1C','B01','C01')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('B1I','B02','C02')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L2')
      CASE('B2A','B2a','B05','C05')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
      CASE('B3I','B03','B06','C03','C06')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L6')
      CASE('B2I','B2B','B2b','B07','C07')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L7')
      CASE('B2X','B2x','B08','C08')
        isys=INDEX(SYS,'C')
        k = pointer_string(nfreq(isys),freq(:,isys),'L8')
      CASE('J01')
        isys=INDEX(SYS,'J')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('J02')
        isys=INDEX(SYS,'J')
        k = pointer_string(nfreq(isys),freq(:,isys),'L2')
      CASE('J05')
        isys=INDEX(SYS,'J')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
      CASE('J06')
        isys=INDEX(SYS,'J')
        k = pointer_string(nfreq(isys),freq(:,isys),'L6')
      CASE('I05')
        isys=INDEX(SYS,'I')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
      CASE('I09')
        isys=INDEX(SYS,'I')
        k = pointer_string(nfreq(isys),freq(:,isys),'L9')
      CASE('L01')
        isys=INDEX(SYS,'L')
        k = pointer_string(nfreq(isys),freq(:,isys),'L1')
      CASE('L05')
        isys=INDEX(SYS,'L')
        k = pointer_string(nfreq(isys),freq(:,isys),'L5')
    END SELECT

    IF (isys.EQ.0 .OR. k.EQ.0) CYCLE

    AT.freq(k,isys)=line(4:6) !this is null for mutl-frequency which no pco

    !! reading
    irow=0
    DO WHILE(.TRUE.)
      READ(lfn,'(A)') line
      IF (line(61:76) .EQ. 'END OF FREQUENCY') EXIT
      IF (line(61:77) .EQ. 'NORTH / EAST / UP') THEN
        READ(line,*,ERR=200) (AT.neu(j,k,isys),j=1,3)
        !! convert unit from mm to m
        DO j=1,3
          AT.neu(j,k,isys)=AT.neu(j,k,isys)*1.d-3
        END DO
      ELSE IF (line(4:8) .EQ. 'NOAZI') THEN
        READ(line(9:),*,ERR=200) (AT.pcv(icol,0,k,isys),icol=1,ncol)
      ELSE
        irow=irow+1
        !IF (irow .GT. 80) WRITE(*,*)irow
        READ(line(9:),*,ERR=200) (AT.pcv(icol,irow,k,isys),icol=1,ncol)
      END IF
    END DO
    IF (irow .NE. nrow) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(rdatx): pcv lost '
      CALL exit(1)
    END IF

    !! convert unit from mm to m
    DO irow=0,nrow
      DO icol=1,ncol
        AT.pcv(icol,irow,k,isys)=AT.pcv(icol,irow,k,isys)*1.d-3
      END DO
    END DO
  END DO

  ! for mult-frequency, if there is no PCC information, using the first, but the error 
  ! will be larger than using the second, so PCC for other frequency is needed
  DO isys=1, MAXSYS
    IF (nfreq(isys) .EQ. 0) CYCLE
    DO k=1, nfreq(isys)
      IF (COUNT(AT.neu(1:3,k,isys).EQ.0.d0) .EQ. 3) THEN
        icol = pointer_string(nfreq(isys),freq(:,isys),'L1')
        !! xsy
        IF (icol .EQ. 0) THEN
          ! icol = pointer_string(nfreq(isys),freq(:,isys),'L2')
          ! IF (icol .EQ. 0) THEN
          !   icol = pointer_string(nfreq(isys),freq(:,isys),'L5')
          ! ENDIF
          ! IF (icol .EQ. 0) CYCLE
          CYCLE
        ENDIF
        AT.neu(1:3,k,isys)=AT.neu(1:3,icol,isys)
        AT.pcv(1:50,0:80,k,isys)=AT.pcv(1:50,0:80,icol,isys)
      END IF
    END DO
  END DO

  RETURN

200 CONTINUE
  WRITE(ERROR_UNIT,'(2A)') '***ERROR(rdatx): read file ', trim(line)
  CALL exit(1)

END SUBROUTINE

!*
SUBROUTINE atxwk()
!!
!*
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: lfn,ierr
  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_FILENAME) :: atxfile

  TYPE(T_FILETABLE) :: FT

  !*
  ! The function called
  !!------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  atxfile = f_tableFileName("antpcv")
  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=atxfile,STATUS='OLD',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A,I4)') '***ERROR(atxwk): open antenna file abs_igs.atx',ierr
    CALL exit(1)
  END IF

  !! read header
  line=' '
  DO WHILE(INDEX(line,'END OF HEADER') .NE. 61)
    READ(lfn,'(A)',END=100) line
    IF (line(1:8) .EQ. 'Changes:') THEN
      READ(lfn,'(A)',END=100) line
      IF (line(3:6).EQ.'week' .AND. line(61:67).EQ.'COMMENT') THEN
        ATXWEEK(7:10)=line(8:11)
        EXIT
      END IF
    END IF
  END DO

100 CONTINUE

  CLOSE(lfn)

  RETURN

END SUBROUTINE
