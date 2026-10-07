!*
SUBROUTINE read_rnxobs_head(lfn,HD,ierr)
!!
!*
USE par
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: lfn,ierr
TYPE(RNXHEAD) :: HD

!*
! The Local variable
!!-------------------------
INTEGER(IT) :: ioerr,i,j,id,nline
CHARACTER(LEN_STRING) :: msg,line,keyword*20
REAL(RL) :: sec

  !*
  ! Start of executable code
  !!------------------------
  ierr = 0

  HD.nsys = 0
  HD.nobstype = 0
  HD.obstype = ''
  DO WHILE(.TRUE.)

    msg='   '
    READ(lfn,'(A)',END=200) line
    keyword = line(61:80)

    !! END of header.
    IF(INDEX(keyword,'END OF HEADER') .NE. 0) RETURN
    IF( (INDEX(keyword,'END OF HEADER').NE.0 .AND. HD.ver.GE.2.d0) &
        .OR. (LEN_TRIM(keyword).EQ.0 .AND. HD.ver.LT.2.d0)) RETURN

    !! RINEX  version
    IF(INDEX(keyword,'RINEX VERSION') .NE. 0) THEN
      READ(line,'(I9,31X,A1)',iostat=ioerr)  i,HD.sys
      IF (ioerr .NE. 0) THEN
        READ(line,'(F9.2,31X,A1)',iostat=ioerr)  HD.ver,HD.sys
      ELSE
        HD.ver=DBLE(i)
      END IF
      IF(HD.sys(1:1) .EQ. ' ') HD.sys='G'
      IF(ioerr .NE. 0) msg = 'READ RINEX VERSION error.'
      !! Raw code, do not support for rnx 4.xx
      ! IF(HD.ver.LT.1.d0 .OR. HD.ver.GT.4.d0) WRITE(msg,'(a,f3.2)') 'invalid RINEX VERSION ',HD.ver
      !! Be competition with rnx 4.xx
      IF(HD.ver.LT.1.d0) WRITE(msg,'(a,f3.2)') 'invalid RINEX VERSION ',HD.ver

    !! site name
    ELSE IF(INDEX(keyword,'MARKER NAME') .NE. 0) THEN
      READ(line,'(a4)') HD.mark

    !! receiver number and type
    ELSE IF(INDEX(keyword,'REC #') .NE. 0) THEN
      HD.recnum=line(1:20)
      HD.rectype=line(21:40)
      HD%recvers=line(41:60)

    !! antenna number and type
    ELSE IF(INDEX(keyword,'ANT #') .NE. 0) THEN
      HD.antnum=line(1:20)
      HD.anttype=line(21:40)

    !! station coordinate
    ELSE IF(INDEX(keyword,'APPROX POSITION').ne.0) THEN
      !READ(line,'(3F14.4)',iostat=ioerr) HD.x,HD.y,HD.z
      READ(line(1:60),*,iostat=ioerr) HD.x,HD.y,HD.z
      IF(ioerr .NE. 0) msg = 'READ APPROX POSITION error'

    !! antenna offset
    ELSE IF(INDEX(keyword,'ANTENNA: DELTA').ne.0) THEN
      !READ(line,'(3f14.4)',iostat=ioerr) HD.h,HD.e,HD.n
      READ(line(1:60),*,iostat=ioerr) HD.h,HD.e,HD.n
      IF(ioerr .NE. 0) msg = 'READ ANTENNA: DELTA error.'

    !! wavelength fact
    ELSE IF(INDEX(keyword,'WAVELENGTH FACT').ne.0) THEN
      READ(line,'(2i6)',iostat=ioerr) HD.fact1,HD.fact2
      IF(ioerr.ne.0) msg = 'READ WAVELENGTH FACT error'

    ELSE IF (INDEX(keyword,'SYS / # / OBS TYPES')) THEN
      !@ CMT BY XSY: W/X/0 > L for LEO
      IF (line(1:1) .EQ. 'X') line(1:1) = 'L'
      IF (line(1:1) .EQ. 'W') line(1:1) = 'L'
      IF (line(1:1) .EQ. '0') line(1:1) = 'L'

      id = INDEX(SYS,line(1:1))
      IF (id .EQ. 0) CYCLE
      ! sinxe 3.03, B1I is set as L2, for previous, we change it to L2
      IF (HD.ver.LE.3.02d0 .AND. id.EQ.INDEX(SYS,'C')) THEN
        j=INDEX(line,'C1')
        IF (j.NE.0) line(j:j+1)='C2'
        j=INDEX(line,'L1')
        IF (j.NE.0) line(j:j+1)='L2'
        j=INDEX(line,'S1')
        IF (j.NE.0) line(j:j+1)='S2'
        j=INDEX(line,'D1')
        IF (j.NE.0) line(j:j+1)='D2'
      END IF

      READ(line,'(3X,I3)',iostat=ioerr) HD.nobstype(id)

      IF (HD.nobstype(id) .GT. MAXOBSTYP) THEN
        WRITE(ERROR_UNIT,*) '***ERROR(read_rnxobs_head): the number of observation is greater than the maximum ',line(1:1),HD.nobstype(id),MAXOBSTYP
        CALL exit(1)
      END IF

      IF (HD.nobstype(id) .EQ. 0) CYCLE
      IF (MOD(HD.nobstype(id),13) .NE. 0) THEN
        nline = INT(HD.nobstype(id)/13)+1
      ELSE
        nline = INT(HD.nobstype(id)/13)
      END IF

      BACKSPACE(lfn)
      DO i=1, nline
        READ(lfn,'(a)',END=200) line
         IF (HD.ver.LE.3.02d0 .AND. id.EQ.INDEX(SYS,'C')) THEN
           j=INDEX(line,'C1')
           IF (j.NE.0) line(j:j+1)='C2'
           j=INDEX(line,'L1')
           IF (j.NE.0) line(j:j+1)='L2'
           j=INDEX(line,'S1')
           IF (j.NE.0) line(j:j+1)='S2'
           j=INDEX(line,'D1')
           IF (j.NE.0) line(j:j+1)='D2'
         END IF
        READ(line,'(6X,13(1x,a3))',iostat=ioerr) &
         (HD.obstype(j,id),j=(i-1)*13+1,MIN(HD.nobstype(id),i*13))
        IF (ioerr .NE. 0) THEN
          msg = 'READ TYPES OF OBSERV error'
          EXIT
        END IF
      END DO
      HD.obstype(HD.nobstype(id)+1:MAXOBSTYP,id)='  '

    ! type of observations
    ELSE IF(INDEX(keyword,'TYPES OF OBSERV').ne.0) THEN
      READ(line,'(i6)',iostat=ioerr) HD.nobstype(1)

      IF (HD.nobstype(1) .GT. MAXOBSTYP) THEN
        WRITE(ERROR_UNIT,*) '***ERROR(read_rnxobs_head): the number of observation is greater than the maximum',HD.nobstype(1),MAXOBSTYP
        CALL exit(1)
      END IF

      IF (HD.nobstype(1) .EQ. 0) CYCLE
      IF (MOD(HD.nobstype(1),9) .NE. 0) THEN
        nline = INT(HD.nobstype(1)/9)+1
      ELSE
        nline = INT(HD.nobstype(1)/9)
      END IF

      BACKSPACE(lfn)
      DO id=1, nline
        READ(lfn,'(a)',END=200) line
        READ(line,'(6X,10(4x,a2))',iostat=ioerr) &
          (HD.obstype(i,1),i=(id-1)*9+1,MIN(HD.nobstype(1),id*9))
        IF(ioerr .NE. 0) THEN
          msg = 'READ TYPES OF OBSERV error'
          EXIT
        END IF
      END DO

      HD.obstype(HD.nobstype(1)+1:MAXOBSTYP,1)='  '

      DO i=2, MAXSYS
        HD.nobstype(i) = HD.nobstype(1)
        HD.obstype(:,i) = HD.obstype(:,1)
      END DO

      !! For WHU's tracking stations, the C1 is C2, and C2 is C7
      i=INDEX(SYS,'C')
      DO j=1, HD.nobstype(i)
        IF (HD.obstype(j,i)(2:2) .EQ. '2') HD.obstype(j,i)(2:2)='7'
        IF (HD.obstype(j,i)(2:2) .EQ. '1') HD.obstype(j,i)(2:2)='2'
      END DO

    !! interval
    ELSE IF(INDEX(keyword,'INTERVAL').ne.0) THEN
      READ(line,*,iostat=ioerr) HD.intv
      IF(ioerr .NE. 0) msg = 'READ INTERVAL error'

    !! start time
    ELSE IF(INDEX(keyword,'TIME OF FIRST OBS').ne.0) THEN
      IF (HD.ver.LT.3.d0) THEN
        READ(line,'(5i6,f13.7)',iostat=ioerr) (HD.t0(i),i=1,5),sec
             HD.t0(6) = nint(sec)
      ELSE
        READ(line,'(5i6,f13.7,5x,a3)',iostat=ioerr) (HD.t0(i),i=1,5),sec,HD.tsys
             HD.t0(6) = nint(sec)
        !!IF((index(HD.tsys,'G').eq.0) .and. (index(HD.tsys,'B').eq.0)) THEN
        !!  WRITE(*,*) 'Unkown Time system, but ',HD.tsys
        !!  CALL exit(1)
        !!ENDIF
      ENDIF
      IF(ioerr.ne.0) msg = ' READ TIME OF FIRST OBS error'

    !! stop  time
    ELSE IF(INDEX(keyword,'TIME OF LAST OBS').ne.0) THEN
      HD.t0(1:5)=0
      sec=0.d0
      IF(HD.ver.LT.3.d0)THEN
        READ(line,'(5i6,f13.7)',iostat=ioerr) (HD%t1(i),i=1,5),sec
        HD.t1(6) = nint(sec)
      ELSE
        READ(line,'(5i6,f13.7,5x,a3)',iostat=ioerr) (HD%t1(i),i=1,5),sec,HD%tsys
             HD.t1(6) = nint(sec)
      END IF
      IF(ioerr.NE.0) msg = ' READ TIME OF LAST OBS error'
    ENDIF
    IF(LEN_TRIM(msg).ne.0) GOTO 100
  ENDDO

  IF (HD.ver .GE. 3.d0) THEN
    DO i=1, MAXSYS
      IF (HD.nobstype(i) .NE. 0) HD.nsys = HD.nsys+1
    END DO
  END IF

!
!! error
100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxobs_head): '//msg
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxobs_head): '//line
! call exit(1)
  ierr = 1
  RETURN

!! END of file
200 CONTINUE
  ierr = 1
  RETURN

END SUBROUTINE
