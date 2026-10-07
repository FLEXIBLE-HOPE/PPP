!*
SUBROUTINE read_rnxnav(csys,flnbrd,mjd0,mjd1,hd,neph,ephm,ephg)
!!
!*
USE brdeph
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
CHARACTER :: csys
CHARACTER(LEN_FILENAME) :: flnbrd
REAL(RL) :: mjd0,mjd1
INTEGER(IT) :: neph(1:*)
TYPE(GPS_BRDEPH) :: ephm(MAXEPH,1:*)
TYPE(GLONASS_BRDEPH) :: ephg(MAXEPH)
TYPE(BRDHEAD) :: HD

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: lfn
  INTEGER(IT) :: i,k,isys
  INTEGER(IT) :: iy,im,id,ih,imin,isec

  TYPE(GPS_BRDEPH) :: eph
  TYPE(GLONASS_BRDEPH) :: ephr
  CHARACTER(LEN_STRING) :: line,fmt

  INTEGER(IT) :: mjd
  REAL(RL) :: dt0,dt1,sod

  LOGICAL(LG) :: lalready

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!----------------------------

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flnbrd)

  line=''
  DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
    READ(lfn,'(A)',END=100,ERR=200) line

    SELECT CASE(line(61:LEN_TRIM(line)))
      ! RINEX 2.00, 2.10, 2.11, 3.00, 3.01, 3.02, 3.03
      CASE('RINEX VERSION / TYPE')
        READ(line,'(F9.2)') HD.ver
        IF (HD.ver .GE. 3.d0) THEN
          IF (csys(1:1).EQ.'M' .AND. line(41:41).NE.'M') THEN
            WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(read_rnxnav): there are only broadcast messages for '//line(41:41)
          END IF
        END IF
      ! RINEX 2.00, 2.10, 2.11, 3.00, 3.01, 3.02, 3.03
      CASE('PGM / RUN BY / DATE')
      CASE('COMMENT')
      ! RINEX 2.00, 2.10, 2.11 for GPS
      CASE('ION ALPHA')
        i=INDEX(SYS,csys)
        IF (i .EQ. 0) THEN
          WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_rnxnav): ION ALPHA is only valid for single system'
          CYCLE
        END IF
        READ(line,'(2X,4D12.4)') HD.ion(1:4,1,i)
      CASE('ION BETA')
        i=INDEX(SYS,csys)
        IF (i .EQ. 0) THEN
          WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_rnxnav): ION BETA is only valid for single system'
          CYCLE
        END IF
        READ(line,'(2X,4D12.4)') HD.ion(1:4,2,i)
      ! RINEX 3.00, 3.01, 3.02, 3.03
      CASE('IONOSPHERIC CORR')
        SELECT CASE(line(1:4))
          CASE('GPS ','GPSA')
            i=INDEX(SYS,'G')
            k=1
          CASE('GPSB')
            i=INDEX(SYS,'G')
            k=2
          CASE('GAL')
            i=INDEX(SYS,'E')
            k=1
          CASE('BDS ','BDSA')
            i=INDEX(SYS,'C')
            k=1
          CASE('BDSB')
            i=INDEX(SYS,'C')
            k=2
          ! 3.03
          CASE('QZS ','QZSA')
            i=INDEX(SYS,'J')
            k=1
          CASE('QZSB')
            i=INDEX(SYS,'J')
            k=2
          ! 3.03
          CASE('IRN ','IRNA')
            i=INDEX(SYS,'I')
            k=1
          CASE('IRNB')
            i=INDEX(SYS,'I')
            k=2
          CASE DEFAULT
            WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_rnxnav): unknown IONOSPHERIC CORR '//TRIM(line(1:4))
            CYCLE
        END SELECT
        HD.ionc(k,i)=line(1:4)
        READ(line,'(5X,4D12.4,1X,A1,1X,A1)') HD.ion(1:4,k,i),HD.tmark(k,i),HD.svid(k,i)
      ! 2.00, 2.11, 2.12
      CASE('DELTA-UTC: A0,A1,T,W')
        i=INDEX(SYS,csys)
        IF (i .EQ. 0) THEN
          WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_rnxnav): DELTA-UTC: A0,A1,T,W is only valid for single system'
          CYCLE
        END IF
        READ(line,'(3X,2D19.12,2F9.0)') HD.tim(1:4,1,i)
      CASE('TIME SYSTEM CORR')
        SELECT CASE(line(1:4))
          CASE('GPUT')
            k=1
            i=INDEX(SYS,'G')
          CASE('GLUT')
            k=1
            i=INDEX(SYS,'R')
          CASE('GAUT')
            k=1
            i=INDEX(SYS,'E')
          CASE('BDUT')
            k=1
            i=INDEX(SYS,'C')
          CASE('QZUT')
            k=1
            i=INDEX(SYS,'J')
          CASE('IRUT')
            k=1
            i=INDEX(SYS,'I')
          CASE('SBUT')
            k=1
            i=INDEX(SYS,'S')
          CASE('GPGA','GAGP')!xsy add :GAGP
            k=2
            i=INDEX(SYS,'G')
          CASE('GLGP')
            k=2
            i=INDEX(SYS,'R')
          CASE('QZGP')
            k=2
            i=INDEX(SYS,'J')
          CASE('IRGP')
            k=2
            i=INDEX(SYS,'I')
          CASE DEFAULT
            WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_rnxnav): unknown TIME SYSTEM CORR '//TRIM(line(1:4))
            CYCLE
        END SELECT
        HD.timc(k,i)=line(1:4)
        READ(line,'(5X,D17.10,D16.9,F7.0,F5.0,1X,A5,1X,I2)') HD.tim(1:4,k,i),HD.satb(k,i),HD.utcid(k,i)
      CASE('LEAP SECONDS')
        IF (HD.ver .LT. 3.0) THEN
          READ(line,'(I6)') HD.leap(1)
        ! RINEX 3.00, 3.01, 3.02, 3.03
        ELSE
          READ(line,'(4I6,A3)') HD.leap,HD.lpts
        END IF
    END SELECT
  END DO

  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=100,ERR=200) line
    BACKSPACE(lfn)

    IF (HD.ver .LT. 3.d0) THEN
      line=csys//TRIM(line)
    END IF

    SELECT CASE(line(1:1))
      CASE('G','E','C','J','I')

        isys=INDEX(SYS,line(1:1))

        IF (HD.ver .LT. 3.d0) THEN
          fmt='(I2,5I3,F5.1,3D19.12/6(3X,4D19.12/),3X,4D19.12)'
        ELSE
          fmt='(A3,1X,I4,5(1X,I2.2),3D19.12/6(4X,4D19.12/),4X,4D19.12)'
        END IF
        !IF (line(1:1).EQ.'J' .AND. csys(1:1).EQ.'J') fmt='(1X,I2.2,5I3,F5.1,3D19.12/6(4X,4D19.12/),4X,4D19.12)'

        IF (HD.ver .LT. 3.d0) THEN
          READ(lfn,fmt,END=100,ERR=300) i,iy,im,id,ih,imin,isec,eph.a0,eph.a1,eph.a2     &
           ,eph.aode,eph.crs,eph.dn,eph.m0,eph.cuc,eph.e,eph.cus,eph.roota,eph.toe,eph.cic   &
           ,eph.omega0,eph.cis,eph.i0,eph.crc,eph.omega,eph.omegadot,eph.idot,eph.resvd0     &
           ,eph.week,eph.resvd1,eph.accu,eph.hlth,eph.tgd,eph.aodc,eph.tom,eph.fih
          WRITE(eph.cprn,'(A1,I2.2)') csys,i

        ELSE
          READ(lfn,fmt,END=100,ERR=300) eph.cprn,iy,im,id,ih,imin,isec,eph.a0,eph.a1,eph.a2     &
           ,eph.aode,eph.crs,eph.dn,eph.m0,eph.cuc,eph.e,eph.cus,eph.roota,eph.toe,eph.cic   &
           ,eph.omega0,eph.cis,eph.i0,eph.crc,eph.omega,eph.omegadot,eph.idot,eph.resvd0     &
           ,eph.week,eph.resvd1,eph.accu,eph.hlth,eph.tgd,eph.aodc,eph.tom,eph.fih

        END IF

        ! According to RINEX 3.03, the broadcast navigation message in their own time system
        ! However, except BDS and GLONASS, Galileo, QZSS, and IRNSS are aligned their time
        ! to GPS time system with nanosecond bias

        ! Do not change the time system

        !! Time of clock
        CALL yr2year(iy)
        eph.mjd=modified_julday(id,im,iy)
        eph.sod=ih*3600.d0+imin*60.d0+DBLE(isec)
        ! BDS/GLONASS Time to GPS Time
        !IF (eph.cprn(1:1) .EQ. 'R') CALL brdtime(eph.cprn,eph.mjd,eph.sod)

        !! the week in broadcast file generated by WHU is GPS week,
        !! but that in IGS meraged file is BDS week
        !! the broadcast file download from iGMAS is also in BDS week
        CALL mjd2wksow(eph.mjd,eph.sod,k,dt1)
        ! roundover
        !IF (eph.week .EQ. 0.d0) then
        !  write(*,*) eph.week 
        !  eph.week = k
        !end if
        ! eph.week in BDS week, so convered to GPS week
        IF (DBLE(DBLE(k)-eph.week) .GT. 1350.d0) eph.week=1356+eph.week

        !! Time of ephemeris
        !mjd=INT((eph.week*7+44244)+eph.toe/86400.d0)
        !sod=(eph.week*7+44244-eph.mjd)*86400.d0+eph.toe
        !! BDS/GLONASS Time to GPS Time
        !IF (eph.cprn(1:1) .EQ. 'R') CALL brdtime(eph.cprn,mjd,sod)
        !CALL mjd2wksow(mjd,sod,k,eph.toe)
        !eph.week=DBLE(k)

        !! check time
        dt0=0.d0
        dt1=0.d0
        IF (mjd0 .NE. 0.d0) THEN
          dt0=eph.mjd+eph.sod/86400.d0-mjd0
        END IF
        IF (mjd1 .NE. 0.d0) THEN
          dt1=eph.mjd+eph.sod/86400.d0-mjd1
        END IF
        IF (dt0.LT.-1.d0/24.d0 .OR. dt1.GT.1.d0/24.d0) CYCLE

        lalready=.FALSE.
        DO i=1, neph(isys)
          IF (ephm(i,isys).cprn.EQ.eph.cprn .AND. eph.mjd.EQ.ephm(i,isys).mjd .AND. ephm(i,isys).sod.EQ.eph.sod) THEN
            lalready = .TRUE.
            EXIT
          END IF
        END DO

        IF (lalready .EQ. .FALSE.) THEN
          neph(isys)=neph(isys)+1
          IF (neph(isys) .GT. MAXEPH) THEN
            WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxnav): exceed the maimum of ephemeris records.'
            CALL exit(1)
          END IF
          ephm(neph(isys),isys)=eph
        END IF

      CASE('R')

        isys=INDEX(SYS,line(1:1))

        IF (HD.ver .LT. 3.d0) THEN
          fmt='(I2,5I3,F5.1,3D19.12/2(3X,4D19.12/),3X,4D19.12)'
        ELSE
          fmt='(A3,1X,I4,5(1X,I2.2),3D19.12,/2(4X,4D19.12/),4X,4D19.12)'
        END IF

        IF (HD.ver .LT. 3.d0) THEN
          READ(lfn,fmt,END=100,ERR=300) i,iy,im,id,ih,imin,isec,ephr.tau, ephr.gamma, ephr.tk, &
                       ephr.pos(1), ephr.vel(1), ephr.acc(1), ephr.health, &
                       ephr.pos(2), ephr.vel(2), ephr.acc(2), ephr.frenum, &
                       ephr.pos(3), ephr.vel(3), ephr.acc(3), ephr.age
          WRITE(ephr.cprn,'(A1,I2.2)') csys,i
        ELSE
          READ(lfn,fmt,END=100,ERR=300) ephr.cprn,iy,im,id,ih,imin,isec,ephr.tau, ephr.gamma, ephr.tk, &
                       ephr.pos(1), ephr.vel(1), ephr.acc(1), ephr.health, &
                       ephr.pos(2), ephr.vel(2), ephr.acc(2), ephr.frenum, &
                       ephr.pos(3), ephr.vel(3), ephr.acc(3), ephr.age
        END IF

        ! Only the health GLONASS satellites
        IF (ephr.health .NE. 0.d0) CYCLE

        CALL yr2year(iy)
        ephr.mjd=modified_julday(id,im,iy)
        ephr.sod=ih*3600.d0+imin*60.d0+DBLE(isec)

        ! GLONASS Time to GPS Time
        CALL brdtime(ephr.cprn,ephr.mjd,ephr.sod)

        !* check time
        dt0=0.d0
        dt1=0.d0
        IF (mjd0 .NE. 0.d0) THEN
          dt0=ephr.mjd+ephr.sod/86400.d0-mjd0
        END IF
        IF (mjd1 .NE. 0.d0) THEN
          dt1=ephr.mjd+ephr.sod/86400.d0-mjd1
        END IF
        IF (dt0.LT.-1.d0/24.d0 .OR. dt1.GT.1.d0/24.d0) CYCLE

        lalready=.FALSE.
        DO i=1, neph(isys)
          IF (ephg(i).cprn.EQ.ephr.cprn .AND. ephg(i).mjd.EQ.ephr.mjd .AND. ephg(i).sod.EQ.ephr.sod) THEN
            lalready = .TRUE.
            EXIT
          END IF
        END DO

        IF (lalready .EQ. .FALSE.) THEN
          neph(isys)=neph(isys)+1
          IF (neph(isys) .GT. MAXEPH) THEN
            WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxnav): exceed the maimum of ephemeris records.'
            CALL exit(1)
          END IF
          ephg(neph(isys))=ephr
        END IF

      CASE('S')

        IF (HD.ver .LT. 3.d0) THEN
          fmt='(I2,5I3,F5.1,3D19.12/2(3X,4D19.12/),3X,4D19.12)'
        ELSE
          fmt='(A3,1X,I4,5(1X,I2.2),3D19.12,/2(4X,4D19.12/),4X,4D19.12)'
        END IF

        IF (HD.ver .LT. 3.d0) THEN
          READ(lfn,fmt,END=100,ERR=300) i,iy,im,id,ih,imin,isec,ephr.tau, ephr.gamma, ephr.tk, &
                       ephr.pos(1), ephr.vel(1), ephr.acc(1), ephr.health, &
                       ephr.pos(2), ephr.vel(2), ephr.acc(2), ephr.frenum, &
                       ephr.pos(3), ephr.vel(3), ephr.acc(3), ephr.age
          WRITE(ephr.cprn,'(A1,I2.2)') csys,i
        ELSE
          READ(lfn,fmt,END=100,ERR=300) ephr.cprn,iy,im,id,ih,imin,isec,ephr.tau, ephr.gamma, ephr.tk, &
                       ephr.pos(1), ephr.vel(1), ephr.acc(1), ephr.health, &
                       ephr.pos(2), ephr.vel(2), ephr.acc(2), ephr.frenum, &
                       ephr.pos(3), ephr.vel(3), ephr.acc(3), ephr.age
        END IF

      CASE DEFAULT
        WRITE(OUTPUT_UNIT,'(A)') '%%%WARNING(read_rnxnav): '//TRIM(line)
        DO WHILE(.TRUE.)
          READ(lfn,'(A)',END=100,ERR=200) line
          SELECT CASE(line(1:1))
            CASE('G','R','E','C','J','I')
              BACKSPACE(lfn)
              EXIT
            CASE DEFAULT
              CONTINUE
          END SELECT
        END DO
      END SELECT

300 CONTINUE

  END DO

  CLOSE(lfn)

  RETURN

100 CONTINUE
    CLOSE(lfn)
    RETURN

200 WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rnxnav): read '//TRIM(line)
    CALL exit(1)

END SUBROUTINE
