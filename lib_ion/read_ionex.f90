!*
SUBROUTINE read_ionex(flnion,IOD)
!!
!!
!*
USE const
USE ion
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
CHARACTER(LEN=*) :: flnion
TYPE(IONEX) :: IOD


  !*
  ! The variables
  !!-----------------------
  INTEGER(IT) :: i,l,k,m,ilat,ihgt
  INTEGER(IT) :: nline,lfn
  
  INTEGER(IT) :: iy,im,id,ih,imin,isec
  
  REAL(RL) :: lat,lon1,lon2,dlon,hgt
  
  CHARACTER(LEN_STRING) :: line
  
  LOGICAL(LG) :: lexist

  !*
  ! The function called
  !!-----------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: get_valid_unit
  
  !*
  ! Start the exectuable code
  !!----------------------------
  
  INQUIRE(FILE=flnion,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(OUTPUT_UNIT,'(A)') '###WARNING(read_ionex): '//TRIM(flnion)//' is not exist!'
    RETURN
  END IF

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=flnion)
  
  line=''
  DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
    READ(lfn,'(A)',ERR=100,END=100) line
    SELECT CASE(line(61:LEN_TRIM(line)))
      CASE('IONEX VERSION / TYPE')
        READ(line,'(F8.1,12X,A1,19X,A3,17X)') IOD.hd.ver,IOD.hd.typ,IOD.hd.sys
      CASE('PGM / RUN BY / DATE')
        READ(line,'(3A20)') IOD.hd.pgm,IOD.hd.run,IOD.hd.dat
      CASE('EPOCH OF FIRST MAP')
        ! READ(line,'(6I6)') iy,im,id,ih,imin,isec
        READ(line,'(5I6,F6.0)') iy,im,id,ih,imin,isec
        IOD.hd.mjd0=modified_julday(id,im,iy)
        IOD.hd.sod0=ih*3600.d0+imin*60.d0+isec
      CASE('EPOCH OF LAST MAP')
        ! READ(line,'(6I6)') iy,im,id,ih,imin,isec
        READ(line,'(5I6,F6.0)') iy,im,id,ih,imin,isec
        IOD.hd.mjd1=modified_julday(id,im,iy)
        IOD.hd.sod1=ih*3600.d0+imin*60.d0+isec
      CASE('INTERVAL')
        READ(line,'(I6)') IOD.hd.intv
      CASE('# OF MAPS IN FILE')
        READ(line,'(I6)') IOD.hd.nmaps
      CASE('MAPPING FUNCTION')
        READ(line,'(2X,A4)') IOD.hd.mapfunc
      CASE('ELEVATION CUTOFF')
        READ(line,'(F8.1)') IOD.hd.elev
      CASE('# OF STATIONS')
        READ(line,'(I6)') IOD.hd.nsit
      CASE('# OF SATELLITES')
        READ(line,'(I6)') IOD.hd.nsat
      CASE('RADIUS')
        READ(line,'(F8.1)') IOD.hd.radius
      CASE('MAP DIMENSION')
        READ(line,'(I6)') IOD.hd.ndim
      CASE('HGT1 / HGT2 / DHGT')
        READ(line,'(2X,3F6.1)') IOD.hd.hgt1,IOD.hd.hgt2,IOD.hd.dhgt
      CASE('LAT1 / LAT2 / DLAT')
        READ(line,'(2X,3F6.1)') IOD.hd.lat1,IOD.hd.lat2,IOD.hd.dlat
      CASE('LON1 / LON2 / DLON')
        READ(line,'(2X,3F6.1)') IOD.hd.lon1,IOD.hd.lon2,IOD.hd.dlon
      CASE('EXPONENT')
        READ(line,'(I6)') IOD.hd.unt
      CASE('PRN / BIAS / RMS')
        i=0
        DO WHILE(INDEX(line,'PRN / BIAS / RMS') .NE. 0)
          i=i+1
          READ(line,'(3X,A3,2F10.3)') IOD.hd.cprn(i),IOD.hd.dcb(i),IOD.hd.rms(i)
          IF (IOD.hd.cprn(i)(1:1) .EQ. ' ') IOD.hd.cprn(i)(1:1)='G'
          READ(lfn,'(A)',ERR=100,END=100) line
        END DO
        IOD.hd.nprn=i
        BACKSPACE(lfn)
      CASE DEFAULT
        CONTINUE
    END SELECT
  END DO
  
  ! Check the consistent
  IF (INT(timdif(IOD.hd.mjd1,IOD.hd.sod1,IOD.hd.mjd0,IOD.hd.sod0)/IOD.hd.intv)+1 .NE. IOD.hd.nmaps) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): The number of maps is not consisted with internal and start/end epochs'
    CALL exit(1)  
  END IF

  IF (IOD.hd.ndim.EQ.2 .AND. IOD.hd.dhgt.NE.0.d0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): dght do not equal 0 for 2D ION maps.'
    CALL exit(1)  
  END IF

  IF (IOD.hd.nmaps .GT. MAXIONMAPS) THEN
    WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): The number of maps is greater then the maximum ', MAXIONMAPS
    CALL exit(1)
  END IF
  
  IOD.nlat=INT(DABS((IOD.hd.lat2-IOD.hd.lat1)/IOD.hd.dlat))+1
  IOD.nlon=INT(DABS((IOD.hd.lon2-IOD.hd.lon1)/IOD.hd.dlon))+1
  IOD.nhgt=INT(DABS((IOD.hd.hgt2-IOD.hd.hgt1)/IOD.hd.dhgt))+1
  
  IF (IOD.nlat .GT. MAXIONLAT) THEN
    WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): The number of latitude is greater then the maximum ', MAXIONLAT
    CALL exit(1)
  END IF

  IF (IOD.nlon .GT. MAXIONLON) THEN
    WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): The number of longtitude is greater then the maximum ', MAXIONLON
    CALL exit(1)
  END IF

  IF (IOD.nhgt .GT. MAXIONHGT) THEN
    WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): The number of height is greater then the maximum ', MAXIONHGT
    CALL exit(1)
  END IF  
  
  
  DO WHILE(INDEX(line,'START OF TEC MAP') .EQ. 0)
    READ(lfn,'(A)',ERR=100,END=100) line
  END DO
  BACKSPACE(lfn)
  
  ! Read the TECU
  DO i=1, IOD.hd.nmaps
    ! Start of TEC MAP
    READ(lfn,'(A)',ERR=100,END=100) line
    READ(line,'(I6)') IOD.iondat(i).ind
    IF (IOD.iondat(i).ind.EQ.0 .OR. IOD.iondat(i).ind.NE.i) THEN
      WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): START OF TEC MAP is 0 or do not equal ',i
      CALL exit(1)
    END IF

    ! Epoch of current MAP
    READ(lfn,'(A)',ERR=100,END=100) line
    ! READ(line,'(6I6)') iy,im,id,ih,imin,isec
    READ(line,'(5I6,F6.0)') iy,im,id,ih,imin,isec
    IOD.iondat(i).mjd=modified_julday(id,im,iy)
    IOD.iondat(i).sod=ih*3600.d0+imin*60.d0+isec
    IF (DABS(timdif(IOD.hd.mjd0,IOD.hd.sod0+(i-1)*IOD.hd.intv,IOD.iondat(i).mjd,IOD.iondat(i).sod)) .GT. MAXWND) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): epoch is not correct.'
      CALL exit(1)
    END IF

    DO l=1, IOD.nlat, 1
      ! LAT/LON/LON2/DLON/H
      READ(lfn,'(A)',ERR=100,END=100) line
      READ(line,'(2X,5F6.1)') lat,lon1,lon2,dlon,hgt
      IF (lon1.NE.IOD.hd.lon1 .OR. lon2.NE.IOD.hd.lon2 .OR. dlon.NE.IOD.hd.dlon) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): longtitude is not consistent.'
        WRITE(ERROR_UNIT,*) lon1,IOD.hd.lon1,lon2,IOD.hd.lon2,dlon,IOD.hd.dlon
        CALL exit(1)
      END IF
  
      IF ((l-1.d0)*IOD.hd.dlat+IOD.hd.lat1 .NE. lat) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): latitude is not consistent.'
        CALL exit(1)
      END IF

      IF (IOD.hd.dlat .NE. 0) THEN
        ilat=INT((lat-IOD.hd.lat1)/IOD.hd.dlat)+1
      ELSE
        ilat=1
      END IF
      IF (IOD.hd.dhgt .NE. 0) THEN
        ihgt=INT((hgt-IOD.hd.hgt1)/IOD.hd.dhgt)+1
      ELSE
        ihgt=1
      END IF
      IOD.iondat(i).ionlat(ilat,ihgt).lat=lat
      IOD.iondat(i).ionlat(ilat,ihgt).hgt=hgt

      nline=INT(IOD.nlon/16)+1
      DO k=1, nline-1, 1
        READ(lfn,'(A)',ERR=100,END=100) line
        READ(line,'(16I5)') (IOD.iondat(i).ionlat(ilat,ihgt).tec(m),m=16*(k-1)+1,16*k)
      END DO
      READ(lfn,'(A)',ERR=100,END=100) line
      READ(line,'(16I5)') (IOD.iondat(i).ionlat(ilat,ihgt).tec(m),m=16*(nline-1)+1,IOD.nlon)

    END DO

    READ(lfn,'(A)',ERR=100,END=100) line
  END DO
  
  DO WHILE(INDEX(line,'START OF RMS MAP') .EQ. 0)
    READ(lfn,'(A)',ERR=100,END=100) line
  END DO
  BACKSPACE(lfn)

  ! Read the RMS
  DO i=1, IOD.hd.nmaps
    ! Start of RMS MAP
    READ(lfn,'(A)',ERR=100,END=100) line
    READ(line,'(I6)') l
    IF (IOD.iondat(i).ind .NE. l) THEN
      WRITE(ERROR_UNIT,'(A,1X,I5)') '***ERROR(read_ionex): START OF TEC MAP is 0 or do not equal ',i
      CALL exit(1)
    END IF

    ! Epoch of current MAP
    READ(lfn,'(A)',ERR=100,END=100) line    
    ! READ(line,'(6I6)') iy,im,id,ih,imin,isec
    READ(line,'(5I6,F6.0)') iy,im,id,ih,imin,isec
    IF (DABS(timdif(IOD.hd.mjd0,IOD.hd.sod0+(i-1)*IOD.hd.intv,IOD.iondat(i).mjd,IOD.iondat(i).sod)) .GT. MAXWND) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): epoch is not correct.'
      CALL exit(1)
    END IF

    DO l=1, IOD.nlat, 1
      ! LAT/LON/LON2/DLON/H
      READ(lfn,'(A)',ERR=100,END=100) line
      READ(line,'(2X,5F6.1)') lat,lon1,lon2,dlon,hgt
      IF (lon1.NE.IOD.hd.lon1 .OR. lon2.NE.IOD.hd.lon2 .OR. dlon.NE.IOD.hd.dlon) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): longtitude is not consistent.'
        CALL exit(1)
      END IF

      IF ((l-1.d0)*IOD.hd.dlat+IOD.hd.lat1 .NE. lat) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(read_ionex): latitude is not consistent.'
        CALL exit(1)
      END IF
  
      IF (IOD.hd.dlat.NE.0) THEN
        ilat=INT((lat-IOD.hd.lat1)/IOD.hd.dlat)+1
      ELSE
        ilat=1
      END IF
      IF (IOD.hd.dhgt.NE.0) THEN
        ihgt=INT((hgt-IOD.hd.hgt1)/IOD.hd.dhgt)+1
      ELSE
        ihgt=1
       END IF

      nline=INT(IOD.nlon/16)+1
      DO k=1, nline-1, 1
        READ(lfn,'(A)',ERR=100,END=100) line
        READ(line,'(16I5)') (IOD.iondat(i).ionlat(ilat,ihgt).rms(m),m=16*(k-1)+1,16*k)
      END DO
      READ(lfn,'(A)',ERR=100,END=100) line
      READ(line,'(16I5)') (IOD.iondat(i).ionlat(ilat,ihgt).rms(m),m=16*(nline-1)+1,IOD.nlon)

    END DO

    READ(lfn,'(A)',ERR=100,END=100) line
  END DO
  
 100 CONTINUE 
 
  CLOSE(lfn)

  RETURN

END SUBROUTINE
