!*
SUBROUTINE gim(mjd,sod,lat,lon,elev,azim,IOD,dion1,drms1)
!!
!!
!! lat : Geodetic latitude of receiver          (rad)
!! lon : Geodetic longitude of receiver         (rad)
!! elev: Elevation angle of satellite           (rad)
!! azim: Geodetic azimuth of satellite          (rad)
!!
!*
USE const
USE ion
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod,lat,lon,elev,azim,dion1,drms1
TYPE(IONEX) :: IOD

  !*
  ! The variables
  !!-----------------------
  REAL(RL) :: dion2,drms2
  INTEGER(IT) :: idx
  REAL(RL) :: lat_i,lon_i,lon_i1,lon_i2

  !*
  ! The funcation called
  !!----------------------------
  REAL(RL) :: timdif
  
  !*
  ! Start the exectuable code
  !!----------------------------
  
  ! Check the time
  IF (timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0).LT.0.d0 .OR. timdif(mjd,sod,IOD.hd.mjd1,IOD.hd.sod1).GT.0.d0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(gim): the required epoch is not in TECU file.'
    CALL exit(1)
  END IF
  
  CALL ipp(lat,lon,elev,azim,lat_i,lon_i)
  
  idx=INT(timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0)/IOD.hd.intv)+1
  
  lon_i1=lon_i+timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0+(idx-1)*IOD.hd.intv)*2*PI/86400.d0
  lon_i2=lon_i+timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0+idx*IOD.hd.intv)*2*PI/86400.d0

  ! [-pi, +pi]
  DO WHILE(lon_i1 .GT. PI)
    lon_i1=lon_i1-2*PI
  END DO
  DO WHILE(lon_i1 .LT. -PI)
    lon_i1=lon_i1+2*PI
  END DO 
  DO WHILE(lon_i2 .GT. PI)
    lon_i2=lon_i2-2*PI
  END DO
  DO WHILE(lon_i2 .LT. -PI)
    lon_i2=lon_i2+2*PI
  END DO
  
  IF (idx .GE. IOD.hd.nmaps) THEN
  
    CALL gettec(lat_i,lon_i1,idx,IOD,dion1,drms1)


  ELSE
    CALL gettec(lat_i,lon_i1,idx  ,IOD,dion1,drms1)
    CALL gettec(lat_i,lon_i2,idx+1,IOD,dion2,drms2)

    dion1 = timdif(IOD.hd.mjd0,IOD.hd.sod0+idx*IOD.hd.intv,mjd,sod)/IOD.hd.intv*dion1- &
            timdif(IOD.hd.mjd0,IOD.hd.sod0+(idx-1)*IOD.hd.intv,mjd,sod)/IOD.hd.intv*dion2

    drms1 = DSQRT((timdif(IOD.hd.mjd0,IOD.hd.sod0+idx*IOD.hd.intv,mjd,sod)/IOD.hd.intv*drms1)**2+ &
            (timdif(IOD.hd.mjd0,IOD.hd.sod0+(idx-1)*IOD.hd.intv,mjd,sod)/IOD.hd.intv*drms2)**2)
    
  END IF

  ! MSLM: VTEC to STEC
  ! IF (dion1 .NE. 9999.d0) dion1=dion1*0.1d0*6371.d0/(6371.d0+506.7d0)*DSIN(0.9782*(PI/2.d0-elev))
  ! IF (drms1 .NE.    0.d0) drms1=drms1*0.1d0*6371.d0/(6371.d0+506.7d0)*DSIN(0.9782*(PI/2.d0-elev))
  ! SLM: Differential Code Bias Estimation using MultiGNSS Observations and Global Ionosphere Maps
  IF (dion1 .NE. 9999.d0) dion1=dion1*0.1d0/(DSQRT(1-(6371.d0/(6371.d0+506.7d0)*DCOS(elev))**2))
  IF (drms1 .NE.    0.d0) drms1=drms1*0.1d0/(DSQRT(1-(6371.d0/(6371.d0+506.7d0)*DCOS(elev))**2))

  RETURN

END SUBROUTINE




!*
SUBROUTINE gim2(mjd,sod,lat_i,lon_i,IOD,dion1)
!!
!!
!! lat : Geodetic latitude of receiver          (rad)
!! lon : Geodetic longitude of receiver         (rad)
!! elev: Elevation angle of satellite           (rad)
!! azim: Geodetic azimuth of satellite          (rad)
!!
!*
USE const
USE ion
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod,lat,lon,elev,azim,dion1
TYPE(IONEX) :: IOD

  !*
  ! The variables
  !!-----------------------
  REAL(RL) :: dion2
  INTEGER(IT) :: idx
  REAL(RL) :: lat_i,lon_i,lon_i1,lon_i2

  !*
  ! The funcation called
  !!----------------------------
  REAL(RL) :: timdif
  
  !*
  ! Start the exectuable code
  !!----------------------------
  
  ! Check the time
  IF (timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0).LT.0.d0 .OR. timdif(mjd,sod,IOD.hd.mjd1,IOD.hd.sod1).GT.0.d0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(gim): the required epoch is not in TECU file.'
    CALL exit(1)
  END IF
  
  idx=INT(timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0)/IOD.hd.intv)+1

  lon_i1=lon_i+timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0+(idx-1)*IOD.hd.intv)*2*PI/86400.d0
  lon_i2=lon_i+timdif(mjd,sod,IOD.hd.mjd0,IOD.hd.sod0+idx*IOD.hd.intv)*2*PI/86400.d0

  ! [-pi, +pi]
  DO WHILE(lon_i1 .GT. PI)
    lon_i1=lon_i1-2*PI
  END DO
  DO WHILE(lon_i1 .LT. -PI)
    lon_i1=lon_i1+2*PI
  END DO 
  DO WHILE(lon_i2 .GT. PI)
    lon_i2=lon_i2-2*PI
  END DO
  DO WHILE(lon_i2 .LT. -PI)
    lon_i2=lon_i2+2*PI
  END DO
  
  IF (idx .GE. IOD.hd.nmaps) THEN
  
    CALL gettec(lat_i,lon_i1,idx,IOD,dion1)


  ELSE
    CALL gettec(lat_i,lon_i1,idx  ,IOD,dion1)
    CALL gettec(lat_i,lon_i2,idx+1,IOD,dion2)

    dion1 = timdif(IOD.hd.mjd0,IOD.hd.sod0+idx*IOD.hd.intv,mjd,sod)/IOD.hd.intv*dion1- &
            timdif(IOD.hd.mjd0,IOD.hd.sod0+(idx-1)*IOD.hd.intv,mjd,sod)/IOD.hd.intv*dion2
    
  END IF

  ! MSLM: VTEC to STEC
  IF (dion1 .NE. 9999.d0) dion1=dion1*0.1d0*6371.d0/(6371.d0+506.7d0)*DSIN(0.9782*(PI/2.d0-elev))

  RETURN

END SUBROUTINE

SUBROUTINE gettec(lat,lon,idx,IOD,dion1,drms1)
!!
!! lat : Geodetic latitude of receiver          (rad)
!! lon : Geodetic longitude of receiver         (rad)
!! elev: Elevation angle of satellite           (rad)
!! azim: Geodetic azimuth of satellite          (rad)
!!
!*
USE const
USE ion
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: idx
REAL(RL) :: lat,lon,elev,azim,dion1,drms1
TYPE(IONEX) :: IOD

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: ilat1,ilon1,ilat2,ilon2
  
  REAL(RL) :: p,q,a11,a21,a12,a22
  REAL(RL) ::     b11,b21,b12,b22
  
  !*
  ! The funcation called
  !!----------------------------
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lon*RAD2DEG.LT.-180.d0 .OR. lon*RAD2DEG.GT.180.d0) THEN
    WRITE(OUTPUT_UNIT,'(A)') '###WARNING(gettec): longitude is out of [-180,180]'
    dion1=9999.d0
    drms1=0.d0
  END IF
  
  IF (lat*RAD2DEG.LT.-(90.d0-DABS(IOD.hd.dlat)) .OR. lat*RAD2DEG.GT.(90.d0-DABS(IOD.hd.dlat))) THEN
    WRITE(OUTPUT_UNIT,'(A)') '###WARNING(gettec): latitude is out of [-87.5,87.5].'
    dion1=9999.d0
    drms1=0.d0
  END IF
  
  IF (idx .GT. IOD.hd.nmaps) THEN
    WRITE(OUTPUT_UNIT,'(A)') '###WARNING(gettec): number of maps is greater then tha maximum.'
    dion1=9999.d0  
    drms1=0.d0  
  END IF
  
  ! LON1/LON2 in [0, 360]
  IF (IOD.hd.lon2 .GT. 180.d0) THEN
    IF (lon*RAD2DEG .LT. 0.d0) THEN
      ilon1=INT((lon*RAD2DEG+360.d0-IOD.hd.lon1)/IOD.hd.dlon)+1
      p=(lon*RAD2DEG+360.d0-IOD.hd.lon1-(ilon1-1)*IOD.hd.dlon)/IOD.hd.dlon
    ELSE
      ilon1=INT((lon*RAD2DEG-IOD.hd.lon1)/IOD.hd.dlon)+1
      p=(lon*RAD2DEG-IOD.hd.lon1-(ilon1-1)*IOD.hd.dlon)/IOD.hd.dlon
    END IF
  ELSE
    ilon1=INT((lon*RAD2DEG-IOD.hd.lon1)/IOD.hd.dlon)+1
    p=(lon*RAD2DEG-IOD.hd.lon1-(ilon1-1)*IOD.hd.dlon)/IOD.hd.dlon
  END IF
  ilon2=ilon1+1

  ilat1=INT((lat*RAD2DEG-IOD.hd.lat1)/IOD.hd.dlat)+1
  ilat2=ilat1+1

  q=(lat*RAD2DEG-IOD.hd.lat1-(ilat1-1)*IOD.hd.dlat)/IOD.hd.dlat
  
  a11=IOD.iondat(idx).ionlat(ilat1,1).tec(ilon1)
  b11=IOD.iondat(idx).ionlat(ilat1,1).rms(ilon1)
  IF (ilon2 .GE. IOD.nlon) THEN
    a21=9999.d0
    b21=0.d0
  ELSE
    a21=IOD.iondat(idx).ionlat(ilat1,1).tec(ilon2)
    b21=IOD.iondat(idx).ionlat(ilat1,1).rms(ilon2)
  END IF
  
  IF (ilat2 .GE. IOD.nlat) THEN
    a12=9999.d0
    b12=0.d0
  ELSE
    a12=IOD.iondat(idx).ionlat(ilat2,1).tec(ilon1)
    b12=IOD.iondat(idx).ionlat(ilat2,1).rms(ilon1)
  END IF
  IF (a12.EQ.9999.d0 .OR. a21.EQ.9999.d0) THEN
    a22=9999.d0
  ELSE
    a22=IOD.iondat(idx).ionlat(ilat2,1).tec(ilon2)  
  END IF
  IF (b12.EQ.0.d0 .OR. b21.EQ.0.d0) THEN
    b22=0.d0
  ELSE
    b22=IOD.iondat(idx).ionlat(ilat2,1).rms(ilon2)  
  END IF
    
  IF (a12.EQ.9999.d0 .OR. a21.EQ.9999.d0 .OR. a21.EQ.9999.d0 .OR. a22.EQ.9999.d0) THEN
    dion1=9999.d0
  ELSE
    dion1=(1-p)*(1-q)*a11+(1-p)*q*a12+p*(1-q)*a21+p*q*a22
  END IF

  IF (b12.EQ.0.d0 .OR. b21.EQ.0.d0 .OR. b21.EQ.0.d0 .OR. b22.EQ.0.d0) THEN
    drms1=0.d0
  ELSE
    drms1=DSQRT(((1-p)*(1-q)*b11)**2+((1-p)*q*b12)**2+(p*(1-q)*b21)**2+(p*q*b22)**2)
  END IF

  RETURN

END SUBROUTINE
