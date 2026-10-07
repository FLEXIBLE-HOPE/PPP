!*
SUBROUTINE oceanload_coef(lat,lon,olc)
!!
!! purpose  : get oceanload coefficients
!! parameter:
!!    input : lat,lon -- geodetic position of station
!!    output: olc     -- oceanload coefficients
!! author   : Geng J
!!
!*
USE par
USE const
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------------
REAL(RL) :: olc(11,6),lat,lon

  !*
  ! The local variables
  !!--------------------------
  LOGICAL(lg) :: lfirst
  INTEGER(IT) :: i,j,lfn,ierr
  CHARACTER(LEN_STRING) :: line
  REAL(RL) :: lat_deg,lon_deg,tol_d,tol_lat,tol_lon,dlat,dlon

  CHARACTER(LEN_FILENAME) :: oceanload

  DATA lfirst/.TRUE./
  SAVE lfirst,lfn,tol_lat

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit


  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst) THEN
    lfirst=.FALSE.
    oceanload = f_tablefilename('ocload')
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=oceanload,STATUS='old',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oceanload_coef): open file oceanload '
      CALL exit(1)
    END IF
    tol_d=1.d4
    tol_lat=tol_d/E_MAJAXIS
  END IF
  tol_lon=tol_lat/dcos(lat)

  !! look for a new one
  DO j=1,6
    DO i=1,11
      olc(i,j)=0.d0
    END DO
  END DO

  REWIND(lfn)

  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=100) line
    i=INDEX(line(1:LEN_TRIM(line)),'lon/lat:')
    IF (i .NE. 0) THEN
      READ(line(i+8:),*) lon_deg,lat_deg
      IF (lon_deg .LT. 0.d0) lon_deg=lon_deg+360.d0
      dlon=dabs(lon_deg*DEG2RAD-lon)
      dlat=dabs(lat_deg*DEG2RAD-lat)
      IF (dlat.LE.tol_lat .AND. dlon.LE.tol_lon) THEN
        DO j=1,6
          READ(lfn,*,iostat=ierr) (olc(i,j),i=1,11)
          IF (ierr .NE. 0) THEN
            WRITE(ERROR_UNIT,'(A)') '***ERROR(oceanload_coef): read file oceanload '
            CALL exit(1)
          END IF

          IF (j .GE. 4) THEN
            DO i=1,11
              olc(i,j)=olc(i,j)*deg2rad
            END DO
          END IF
        END DO
        RETURN
      END IF
    END IF
  END DO

100 WRITE(OUTPUT_UNIT,'(A)') '###WARNING(oceanload_coef): no oceanload coefficients'

  RETURN

END SUBROUTINE
