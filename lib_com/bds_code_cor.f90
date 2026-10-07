!*
SUBROUTINE bds_code_cor(ctype,nfreq,cfreq,elev,bias)
!!
!! Correction of code pseudorange measurements based on the model presented
!! Lambert Wanninger, Susanne Beer. BeiDou satellite-induced code pseudorange
!! variations: diagnosis and therapy. GPS Solut DOI 10.1007/s10291-014-0423-3
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
CHARACTER(LEN=*) :: ctype
CHARACTER(LEN_FREQ) :: cfreq(MAXFREQ)
INTEGER(IT) :: nfreq
REAL(RL) :: elev, bias(1:*)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: isat,ifreq(3),i,idx

  REAL(RL) :: coef(10,3,2),alpha
  DATA coef  /-0.55d0, -0.40d0, -0.34d0, -0.23d0, -0.15d0, &
              -0.04d0,  0.09d0,  0.19d0,  0.27d0,  0.35d0, &
              -0.71d0, -0.36d0, -0.33d0, -0.19d0, -0.14d0, &
              -0.03d0,  0.08d0,  0.17d0,  0.24d0,  0.33d0, &
              -0.27d0, -0.23d0, -0.21d0, -0.15d0, -0.11d0, &
              -0.04d0,  0.05d0,  0.14d0,  0.19d0,  0.32d0, &
              -0.47d0, -0.38d0, -0.32d0, -0.23d0, -0.11d0, &
               0.06d0,  0.34d0,  0.69d0,  0.97d0,  1.05d0, &
              -0.40d0, -0.31d0, -0.26d0, -0.18d0, -0.06d0, &
               0.09d0,  0.28d0,  0.48d0,  0.64d0,  0.69d0, &
              -0.22d0, -0.15d0, -0.13d0, -0.10d0, -0.04d0, &
               0.05d0,  0.14d0,  0.27d0,  0.36d0,  0.47d0 /

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (nfreq .GT. 3) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(bds_code_cor): the number of frequencies is greater than 3 for BeiDou'
    CALL exit(1)
  END IF

  SELECT CASE(ctype)
    CASE('BEIDOU-2I')
      isat=1
    CASE('BEIDOU-2M')
      isat=2
    CASE DEFAULT
      RETURN
  END SELECT

  ifreq=0
  DO i=1, nfreq
    SELECT CASE(cfreq(i))
      CASE('L2')
        ifreq(i)=1
      CASE('L7')
        ifreq(i)=2
      CASE('L6')
        ifreq(i)=3
    END SELECT
  END DO

  IF (elev .LE. 0.d0) THEN
    DO i=1, nfreq
      bias(i)=coef(1,ifreq(i),isat)
    END DO
  ELSE IF (elev .GE. 90.d0) THEN
    DO i=1, nfreq
      bias(i)=coef(10,ifreq(i),isat)
    END DO
  ELSE
    idx=INT(elev/10.d0)+1
    DO i=1, nfreq
      alpha=(coef(idx+1,ifreq(i),isat)-coef(idx,ifreq(i),isat))/10.d0
      bias(i)=alpha*(elev-(idx-1)*10.d0)+coef(idx,ifreq(i),isat)
    END DO
  END IF

  !*
  ! Return
  !!-----------------------------

  RETURN

END SUBROUTINE
