!*
SUBROUTINE mean_pole(conv,mjd,xpm,ypm)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------
CHARACTER(LEN=*) :: conv
REAL(RL) :: mjd ! in TT
REAL(RL) :: xpm,ypm

  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: n,m
  REAL(RL) :: doy

  REAL(RL) :: xcof(4,2),ycof(4,2)
  DATA xcof / 55.974d0, 1.8243d0, 0.18413d0, 0.007024d0, &
              23.513d0, 7.6141d0, 0.00000d0, 0.000000d0 /
  DATA ycof /346.346d0, 1.7896d0,-0.10729d0,-0.000908d0, &
             358.891d0,-0.6287d0, 0.00000d0, 0.000000d0 /

  !*
  ! Start the exectuable code
  !!---------------------------

  xpm=0.d0
  ypm=0.d0

  !! J2000.0: http://baike.baidu.com/view/2171240.html
  !! Julian day: http://baike.baidu.com/view/37497.html
  SELECT CASE(TRIM(conv))
    CASE('IERS2003')
      xpm=0.054d0+0.00083d0*(mjd-51544.5d0)/365.25d0
      ypm=0.357d0+0.00395d0*(mjd-51544.5d0)/365.25d0
    CASE('IERS2010')
      n=1
      IF (mjd .GE. 55197.d0) n=2
      doy=(mjd-51544.5d0)/365.25d0
      DO m=1, 4
        xpm=xpm+(doy**(m-1))*xcof(m,n)
        ypm=ypm+(doy**(m-1))*ycof(m,n)
      END DO
      ! mas to arcsec
      xpm=xpm*1.d-3
      ypm=ypm*1.d-3
      ! linear mean pole for updated IERS2010
!      doy=(mjd-51544.5d0)/365.25d0
!      xpm= 55.0+1.677*doy
!      ypm=320.5+3.460*doy
!      ! mas to arcsec
!      xpm=xpm*1.d-3
!      ypm=ypm*1.d-3
    CASE DEFAULT
      WRITE(ERROR_UNIT,'(A)') '***ERROR(mean_pole): unknown conventions '//TRIM(conv)
      CALL exit(1)
  END SELECT

  RETURN

END SUBROUTINE
