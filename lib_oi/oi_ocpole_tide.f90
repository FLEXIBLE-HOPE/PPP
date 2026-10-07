!*
SUBROUTINE oi_ocpole_tide(mjd,xpole,ypole,ndegree,dc,ds)
!!
!! To compute the effect of the ocean pole tide according to
!! IERS 2010 P94
!*
USE tables
IMPLICIT NONE

!*
! The arguments
!!--------------------------
INTEGER(IT) :: ndegree
REAL(RL) :: dc(MAXOPTDEG,0:MAXOPTDEG)
REAL(RL) :: ds(MAXOPTDEG,0:MAXOPTDEG)
REAL(RL) :: mjd,xpole,ypole

  !*
  ! The local variables
  !!----------------------------

  INTEGER(IT) :: lfn,n,m

  REAL(RL), PARAMETER :: MAXOPDEG = 360
  REAL(RL) :: kn(MAXOPDEG)
  REAL(RL) :: Re_Anm(MAXOPDEG,0:MAXOPDEG)
  REAL(RL) :: Im_Anm(MAXOPDEG,0:MAXOPDEG)
  REAL(RL) :: Re_Bnm(MAXOPDEG,0:MAXOPDEG)
  REAL(RL) :: Im_Bnm(MAXOPDEG,0:MAXOPDEG)

  REAL(RL) :: fac1,fac2
  REAL(RL) :: xpm,ypm,fac,rn

  CHARACTER(LEN_STRING) :: line

  LOGICAL(LG) :: lfirst,lexist
  DATA lfirst /.TRUE./

  SAVE lfirst,Re_Anm,Im_Anm,Re_Bnm,Im_Bnm,kn,fac

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    fac=(E_ROTATE**2)*(E_MAJAXIS**4)/GME*4*PI*6.67428D-11*1.025D+03/9.7803278D0

    kn=0.d0
    kn(2)=-0.3075d0
    kn(3)=-0.1950d0
    kn(4)=-0.1320d0
    kn(5)=-0.1032d0
    kn(6)=-0.0892d0

    line=f_tableFileName("optide")
    lfn=get_valid_unit(10)
    INQUIRE(FILE=line,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(ocean_pole_tide): the ocean tide file is not exist.'
      CALL exit
    END IF
    OPEN(UNIT=lfn,FILE=line,STATUS='OLD')
    ! skype the first line
    READ(lfn,*)
    DO WHILE(.TRUE.)
      READ(lfn,'(A)',END=100) line
      READ(line,*) n,m
      IF (n.GE.MAXOPDEG .OR. m.GT.MAXOPDEG) CYCLE
      READ(line,*) n,m,Re_Anm(n,m),Re_Bnm(n,m),Im_Anm(n,m),Im_Bnm(n,m)
    END DO
100 CONTINUE
    CLOSE(lfn)
  END IF

  CALL mean_pole('IERS2010',mjd,xpm,ypm)
  xpm=xpole-xpm
  ypm=-(ypole-ypm)

  !! arcsec to radian
  xpm=xpm*ARCSEC2RAD
  ypm=ypm*ARCSEC2RAD

  fac1=xpm*0.6870d0+ypm*0.0036d0
  fac2=ypm*0.6870d0-xpm*0.0036d0

  DO n=1, ndegree
    rn=fac*(1.0d0+kn(n))/DBLE(2*n+1)
    DO m=0, n
      dc(n,m)=rn*(Re_Anm(n,m)*fac1+Im_Anm(n,m)*fac2)
      ds(n,m)=rn*(Re_Bnm(n,m)*fac1+Im_Bnm(n,m)*fac2)
    END DO
  END DO

  RETURN

END SUBROUTINE
