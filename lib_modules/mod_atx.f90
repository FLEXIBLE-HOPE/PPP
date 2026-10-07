!*
MODULE atx
!*
USE par
IMPLICIT NONE

  TYPE antatx
    INTEGER(IT) :: nfreq = 0
    CHARACTER(LEN_FREQ) :: freq(MAXFREQ,MAXSYS) = ''
    CHARACTER(LEN_ANTENNA) :: antnam = ''
    CHARACTER(LEN_ANTENNA) :: antnum = ''
    REAL(RL) :: mjd1 = 0.D0
    REAL(RL) :: mjd2 = 0.D0
    REAL(RL) :: zen1 = 0.D0
    REAL(RL) :: zen2 = 0.D0
    REAL(RL) :: dzen = 0.D0
    REAL(RL) :: dazi = 0.D0
    REAL(RL) :: neu(3,MAXFREQ,MAXSYS) = 0.D0
    !REAL(RL) :: pcv(50,0:200,MAXFREQ,MAXSYS) = 0.D0
    ! For GALILEO
    REAL(RL) :: pcv(50,0:200,MAXFREQ,MAXSYS) = 0.D0
    !REAL(RL) :: pcv(50,0:80,MAXFREQ,MAXSYS) = 0.D0
    !REAL(RL) :: pcv(30,0:80,MAXFREQ,MAXSYS) = 0.D0
  END TYPE

END MODULE
