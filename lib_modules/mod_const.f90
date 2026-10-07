!*
MODULE const
!!
!! POSITION AND NAVIGATION DATA ANALYST (PANDA) SOFTWARE HAS BEEN DEVELOPED BY
!! GNSS RESEARCH CENTER, WUHAN UNIVERSITY SINCE 2000
!!
!! THIS VERSION HAS BEEN MAINTAINED BY JING GUO SINCE 2013, SO IT IS
!! JING GUO'S POSITION AND NAVIGATION DATA ANALYST (gPANDA).
!!
!! VERSION 1.0 (2013/07/09): INITIAL VERSION IS FOCUSED ON PRECISE ORBIT DETERMINATION
!!   FOR MULTI-GNSS SATELLITES (GPS, GLONASS, BEIDOU, GALILEO, QZSS, AND SO ON)
!!
!! VERSION 1.2 (2015/10/17): SUPPORT THE LEO SATELLITES, AND LEOs COULD ALSO USED
!!   AS NAVIGATION SATELLITES
!!
!! VERSION 1.3 (2016/08/24): ADD BEIDOU 3
!!
!! VERSION 1.31 (2017/02/18): UPDATE CKDCTRL FOR IONOSPHERE

!! CHANGES LOG:
!!   April 18, 2014: Add yaw attitude for GLONASS-M, the test results show
!!     improved orbit.
!!
!!   April 19, 2014: Estimated IFB for each GLONASS satellites, the test
!!     results demonstrate the orbit improved a bit once the the IFB
!!     estimated at last.
!!
!!   April 20-23, 2014: fixed bug to ISB estimation, due to add the 'sig'
!!     to the relitavity corrections.
!!
!!   April 24, 2014: there is noticable periodic signal in Galielo yaw angle,
!!     maybe related to clock correction, but still unsolved.
!!     The tests shows that the GPS orbit and clock products achieve same
!!     accuracy as CODE, and GFZ for GPS.
!!
!!   April 25, 2014: As to discussion with Bin Wang, the periodical signal
!!     in satellite clock is related to thermal effect, so the problem could
!!     be omitted.
!!
!!   May 28, 2014: fix minor bug in sp3orb and orb2sp3
!!
!!   4 Jun, 2014: fix minor bug in mkclk
!!      update the cc2noncc to support RINEX 3 file
!!
!!   12 Jun, 2014: extend the SLR program to estimate the orbit
!!
!!   13 Jun, 2014: add the PCV estimation in LSQ
!!      and also update the PCO and PCV estimation for each frequency
!!
!!   20 Jun, 2014: update the ambfix program to fix the ambiguity for
!!      raw observations
!!
!!   31 Oct, 2015: add the LEO for data processing, fixed the bug related to
!!      compute the mean pole. The transformation from ITRF to GCRF needs to
!!      be revised for velocity. The eclipse has been modified for QZSS
!!
!!   01 Jun, 2018: changed svnav and atx to IGS convension
!!      Add IRNSS and more QZSS satellites
!*
USE par
IMPLICIT NONE

  CHARACTER(LEN=50) :: VERSION = " JING GUO's PANDA(@gPANDA) Version: 1.3 (20170218)"

  REAL(RL), PARAMETER :: PI = 3.14159265358979323846D0
  REAL(RL), PARAMETER :: RAD2DEG = 180.D0/PI
  REAL(RL), PARAMETER :: DEG2RAD = PI/180.D0
  REAL(RL), PARAMETER :: ARCSEC2RAD = PI/(3600.D0*180.D0)
  REAL(RL), PARAMETER :: SEC2RAD = PI/43200.D0

  !! Groten E (2004) Fundamental parameters and current (2004) best estimates of the parameters
  !! of common relevance to astronomy, geodesy, and geodynamics. Journal of Geodesy
  !! IERS 2010, TT-compatible, not TCG, TDG, TCB, TGD
  REAL(RL), PARAMETER :: VEL_LIGHT = 299792458.D0 ! defining, same for all time system
  REAL(RL), PARAMETER :: GME = 3.986004415D14     ! TT-compatible, same as CODE
  REAL(RL), PARAMETER :: GMS = 1.327124420065D20  ! TT-compatible, CODE 1.3271245D20
  REAL(RL), PARAMETER :: E_ROTATE = 7.2921151467D-5 ! defining, 1.00273781191135448 rev/UT1day
  REAL(RL), PARAMETER :: E_MAJAXIS = 6378136.3D0  ! 6378136.3 (6378136.6 old) TT-compatible, and most geopotential model used
  REAL(RL), PARAMETER :: MEMR = 0.0123000371d0    ! MOON-Earth Mass Ratio, the difference between TT and TCG is quite small

  REAL(RL), PARAMETER :: OFF_GPS2TAI = 19.D0
  REAL(RL), PARAMETER :: OFF_TAI2TT = 32.184D0
  REAL(RL), PARAMETER :: OFF_GPS2TT = OFF_GPS2TAI+OFF_TAI2TT
  REAL(RL), PARAMETER :: OFF_MJD2JD = 2400000.5D0

  CHARACTER(MAXSYS) :: SYS = 'GRECSJIL'
  !! for iGMAS due to W could be tracked by much more sites, however the W is too bad
  CHARACTER(LEN=16) :: OBSTYPE = 'PWCIXSAQLDBYMZN '
  !CHARACTER(LEN=16) :: OBSTYPE = 'PWCDYMNABCIQSLXZ'
  ! CHARACTER(LEN=16) :: OBSTYPE = 'CPWIXSAQLDBYMZN '

  ! Frequency for GNSS satellites
  REAL(RL), PARAMETER :: GPS_FREQ = 10230000.D0
  REAL(RL), PARAMETER :: GPS_L1 = 154*GPS_FREQ
  REAL(RL), PARAMETER :: GPS_L2 = 120*GPS_FREQ
  REAL(RL), PARAMETER :: GPS_L5 = 115*GPS_FREQ

  REAL(RL), PARAMETER :: GLS_FREQ = 178000000.D0
  REAL(RL), PARAMETER :: GLS_L1 = 9*GLS_FREQ
  REAL(RL), PARAMETER :: GLS_L2 = 7*GLS_FREQ
  REAL(RL), PARAMETER :: GLS_DL1 = 562500.D0
  REAL(RL), PARAMETER :: GLS_DL2 = 437500.D0
  REAL(RL) :: GLS_FAC(2)
  DATA GLS_FAC /0.187136366d0, 0.240603899d0/

  REAL(RL), PARAMETER :: GAL_E1  = GPS_L1
  REAL(RL), PARAMETER :: GAL_E5  = 1191795000.D0
  REAL(RL), PARAMETER :: GAL_E5A = GPS_L5
  REAL(RL), PARAMETER :: GAL_E5B = 1207140000.D0
  REAL(RL), PARAMETER :: GAL_E6  = 1278750000.D0

  REAL(RL), PARAMETER :: BDS_B1 = 1561098000.D0
  REAL(RL), PARAMETER :: BDS_B2 = GAL_E5B
  REAL(RL), PARAMETER :: BDS_B3 = 1268520000.D0

  REAL(RL), PARAMETER :: BDS_B1C = GPS_L1
  REAL(RL), PARAMETER :: BDS_B1I = BDS_B1
  REAL(RL), PARAMETER :: BDS_B2A = GPS_L5
  REAL(RL), PARAMETER :: BDS_B2B = GAL_E5B       !B2B
  REAL(RL), PARAMETER :: BDS_B2C = GAL_E5        !B2a+b
  REAL(RL), PARAMETER :: BDS_B3I = BDS_B3
  REAL(RL), PARAMETER :: BDS_BS  = 2492028000.D0

  REAL(RL), PARAMETER :: QZS_L1  = GPS_L1
  REAL(RL), PARAMETER :: QZS_L2  = GPS_L2
  REAL(RL), PARAMETER :: QZS_L5  = GPS_L5
  REAL(RL), PARAMETER :: QZS_LEX = 1278750000.D0

  REAL(RL), PARAMETER :: IRS_L5  = GPS_L5
  REAL(RL), PARAMETER :: IRS_S   = 2492028000.D0

  REAL(RL), PARAMETER :: LEO_L1  = GPS_L1
  REAL(RL), PARAMETER :: LEO_L2  = GPS_L2
  REAL(RL), PARAMETER :: LEO_L5  = GPS_L5

END MODULE
