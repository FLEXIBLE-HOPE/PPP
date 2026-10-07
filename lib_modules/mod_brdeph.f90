!*
MODULE brdeph
!*
USE const
USE par
IMPLICIT NONE

  ! GNSS NAVIGATION MESSAGE FILE - HEADER SECTION
  TYPE BRDHEAD
    REAL(RL) :: ver = 0.d0
    ! Ionospheric correction parameters
    ! correction type, alpha 1, beta 2
    CHARACTER(LEN_EPHION) :: ionc(2,MAXSYS) = ''
    ! parameters
    REAL(RL) :: ion(4,2,MAXSYS) = 0.d0
    ! The following two parameters are mandatory for BDS
    ! time mark, transmission time (seconds of week)
    CHARACTER :: tmark(2,MAXSYS)
    ! SV ID, identify which satellite provided the ionospheric parameters
    CHARACTER :: svid(2,MAXSYS)
    ! Corrections to transform the system time to UTC or other time systems
    ! correction type, UTC 1, other system 2
    CHARACTER(LEN_EPHION) :: timc(2,MAXSYS) = ''
    ! parameters, 1-2 parameters, 3-4 reference time for polynomial
    REAL(RL) :: tim(4,2,MAXSYS) = 0.d0
    ! satellite broadcasting
    CHARACTER(5) :: satb(2,MAXSYS)
    ! UTC identifier
    INTEGER(IT) :: utcid(2,MAXSYS)
    ! Number of leap seconds since 6-Jan-1980
    INTEGER(IT) :: leap(4) = 0
    ! Time system identifier, only GPS or BDS
    CHARACTER(3) :: lpts
  END TYPE

  ! Brodcast ephemeris type for GPS, Galieo, Beidou, QZSS, and IRNSS
  TYPE GPS_BRDEPH
  !SEQUENCE
    CHARACTER(LEN_PRN) :: cprn = ''
    ! jd, sod : Time of Clock (Toc)
    INTEGER(IT) :: mjd = 0
    REAL(RL)    :: sod = 0.D0
    ! a0: SV clock offset
    ! a1: SV clock drift
    ! a2: SV clock drift rate
    REAL(RL) :: a0 = 0.D0,a1 = 0.D0, a2 = 0.D0
    ! ORBIT-1
    ! aode: age of ephemeris upload (IODE) 
    ! crs: Orbtital radius correction
    ! dn: Mean motion difference
    ! m0: Mean anomaly at reference epoch
    REAL(RL) :: aode = 0.D0, crs = 0.D0, dn = 0.D0, m0=0.d0
    ! ORBIT-2
    ! e: Eccentricity
    ! cuc, cus: Latitude argument correction
    ! roota: Square root of semi-major axis
    REAL(RL) :: cuc = 0.D0, e = 0.D0, cus = 0.D0, roota = 0.D0
    ! ORBIT-3
    ! toe, week: Ephemerides reference epoch in seconds with the week
    ! cis, cic: Inclination correction
    ! omega0: Longtitude of ascending node at the begining of the week
    REAL(RL) :: toe = 0.D0, cic = 0.D0, omega0 = 0.D0, cis = 0.D0
    ! ORBIT-4
    ! i0: Inclination at reference epoch
    ! crc: Orbtital radius correction
    ! omega: Argument of perigee
    ! omegadot: Rate of node's right ascension
    REAL(RL) :: i0 = 0.D0, crc = 0.D0, omega = 0.D0, omegadot = 0.D0
    ! ORBIT-5
    ! idot: Rate of inclination angle
    ! sesvd0: Codes on L2 channel (GPS,QZSS), Data sources (Galileo, see RINEX3.03), Spare (BDS,IRNSS)
    ! week: GPS Week #(to go with ToE)
    ! resvd1: L2P data flag (GPS,QZSS), spare (Galileo, BDS, IRNSS)
    REAL(RL) :: idot = 0.D0, resvd0 = 0.D0, week = 0.D0, resvd1 = 0.D0
    ! ORBIT-6
    ! accu: SV accuracy (GPS,QZSS), SISA (signal in space accuracy, Galileo), User Range Accuray (URA, IRNSS)
    ! hlth: SV health, for Galileo should converted to INTEGER
    ! tgd: Time group delay (GPS,QZSS), BGD E5a/E1 (Galileo), TGD1 B1/B3 (BDS)
    ! aodc: Age of clock parameter upload (GPS,QZSS), BDG E5b/E1 (Galileo), TGD2 B2/B3 (BDS), black (IRNSS)
    REAL(RL) :: accu = 0.D0, hlth = 0.D0, tgd = 0.D0, aodc = 0.D0
    ! ORBIT-7
    ! tom: Transmission time of message
    ! fih: Fit interval in hours (GPS,QZSS), spare (Galileo, IRNSS), AODC age of data clock (BDS)
    REAL(RL) :: tom = 0.D0
    REAL(RL) :: fih = 0.d0
  END TYPE GPS_BRDEPH
  
  ! Brodcast ephemeris type for GLONASS, SBAS(maybe)
  TYPE GLONASS_BRDEPH
  !SEQUENCE
    ! Slot. number in satellite constellation
    CHARACTER(LEN_PRN) :: cprn = ''
    ! jd, sod : Time of Clock
    INTEGER(IT) :: mjd = 0
    REAL(RL)    :: sod = 0.D0
    ! tau: SV clock bias
    REAL(RL) :: tau = 0.D0
    ! SV relative frequency bias
    REAL(RL) :: gamma = 0.D0
    ! Message frame time (tk+nd*86400) in seconds of the UTC week (GLONASS)
    ! Transmission time of message (start of the message) in GPS seconds of the week (SBAS)
    REAL(RL) :: tk = 0.D0
    ! ORBIT-1
    ! Coordinate at ephemerides reference epoch in PZ-90
    REAL(RL) :: pos(3) = 0.D0
    ! SV health
    REAL(RL) :: health = 0.D0
    ! ORBIT-2
    ! Velocity at ephemerides reference epoch in PZ-90
    REAL(RL) :: vel(3) = 0.D0
    ! frequency number (GLONASS)
    ! Accuracy code (URA, SBAS)
    REAL(RL) :: frenum = 0.D0
    ! ORBIT-3
    ! Acceleration at ephemerides reference epoch in PZ-90
    REAL(RL) :: acc(3) = 0.D0
    ! age of operation information
    REAL(RL) :: age
  END TYPE GLONASS_BRDEPH

  ! Clock coefficients of navigation message for RT clock determination
  TYPE BRDCK
    ! The reference time of clock coefficients
    INTEGER(IT) :: mjd
    REAL(RL) :: sod
    !  The clock offset, rate, and acceleration
    REAL(RL) :: a(0:2)
  END TYPE BRDCK

END MODULE brdeph
