!*
MODULE station
USE CONST
IMPLICIT NONE

  TYPE site
    LOGICAL(LG) :: lfcb = .FALSE.
    CHARACTER(LEN_SITENAME) :: name = ''
    CHARACTER(LEN_SITETYPE) :: skd = ''
    CHARACTER(LEN_CLKTYPE) :: clk = ''
    CHARACTER(LEN_PCVTYPE) :: pcv = ''
    CHARACTER(LEN_PRN) :: cprn=''

    REAL(RL) :: x(3) = 0.D0
    REAL(RL) :: dx0(9) = 0.D0
    REAL(RL) :: qx(3) = 0.D0
    REAL(RL) :: geod(3) = 0.D0

    REAL(RL) :: rot_l2f(3,3) = 0.D0
    REAL(RL) :: enu(3,MAXFREQ,MAXSYS) = 0.D0 !SITE PCO,站心系
    REAL(RL) :: enu0(3) = 0.D0               !SITE ENU
    REAL(RL) :: denu0(3) = 0.D0
    ! north orientation correction angle
    REAL(RL) :: noca=0.d0

    REAL(RL) :: rnxver = 0.D0 ! 3.04
    CHARACTER(LEN_ANTENNA) :: rectyp = ''
    CHARACTER(LEN_ANTENNA) :: recver = ''
    CHARACTER(LEN_ANTENNA) :: anttyp = ''

    INTEGER(IT) :: ileo = 0
    INTEGER(IT) :: iptatx = 0

    INTEGER(IT) :: jdlog = 0
    REAL(RL) :: sodlog = 0
    REAL(RL) :: intvlog = 0
    INTEGER(IT) :: iunit = 0
    CHARACTER(LEN_FILENAME) :: obsfile = ''
    INTEGER(IT) :: lfnlog = 0
    CHARACTER(LEN_FILENAME) :: logfile = ''

    REAL(RL) :: inclk = 0.D0
    REAL(RL) :: rclock(MAXSYS) = 0.D0
    REAL(RL) :: dclk0(MAXSYS) = 0.D0
    REAL(RL) :: qclk(MAXSYS) = 0.D0

    REAL(RL) :: sigr(MAXSYS) = 0.D0
    REAL(RL) :: sigp(MAXSYS) = 0.D0
    REAL(RL) :: cutoff = 0.D0

    LOGICAL(LG) :: first(MAXSAT) = .TRUE.
    REAL(RL) :: prephi(MAXSAT) = 0.D0

    REAL(RL) :: olc(11,6) = 0.D0

    CHARACTER(LEN_ZTDMAP) :: map = ''
    REAL(RL) :: ztdcor = 0.D0  ! the residual of ZWD, 估计湿延迟残余量
    REAL(RL) :: dztd0 = 0.D0
    REAL(RL) :: qztd = 0.D0
    REAL(RL) :: p0 = 0.D0
    REAL(RL) :: t0 = 0.D0
    REAL(RL) :: dt0 = 0.D0
    REAL(RL) :: hr0 = 0.D0
    REAL(RL) :: undu = 0.D0

    REAL(RL) :: zcor = 0.D0
    REAL(RL) :: zdd = 0.D0
    REAL(RL) :: zwd = 0.D0

    REAL(RL) :: dgrd0 = 0.D0
    REAL(RL) :: qgrd = 0.D0
    REAL(RL) :: grd(2) = 0.D0
    REAL(RL) :: gcor(2) = 0.D0

    REAL(RL) :: dion0 = 0.d0
    REAL(RL) :: qion = 0.d0
    REAL(RL) :: ion(MAXSAT) = 0.d0
    REAL(RL) :: sisre(MAXSAT) = 0.d0
    REAL(RL) :: iondel(MAXSAT,MAXFREQ) = 0.d0
    REAL(RL) :: trpdel(MAXSAT) = 0.d0
    !xsy
    REAL(RL) :: dop(3) = 0.d0
    !xsy-2023-12-14: 测站每个系统的每个频点的测距码(C1W)
    CHARACTER(LEN_PRN) :: freqChannel(MAXSYS,MAXFREQ) = ''
    !xsy-2024-04-09
    INTEGER(IT) :: nsat = 0

    INTEGER(IT) :: nfixnl = 0

  END TYPE

END MODULE
