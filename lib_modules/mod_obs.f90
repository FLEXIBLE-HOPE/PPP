!*
MODULE observation
!*
USE par
USE const
IMPLICIT NONE

  !*
  ! Rinex header information
  !!------------------------
  TYPE rnxhead
    REAL(RL)  :: ver
    CHARACTER(4) :: sys
    CHARACTER(4) :: mark
    CHARACTER(LEN_ANTENNA) :: rectype,anttype,recvers
    CHARACTER(LEN_ANTENNA) :: recnum,antnum
    REAL(RL) :: x,y,z,h,e,n
    INTEGER(IT) :: fact1,fact2,nobstype(MAXSYS)
    CHARACTER(LEN_OBSTYPE) :: obstype(MAXOBSTYP,MAXSYS)
    REAL(RL) :: intv
    INTEGER(IT) :: nprn
    CHARACTER(LEN_PRN) :: cprn(MAXSAT)
    INTEGER(IT) :: T0(6),T1(6)
    ! the number of GNSS systems in the observation file
    INTEGER(IT) :: nsys
    CHARACTER(8) :: tsys                           !for rinex 3.0 time system
    CHARACTER :: usetype(MAXFREQ,MAXSYS)=''
  END TYPE rnxhead

  TYPE rnxobs

    ! be mind for the HISI device, used for P1-P2 thred of read_rnxobs.f90
    LOGICAL(LG) :: lhisi = .FALSE.

    INTEGER(IT) :: nprn = 0
    CHARACTER(LEN_PRN) :: cprn(MAXSAT) = ''

    !NOLS DETECTER,    0:NLOS,  1:LOS
    INTEGER(IT) :: nlosflag(MAXSAT) = 0
    ! xsy-2024-02-19: LLI
    INTEGER(IT) :: lli(MAXSAT,MAXFREQ) = 0

    INTEGER(IT) :: ntyp = 0
    CHARACTER(LEN_OBSTYPE) :: obstyp(MAXOBSTYP,MAXSYS) = ''

    INTEGER(IT) :: jd = 0
    REAL(RL) :: tsec = 0.d0
    REAL(RL) :: dtrcv = 0.d0

    REAL(RL) :: obs(MAXSAT,4*MAXFREQ) = 0.d0                !xsy: Phase Code Doppler SNR
    CHARACTER(LEN_OBSTYPE) :: fob(MAXSAT,4*MAXFREQ) = ''
   
    !!xsy-22-9-18: 存储观测文件各个卫星系统采用的频率和通道,用来匹配OSB
    CHARACTER(LEN_OBSTYPE) ::  code_type(MAXFREQ,MAXSYS) =''
    CHARACTER(LEN_OBSTYPE) :: phase_type(MAXFREQ,MAXSYS) = ''

    REAL(RL) :: lifamb(MAXSAT,MAXFREQ,2) = 0.d0
    REAL(RL) :: spndel(MAXSAT,MAXFREQ,2) = 0.d0
    INTEGER(IT) :: flag(MAXSAT,MAXFREQ) = 0

    INTEGER(IT) :: npar = 0
    INTEGER(IT) :: ltog(MAXPARSIT,MAXSAT) = 0

    REAL(RL) :: azim(MAXSAT) = 0.d0
    REAL(RL) :: elev(MAXSAT) = 0.d0
    REAL(RL) :: mltperr(MAXSAT,MAXFREQ) = 0.d0
    REAL(RL) :: nadir(MAXSAT) = 0.d0
    REAL(RL) :: azum(MAXSAT) = 0.d0
    REAL(RL) :: beta(MAXSAT) = 0.d0
    REAL(RL) :: mu(MAXSAT) = 0.d0
    REAL(RL) :: delay(MAXSAT) = 0.d0
    REAL(RL) :: omc(MAXSAT,2*MAXFREQ) = 0.d0
    REAL(RL) :: var(MAXSAT,2*MAXFREQ) = 0.d0
    REAL(RL) :: amat(MAXPARSIT,MAXSAT) = 0.d0
    CHARACTER(LEN_PARNAME) :: pname(MAXPARSIT) = ''

    REAL(RL) :: zmap(MAXSAT) = 0.d0
    REAL(RL) :: clkjmp = 0.d0

  END TYPE

  TYPE slrobs
    CHARACTER(LEN_PRN) :: cprn = ''
    CHARACTER(LEN_SITENAME) :: snam = ''

    INTEGER(IT) :: mjd = 0
    REAL(RL) :: sod = 0.d0
    REAL(RL) :: delay = 0.d0

    REAL(RL) :: p0 = 0.d0
    REAL(RL) :: t0 = 0.d0
    REAL(RL) :: hr0 = 0.d0
    REAL(RL) :: wave_length = 0.d0
    REAL(RL) :: sys_delay = 0.d0

    REAL(RL) :: omc=0.d0
    REAL(RL) :: azim=0.d0
    REAL(RL) :: elev=0.d0
    REAL(RL) :: var = 0.d0
    REAL(RL) :: amat(MAXPARSIT) = 0.d0
    CHARACTER(LEN_PARNAME) :: pname(MAXPARSIT) = ''

    INTEGER(IT) :: npar = 0
    INTEGER(IT) :: ltog(MAXPARSIT,MAXSAT) = 0

    LOGICAL(LG) :: ltro = .FALSE.
    LOGICAL(LG) :: lcom = .FALSE.
    LOGICAL(LG) :: lsitd = .FALSE.
    LOGICAL(LG) :: lsatd = .FALSE.
  END TYPE

  TYPE islobs
    CHARACTER(LEN_PRN) :: cprn(2)
    INTEGER(IT) :: mjd(2) = 0
    REAL(RL) :: sod(2) = 0.d0
    REAL(RL) :: isl(2) = 0.d0

    INTEGER(IT) :: flag = 0

    INTEGER(IT) :: nus
    REAL(RL) :: peod = 0.d0
    CHARACTER(LEN_PARNAME) :: uname(10) = ''

    REAL(RL) :: omc = 0.d0
    REAL(RL) :: sig = 0.d0
    REAL(RL) :: var = 0.d0
    REAL(RL) :: azim=0.d0
    REAL(RL) :: elev=0.d0
    REAL(RL) :: nadir = 0.d0
    REAL(RL) :: azum = 0.d0

    REAL(RL) :: amat(MAXPARSIT,MAXSAT) = 0.d0
    CHARACTER(LEN_PARNAME) :: pname(MAXPARSIT) = ''

    INTEGER(IT) :: npar = 0
    INTEGER(IT) :: ltog(MAXPARSIT,MAXSAT) = 0
  END TYPE

END MODULE
