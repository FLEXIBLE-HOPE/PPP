!*
MODULE ckdctrl
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  TYPE ckdcfg
    ! socket number and port
    INTEGER(IT) :: sock = 0
    CHARACTER(LEN_PORT) :: port = ''
    LOGICAL(LG) :: ldebug = .FALSE.
    LOGICAL(LG) :: lpost = .FALSE.
    LOGICAL(LG) :: llog = .FALSE.
    ! lsq orbit determination
    LOGICAL(LG) :: lpod = .FALSE.
    ! lsq fixing ambiguity
    LOGICAL(LG) :: lamb = .FALSE.
    ! using inter-satellite link data
    LOGICAL(LG) :: lisl = .FALSE.
    ! tread BDS-2 and BDS-3 as indepedent system
    LOGICAL(LG) :: lbds = .FALSE.
    ! estimate receiver clock(FALSE) or ISB(TRUE)
    LOGICAL(LG) :: lisb = .TRUE.
    ! estimate the RECDCB for the third frequency code observations(FALSE)
    LOGICAL(LG) :: lrecdcb = .FALSE.    
    ! estimate LEO time synchronization bias (FALSE) or IFB(TRUE)
    LOGICAL(LG) :: lleoifb = .TRUE.
    ! correct the consistency when processing HISI data
    LOGICAL(LG) :: lhisi = .FALSE.

    ! Observation used
    CHARACTER(LEN=20) :: uobs = 'CODE PHASE'
    CHARACTER(LEN=20) :: cobs = 'IF'
    CHARACTER(LEN=20) :: dobs = 'UD'
    INTEGER(IT) :: nfq(MAXSYS) = 0

    ! The GNSS system
    INTEGER(IT) :: iref = 0
    INTEGER(IT) :: nsys = 0
    CHARACTER(MAXSYS) :: system = ''
    INTEGER(IT) :: nfreq(MAXSYS) = 0                      !每个系统频率数,支持混频
    CHARACTER(LEN_FREQ) :: freq(MAXFREQ,MAXSYS) = ''
    
    INTEGER(IT) :: refnprn(MAXSYS) = 0                    !各系统参考星
    CHARACTER(LEN_PRN) :: refcprn(MAXSYS) = ''
    INTEGER(IT) :: refsat = 0                             !不分系统时的参考星
    CHARACTER(LEN_PRN) :: refprn = ''

    !@CMT BY XSY: reference satellite and station for sdb program
    INTEGER(IT) :: ind_rsat = 0
    INTEGER(IT) :: ind_rsit = 0
    CHARACTER(LEN_PRN) :: rsat = ''
    CHARACTER(LEN_SITENAME) :: rsit = ''

    CHARACTER(LEN=10) :: ArBiasMode = ''
    !@CMT BY XSY: DCB mode: osb or cc2nocc
    CHARACTER(LEN=10) :: codebias = ''    
    !@CMT BY XSY: OSB mode
    CHARACTER(LEN=10) :: osbMode = ''

    ! Time of epoch
    INTEGER(IT) :: mjd = 0
    INTEGER(IT) :: totepc = 0
    REAL(RL) :: sod = 0.D0
    REAL(RL) :: dintv = 0.D0

    ! the begining and end time
    INTEGER(IT) :: mjd0 = 0
    INTEGER(IT) :: mjd1 = 0
    REAL(RL) :: sod0 = 0.d0
    REAL(RL) :: sod1 = 0.d0

    ! stations and satellites
    INTEGER(IT) :: nsit = 0
    INTEGER(IT) :: nprn = 0
    INTEGER(IT) :: nleo = 0
    CHARACTER(LEN_PRN) :: cprn(MAXSAT) = ''

    ! carrier-phase smooth the code
    ! NONE/NO
    ! PHASE
    ! DOPPLER
    CHARACTER(LEN_ZTDMOD) :: pscmod = ''
    INTEGER(IT) :: pscwnd = 0
    ! troposphere dealy estimation approach
    ! epoch wise or piece-wise constant
    CHARACTER(LEN_ZTDMOD) :: ztdmod = ''
    ! troposphere gradient dealy estimation approach
    ! epoch wise (EW) or piece-wise constant (PWC)
    CHARACTER(LEN_ZTDMOD) :: grdmod = ''
    ! ionosphere delay estimation approach
    ! option:
    ! 1) SID: slant ionosphere delay, EW/PWC
    ! 2) ZTD: zenith ionosphere delay, EW/PWC
    ! 3) KLO: correction by klobuchar model
    ! 4) NEQ: correction by NeQuick model
    ! 5) GIM: correction by GIM model
    CHARACTER(LEN_IONMOD) :: ionmod = ''
    !@CMT BY XSY: PPPRTK
    REAL(RL) :: ionConstrians = 10.d0
    REAL(RL) :: ztdConstrians = 10.d0
    !@CMT BY XSY: PPPRTK ionsphere product model
    ! 1) PDUD: undifference products from PANDA intion
    ! 2) PDSD: satellite difference products from PANDA intion
    ! 3) SD  : satellite difference products from huixia chen
    ! 4) GRID: grid products from huixia chen
    CHARACTER(LEN_IONMOD) :: ionpro = 'PDUD'
    !@CMT BY XSY: SIT constrians
    LOGICAL(LG) :: lsitcons = .FALSE.
    REAL(RL) :: sitConstrians(3) = 10.d0
    !@CMT BY XSY: set the lockout span
    REAL(RL) :: lockout(2) = 0.d0

    ! data pre-processing
    INTEGER(IT) :: gap = 0
    REAL(RL) :: lg = 0.D0
    REAL(RL) :: lw = 0.D0
    !@CMT BY XSY: LLI for cycle slip detection
    LOGICAL(LG) :: lli_flag = .FALSE.

    ! Gravity field modeling
    LOGICAL(IT) :: lgfm = .FALSE.
    INTEGER(IT) :: mindeg = 0
    INTEGER(IT) :: maxdeg = 0
    INTEGER(IT) :: ngc = 0
    INTEGER(IT) :: ltog(0:MAXGRADEG,0:MAXGRADEG,2) = 0

    !@CMT BY XSY: NLOS factor
    LOGICAL(LG) :: nlos = .FALSE. !true:直接删除nlos观测值
    INTEGER(IT) :: nlosfac = 1
    !@CMT BY XSY: SNR detele
    INTEGER(IT) :: SnrLimit = 32
    !@CMT BY XSY: IGG3 decision
    LOGICAL(LG) :: Igg3 = .FALSE. 
    REAL(RL) :: Igg3_k0 = 0.d0
    REAL(RL) :: Igg3_k1 = 0.d0
    !@CMT BY XSY: repair cycle slip or not
    LOGICAL(LG) :: repslip = .FALSE.
    !@CMT BY XSY: DIA detection thershold
    REAL(RL) :: DiaSig = 0.d0
    !@CMT BY XSY: Reconvergence time
    INTEGER(IT) :: ReConvTime = 0
    !@CMT BY XSY: AR mode
    CHARACTER(LEN=20) :: ArMode = ''
    !@CMT BY XSY: Whether map FCB to OSB
    LOGICAL(LG) :: FcbMapToOsb = .FALSE.

    CHARACTER(LEN_STRING) :: freq_used = ''
    CHARACTER(LEN_STRING) :: signal_used = ''
   
    ! Files
    CHARACTER(LEN_FILENAME) :: flnorb = ''
    CHARACTER(LEN_FILENAME) :: flnleo = ''
    CHARACTER(LEN_FILENAME) :: flnerp = ''
    CHARACTER(LEN_FILENAME) :: flnrck = ''
    CHARACTER(LEN_FILENAME) :: flnztd = ''
    CHARACTER(LEN_FILENAME) :: flndop = ''
    CHARACTER(LEN_FILENAME) :: flnres = ''
    CHARACTER(LEN_FILENAME) :: flncck = ''
    CHARACTER(LEN_FILENAME) :: flnamb = ''
    CHARACTER(LEN_FILENAME) :: flnsnx = ''
    CHARACTER(LEN_FILENAME) :: flnfnl = ''
    CHARACTER(LEN_FILENAME) :: flnfwl = ''
    CHARACTER(LEN_FILENAME) :: flnfewl = ''
    CHARACTER(LEN_FILENAME) :: flnfeewl = ''
    CHARACTER(LEN_FILENAME) :: flnfhewl = ''
    CHARACTER(LEN_FILENAME) :: flnpos = ''
    CHARACTER(LEN_FILENAME) :: flnkin = ''
    CHARACTER(LEN_FILENAME) :: flnslr = ''
    CHARACTER(LEN_FILENAME) :: flnisl = ''
    CHARACTER(LEN_FILENAME) :: flnion = ''
    CHARACTER(LEN_FILENAME) :: flnbrd = ''
    CHARACTER(LEN_FILENAME) :: flnifcb = ''
    CHARACTER(LEN_FILENAME) :: flnmw = ''
    ! For lsq
    CHARACTER(LEN_FILENAME) :: flnpar = ''
    CHARACTER(LEN_FILENAME) :: flnnom = ''
    !@CMT BY XSY:
    CHARACTER(LEN_FILENAME) :: flnosb = ''
    CHARACTER(LEN_FILENAME) :: flnosbcode = ''  !CodeOsb file for FCB as input file
    CHARACTER(LEN_FILENAME) :: flnosbupd = ''   !OSB file for FCB by my own FCB, Phase OSB and Code osb
    CHARACTER(LEN_FILENAME) :: flnpco = ''
    CHARACTER(LEN_FILENAME) :: flnrbias = ''    !FPA bias file between different receiver types, including NL/WL/EWL/EEWL/HEWL
    CHARACTER(LEN_FILENAME) :: flnsdb = '' 
    CHARACTER(LEN_FILENAME) :: flnomc = ''
    CHARACTER(LEN_FILENAME) :: flngfif = ''
    CHARACTER(LEN_FILENAME) :: flndcb = ''
    CHARACTER(LEN_FILENAME) :: flngf = ''   

    ! Temporal files for smooth
    INTEGER(IT) :: lfnobs = ERROR_UNIT
    INTEGER(IT) :: lfncid = ERROR_UNIT
    INTEGER(IT) :: lfnrem = ERROR_UNIT


    ! The options for GLONASS ISB
    LOGICAL(LG) :: lrisb = .FALSE.
    INTEGER(IT) :: imw = -1
    ! The configration options for ambiguity
    ! Integer Ambiguity resolution or not
    LOGICAL(LG) :: liar = .TRUE.
    ! The minmam common period for satellite paries in second
    REAL(RL) :: minsec_common = 1200.D0
    !@CMT BY XSY: whether to correct the MW PCO for WL ambiguity fixing
    LOGICAL(LG) :: if_pco_corr = .FALSE.
    ! The decision parameters for wide-lane ambiguity resolution in cycle
    REAL(RL) :: wl_maxdev = 0.1D0
    REAL(RL) :: wl_maxsig = 0.1D0
    REAL(RL) :: wl_alpha = 1.D3
    REAL(RL) :: wl_dintv = 30.d0
    ! The decision parameters for narrow-lane ambiguity resolution in cycle
    ! For ambfix
    REAL(RL) :: nl_maxdev = 0.1D0
    REAL(RL) :: nl_maxsig = 0.1D0
    REAL(RL) :: nl_alpha = 1.D3
    REAL(RL) :: bl_limit = 5000.d0
    ! The cut-off angle of NL ambiguities
    REAL(RL) :: cutoff = 10.d0
    ! The cut-off angle of WL ambiguities
    REAL(RL) :: wlcutoff = 10.d0
    ! The length for updating Melbourne-Wubbean ambiguity
    REAL(RL) :: updwl = 86400.D0
    ! Maximum number of excluded ambiguities
    INTEGER(IT) :: nl_maxdel = 0
    ! Minimum number of reserved ambiguities
    INTEGER(IT) :: nl_minsav = 0
    ! Threshold values for chi- and ratio tests for NL and WL
    REAL(RL) :: nl_ratio = 0
    REAL(RL) :: wl_ratio = 0
    !@CMT BY XSY: used for ppp, no fixing satellite
    CHARACTER(LEN_STRING) :: nofixsat = ''
    !@CMT BY XSY: used for fcb, the minmum common stations
    INTEGER(IT) :: commonsit = 5

    !@CMT BY XSY: used for sdb, the MW minimu arc length for sdb program
    REAL(RL) :: mwarc_gnss2leo = 600.D0
    REAL(RL) :: mwarc_gnss2sit = 600.D0
    REAL(RL) :: mwarc_leo2sit  = 600.D0
    !@CMT BY XSY: used for sdb, the input MW SDB as init
    CHARACTER(LEN_FILENAME) :: flnmw_init = ''   
    !@CMT BY XSY: delete the bad station for -net model of sdb and dcb program
    CHARACTER(LEN_SITENAME) :: delsit(MAXSIT) = ''

    !@CMT BY XSY: used for rnxcmp
    CHARACTER(LEN=20) :: done = ''
    LOGICAL(LG) :: lint = .FALSE.
    LOGICAL(LG) :: lbroad = .FALSE.
    LOGICAL(LG) :: litrs = .FALSE.
    LOGICAL(LG) :: lsp3 = .FALSE.
    CHARACTER(LEN=20) :: brdc(MAXSYS) = ''
    CHARACTER(LEN_FILENAME) :: flntec = ''
  END TYPE

  TYPE EOPMOD
    CHARACTER(LEN_STRING):: model
    REAL(RL) :: x(6) = 0.d0
    REAL(RL) :: dx0(3) = 0.d0
    REAL(RL) :: q(3) = 0.d0
  END TYPE EOPMOD

  !! solution report
  TYPE SOL
    INTEGER(IT) :: ncad = 0         !xsy: 即将固定的窄巷数量
    INTEGER(IT) :: nfix = 0         !xsy: 固定的窄巷数量
    REAL(RL) :: refx(3) = 0.D0      !xsy: 参考坐标，来自gnss.sit
    REAL(RL) :: ratio = 0.D0        !xsy：可靠性检验
    REAL(RL) :: adop = 0.D0
    REAL(RL) :: fpos(3) = 0.D0      !xsy：浮点解
    REAL(RL) :: xpos(3) = 0.D0      !xsy: 固定解
    REAL(RL) :: rot(3,3) = 0.D0     !xsy: rot_l2f(enu2xyz)

    INTEGER(IT) :: fix_ewl = 0      !xsy: 是否成功固定，0：NO   1：YES
    INTEGER(IT) :: fixnum_ewl = 0
    INTEGER(IT) :: fix_wl = 0       !xsy: 是否成功固定，0：NO   1：YES
    INTEGER(IT) :: fixnum_wl = 0
    INTEGER(IT) :: fix_nl = 0       !xsy: 是否成功固定，0：NO   1：YES
    INTEGER(IT) :: fixnum_nl = 0
  END TYPE

  !! receiver information
  TYPE RECINFO
      LOGICAL(LG) :: lstation = .FALSE.
      CHARACTER(LEN_SITENAME) :: name = ''   
      CHARACTER(LEN_STRING) :: rectype = ''
      CHARACTER(LEN_STRING) :: recvers = ''
      CHARACTER(LEN_STRING) :: anttype = ''
  END TYPE

END MODULE
