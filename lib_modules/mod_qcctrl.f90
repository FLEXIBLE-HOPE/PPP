!*
MODULE qcctrl
!*
USE par
IMPLICIT NONE

   TYPE qccfg

    ! the begining and end time
    INTEGER(IT) :: mjd = 0
    INTEGER(IT) :: mjd0 = 0
    INTEGER(IT) :: mjd1 = 0
    REAL(RL) :: sod = 0.d0
    REAL(RL) :: sod0 = 0.d0
    REAL(RL) :: sod1 = 0.d0
    REAL(RL) :: dintv = 0.d0
    REAL(RL) :: pintv = 0.d0

    ! stations and satellites
    INTEGER(IT) :: nsit = 0
    INTEGER(IT) :: nprn = 0
    INTEGER(IT) :: nleo = 0
    CHARACTER(LEN_PRN) :: cprn(MAXSAT) = ''

    ! Observation used
    CHARACTER(LEN=20) :: uobs = 'CODE PHASE'
    CHARACTER(LEN=20) :: cobs = 'IF'
    CHARACTER(LEN=20) :: dobs = 'UD'
    CHARACTER(LEN_ZTDMOD) :: ztdmod = ''
    INTEGER(IT) :: nfq(MAXSYS) = 0

    ! statistics information
    !! For edtres
    INTEGER(IT) :: nobs(MAXSYS+1) = 0
    INTEGER(IT) :: pobs(MAXSYS+1) = 0
    INTEGER(IT) :: nele(MAXSYS+1) = 0

    ! The GNSS system
    INTEGER(IT) :: iref = 0
    INTEGER(IT) :: nsys = 0
    CHARACTER(MAXSYS) :: system = ''
    INTEGER(IT) :: nfreq(MAXSYS) = 0
    CHARACTER(LEN_FREQ) :: freq(MAXFREQ,MAXSYS) = ''

    CHARACTER(LEN_FILENAME) :: flnorb = ''
    CHARACTER(LEN_FILENAME) :: flnleo = ''
    CHARACTER(LEN_FILENAME) :: flncck = ''
    CHARACTER(LEN_FILENAME) :: flnrck = ''
    CHARACTER(LEN_FILENAME) :: flnztd = ''
    CHARACTER(LEN_FILENAME) :: flnkin = ''
    CHARACTER(LEN_FILENAME) :: flnion = ''
  
    !! For edtres
    INTEGER(IT) :: lfnsum
    INTEGER(IT) :: lfnres

    ! The control parameters
    LOGICAL(LG) :: lrenew = .FALSE.
    LOGICAL(LG) :: llc_check = .FALSE.
    LOGICAL(LG) :: lpc_check = .FALSE.
    LOGICAL(LG) :: lmw_check = .FALSE.
    LOGICAL(LG) :: llg_check = .FALSE.
    LOGICAL(LG) :: lorb_use = .FALSE.
    LOGICAL(LG) :: lrinex_flag_use = .FALSE.
    INTEGER(IT) :: len_gap = 0
    INTEGER(IT) :: len_short = 0
    INTEGER(IT) :: nsat = 0
    REAL(RL) :: cutoff = 0.d0
    REAL(RL) :: lc_limit = 0.d0
    REAL(RL) :: pc_limit = 0.d0

    REAL(RL) :: lw_limit = 0.d0
    REAL(RL) :: lg_limit = 0.d0
    REAL(RL) :: lg_rms_limit = 0.d0
    INTEGER(IT) :: nsec_per_degree_lg = 0
    INTEGER(IT) :: niter_lg = 0

    REAL(RL) :: max_mean_namb = 0.d0
    REAL(RL) :: min_percent = 0.d0
    INTEGER(IT) :: min_mean_nprn = 0

    !! For edtres
    LOGICAL(LG) :: lcode = .FALSE.
    REAL(RL) :: bad_limit = 0.d0
    REAL(RL) :: jump_limit = 0.d0

    !! For extclk
    LOGICAL(LG) :: lsat = .FALSE.
    LOGICAL(LG) :: lrec = .FALSE.

  END TYPE qccfg

END MODULE qcctrl
