!*
SUBROUTINE get_oi_args(CKF,SAT)
!!
!*
USE const
USE orbit
USE tables
USE station
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
TYPE(ORBCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)

 !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i = 0, j = 0
  INTEGER(IT) :: k = 0, l = 0
  INTEGER(IT) :: narg = 0
  INTEGER(IT) :: lfn = 0
  INTEGER(IT) :: lfnobj = 0
  INTEGER(IT) :: ierr = 0
  INTEGER(IT) :: norder_adams = 0
  REAL(RL) :: int_step = 0.d0
  REAL(RL) :: out_step = 0.d0

  CHARACTER(LEN_FILENAME) :: cfg
  CHARACTER(LEN_STRING) :: msg = ''
  CHARACTER(LEN_STRING) :: key = ''
  CHARACTER(LEN_STRING) :: fmt = ''
  CHARACTER(LEN_STRING) :: bracket = ''

  INTEGER(IT) :: nprn
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  INTEGER(IT) :: nword = 0
  CHARACTER(LEN_ORBPAR) :: word(30) = ''

  REAL(RL) :: isec
  INTEGER(IT) :: iy,im,id,ih,imi,seslen

  LOGICAL(LG) :: lexist,lfound
  TYPE(T_FILETABLE) :: FT
  TYPE(SITE) :: SIT

  !*
  ! The function called
  !!----------------------
  REAL(RL) :: timdif
  CHARACTER(LEN_STRING) :: findkey
  CHARACTER(LEN_STRING) :: lower_string
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: modified_julday
  CHARACTER :: cprn2sys

  !*
  ! Start the exectuable code
  !!---------------------------

  !! Read arguements
  narg=iargc()
  IF (narg .NE. 2) THEN
    WRITE(OUTPUT_UNIT,'(A)') ' Position And Navigation Data Analyst (@PANDA) software'
    WRITE(OUTPUT_UNIT,'(A)') TRIM(VERSION)
    WRITE(OUTPUT_UNIT,'(A)') ' Orbit Integration (OI)'
    WRITE(OUTPUT_UNIT,'(A)') ''
    WRITE(OUTPUT_UNIT,'(A)') ' USAGE: oi file_table flnics'
    CALL exit(1)
  END IF

  CALL getarg(1,cfg)

  INQUIRE(FILE=cfg,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): '//TRIM(cfg)//' is not exist'
    CALL exit(1)
  END IF
  CALL read_filetable(cfg,FT)

  CALL getarg(2,CKF.flnics)
  INQUIRE(FILE=CKF.flnics,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): '//TRIM(CKF.flnics)//' is not exist'
    CALL exit(1)
  END IF
  CKF.flnorb='orb'//TRIM(CKF.flnics(4:))

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=CKF.flnics,STATUS='OLD',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): open '//TRIM(CKF.flnics)
    CALL exit(1)
  END IF

  ! scan the ics file to get the satellites
  nprn=0
  READ(lfn,'(A/A)') msg,msg
  i=INDEX(msg,':')
  READ(msg(i+1:),*) CKF.rmjd,CKF.rsod
  DO WHILE(.TRUE.)
    READ(lfn, '(A)') msg
    IF (INDEX(msg, 'END of FILE') .NE. 0) EXIT
    nprn=nprn+1
    READ(msg(11:), *) cprn(nprn)
    DO WHILE (INDEX(msg,'END of SAT') .EQ. 0)
      READ(lfn,'(A)') msg
    END DO
  END DO
  CLOSE(lfn)
  IF (nprn .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): there are no satellites in '//TRIM(CKF.flnics)
    CALL exit(1)
  END IF

  cfg=f_tablefilename('config')
  INQUIRE(FILE=cfg,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): '//TRIM(cfg)//' is not exist'
    CALL exit(1)
  END IF

  lfn=get_valid_unit(10)
  OPEN(UNIT=lfn,FILE=cfg,STATUS='OLD',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): open '//TRIM(cfg)
    CALL exit(1)
  END IF

  msg = 'Start time&session length'
  key = findkey(lfn,msg,'')
  IF (key(1:5) .EQ. 'EMPTY') GOTO 100
  READ (key,*,err=200) iy,im,id,ih,imi,isec,seslen
  CALL yr2year(iy)
  CKF.mjd0 = modified_julday(id,im,iy)
  CKF.sod0 = ih*3600.d0+imi*60.d0+isec
  CKF.mjd1 = INT(CKF.mjd0+(CKF.sod0+seslen)/86400.d0)
  CKF.sod1 = CKF.sod0+seslen-(CKF.mjd1-CKF.mjd0)*86400.d0

  msg = 'Gravity field modeling'
  key = findkey(lfn,msg,'')
  IF (key(1:5) .EQ. 'EMPTY') GOTO 100
  IF (key(1:2).NE.'NO' .AND. key(1:2).NE.'no') THEN
    CKF.lgfm = .TRUE.
    READ (key,*,err=200) CKF.mindeg, CKF.maxdeg
    CKF.ngc=(CKF.maxdeg+1)**2-CKF.mindeg**2
    CKF.ltog=0
    CALL shcoef_order_wise(CKF.maxdeg,CKF.mindeg,CKF.ltog)
  ELSE
    CKF.lgfm = .FALSE.
    CKF.ngc = 0
    CKF.ltog=0
  END IF

  cfg=f_tablefilename('object')
  INQUIRE(FILE=cfg,EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): '//TRIM(cfg)//' is not exist'
    CALL exit(1)
  END IF

  lfnobj=get_valid_unit(10)
  OPEN(UNIT=lfnobj,FILE=cfg,STATUS='OLD',IOSTAT=ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): open '//TRIM(cfg)
    CALL exit(1)
  END IF

  msg='+Satellite used GNSS'
  key=findkey(lfnobj,msg,'')
  IF(key(1:5) .EQ. 'EMPTY') GOTO 100
  fmt=TRIM(key)

  i=0
  DO WHILE(INDEX(key,'-Satellite used GNSS') .NE. 1)
    READ(lfnobj,'(A)',END=100) key
    IF (key(1:1) .NE. ' ') CYCLE

    i=i+1
    IF (i .GT. MAXSAT) THEN
      WRITE(OUTPUT_UNIT,'(A,I5)') '%%%MESSAGE(get_oi_args): Exceeding the MAXSAT ', MAXSAT
      EXIT
    END IF

    READ(key,fmt,ERR=200) SAT(i).cprn,SAT(i).type,SAT(i).pcv,SAT(i).clk,SAT(i).dclk0, &
       SAT(i).qclk, (SAT(i).dx0(j),j=1,MIN(21,MAXICS))
    SAT(i).dx0(1:6)=SAT(i).dx0(1:6)*1.d-3

    IF (pointer_string(nprn,cprn,SAT(i).cprn) .EQ. 0) THEN
      i=i-1
    ELSE
      CKF.cprn(i) = SAT(i).cprn
      CALL read_svnav(CKF.mjd0+CKF.sod0/86400.d0,SAT(i))
    END IF

  END DO

  msg='+Station used'
  key=findkey(lfnobj,msg,'')
  IF (key(1:5) .EQ. 'EMPTY') GOTO 100
  fmt=TRIM(key)
  DO WHILE(INDEX(key,'-Station used') .NE. 1)
    READ(lfnobj,'(A)',END=100) key
    IF(key(1:1) .NE.' ') CYCLE

    READ(key,fmt,ERR=200) SIT.name, SIT.skd, SIT.pcv, SIT.clk, &
      (SIT.dclk0(j), SIT.qclk(j),j=1,MAXSYS), &
      SIT.cutoff, SIT.map, SIT.dztd0, SIT.qztd, SIT.dgrd0, &
      SIT.qgrd, SIT.dion0, SIT.qion, (SIT.sigr(j),SIT.sigp(j),j=1,MAXSYS),SIT.cprn

    IF (SIT.skd(1:1) .EQ. 'D') THEN
      i=i+1
      IF (pointer_string(nprn,cprn,SIT.cprn) .EQ. 0) THEN
        i=i-1
      ELSE
        SAT(i).cprn=SIT.cprn
        CKF.cprn(i)=SIT.cprn
        CALL read_svnav(CKF.mjd0+CKF.sod0/86400.d0,SAT(i))

        msg=lower_string(SAT(i).type)
        CALL get_file_name(.FALSE.,'att','SATNAM='//TRIM(msg),iy,im,id,ih,SAT(i).flnatt)
        INQUIRE(FILE=SAT(i).flnatt,EXIST=lexist)
        IF (lexist .EQ. .FALSE.) THEN
          WRITE(OUTPUT_UNIT,'(A)') '###MESSAGE(get_oi_args): '//TRIM(SAT(i).flnatt)//' is not exist'
        END IF

        CALL get_file_name(.FALSE.,'acc','SATNAM='//TRIM(msg),iy,im,id,ih,SAT(i).flnacc)
        INQUIRE(FILE=SAT(i).flnacc,EXIST=lexist)
        IF (lexist .EQ. .FALSE.) THEN
          WRITE(OUTPUT_UNIT,'(A)') '###MESSAGE(get_oi_args): '//TRIM(SAT(i).flnacc)//' is not exist'
        END IF
      END IF
    END IF

  END DO
  CKF.nprn=i
  CLOSE(lfnobj)

  CKF.int_step=60.0d0
  CKF.out_step=300.d0
  CKF.norder_adams=11
  DO i=1, CKF.nprn

    SELECT CASE(CKF.cprn(i)(1:1))
      CASE('G')
        bracket='Force model GPS'
      CASE('R')
        bracket='Force model GLONASS'
      CASE('E')
        bracket='Force model GALILEO'
      CASE('C')
        bracket='Force model BEIDOU'
      CASE('S')
        bracket='Force model SBAS'
      CASE('J')
        bracket='Force model QZSS'
      CASE('I')
        bracket='Force model IRNSS'
      CASE('L')
        bracket='Force model '//TRIM(SAT(i).type)
      CASE DEFAULT
        WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): unknown GNSS system '//cprn2sys(CKF.cprn(i))
        CALL exit(1)
    END SELECT

    msg='Integration step etc.'
    key=findkey(lfn,msg,bracket)
    IF (key(1:5) .NE. 'EMPTY' ) THEN
      READ(key, *, ERR=200) int_step,out_step,norder_adams
      IF (int_step .LT. CKF.int_step) CKF.int_step=int_step
      IF (out_step .LT. CKF.out_step) CKF.out_step=out_step
      IF (norder_adams .GT. CKF.norder_adams) CKF.norder_adams=norder_adams
    END IF

    k=0
    DO l=1, MAXFORCE
      j=LEN_TRIM(FORCE_ALL(l))
      IF (j .EQ. 0) CYCLE
      msg=findkey(lfn,FORCE_ALL(l)(1:j),bracket)
      IF (msg(1:5).NE.'EMPTY' .AND. msg(1:4).NE.'NONE') THEN
        k=k+1
        IF (k .GT. MAXFORCE) THEN
          WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): too many forces'
          CALL exit(1)
        END IF
        SAT(i).force_name(k) =TRIM(FORCE_ALL(l))
        SAT(i).force_model(k)=TRIM(msg)
      END IF
    END DO
    IF (k .EQ. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): there are no active force models'
      CALL exit(1)
    ENDIF
    SAT(i).nforce=k

    msg='Variation'
    key=findkey(lfn,msg,bracket)
    IF (key(1:2) .NE. 'NO') THEN
      CKF.lpart=.TRUE.
    END IF

    SAT(i).npar=6
    SAT(i).pname(1)='PXSAT'
    SAT(i).pname(2)='PYSAT'
    SAT(i).pname(3)='PZSAT'
    SAT(i).pname(4)='VXSAT'
    SAT(i).pname(5)='VYSAT'
    SAT(i).pname(6)='VZSAT'
    DO l=1, SAT(i).nforce
      msg=TRIM(SAT(i).force_model(l))
      CALL split_string(.TRUE., msg,'[', ']',' ',nword,word)
      DO j=1, nword
        IF (INDEX(SAT(i).force_name(l),'Customer').NE.0 .AND. INDEX(word(j),'_').EQ.0) THEN
          word(j) = 'EMP_'//word(j)
        END IF
        lfound=.FALSE.
        DO k=6, SAT(i).npar
          IF (SAT(i).pname(k) .EQ. word(j)) THEN
            lfound=.TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. lfound) THEN
          SAT(i).npar=SAT(i).npar+1
          IF (SAT(i).npar .GT. MAXICS) THEN
            WRITE(ERROR_UNIT,'(A)') ' ***ERROR(get_oi_args): there are too many parameters for satellite '//SAT(i).cprn
            CALL exit(1)
          END IF
          SAT(i).pname(SAT(i).npar)=word(j)
        END IF
      END DO
    END DO

   ! Check it with parameter names
    DO l=7, SAT(i).npar
      lfound=.FALSE.
      DO j=1, MAXORBPAR
        IF (SAT(i).pname(l) .EQ. PARAM_ALL(j)) THEN
          lfound=.TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. lfound) THEN
        WRITE(ERROR_UNIT,'(A)') ' ***ERROR(get_oi_args): unknown force model parameter '//TRIM(SAT(i).pname(l))
        CALL exit(1)
      END IF
    END DO

    IF (CKF.lgfm .EQ. .TRUE.) THEN
      IF (SAT(i).npar+CKF.ngc .GT. MAXICS) THEN
        WRITE(ERROR_UNIT,'(A)') ' ***ERROR(get_oi_args): there are too many parameters for satellite '//SAT(i).cprn
        CALL exit(1)
      END IF
      DO j=CKF.mindeg, CKF.maxdeg
        DO l=0, j
          WRITE(SAT(i).pname(SAT(i).npar+CKF.ltog(j,l,1)),"(A3,2I3.3)") 'CNM',j,l
          IF (l .NE. 0) THEN
            WRITE(SAT(i).pname(SAT(i).npar+CKF.ltog(j,l,2)),"(A3,2I3.3)") 'SNM',j,l
          END IF
        END DO
      END DO
      SAT(i).npar=SAT(i).npar+CKF.ngc
    END IF

    k=6*(SAT(i).npar-6+1+6)
    IF (k .GT. CKF.nequ) CKF.nequ=k

  END DO

  CKF.dintv=CKF.out_step
  CKF.sod0 =NINT(CKF.sod0/CKF.dintv)*CKF.dintv
  CKF.sod1 =NINT(CKF.sod1/CKF.dintv)*CKF.dintv
  CALL timinc(CKF.mjd0,CKF.sod0,-8.d0*CKF.dintv,CKF.mjd0,CKF.sod0)
  CALL timinc(CKF.mjd1,CKF.sod1, 8.d0*CKF.dintv,CKF.mjd1,CKF.sod1)

  IF (CKF.lpart .EQ. .FALSE.) CKF.nequ=6

  RETURN


100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): find option '//TRIM(msg)//TRIM(key)
  CALL exit(1)

200 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(get_oi_args): read option '//TRIM(msg)//TRIM(key)
  CALL exit(1)

END SUBROUTINE
