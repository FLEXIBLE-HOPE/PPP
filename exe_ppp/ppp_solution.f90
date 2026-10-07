!
!! purpose  : output solutions
!! parameter:
!!    input : SCF -- sri configuration
!!            SITE -- station struct
!!            OB -- rinex observations
!!            NM,PM -- parameter & information matrix
!!    output:
!! author   : Geng J
!! created  : Nov. 14, 2007
!
SUBROUTINE ppp_solution(CKF,SIT,SAT,OB,NM,PM,AM,SL,Qxyz)
!*
USE info
USE const
USE ckdctrl
USE station
USE satellite
USE ambiguity
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
TYPE(RNXOBS) :: OB(1:*)
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT(1:*)
TYPE(CKDCFG) :: CKF
TYPE(INFM) :: NM(1:*)
TYPE(PRMT) :: PM(CKF.nsys+1+2+3+2*CKF.nprn+(MAXFREQ-2)*MAXSAT,1:*)
TYPE(AMBT) :: AM(MAXFREQ*MAXSAT,1:*)
TYPE(SOL) :: SL(1:*)
REAL(RL) :: Qxyz(3,3)

  !*
  ! The local variables
  !!--------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,j,isit,isys,ipar,isat,iy,imon,id,ih,im,ifreq,  gpsweek,wd
  INTEGER(IT) :: lfnrck,lfnpos,lfnamb,lfnztd,lfnres,lfngm,nobs
  REAL(RL) :: sec,phase,pseud,rwl,swl,geo(3), sow,rot_l2f_T(3,3),rot_f2l(3,3),det
  CHARACTER :: cid
  CHARACTER(LEN=2) :: ab
  CHARACTER(LEN_STRING) :: line,contype

  INTEGER(IT) :: nlfix(CKF%nprn,CKF%nsit)

  DATA lfirst,lfnrck,lfnpos,lfnamb,lfnztd,lfnres,lfngm /.TRUE.,6*0/
  SAVE lfirst,lfnrck,lfnpos,lfnamb,lfnztd,lfnres,lfngm,nobs

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  CHARACTER(LEN=25) :: run_tim

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    !! receiver clock
    lfnrck=get_valid_unit(10)
    OPEN(UNIT=lfnrck,file=CKF.flnrck)
    !! zenith troposphere delay
    lfnztd=get_valid_unit(10)
    OPEN(UNIT=lfnztd,file=CKF.flnztd)
    !! ambiguity solution
    lfnamb=get_valid_unit(10)
    OPEN(UNIT=lfnamb,file=CKF.flnamb)
    CALL ppp_wt_cfg(CKF,SIT,lfnamb)
    !! position
    lfnpos=get_valid_unit(10)
    OPEN(UNIT=lfnpos,file=CKF.flnpos)
    CALL ppp_wt_cfg(CKF,SIT,lfnpos)
    !! residuals
    lfnres=get_valid_unit(10)
    OPEN(UNIT=lfnres,file=CKF.flnres)
    CALL ppp_wt_cfg(CKF,SIT,lfnres)
    WRITE(lfnres,'(A)')'#SITE PRN Ifreq   RES_phase        RES_code    Weight_phase    Wight_code  qcflag  ELEV  AZIM  NLOSflag  SNR  Sig0'

    IF (CKF.cobs(1:2) .EQ. 'IF') THEN
      nobs=2
    ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
      nobs=MAXVAL(CKF.nfq)*2
    ELSE IF (CKF.cobs(1:7) .EQ. 'GRAPHIC') THEN
      nobs=2
    END IF
  END IF

  !! time tag
  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)

  !! write receiver clocks
  IF (lfnrck .NE. 0) THEN
    DO isys=1, CKF.nsys
      DO isit=1, CKF.nsit
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECCLK'//CKF.system(isys:isys))
        IF (CKF.system(isys:isys).EQ.'R') THEN
          DO isat=1, CKF.nprn
            IF (CKF.cprn(isat)(1:1) .NE. 'R') CYCLE
            DO ipar=1, OB(isit).npar
              IF (INDEX(OB(isit).pname(ipar),'RECCLK'//CKF.cprn(isat)).NE.0) EXIT
            END DO
            IF (ipar .LE. OB(isit).npar) THEN
              ipar=OB(isit).ltog(ipar,isat)
              IF (PM(ipar,isit).iobs .GT. 0) THEN
                WRITE(lfnrck,'(A4,1X,A15,I7,F10.2,F17.6,F14.6)') SIT(isit).name,PM(ipar,isit).pname,CKF.mjd,CKF.sod,&
                                        PM(ipar,isit).xini,PM(ipar,isit).xcor
              END IF
            ELSE
              ipar=0
            END IF
          END DO
        END IF
        
        !@CMT BY XSY: time synchronization bias for LEO satellites
        IF (CKF%lleoifb) THEN
          IF (CKF.system(isys:isys).EQ.'L') THEN
            DO isat=1, CKF.nprn
              IF (CKF.cprn(isat)(1:1) .NE. 'L') CYCLE
              DO ipar=1, OB(isit).npar
                IF (INDEX(OB(isit).pname(ipar),'RECCLK'//CKF.cprn(isat)).NE.0) EXIT
              END DO
              IF (ipar .LE. OB(isit).npar) THEN
                ipar=OB(isit).ltog(ipar,isat)
                IF (PM(ipar,isit).iobs .GT. 0) THEN
                  WRITE(lfnrck,'(A4,1X,A15,I7,F10.2,F17.6,F14.6)') SIT(isit).name,PM(ipar,isit).pname,CKF.mjd,CKF.sod,&
                                          PM(ipar,isit).xini,PM(ipar,isit).xcor
                END IF
              ELSE
                ipar=0
              END IF
            END DO
          END IF
        END IF

        IF (ipar .EQ. 0) CYCLE
        ipar=OB(isit).ltog(ipar,1)
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          WRITE(lfnrck,'(A4,1X,A15,I7,F10.2,F17.6,F14.6)') SIT(isit).name,PM(ipar,isit).pname,CKF.mjd,CKF.sod,&
                                        PM(ipar,isit).xini,PM(ipar,isit).xcor
        END IF

      END DO
    END DO
  END IF

  nlfix=0
  !! write ambiguities
  IF (lfnamb .NE. 0) THEN
    WRITE(lfnamb,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
    DO isit=1, CKF.nsit
      DO ipar=1, NM(isit).ns
        isat =AM(ipar,isit).psat
        ifreq=AM(ipar,isit).ifreq
        IF (ifreq.EQ.1) THEN
          nlfix(isat,isit)=AM(ipar,isit)%ifab
          DO isys=1,CKF%nsys
            IF (CKF%cprn(isat)(1:1).EQ.CKF%system(isys:isys)) THEN
              IF (CKF%refcprn(isys).EQ.CKF%cprn(isat)) THEN
                nlfix(isat,isit)=10
              END IF
            END IF
          END DO
        END IF        
        IF (AM(ipar,isit).iobs .GT. 0) THEN
          IF (CKF.liar) THEN
            !@ CMT BY XSY: FIVE FREQUENCY
            WRITE(lfnamb,'(A4,1X,I1,1X,A3,1X,I1,1X,8F18.6,2F18.10,F9.4,F6.1)') SIT(isit).name,AM(ipar,isit).ifab,CKF.cprn(AM(ipar,isit).psat),&
              ifreq,AM(ipar,isit).xini,AM(ipar,isit).xcor,AM(ipar,isit).xini+AM(ipar,isit).xcor,AM(ipar,isit).abhewl,AM(ipar,isit).abeewl,AM(ipar,isit).abewl,AM(ipar,isit).abwl,AM(ipar,isit).abnl,AM(ipar,isit).ptime(1:2),&
              AM(ipar,isit).xsig,AM(ipar,isit).elev/AM(ipar,isit).iobs  !AM(ipar,isit).xini+AM(ipar,isit).xcor:浮点模糊度   AM(ipar,isit).ifab:固定类型
          ELSE
            IF (CKF.cobs(1:2) .EQ. 'IF') THEN
              rwl=AM(ipar,isit).xrwl/AM(ipar,isit).weig+AM(ipar,isit).abin
              swl=0.d0
              IF (AM(ipar,isit).iobs .GT. 1) THEN
                swl=DSQRT((AM(ipar,isit).xswl-AM(ipar,isit).xrwl**2/AM(ipar,isit).weig)/AM(ipar,isit).weig/(AM(ipar,isit).iobs-1))
              END IF
            ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
              rwl=0.d0
              swl=0.d0
            END IF
            WRITE(lfnamb,'(A4,1X,A1,1X,A3,1X,I1,1X,4F18.6,2F18.10,2F9.4,F6.1)') SIT(isit).name,'0',CKF.cprn(isat),ifreq,AM(ipar,isit).xini,AM(ipar,isit).xcor,AM(ipar,isit).xini+AM(ipar,isit).xcor,rwl,&
                                   AM(ipar,isit).ptime(1:2),AM(ipar,isit).xsig,swl,AM(ipar,isit).elev/AM(ipar,isit).iobs
          END IF
        END IF
      END DO
    END DO
  END IF

  !! write zenith troposphere delay
  IF (lfnztd .NE. 0 .AND. CKF%ztdmod(1:3).NE.'FIX') THEN
    DO isit=1, CKF.nsit
      IF (CKF.ztdmod(1:4).EQ.'NONE' .OR. CKF.ztdmod(1:3).EQ.'FIX' .OR. SIT(isit)%skd(1:1).EQ.'D') THEN
      ELSE
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,CKF.ztdmod)
        ipar=OB(isit).ltog(ipar,1)
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,3F11.6)') 'ZTD   ',SIT(isit).name,CKF.mjd,CKF.sod,SIT(isit).zdd,&
                      SIT(isit).zwd+SIT(isit).ztdcor,PM(ipar,isit).xcor !xcor:湿延迟残余参数的改正数,ztdcor:残余量初值
        END IF

        IF (CKF.grdmod(1:4) .NE. 'NONE') THEN
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'N'//CKF.grdmod)
          ipar=OB(isit).ltog(ipar,1)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,F11.6)') 'NGRD  ', SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest
          END IF
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'S'//CKF.grdmod)
          ipar=OB(isit).ltog(ipar,1)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,F11.6)') 'SGRD  ', SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest
          END IF
        END IF
      END IF

      IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
        DO isat=1, CKF.nprn
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'ION'//CKF.cprn(isat))
          ipar=OB(isit).ltog(ipar,isat)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,F11.6,F6.1,1X,I6)') 'ION'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest,OB(isit)%elev(isat)*RAD2DEG,nlfix(isat,isit)
          END IF
        END DO
      END IF

      IF (CKF%lbds) THEN
        DO isat=1, CKF.nprn
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'SISRE'//CKF.cprn(isat))
          ipar=OB(isit).ltog(ipar,isat)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A8,1X,A4,I7,F10.2,2F11.3,F6.1)') 'SISRE'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xini,PM(ipar,isit).xest,OB(isit)%elev(isat)*RAD2DEG
          END IF
        END DO
      END IF

      ! IF (SIT(isit)%skd(1:2).EQ.'DE') THEN
      !   ipar=pointer_string(OB(isit).npar,OB(isit).pname,'DRAG_c')
      !   ipar=OB(isit).ltog(ipar,1)
      !   IF (PM(ipar,isit).iobs .GT. 0) THEN
      !     WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,2F11.6)') 'DRAG_c',SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xini,PM(ipar,isit).xest
      !   END IF      
      ! END IF
    END DO
   
    !@CMT BY XSY: [RECDCB] The third frequency code observation is a useless contribution because its' low wight
    IF (CKF%lrecdcb) THEN
      IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
        DO isit=1, CKF.nsit
          DO isat=1, CKF.nprn
            isys=INDEX(SYS,CKF.cprn(isat)(1:1))
            IF (CKF%nfq(isys) .GE. 3) THEN
              IF (OB(isit)%elev(isat) .LT. SIT(isit)%cutoff) CYCLE
              
              ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECDCBL3'//CKF.cprn(isat))
              IF (ipar .EQ. 0) CYCLE
              ipar=OB(isit).ltog(ipar,isat)
              IF (PM(ipar,isit).iobs .GT. 0) THEN
                WRITE(lfnztd,'(A11,1X,A4,I7,F10.2,F11.6,F11.1)') 'RECDCBL3'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest,OB(isit)%elev(isat)*RAD2DEG
              END IF

              IF (CKF%nfq(isys) .GE. 4) THEN
                ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECDCBL4'//CKF.cprn(isat))
                IF (ipar .EQ. 0) CYCLE
                ipar=OB(isit).ltog(ipar,isat)
                IF (PM(ipar,isit).iobs .GT. 0) THEN
                  WRITE(lfnztd,'(A11,1X,A4,I7,F10.2,F11.6,F11.1)') 'RECDCBL4'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest,OB(isit)%elev(isat)*RAD2DEG
                END IF

                IF (CKF%nfq(isys) .GE. 5) THEN
                  ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECDCBL5'//CKF.cprn(isat))
                  IF (ipar .EQ. 0) CYCLE
                  ipar=OB(isit).ltog(ipar,isat)
                  IF (PM(ipar,isit).iobs .GT. 0) THEN
                    WRITE(lfnztd,'(A11,1X,A4,I7,F10.2,F11.6,F11.1)') 'RECDCBL5'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest,OB(isit)%elev(isat)*RAD2DEG
                  END IF                 
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END IF

  !! static & kinematic position estimates
  IF (lfnpos .NE. 0) THEN
    DO isit=1, CKF.nsit
      IF (SIT(isit).skd(1:2) .EQ. 'DE') THEN
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'PXSAT')
        IF (ipar .EQ. 0) THEN
          WRITE(OUTPUT_UNIT,'(A,A,I6,F9.2,A2,64X,"  (",I2,F6.1,")")') run_tim(), SIT(isit).name, &
               CKF.mjd,CKF.sod,ab,NM(isit).nobs/2,NM(isit).esig
          CYCLE
        END IF
        ipar=OB(isit).ltog(ipar,1)
        ab=' '
        IF (PM(ipar,isit).iobs .LT. 4) THEN
          ab='* '
        ELSE IF (SIT(isit).nfixnl .GE. 4) THEN
          ab='x '
        END IF

        SL(isit).xpos(1)=PM(ipar,isit).xest-SL(isit).refx(1)
        SL(isit).xpos(2)=PM(ipar+1,isit).xest-SL(isit).refx(2)
        SL(isit).xpos(3)=PM(ipar+2,isit).xest-SL(isit).refx(3)
        SL(isit).xpos(1:3)=SL(isit).xpos(1:3)*1.d5
        CALL matinv(SIT(isit).rot_l2f,3,3,swl)
        CALL matmpy(SIT(isit).rot_l2f,SL(isit).xpos,geo,3,3,1)
        SL(isit).fpos(1)=SL(isit).fpos(1)-SL(isit).refx(1)
        SL(isit).fpos(2)=SL(isit).fpos(2)-SL(isit).refx(2)
        SL(isit).fpos(3)=SL(isit).fpos(3)-SL(isit).refx(3)
        SL(isit).fpos(1:3)=SL(isit).fpos(1:3)*1.d5
        WRITE(OUTPUT_UNIT,'(A,A,I5,4I3,F5.1,A2,3F8.1,"  [",I3,I3,F6.1,"]",3F8.1,"  (",I2,F6.2,")",6I3)') run_tim(), SIT(isit).name, &
              iy,imon,id,ih,im,sec,ab,SL(isit).xpos(1:3),SL(isit).ncad,SL(isit).nfix,SL(isit).ratio,SL(isit).fpos(1:3), &
              SIT(isit)%nsat,NM(isit).esig,SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl, &
              SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl        
        ! WRITE(lfnpos,'(I4,4I3,F5.1,A2,A4,3F13.3,3F9.3,3F13.3,3F23.12,I5)') iy,imon,id,ih,im,sec,ab,SIT(isit).name, &
        !   PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,PM(ipar,isit).xsig,PM(ipar+1,isit).xsig, &
        !   PM(ipar+2,isit).xsig,SL(isit).xpos(1:3),geo(1)*RAD2DEG,geo(2)*RAD2DEG,geo(3),NM(isit).nobs/nobs    
        WRITE(lfnpos,'(I4,4I3,F5.1,A2,A4,3F13.3,3F9.3,3F13.3,3F23.12,I5,1X,6I3,F10.3)') iy,imon,id,ih,im,sec,ab,SIT(isit).name, &
              PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,PM(ipar,isit).xsig,PM(ipar+1,isit).xsig, &
              PM(ipar+2,isit).xsig,SL(isit).xpos(1:3),geo(1)*RAD2DEG,geo(2)*RAD2DEG,geo(3),SIT(isit)%nsat,&
              SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl,SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl,SIT(isit).dop(1) !PDOP              
      ELSE
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'STAPX')
        IF (ipar .EQ. 0) THEN
          WRITE(OUTPUT_UNIT,'(A,A,I6,F9.2,A2,64X,"  (",I2,F6.1,")")') run_tim(), SIT(isit).name, &
              CKF.mjd,CKF.sod,ab,SIT(isit)%nsat,NM(isit).esig
          CYCLE
        END IF
        ipar=OB(isit).ltog(ipar,1)
        ab=' '
        ! IF (PM(ipar,isit).iobs .LT. 5) THEN   !和位置参数有关的观测方程数,对于RAW模式会出错
        IF (SIT(isit)%nsat .LT. 5) THEN !根据有效卫星判断最合适
          ab='* '
          SL(isit)%fix_ewl=0
          SL(isit)%fix_wl=0
          SL(isit)%fix_nl=0
          SL(isit)%fixnum_ewl=0
          SL(isit)%fixnum_wl=0
          SL(isit)%fixnum_nl=0
        ELSE IF (SIT(isit).nfixnl .GT. 4) THEN
          ab='x '
        END IF
        !IF (((SL(isit).ncad-SL(isit).nfix).LE.2) .AND. (2*nobs*SL(isit).nfix.GE.NM(isit).nobs)) CKF.lamb=.TRUE.
        !IF (SL(isit).nfix+CKF.nl_maxdel.GE.NM(isit).nobs/nobs) CKF.lamb=.TRUE.
        SL(isit).xpos(1)=PM(ipar,isit).xest-SL(isit).refx(1)
        SL(isit).xpos(2)=PM(ipar+1,isit).xest-SL(isit).refx(2)
        SL(isit).xpos(3)=PM(ipar+2,isit).xest-SL(isit).refx(3)
        SL(isit).xpos(1:3)=SL(isit).xpos(1:3)*1.d2
        CALL matmpy(SL(isit).xpos,SL(isit).rot,SL(isit).xpos,1,3,3)
        ! IF (SIT(isit).skd(1:2) .NE. 'DP') THEN
        !   CALL matmpy(SL(isit).xpos,SL(isit).rot,SL(isit).xpos,1,3,3)
        ! END IF
        SL(isit).fpos(1)=SL(isit).fpos(1)-SL(isit).refx(1)
        SL(isit).fpos(2)=SL(isit).fpos(2)-SL(isit).refx(2)
        SL(isit).fpos(3)=SL(isit).fpos(3)-SL(isit).refx(3)
        SL(isit).fpos(1:3)=SL(isit).fpos(1:3)*1.d2
        !WRITE(1005,'(I5,4I3,F11.7,2X,A4,3F13.3)')iy,imon,id,ih,im,sec,SIT(isit).name,SL(isit).refx(1),SL(isit).refx(2),SL(isit).refx(3)
        CALL matmpy(SL(isit).fpos,SL(isit).rot,SL(isit).fpos,1,3,3)
        CALL xyzblh(PM(ipar:ipar+2,isit).xest,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,geo)
        IF (ab.EQ.'* ') THEN
          WRITE(OUTPUT_UNIT,'(A,A,I5,4I3,F5.1,A2,3F8.1,"  [",I3,I3,F6.1,"]",3F8.1,"  (",I2,F6.2,")",6I3)') run_tim(), SIT(isit).name, &
                iy,imon,id,ih,im,sec,ab,SL(isit).xpos(1:3),SL(isit).ncad,SL(isit).nfix,SL(isit).ratio,SL(isit).fpos(1:3), &
                SIT(isit)%nsat,NM(isit).esig,SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl, &
                SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl
        ELSE
          WRITE(OUTPUT_UNIT,'(A,A,I5,4I3,F5.1,A2,3F8.1,"  [",I3,I3,F6.1,"]",3F8.1,"  (",I2,F6.2,")",6I3)') run_tim(), SIT(isit).name, &
                iy,imon,id,ih,im,sec,ab,SL(isit).xpos(1:3),SL(isit).ncad,SL(isit).nfix,SL(isit).ratio,SL(isit).fpos(1:3), &
                ! (NM(isit).nobs-NM(isit)%ncstr)/nobs,NM(isit).esig,SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl, &
                SIT(isit)%nsat,NM(isit).esig,SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl, &
                SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl
        END IF
        IF (ab.EQ.'* ') THEN
          WRITE(lfnpos,'(I4,4I3,F5.1,A2,A4,3F13.3,3F9.3,3F13.3,3F23.12,I5,1X,6I3,F10.3)') iy,imon,id,ih,im,sec,ab,SIT(isit).name, &
          PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,PM(ipar,isit).xsig,PM(ipar+1,isit).xsig, &
          PM(ipar+2,isit).xsig,SL(isit).xpos(1:3),geo(1)*RAD2DEG,geo(2)*RAD2DEG,geo(3),SIT(isit)%nsat,&
          SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl,SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl,SIT(isit).dop(1) !PDOP        
        ELSE
          WRITE(lfnpos,'(I4,4I3,F5.1,A2,A4,3F13.3,3F9.3,3F13.3,3F23.12,I5,1X,6I3,F10.3)') iy,imon,id,ih,im,sec,ab,SIT(isit).name, &
          PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,PM(ipar,isit).xsig,PM(ipar+1,isit).xsig, &
          ! PM(ipar+2,isit).xsig,SL(isit).xpos(1:3),geo(1)*RAD2DEG,geo(2)*RAD2DEG,geo(3),(NM(isit).nobs-NM(isit)%ncstr)/nobs,&
          PM(ipar+2,isit).xsig,SL(isit).xpos(1:3),geo(1)*RAD2DEG,geo(2)*RAD2DEG,geo(3),SIT(isit)%nsat,&
          SL(isit).fix_ewl,SL(isit).fix_wl,SL(isit).fix_nl,SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl,SIT(isit).dop(1) !PDOP
        END IF

        !! xsy 2024-02-21: 输出成GPS周、周内秒
        ! gpsweek=(CKF.mjd-44244)/7
        ! wd=CKF.mjd-44244-gpsweek*7
        ! sow=wd*86400+ih*3600+im*60+sec
        ! IF (ab.EQ.'* ') THEN
        !   WRITE(lfnpos,'(I4,2X,F12.3,2X,A4,2X,3F13.3,3F9.3,2X,I5,F10.3,3I6,F8.1,F10.3)')gpsweek,sow,SIT(isit).name,&
        !       PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,&
        !       PM(ipar,isit).xsig,PM(ipar+1,isit).xsig,PM(ipar+2,isit).xsig,&
        !       SIT(isit)%nsat,SIT(isit).dop(1),SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl,SL(isit).ratio,SL(isit).adop        
        ! ELSE
        !   WRITE(lfnpos,'(I4,2X,F12.3,2X,A4,2X,3F13.3,3F9.3,2X,I5,F10.3,3I6,F8.1,F10.3)')gpsweek,sow,SIT(isit).name,&
        !       PM(ipar,isit).xest,PM(ipar+1,isit).xest,PM(ipar+2,isit).xest,&
        !       PM(ipar,isit).xsig,PM(ipar+1,isit).xsig,PM(ipar+2,isit).xsig,&
        !       SIT(isit)%nsat,SIT(isit).dop(1),SL(isit).fixnum_ewl,SL(isit).fixnum_wl,SL(isit).fixnum_nl,SL(isit).ratio,SL(isit).adop
        ! END IF
      END IF
    END DO
  END IF

  !! residuals
  IF (lfnres .NE. 0) THEN
    WRITE(lfnres,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
    DO isit=1, CKF.nsit
      i=1
      DO WHILE(i .LE. NM(isit).nobs)
        isat =NM(isit).ipob(i,1)
        ifreq=NM(isit).ipob(i,4)
        ! 观测值类型:1[相位],2[伪距],3[电离层约束],4[对流层约束],5[Station]
        IF (NM(isit).ipob(i,2).GT.2) THEN
          ! 约束信息不输出
          phase=NM(isit).resi(i)/NM(isit).weig(i)
          contype = ''
          IF (NM(isit).ipob(i,2).EQ.5) contype='STA'
          IF (NM(isit).ipob(i,2).EQ.4) contype='ZWD'
          IF (NM(isit).ipob(i,2).EQ.3) THEN
            IF (isat.EQ.0) THEN
              i=i+1
              CYCLE
            ELSE
              contype='ION'//CKF%cprn(isat)
            END IF
          END IF
          WRITE(lfnres,'(A4,1X,A6,1X,I1,1X,2D16.8,2F10.3)') SIT(isit).name,contype,NM(isit).ipob(i,2),&
              phase,NM(isit).weig(i),OB(isit).elev(isat)*RAD2DEG,OB(isit).azim(isat)*RAD2DEG  
          i=i+1
        ELSE
          IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs, 'CODE').NE.0) THEN
            phase=NM(isit).resi(i)/NM(isit).weig(i)
            pseud=NM(isit).resi(i+1)/NM(isit).weig(i+1)
          ELSE IF (INDEX(CKF.uobs,'PHASE') .NE. 0) THEN
            phase=NM(isit).resi(i)/NM(isit).weig(i)
            pseud=0.d0
          ELSE
            phase=0.d0
            pseud=NM(isit).resi(i)/NM(isit).weig(i)
          END IF
          IF (i .EQ. 1) THEN
            WRITE(lfnres,'(A4,1X,A3,1X,I1,1X,4D16.8,I3,F8.3,F9.3,I3,2F8.3)') SIT(isit).name,CKF.cprn(isat),ifreq,&
                phase,pseud,NM(isit).weig(i),NM(isit).weig(i+1),OB(isit).flag(isat,1),&
                OB(isit).elev(isat)*RAD2DEG,OB(isit).azim(isat)*RAD2DEG,OB(isit).nlosflag(isat),OB(isit).obs(isat,ifreq+3*MAXFREQ),NM(isit).esig
          ELSE
            WRITE(lfnres,'(A4,1X,A3,1X,I1,1X,4D16.8,I3,F8.3,F9.3,I3,F8.3)') SIT(isit).name,CKF.cprn(isat),ifreq,&
                phase,pseud,NM(isit).weig(i),NM(isit).weig(i+1),OB(isit).flag(isat,1),&
                OB(isit).elev(isat)*RAD2DEG,OB(isit).azim(isat)*RAD2DEG,OB(isit).nlosflag(isat),OB(isit).obs(isat,ifreq+3*MAXFREQ)
          
          END IF
          IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs, 'CODE').NE.0) THEN
            i=i+2
          ELSE
            i=i+1
          END IF
        END IF
      END DO
    END DO
  END IF

!  OPEN(UNIT=lfngm,file='gm.html')
!  DO WHILE(.TRUE.)
!   READ(lfngm,'(A)',END=100) line
!   IF (INDEX(line,'</body>') .NE. 0) EXIT
!  END DO
!  DO isit=1, CKF.nsit
!    ipar=pointer_string(OB(isit).npar,OB(isit).pname,'STAPX')
!    IF (ipar .NE. 0) THEN
!      ipar=OB(isit).ltog(ipar,1)
!      CALL xyzblh(PM(ipar:ipar+2,isit).xest,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,geo)
!      WRITE(lfngm,'(A,2F15.8,A)') '<script>gotoLocation(',geo(1)*RAD2DEG,geo(2)*RAD2DEG,')</script>'
!    END IF
!  END DO
!  WRITE(lfngm,'(A)') '</body>'
!  WRITE(lfngm,'(A)') '</html>'
 
!100 CLOSE(lfngm)

  RETURN

END SUBROUTINE
