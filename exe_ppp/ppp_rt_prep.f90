!!
!! purpose  : detect bad data/cycle slips by checking change of WL and LG  for GNSS
!! parameter: 
!!    input : iepo -- epoch number
!!            CKF  -- configure information
!!    output: RP   -- saved WL & LG
!*
SUBROUTINE ppp_rt_prep(CKF,OB,SAT)
!!
!*
USE const
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! Start the exectuable
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB(1:*)
TYPE(SATE) :: SAT(MAXSAT)

  !*
  ! The local variables
  !!------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: isit,isat,isys,ljd(MAXSAT,MAXSIT), iy,imon,id,ih,im,interrupt_flag(MAXSIT),nx,fg(MAXSAT),k,j,ind
  REAL(RL) :: ltm(MAXSAT,MAXSIT),lsg(MAXSAT,MAXSIT),dlsg(MAXSAT,MAXSIT),lsw(MAXSAT,MAXSIT),lsewl(MAXSAT,MAXSIT),isec,lseewl(MAXSAT,MAXSIT),lshewl(MAXSAT,MAXSIT),inlg,inlw
  REAL(RL) :: lcode(MAXSAT,MAXSIT),lphase(MAXSAT,MAXSIT),ldoppler(MAXSAT,MAXSIT)

  REAL(RL) :: dt,lgr,dlgr,lw,lew,leewl,lhewl,wwl,lp,lc,ljump,rx(MAXSAT),wx(MAXSAT),mean,rms,sig

  REAL(RL) :: lgfif1(MAXSAT,MAXSIT),lgfif2(MAXSAT,MAXSIT),lgfif3(MAXSAT,MAXSIT)
  REAL(RL) :: gfif1,gfif2,gfif3,thre

  DATA lfirst /.TRUE./
  SAVE lfirst,ljd,ltm,lsg,dlsg,lsw,lsewl,lseewl,lshewl,lcode,lphase,ldoppler,lgfif1,lgfif2,lgfif3

  !*
  ! The function called
  !!------------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: pointer_string


  !*
  ! Start the exectuable code
  !!----------------------------
  
  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,isec)
  WRITE(1006,'((A),I5,4I3,F11.7,I7,F10.2,(A))')'TIM ',iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' [-> isit prn Dgf Dmw Demw Deemw Dhemw elev GF MW EMW EEMW HEMW <-]'
  IF (MOD(INT((CKF.mjd-CKF.mjd0)*86400.0+CKF.sod-CKF.sod0),CKF.ReConvTime) .EQ. 0) THEN
     lfirst=.TRUE.
  END IF

  IF (lfirst) THEN
    lfirst=.FALSE.
    DO isit=1, MAXSIT
      DO isat=1, CKF.nprn
        ljd(isat,isit)=0
        ltm(isat,isit)=0.d0
        lsg(isat,isit)=0.d0
        dlsg(isat,isit)=0.d0
        lsw(isat,isit)=0.d0
        lsewl(isat,isit)=0.d0
        lseewl(isat,isit)=0.d0
        lshewl(isat,isit)=0.d0

        lgfif1(isat,isit)=0.d0
        lgfif2(isat,isit)=0.d0
        lgfif3(isat,isit)=0.d0

        lcode(isat,isit)=0.d0
        lphase(isat,isit)=0.d0
        ldoppler(isat,isit)=0.d0
      END DO
    END DO
  END IF

  interrupt_flag=0
  !@ CMT BY XSY: 如果当前历元有效卫星数小于4，那么继续保存上一历元的LG、LW、LEW,否则会对周条探测和修复造成误区
  DO isit=1, CKF.nsit 
    IF (COUNT(OB(isit).obs(1:CKF.nprn,MAXFREQ+1).NE.0.d0) .LT. 4) THEN
      interrupt_flag(isit)=1
    END IF
  END DO 

  !! real-time preprocessing
  DO isat=1, CKF.nprn
    isys=INDEX(SYS,CKF.cprn(isat)(1:1))

    DO isit=1, CKF.nsit
      IF (interrupt_flag(isit) .EQ. 1) CYCLE
      IF (OB(isit).obs(isat,1) .EQ. 0.D0) CYCLE
      
      !@ CMT BY XSY: GAP is so large that a new ambiguity has to be set
      IF (ljd(isat,isit) .NE. 0) THEN
        dt=timdif(CKF.mjd,CKF.sod,ljd(isat,isit),ltm(isat,isit))
        !@ CMT BY XSY: for LEO, should be shorter
        IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
          IF (dt .GT. 60) THEN
            OB(isit).flag(isat,1:MAXFREQ)=1
            WRITE(*,'((A))') ' ... A big gap for : '//CKF.cprn(isat)
          END IF
        ELSE
          IF (dt .GT. CKF.gap) THEN
            OB(isit).flag(isat,1:MAXFREQ)=1
            WRITE(*,'((A))') ' ... A big gap for : '//CKF.cprn(isat)
          END IF
        END IF
      END IF

      lgr=0.d0
      dlgr=0.d0
      lw=0.d0
      lew=0.d0
      leewl=0.d0
      lhewl=0.d0
      gfif1=0.d0
      gfif2=0.d0
      gfif3=0.d0
      lp=0.d0
      lc=0.d0 

      !@ CMT BY XSY: THE BIASC FREQUENCY MUAT BE EXSIT
      IF (CKF%nfreq(isys).GE.2) THEN
        IF (OB(isit)%obs(isat,1).EQ.0.D0 .OR. OB(isit)%obs(isat,2).EQ.0.D0) CYCLE
        !! LG AND MW
        lgr=(SAT(isat).lamda(1)*OB(isit).obs(isat,1)-SAT(isat).lamda(2)*OB(isit).obs(isat,2))/(SAT(isat).lamda(2)-SAT(isat).lamda(1))
        lw=OB(isit).obs(isat,1)-OB(isit).obs(isat,2)-&
                (SAT(isat).g*OB(isit).obs(isat,MAXFREQ+1)+OB(isit).obs(isat,MAXFREQ+2))/(1.d0+SAT(isat).g)/SAT(isat).lamdw                    
      END IF

      IF (CKF%lhisi) THEN
        ind=2 ! WL:1 EWL:2
      ELSE
        ind=1
      END IF
      !@ CMT BY XSY: THE THIRD FREQUENCY IS NOT MUST BE EXSIT        
      IF (CKF.nfreq(isys) .GE. 3 .AND. OB(isit)%obs(isat,3).NE.0.D0) THEN
        wwl=SAT(isat).freq(ind)/SAT(isat).freq(3)
        lew=OB(isit).obs(isat,ind)-OB(isit).obs(isat,3)- &
                (wwl*OB(isit).obs(isat,MAXFREQ+ind)+OB(isit).obs(isat,MAXFREQ+3))/(1.d0+wwl)/VEL_LIGHT*(SAT(isat).freq(ind)-SAT(isat).freq(3))
        gfif1=(SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2))- &
               SAT(isat).lamda(2)*OB(isit).obs(isat,2)*(SAT(isat).freq(2)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2)))- &
              (SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(3)**2))- &
               SAT(isat).lamda(3)*OB(isit).obs(isat,3)*(SAT(isat).freq(3)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(3)**2)))        
      END IF
      IF (CKF.nfreq(isys) .GE. 4 .AND. OB(isit)%obs(isat,4).NE.0.D0) THEN        
        wwl=SAT(isat).freq(ind)/SAT(isat).freq(4)
        leewl=OB(isit).obs(isat,ind)-OB(isit).obs(isat,4)- &
                (wwl*OB(isit).obs(isat,MAXFREQ+ind)+OB(isit).obs(isat,MAXFREQ+4))/(1.d0+wwl)/VEL_LIGHT*(SAT(isat).freq(ind)-SAT(isat).freq(4))
        gfif2=(SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2))- &
               SAT(isat).lamda(2)*OB(isit).obs(isat,2)*(SAT(isat).freq(2)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2)))- &
              (SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(4)**2))- &
               SAT(isat).lamda(4)*OB(isit).obs(isat,4)*(SAT(isat).freq(4)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(4)**2)))                 
      END IF
      IF (CKF.nfreq(isys) .GE. 5 .AND. OB(isit)%obs(isat,5).NE.0.D0) THEN
        wwl=SAT(isat).freq(ind)/SAT(isat).freq(5)
        lhewl=OB(isit).obs(isat,ind)-OB(isit).obs(isat,5)- &
                (wwl*OB(isit).obs(isat,MAXFREQ+ind)+OB(isit).obs(isat,MAXFREQ+5))/(1.d0+wwl)/VEL_LIGHT*(SAT(isat).freq(ind)-SAT(isat).freq(5))
        gfif3=(SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2))- &
               SAT(isat).lamda(2)*OB(isit).obs(isat,2)*(SAT(isat).freq(2)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(2)**2)))- &
              (SAT(isat).lamda(1)*OB(isit).obs(isat,1)*(SAT(isat).freq(1)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(5)**2))- &
               SAT(isat).lamda(5)*OB(isit).obs(isat,5)*(SAT(isat).freq(5)**2/(SAT(isat).freq(1)**2-SAT(isat).freq(5)**2)))                  
      END IF

      !@ CMT BY XSY: FOR THE SINGLE FREQUENCY DETECT
      IF (CKF%nfreq(isys).EQ.1) THEN
        IF (OB(isit).obs(isat,2*MAXFREQ+1).NE.0.d0 .AND. OB(isit).obs(isat,1).NE.0.d0 .AND. lphase(isat,isit).NE.0.d0 .AND. ldoppler(isat,isit).NE.0.d0) THEN
          lp=OB(isit).obs(isat,1)-lphase(isat,isit)+(OB(isit).obs(isat,2*MAXFREQ+1)+ldoppler(isat,isit))*dt*0.5d0
          lc=(OB(isit).obs(isat,1+MAXFREQ)-lcode(isat,isit))+(OB(isit).obs(isat,2*MAXFREQ+1)+ldoppler(isat,isit))*dt*0.5d0*SAT(isat).lamda(1)
        END IF
      END IF

      !@ CMT BY XSY: CYCLE SLIP DECISION
      ! IF (CKF%cprn(isat)(1:1).NE.'L') THEN
        IF(CKF.nfreq(isys) .GE. 2)THEN
          !! new ambiguity cycle slip
          IF(ljd(isat,isit) .NE. 0) THEN
            WRITE(1006,'(I4,A5,11F14.3,A,3F14.3)')isit,CKF.cprn(isat),(lgr-lsg(isat,isit)),(lw-lsw(isat,isit)),(lew-lsewl(isat,isit)),&
                      (leewl-lseewl(isat,isit)),(lhewl-lshewl(isat,isit)),OB(isit)%elev(isat)*RAD2DEG,lgr,lw,lew,leewl,lhewl,'  GFIF:',&
                      (gfif1-lgfif1(isat,isit)),(gfif2-lgfif2(isat,isit)),(gfif3-lgfif3(isat,isit)) !,gfif1,gfif2,gfif3
            inlg=CKF%lg
            inlw=CKF%lw
            IF (DABS(lgr-lsg(isat,isit)).GT.inlg .OR. DABS(lw-lsw(isat,isit)).GT.inlw .OR. DABS(lew-lsewl(isat,isit)).GT.inlw .OR. &
                DABS(leewl-lseewl(isat,isit)).GT.inlw .OR. DABS(lhewl-lshewl(isat,isit)).GT.inlw) THEN              
                IF (.NOT.CKF%lhisi) THEN
                  !@ CMT BY XSY: Preventing the misdetection of a large number of cycle slip due to intense solar activity
                  IF (DABS(lgr-lsg(isat,isit)).GT.inlg .AND. DABS(lw-lsw(isat,isit)).LT.0.9D0 .AND. DABS(lew-lsewl(isat,isit)).LT.0.9D0 .AND. &
                      DABS(leewl-lseewl(isat,isit)).LT.0.9D0 .AND. DABS(lhewl-lshewl(isat,isit)).LT.0.9D0) THEN                      
                    ! WRITE(*,'((A))') ' ... BE MIND THE LARGE LG: '//CKF.cprn(isat)                                     
                  ELSE
                    OB(isit).flag(isat,1:MAXFREQ)=1
                    WRITE(*,'((A))') ' ... new amb (CYCLE CLIP): '//CKF.cprn(isat)                   
                  END IF
                ELSE
                  OB(isit).flag(isat,1:MAXFREQ)=1
                  WRITE(*,'((A))') ' ... new amb (CYCLE CLIP): '//CKF.cprn(isat)                       
                END IF
            END IF
            !@CMT BY XSY: GFIF [m] to detect the cycle slip, except for BDS-3 B1C and B1I
            IF (CKF%cprn(isat)(1:1).EQ.'C' .AND. pointer_string(CKF%nfreq(isys),CKF%freq(:,isys),'L1').NE.0 .AND. pointer_string(CKF%nfreq(isys),CKF%freq(:,isys),'L2').NE.0) THEN
              thre=2.0d0
            ELSE
              thre=0.1d0
            END IF
            !@CMT BY XSY: please be mind the 0.1
            ! question 1. For IGS station, it may be even smaller such as 0.08. 
            ! question 2. For kinematic data with MP error, this may be larger such as 0.2.
            ! question 3. For GPS L5, if no IFCB products, gfif combination will introduce 0.2 bias [1.3*0.15], it must be set lager in this suitation.
            IF (OB(isit).flag(isat,1).EQ.0 .AND. &
                (DABS(gfif1-lgfif1(isat,isit)).GT.thre .OR. &
                DABS(gfif2-lgfif2(isat,isit)).GT.thre .OR. &
                DABS(gfif3-lgfif3(isat,isit)).GT.thre)) THEN
                OB(isit).flag(isat,1:MAXFREQ)=1
                WRITE(*,'((A))') ' ... new amb (CYCLE CLIP): '//CKF.cprn(isat)//' GFIF'                
            END IF
          !! first ambiguity for this satellite
          ELSE
            OB(isit).flag(isat,1:MAXFREQ)=1
            ! WRITE(*,'((A))') ' ... new amb (NEW SATELLITE): '//CKF.cprn(isat)    
          END IF
        END IF

        IF(CKF.nfreq(isys) .EQ. 1)THEN
          !! new ambiguity cycle slip
          IF (ljd(isat,isit) .NE. 0) THEN
            IF ( DABS(lp).GT.CKF.lw) THEN
                WRITE(*,*) 'Single-freq Doppler CYCLE SLIP on '//CKF.cprn(isat),isit,lp
                OB(isit).flag(isat,1:MAXFREQ)=1
            END IF
          !! first ambiguity for this satellite
          ELSE
            OB(isit).flag(isat,1:MAXFREQ)=1   
            WRITE(*,'((A))') ' ... new amb (NEW SATELLITE): '//CKF.cprn(isat)    
          END IF
        END IF

      ! !@ CMT BY XSY: DOUBLE-DIFFERENCE LG CYCLE SLIP FOR LEO
      ! ELSE
      !   IF (CKF%nfreq(isys) .GE. 2) THEN
      !     dlgr=lgr-lsg(isat,isit)
      !   END IF
      !   IF (ljd(isat,isit) .NE. 0) THEN
      !       WRITE(1006,'(I4,A5,11F14.3)')isit,CKF.cprn(isat),(dlgr-dlsg(isat,isit)),(lw-lsw(isat,isit)),(lew-lsewl(isat,isit)),&
      !                 (leewl-lseewl(isat,isit)),(lhewl-lshewl(isat,isit)),OB(isit)%elev(isat)*RAD2DEG,dlgr,lw,lew,leewl,lhewl    
      !       inlg=0.2d0*CKF%dintv
      !       inlw=CKF%lw
      !       IF(DABS(dlgr-dlsg(isat,isit)).GT.inlg .OR. DABS(lw-lsw(isat,isit)).GT.inlw) THEN
      !           OB(isit)%flag(isat,1:MAXFREQ)=1
      !           WRITE(*,'((A))') ' ... new amb (CYCLE CLIP): '//CKF%cprn(isat)    
      !       END IF                            
      !   ELSE
      !     OB(isit)%flag(isat,1:MAXFREQ)=1
      !     ! WRITE(*,'((A))') ' ... new amb (NEW SATELLITE): '//CKF.cprn(isat)              
      !   END IF
      ! END IF

      !@ CMT BY XSY: LLI TO DETECT THE CYCLE SLIP
      IF (CKF%lli_flag) THEN
        IF (OB(isit)%flag(isat,1).EQ.0) THEN
          IF (COUNT(OB(isit)%lli(isat,:).EQ.1).NE.0) THEN
            OB(isit)%flag(isat,1:MAXFREQ)=1
            WRITE(*,'((A))') ' ... new amb (LLI): '//CKF.cprn(isat)    
          END IF
        END IF
      END IF

      !! reserve obs
      ljd(isat,isit)=CKF.mjd
      ltm(isat,isit)=CKF.sod

      lsg(isat,isit)=lgr
      dlsg(isat,isit)=dlgr
      lsw(isat,isit)=lw
      lsewl(isat,isit)=lew
      lseewl(isat,isit)=leewl
      lshewl(isat,isit)=lhewl

      lgfif1(isat,isit)=gfif1
      lgfif2(isat,isit)=gfif2
      lgfif3(isat,isit)=gfif3     

      lphase(isat,isit)=OB(isit).obs(isat,1)
      lcode(isat,isit)=OB(isit).obs(isat,1+MAXFREQ)
      ldoppler(isat,isit)=OB(isit).obs(isat,2*MAXFREQ+1)

    END DO
  END DO

  ! clock jump: code jump ,phase not,需要判断是伪距跳、还是相位跳，钟跳影响定位实际上是导致码相不一致，造成重收敛
  ! 1. 伪距跳-相位不跳: 影响定位，但是可以通过反向修复相位，使得相位与伪距一致
  ! 2. 相位跳-伪距不跳: 不存在
  ! 3. 伪距跳-相位跳  : 不影响定位
  ! 4. 伪距跳-相位不跳时，如果根据RTCM协议调了相位，可能会造成问题
  IF(CKF%nfreq(CKF%iref) .EQ. 1)THEN
     DO isit=1, CKF%nsit
       nx=0
       !OB(isit)%clkjmp=0.d0
       DO isat=1, CKF%nprn

         isys=INDEX(SYS,CKF%cprn(isat)(1:1))

         IF(OB(isit)%obs(isat,1) .NE. 0.d0) THEN

            ljump=(OB(isit)%obs(isat,1+MAXFREQ)-lcode(isat,isit))-((OB(isit)%obs(isat,1)-lphase(isat,isit)))*SAT(isat)%lamda(1)

            IF(OB(isit)%flag(isat,1) .NE. 1) THEN!将钟跳和周跳区分开
               nx=nx+1
               rx(nx)=ljump
               fg(nx)=0
               wx(nx)=1.d0
            END IF

            !! reserve obs
            ljd(isat,isit)=CKF%mjd
            ltm(isat,isit)=CKF%sod
            lphase(isat,isit)=OB(isit)%obs(isat,1)
            lcode(isat,isit)=OB(isit)%obs(isat,1+MAXFREQ)
            ldoppler(isat,isit)=OB(isit)%obs(isat,1+2*MAXFREQ)
         END IF
       END DO

       IF(nx .GT. 0) THEN
          CALL get_wgt_mean(.TRUE.,rx,fg,wx,nx,k,mean,rms,sig)!|rx-mean|>=3rms,flg=2,剔除里群值
          DO WHILE(nx-k.GT.2 .AND. rms.GT.30.d0)
            j=k
            CALL sign_robust(nx,rx,fg,10.d0,k)
            IF(k .EQ. j) EXIT
            CALL get_wgt_mean(.FALSE.,rx,fg,wx,nx,k,mean,rms,sig)
          END DO
          !!! 100.0 is not the optimal,
          IF(dabs(mean) .gt. 100.d0) then!码相历元间差分观测值rx的均值，如果均值过大，则认为发生钟跳
            write(*,*) 'CLOCK JUMP on '//CKF%cprn(isat),isit,mean
            OB(isit)%clkjmp=mean+OB(isit)%clkjmp
          END IF
       END IF

       ! correct the phase with accumulated clock jumps
       DO isat=1, CKF%nprn
         ! add ambiguity, but not the best approach, it will result reconvergence
         ! OB(isit).flag(isat,1:MAXFREQ)=1
         IF(OB(isit)%obs(isat,1) .NE. 0.d0) then
            OB(isit)%obs(isat,1)=OB(isit)%obs(isat,1)+OB(isit)%clkjmp/SAT(isat)%lamda(1)
         END IF
       END DO
     END DO
  END IF

  RETURN

END SUBROUTINE
