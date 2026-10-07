!
!! purpose   : add observation equations to primary estimator
!! paraemters: 
!!             SCF    -- SRIF control parameters
!!             OB     -- Observation struct
!!             PM,NM  -- normal matrix & PAR table
!! author    : Geng J
!! created   : Nov. 11, 2007
!
SUBROUTINE ppp_addob_pri(CKF,SAT,OB,SIT,NM,PM,AM,infs,sion)
!!
!*
USE info
USE const
USE ckdctrl
USE ambiguity
USE satellite
USE observation
USE station
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT
TYPE(RNXOBS) :: OB
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
TYPE(AMBT) :: AM(1:*)
REAL(RL) :: infs(1:*)
REAL(RL) :: sion(MAXSAT)

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: ind,isat,ipar,i,j
  REAL(RL) :: phase,range,ewl,wl
  INTEGER(IT) :: ifreq,isys,iamb
  CHARACTER(LEN=5) :: camb
  CHARACTER(LEN=11) :: crfcb
  INTEGER(IT) :: pamb(MAXFREQ),pointer_string

  REAL(RL) :: sion_ele(MAXSYS)
  INTEGER(IT) :: sion_ref(MAXSYS),ind_ref(MAXSYS)

!  INTEGER(IT) :: k,nx,fg(MAXSAT),itg(MAXSAT)
!  REAL(RL) :: rx(MAXSAT),wx(MAXSAT),mean,sig,rms

  !*
  ! Start the exectuable code
  !!---------------------------

  i = COUNT(OB.omc(1:CKF.nprn,1) .NE. 0.D0)*2*MAXVAL(CKF.nfq)
  IF (NM.imtx+i+1 .GT. NM.nmtx) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_addob_pri): S matrix too small '
    CALL exit(1)
  END IF

  DO isat=1, CKF.nprn

    isys=INDEX(SYS,CKF.cprn(isat)(1:1))

    pamb=0
    DO ifreq=1, CKF.nfq(isys)

      WRITE(camb,'(A4,I1)') 'AMBL',ifreq
      WRITE(crfcb,'(A7,I1,A3)') 'RECDCBL',ifreq,CKF.cprn(isat)

      IF (OB.omc(isat,ifreq).EQ.0.d0 .OR. OB.omc(isat,MAXFREQ+ifreq).EQ.0.d0) CYCLE
      
      !! Melbourne-Wubbean combination observable for IF
      IF (CKF.cobs(1:2).EQ.'IF' .AND. CKF.nfreq(isys).GE.3) THEN
        iamb=pointer_string(OB.npar,OB.pname,'AMBL1')
        iamb=OB.ltog(iamb,isat)-NM.npc
        IF (AM(iamb).fewl .EQ.0.d0 .OR. AM(iamb).fwl .EQ.0.d0) GOTO 100
      END IF

      range=0.d0
      phase=0.d0

      !! right side & observation index
      IF (INDEX(CKF.uobs,'PHASE') .NE. 0) THEN
        NM.nobs=NM.nobs+1
        IF (CKF.cobs(1:2) .EQ. 'IF') THEN
          IF (CKF.nfreq(isys) .EQ. 2) THEN
            phase=OB.omc(isat,1)*SAT(isat).fac(1)-OB.omc(isat,2)*SAT(isat).fac(2)
            NM.weig(NM.nobs)=1.d0/DSQRT(SAT(isat).fac(1)**2*OB.var(isat,1)+SAT(isat).fac(2)**2*OB.var(isat,2))
          END IF
          ! IF (CKF.nfreq(isys) .GE. 3) THEN
          !   phase=SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(3))* &
          !         (SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(2))*OB.omc(isat,1) -&
          !          SAT(isat).freq(2)/(SAT(isat).freq(1)-SAT(isat).freq(2))*OB.omc(isat,2))-&
          !         SAT(isat).freq(3)/(SAT(isat).freq(1)-SAT(isat).freq(3))* &
          !         (SAT(isat).freq(2)/(SAT(isat).freq(2)-SAT(isat).freq(3))*OB.omc(isat,2) -&
          !          SAT(isat).freq(3)/(SAT(isat).freq(2)-SAT(isat).freq(3))*OB.omc(isat,3))
          !   NM.weig(NM.nobs)=SAT(isat).freq(1)**2/(SAT(isat).freq(1)-SAT(isat).freq(3))**2* &
          !          (SAT(isat).freq(1)**2/(SAT(isat).freq(1)-SAT(isat).freq(2))**2*OB.var(isat,1) -&
          !          SAT(isat).freq(2)/(SAT(isat).freq(1)-SAT(isat).freq(2))*OB.var(isat,2))+&
          !         SAT(isat).freq(3)**2/(SAT(isat).freq(1)-SAT(isat).freq(3))**2* &
          !         (SAT(isat).freq(2)**2/(SAT(isat).freq(2)-SAT(isat).freq(3))**2*OB.var(isat,2) -&
          !          SAT(isat).freq(3)**2/(SAT(isat).freq(2)-SAT(isat).freq(3))**2*OB.var(isat,3))
            
          !   phase=VEL_LIGHT/(SAT(isat).freq(1)-SAT(isat).freq(2))*(OB.obs(isat,1)-OB.obs(isat,2)-AM(ind).abewl)-&
          !        VEL_LIGHT/(SAT(isat).freq(2)-SAT(isat).freq(3))*(OB.obs(isat,2)-OB.obs(isat,3)-AM(ind).abwl)
          !   NM.weig(NM.nobs)=1.d0/DSQRT(NM.weig(NM.nobs))
          ! END IF
        ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
          phase=OB.omc(isat,ifreq)
          NM.weig(NM.nobs)=1.d0/DSQRT(OB.var(isat,ifreq))
        ELSE IF (CKF.cobs(1:7) .EQ. 'GRAPHIC') THEN
          phase=(OB.omc(isat,ifreq)+OB.omc(isat,MAXFREQ+ifreq))/2.d0
          NM.weig(NM.nobs)=2.d0/DSQRT(OB.var(isat,ifreq)+OB.var(isat,MAXFREQ+ifreq))
        END IF

        NM.ipob(NM.nobs,1)=isat
        NM.ipob(NM.nobs,2)=1    ! carrier-phase
        NM.ipob(NM.nobs,4)=ifreq
      END IF

      IF (INDEX(CKF.uobs, 'CODE') .NE. 0) THEN
        NM.nobs=NM.nobs+1
        IF (CKF.cobs(1:2) .EQ. 'IF') THEN
          IF (CKF.nfreq(isys) .EQ. 2) THEN
            range=OB.omc(isat,MAXFREQ+1)*SAT(isat).fac(1)-OB.omc(isat,MAXFREQ+2)*SAT(isat).fac(2)
            NM.weig(NM.nobs)=1.d0/DSQRT(SAT(isat).fac(1)**2*OB.var(isat,MAXFREQ+1)+SAT(isat).fac(2)**2*OB.var(isat,MAXFREQ+2))
          END IF
          ! ! need stable bias, otherwise it is easily effect by the bias
          ! IF (CKF.nfreq(isys) .GE. 3) THEN
          !   iamb=pointer_string(OB.npar,OB.pname,'AMBL1')
          !   iamb=OB.ltog(iamb,isat)-NM.npc
          !   IF (AM(iamb).fewl.NE.0.d0 .AND. AM(iamb).fwl.NE.0.d0 .AND. AM(iamb).ifab.GT.1) THEN
          !     wl =SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(2))*OB.omc(isat,1)-&
          !         SAT(isat).freq(2)/(SAT(isat).freq(1)-SAT(isat).freq(2))*OB.omc(isat,2)-AM(iamb).fwl*VEL_LIGHT/(SAT(isat).freq(1)-SAT(isat).freq(2))
          !     ewl=SAT(isat).freq(2)/(SAT(isat).freq(2)-SAT(isat).freq(3))*OB.omc(isat,2)-&
          !         SAT(isat).freq(3)/(SAT(isat).freq(2)-SAT(isat).freq(3))*OB.omc(isat,3)-AM(iamb).fewl*VEL_LIGHT/(SAT(isat).freq(2)-SAT(isat).freq(3))
          !     range=SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(3))*wl-&
          !           SAT(isat).freq(3)/(SAT(isat).freq(1)-SAT(isat).freq(3))*ewl
                    
          !     wl=(SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(2)))**2*OB.var(isat,1)+&
          !         (SAT(isat).freq(2)/(SAT(isat).freq(1)-SAT(isat).freq(2)))**2*OB.var(isat,2)
          !     ewl=(SAT(isat).freq(2)/(SAT(isat).freq(2)-SAT(isat).freq(3)))**2*OB.var(isat,2)+&
          !         (SAT(isat).freq(3)/(SAT(isat).freq(2)-SAT(isat).freq(3)))**2*OB.var(isat,3)
          !     !NM.weig(NM.nobs)=1.d0/DSQRT(wl+ewl)
          !     WRITE(*,*) CKF.cprn(isat),range,NM.weig(NM.nobs)!,AM(iamb).abewl,AM(iamb).abwl
          !   END IF
          ! END IF
        ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
          range=OB.omc(isat,MAXFREQ+ifreq)
          NM.weig(NM.nobs)=1.d0/DSQRT(OB.var(isat,MAXFREQ+ifreq))
          IF (.NOT.CKF%lrecdcb) THEN
            IF (ifreq.GE.3) NM.weig(NM.nobs)=1.D-12
          END IF
        ELSE IF (CKF.cobs(1:7) .EQ. 'GRAPHIC') THEN
          range=OB.omc(isat,MAXFREQ+ifreq)
          NM.weig(NM.nobs)=1.d0/DSQRT(OB.var(isat,MAXFREQ+ifreq))
        END IF

        NM.ipob(NM.nobs,1)=isat
        NM.ipob(NM.nobs,2)=2    ! code range
        NM.ipob(NM.nobs,4)=ifreq
      END IF

      !! check size
      IF (NM.imtx+NM.nobs+1 .GT. NM.nmtx) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_addob_pri): information matrix in primary filter'
        CALL exit(1)
      END IF

      !! clean
      DO j=1,NM.imtx+1
        IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs,'CODE').NE.0) THEN
          infs(NM.iptx(j)+NM.imtx+NM.nobs-1)=0.d0
        END IF
        infs(NM.iptx(j)+NM.imtx+NM.nobs)=0.d0
      END DO

      !! observation equations
      DO ipar=1,OB.npar
        IF (OB.ltog(ipar,isat) .EQ. 0) CYCLE

        IF (OB.pname(ipar)(1:7) .EQ. 'RECDCBL') THEN
          IF (OB.pname(ipar)(1:11) .NE. crfcb) CYCLE
        END IF

        ind=OB.ltog(ipar,isat)

        !@ CMT BY XSY: 各系统钟差
        ! IF (OB.pname(ipar)(1:6).EQ.'RECCLK' .AND. OB.pname(ipar)(7:7).NE.CKF.cprn(isat)(1:1)) CYCLE
        IF (OB.pname(ipar)(1:6).EQ.'RECCLK') THEN
          !@ CMT BY XSY: 参考系统:REC, 其他系统:REC+ISB
          IF (OB.pname(ipar)(7:7).NE.SYS(CKF.iref:CKF.iref) .AND. OB.pname(ipar)(7:7).NE.CKF.cprn(isat)(1:1)) CYCLE
        END IF

        !! initial value for the new ambiguites
        IF (OB.pname(ipar)(1:5) .EQ. TRIM(camb)) THEN
          pamb(ifreq)=ind
          AM(ind-NM.npc).iobs=AM(ind-NM.npc).iobs+1
          ! already added in ppp_abfix_emw
          IF (CKF.cobs(1:2).EQ.'IF' .AND. MAXVAL(CKF.nfreq).GE.3) AM(ind-NM.npc).iobs=AM(ind-NM.npc).iobs-1
          AM(ind-NM.npc).elev=AM(ind-NM.npc).elev+OB.elev(isat)*RAD2DEG
          IF (OB.flag(isat,ifreq) .NE. 0) THEN
            IF (CKF.cobs(1:7) .NE. 'GRAPHIC') THEN
              ! a priori iono-free ambiguit 
              AM(ind-NM.npc).xini=phase-range
              !@ CMT BY XSY: BE CARE OF THIS
              ! for raw, we set it as zero to deal with the multi-frequency ambiguities
              IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
                AM(ind-NM.npc).xini=0.d0
              END IF
            ELSE
              AM(ind-NM.npc).xini=(OB.omc(isat,ifreq)-OB.omc(isat,MAXFREQ+ifreq))/2.d0
            END IF
          END IF
          phase=phase-OB.amat(ipar,isat)*AM(ind-NM.npc).xini
         
          ! for raw obs, only ambigtuity term is estimated
          IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
            OB.amat(ipar,isat)=SAT(isat).lamda(ifreq)
          END IF
        ELSE IF (OB.pname(ipar)(1:3) .EQ. 'AMB') THEN
          CYCLE
        ELSE ! crfcb only for code
          IF (INDEX(CKF.uobs,'PHASE') .NE. 0) PM(ind).iobs=PM(ind).iobs+1
          IF (INDEX(CKF.uobs, 'CODE') .NE. 0) PM(ind).iobs=PM(ind).iobs+1
          ! IF (OB.pname(ipar)(1:6).EQ.'RECCLK') WRITE(*,*)SIT%name,ind,TRIM(PM(ind)%pname),PM(ind).iobs
        END IF

        ! Ionosphere delay on L1
        IF (OB.pname(ipar)(1:6) .EQ. 'ION'//CKF.cprn(isat)) THEN
          OB.amat(ipar,isat)=(SAT(isat).freq(1)/SAT(isat).freq(ifreq))**2
        END IF

        !! add coefficients
        !IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs,'CODE').NE.0 .AND. CKF.cobs(1:7).NE.'GRAPHIC') THEN
        IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs,'CODE').NE.0) THEN
          ! all parameter, expect for reccode
          IF (OB.pname(ipar)(1:11) .NE. crfcb) THEN
            infs(NM.iptx(ind)+NM.imtx+NM.nobs-1)=OB.amat(ipar,isat)*NM.weig(NM.nobs-1)
          END IF
          ! ionosphere, change rate
          IF (OB.pname(ipar)(1:6) .EQ. 'ION'//CKF.cprn(isat)) THEN
            infs(NM.iptx(ind)+NM.imtx+NM.nobs-1)=-OB.amat(ipar,isat)*NM.weig(NM.nobs-1)
          END IF
          ! amb only for phase
          IF (OB.pname(ipar)(1:5) .NE. camb) THEN
            infs(NM.iptx(ind)+NM.imtx+NM.nobs)=OB.amat(ipar,isat)*NM.weig(NM.nobs)
          END IF
          ! WRITE(5000,'(A3,1X,A20,I5,2F23.12)') CKF.cprn(isat),TRIM(OB.pname(ipar)),OB.ltog(ipar,isat), &
          !                                    infs(NM.iptx(ind)+NM.imtx+NM.nobs-1)/NM.weig(NM.nobs-1), &
          !                                    infs(NM.iptx(ind)+NM.imtx+NM.nobs)/NM.weig(NM.nobs)
        ELSE
          ! phase or code
          infs(NM.iptx(ind)+NM.imtx+NM.nobs)=OB.amat(ipar,isat)*NM.weig(NM.nobs)
          ! phase, ion
          IF (OB.pname(ipar)(1:6).EQ.'ION'//CKF.cprn(isat) .AND. INDEX(CKF.uobs,'PHASE').NE.0) THEN
            infs(NM.iptx(ind)+NM.imtx+NM.nobs)=-infs(NM.iptx(ind)+NM.imtx+NM.nobs)
          END IF
          ! phase, recdcb
          IF (OB.pname(ipar)(1:11).EQ.crfcb .AND. INDEX(CKF.uobs,'PHASE').NE.0) THEN
            infs(NM.iptx(ind)+NM.imtx+NM.nobs)=0.d0
          END IF
        END IF
      END DO

      !! add right-hand side
      !IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs,'CODE').NE.0 .AND. CKF.cobs(1:7).NE.'GRAPHIC') THEN
      IF (INDEX(CKF.uobs,'PHASE').NE.0 .AND. INDEX(CKF.uobs,'CODE').NE.0) THEN
        infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs-1)=phase*NM.weig(NM.nobs-1)
        infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs)  =range*NM.weig(NM.nobs)
      !! graphic observation is also in this condition
      ELSE IF (INDEX(CKF.uobs,'PHASE') .NE. 0) THEN
        infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs)  =phase*NM.weig(NM.nobs)
      ELSE
        infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs)  =range*NM.weig(NM.nobs)
      END IF

    END DO


 100 CONTINUE

    !! next observations
  END DO
  
  NM%ncstr=0

  ! Station constrains
  IF (CKF%lsitcons .EQ. .TRUE.) THEN
    DO ipar=1, OB%npar
      IF (OB%ltog(ipar,1) .EQ. 0) CYCLE
      IF (OB%pname(ipar)(1:4).NE.'STAP') CYCLE
      ind = OB%ltog(ipar,1)
      NM%nobs=NM%nobs+1
      NM%ncstr=NM%ncstr+1
      NM.ipob(NM%nobs,1)=0    !第几颗卫星
      NM.ipob(NM%nobs,2)=5    !观测值类型:1[相位],2[伪距],3[电离层约束],4[对流层约束],5[Station]
      NM.ipob(NM%nobs,3)=0    !不确定
      NM.ipob(NM%nobs,4)=0    !频点
      ! NM%weig(NM%nobs)=1.d0/CKF%sitConstrians(1)
      IF (OB%pname(ipar)(1:5).EQ.'STAPX') NM%weig(NM%nobs)=1.d0/CKF%sitConstrians(1)
      IF (OB%pname(ipar)(1:5).EQ.'STAPY') NM%weig(NM%nobs)=1.d0/CKF%sitConstrians(2)
      IF (OB%pname(ipar)(1:5).EQ.'STAPZ') NM%weig(NM%nobs)=1.d0/CKF%sitConstrians(3)
      infs(NM%iptx(ind)+NM%imtx+NM%nobs) = 1.d0*NM%weig(NM%nobs)  !add coffeicens
      infs(NM%iptx(NM%imtx+1)+NM%imtx+NM%nobs) = 0.d0          
    END DO
  END IF

  ! ZWD constrains
  IF (CKF%ztdmod(1:3).EQ.'FIX' .AND. SIT%ztdcor.NE.0.D0) THEN
    DO ipar=1,OB%npar
      IF (OB%ltog(ipar,1).EQ.0) CYCLE
      IF (OB%pname(ipar)(1:3).NE.'FIX') CYCLE
      ind = OB%ltog(ipar,1)
      NM%nobs=NM%nobs+1
      NM%ncstr=NM%ncstr+1
      NM%ipob(NM%nobs,1)=0
      NM%ipob(NM%nobs,2)=4
      NM%ipob(NM%nobs,3)=0
      NM%ipob(NM%nobs,4)=0
      NM%weig(NM.nobs)=1.d0/CKF%ztdConstrians
      infs(NM%iptx(ind)+NM%imtx+NM%nobs) = 1.d0*NM%weig(NM%nobs)
      infs(NM%iptx(NM%imtx+1)+NM%imtx+NM%nobs) = 0.d0
    END DO
  END IF

  ! Ionsphere constrains by satellite difference
  IF (CKF.ionmod(1:5) .EQ. 'FIXSD') THEN
    sion_ele = 0.d0
    sion_ref = 0
    ! select the refsat
    DO isys=1,CKF.nsys
      DO isat=1, CKF.nprn
        IF (CKF.cprn(isat)(1:1).NE.CKF.system(isys:isys)) CYCLE
        IF (OB.omc(isat,1).EQ.0.d0 .or. OB.omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE

        IF (sion(isat) .NE. 0.d0) THEN
          IF (OB.elev(isat) .GT. sion_ele(isys)) THEN
            sion_ele(isys) = OB.elev(isat)
            sion_ref(isys) = isat
          END IF
        END IF
      END DO
    END DO
    ! get the refsat ION flag
    DO isys=1,CKF.nsys
      DO isat=1,CKF.nprn
        IF (CKF.cprn(isat)(1:1).NE.CKF.system(isys:isys)) CYCLE
        IF (isat.NE.sion_ref(isys)) CYCLE
        DO ipar=1, OB.npar
          IF (OB.ltog(ipar,isat) .EQ. 0) CYCLE
          IF (OB.pname(ipar)(1:6).NE. 'ION'//CKF.cprn(isat)) CYCLE
          ind_ref(isys) = OB.ltog(ipar,isat)
        END DO
      END DO
    END DO
    ! add constrains
    DO isys=1,CKF.nsys
      DO isat=1,CKF.nprn
        IF (CKF.cprn(isat)(1:1).NE.CKF.system(isys:isys)) CYCLE
        IF (isat.EQ.sion_ref(isys)) CYCLE

        DO ipar=1, OB.npar
          IF (OB.ltog(ipar,isat) .EQ. 0) CYCLE
          IF (OB.pname(ipar)(1:6).NE. 'ION'//CKF.cprn(isat)) CYCLE
          ind = OB.ltog(ipar,isat)
          IF (PM(ind).iobs .LE. 0) CYCLE
          
          IF (sion(isat).NE.0.d0) THEN
            NM.nobs=NM.nobs+1
            NM%ncstr=NM%ncstr+1
            NM.ipob(NM.nobs,1)=isat         !第几颗卫星
            NM.ipob(NM.nobs,2)=3            !观测值类型:1[相位],2[伪距],3[电离层约束],4[对流层约束]
            NM.ipob(NM.nobs,3)=3            !为什么3?
            NM.ipob(NM.nobs,4)=MAXFREQ+1    !频点
            
            IF ((OB%elev(isat)+OB%elev(sion_ref(isys)))/2.d0*RAD2DEG .LE. 30.d0) THEN
              NM%weig(NM.nobs)=(1.d0/CKF%ionConstrians)*(2.d0*dsin((OB%elev(isat)+OB%elev(sion_ref(isys)))/2.d0))
            ELSE
              NM%weig(NM.nobs)=(1.d0/CKF%ionConstrians)
            END IF

            infs(NM.iptx(ind)+NM.imtx+NM.nobs) = 1.d0*NM.weig(NM.nobs)  !add coffeicens
            infs(NM.iptx(ind_ref(isys))+NM.imtx+NM.nobs) = -1.d0*NM.weig(NM.nobs)
            infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs) = 0.d0             !add right-hand side          
          END IF
        END DO
      END DO
    END DO
  END IF

  ! UnDifference constrains
  IF (CKF.ionmod(1:3) .EQ. 'FIX' .AND. CKF%ionmod(1:5) .NE. 'FIXSD') THEN

    DO isat=1, CKF.nprn
      IF (OB.omc(isat,1).EQ.0.d0 .or. OB.omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE
      DO ipar=1, OB.npar
        IF (OB.ltog(ipar,isat) .EQ. 0) CYCLE
        IF (OB.pname(ipar)(1:6).NE. 'ION'//CKF.cprn(isat)) CYCLE
        ind = OB.ltog(ipar,isat)

        IF (PM(ind).iobs .LE. 0) CYCLE
        !电离层产品存在则约束
        IF (sion(isat).NE.0.d0) THEN
          NM.nobs=NM.nobs+1
          NM%ncstr=NM%ncstr+1
          NM.ipob(NM.nobs,1)=isat         !第几颗卫星
          NM.ipob(NM.nobs,2)=3            !观测值类型:1[相位],2[伪距],3[电离层约束],4[对流层约束]
          NM.ipob(NM.nobs,3)=3            !不确定
          NM.ipob(NM.nobs,4)=MAXFREQ+1    !频点
          
          IF ((OB%elev(isat)+OB%elev(sion_ref(isys)))/2.d0*RAD2DEG .LE. 30.d0) THEN
            NM%weig(NM.nobs)=(1.d0/CKF%ionConstrians)*(2.d0*dsin((OB%elev(isat)+OB%elev(sion_ref(isys)))/2.d0))
          ELSE
            NM%weig(NM.nobs)=(1.d0/CKF%ionConstrians)
          END IF

          infs(NM.iptx(ind)+NM.imtx+NM.nobs) = 1.d0*NM.weig(NM.nobs)  !add coffeicens
          infs(NM.iptx(NM.imtx+1)+NM.imtx+NM.nobs) = 0.d0             !add right-hand side          
        END IF
      END DO
    END DO
  END IF

  ! WRITE(6000,*)CKF.mjd,CKF.sod
  ! DO i=1,NM.imtx+1
  !   WRITE(6000,'(<NM.imtx>(F8.3,2X))')(infs(NM.iptx(i)+j),j=1,NM.imtx)
  ! END DO
  RETURN

END SUBROUTINE
