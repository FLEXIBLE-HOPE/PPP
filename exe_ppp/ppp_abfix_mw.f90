!!
SUBROUTINE ppp_abfix_mw(CKF,isit,SAT,OB,NM,AM,UPD,SL,SIT)
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE satellite
USE observation
USE station
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SATE) :: SAT(1:*)
TYPE(INFM) :: NM
TYPE(AMBT) :: AM(1:*)
TYPE(FCB) :: UPD
TYPE(SOL) :: SL
TYPE(SITE) :: SIT
INTEGER(IT) :: isit

  !*
  ! The local variables
  !!-------------------------
  LOGICAL(LG) :: lfirst(MAXSIT)
  INTEGER(IT) :: i,ind,iamb,isat,nxl,ndl,isys,ifg(MAXSAT),ipt(MAXSAT)
  REAL(RL) :: rwl,wwl,rxl(MAXSAT),pxl(MAXSAT),wgt(MAXSAT),fxl(MAXSYS,MAXSIT),vxl,sxl,alpha
  INTEGER(IT) :: irectype

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_string

  DATA lfirst/MAXSIT*.TRUE./
  SAVE lfirst,fxl

  !*
  ! Start the exectuable code
  !!----------------------------
  irectype = pointer_string(MAXRECTYPE,UPD%rectype,TRIM(SIT%rectyp))

  WRITE(1004,'((A),I7,F10.2,(A))') 'TIM',CKF%mjd,CKF%sod,'----------------------------------------------------------------IF'
  WRITE(1004,'(A)')'---PRN--------------ELE-------FLT-------SIG-------FIX-------UPD-----rBIAS--------CMT'
  IF (lfirst(isit) .EQ. .TRUE.) then
    lfirst(isit)=.FALSE.
    fxl(:,isit)=10.d0
  END IF
  SL.fixnum_wl = 0
  !! accumulate real-valued Melbourne-Wubbena combination observable
  iamb=pointer_string(OB.npar,OB.pname,'AMBL1')
  DO isat=1,CKF.nprn
    IF (OB.omc(isat,1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+1).NE.0.d0) THEN

      ind=OB.ltog(iamb,isat)-NM.npc

      isys=INDEX(SYS,CKF.cprn(isat)(1:1))

      IF (CKF.if_pco_corr) THEN
          !! xsy: only the Z PCO of satellite
          ! rwl=(OB.obs(isat,1) + SAT(isat).xyz(3,1)*SAT(isat).freq(1)/VEL_LIGHT)- &
          !     (OB.obs(isat,2) + SAT(isat).xyz(3,2)*SAT(isat).freq(2)/VEL_LIGHT)- &
          !     (SAT(isat).g*(OB.obs(isat,MAXFREQ+1)+SAT(isat).xyz(3,1)) + (OB.obs(isat,MAXFREQ+2)+SAT(isat).xyz(3,2)))/(1.d0+SAT(isat).g)/SAT(isat).lamdw
          !! xsy: the Z PCO of satellite and SIT PCO
          rwl=(OB.obs(isat,1) + (SIT.enu(3,1,isys)*dsin(OB.elev(isat)) + SAT(isat).xyz(3,1)*dcos(OB.nadir(isat)))*SAT(isat).freq(1)/VEL_LIGHT)- &
              (OB.obs(isat,2) + (SIT.enu(3,2,isys)*dsin(OB.elev(isat)) + SAT(isat).xyz(3,2)*dcos(OB.nadir(isat)))*SAT(isat).freq(2)/VEL_LIGHT)- &
              (SAT(isat).g*(OB.obs(isat,MAXFREQ+1)+SIT.enu(3,1,isys)*dsin(OB.elev(isat)) + SAT(isat).xyz(3,1)*dcos(OB.nadir(isat))) + &
              (OB.obs(isat,MAXFREQ+2)+SIT.enu(3,2,isys)*dsin(OB.elev(isat)) + SAT(isat).xyz(3,2)*dcos(OB.nadir(isat))))/(1.d0+SAT(isat).g)/SAT(isat).lamdw
      ELSE
          rwl=OB.obs(isat,1)-OB.obs(isat,2)- &
              (SAT(isat).g*OB.obs(isat,MAXFREQ+1)+OB.obs(isat,MAXFREQ+2))/(1.d0+SAT(isat).g)/SAT(isat).lamdw
      END IF

      IF (CKF%ArBiasMode .EQ. 'FCB') THEN
        IF (UPD%wfcb(isat) .NE. 10.d0) THEN
          rwl = rwl - UPD%wfcb(isat)
        END IF
      END IF
      !! xsy: correct the FPA(receiver-specific) bias
      IF (irectype .NE. 0) THEN
        rwl = rwl - UPD%recfcb(isat,2,irectype)
      END IF

      wwl=1.d0
      IF (OB.elev(isat)*RAD2DEG .LE. 30.d0) wwl=wwl*2.d0*DSIN(OB.elev(isat))
      IF (OB.flag(isat,1) .NE. 0) THEN
        AM(ind).abin=NINT(rwl)                ! a priori value, must be integer
      END IF
      rwl=rwl-AM(ind).abin
      AM(ind).xrwl=AM(ind).xrwl+wwl*rwl       ! mean
      AM(ind).weig=AM(ind).weig+wwl           ! weight
      AM(ind).xswl=AM(ind).xswl+wwl*rwl**2    ! sigma
    END IF
  END DO

  !! resolve ambiguities or not
  IF (CKF.liar) then

    DO isys=1, CKF.nsys
      !! correct for wide-lane FCBs and collect eligible widelane ambiguity estimates
      nxl=0
      ipt=-1
      DO ind=1, NM.ns
        isat=AM(ind).psat
        IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.2)) AM(ind).abwl=0.5d0
        ! only one epoch
        !IF (AM(ind).iobs .LE. 1) CYCLE
        IF (AM(ind).iobs .LT. 1) CYCLE
        ! mean elevation angle
        !IF (AM(ind).elev/AM(ind).iobs .LE. CKF.cutoff) CYCLE
        ! no enough period
        ! IF ((AM(ind).ptime(2)-AM(ind).ptime(1))*86400.d0 .LT. CKF.minsec_common) CYCLE

        IF (CKF%ArBiasMode .EQ. 'FCB') THEN
          IF (UPD%wfcb(isat) .EQ. 10.d0) THEN
            WRITE(1004,'((A),A4,A)')'WL',CKF%cprn(isat),'---------NONE WL FCB'
            CYCLE
          END IF
        END IF

        nxl=nxl+1
        rxl(nxl)=AM(ind).xrwl/AM(ind).weig+AM(ind).abin
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.2)) AM(ind).abwl=rxl(nxl)
        IF (CKF%ArBiasMode .EQ. 'OSB') THEN
          pxl(nxl)=DSQRT((AM(ind).xswl-AM(ind).weig*(AM(ind).xrwl/AM(ind).weig)**2)/AM(ind).weig/(AM(ind).iobs)+SAT(isat).posbstd(1)**2+SAT(isat).posbstd(2)**2)
        ELSE
          pxl(nxl)=DSQRT((AM(ind).xswl-AM(ind).weig*(AM(ind).xrwl/AM(ind).weig)**2)/AM(ind).weig/(AM(ind).iobs))
        END IF
        ifg(nxl)=0
        wgt(nxl)=1.d0
        ipt(ind)=nxl
      END DO
      IF (nxl .EQ. 0) CYCLE

      !! estimate receiver-specific FCB
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl(isys,isit),vxl,sxl)

      !! fix undifferenced widelane ambiguities
      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.25d0) THEN
        DO ind=1,NM.ns
          isat=AM(ind).psat
          IF (ipt(ind) .EQ. -1) CYCLE
          IF (ifg(ipt(ind)) .NE. 0) THEN
            WRITE(1004,'((A),A4,I3,1X,I3,F10.1,2X,A)')'WL',CKF%cprn(isat),ind,nxl,OB.elev(isat)*RAD2DEG,'---------BIG rBIAS'
            CYCLE
          END IF
          rwl=AM(ind).abwl-fxl(isys,isit)
          ! WRITE(1004,'((A),A4,I3,1X,I3,5F10.3)')'WL',CKF%cprn(isat),ind,nxl,AM(ind).abwl,fxl(isys,isit),rwl,pxl(ipt(ind)),UPD%wfcb(isat)
          CALL prob_resol(rwl,pxl(ipt(ind)),1,CKF.wl_maxdev,CKF.wl_maxsig,alpha)

          IF (CKF.lamb .EQ. .TRUE.) THEN
            IF (AM(ind).ifab .LT. 2) THEN
              IF (alpha .GT. CKF.wl_alpha) THEN
                AM(ind).ifab=2
                AM(ind).abwl=NINT(rwl)
                AM(ind).fwl=NINT(rwl)+fxl(isys,isit)
              ELSE
                AM(ind).abwl=rwl
                AM(ind).fwl=0.d0 !AM(ind).abwl
              END IF
            END IF
          ELSE
            IF (alpha .GT. CKF.wl_alpha) THEN
              AM(ind).ifab=2
              AM(ind).abwl=NINT(rwl)
              AM(ind).fwl=NINT(rwl)+fxl(isys,isit)
              SL.fixnum_wl=SL.fixnum_wl+1
            ELSE
              AM(ind).abwl=rwl
              AM(ind).fwl=0.d0 !AM(ind).abwl
            END IF
          END IF
          IF (AM(ind).ifab .GE. 2) THEN
            IF (irectype .NE. 0) THEN
              WRITE(1004,'((A),A4,I3,1X,I3,F10.1,5F10.3,5X,A)')'WL',CKF%cprn(isat),ind,nxl,OB.elev(isat)*RAD2DEG,rwl,pxl(ipt(ind)),AM(ind).abwl,UPD%wfcb(isat),UPD%recfcb(isat,2,irectype),'FIX_WL'
            ELSE
              WRITE(1004,'((A),A4,I3,1X,I3,F10.1,5F10.3,5X,A)')'WL',CKF%cprn(isat),ind,nxl,OB.elev(isat)*RAD2DEG,rwl,pxl(ipt(ind)),AM(ind).abwl,UPD%wfcb(isat),0.d0,'FIX_WL'
            END IF
          ELSE
            IF (irectype .NE. 0) THEN
              WRITE(1004,'((A),A4,I3,1X,I3,F10.1,5F10.3)')'WL',CKF%cprn(isat),ind,nxl,OB.elev(isat)*RAD2DEG,rwl,pxl(ipt(ind)),AM(ind).abwl,UPD%wfcb(isat),UPD%recfcb(isat,2,irectype)
            ELSE
              WRITE(1004,'((A),A4,I3,1X,I3,F10.1,5F10.3)')'WL',CKF%cprn(isat),ind,nxl,OB.elev(isat)*RAD2DEG,rwl,pxl(ipt(ind)),AM(ind).abwl,UPD%wfcb(isat),0.d0
            END IF
          END IF

        END DO
      END IF
    END DO
  END IF

  SL.fix_wl = 0
  IF (SL.fixnum_wl .GT. 4) THEN
    SL.fix_wl = 1
  END IF

  RETURN

END SUBROUTINE
