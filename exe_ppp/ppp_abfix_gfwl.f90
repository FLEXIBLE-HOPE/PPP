!!
SUBROUTINE ppp_abfix_gfwl(CKF,isit,SAT,OB,NM,AM,UPD)
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE satellite
USE observation
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
INTEGER(IT) :: isit

  !*
  ! The local variables
  !!-------------------------
  LOGICAL(LG) :: lfirst(MAXSIT)
  INTEGER(IT) :: i,ind,iamb,isat,nxl,ndl,isys,ifg(MAXSAT),ipt(MAXSAT)
  REAL(RL) :: rwl,wwl,rxl(MAXSAT),pxl(MAXSAT),wgt(MAXSAT),fxl(MAXSYS,MAXSIT),vxl,sxl,alpha

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_string

  DATA lfirst/MAXSIT*.TRUE./
  SAVE lfirst,fxl

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst(isit) .EQ. .TRUE.) then
    lfirst(isit)=.FALSE.
    fxl(:,isit)=10.d0
  END IF

  !! accumulate real-valued Melbourne-Wubbena combination observable
  iamb=pointer_string(OB.npar,OB.pname,'AMBL2')
  DO isat=1,CKF.nprn
    IF (OB.omc(isat,1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+1).NE.0.d0) THEN

      ind=OB.ltog(iamb,isat)-NM.npc

      wwl=SAT(isat).freq(1)/SAT(isat).freq(2)
      rwl=OB.obs(isat,1)-OB.obs(isat,2)- &
          (wwl*OB.obs(isat,MAXFREQ+1)+OB.obs(isat,MAXFREQ+2))/(1.d0+wwl)/VEL_LIGHT*(SAT(isat).freq(1)-SAT(isat).freq(2))
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
        IF (INDEX(AM(ind).pname,'AMBL2') .EQ. 0) CYCLE
        isat=AM(ind).psat
        IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.2)) AM(ind).abwl=0.5d0
        ! only one epoch
        IF (AM(ind).iobs .LE. 1) CYCLE
        ! no widelane FCB
        IF (UPD.wfcb(isat) .EQ. 10.d0) CYCLE
        ! mean elevation angle
        IF (AM(ind).elev/AM(ind).iobs .LE. CKF.cutoff) CYCLE
        ! no enough period
        ! IF ((AM(ind).ptime(2)-AM(ind).ptime(1))*86400.d0 .LT. CKF.minsec_common) CYCLE
        nxl=nxl+1
        rxl(nxl)=AM(ind).xrwl/AM(ind).weig+AM(ind).abin-UPD.wfcb(isat)
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.2)) AM(ind).abwl=rxl(nxl)
        pxl(nxl)=DSQRT((AM(ind).xswl-AM(ind).weig*(AM(ind).xrwl/AM(ind).weig)**2)/AM(ind).weig/&
                       (AM(ind).iobs-1)+UPD.wsl(isat)**2)
        ifg(nxl)=0
        wgt(nxl)=1.d0
        ipt(ind)=nxl
      END DO
      IF (nxl .EQ. 0) RETURN

      !! estimate receiver-specific FCB
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl(isys,isit),vxl,sxl)

      !! fix undifferenced widelane ambiguities
      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.25d0) THEN
        DO ind=1,NM.ns
          IF (ipt(ind) .EQ. -1) CYCLE
          IF (ifg(ipt(ind)) .NE. 0) CYCLE
          rwl=AM(ind).abwl-fxl(isys,isit)
          CALL prob_resol(rwl,pxl(ipt(ind)),1,CKF.wl_maxdev,CKF.wl_maxsig,alpha)


          IF (CKF.lamb .EQ. .TRUE.) THEN
            IF (AM(ind).ifab .LT. 2) THEN
              IF (alpha .GT. CKF.wl_alpha) THEN
                AM(ind).ifab=2
                AM(ind).abwl=NINT(rwl)
              ELSE
                AM(ind).abwl=rwl
              END IF
            END IF
          ELSE

            IF (alpha .GT. CKF.wl_alpha) THEN
              AM(ind).ifab=2
              AM(ind).abwl=NINT(rwl) !+fxl(isys,isit)-UPD.wfcb(isat)
            ELSE
              AM(ind).abwl=rwl
            END IF

          END IF
        END DO
      END IF
    END DO
  END IF

  RETURN

END SUBROUTINE
