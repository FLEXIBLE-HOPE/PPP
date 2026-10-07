!!
SUBROUTINE ppp_abfix_emw(CKF,isit,SAT,OB,NM,AM,UPD)
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
  iamb=pointer_string(OB.npar,OB.pname,'AMBL1')
  DO isat=1,CKF.nprn
    IF (OB.omc(isat,1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+1).NE.0.d0) THEN

      ind=OB.ltog(iamb,isat)-NM.npc
      
      wwl=SAT(isat).freq(2)/SAT(isat).freq(3)
      rwl=OB.obs(isat,2)-OB.obs(isat,3)-(wwl*OB.obs(isat,MAXFREQ+2)+ &
          OB.obs(isat,MAXFREQ+3))/(1.d0+wwl)/VEL_LIGHT*(SAT(isat).freq(2)-SAT(isat).freq(3))
      wwl=1.d0
      IF (OB.elev(isat)*RAD2DEG .LE. 30.d0) wwl=wwl*2.d0*DSIN(OB.elev(isat))
      IF (OB.flag(isat,1) .NE. 0) THEN
        AM(ind).abein=NINT(rwl)                ! a priori value, must be integer
      END IF
      rwl=rwl-AM(ind).abein
      AM(ind).xrewl=AM(ind).xrewl+wwl*rwl       ! mean
      AM(ind).eweig=AM(ind).eweig+wwl
      AM(ind).xsewl=AM(ind).xsewl+wwl*rwl**2    ! sigma
 
      AM(ind).iobs=AM(ind).iobs+1 
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
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.1)) AM(ind).abewl=0.5d0
        ! only one epoch
        IF (AM(ind).iobs .LT. 1) CYCLE
        ! no widelane FCB
        IF (UPD.ewfcb(isat) .EQ. 10.d0) CYCLE
        ! mean elevation angle
        !IF (AM(ind).elev/AM(ind).iobs .LE. CKF.cutoff) CYCLE
        ! no enough period
        ! IF ((AM(ind).ptime(2)-AM(ind).ptime(1))*86400.d0 .LT. CKF.minsec_common) CYCLE
        nxl=nxl+1
        rxl(nxl)=AM(ind).xrewl/AM(ind).eweig+AM(ind).abein-UPD.ewfcb(isat)
        IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.LT.1)) AM(ind).abewl=rxl(nxl)
        pxl(nxl)=DSQRT((AM(ind).xsewl-AM(ind).eweig*(AM(ind).xrewl/AM(ind).eweig)**2)/AM(ind).eweig/ &
                       !(AM(ind).iobs-1)+UPD.ewsl(isat)**2)
                       (AM(ind).iobs)+UPD.ewsl(isat)**2)
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
          rwl=AM(ind).abewl-fxl(isys,isit)
          isat=AM(ind).psat
          CALL prob_resol(rwl,pxl(ipt(ind)),1,CKF.wl_maxdev,CKF.wl_maxsig,alpha)

          IF (CKF.lamb .EQ. .TRUE.) THEN 
            IF (AM(ind).ifab .LT. 1) THEN
              IF (alpha .GT. CKF.wl_alpha) THEN
                AM(ind).ifab=1
                AM(ind).abewl=NINT(rwl)
                AM(ind).fewl=NINT(rwl)+fxl(isys,isit)+UPD.ewfcb(isat)
              ELSE
                AM(ind).abewl=rwl
                AM(ind).fewl=0.d0 !AM(ind).abewl
              END IF
            END IF
          ELSE

            IF (alpha .GT. CKF.wl_alpha) THEN
              AM(ind).ifab=1
              AM(ind).abewl=NINT(rwl) !+fxl(isys,isit)-UPD.ewfcb(isat)
              AM(ind).fewl=NINT(rwl)+fxl(isys,isit)+UPD.ewfcb(isat)
            ELSE
              AM(ind).abewl=rwl
              AM(ind).fewl=0.d0 !AM(ind).abewl
            END IF

          END IF
        END DO
      END IF
    END DO
  END IF

  RETURN

END SUBROUTINE
