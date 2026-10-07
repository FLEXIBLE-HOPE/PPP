!*
SUBROUTINE ppp_abfix_fcb_pri(CKF,SIT,SAT,UPD,NM,PM,AM,QM,SL)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
TYPE(AMBT) :: AM(1:*)
TYPE(INVM) :: QM
TYPE(SOL) :: SL
TYPE(FCB) :: UPD
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT

  !*
  ! The local variables
  !!-------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,isat,kpt,nfix(MAXSYS),nxl,ndl,maxtim,ifg(MAXSAT),ipt(MAXSAT),isys
  REAL(RL) :: rxl(MAXSAT),wgt(MAXSAT),fxl,vxl,sxl,maxele
  TYPE(AMBD), POINTER :: AB(:)

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! prepare auxiliary AB
  IF (NM.ns .GT. 0) THEN
    ALLOCATE(AB(NM.ns))
  ELSE
    RETURN
  END IF
  DO i=1, NM.ns
    IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(i).ifab.LT.4)) AM(i).abnl=0.5d0
    AB(i).pab=i
  END DO

  DO isys=1, CKF.nsys

    nfix(isys)=0
    nxl=0
    ipt=0

    DO i=1, NM.ns
      isat=AM(i).psat
      IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

      IF (AM(i).iobs.EQ.0 .OR. AM(i).ifab.LT.2 .OR. UPD.nfcb(isat).EQ.10.d0 .OR.(abs(int(AM(i).abwl)-AM(i).abwl).NE.0.d0) .OR.AM(i).elev/AM(i).iobs.LE.CKF.cutoff) THEN
        ! LN estimate
        AB(i).abst=-1
        AB(i).abfr=AM(i).xini+AM(i).xcor
      ELSE
        ! NL(N1)
        AB(i).abst=0
        AB(i).abfr=(AM(i).xini+AM(i).xcor)/SAT(isat).lamdn-AM(i).abwl/(SAT(isat).g-1.d0)-UPD.nfcb(isat)
        nxl=nxl+1
        rxl(nxl)=AB(i).abfr
        wgt(nxl)=1.d0
        ifg(nxl)=0
        ipt(nxl)=i  ! pointer to AM
        IF (AM(i).ifab .EQ. 4) THEN
          nfix(isys)=nfix(isys)+1
          AB(i).abst=1     ! narrowlane fixed
          AB(i).abfx=AM(i).abnl
        END IF
      END IF
    END DO
 
    !! estimate receiver-specific FCB
    IF (nfix(isys).EQ.0 .AND. nxl.GT.0) THEN
      fxl=10.d0
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl,vxl,sxl)
      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.4d0) THEN

        !! choose a reference satellite
        maxtim=0
        maxele=0.d0
        DO i=1,nxl
          IF (ifg(i).NE.0) CYCLE
          IF (AM(ipt(i)).iobs .GT. maxtim) THEN
            kpt=ipt(i)  ! which AM or AB
            maxtim=AM(ipt(i)).iobs
          ELSE IF (AM(ipt(i)).iobs.EQ.maxtim .AND. &
                AM(ipt(i)).elev/AM(ipt(i)).iobs.GT.maxele) THEN
            kpt=ipt(i)
            maxele=AM(ipt(i)).elev/AM(ipt(i)).iobs
          END IF
        END DO
        nfix(isys)=1    ! fixed as datum, so-called presdu-fixed
        AB(kpt).abst=2
        AB(kpt).abfx=NINT(AB(kpt).abfr-fxl)
      END IF
      ! DO i=1, NM.ns
      !   isat=AM(i).psat
      !   IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE
      !   IF (AB(i).abst .NE. -1) THEN
      !     AB(i).abfr = AB(i).abfr-fxl
      !   END IF
      ! END DO
    END IF
  END DO

  !! narrow-lane ambiguity fixing
  SL.ncad =0
  SL.nfix =0
  SL.ratio=0.d0
  IF (SUM(nfix) .GE. 0) THEN
    !! The following code is to get the value and variance for the
    !! unfixed ambiguities if there are some ambiguities are fixed
    !! so

    !! place fixed ambiguities at the end
    CALL ppp_plc_fixed(AB,QM,QM.invx)

    !! transform invx from LC to L1
    CALL ppp_map_invx(AM,SAT,AB,QM,QM.invx)

    !! impose constraints from already fixed ambiguities
    CALL ppp_add_ambcon(PM,AB,QM,QM.invx)

    !! place precise ambiguities at the end
    CALL ppp_plc_float(AB,QM,QM.invx)

    !! ambiguity resolution using LAMBDA method
    SL.ncad=QM.ndam+QM.nfix
    IF (QM.ndam .GT. 0) THEN

      CALL ppp_abfix_lambda(AB,QM,QM.invx,CKF.nl_maxdel,CKF.nl_minsav,CKF.nl_ratio,SL.ratio)

      !! store fixed widelane
      IF (QM.ncad .GT. 0) THEN
        SL.nfix=QM.ncad+QM.nfix

        !! sort inverted normal matrix
        CALL ppp_plc_fixed(AB,QM,QM.invx)

        !! further apply newly fixed widelane
        CALL ppp_add_ambcon(PM,AB,QM,QM.invx)
      END IF

      !! record ambiguity information
      DO i=1, NM.ns
        IF (AB(i).abst.EQ.1 .OR. AB(i).abst.EQ.2 .AND. QM.ncad.GT.0) THEN
          SIT.nfixnl=SIT.nfixnl+1
          IF (CKF.lamb.EQ..TRUE. .AND. AM(AB(i).pab).ifab.EQ.4) CYCLE
          AM(AB(i).pab).ifab=4
          AM(AB(i).pab).abnl=AB(i).abfx
          !write(*,*) CKF.cprn(AM(AB(i).pab).psat),AM(AB(i).pab).abnl
        ELSE
          AM(AB(i).pab).abnl=AB(i).abfr
        END IF
      END DO
    ELSE
      IF (CKF.lamb .EQ. .TRUE.) THEN
        !! JG: output issues, to be slove
        SL.nfix=QM.nfix
        DO i=1, NM.ns
          IF (AB(i).abst.EQ.1 .OR. AB(i).abst.EQ.2) THEN
            SIT.nfixnl=SIT.nfixnl+1
          END IF
        END DO
      END IF
    END IF
  END IF

  SL.fix_nl = 0
  SL.fixnum_nl = 0
  IF (SIT.nfixnl .GT. 4) THEN
    SL.fix_nl = 1
    SL.fixnum_nl = SIT.nfixnl
  END IF

  !! clean memory
  IF (NM.ns .GT. 0) DEALLOCATE(AB)

  RETURN

END SUBROUTINE