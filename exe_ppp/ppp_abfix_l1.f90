!*
SUBROUTINE ppp_abfix_l1(CKF,SIT,SAT,UPD,NM,PM,AM,QM,SL)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE ambiguity
USE observation
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
  INTEGER(IT) :: i,j,k,isat,kpt,nfix(MAXSYS),nxl,ndl,maxtim,ifg(MAXSAT),ipt(MAXSAT),isys,ind,ntot
  REAL(RL) :: rxl(MAXSAT),wgt(MAXSAT),pxl(MAXSAT),fxl,vxl,sxl,maxele,alpha,rwl,wwl,g1,g2
  REAL(RL) :: rfcb(3,MAXSYS),disall(2),q22(3,3)
  TYPE(AMBD), POINTER :: AB(:)
  REAL(RL), ALLOCATABLE :: map(:,:),mapt(:,:),maptx(:,:),bias(:),qxx(:),invx(:,:)

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! prepare auxiliary AB
  IF (NM.ns .GT. 0) ALLOCATE(AB(NM.ns))

  ALLOCATE(map(QM.ntot,QM.ntot))
  ALLOCATE(mapt(QM.ntot,QM.ntot))
  ALLOCATE(maptx(QM.ntot,QM.ntot))
  ALLOCATE(invx(QM.ntot,QM.ntot))
  ALLOCATE(bias(QM.ntot))
  ALLOCATE(qxx(QM.ntot*(QM.ntot+1)/2))

  DO i=1, NM.ns
    ! do not hold, or hold but not fix
    IF (CKF.lamb.EQ..FALSE. .OR. (CKF.lamb.EQ..TRUE. .AND. AM(i).ifab.LT.2)) THEN
      AM(i).abewl=0.5d0
      AM(i).abwl=0.5d0
      AM(i).abnl=0.5d0
    END IF
    AB(i).pab=i
  END DO

  map=0.d0
  mapt=0.d0
  invx=0.d0
  ntot=QM.ntot
  DO i=1, QM.ntot
    map(i,i)=1.d0
    mapt(i,i)=1.d0
    DO j=i, QM.ntot
      invx(j,i)=QM.invx(j,i)
      invx(i,j)=QM.invx(j,i)
    END DO
  END DO

  ! map ambigtuiry for N1,N2,N5 to N1,Newl,Nwl
  DO i=QM.nxyz+1,QM.ntot
    IF (TRIM(AM(i-QM.nxyz).pname) .NE. 'AMBL1') CYCLE

    map(i,i)    = 1.d0
    map(i,i+1)  = 0.d0
    map(i,i+2)  = 0.d0
    map(i+1,i)  = 1.d0
    map(i+1,i+1)=-1.d0
    map(i+1,i+2)= 0.d0
    map(i+2,i)  = 0.d0
    map(i+2,i+1)= 1.d0
    map(i+2,i+2)=-1.d0

    mapt(i,i)    = 1.d0
    mapt(i+1,i)  = 0.d0
    mapt(i+2,i)  = 0.d0
    mapt(i,i+1)  = 1.d0
    mapt(i+1,i+1)=-1.d0
    mapt(i+2,i+1)= 0.d0
    mapt(i,i+2)  = 0.d0
    mapt(i+1,i+2)= 1.d0
    mapt(i+2,i+2)=-1.d0

  END DO

  DO i=1, QM.ntot
    DO j=i+1, QM.ntot
      QM.invx(i,j)=QM.invx(j,i)
    END DO
  END DO

  ! mapping raw ambiguity to L1, WL, EWL combination
  CALL matmpy(map,QM.invx(1:QM.ntot,1:QM.ntot),maptx,QM.ntot,QM.ntot,QM.ntot)
  CALL matmpy(maptx,mapt,QM.invx(1:QM.ntot,1:QM.ntot),QM.ntot,QM.ntot,QM.ntot)

  rfcb=10.d0
  ! To obtain the receiver FCB for EWL
  DO isys=1, CKF.nsys

    nfix(isys)=0
    nxl=0
    ipt=0

    DO i=1, NM.ns
      isat=AM(i).psat
      IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

      ! EWL ambiguities
      IF (TRIM(AM(i).pname) .NE. 'AMBL3') CYCLE

      ! no observation or FCB corrections, we set abst as no information
      IF (AM(i).iobs.EQ.0 .OR. UPD.ewfcb(isat).EQ.10.d0 .OR. AM(i).elev/AM(i).iobs.LE.CKF.cutoff) THEN
        AB(i).abst=-1
        AB(i).abfr=(AM(i-1).xini+AM(i-1).xcor)-(AM(i).xini+AM(i).xcor)
      ! with etimations as well as FCB corrections
      ELSE
        AB(i).abst=0
        AB(i).abfr=(AM(i-1).xini+AM(i-1).xcor)-(AM(i).xini+AM(i).xcor)-UPD.ewfcb(isat)
        nxl=nxl+1
        rxl(nxl)=AB(i).abfr-NINT(AB(i).abfr)
        wgt(nxl)=1.d0
        ifg(nxl)=0
        ipt(nxl)=i
      END IF
    END DO
 
    !! estimate receiver-specific FCB
    IF (nfix(isys).EQ.0 .AND. nxl.GT.0) THEN
      fxl=10.d0
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl,vxl,sxl)

      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.4d0) THEN
        rfcb(2,isys)=fxl
      END IF
    END IF
  END DO


  ! To fix the WL measurements
  DO isys=1, CKF.nsys

    nfix(isys)=0
    nxl=0
    ipt=0

    DO i=1, NM.ns
      isat=AM(i).psat
      IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

      ! EWL ambiguities
      IF (TRIM(AM(i).pname) .NE. 'AMBL2') CYCLE

      ! no observation or FCB corrections, we set abst as no information
      IF (AM(i).iobs.EQ.0 .OR. UPD.wfcb(isat).EQ.10.d0) THEN ! .OR. AM(i).elev/AM(i).iobs.LE.CKF.cutoff) THEN
        AB(i).abst=-1
        ! 2-3
        ! ambiguity in normal equation already        
        AB(i).abfr=(AM(i-1).xini+AM(i-1).xcor)-(AM(i).xini+AM(i).xcor)
      ! with etimations as well as FCB corrections
      ELSE
        AB(i).abst=0
        AB(i).abfr=(AM(i-2).xini+AM(i-2).xcor)-(AM(i-1).xini+AM(i-1).xcor)-UPD.wfcb(isat)
        nxl=nxl+1
        rxl(nxl)=AB(i).abfr-NINT(AB(i).abfr)
        wgt(nxl)=1.d0
        ifg(nxl)=0
        ipt(nxl)=i
        IF (AM(i).ifab .EQ. 2) THEN
          nfix(isys)=nfix(isys)+1
          AB(i).abst=1     ! WL fixed
          AB(i).abfx=AM(i).abwl
          rfcb(3,isys)=AB(i).abfr-AB(i).abfx
        END IF
      END IF
    END DO
 
    !! estimate receiver-specific FCB
    IF (nfix(isys).EQ.0 .AND. nxl.GT.0) THEN
      fxl=10.d0
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl,vxl,sxl)

      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.4d0) THEN
        rfcb(3,isys)=fxl
      END IF
    END IF
  END DO

  ! to fix the (extra-) wide-lane for narrow lane fixing
  nfix=0
  DO i=1, NM.ns

    IF (TRIM(AM(i).pname) .NE. 'AMBL1') CYCLE
    isat=AM(i).psat
    isys=INDEX(CKF.system,CKF.cprn(isat)(1:1))

    IF (AM(i+1).ifab.EQ.2 .AND. AM(i+2).ifab.EQ.2) THEN
      nfix(isys)=nfix(isys)+2
      CYCLE
    END IF

    IF (AB(i+1).abst.EQ.0 .AND. rfcb(2,isys).NE.10.d0) THEN 
      bias(2)=AB(i+1).abfr-rfcb(2,isys)
    ELSE
      bias(2)=AB(i+1).abfr
    END IF

    IF (AB(i+2).abst.EQ.0 .AND. rfcb(3,isys).NE.10.d0) THEN 
      bias(3)=AB(i+2).abfr-rfcb(3,isys)
    ELSE
      bias(3)=AB(i+2).abfr
    END IF

    qxx(1)=QM.invx(QM.nxyz+i+1,QM.nxyz+i+1)
    qxx(2)=QM.invx(QM.nxyz+i+2,QM.nxyz+i+1)
    qxx(3)=QM.invx(QM.nxyz+i+2,QM.nxyz+i+2)

    CALL ambslv(2,qxx(1:3),bias(2:3),disall)

    !write(*,*) CKF.cprn(isat),AM(i+1).ifab,AM(i+2).ifab,bias(2:3)
    IF (disall(2)/disall(1) .GT. CKF.nl_ratio) THEN
      nfix(isys)=nfix(isys)+1
      AB(i+1).abfx=bias(2)
      AB(i+1).abst=1
      AM(i+1).ifab=2
      AM(i+1).abewl=bias(2)

      nfix(isys)=nfix(isys)+1
      AB(i+2).abfx=bias(3)
      AB(i+2).abst=1
      AM(i+2).ifab=2
      AM(i+2).abwl =bias(3)
    END IF
  END DO


  DO i=1, NM.ns
    IF (TRIM(AM(i).pname) .NE. 'AMBL1') CYCLE
    AB(i).abst=-1
    AB(i).abfr=AM(i).xini+AM(i).xcor
  END DO

  IF (SUM(nfix) .GT. 0) THEN
    !! The following code is to get the value and variance for the
    !! unfixed ambiguities if there are some ambiguities are fixed

    !! place fixed ambiguities at the end
    CALL ppp_plc_fixed(AB,QM,QM.invx)

    !! impose constraints from already fixed ambiguities
    CALL ppp_add_ambcon(PM,AB,QM,QM.invx)

    !! place precise ambiguities at the end
    CALL ppp_plc_float(AB,QM,QM.invx)

      !! record ambiguity information
    DO i=1, NM.ns
      IF (AB(i).abst.EQ.1 .OR. AB(i).abst.EQ.2) THEN
        AM(AB(i).pab).ifab=2
        SELECT CASE(TRIM(AM(AB(i).pab).pname))
          CASE('AMBL1')
            AM(AB(i).pab).abnl=AB(i).abfx
          CASE('AMBL2')
            AM(AB(i).pab).abewl=AB(i).abfx
          CASE('AMBL3')
            AM(AB(i).pab).abwl=AB(i).abfx
        END SELECT
      ELSE
        SELECT CASE(TRIM(AM(AB(i).pab).pname))
          CASE('AMBL1')
            AM(AB(i).pab).abnl=AB(i).abfr
          CASE('AMBL2')
            AM(AB(i).pab).abewl=AB(i).abfr
          CASE('AMBL3')
            AM(AB(i).pab).abwl=AB(i).abfr
        END SELECT
      END IF
    END DO
  END IF



  !! to obtain the narrow lane
  DO isys=1, CKF.nsys

    nfix(isys)=0
    nxl=0
    ipt=0

    DO i=1, QM.ntot-QM.nxyz

      j=AB(i).pab

      isat=AM(j).psat
      IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

      ! NL ambiguities
      IF (TRIM(AM(j).pname) .NE. 'AMBL1') CYCLE

      g1=SAT(isat).freq(1)/(SAT(isat).freq(1)-SAT(isat).freq(2))
      g2=SAT(isat).freq(2)/(SAT(isat).freq(1)-SAT(isat).freq(2))

      ! no observation or FCB corrections, we set abst as no information
      IF (AM(j).iobs.EQ.0 .OR. UPD.nfcb(isat).EQ.10.d0 .OR. AM(j+2).ifab.LT.2 .OR. AM(j).elev/AM(j).iobs.LE.CKF.cutoff) THEN
        AB(i).abst=-1
        ! using the updated L1
        AB(i).abfr=AB(i).abfr
      ELSE
        AB(i).abst=0
        AB(i).abfr=AB(i).abfr+g2*UPD.wfcb(isat)-UPD.nfcb(isat)
        nxl=nxl+1
        rxl(nxl)=AB(i).abfr
        wgt(nxl)=1.d0
        ifg(nxl)=0
        ipt(nxl)=i
        IF (AM(j).ifab .EQ. 2) THEN
          nfix(isys)=nfix(isys)+1
          AB(i).abst=1 
          AB(i).abfx=AM(j).abnl
        END IF
      END IF

    END DO
 
    !! estimate receiver-specific FCB
    IF (nfix(isys).EQ.0 .AND. nxl.GT.0) THEN
      fxl=10.d0
      CALL proc_xl_dirc(nxl,rxl,wgt,ifg,ndl,fxl,vxl,sxl)

      IF ((nxl-ndl)*1.d0/nxl.GE.0.6d0 .AND. vxl.LE.0.4d0) THEN
        rfcb(1,isys)=fxl

        !! choose a reference satellite
        maxtim=0
        maxele=0.d0
        DO i=1,nxl
          IF (ifg(i).NE.0) CYCLE
          IF (AM(ipt(i)).iobs .GT. maxtim) THEN
            kpt=ipt(i)  ! which AM or AB
            maxtim=AM(ipt(i)).iobs
          ELSE IF (AM(ipt(i)).iobs.EQ.maxtim .AND. AM(ipt(i)).elev/AM(ipt(i)).iobs.GT.maxele) THEN
            kpt=ipt(i)
            maxele=AM(ipt(i)).elev/AM(ipt(i)).iobs
          END IF
        END DO
        nfix(isys)=1    ! fixed as datum, so-called presdu-fixed
        AB(kpt).abst=2
        AB(kpt).abfx=NINT(AB(kpt).abfr-fxl)
      END IF
    END IF
  END DO

  !! narrow-lane ambiguity fixing
  SL.ncad =0
  SL.nfix =0
  SL.ratio=0.d0
  IF (SUM(nfix) .GT. 0) THEN
    !! The following code is to get the value and variance for the
    !! unfixed ambiguities if there are some ambiguities are fixed

    !! place fixed ambiguities at the end
    CALL ppp_plc_fixed(AB,QM,QM.invx)

    !! impose constraints from already fixed ambiguities
    CALL ppp_add_ambcon(PM,AB,QM,QM.invx)

    !! place precise ambiguities at the end
    CALL ppp_plc_float(AB,QM,QM.invx)

    !! ambiguity resolution using LAMBDA method
    SL.ncad=QM.ndam+QM.nfix
    IF (QM.ndam .GT. 0) THEN

      CALL ppp_abfix_lambda(AB,QM,QM.invx,CKF.nl_maxdel,CKF.nl_minsav,CKF.nl_ratio,SL.ratio)
      !CALL ppp_abfix_lambda(AB,QM,QM.invx,INT((QM.ntot-QM.nxyz)/2),CKF.nl_minsav,CKF.nl_ratio,SL.ratio)

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
          IF (TRIM(AM(AB(i).pab).pname) .EQ. 'AMBL1') SIT.nfixnl=SIT.nfixnl+1
          IF (CKF.lamb.EQ..TRUE. .AND. AM(AB(i).pab).ifab.EQ.2) CYCLE
          AM(AB(i).pab).ifab=2
          SELECT CASE(TRIM(AM(AB(i).pab).pname))
            CASE('AMBL1')
              AM(AB(i).pab).abnl=AB(i).abfx
            CASE('AMBL2')
              AM(AB(i).pab).abewl=AB(i).abfx
            CASE('AMBL3')
              AM(AB(i).pab).abwl=AB(i).abfx
          END SELECT
        ELSE
          SELECT CASE(TRIM(AM(AB(i).pab).pname))
            CASE('AMBL1')
              AM(AB(i).pab).abnl=AB(i).abfr
            CASE('AMBL2')
              AM(AB(i).pab).abewl=AB(i).abfr
            CASE('AMBL3')
              AM(AB(i).pab).abwl=AB(i).abfr
          END SELECT
        END IF
      END DO
    ELSE
      IF (CKF.lamb .EQ. .TRUE.) THEN
        !! JG: output issues, to be slove
        SL.nfix=QM.nfix
        DO i=1, NM.ns
          IF (AB(i).abst.EQ.1 .OR. AB(i).abst.EQ.2) THEN
            IF (TRIM(AM(AB(i).pab).pname) .EQ. 'AMBL1') SIT.nfixnl=SIT.nfixnl+1
          END IF
        END DO
      END IF
    END IF
  END IF

  !! clean memory
  IF (NM.ns .GT. 0) DEALLOCATE(AB)

  DEALLOCATE(map)
  DEALLOCATE(mapt)
  DEALLOCATE(maptx)
  DEALLOCATE(bias)
  DEALLOCATE(qxx)
  DEALLOCATE(invx)

  RETURN

END SUBROUTINE
