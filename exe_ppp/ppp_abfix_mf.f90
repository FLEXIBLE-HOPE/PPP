!*
SUBROUTINE ppp_abfix_mf(CKF,SIT,SAT,UPD,NM,PM,AM,QM,SL)
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
  INTEGER(IT) :: i,j,k,m,isat,kpt,nfix(MAXSYS),nxl,ndl,maxtim,ifg(MAXSAT),ipt(MAXSAT),isys,ind,ntot,jsat
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
  ! To fix the EWL measurements
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
        rfcb(3,isys)=fxl
        
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

  ! L1 and WL
  DO i=1, NM.ns
    isat=AM(i).psat
    IF (TRIM(AM(i).pname) .EQ. 'AMBL3') CYCLE

    ! do not fix L1 and WL
    AB(i).abst=-1
    
    IF (TRIM(AM(i).pname) .EQ. 'AMBL2') then
      AB(i).abfr=AM(i-1).xini+AM(i-1).xcor-(AM(i).xini+AM(i).xcor)
    ELSE 
      AB(i).abfr=AM(i).xini+AM(i).xcor
    END IF
  END DO

  write(*,*) QM.ntot,QM.nxyz,QM.ndam,QM.ncad,QM.nfix
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

      write(*,*) QM.ntot,QM.nxyz,QM.ndam,QM.ncad,QM.nfix
      CALL ppp_abfix_lambda(AB,QM,QM.invx,CKF.nl_maxdel,CKF.nl_minsav,CKF.nl_ratio,SL.ratio)
      write(*,*) QM.ntot,QM.nxyz,QM.ndam,QM.ncad,QM.nfix

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
          IF (TRIM(AM(AB(i).pab).pname).NE.'AMBL1') CYCLE
          SIT.nfixnl=SIT.nfixnl+1
          IF (CKF.lamb.EQ..TRUE. .AND. AM(AB(i).pab).ifab.EQ.2) CYCLE
          AM(AB(i).pab).ifab=2
          AM(AB(i).pab).abnl=AB(i).abfx
        ELSE
          AM(AB(i).pab).abnl=AB(i).abfr
        END IF
      END DO
    ELSE
      IF (CKF.lamb .EQ. .TRUE.) THEN
        !! JG: output issues, to be sloved
        SL.nfix=QM.nfix
        DO i=1, NM.ns
          IF (AB(i).abst.EQ.1 .OR. AB(i).abst.EQ.2) THEN
            SIT.nfixnl=SIT.nfixnl+1
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
