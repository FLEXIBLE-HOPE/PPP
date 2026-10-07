!!
!! purpose  : intialize PM & NM for primary filter
!! parameter:
!!    input : SCF   -- configuration information
!!            SITE  -- station information
!!            OB    -- rinex observation
!!    output: NM,PM -- information matrix & parameters
!! author   : Geng J
!! created  : Nov. 10, 2007
!!
!*
SUBROUTINE ppp_init_pri(CKF,SAT,SIT,OB,NM,PM,isit)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
TYPE(CKDCFG) :: CKF
TYPE(SITE) :: SIT
TYPE(SATE) :: SAT(1:*)
TYPE(RNXOBS) :: OB
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
INTEGER(IT) :: isit
  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,l,ip,ic,is,isat,idx
  INTEGER(IT) :: ipar,minut,isys,ifreq,ksat
  REAL(RL) :: val
  CHARACTER(LEN=3) xyz
  CHARACTER(LEN=2) ns

  DATA xyz,ns /'XYZ','NS'/
  SAVE xyz,ns

  !*
  ! Start the exectuable code
  !!--------------------------

  !! initialization of parameters & information matrix
  ip=0
  ic=0
  is=0
  idx=0

  !! station dependent
  OB.npar=0
  OB.ltog=0
  IF (SIT.skd(1:1) .EQ. 'S') THEN
    DO i=1,3
      ic=ic+1
      ipar=NM.np+ic
      PM(ipar).pname='STAP'//xyz(i:i)
      PM(ipar).xini =SIT.x(i)
      OB.npar=OB.npar+1
      OB.pname(OB.npar) =PM(ipar).pname
      DO isat=1, CKF.nprn
        OB.ltog(OB.npar,isat)=ipar
      END DO
      NM.infs(ipar,ipar)=1.d0/SIT.dx0(i)
    END DO
  ELSE IF (SIT.skd(1:1).EQ.'K' .OR. SIT.skd(1:2).EQ.'DP' .OR. SIT.skd(1:2).EQ.'DK') THEN
    DO i=1,3
      ip=ip+1
      ipar=ip
      PM(ipar).pname='STAP'//xyz(i:i)
      PM(ipar).xini=SIT.x(i)
      ! White noise
      IF (SIT.skd(1:2) .EQ. 'DP') THEN
        PM(ipar).map=0.d0
        PM(ipar).rw=1.d0/SIT.qx(i)
      ! Random walk
      ELSE
        PM(ipar).map=1.d0
        PM(ipar).rw=1.d0/(SIT.qx(i)*DSQRT(CKF.dintv))
      END IF
      OB.npar=OB.npar+1
      OB.pname(OB.npar)=PM(ipar).pname
      DO isat=1, CKF.nprn
        OB.ltog(OB.npar,isat)=ipar
      END DO
      NM.infs(ipar,ipar)=1.d0/SIT.dx0(i)
    END DO
  !! reduced-dynamic orbit solution
  ELSE IF (SIT.skd(1:2) .EQ. 'DE') THEN
    isat=CKF.nprn+SIT.ileo
    DO i=1, SAT(isat).npar
      ip=ip+1
      ipar=ip
      PM(ipar).pname=TRIM(SAT(isat).pname(i))
      PM(ipar).ptime=CKF.mjd0+CKF.sod0/86400.D0
     
      !! we use the radam walk model for velocity and dynamic model
      PM(ipar).xini=SAT(isat).x(i)
      PM(ipar).map=1.d0
     
      ! PM(ipar).rw=1.d0/(SAT(isat).qx(i)*DSQRT(CKF.dintv))
      OB.npar=OB.npar+1
      OB.pname(OB.npar)=PM(ipar).pname
      DO ksat=1, CKF.nprn
        OB.ltog(OB.npar,ksat)=ipar
      END DO
      !! position: meter-level
      !! velocity: centimeter-level
      !! radiation pressure: 1
      !! drag: 2
      !! customer acceleration: 10-9
      SELECT CASE(TRIM(PM(ipar).pname))
        CASE('PXSAT','PYSAT','PZSAT')
          !! process noise: 0.5 m
          PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*5.d-4)
          !! initial variance in 10 m
          NM.infs(ipar,ipar)=1.d0/1.d-2
          ! PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*5.d-3)
          ! NM.infs(ipar,ipar)=1.d0/1.d-2
        CASE('VXSAT','VYSAT','VZSAT')
          !! process noise: 1 cm/s
          PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d-5)
          !! initial variance in 1 m/s
          NM.infs(ipar,ipar)=1.d0/1.d-3
          ! PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d-4)
          ! NM.infs(ipar,ipar)=2.d0/1.d-3
        CASE('DRAG_c')
          ! 2.0
          PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d-1)
          NM.infs(ipar,ipar)=1.d0/2.d0
          ! PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d0)
          ! NM.infs(ipar,ipar)=1.d0/1.d2
        CASE('SR_scale')
          ! 1.0
          PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d0)
          NM.infs(ipar,ipar)=1.d0/1.d0
          ! PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d0)
          ! NM.infs(ipar,ipar)=1.d0/1.d2
        CASE('EMP_Sa1','EMP_Ca1','EMP_Sc1','EMP_Cc1')
          ! 1.d-2
          PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d-2)
          NM.infs(ipar,ipar)=1.d0/1.d-2
          ! PM(ipar).rw=1.d0/(DSQRT(CKF.dintv)*1.d0)
          ! NM.infs(ipar,ipar)=1.d0/1.d2
        CASE DEFAULT
          PM(ipar).rw=1.d0/(SAT(isat).qx(i)*DSQRT(CKF.dintv))
          NM.infs(ipar,ipar)=1.d0/SAT(isat).dx0(i)
      END SELECT
      !write(*,*) PM(ipar).pname,PM(ipar).map,PM(ipar).rw,NM.infs(ipar,ipar)
    END DO
  ELSE IF (SIT.skd(1:1).NE.'F' .AND. SIT.skd(1:2).NE.'DF') THEN
    WRITE(ERROR_UNIT,'(2A)') '***ERROR(ppp_init_pri): SITE types unknown ', SIT.skd
    CALL exit(1)
  END IF

  !! receiver clock for ionosphere-free observable
  !! please note the ISB and IFB parameters maybe should be estimated as
  !! random parameters due to that the run
  DO i=1, CKF.nsys

    !! There must be GPS
    IF (CKF.system(i:i).EQ.'R' .AND. CKF.iref .EQ.INDEX(SYS,'G')) CYCLE

    isys=INDEX(SYS,CKF.system(i:i))

    IF (CKF%lisb) THEN
      IF (isys .EQ. CKF.iref) THEN
        ip=ip+1
        ipar=ip
        PM(ipar).pname='RECCLK'//CKF.system(i:i)
        ! White noise
        !IF (SIT.qclk(isys) .EQ. 0.d0) THEN
          PM(ipar).map=0.d0
          PM(ipar).rw=1.d0/SIT.dclk0(isys)
        ! Random walk
        !ELSE
        !  PM(ipar).map=1.d0
        !  PM(ipar).rw=1.d0/(SIT.qclk(isys)*DSQRT(CKF.dintv))
        !END IF
        PM(ipar).xini =0.d0
        OB.npar=OB.npar+1
        OB.pname(OB.npar) =PM(ipar).pname
        DO isat=1, CKF.nprn
          OB.ltog(OB.npar,isat)=ipar
        END DO
        NM.infs(ipar,ipar)=1.d0/SIT.dclk0(isys)
      ELSE
        ic=ic+1
        ipar=NM.np+ic
        PM(ipar).pname='RECCLK'//CKF.system(i:i)
        PM(ipar).map=0.d0
        PM(ipar).rw=0.d0

        PM(ipar).xini =0.d0
        OB.npar=OB.npar+1
        OB.pname(OB.npar) =PM(ipar).pname
        DO isat=1, CKF.nprn
          OB.ltog(OB.npar,isat)=ipar
        END DO
        NM.infs(ipar,ipar)=1.d0/SIT.dclk0(isys)
      END IF
    ELSE
      ip=ip+1
      ipar=ip
      PM(ipar).pname='RECCLK'//CKF.system(i:i)
      ! White noise
      !IF (SIT.qclk(isys) .EQ. 0.d0) THEN
        PM(ipar).rw=1.d0/SIT.dclk0(isys)
        PM(ipar).map=0.d0
      ! Random walk
      !ELSE
      !  PM(ipar).map=1.d0
      !  PM(ipar).rw=1.d0/(SIT.qclk(isys)*DSQRT(CKF.dintv))
      !END IF
      PM(ipar).xini =0.d0
      OB.npar=OB.npar+1
      OB.pname(OB.npar) =PM(ipar).pname
      DO isat=1, CKF.nprn
        OB.ltog(OB.npar,isat)=ipar
      END DO
      NM.infs(ipar,ipar)=1.d0/SIT.dclk0(isys)      
    END IF
  END DO

  !! IFB parameters for GLONASS
  IF (INDEX(CKF.system,'R').NE.0 .AND.  CKF.iref .EQ.INDEX(SYS,'G')) THEN
    DO isat=1, CKF.nprn
      IF (CKF.cprn(isat)(1:1) .EQ. 'R') THEN
        ic=ic+1
        ipar=NM.np+ic
        WRITE(PM(ipar).pname,'(A9,SP,I2)') 'RECCLK'//CKF.cprn(isat),SAT(isat).ifreq
        PM(ipar).map=0.d0
        PM(ipar).rw=0.d0
        PM(ipar).xini =0.d0
        OB.npar=OB.npar+1
        OB.pname(OB.npar)=PM(ipar).pname
        OB.ltog(OB.npar,isat)=ipar
        NM.infs(ipar,ipar)=1.d0/SIT.dclk0(INDEX(SYS,'R'))
      END IF
    END DO
  ELSE IF (CKF.iref .EQ. INDEX(SYS,'R')) THEN
    DO isat=1, CKF.nprn
      IF (CKF.cprn(isat)(1:1) .EQ. 'R') THEN
        IF (idx .EQ. 0) idx=isat
        IF (isat .EQ. idx) CYCLE
        ic=ic+1
        ipar=NM.np+ic
        WRITE(PM(ipar).pname,'(A9,SP,I2)') 'RECCLK'//CKF.cprn(isat),SAT(isat).ifreq
        PM(ipar).map=0.d0
        PM(ipar).rw=0.d0
        PM(ipar).xini =0.d0
        OB.npar=OB.npar+1
        OB.pname(OB.npar)=PM(ipar).pname
        OB.ltog(OB.npar,isat)=ipar
        NM.infs(ipar,ipar)=1.d0/SIT.dclk0(INDEX(SYS,'R'))
      END IF
    END DO
  END IF

  !@CMT BY XSY: [RECDCB] The third frequency code observation is a useless contribution because its' low wight
  IF (CKF%lrecdcb) THEN
    IF (TRIM(CKF.cobs) .EQ. 'RAW') THEN
      DO isat=1, CKF.nprn
        isys=INDEX(SYS,CKF.cprn(isat)(1:1))
        IF (CKF.nfq(isys) .GE. 3) THEN
          DO i=3,CKF.nfq(isys)
            ic=ic+1
            ipar=NM.np+ic
            WRITE(PM(ipar).pname,'(A7,I1,A3)') 'RECDCBL',i,CKF.cprn(isat)
            PM(ipar).xini=0.d0
            OB.npar=OB.npar+1
            OB.pname(OB.npar)=PM(ipar).pname
            OB.ltog(OB.npar,isat)=ipar
            NM.infs(ipar,ipar)=1.d-4

            ! ip=ip+1
            ! ipar=ip
            ! WRITE(PM(ipar).pname,'(A7,I1,A3)') 'RECDCBL',i,CKF.cprn(isat)
            ! PM(ipar)%pcode(1)=isit
            ! PM(ipar)%pcode(2)=isat
            ! PM(ipar)%pcode(3)=i

            ! PM(ipar).xini=0.d0
            ! ! Random walk
            ! PM(ipar).map=1.d0
            ! PM(ipar).rw=1.d0/0.01D0
            ! OB.npar=OB.npar+1
            ! OB.pname(OB.npar)=PM(ipar).pname
            ! OB.ltog(OB.npar,isat)=ipar
            ! NM.infs(ipar,ipar)=1.d-4
          END DO      
        END IF
      END DO
    END IF
  END IF
  
  !! zenith atmospheric delay
  ! IF (CKF%ztdmod(1:4).NE.'NONE' .AND. CKF%ztdmod(1:3).NE.'FIX' .AND. SIT%skd(1:1).NE.'D') THEN
  IF (CKF%ztdmod(1:4).NE.'NONE' .AND. SIT%skd(1:1).NE.'D') THEN
    ip=ip+1
    ipar=ip
    PM(ipar).pname=CKF.ztdmod
    PM(ipar).map  =1.d0
    PM(ipar).xini =0.d0
    i=INDEX(CKF.ztdmod,':')
    IF (i .NE. 0) READ(CKF.ztdmod(i+1:),*) minut
    IF (i .EQ. 0) THEN
      PM(ipar).rw =1.d0/(SIT.qztd*DSQRT(CKF.dintv/3600.d0))
    ELSE
      PM(ipar).rw =1.d0/(SIT.qztd*DSQRT(minut/60.d0))
    END IF
    OB.npar=OB.npar+1
    OB.pname(OB.npar)=PM(ipar).pname
    DO isat=1, CKF.nprn
      OB.ltog(OB.npar,isat)=ipar
    END DO
    NM.infs(ipar,ipar)=1.d0/SIT.dztd0
  END IF

  !! Gradient troposphere model
  IF (CKF.grdmod(1:4) .NE. 'NONE' .AND. SIT.skd(1:1).NE.'D') THEN
    DO j=1, 2
      ip=ip+1
      ipar=ip
      PM(ipar).pname=ns(j:j)//CKF.grdmod
      PM(ipar).map  =1.d0
      PM(ipar).xini =0.d0
      i=INDEX(CKF.grdmod,':')
      IF (i .NE. 0) READ(CKF.grdmod(i+1:),*) minut
      IF (i .EQ. 0) THEN
        PM(ipar).rw =1.d0/(SIT.qgrd*DSQRT(CKF.dintv/3600.d0))
      ELSE
        PM(ipar).rw =1.d0/(SIT.qgrd*DSQRT(minut/60.d0))
      END IF
      OB.npar=OB.npar+1
      OB.pname(OB.npar)=PM(ipar).pname
      DO isat=1, CKF.nprn
        OB.ltog(OB.npar,isat)=ipar
      END DO
      NM.infs(ipar,ipar)=1.d0/SIT.dgrd0
    END DO
  END IF

  IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
    DO isat=1, CKF.nprn
      isys=INDEX(SYS,CKF%cprn(isat)(1:1))
      IF (CKF%nfreq(isys).LT.2) CYCLE
      
      ip=ip+1
      ipar=ip
      PM(ipar).pname='ION'//CKF.cprn(isat)
      PM(ipar).pcode(1)=isit
      PM(ipar).pcode(2)=isat

      PM(ipar).map  =1.d0
      PM(ipar).xini =0.d0
      ! White noise
      IF (SIT.qion .EQ. 0.d0) THEN
        PM(ipar).rw=1.d0/SIT.dion0
        PM(ipar).map=0.d0
      ! Random walk
      ELSE
        PM(ipar).map=1.d0
        PM(ipar).rw=1.d0/(SIT.qion*DSQRT(CKF.dintv))
      END IF
      OB.npar=OB.npar+1
      OB.pname(OB.npar)=PM(ipar).pname
      OB.ltog(OB.npar,isat)=ipar
      NM.infs(ipar,ipar)=1.d0/SIT.dion0
    END DO
  END IF

  ! the broadcast ephemeris errors absorbing parameter 
  IF (CKF.lbds .EQ. .TRUE.) THEN
    DO isat=1, CKF.nprn
      ip=ip+1
      ipar=ip
      PM(ipar).pname='SISRE'//CKF.cprn(isat)
      PM(ipar).xini =0.d0
      ! White noise
      IF (SIT.qion .EQ. 0.d0) THEN
        PM(ipar).rw=1.d0/SIT.dion0
        PM(ipar).map=0.d0
      ! Random walk
      ELSE
        PM(ipar).map=1.d0
        PM(ipar).rw=1.d0/(SIT.qion*DSQRT(CKF.dintv))
      END IF
      OB.npar=OB.npar+1
      OB.pname(OB.npar)=PM(ipar).pname
      OB.ltog(OB.npar,isat)=ipar
      NM.infs(ipar,ipar)=1.d0/SIT.dion0
    END DO
  END IF

  !! estimate the time synchronization bias for LEO satellites as IFB, similar to GLONASS IFB
  IF (CKF%lleoifb) THEN
    DO isat=1, CKF%nprn
      IF (CKF%cprn(isat)(1:1) .EQ. 'L') THEN
        ic=ic+1
        ipar=NM%np+ic
        PM(ipar)%pname='RECCLK'//CKF%cprn(isat)
        PM(ipar)%map=0.d0
        PM(ipar)%rw=0.d0
        PM(ipar)%xini =0.d0
        OB%npar=OB%npar+1
        OB%pname(OB%npar)=PM(ipar)%pname
        OB%ltog(OB%npar,isat)=ipar
        NM%infs(ipar,ipar)=1.d0/1000.d0
      END IF
    END DO
  END IF

  !! ambiguity (dynamic), not saved in PM but in AM
  IF (INDEX(CKF.uobs,'PHASE') .NE. 0) THEN
    DO ifreq=1, MAXVAL(CKF.nfq)
      OB.npar=OB.npar+1
      WRITE(OB.pname(OB.npar),'(A4,I1)') 'AMBL',ifreq   !L1-L2-L3
    END DO
    DO isat=1, CKF.nprn
      !! revise kinematic PPP, no ambiguity parameters
      isys=INDEX(SYS,CKF.cprn(isat)(1:1))
      DO ifreq=1, CKF.nfq(isys)
        OB.ltog(OB.npar-MAXVAL(CKF.nfq)+ifreq,isat)=0
        IF (CKF.llog .EQ. .FALSE.) THEN
          OB.lifamb(isat,ifreq,1)=0.d0
          OB.lifamb(isat,ifreq,2)=0.d0
        END IF
      END DO
    END DO
  END IF

  !! output to screen
  ! WRITE(OUTPUT_UNIT,'(A)') SIT.name
  ! WRITE(OUTPUT_UNIT,'(A)') 'PARAMETERS NO.   NAME          INITIAL VARIANCE'
  ! DO ipar=1,NM.imtx
  !   val=NM.infs(ipar,ipar)
  !   WRITE(OUTPUT_UNIT,'(I14,2X,A20,2X,2(E20.14,1X))') ipar,PM(ipar).pname,val,PM(ipar).xini
  ! END DO
  ! STOP
  RETURN

END SUBROUTINE
