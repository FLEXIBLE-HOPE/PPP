!*
SUBROUTINE ppp_gps_model(CKF,SIT,OB,SAT,BHD,IOD,ifcb,sion)
!!
!*
USE const
USE ion
USE brdeph
USE station
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(SATE) :: SAT(MAXSAT),SATX
TYPE(SITE) :: SIT
TYPE(RNXOBS) :: OB
TYPE(CKDCFG) :: CKF
TYPE(BRDHEAD) :: BHD
TYPE(IONEX) :: IOD
REAL(RL) :: ifcb(MAXSAT)
REAL(RL) :: sion(MAXSAT)

  !*
  ! The local variables
  !!-----------------------
  LOGICAL(LG) :: flag(MAXSAT),lcont,lfind
  INTEGER(IT) :: i,j,k,ite,isat,jd_send,jd_recv,ierr,nx,fg(MAXSAT),itg(MAXSAT),isys,l
  REAL(RL) :: sod_send,sod_recv,drecclk(MAXSYS),dsatclk,mean,rms,sig,thred(2)
  REAL(RL) :: fjd,drate,tdelay,ddelay,dmap,ah,aw
  REAL(RL) :: sitrad(MAXFREQ),satrad(MAXFREQ),reldel(MAXFREQ),nadir(MAXFREQ)
  REAL(RL) :: r1(3,MAXFREQ),r2(3),r1leng(MAXFREQ),phase1(MAXFREQ),range1(MAXFREQ)
  REAL(RL) :: rx(MAXSAT),wx(MAXSAT),xsun(6),xlun(6),dx(3),iondel(MAXFREQ),dpcv(0:MAXPCVDEG,MAXFREQ,MAXSYS)
  REAL(RL) :: gmst,xpole,ypole,trpdel,ztdpart,dphwp(MAXFREQ),pcv(MAXFREQ,MAXSYS)
  REAL(RL) :: xant_f(6,MAXFREQ,MAXSYS),xant_j(6,MAXFREQ,MAXSYS),xsat_j(6,MAXFREQ),dump(3),dloudx(3)
  REAL(RL) :: scal,utcut1r,xhelp(2),rot_f2j(3,3),rot_rat(3,3),dxmat(3,3),dymat(3,3),sodut1,sec,sod0
  INTEGER(IT) :: mjdut1,npwc(MAXICS),lnpwc(MAXICS),ileo,lnpar, iy,imon,id,ih,im,mjd0
  REAL(RL) :: lpart(MAXPWC*MAXPWCEQUS),dlpart(MAXPWC*MAXPWCEQUS)
  REAL(RL) :: qm(4),rot_s2j(3,3),acr(3,3)
  REAL(RL) :: cbias(MAXFREQ),ecefsat(6),lmc,IFcode,IFphase,Rsat2sit(MAXSAT)
  REAL(RL) :: phi1(6,MAXICS),phi2(6,MAXICS)
  CHARACTER(LEN=1) :: conion
  SAVE sod0
  ! 用来计算卫星信号发射的惯性系和地固系的旋转矩阵
  REAL(RL) ::rot_f2j1(3,3),rot_rat1(3,3),dxmat1(3,3),dymat1(3,3),gmst1,xpole1,ypole1,sodut11,p1p2
  INTEGER(IT) :: mjdut11
  
  !*
  ! function called
  !!---------------------
  REAL(RL) :: dot, wetpp, timdif

  !*
  ! Start the exectuable code
  !!----------------------------

  OB.amat=0.d0
  OB.zmap=0.d0
  OB.azim=0.d0
  OB.elev=0.d0
  OB.azum=0.d0
  OB.nadir=0.d0
  OB.beta=0.d0
  OB.mu=0.d0
  OB.var =0.d0
  OB.omc =0.d0

  !! convert antenna offset from enu to xyz in earth fixed system, the unit is meter
  IF (SIT.skd(1:1) .NE. 'D') THEN
    DO k=1, CKF.nsys
      isys=INDEX(SYS,CKF.system(k:k))
      DO j=1, CKF.nfreq(isys)
        dump(1:3)=SIT.enu0(1:3)+SIT.enu(1:3,j,isys)
        CALL matmpy(SIT.rot_l2f,dump,dump,3,3,1)
        DO i=1,3
          xant_f(i,j,isys)=SIT.x(i)+dump(i)
          xant_f(i+3,j,isys)=0.d0
        END DO
      END DO
    END DO

    !! atmosphere correction
    !@ CMT BY XSY: 计算气象参数
    !CALL gpt2(CKF.mjd+CKF.sod/86400.D0,SIT.geod(1),SIT.geod(2),SIT.geod(3),1,0, &
    !          SIT.p0,SIT.t0,SIT.dt0,SIT.hr0,ah,aw,SIT.undu)
    CALL gpt2_1w(CKF.mjd+CKF.sod/86400.D0,SIT.geod(1),SIT.geod(2),SIT.geod(3),1,0, &
              SIT.p0,SIT.t0,SIT.dt0,dump(1),SIT.hr0,ah,aw,dump(2),SIT.undu)
    CALL drysaas(SIT.p0,SIT.geod(1),SIT.geod(3),SIT.zdd)
    SIT.hr0=1.d0/wetpp(SIT.t0,1.d0/SIT.hr0 )
    CALL wetsaas(SIT.hr0,SIT.t0,SIT.geod(1),SIT.geod(3),SIT.zwd)

    lnpar=0
    ileo=0
  ELSE
    ileo=SIT.ileo+CKF.nprn
    IF (SIT.skd(1:2) .EQ.'DE') THEN
      lnpar=SAT(ileo).npar
    ELSE
      lnpar=0
    END IF
    !! save the initial position and velocity
    SATX=SAT(ileo)
  END IF

  !! receiver clock correction
  drecclk=SIT.rclock/VEL_LIGHT

  !! if receiver clock is too bad, we repeat modeling
  ite = 0
  flag = .TRUE.

100 CONTINUE
  ite=ite+1

  DO k=1, CKF.nsys

    isys=INDEX(SYS,CKF.system(k:k))
    !CALL timinc(CKF.mjd,CKF.sod,-drecclk(isys),jd_recv,sod_recv)
    CALL timinc(OB.jd,OB.tsec,-drecclk(CKF.iref),jd_recv,sod_recv)

    !! Compute transformation matrix from earth-fixed to inertial system
    IF (CKF.lpost .EQ. .FALSE.) THEN
      CALL pth_orbit_igserp(CKF.flnerp,jd_recv,sod_recv,utcut1r,xhelp)
    ELSE
      !! IERS ERP parameters in UTC for EPO C04 products
      CALL timinc(jd_recv,sod_recv,OFF_GPS2TAI,jd_send,sod_send)
      CALL taiutc(jd_send,sod_send,jd_send,sod_send)

      CALL table_linear_interpolate('poleut1',.FALSE.,jd_send+sod_send/86400.d0,dump)
      xhelp(1:2)=dump(1:2)
      utcut1r=dump(3)
    END IF
    !@ CMT BY XSY: 需要用外部SOFA库，cmake在Release优化模式下编译的结果和Debug模式有差异
    CALL itrs2gcrs('IERS2010',jd_recv,sod_recv,utcut1r,xhelp,rot_f2j,rot_rat,dxmat,dymat,mjdut1,sodut1,gmst,xpole,ypole)

    IF (SIT.skd(1:1) .NE. 'D') THEN

      DO j=1, CKF.nfreq(isys)
        CALL matmpy(rot_f2j,xant_f(1,j,isys),xant_j(1,j,isys),3,3,1)
        CALL matmpy(rot_f2j,xant_f(4,j,isys),dx         ,3,3,1)
        CALL matmpy(rot_rat,xant_f(1,j,isys),xant_j(4,j,isys),3,3,1)
        xant_j(4:6,j,isys)=xant_j(4:6,j,isys)+dx(1:3)
      END DO

      !! read solar and lunar tables jd_send sod_send (AT)
      fjd=OFF_MJD2JD+jd_recv+(sod_recv-0.075d0+OFF_GPS2TT)/86400.d0
      CALL pleph(fjd,11,3,xsun)
      CALL pleph(fjd,10,3,xlun)
      xsun=xsun*1.d3
      xlun=xlun*1.d3

      DO j=1, CKF.nfreq(isys)
        !! tide dpsplacement
        CALL tide_displace(jd_recv,sod_recv,mjdut1,sodut1,xant_f(1,j,isys),xsun,xlun,&
            rot_f2j,SIT.rot_l2f,SIT.geod(1),SIT.geod(2),xpole,ypole,SIT.olc,dx)
        DO i=1,3
          xant_j(i,j,isys)=xant_j(i,j,isys)+dx(i)
        END DO
      END DO

    !@ CMT BY XSY: 几何学定轨,初值从orb文件输入
    ELSE IF (SIT.skd(1:2) .EQ. 'DP') THEN

      CALL leo_interpolate_orbit(CKF.flnleo,jd_recv,sod_recv,SIT.ileo,SAT(ileo).cprn,SAT(ileo).npar, &
                     SAT(ileo).pname,.TRUE.,.TRUE.,.FALSE.,xsat_j(1,1),xsat_j(4,1),lnpwc,lpart,dlpart)
      xsat_j(1:6,1) = xsat_j(1:6,1)*1.D3 !km2m

      CALL read_sat_asc(SAT(ileo).cprn,SAT(ileo).flnatt,jd_recv+sod_recv/86400.d0,qm,lfind)
      IF (lfind .EQ. .FALSE.) THEN
        CALL leo_att_non(SAT(ileo).type,xsat_j(1:6,1),rot_s2j)
      ELSE
        CALL q2rotmat(SAT(ileo).type,qm,rot_s2j)
        IF (SAT(ileo).type(1:5) .EQ. 'SWARM') THEN
          CALL matmpy(rot_f2j,rot_s2j,rot_s2j,3,3,3)
        END IF
      END IF

      ! Transform to along-track across-track and radial directions
      DO i=1, 3
        acr(1,1:3)=xsat_j(4:6,1)
        acr(3,1:3)=xsat_j(1:3,1)
      END DO
      CALL unit_vector(3,acr(1,1:3),acr(1,1:3),sig)
      CALL unit_vector(3,acr(3,1:3),acr(3,1:3),sig)
      CALL cross(acr(3,1:3),acr(1,1:3),acr(2,1:3))
      CALL unit_vector(3,acr(2,1:3),acr(2,1:3),sig)
      CALL matmpy(acr,rot_f2j,SIT.rot_l2f,3,3,3)
      
      ! CALL matmpy(rot_f2j,acr,SIT.rot_l2f,3,3,3)
      SAT(ileo).x=0.d0
      DO i=1,3
        DO j=1,3
          SAT(ileo).x(i)=SAT(ileo).x(i)+rot_f2j(j,i)*xsat_j(j,1)
        END DO
      END DO

      DO j=1, CKF.nfreq(isys)
        CALL matmpy(SAT(ileo).rant,SIT.enu(1:3,j,isys),dump,3,3,1) ! PCO
        dump=SAT(ileo).xyz0(1:3)+dump ! svnav
        CALL matmpy(rot_s2j,dump,dump,3,3,1)
        DO i=1,3
          xant_j(i,j,isys)=xsat_j(i,1)+dump(i)
          xant_j(i+3,j,isys)=xsat_j(i+3,1)
        END DO
      END DO
      !! orbit in ITRF
      SIT.x=0.d0
      DO i=1,3
        DO j=1,3
           SIT.x(i)=SIT.x(i)+rot_f2j(j,i)*xsat_j(j,1)
         END DO
      END DO

    !@ CMT BY XSY: 几何学定轨,初值从kin文件输入
    ELSE IF (SIT.skd(1:2) .EQ. 'DK') THEN

      CALL read_sat_asc(SAT(ileo).cprn,SAT(ileo).flnatt,jd_recv+sod_recv/86400.d0,qm,lfind)
      IF (lfind .EQ. .FALSE.) THEN
        rot_s2j=0.d0
        DO i=1, 3
          rot_s2j(i,i)=1.d0
        END DO
      ELSE
        CALL q2rotmat(SAT(ileo).type,qm,rot_s2j)
        IF (SAT(ileo).type(1:5) .EQ. 'SWARM') THEN
          CALL matmpy(rot_f2j,rot_s2j,rot_s2j,3,3,3)
        END IF
      END IF

      ! pos in icrf
      xant_j(:,:,isys)=0.d0
      DO j=1, CKF.nfreq(isys)
        CALL matmpy(rot_f2j,SIT.x(1),xant_j(1,j,isys),3,3,1)
        !IF (xsat_j(4,1).NE.0.d0) xant_j(4:6,j,isys)=xsat_j(4:6,1)
      END DO

      DO j=1, CKF.nfreq(isys)
        CALL matmpy(SAT(ileo).rant,SIT.enu(1:3,j,isys),dump,3,3,1)
        dump=SAT(ileo).xyz0(1:3)+dump
        CALL matmpy(rot_s2j,dump,dump,3,3,1)
        DO i=1,3
          xant_j(i,j,isys)=xant_j(i,j,isys)+dump(i)
        END DO
      END DO

    !@ CMT BY XSY: 简化动力学定轨,初值从orb文件输入,估计力学参数
    ELSE IF (SIT.skd(1:2) .EQ. 'DE') THEN
      ! Reference value for evaluation
      CALL leo_interpolate_orbit(CKF.flnleo,CKF.mjd,CKF.sod,SIT.ileo,SAT(ileo).cprn,SAT(ileo).npar, &
                    SAT(ileo).pname,.TRUE.,.TRUE.,.FALSE.,xsat_j(1,2),xsat_j(4,2),lnpwc,lpart,dlpart)
      SIT.dx0(1:6)=xsat_j(1:6,2)
      ! Transform to along-track across-track and radial directions
      DO i=1, 3
        acr(1,1:3)=xsat_j(4:6,2)
        acr(3,1:3)=xsat_j(1:3,2)
      END DO
      CALL unit_vector(3,acr(1,1:3),acr(1,1:3),sig)
      CALL unit_vector(3,acr(3,1:3),acr(3,1:3),sig)
      CALL cross(acr(3,1:3),acr(1,1:3),acr(2,1:3))
      CALL unit_vector(3,acr(2,1:3),acr(2,1:3),sig)
      SIT.rot_l2f=acr ! ACR to j2000

      CALL ppp_oi(CKF.mjd,CKF.sod,timdif(jd_recv,sod_recv,CKF.mjd,CKF.sod),SAT(ileo),xsat_j(1:6,1))
      xsat_j(1:6,1) = xsat_j(1:6,1)*1.d3

      l=0
      lpart=0.d0
      DO i=0, SAT(ileo).npar-1
        DO j=1, 3
          l=l+1
          lpart(l)=SAT(ileo).phi(i*6+j)
        END DO
      END DO
      ! write(*,*)lpart(1:SAT(ileo).npar*3)
      CALL read_sat_asc(SAT(ileo).cprn,SAT(ileo).flnatt,jd_recv+sod_recv/86400.d0,qm,lfind)
      IF (lfind .EQ. .FALSE.) THEN
        CALL leo_att_non(SAT(ileo).type,xsat_j(1:6,1),rot_s2j)
      ELSE
        CALL q2rotmat(SAT(ileo).type,qm,rot_s2j)
        IF (SAT(ileo).type(1:5) .EQ. 'SWARM') THEN
          CALL matmpy(rot_f2j,rot_s2j,rot_s2j,3,3,3)
        END IF
      END IF

      DO j=1, CKF.nfreq(isys)
        CALL matmpy(SAT(ileo).rant,SIT.enu(1:3,j,isys),dump,3,3,1)
        dump=SAT(ileo).xyz0(1:3)+dump
        CALL matmpy(rot_s2j,dump,dump,3,3,1)
        DO i=1,3
          xant_j(i,j,isys)=xsat_j(i,1)+dump(i)
          xant_j(i+3,j,isys)=xsat_j(i+3,1)
        END DO
      END DO
    END IF
  END DO

  !! loop over all satellites
  DO isat=1, CKF.nprn
    IF(OB.obs(isat,1).EQ.0.d0 .OR. OB.obs(isat,MAXFREQ+1).EQ.0.d0 .OR. .NOT.flag(isat)) CYCLE

    isys = INDEX(SYS,CKF.cprn(isat)(1:1))

    !### BLOCK 1
    !*** ITERATION OF SEND TIME and GET THE GEOMETRIC DISTANCE
    ddelay=0.1d0
    DO WHILE(dabs(ddelay) .GT. 1.d-9)

      !CALL timinc(CKF.mjd,CKF.sod,-drecclk(isys)-OB.delay(isat),jd_send,sod_send)
      CALL timinc(OB.jd,OB.tsec,-drecclk(CKF.iref)-OB.delay(isat),jd_send,sod_send)

      !! compute the satellite position and velocity coordinates
      IF (CKF.lpost .EQ. .TRUE.) THEN
        CALL everett_interp_orbit_post(CKF.flnorb,.TRUE.,.TRUE.,jd_send,sod_send,CKF.cprn(isat),xsat_j(1,1),xsat_j(4,1))
      ELSE
        CALL everett_interp_orbit(CKF.flnorb,.TRUE.,.TRUE.,jd_send,sod_send,CKF.cprn(isat),xsat_j(1,1),xsat_j(4,1))
      END IF

      IF (ALL(xsat_j(1:6,1) .EQ. 1.d15)) EXIT

      !! from km to m
      xsat_j(1:6,1) = xsat_j(1:6,1)*1.D3
      DO j=2, CKF.nfreq(isys)
        xsat_j(1:6,j) = xsat_j(1:6,1)
      END DO
    
      !! ECEF SATpos 参考orb2sp3,计算惯性系卫星位置旋转矩阵需要在卫星信号发射时刻计算
      ! IF (CKF.cprn(isat).EQ.'G04') THEN
      !   CALL itrs2gcrs('IERS2010',jd_send,sod_send,utcut1r,xhelp,rot_f2j1,rot_rat1,dxmat1,dymat1,mjdut11,sodut11,gmst1,xpole1,ypole1)
      !   DO i=1, 3
      !     ecefsat(i)=0.d0
      !     ecefsat(i+3)=0.d0
      !     DO j=1, 3
      !       ecefsat(i)=ecefsat(i)+rot_f2j1(j,i)*xsat_j(j,1)
      !       ! for nsat, rotation of earth need to be considered
      !       ecefsat(i+3)=ecefsat(i+3)+rot_f2j1(j,i)*xsat_j(j+3,1)+rot_rat1(j,i)*xsat_j(j,1) 
      !     END DO
      !   END DO
      !   WRITE(*,'(A4,I8,F14.8,6(F16.6,2X))')CKF.cprn(isat),jd_send,sod_send,ecefsat(1:6)
      ! END IF
      SAT(isat).satpos(1:6)=xsat_j(1:6,1)

      !! the satellite unit vectors (should be corrected for yaw error)
      fjd=OFF_MJD2JD+jd_send+(sod_send+OFF_GPS2TT)/86400.d0
      CALL pleph(fjd,11,3,xsun)
      xsun=xsun*1.d3
      CALL rot_scfix2j2000(jd_send,sod_send,CKF.cprn(isat),SAT(isat).csvn,SAT(isat).type,xsat_j(1,1),xsun, &
                           SAT(isat).xscf,SAT(isat).yscf,SAT(isat).zscf)

      DO j=1, CKF.nfreq(isys)
        DO i=1,3
          xsat_j(i,j)=xsat_j(i,j)+(SAT(isat).xyz(1,j)*SAT(isat).xscf(i)+&
                                   SAT(isat).xyz(2,j)*SAT(isat).yscf(i)+&
                                   SAT(isat).xyz(3,j)*SAT(isat).zscf(i))
        END DO
      END DO

      !! geometric distance
      DO j=1, CKF.nfreq(isys)
        DO i=1,3
          r1(i,j)=xsat_j(i,j)-xant_j(i,j,isys)
        END DO
        r1leng(j)=dsqrt(dot(3,r1(1,j),r1(1,j)))
      END DO
      Rsat2sit(isat) = r1leng(1)
      tdelay=r1leng(1)/VEL_LIGHT

      DO i=1,3
        dloudx(i)=r1(i,1)/r1leng(1)
      END DO

      !! decide whether to iterate the delay calculation
      ddelay=tdelay-OB.delay(isat)
      OB.delay(isat)=tdelay
    END DO

    ! satellite missing
    IF (ALL(xsat_j(1:6,1) .EQ. 1.d15)) CYCLE

    !### END OF BLOCK 1
    !###
    !### BLOCK 2
    !*** COMPUTE CORRECTON

    !! Compute the atmospheric corrections to the delay
    ! compute the station-satellite elevation angle and azimuth from north
    IF (SIT.skd(1:1) .NE. 'D') THEN
      CALL matmpy(r1(1,1),rot_f2j,r2,1,3,3)
      CALL matmpy(r2,SIT.rot_l2f,dump,1,3,3)
      OB.azim(isat)=DATAN2(dump(1),dump(2))
      OB.elev(isat)=DATAN(dump(3)/DSQRT(dump(1)**2+dump(2)**2))
    ELSE
      CALL matmpy(rot_s2j,SAT(ileo).bvec,r2,3,3,1)
      CALL unit_vector(3,r2,r2,sig)
      sig=dot(3,dloudx,r2)
      OB.elev(isat)=PI/2.d0-DACOS(sig)
      DO k = 1,3
        dx(k)=dloudx(k)-sig*r2(k)
      END DO
      CALL unit_vector(3,dx,dx,sig)
      CALL matmpy(rot_s2j,SAT(ileo).avec,dump,3,3,1)
      CALL unit_vector(3,dump,dump,sig)
      CALL cross(dump,r2,r2)
      OB.azim(isat)=DATAN2(dot(3,dx,r2),dot(3,dx,dump))
    END IF
    
    !! nadir angle for GPS satellite
    DO j=1, CKF.nfreq(isys)
      nadir(j)=dot(3,xsat_j(1,j),r1(1,j))/DSQRT(dot(3,xsat_j(1,j),xsat_j(1,j)))/DSQRT(dot(3,r1(1,j),r1(1,j)))
      nadir(j)=DACOS(nadir(j))
    END DO
    OB.nadir(isat)=nadir(1)
    DO j=1, 3
      dump(j)=-dloudx(j)-DCOS(nadir(1))*SAT(isat).zscf(j)
    END DO
    CALL unit_vector(3,dump,dump,OB.azum(isat))
    !OB.azum(isat)=dot(3,SAT(isat).xscf,dump)/DSQRT(dot(3,SAT(isat).xscf,SAT(isat).xscf))/DSQRT(dot(3,dump,dump))
    !OB.azum(isat)=DACOS(OB.azum(isat))
    CALL cross(SAT(isat).yscf,-dloudx,r2)
    OB.azum(isat)=DATAN2(dot(3,dump,r2),dot(3,dump,SAT(isat).yscf))
    CALL betau(xsat_j(1,1),xsun,OB.beta(isat),OB.mu(isat))

    IF (SIT.skd(1:1) .NE. 'D') THEN
      CALL VMF1_HT(ah,aw,CKF.mjd+CKF.sod/86400.D0,SIT.geod(1),SIT.geod(3),PI/2.D0-OB.elev(isat),dmap,ztdpart)
      trpdel=dmap*SIT.zdd+ztdpart*(SIT.zwd+SIT.ztdcor)
      !@ CMT BY XSY: runPPPtoION.py提取完整的对流层湿延迟[zwd],而lsq中仅提取湿延迟残余[ztdcor],如果跟LSQ保持一致,则注释下面三行
      IF (CKF%ztdmod(1:3).EQ.'FIX' .AND. SIT%ztdcor.NE.0.D0) THEN
        trpdel=dmap*SIT.zdd+ztdpart*SIT.ztdcor
      END IF
      OB.zmap(isat)=ztdpart

      SIT.grd(1) = -(SIT.zwd+SIT.ztdcor)/DTAN(OB.elev(isat))*DCOS(OB.azim(isat))   ! north
      SIT.grd(2) = -(SIT.zwd+SIT.ztdcor)/DTAN(OB.elev(isat))*DSIN(OB.azim(isat))   ! east
    ELSE
      trpdel=0.d0
    END IF
    SIT%trpdel(isat)=trpdel

    ! slant ionosphere delay
    IF (CKF.ionmod(1:3) .EQ. 'SID') THEN
      !CALL klobuchar(SIT.geod(1)*RAD2DEG,SIT.geod(2)*RAD2DEG,OB.elev(isat)*RAD2DEG,OB.azim(isat)*RAD2DEG,CKF.sod,BHD.ion(1:4,1,isys),BHD.ion(1:4,2,isys),drate)
      !SIT.ion(isat)=drate
      !DO j=1, CKF.nfreq(isys)
      !  iondel(j)=(GPS_L1/SAT(isat).freq(j))**2*drate
      !END DO

      ! slant ionosphere from GIM
!      CALL gim(CKF.mjd,CKF.sod,SIT.geod(1),SIT.geod(2),OB.elev(isat),OB.azim(isat),IOD,drate,rms)
!      IF (drate .NE. 9999.d0) THEN
!        SIT.ion(isat)=drate*40.28*1.d16/(GPS_L1)**2 !+SAT(isat).fac(2)*dcb(isat)*1.d-9*VEL_LIGHT
!        DO j=1, CKF.nfreq(isys)
!          iondel(j)=SAT(isat).freq(1)**2/SAT(isat).freq(j)**2*SIT.ion(isat)                    !drate*40.28*1.d16/SAT(isat).freq(j)**2
!        END DO
!      ELSE
!        SIT.ion(isat)=0.d0
!        iondel=0.d0
!      END IF

      ! slant ionosphere from other nearest stations
!      CALL read_ion(CKF.mjd,CKF.sod,CKF,sion)
!      IF (sion(isat) .NE. 0.d0) THEN
!        SIT.ion(isat)=sion(isat) !+SAT(isat).fac(2)*dcb(isat)*1.d-9*VEL_LIGHT
!        DO j=1, CKF.nfreq(isys)
!          iondel(j)=SAT(isat).freq(1)**2/SAT(isat).freq(j)**2*SIT.ion(isat)
!        END DO
!      ELSE
!        iondel=0.d0
!      END IF

      ! the initial value should be different for double and signle frequency
      IF (CKF.nfreq(isys) .GE. 2) THEN
        SIT.ion(isat)=(OB.obs(isat,MAXFREQ+1)-OB.obs(isat,MAXFREQ+2))/(1.d0-(SAT(isat).freq(1)/SAT(isat).freq(2))**2)
        DO j=1, CKF.nfreq(isys)
          iondel(j)=(SAT(isat).freq(1)/SAT(isat).freq(j))**2*SIT.ion(isat)
        END DO
      ELSE
        iondel=0.d0
      END IF
    ! zenith ionosphere delay
    ELSE IF (CKF.ionmod(1:3) .EQ. 'ZID') THEN
      iondel=0.d0
    ! model correction
    ELSE IF (CKF.ionmod(1:3) .EQ. 'BRD') THEN
      iondel=0.d0
      SELECT CASE(CKF.cprn(isat)(1:1))
        CASE('G','C','J','I')
          ! Klobuchar model: GPS
          !! The input sow should be in UTC, for simplicity, the GPST is used
          !! alpha and beta is not fixed, it should be changed based on time
          CALL klobuchar(SIT.geod(1)*RAD2DEG,SIT.geod(2)*RAD2DEG,OB.elev(isat)*RAD2DEG,OB.azim(isat)*RAD2DEG,CKF.sod,BHD.ion(1:4,1,isys),BHD.ion(1:4,2,isys),drate)
          IF (CKF.cprn(isat)(1:1) .EQ. 'I') THEN
            DO j=1, CKF.nfreq(isys)
              iondel(j)=(GPS_L5/SAT(isat).freq(j))**2*drate
            END DO
          ELSE
            DO j=1, CKF.nfreq(isys)
              iondel(j)=(GPS_L1/SAT(isat).freq(j))**2*drate
            END DO
          END IF
        CASE('E')
          ! NeQuick: Galileo
      END SELECT
    ! Global Ionosphere Model
    ELSE IF (CKF.ionmod(1:3) .EQ. 'GIM') THEN
      CALL gim(CKF.mjd,CKF.sod,SIT.geod(1),SIT.geod(2),OB.elev(isat),OB.azim(isat),IOD,drate,rms)
      IF (drate .NE. 9999.d0) THEN
        SIT.ion(isat)=drate*40.28*1.d16/(GPS_L1)**2 !+SAT(isat).fac(2)*dcb(isat)*1.d-9*VEL_LIGHT
        DO j=1, CKF.nfreq(isys)
          iondel(j)=drate*40.28*1.d16/SAT(isat).freq(j)**2
        END DO
      ELSE
        SIT.ion(isat)=0.d0
        iondel=0.d0
      END IF
    ELSE IF (CKF.ionmod(1:3) .EQ. 'FIX') THEN
      ! CALL read_ion(CKF.mjd,CKF.sod,CKF,sion)
      !@ CMT BY XSY: UD电离层产品,若是单差,此处需要修改
      IF (sion(isat) .NE. 0.d0) THEN
        SIT.ion(isat)=sion(isat)
        DO j=1, CKF.nfreq(isys)
          iondel(j)=SAT(isat).freq(1)**2/SAT(isat).freq(j)**2*SIT.ion(isat)
        END DO
      !@ CMT BY XSY: 没有电离层产品则用双频P码改正
      ELSE
        IF (CKF.nfreq(isys) .GE. 2) THEN
          SIT.ion(isat)=(OB.obs(isat,MAXFREQ+1)-OB.obs(isat,MAXFREQ+2))/(1.d0-(SAT(isat).freq(1)/SAT(isat).freq(2))**2)
          DO j=1, CKF.nfreq(isys)
            iondel(j)=(SAT(isat).freq(1)/SAT(isat).freq(j))**2*SIT.ion(isat)
          END DO
        END IF
      END IF
    ELSE
      iondel=0.d0
    END IF
    SIT.iondel(isat,:)=iondel(:)

    !! compute antenna / transmitter orientation dependent phase (cycle)
    !! corrections for Right Circularly Polarized electro magnetic waves
    IF (SIT.skd(1:1) .NE. 'D') THEN
      DO j=1, CKF.nfreq(isys)
        CALL phase_windup(.TRUE.,SIT.first(isat),rot_f2j,SIT.rot_l2f,SAT(isat).xscf,SAT(isat).yscf,&
                                                 SAT(isat).zscf,r1(1,j),SIT.prephi(isat),dphwp(j))
      END DO
    ELSE
      CALL phase_windup(.FALSE.,SIT.first(isat),rot_s2j,SAT(ileo).rant,SAT(isat).xscf,SAT(isat).yscf,&
                                                 SAT(isat).zscf,r1(1:3,1),SIT.prephi(isat),dphwp(1))
      dphwp(2:CKF.nfreq(isys))=dphwp(1)
    END IF

    !@ CMT BY XSY: FOR LEO ONBOARD CLOCK,THE RELETIVISTIC TIME DELAY IS INCLUDED
    ! IF (CKF%cprn(isat)(1:1).NE.'L') THEN
      !! general relativistic time delay due to the Earth gravity (meter)
      DO j=1, CKF.nfreq(isys)
        reldel(j)=2.d0*dot(3,xsat_j(1,j),xsat_j(4,j))/VEL_LIGHT
        ! reldel(j)=reldel(j)-2.d0*dot(3,xant_j(1,j,isys),xant_j(4,j,isys))/VEL_LIGHT
      END DO
    ! ELSE
    !   reldel=0.d0
    ! END IF

    !! gravitional change effects on the satellite oscillators.
    DO j=1, CKF.nfreq(isys)
      sitrad(j)=dsqrt(dot(3,xant_j(1,j,isys),xant_j(1,j,isys)))
      satrad(j)=dsqrt(dot(3,xsat_j(1,j),xsat_j(1,j)))
      reldel(j)=reldel(j)+2.d0*GME/(VEL_LIGHT)**2*LOG((sitrad(j)+satrad(j)+r1leng(j))/(sitrad(j)+satrad(j)-r1leng(j)))
    END DO      

    !! pcv correction for satellite and receiver antenna (meter)
    !! dpcv only contain the PCV for satellites
    CALL get_ant_pcv(SIT.iptatx,SAT(isat).iptatx,PI/2.d0-OB.elev(isat),OB.azim(isat),nadir,pcv,dpcv)
    
    OB.delay(isat)=(r1leng(1)+trpdel+reldel(1)+pcv(1,isys)+iondel(1)+SIT.sisre(isat))/VEL_LIGHT+dphwp(1)/SAT(isat).freq(1)
    DO j=1, CKF.nfreq(isys)
      !range1(j)=r1leng(j)+VEL_LIGHT*(drecclk(isys)-SAT(isat).sclock)+trpdel+reldel(j)+iondel(j)+pcv(j,isys)
      range1(j)=r1leng(j)+VEL_LIGHT*(drecclk(CKF.iref)-SAT(isat).sclock)+trpdel+reldel(j)+iondel(j)+pcv(j,isys)+SIT%sisre(isat)
      phase1(j)=range1(j)+dphwp(j)*VEL_LIGHT/SAT(isat).freq(j)-2.d0*iondel(j)
    END DO

    !! Finally, form observed minus calculated, in meter
    DO j=1, CKF.nfreq(isys)
      IF (OB.obs(isat,MAXFREQ+j).EQ.0.D0) CYCLE
      SIT%freqChannel(isys,j) = OB%fob(isat,j)

      OB.omc(isat,MAXFREQ+j)=OB.obs(isat,MAXFREQ+j)-range1(j)
      OB.omc(isat,j)=OB.obs(isat,j)*VEL_LIGHT/SAT(isat).freq(j)-phase1(j)
      
      !@ CMT BY XSY: GPS L5 AND BDS-2 B2
      IF (CKF%liar .AND. CKF%ArBiasMode .EQ. 'FCB') THEN
        IF ((CKF%cprn(isat)(1:1).EQ.'G' .AND. CKF%freq(j,isys).EQ.'L5') .OR. (INDEX(TRIM(SAT(isat)%type),'BEIDOU-2').NE.0 .AND. CKF%freq(j,isys).EQ.'L7')) THEN
          OB.omc(isat,j)=OB.omc(isat,j)+(SAT(isat).freq(1)**2/SAT(isat).freq(j)**2-1)*ifcb(isat) !ifcb in meter
        END IF
      END IF
      !@ CMT BY XSY: FOR LEO BIAS
      ! IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
      !   OB.omc(isat,MAXFREQ+j)=OB.omc(isat,MAXFREQ+j)-680.D0
      !   OB.omc(isat,j)=OB.omc(isat,j)-680.D0  
      ! END IF
    END DO

    !! correction the code bias for BeiDou observations
    !! for LEO satellites (at least FY3C), the systematic errors are not observed
    IF (INDEX(SAT(isat).type,'BEIDOU-2').NE.0 .AND. SIT.skd(1:1).NE.'D') THEN
      isys=INDEX(SYS,'C')
      cbias=0.d0
      CALL bds_code_cor(SAT(isat).type,CKF.nfreq(isys),CKF.freq(:,isys),OB.elev(isat)*RAD2DEG,cbias)
      DO j=1, CKF.nfreq(isys)
        OB.omc(isat,MAXFREQ+j)=OB.omc(isat,MAXFREQ+j)+cbias(j)
      END DO
    END IF

    !! weight of observations
    DO j=1, CKF.nfreq(isys)
      !@ CMT BY XSY: NLOS DOWN WEIGHT
      IF (OB.nlosflag(isat).EQ.0)THEN
        OB.var(isat,j)=(SIT.sigp(isys)*CKF.nlosfac)**2
        OB.var(isat,MAXFREQ+j)=(SIT.sigr(isys)*CKF.nlosfac)**2
      ELSE
        OB.var(isat,j)=(SIT.sigp(isys))**2
        OB.var(isat,MAXFREQ+j)= (SIT.sigr(isys))**2
      END IF
      !@ CMT BY XSY: FOR THE THIRD FREQUENCY OBSERVATIONS, DOWN THE WEIGHT
      ! IF (j.GE.3) THEN
      !   OB.var(isat,j)=OB.var(isat,j)/0.9d0
      !   OB.var(isat,MAXFREQ+j)=OB.var(isat,MAXFREQ+j)/0.9d0
      ! END IF
      !@ CMT BY XSY: DOWN THE WIGHT FOR THE POINYER SATELLITES
      ! IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
      !   OB.var(isat,j)=OB.var(isat,j)/0.9d0
      !   OB.var(isat,MAXFREQ+j)=OB.var(isat,MAXFREQ+j)/0.9d0        
      ! END IF     
    END DO

    IF (OB.elev(isat)*RAD2DEG .LE. 30.d0) THEN
      scal = 2.d0*dsin(OB.elev(isat))
      DO i=1,2*MAXFREQ
        OB.var(isat,i)=OB.var(isat,i)/scal**2
      END DO
    END IF

    ! IF (flag(isat).EQ..FALSE.) OB.var(isat,1:2*MAXFREQ)=OB.var(isat,1:2*MAXFREQ)*1.d10
    ! IF (flag(isat).EQ..FALSE. .AND. TRIM(CKF.uobs).EQ.'CODE') OB.omc(isat,1:2*MAXFREQ)=0.d0
    IF (TRIM(SAT(isat).type).EQ.'BEIDOU-2G' .OR. INDEX(SAT(isat).type,'BEIDOU-3G').NE.0) OB.var(isat,1:2*MAXFREQ)=OB.var(isat,1:2*MAXFREQ)*4

    !! Compute the delay rate and the predicted delay
    ! DO i=1,3
    !   dump(i)=xsat_j(i+3,1)-xant_j(i+3,1,isys)
    ! END DO
    ! drate=dot(3,dump,r1(1,1))/(VEL_LIGHT*r1leng(1))
    ! OB.delay(isat)=OB.delay(isat)+drate*CKF.dintv

    !! observation equation
    CALL ppp_partial(OB.npar,dloudx,drate,OB.pname,rot_f2j,ztdpart,SIT.grd,lnpar,lpart,OB.amat(1,isat))
  !! next satellite
  END DO

  !! offset of receiver clock
  IF (ite .LE. 10) THEN

    lcont = .FALSE.

    DO l=1, CKF.nsys

      isys=INDEX(SYS,CKF.system(l:l))

      nx=0
      DO isat=1, CKF.nprn
        IF (OB.omc(isat,MAXFREQ+1).NE.0.d0 .AND. flag(isat) .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(l:l)  &
            .AND. OB.elev(isat).GE.SIT.cutoff) THEN

          !@ CMT BY XSY: NLOS卫星会拉坏omc均值，影响判断，所以不参与接收机钟差迭代计算，对于没有经过SM算法探测NLOS的观测文件，不需要这一行
          ! IF (OB.nlosflag(isat) .EQ. 0) CYCLE

          nx=nx+1
          IF (CKF.nfreq(isys) .GE. 2) THEN
            rx(nx)=OB.omc(isat,MAXFREQ+1)*SAT(isat).fac(1)-OB.omc(isat,MAXFREQ+2)*SAT(isat).fac(2)
          ELSE
            !! contain the ionosphere
            rx(nx)=OB.omc(isat,MAXFREQ+1)
          END IF
          fg(nx)=0
          wx(nx)=1.d0
          itg(nx)=isat
        END IF
      END DO

      IF (nx .EQ. 0) CYCLE

      IF (nx .GT. 0) THEN
        CALL get_wgt_mean(.TRUE.,rx,fg,wx,nx,k,mean,rms,sig)
        ! IF (sig.GT.1.0D0) THEN
        !   thred(1)=4.d0
        !   thred(2)=4.d0
        ! ELSE
        !   thred(1)=30.d0
        !   thred(2)=10.d0
        ! END IF     
        DO WHILE(nx-k.GT.2 .AND. rms.GT.30.d0)
          j=k
          CALL sign_robust(nx,rx,fg,10.d0,k)
          IF (k .EQ. j) EXIT
          CALL get_wgt_mean(.FALSE.,rx,fg,wx,nx,k,mean,rms,sig)
        END DO
        IF (k .GT. 0) THEN
          DO i=1, nx
            IF(fg(i) .NE. 0) THEN
              flag(itg(i))=.FALSE.
              OB.omc(itg(i),1:2*MAXFREQ)=0.d0
              !! for kinematic position, this will remove too much observations, as the
              !! a prior station is not good enough
              ! IF (SIT.skd(1:2).EQ.'DK' .OR. SIT.skd(1:2).EQ.'K' .OR. SIT.skd(1:2).EQ.'DE') THEN
              !   IF (ite .GT. 5) THEN
              !     flag(itg(i))=.FALSE.
              !     OB.omc(itg(i),1:2*MAXFREQ)=0.d0
              !   END IF
              ! ELSE
              !   flag(itg(i))=.FALSE.
              !   OB.omc(itg(i),1:2*MAXFREQ)=0.d0
              ! END IF              
              IF (fg(i) .EQ. 2) THEN
                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4,I3)') '*** REMOVE(range_wgt ) ',&
                     CKF.mjd,CKF.sod,SIT.name,CKF.cprn(itg(i)),rx(i)-mean,mean,rms,OB.nlosflag(itg(i))
              ELSE
                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4,I3)') '*** REMOVE(range_sign) ',&
                     CKF.mjd,CKF.sod,SIT.name,CKF.cprn(itg(i)),rx(i)-mean,mean,rms,OB.nlosflag(itg(i))
              END IF
            END IF
          END DO
        END IF
        IF (isys .EQ. CKF.iref) THEN
          IF (dabs(mean/VEL_LIGHT).GT.1.d-7 .OR. SIT.rclock(isys).EQ.0.d0) THEN
            lcont=.TRUE.
            SIT.rclock(isys)=SIT.rclock(isys)+mean
            drecclk(isys)=SIT.rclock(isys)/VEL_LIGHT
          END IF
        ELSE
          ! ISB
          SIT.rclock(isys)=mean
        END IF
      END IF

    END DO

    IF (lcont .EQ. .TRUE.) GOTO 100

  END IF

  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)
  IF (CKF.sod .GT. sod0) THEN

      WRITE(1010,'(A3,1X,I4,4I3,F5.1,2X,(A))')'TIM',iy,imon,id,ih,im,sec,'--PHASE------------CODE---------DOPPLER-------------SNR-------IONpri-------TROPpri&
            --------ELEV--------AZIM---SlipFlag-----MultP1--------MultP2--LMC+2*GIM[m]-------IFphase--------IFcode---------IFlmc---------P1-P2------R1length'

      DO isat=1,CKF.nprn
         IF(OB.obs(isat,1).EQ.0.d0 .OR. OB.obs(isat,MAXFREQ+1).EQ.0.d0) CYCLE
         isys = INDEX(SYS,CKF.cprn(isat)(1:1))

         IF(CKF.nfreq(isys) .GE. 2)THEN
            OB.mltperr(isat,1) = OB.obs(isat,MAXFREQ+1) - (SAT(isat).g2+1)/(SAT(isat).g2-1) * OB.obs(isat,1) * VEL_LIGHT/SAT(isat).freq(1) + &
                                                                       2.0/(SAT(isat).g2-1) * OB.obs(isat,2) * VEL_LIGHT/SAT(isat).freq(2)
            OB.mltperr(isat,2) = OB.obs(isat,MAXFREQ+2) + (SAT(isat).g2+1)/(SAT(isat).g2-1) * OB.obs(isat,2) * VEL_LIGHT/SAT(isat).freq(2) - &
                                                          SAT(isat).g2*2.0/(SAT(isat).g2-1) * OB.obs(isat,1) * VEL_LIGHT/SAT(isat).freq(1)
         END IF

         IFphase = OB.obs(isat,1)*SAT(isat).lamda(1)*SAT(isat).fac(1) - OB.obs(isat,2)*SAT(isat).lamda(2)*SAT(isat).fac(2)
         IFcode  = OB.obs(isat,1+MAXFREQ)           *SAT(isat).fac(1) - OB.obs(isat,2+MAXFREQ)           *SAT(isat).fac(2)
         conion=''
         IF (sion(isat).NE.0.d0) conion='Y'
         
         DO i=1,CKF.nfreq(isys)
            IF (OB%obs(isat,i).EQ.0.D0 .OR. OB%obs(isat,MAXFREQ+i).EQ.0.D0) CYCLE

            lmc = VEL_LIGHT/SAT(isat).freq(i)*OB.obs(isat,i)-OB.obs(isat,i+MAXFREQ)+2*SIT.iondel(isat,i)
            p1p2 = OB.obs(isat,1+MAXFREQ)-OB.obs(isat,2+MAXFREQ)
            WRITE(1010,'(1X,A4,2X,A3,2X,I3,2X,A3,4(F14.3,2X),1X,F10.3,A2,2X,F10.3,2X,F10.3,2X,F10.3,6X,I2,7F14.3,F14.3)')SIT.name,CKF.cprn(isat),i,OB.fob(isat,i), &
                        OB.obs(isat,i),OB.obs(isat,i+MAXFREQ),OB.obs(isat,i+2*MAXFREQ),OB.obs(isat,i+3*MAXFREQ), &
                        SIT.iondel(isat,i),conion,SIT%trpdel(isat),OB.elev(isat)*RAD2DEG,OB.azim(isat)*RAD2DEG,OB.flag(isat,1),OB.mltperr(isat,1),OB.mltperr(isat,2), &
                        lmc,IFphase,IFcode,IFcode-IFphase,p1p2,Rsat2sit(isat)
         END DO
      END DO
  END IF
  mjd0=CKF.mjd
  sod0=CKF.sod

  RETURN

END SUBROUTINE
