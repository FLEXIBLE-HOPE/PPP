!*
SUBROUTINE fcb_gps_model(CKF,SIT,OB,SAT,BHD,IOD,ifcb)
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
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT
TYPE(RNXOBS) :: OB
TYPE(CKDCFG) :: CKF
TYPE(BRDHEAD) :: BHD
TYPE(IONEX) :: IOD
REAL(RL) :: ifcb(MAXSAT)

  !*
  ! The local variables
  !!-----------------------
  LOGICAL(LG) :: flag(MAXSAT),lcont,lfind
  INTEGER(IT) :: i,j,k,ite,isat,jd_send,jd_recv,ierr,nx,fg(MAXSAT),itg(MAXSAT),isys,l
  REAL(RL) :: sod_send,sod_recv,drecclk(MAXSYS),dsatclk,mean,rms,sig
  REAL(RL) :: fjd,drate,tdelay,ddelay,dmap,ah,aw
  REAL(RL) :: sitrad(MAXFREQ),satrad(MAXFREQ),reldel(MAXFREQ),nadir(MAXFREQ)
  REAL(RL) :: r1(3,MAXFREQ),r2(3),r1leng(MAXFREQ),phase1(MAXFREQ),range1(MAXFREQ)
  REAL(RL) :: rx(MAXSAT),wx(MAXSAT),xsun(6),xlun(6),dx(3),iondel(MAXSYS),dpcv(0:MAXPCVDEG,MAXFREQ,MAXSYS)
  REAL(RL) :: gmst,xpole,ypole,trpdel,ztdpart,dphwp(MAXFREQ),pcv(MAXFREQ,MAXSYS)
  REAL(RL) :: xant_f(6,MAXFREQ,MAXSYS),xant_j(6,MAXFREQ,MAXSYS),xsat_j(6,MAXFREQ),dump(3),dloudx(3)
  REAL(RL) :: scal,utcut1r,xhelp(2),rot_f2j(3,3),rot_rat(3,3),dxmat(3,3),dymat(3,3),sodut1
  INTEGER(IT) :: mjdut1,npwc(MAXICS),lnpwc(MAXICS),ileo,lnpar
  REAL(RL) :: lpart(MAXPWC*MAXPWCEQUS),dlpart(MAXPWC*MAXPWCEQUS)
  REAL(RL) :: qm(4),rot_s2j(3,3)
  REAL(RL) :: cbias(MAXFREQ),sion(MAXSAT)

  !*
  ! function called
  !!---------------------
  REAL(RL) :: dot, wetpp

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
    IF (SIT.skd(1:2) .EQ.'DD') THEN
      lnpar=SAT(ileo).npar
    ELSE
      lnpar=0
    END IF
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
            rot_f2j,SIT.rot_l2f,SIT.geod(1),SIT.geod(2),xpole,ypole,SIT.olc,dx) ! in km
        DO i=1,3
          xant_j(i,j,isys)=xant_j(i,j,isys)+dx(i)
        END DO
      END DO
    ! LEO [DP]
    ELSE
      
      CALL leo_interpolate_orbit(CKF.flnleo,jd_recv,sod_recv,SIT.ileo,SAT(ileo).cprn,SAT(ileo).npar, &
                     SAT(ileo).pname,.TRUE.,.TRUE.,.FALSE.,xsat_j(1,1),xsat_j(4,1),lnpwc,lpart,dlpart)
      xsat_j(1:6,1) = xsat_j(1:6,1)*1.D3

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
      SIT.x=0.d0
      DO i=1,3
        DO j=1,3
           SIT.x(i)=SIT.x(i)+rot_f2j(j,i)*xsat_j(j,1)
         END DO
      END DO
        
    END IF

  END DO

  !! loop over all satellites
  DO isat=1, CKF.nprn

    IF(OB.obs(isat,1).EQ.0.d0 .OR. OB.obs(isat,MAXFREQ+1).EQ.0.d0 .OR. .NOT.flag(isat)) CYCLE

    isys = INDEX(SYS,CKF.cprn(isat)(1:1))

    !### BLOCK 1
    !### ITERATION OF SEND TIME and GET THE GEOMETRIC DISTANCE
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
    !### COMPUTE CORRECTON

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
      OB.zmap(isat)=ztdpart

      SIT.grd(1) = -(SIT.zwd+SIT.ztdcor)/DTAN(OB.elev(isat))*DCOS(OB.azim(isat))   ! north
      SIT.grd(2) = -(SIT.zwd+SIT.ztdcor)/DTAN(OB.elev(isat))*DSIN(OB.azim(isat))   ! east
    ELSE
      trpdel=0.d0
    END IF

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
      CALL read_ion(CKF.mjd,CKF.sod,CKF,sion)
      IF (sion(isat) .NE. 0.d0) THEN
        SIT.ion(isat)=sion(isat)
        DO j=1, CKF.nfreq(isys)
          iondel(j)=SAT(isat).freq(1)**2/SAT(isat).freq(j)**2*SIT.ion(isat)
        END DO
      ELSE
        iondel=0.d0
      END IF
    ELSE
      iondel=0.d0
    END IF

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

    !! general relativistic time delay due to the Earth gravity (meter)
    DO j=1, CKF.nfreq(isys)
      reldel(j)=2.d0*dot(3,xsat_j(1,j),xsat_j(4,j))/VEL_LIGHT
      ! reldel(j)=reldel(j)-2.d0*dot(3,xant_j(1,j,isys),xant_j(4,j,isys))/VEL_LIGHT
    END DO

    !! gravitional change effects on the satellite oscillators.
    DO j=1, CKF.nfreq(isys)
      sitrad(j)=dsqrt(dot(3,xant_j(1,j,isys),xant_j(1,j,isys)))
      satrad(j)=dsqrt(dot(3,xsat_j(1,j),xsat_j(1,j)))
      reldel(j)=reldel(j)+2.d0*GME/(VEL_LIGHT)**2*LOG((sitrad(j)+satrad(j)+r1leng(j))/(sitrad(j)+satrad(j)-r1leng(j)))
    END DO

    !! pcv correction for satellite and receiver antenna (meter)
    CALL get_ant_pcv(SIT.iptatx,SAT(isat).iptatx,PI/2.d0-OB.elev(isat),OB.azim(isat),nadir,pcv,dpcv)

    OB.delay(isat)=(r1leng(1)+trpdel+reldel(1)+pcv(1,isys)+iondel(1))/VEL_LIGHT+dphwp(1)/SAT(isat).freq(1)
    DO j=1, CKF.nfreq(isys)
      !range1(j)=r1leng(j)+VEL_LIGHT*(drecclk(isys)-SAT(isat).sclock)+trpdel+reldel(j)+iondel(j)+pcv(j,isys)
      range1(j)=r1leng(j)+VEL_LIGHT*(drecclk(CKF.iref)-SAT(isat).sclock)+trpdel+reldel(j)+iondel(j)+pcv(j,isys)
      phase1(j)=range1(j)+dphwp(j)*VEL_LIGHT/SAT(isat).freq(j)-2.d0*iondel(j)
    END DO

    !! Finally, form observed minus calculated, in meter
    DO j=1, CKF.nfreq(isys)
      IF (OB.obs(isat,MAXFREQ+j).EQ.0.D0) CYCLE

      OB.omc(isat,j)=OB.obs(isat,j)*VEL_LIGHT/SAT(isat).freq(j)-phase1(j)
      OB.omc(isat,MAXFREQ+j)=OB.obs(isat,MAXFREQ+j)-range1(j)
      
      !@ CMT BY XSY: GPS L5 AND BDS-2 B2
      IF ((CKF%cprn(isat)(1:1).EQ.'G' .AND. CKF%freq(j,isys).EQ.'L5') .OR. (INDEX(TRIM(SAT(isat)%type),'BEIDOU-2').NE.0 .AND. CKF%freq(j,isys).EQ.'L7')) THEN
        OB.omc(isat,j)=OB.omc(isat,j)+(SAT(isat).freq(1)**2/SAT(isat).freq(3)**2-1)*ifcb(isat) !ifcb in meter
      END IF      
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
      OB.var(isat,j)=(SIT.sigp(isys))**2
      OB.var(isat,MAXFREQ+j)= SIT.sigr(isys)**2
    END DO
    IF (OB.elev(isat)*RAD2DEG .LE. 30.d0) THEN
      scal = 2.d0*dsin(OB.elev(isat))
      DO i=1,2*MAXFREQ
        OB.var(isat,i)=OB.var(isat,i)/scal**2
      END DO
    END IF
    IF (TRIM(SAT(isat).type).EQ.'BEIDOU-2G' .OR. INDEX(SAT(isat).type,'BEIDOU-3G').NE.0) OB.var(isat,1:2*MAXFREQ)=OB.var(isat,1:2*MAXFREQ)*4

    !! Compute the delay rate and the predicted delay
    !DO i=1,3
    !  dump(i)=xsat_j(i+3,1)-xant_j(i+3,1,isys)
    !END DO
    ! drate=dot(3,dump,r1(1,1))/(VEL_LIGHT*r1leng(1))
    ! OB.delay(isat)=OB.delay(isat)+drate*CKF.dintv

    !! observation equation
    CALL fcb_partial(OB.npar,dloudx,drate,OB.pname,OB.ltog(1,isat),rot_f2j,ztdpart,SIT.grd,OB.amat(1,isat))

  !! next satellite
  END DO

  !! offset of receiver clock
  IF (ite .LE. 10) THEN

    lcont=.FALSE.

    DO l=1, CKF.nsys

      isys=INDEX(SYS,CKF.system(l:l))

      nx=0
      DO isat=1, CKF.nprn
        IF (OB.omc(isat,MAXFREQ+1).NE.0.d0 .AND. flag(isat) .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(l:l)  &
            .AND. OB.elev(isat).GE.SIT.cutoff) THEN
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
              IF (fg(i) .EQ. 2) THEN
                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4)') 'REMOVE(range_wgt) ',&
                     CKF.mjd,CKF.sod,SIT.name,CKF.cprn(itg(i)),rx(i)-mean,mean,rms
              ELSE
                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4)') 'REMOVE(range_sign) ',&
                     CKF.mjd,CKF.sod,SIT.name,CKF.cprn(itg(i)),rx(i)-mean,mean,rms
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

  RETURN

END SUBROUTINE
