!*
SUBROUTINE oi_fright_acc(mjd,sod,funct,acc,amat,bmat,cmat)
!!
!*
USE orbit
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(SATE) :: SAT0
TYPE(ORBCFG) :: CKF0
INTEGER(IT) :: mjd,nequ,lfn
REAL(RL) :: sod,x(1:*),funct(1:*)
REAL(RL) :: acc(1:*),amat(3,3)
REAL(RL) :: bmat(3,3),cmat(1:*)

  !*
  ! The local variables
  !!--------------------------
  TYPE(PLANET_INFO) :: PL
  TYPE(SATEPAR) :: ICS
  TYPE(SATE) :: SAT
  TYPE(SATEPAN) :: PAN
  TYPE(ORBCFG) :: CKF

  INTEGER(IT) :: i,j,k,ii,iequ
  INTEGER(IT) :: otdeg,optdeg,aodeg,ndegree

  REAL(RL) :: period,dt,rot_s2j(3,3),rot_l2c(3,3)
  REAL(RL) :: xics(MAXEQUS),mjdtt

  REAL(RL) :: acc1(3),acc2(3),lambda

  !! solid earth tides and pole tide
  REAL(RL) :: dce(4,0:4),dse(4,0:4)
  !! ocean tides
  REAL(RL) :: dco(MAXOCNDEG,0:MAXOCNDEG)
  REAL(RL) :: dso(MAXOCNDEG,0:MAXOCNDEG)
  !! ocean pole tides
  REAL(RL) :: dcop(MAXOPTDEG,0:MAXOPTDEG)
  REAL(RL) :: dsop(MAXOPTDEG,0:MAXOPTDEG)
  !! atmosphere and ocean dealiasing
  REAL(RL) :: dca(MAXAODDEG,0:MAXAODDEG)
  REAL(RL) :: dsa(MAXAODDEG,0:MAXAODDEG)

  INTEGER(IT) :: mjdutc,mjdtai
  REAL(RL) :: sodutc,sodtai
  INTEGER(IT) :: mjdut1 = 0
  REAL(RL) :: sodut1 = 0.d0
  REAL(RL) :: utcut1r = 0.D0
  REAL(RL) :: xhelp(2) = 0.D0
  REAL(RL) :: mate2j(3,3) = 0.D0
  REAL(RL) :: rmte2j(3,3) = 0.D0
  REAL(RL) :: dxmat(3,3) = 0.D0
  REAL(RL) :: dymat(3,3) = 0.D0
  REAL(RL) :: gmst = 0.D0
  REAL(RL) :: xpole = 0.D0
  REAL(RL) :: ypole = 0.D0
  REAL(RL) :: erp(3),q(4)

  INTEGER(IT) :: nman,nman0
  REAL(RL) :: tman(2,MAXMAN),tman0(2,MAXMAN)

  LOGICAL(LG) :: lfind

  SAVE CKF,SAT,PL,PAN,xics,period,ICS,tman,nman

  !*
  ! Start the exectuable code
  !!---------------------------
  
  otdeg=0
  optdeg=0
  aodeg=0
  ndegree=0

  dce=0.d0
  dse=0.d0

  dco=0.d0
  dso=0.d0

  dcop=0.d0
  dsop=0.d0

  dca=0.d0
  dsa=0.d0

  ! Initialize
  DO i=1, CKF.nequ
    acc(i)=0.d0
  END DO
  DO i=1, 3
    acc1(i)=0.d0
    acc2(i)=0.d0
    DO j=1, 3
      amat(i,j)=0.d0
      bmat(i,j)=0.d0
    END DO
  END DO

  ! find the right model parameter values
  k=0
  dt=mjd+sod/86400.d0
  DO i=1, SAT.npar

    j=0
    ! Before the start time of the first one, use the first one
    IF ((dt-ICS.ptime(1,k+1))*86400.d0 .LT. MAXWND) THEN
      ii=k+1

    ! After the end time of the last one, use the last one
    ELSE IF((dt-ICS.ptime(2,k+ICS.npwc(i)))*86400.d0 .GT. -MAXWND) THEN
      ii = k+ICS.npwc(i)

    ! Search for the right one
    ELSE
      DO j=1, ICS.npwc(i)
        IF (dt.GE.ICS.ptime(1,k+j) .AND. dt.LT.ICS.ptime(2,k+j) ) ii=k+j
      END DO
    END IF
    xics(i)=ICS.val(ii)
    k=k+ICS.npwc(i)

    ! For velocity (impulse)
    IF (j.NE.0 .AND. DABS((dt-ICS.ptime(1,ii))*86400.d0).LE.MAXWND .AND. i.GE.4 .AND. i.LE.6) THEN
      funct(i)=funct(i)+xics(i)
    END IF
  END DO

  CALL timinc(mjd,sod,OFF_GPS2TAI,mjdtai,sodtai)
  CALL taiutc(mjdtai,sodtai,mjdutc,sodutc)
  CALL table_linear_interpolate('poleut1',.FALSE.,mjdutc+sodutc/86400.d0,erp)
  xhelp=erp(1:2)
  utcut1r=erp(3)
  CALL itrs2gcrs('IERS2010',mjd,sod,utcut1r,xhelp,mate2j,rmte2j,dxmat,dymat,mjdut1,sodut1,gmst,xpole,ypole)

  DO i=1, 6
    PL.xj(i,PL.isc)=funct(i)
  END DO
  CALL oi_eph_planet(PL,mjd,sod,mate2j,rmte2j)

  !! Satellite attitude
  IF (SAT.cprn(1:1) .EQ. 'L') THEN
    !! no initial value for flnatt
    CALL read_sat_asc(SAT.cprn,SAT.flnatt,mjd+sod/86400.d0,q,lfind)
    IF (lfind .EQ. .TRUE.) THEN
      CALL q2rotmat(SAT.type,q,rot_s2j)
      IF (SAT.type(1:5) .EQ. 'SWARM') THEN
        CALL matmpy(mate2j,rot_s2j,rot_s2j,3,3,3)
      END IF
    ELSE
      !! maybe we need to define the attitude control mode (to be solved)
      CALL leo_att_non(SAT.type,PL.xj(1,PL.isc),rot_s2j)
    END IF
    !write(*,*) rot_s2j

    !! the solar panel direction
    DO i=1, PAN.npan
      CALL matmpy(rot_s2j,PAN.norm(1:3,i),PAN.normj(1:3,i),3,3,1)
    END DO
  END IF

  CALL scf2ecef(PL.xj(1,PL.isc),rot_l2c)

  mjdtt=mjd+(sod+OFF_GPS2TT)/86400.d0

  DO i=1, SAT.nforce

    IF (SAT.force_model(i)(1:5).EQ.'EMPTY' .OR. SAT.force_model(i)(1:4).EQ.'NONE') CYCLE
    SELECT CASE(TRIM(SAT.force_name(i)))
      CASE('Point mass')
        CALL oi_point_mass(CKF.lpart,PL,SAT.force_model(i),acc,amat)

      CASE('Solid Earth tides')
        CALL oi_solid_tides(PL,SAT.force_model(i),mjdtt,gmst,dce,dse)

      CASE('Ocean tides')
        j=INDEX(SAT.force_model(i),' ')
        READ(SAT.force_model(i)(j+1:),*) otdeg
        IF (otdeg .EQ. 0) CYCLE
        IF (INDEX(SAT.force_model(i),'EOT11a') .NE. 0) THEN
          CALL oi_eot11a_tide(gmst,mjdtt,otdeg,dco,dso)
        ELSE
          CALL oi_fes2014b_tide(gmst,mjdtt,otdeg,dco,dso)
        END IF
        

      CASE('Earth pole tide')
        CALL oi_etpole_tide(mjdtt,xpole,ypole,dce(2,1),dse(2,1))

      CASE('Ocean pole tide')
        READ(SAT.force_model(i),*) optdeg
        IF (optdeg .EQ. 0) CYCLE
        CALL oi_ocpole_tide(mjdtt,xpole,ypole,optdeg,dcop,dsop)

      CASE('ATM and ocean var.')
        READ(SAT.force_model(i),*) aodeg
        IF (aodeg .EQ. 0) CYCLE
        CALL oi_aod(mjdtt,aodeg,dca,dsa)

      CASE('Gravity model')
        READ(SAT.force_model(i),*) ndegree
        CALL oi_gravity_pines(CKF.lpart,mjdtt,PL.xe(1,PL.isc),mate2j, &
             ndegree,dce,dse,otdeg,dco,dso,optdeg,dcop,dsop,aodeg,dca,dsa, &
             CKF.lgfm,CKF.mindeg,CKF.maxdeg,CKF.ltog,SAT.npar-CKF.ngc,acc,amat,cmat)

      CASE('Relativity')
        CALL oi_relativity(PL.gm(PL.icb),PL.gm(PL.isun),PL.xj(1,PL.isc),PL.xj(1,PL.isun),acc)

      CASE('Solar radiation')
        CALL oi_shadow_factor(PL.radius(PL.isun),PL.radius(PL.icb),PL.radius(4),PL.xj(1,PL.isun),PL.xj(1,PL.isc),PL.xj(1,4),lambda)
        IF (SAT.force_model(i)(1:4) .EQ. 'BERN') THEN
          CALL oi_srp_bern(SAT.force_model(i),CKF.lpart,mjd,sod,SAT.cprn,SAT.csvn,SAT.type,PAN,SAT.mass,lambda,SAT.npar,SAT.pname,xics,&
               PL.xj(1,PL.isc),PL.xj(1,PL.isun),acc1,cmat)
        ELSE IF (SAT.force_model(i)(1:4) .EQ. 'BOXW') THEN
          CALL oi_srp_boxw(CKF.lpart,mjd,sod,SAT.cprn,SAT.csvn,SAT.type,PAN,SAT.mass,lambda,SAT.npar,SAT.pname,xics,PL.xj(1,PL.isc), &
               PL.xj(1,PL.isun),acc1,cmat)
        ELSE
          CALL oi_srp_gen(CKF.lpart,PAN,SAT.mass,lambda,SAT.npar,SAT.pname,xics,PL.xj(1,PL.isc), &
               PL.xj(1,PL.isun),acc1,cmat)
        END IF

      CASE('Earth radiation')
        CALL ERPFBOXW(SAT.type,SAT.cprn,SAT.csvn,SAT.mass,mjdtt,1,PL.xj(1,PL.isc),PL.xj(1,PL.isun),acc1)

      CASE('Atmosphere drag')
        CALL oi_atm_drag(SAT.force_model(i),CKF.lpart,mjdutc+sodutc/86400.d0,PAN,SAT.mass,SAT.npar,SAT.pname,xics, &
               mate2j,PL.xe(1,PL.isc),PL.xe(1,PL.isun),acc1,cmat)

      CASE('Accelerator obs')
        ! CALL oi_accelerator()

      CASE('Thrust force')
        CALL oi_thrust_acc(CKF.lpart,SAT.npar,SAT.pname,xics,nman,tman,rot_l2c,mjd+sod/86400.d0,cmat,acc1)

      CASE('Customer model')
        dt=(mjd-CKF.rmjd)*86400.d0+(sod-CKF.rsod)
        CALL oi_customer_acc(CKF.lpart,SAT.cprn,SAT.npar,SAT.pname,xics,period,dt,rot_l2c,cmat,acc1)

      CASE DEFAULT
        WRITE(OUTPUT_UNIT,'(A)') '###WARNING(oi_fright_acc): unknown force model '//TRIM(SAT.force_name(i))
    END SELECT
  END DO

  DO i=1, 3
    acc(i)=acc(i)+acc1(i)+acc2(i)
  END DO

  IF (.NOT. CKF.lpart) RETURN

  ! Rightside of variation equaiton:  AMat*partial(x)/partial(q)+BMat*partial(x)/partial(dq)
  ! iequ==1 is for motion equation, k is for components x,y,z
  DO iequ=2, 7
    i=iequ-1
    DO k=1, 3
      acc(i*3+k)=0.d0
      DO j=1, 3
        ii=6*i+j
        acc(i*3+k)=acc(i*3+k)+amat(k,j)*funct(ii)+bmat(k,j)*funct(ii+3)
      END DO
    END DO
  END DO

  ! Rightside of variation equaiton: AMat*partial(x)/partial(q)+BMat*partial(x)/partial(dq)+CMat
  ! for force model parameters
  DO iequ=8, CKF.nequ/6
    i=iequ-1
    DO k=1, 3
      acc(i*3+k)=cmat((iequ-8)*3+k)
      DO j=1, 3
        ii=6*i+j
        acc(i*3+k)=acc(i*3+k)+amat(k,j)*funct(ii)+bmat(k,j)*funct(ii+3)
      END DO
    END DO
  END DO

  RETURN

ENTRY oi_fright_init(CKF0,SAT0,TMAN0,NMAN0)

  CKF=CKF0
  SAT=SAT0

  CALL read_ics(CKF.flnics,SAT,ICS)

  k=0
  DO i=1, SAT.npar
    xics(i)=ICS.val(k+1)
    k=k+ICS.npwc(i)
  END DO

  PAN.cprn=SAT.cprn
  DO i=1, SAT.nforce
    IF (INDEX(SAT.force_model(i),'PANEL') .NE. 0) THEN
      !! need to be considerated for LEO satellites
      IF (PAN.npan .EQ. 0) CALL oi_sat_panel(SAT.type,PAN)
    END IF
  END DO

  ! get the maneuver information
  NMAN=0
  NMAN0=0
  TMAN=0.d0
  TMAN0=0.d0
  DO i=1, SAT.nforce
    IF (INDEX(SAT.force_name(i),'Thrust force') .EQ. 0) CYCLE
    IF (INDEX(SAT.force_model(i),'YES') .EQ. 0) CYCLE
    CALL read_thrust(CKF.mjd0+CKF.sod0/86400.d0,CKF.mjd1+CKF.sod1/86400.d0,SAT.type,tman)
    DO k=1, MAXMAN
      IF (TMAN(1,k) .EQ. 0.d0) CYCLE
      NMAN=NMAN+1
    END DO
  END DO
  TMAN0=TMAN
  NMAN0=NMAN

  CALL oi_planet_init(PL)
  WRITE(PL.name(1),'(A3)') SAT.cprn

  CALL eci2orb(SAT.cprn,xics(1:6),period)

  RETURN

ENTRY oi_fright_initorb(nequ,x,lfn)

  DO i=1, 6
    x(i)=xics(i)
  END DO

  k=6
  DO i=2, nequ/6
    DO j=1, 6
      k=k+1
      IF (i-1 .EQ. j) THEN
        x(k)=1.d0
      ELSE
        x(k)=0.d0
      END IF
    END DO
  END DO

  IF (lfn .NE. 0) THEN
    WRITE(lfn) ICS
  END IF

  RETURN

END SUBROUTINE
