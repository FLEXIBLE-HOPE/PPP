!*
PROGRAM oi
!!
!*
USE orbit
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  !*
  ! The local variables
  !!---------------------------
  TYPE(ORBCFG) :: CKF
  TYPE(SATE) :: SAT(MAXSAT)

  INTEGER(IT) :: i,j,k,kk,nsign
  INTEGER(IT) :: isat,nepo,iequ
  INTEGER(IT) :: idir,ilast,iepo
  INTEGER(IT) :: lfnt,lfnsat,lfn

  REAL(RL) :: sod0,sod1
  REAL(RL) :: tb,te,t,tdir,hh
  REAL(RL) :: x(MAXEQUS),acc(MAXEQUS)
  REAL(RL) :: amat(3,3),bmat(3,3),cmat(3*MAXICS)
  REAL(RL) :: funct(0:20,MAXEQUS),fright(0:20,MAXEQUS)

  INTEGER(IT) :: nman,iman,ii
  LOGICAL(LG) :: lman
  REAL(RL) :: tman(2,MAXMAN),t0,t1,samp,hs
  LOGICAL(LG) :: loutref,lin,lrkf

  CHARACTER(LEN=2) :: cdir

  !*
  ! The function called
  !!---------------------------
  CHARACTER(LEN_RUNTIM) :: run_tim
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  CALL get_oi_args(CKF,SAT)

  lfnt=get_valid_unit(500)

  DO isat=1, CKF.nprn

    WRITE(OUTPUT_UNIT,'(A,A)') run_tim(), 'Integrating for satellite '//CKF.cprn(isat)

    CALL oi_fright_init(CKF,SAT(isat),tman,nman)

    ! decide if forward and/or backward integeration are needed
    sod0=(CKF.mjd0-CKF.rmjd)*86400.d0+CKF.sod0
    sod1=(CKF.mjd1-CKF.rmjd)*86400.d0+CKF.sod1
    loutref=.FALSE.

    IF (sod0.LE.CKF.rsod .AND. sod1.GE.CKF.rsod) loutref=.TRUE.

    ! Be aware, backward must be before forward as barckward result has to be
    ! reversed and put to the scarth file at first.
    cdir=''
    tb=sod0
    te=sod1
    IF (tb .GT. CKF.rsod) tb=CKF.rsod
    IF (te .GT. CKF.rsod) cdir='F'
    IF (tb .LT. CKF.rsod) cdir='B'//TRIM(cdir)

    lfnsat=lfnt+isat
    OPEN(UNIT=lfnsat,STATUS='SCRATCH',FORM='UNFORMATTED')

    ! maneuver
    lman=.FALSE.
    IF (nman .GT. 0) THEN
      lman=.TRUE.
      DO iman=1, nman
        tman(1,iman)=(tman(1,iman)-CKF.rmjd)*86400.d0
        tman(2,iman)=(tman(2,iman)-CKF.rmjd)*86400.d0
      END DO
    END IF

    nepo=0
    DO idir=1, LEN_TRIM(cdir)

      ! For each satellite one scratch file is opened for the whole arc. for
      ! backward, the result is in a reversed order in time, so we put it in
      ! another scratch file (unit 3), after integration finished we read it
      ! from end to the begin and put it into the scratch file for each sat.
      ! For forward result is put into the same scratch file directly. That
      ! why backward integration must be performed firstly
      IF (cdir(idir:idir) .EQ. 'B') THEN
        lfn=get_valid_unit(10)
        OPEN(UNIT=lfn,STATUS='SCRATCH',FORM='UNFORMATTED')
      ELSE
        lfn=lfnsat
      END IF

      ! The right start and stop time and step sign for backward and forward
      ! tidir1 is the destination
      IF (cdir(idir:idir) .EQ. 'B' ) THEN
        tdir=tb
        nsign=-1
      ELSE
        tdir=te
        nsign=1
      END IF

      IF (idir .EQ. 1) THEN
        CALL oi_fright_initorb(CKF.nequ,x,lfnsat)
      ELSE
        CALL oi_fright_initorb(CKF.nequ,x,0)
      END IF

      ! single step RKF, RKF step is 1/6 of ADAMS step
      t=CKF.rsod
      hh=-CKF.int_step*nsign/30
      DO i=CKF.norder_adams,1,-1
        IF (i .NE. CKF.norder_adams) THEN
          lin=.FALSE.
          IF (lman .EQ. .TRUE.) THEN
            t1=t-CKF.int_step*nsign
            DO iman=1, nman
              IF (nsign .EQ. -1) THEN
                IF (t1.GE.tman(1,iman) .AND. t.LE.tman(2,iman) .AND.(t.LT.tman(1,iman).OR.t1.GT.tman(2,iman))) lin=.TRUE.
              ELSE
                IF (t1.LE.tman(2,iman) .AND. t.GE.tman(1,iman) .AND.(t1.LT.tman(1,iman).OR.t.GT.tman(2,iman))) lin=.TRUE.
              END IF
              IF (lin .EQ. .TRUE.) exit
            END DO
          END IF
          ii=5
          IF (lin .EQ. .TRUE.) THEN
            samp=MIN((tman(2,iman)-tman(1,iman))*0.1d0,0.01d0)
            DO WHILE(CKF.int_step/DBLE(ii).GT.samp)
              ii=ii+1
            END DO
          END IF
          hh=-CKF.int_step*nsign/ii
          DO k=1, ii ! 30
            CALL oi_rkf(CKF.rmjd,t,x,hh,x,CKF.nequ,acc,amat,bmat,cmat)
            t=t+hh
          END DO
        END IF
        CALL oi_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)

        ! Save funct and fright for multi-step method in reverse sequence
        DO j=1, CKF.nequ
          funct(i-1,j)=x(j)
        END DO
        DO j=1, CKF.nequ/6
          kk=(j-1)*6
          DO k=1, 3
            fright(i-1,kk+k)=x(kk+k+3)
            fright(i-1,kk+3+k)=acc((j-1)*3+k)
          END DO
        END DO
      END DO

      ! ADAMS, iepo is integrating epoch
      ilast=CKF.norder_adams
      IF (idir.EQ.1 .AND. loutref) THEN
        WRITE(lfn) CKF.rsod,((funct(ilast-1,(k-1)*6+j),j=1,6),k=1,CKF.nequ/6)
        nepo=nepo+1
      END IF

      iepo=0
      hh=CKF.int_step*nsign
      DO WHILE (t*nsign .LT. tdir*nsign)
        iepo=iepo+1
        t=CKF.rsod+iepo*hh

        lrkf=.FALSE.
        IF (lman .EQ. .TRUE.) THEN
          DO iman=1, nman
            IF (nsign .EQ. 1) THEN
              IF (t.GE.tman(1,iman).AND.t.LE.tman(2,iman)+ilast*CKF.int_step) lrkf=.TRUE.
            ELSE
              IF (t.le.tman(2,iman).and.t.ge.tman(1,iman)-ilast*CKF.int_step) lrkf=.TRUE.
            END IF
            IF (lrkf .EQ. .TRUE.) exit
          END DO
        END IF

        IF (lrkf .EQ. .TRUE.) THEN
          lin=.FALSE.
          t0=t-hh
          DO iman=1, nman
            IF (nsign .EQ. -1) THEN
              IF (t0.GE.tman(1,iman).AND.t.LE.tman(2,iman).AND.(t.LT.tman(1,iman).OR.t0.GT.tman(2,iman))) lin=.TRUE.
            ELSE
              IF (t0.LE.tman(2,iman).AND.t.GE.tman(1,iman).AND.(t0.LT.tman(1,iman).OR.t.GT.tman(2,iman))) lin=.TRUE.
            END IF
            IF (lin .EQ. .TRUE.) exit
          END DO
          ii=5
          IF (lin .EQ. .TRUE.) THEN
            samp=MIN((tman(2,iman)-tman(1,iman))*0.1d0,0.01d0)
            DO WHILE(CKF.int_step/DBLE(ii) .GT. samp)
              ii=ii+1
            END DO
          END IF

          hs = -CKF.int_step*nsign/ii
          t = t0
          DO k=1, ii
            CALL oi_rkf(CKF.rmjd,t,x,-hs,x,CKF.nequ,acc,amat,bmat,cmat)
            t=t-hs
          ENDDO
          CALL oi_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)

          DO j=1, CKF.nequ
            funct(ilast,j)=x(j)
          END DO
          DO j=1, CKF.nequ/6
            kk=(j-1)*6
            DO k=1, 3
              fright(ilast,kk+k)=funct(ilast,kk+k+3)
              fright(ilast,kk+3+k)=acc((j-1)*3+k)
            END DO
          END DO

        ELSE

          ! Adams prediction for iepo
          DO iequ=1, CKF.nequ
            CALL oi_adams('PRE',CKF.norder_adams,funct(0,iequ),fright(0,iequ),hh)
          END DO

          ! f_right at iepo based on prediction
          DO j=1, CKF.nequ
            x(j)=funct(ilast,j)
          END DO
          CALL oi_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)
          DO j=1, CKF.nequ
            funct(ilast,j)=x(j)
          END DO
          DO j=1, CKF.nequ/6
            kk=(j-1)*6
            DO k=1, 3
              fright(ilast,kk+k)=funct(ilast,kk+k+3)
              fright(ilast,kk+3+k)=acc((j-1)*3+k)
            END DO
          END DO

          ! Adams correction
          DO iequ=1, CKF.nequ
            CALL oi_adams('COR',CKF.norder_adams,funct(0,iequ),fright(0,iequ),hh)
          END DO
        END IF

        ! Output to scratch file
        IF (t.GT.sod0-1.d-5 .AND. t.LT.sod1+1.d-5 .AND. &
             DABS(t-CKF.rsod-NINT((t-CKF.rsod)/(nsign*CKF.out_step))*nsign*CKF.out_step).LE.1d-5) THEN
          WRITE(lfn) t,((funct(ilast,(k-1)*6+j),j=1,6),k=1,CKF.nequ/6)
          nepo=nepo+1
        END IF

        ! f_right at iepo based on estimation
        DO j=1, CKF.nequ
          x(j)=funct(ilast,j)
        END DO
        CALL oi_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)

        ! Get back in case that vel. changed
        DO j=1, CKF.nequ
          funct(ilast,j)=x(j)
        END DO
        DO j=1, CKF.nequ/6
          kk=(j-1)*6
          DO k=1, 3
            fright(ilast,kk+k)=funct(ilast,kk+k+3)
            fright(ilast,kk+3+k)=acc((j-1)*3+k)
          END DO
        END DO

        ! Shift funct and fright for the next step
        DO i=0, CKF.norder_adams-1
          DO iequ=1, CKF.nequ
            funct(i,iequ)=funct(i+1,iequ)
            fright(i,iequ)=fright(i+1,iequ)
          END DO
        END DO
      END DO

      ! Reverse backward result and write to the scratch file for this satellite
      IF (cdir(idir:idir) .EQ. 'B') THEN
        DO WHILE(t .LT. CKF.rsod)
          BACKSPACE(lfn)
          READ(lfn) t,((funct(ilast,(k-1)*6+j),j=1,6),k=1,CKF.nequ/6)
          WRITE(lfnsat) t,((funct(ilast,(k-1)*6+j),j=1,6),k=1,CKF.nequ/6)
          BACKSPACE(lfn)
        ENDDO
        CLOSE(lfn)
      END IF

    END DO

    REWIND(lfnsat)

  END DO

  CALL oi_meg_flnorb(CKF,lfnt,nepo)

  CALL close_all_files()

  STOP

END PROGRAM
