!*
SUBROUTINE chsip(CKF,SAT,OB,OBH,AM,NM)
!!
!! Repair cycle slips by epoch difference
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------------
TYPE(SATE) :: SAT(MAXSAT)
TYPE(RNXOBS) :: OB,OBH
TYPE(CKDCFG) :: CKF
TYPE(AMBT) :: AM(1:*)
TYPE(INFM) :: NM

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i,ind,isat,ksat,iamb,nost,nsip,nfix,bsip(MAXSAT),isys,npar
  REAL(RL) :: wcls(MAXSAT),ncls(MAXSAT)
  REAL(RL) :: elev,pomc(2,MAXSAT),comc(2,MAXSAT)

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_int
  INTEGER(IT) :: pointer_string
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!----------------------------

  npar=0
  IF (pointer_string(OB.npar,OB.pname,'STAPX') .NE. 0) npar=3

  DO isys=1, CKF.nsys

    !! whether repair or not
    iamb=pointer_string(OB.npar,OB.pname,'AMBL1')
    DO isat=1,CKF.nprn
      IF (OB.omc(isat,1).NE.0.d0 .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(isys:isys)) THEN
        IF (OB.flag(isat,1).NE.0 .AND. OB.ltog(iamb,isat).NE.0) EXIT
      END IF
    END DO
    IF (isat .GT. CKF.nprn) CYCLE

    !! form epoch difference
    ! wide-lane cycle slip
    wcls=0.5d0
    ! L1/narrow-lane cycle slip
    ncls=0.5d0
    IF (timdif(CKF.mjd,CKF.sod,OBH.jd,OBH.tsec) .LE. 30.d0) THEN
    !IF (timdif(CKF.mjd,CKF.sod,OBH.jd,OBH.tsec) .LE. 5*CKF.dintv) THEN
      ! number of available satellites
      nost=0
      ! number of to-be-repaired cycle slips
      nsip=0
      ! reference satellite
      ksat=0
      ! maximum elevation
      elev=0.d0
      DO isat=1,CKF.nprn
        bsip(isat)=0       ! flag for slip status
        pomc(1:2,isat)=0.d0
        comc(1:2,isat)=0.d0

        IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

        IF (OB.omc(isat,1).NE.0.d0 .AND. OBH.omc(isat,1).NE.0.d0) THEN
          nost=nost+1
          comc(1,isat)=OB.omc(isat,1)-OBH.omc(isat,1) ! mitigate ionospheric delays
          comc(2,isat)=OB.omc(isat,2)-OBH.omc(isat,2)
          pomc(1,isat)=OB.omc(isat,MAXFREQ+1)-OBH.omc(isat,MAXFREQ+1)
          pomc(2,isat)=OB.omc(isat,MAXFREQ+2)-OBH.omc(isat,MAXFREQ+2)
          IF (OB.flag(isat,1) .NE. 0) THEN
            nsip=nsip+1
            bsip(isat)=1       ! slip set
          END IF
          IF (OB.elev(isat) .GT. elev) THEN
            ksat=isat
            elev=OB.elev(isat)
          END IF
        END IF
      END DO

      !! is it possible to repair
      nfix=0  ! number of successfully repaired cycle slips
      IF (nsip.GT.0 .AND. nost.GT.4) THEN

        !! reference ambiguity if all cycle slips
        IF (nsip .EQ. nost) THEN
           bsip(ksat)=-1
           wcls(ksat)=NINT(comc(1,ksat)/SAT(ksat).lamda(1)-comc(2,ksat)/SAT(ksat).lamda(2)-&
                           pomc(1,ksat)/SAT(ksat).lamda(1)+pomc(2,ksat)/SAT(ksat).lamda(2))
           ncls(ksat)=NINT(comc(1,ksat)/SAT(ksat).lamda(1)-pomc(1,ksat)/SAT(ksat).lamda(1))
        END IF

        !! wide-lane cycle-slip resolution
        CALL rep_widelane(SAT,CKF.nprn,npar,ksat,bsip,OB.amat,pomc,comc,wcls)

        nfix=COUNT(wcls(1:CKF.nprn) .NE. 0.5d0)
        ksat=pointer_int(CKF.nprn,bsip,-1)
        IF (nfix.EQ.1 .AND. ksat.NE.0) THEN  ! failed repair
          nfix =0
          wcls(ksat)=0.5d0
          ncls(ksat)=0.5d0
        END IF

        !! narrow-lane cycle-slip resolution
        ! if widelane fixed successfully
        IF (nfix .GT. 0) THEN
          ! number of available satellites
          nost=0
          ksat=0
          elev=0.d0
          DO isat=1, CKF.nprn
            IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

            IF ((bsip(isat).LE.0.AND.comc(1,isat).NE.0.d0) .OR. & ! good satellite
                (bsip(isat).GT.0.AND.wcls(isat).NE.0.5d0)) THEN ! widelane resolved
              nost=nost+1
              IF (OB.elev(isat) .GT. elev) THEN                 ! a new reference satellite
                ksat=isat
                elev=OB.elev(isat)
              END IF
            END IF
            ! failed cycle slip repair
            IF (bsip(isat).GT.0 .AND. wcls(isat).EQ.0.5d0) THEN
              bsip(isat)=0
              pomc(1:2,isat)=0.d0
              comc(1:2,isat)=0.d0
            END IF
          END DO

          !! begin to fix L1 cycle slips
          nfix=0
          IF (nost .GT. 4) THEN
            CALL rep_l1phase(SAT,CKF.nprn,npar,ksat,bsip,OB.amat,pomc,comc,wcls,ncls)
            nfix=COUNT(ncls(1:CKF.nprn) .NE. 0.5d0)
            ksat=pointer_int(CKF.nprn,bsip,-1)
            IF (nfix.EQ.1 .AND. ksat.NE.0) THEN  ! failed repair
              nfix=0
              wcls(ksat)=0.5d0
              ncls(ksat)=0.5d0
            END IF
          END IF
        END IF

        !! update AM and OB
        IF (nfix .GT. 0) THEN
          WRITE(OUTPUT_UNIT,'(A,I7,F9.2,I3)') 'SLIP REP', CKF.mjd, CKF.sod, nfix
          DO isat=1, CKF.nprn
            IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE
            IF (wcls(isat).NE.0.5d0 .AND. ncls(isat).NE.0.5d0) THEN
              ind=OB.ltog(iamb,isat)-NM.npc
              OB.flag(isat,1)=0
              AM(ind).wcls=wcls(isat)
              AM(ind).ncls=ncls(isat)
              AM(ind).abin=AM(ind).abin+wcls(isat)
              AM(ind).xini=AM(ind).xini+SAT(isat).lamdw*wcls(isat)/(1.d0+SAT(isat).g)+SAT(isat).lamdn*ncls(isat)
              IF (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.EQ.2) THEN
                AM(ind).abwl=AM(ind).abwl+wcls(isat)
              ELSE IF (CKF.lamb.EQ..TRUE. .AND. AM(ind).ifab.EQ.4) THEN
                AM(ind).abwl=AM(ind).abwl+wcls(isat)
                AM(ind).abnl=AM(ind).abnl+ncls(isat)
              END IF
              WRITE(OUTPUT_UNIT,'(A,A3,2F8.1)') '---- SAT ',CKF.cprn(isat), AM(ind).wcls, AM(ind).ncls
            END IF
          END DO
        END IF
      END IF
    END IF
  END DO

  RETURN

END SUBROUTINE
