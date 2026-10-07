!*
SUBROUTINE cyslip(CKF,SAT,OB,OBH,AM,NM)
!!
!! REPAIRING CYCLE SLIPS FOR RAW OBSERVATION BY EPOCH DIFFERENCE
!! JING GUO
!! 2014/01/06
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
!!------------------------
TYPE(SATE) :: SAT(MAXSAT)
TYPE(RNXOBS) :: OB,OBH
TYPE(CKDCFG) :: CKF
TYPE(AMBT) :: AM(1:*)
TYPE(INFM) :: NM

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i,k,ind,isat,ksat,iamb,nost,nsip
  INTEGER(IT) :: bsip(MAXSAT),isys,npar,ifreq,nfix
  REAL(RL) :: rcls(MAXSAT)
  REAL(RL) :: elev,pomc(MAXSAT),comc(MAXSAT)
  CHARACTER(LEN=5) :: camb

  !*
  ! The function called
  !!----------------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: pointer_int
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!----------------------------

  npar=0
  IF (pointer_string(OB.npar,OB.pname,'STAPX') .NE. 0) npar=3

  DO isys=1, CKF.nsys

    k=INDEX(SYS,CKF.system(isys:isys))

    DO ifreq=1, CKF.nfq(k)

      WRITE(camb,'(A4,I1)') 'AMBL',ifreq

      !! whether repair or not
      iamb=pointer_string(OB.npar,OB.pname,camb)
      DO isat=1, CKF.nprn
        IF (OB.omc(isat,ifreq).NE.0.d0 .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(isys:isys)) THEN
          IF (OB.flag(isat,ifreq).NE.0 .AND. OB.ltog(iamb,isat).NE.0) EXIT
        END IF
      END DO
      IF (isat .GT. CKF.nprn) CYCLE

      !! form epoch difference
      rcls=0.5d0
      IF (timdif(CKF.mjd,CKF.sod,OBH.jd,OBH.tsec) .LE. 30.d0) THEN

        ! number of available satellites
        nost=0
        ! number of to-be-repaired cycle slips
        nsip=0
        ! reference satellite
        ksat=0
        ! maximum elevation
        elev=0.d0
        DO isat=1, CKF.nprn
          bsip(isat)=0       ! flag for slip status
          pomc(isat)=0.d0
          comc(isat)=0.d0

          IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE

          IF (OB.omc(isat,ifreq).NE.0.d0 .AND. OBH.omc(isat,ifreq).NE.0.d0) THEN
            nost=nost+1
            comc(isat)=OB.omc(isat,ifreq)-OBH.omc(isat,ifreq)
            pomc(isat)=OB.omc(isat,MAXFREQ+ifreq)-OBH.omc(isat,MAXFREQ+ifreq)
            IF (OB.flag(isat,ifreq) .NE. 0) THEN
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
            !! ionosphere prediction is nessary for this
            !! this should be ingestivated further
            rcls(ksat)=NINT((comc(ksat)-pomc(ksat))/SAT(ksat).lamda(ifreq))
          END IF

          CALL rep_rawphase(SAT,CKF.nprn,npar,ifreq,ksat,bsip,OB.amat,pomc,comc,rcls)

          nfix=COUNT(rcls(1:CKF.nprn) .NE. 0.5d0)
          ksat=pointer_int(CKF.nprn,bsip,-1)
          IF (nfix.EQ.1 .AND. ksat.NE.0) THEN  ! failed repair
            nfix=0
            rcls(ksat)=0.5d0
          END IF

          !! update AM and OB
          IF (nfix .GT. 0) THEN
            WRITE(OUTPUT_UNIT,'(A,I7,F9.2,I3,I3)') 'SLIP REP',CKF.mjd,CKF.sod,nfix,ifreq
            DO isat=1, CKF.nprn
              IF (CKF.cprn(isat)(1:1) .NE. CKF.system(isys:isys)) CYCLE
              IF (rcls(isat) .NE. 0.5d0) THEN
                ind=OB.ltog(iamb,isat)-NM.npc
                OB.flag(isat,ifreq)=0
                AM(ind).wcls=rcls(isat)
                AM(ind).abin=AM(ind).abin+rcls(isat)
                AM(ind).xini=AM(ind).xini+SAT(isat).lamda(ifreq)*rcls(isat)
                WRITE(OUTPUT_UNIT,'(A,A3,I3,F8.1)') '---- SAT ',CKF.cprn(isat),ifreq,AM(ind).wcls
              END IF
            END DO
          END IF
        END IF
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE
