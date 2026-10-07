!
!! purpose   : add the new ambiguities to the estimator. Newly occured or detected
!!             is indicated by OB.flag.ne.0. Parameter table must be updated and 
!!             so do the estimator.
!! parameter :
!!    input  : SCF    -- SRIF struct
!!             OB     -- observation struct
!!             NM,AM  -- normal matrix & ambiguity table
!! author    : Geng J
!! created   : Nov. 11, 2007
!
SUBROUTINE ppp_add_newamb(CKF,OB,AM,NM,infs)
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(INFM) :: NM
TYPE(AMBT) :: AM(1:*)
REAL(RL) :: infs(1:*)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,isat,iamb,nbias
  INTEGER(IT) :: ifreq,isys
  CHARACTER(LEN=5) :: camb

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  nbias=0

  DO isat=1, CKF.nprn

    isys=INDEX(SYS,CKF.cprn(isat)(1:1))

    DO ifreq=1, CKF.nfq(isys)

      WRITE(camb,'(A4,I1)') 'AMBL',ifreq
      iamb=pointer_string(OB.npar,OB.pname,camb)

      IF (OB.omc(isat,ifreq) .EQ. 0.d0) CYCLE
      IF (OB.flag(isat,ifreq) .EQ. 1) THEN

        !! add a new ambiguity to primary filter
        nbias=nbias+1
        IF (NM.imtx+nbias+1.GT.NM.nmtx .OR. NM.ns+nbias.GT.NM.namb) THEN
          WRITE(ERROR_UNIT,'(A,I5,F15.7,2I5)') '***ERROR(ppp_add_newamb): matrix too small ',CKF.mjd,CKF.sod
          CALL exit(1)
        END IF
        AM(NM.ns+nbias).pname=camb
        AM(NM.ns+nbias).ifab =0
        AM(NM.ns+nbias).iobs =0
        AM(NM.ns+nbias).elev =0.d0
        AM(NM.ns+nbias).psat =isat
        AM(NM.ns+nbias).ifreq =ifreq
        IF (CKF.llog .EQ. .FALSE.) THEN
          AM(NM.ns+nbias).ptime=CKF.mjd+CKF.sod/86400.d0   ! depends on real observations
        ELSE
          AM(NM.ns+nbias).ptime(1)=OB.lifamb(isat,ifreq,1)
          AM(NM.ns+nbias).ptime(2)=OB.lifamb(isat,ifreq,2)
        END IF
        AM(NM.ns+nbias).xini =0.d0
        AM(NM.ns+nbias).xcor =0.d0
        AM(NM.ns+nbias).xsig =0.d0
        AM(NM.ns+nbias).abin =0.d0
        AM(NM.ns+nbias).abein =0.d0
        AM(NM.ns+nbias).weig =0.d0
        AM(NM.ns+nbias).eweig =0.d0
        AM(NM.ns+nbias).eeweig=0.d0
        AM(NM.ns+nbias).heweig=0.d0 
        AM(NM.ns+nbias).xrwl =0.d0
        AM(NM.ns+nbias).xswl =0.d0
        AM(NM.ns+nbias).xrewl =0.d0
        AM(NM.ns+nbias).xsewl =0.d0
        AM(NM.ns+nbias).xreewl =0.d0
        AM(NM.ns+nbias).xseewl =0.d0
        AM(NM.ns+nbias).xrhewl =0.d0
        AM(NM.ns+nbias).xshewl =0.d0

        AM(NM.ns+nbias).hewcls=0.5d0
        AM(NM.ns+nbias).eewcls=0.5d0
        AM(NM.ns+nbias).ewcls =0.5d0
        AM(NM.ns+nbias).wcls =0.5d0
        AM(NM.ns+nbias).ncls =0.5d0
        OB.ltog(iamb,isat)=NM.imtx+nbias !@ CMT BY XSY: imtx-npc = ns -> AM

        !! add a new ambiguity to auxiliary filter
        IF (CKF.liar) THEN
          AM(NM.ns+nbias).abhewl=0.5d0
          AM(NM.ns+nbias).abeewl=0.5d0
          AM(NM.ns+nbias).abewl=0.5d0
          AM(NM.ns+nbias).abwl=0.5d0
          AM(NM.ns+nbias).abnl=0.5d0
        END IF
      END IF
    END DO
  END DO
  IF (nbias .EQ. 0) RETURN

  !! shift the rightside of S to the last column and clean for primary filter [右移右函数]
  DO i=1,NM.imtx+nbias
    IF (i .GT. NM.imtx) THEN
      infs(NM.iptx(NM.imtx+1+nbias)+i)=0.d0
    ELSE
      infs(NM.iptx(NM.imtx+1+nbias)+i)=infs(NM.iptx(NM.imtx+1)+i)
    END IF
  END DO
  DO j=1,NM.imtx+nbias
    DO i=1,NM.imtx+nbias
      IF (i.LE.NM.imtx .AND. j.LE.NM.imtx) CYCLE
      infs(NM.iptx(j)+i)=0.d0
      IF (i .EQ. j) infs(NM.iptx(j)+i)=1.d-4
      !! For LUTAN1A LEO RAW PPP
      ! IF (CKF%cobs(1:3).EQ.'RAW') THEN
      !   ! IF (i .EQ. j) infs(NM.iptx(j)+i)=1.d-10 !@ CMT BY XSY: FOR LEO ONBOARD RAW PPP, IT SHOULD BE LAGER
      ! ELSE
      !   IF (i .EQ. j) infs(NM.iptx(j)+i)=1.d-4
      ! END IF      
    END DO
  END DO
  NM.ns=NM.ns+nbias
  NM.imtx=NM.imtx+nbias ! PM:np+nc AM:ns

  RETURN

END SUBROUTINE
