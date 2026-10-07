!*
SUBROUTINE corr_osb(CKF,OB,SAT,SIT,isit)
!!
!! Correction for one-day OSB [ICLK]
!!
!*
USE const
USE ckdctrl
USE observation
USE satellite
USE station
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT
TYPE(SATE) :: SAT(MAXSAT)
INTEGER(IT) :: isit

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: isat,ifreq,iy,imon,id,ih,im 
  LOGICAL(LG) :: lfirst(MAXSIT)
  REAL(RL) :: osb_corr(MAXSAT,MAXFREQ,4,MAXSIT),sec
  CHARACTER(LEN=4) :: flag

  DATA lfirst /MAXSIT*.TRUE./
  SAVE osb_corr,lfirst

  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-------------------------------
  !! 读OSB文件
  IF (lfirst(isit) .EQ. .TRUE.) THEN
     CALL read_osb(CKF,SAT,osb_corr,OB,isit)
     lfirst(isit) = .FALSE.
  END IF

  !! time tag
  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)
  WRITE(1008,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
  
  !code osb or cc2nocc
  IF (INDEX(TRIM(CKF.codebias),'osb').NE.0) THEN
    DO isat=1, CKF.nprn
      DO ifreq=1, MAXFREQ
        flag = ''
        IF (OB.obs(isat,ifreq+MAXFREQ) .NE. 0) THEN
          OB.obs(isat,ifreq+MAXFREQ) = OB.obs(isat,ifreq+MAXFREQ) - (osb_corr(isat,ifreq,1,isit)*VEL_LIGHT*1.0d-9)
          IF (osb_corr(isat,ifreq,1,isit) .EQ. 0.d0) THEN
            flag = 'NONE'
            CKF%nofixsat = TRIM(CKF%nofixsat)//' '//CKF.cprn(isat)
          END IF
          WRITE(1008,'(2X,A5,(A),2A4,1X,F14.4,A5,F14.3,A4,F14.3)')SIT.name,'  CODE: ',CKF.cprn(isat),OB.fob(isat,ifreq+MAXFREQ),osb_corr(isat,ifreq,1,isit),&
                flag,OB.obs(isat,ifreq+MAXFREQ)+osb_corr(isat,ifreq,1,isit)*VEL_LIGHT*1.0d-9,' => ',OB.obs(isat,ifreq+MAXFREQ)
        END IF
      END DO
    END DO
  END IF

  IF (CKF.liar .EQ. .TRUE.) THEN
    DO isat=1, CKF.nprn
      DO ifreq=1, MAXFREQ
        flag = ''
        IF (OB.obs(isat,ifreq) .NE. 0) THEN
          OB.obs(isat,ifreq) = OB.obs(isat,ifreq) - (osb_corr(isat,ifreq,3,isit)*VEL_LIGHT*1.0d-9)/(VEL_LIGHT/SAT(isat).freq(ifreq))
          IF (osb_corr(isat,ifreq,3,isit) .EQ. 0.d0) THEN
            flag = 'NONE'
            CKF%nofixsat = TRIM(CKF%nofixsat)//' '//CKF.cprn(isat)
          END IF
          WRITE(1008,'(2X,A5,(A),2A4,1X,F14.4,A5,F14.3,A4,F14.3)')SIT.name,' PHASE: ',CKF.cprn(isat),OB.fob(isat,ifreq),osb_corr(isat,ifreq,3,isit),&
                flag,OB.obs(isat,ifreq)+osb_corr(isat,ifreq,3,isit)*VEL_LIGHT*1.0d-9,' => ',OB.obs(isat,ifreq)
        END IF
      END DO
    END DO 
  END IF

  RETURN

END SUBROUTINE