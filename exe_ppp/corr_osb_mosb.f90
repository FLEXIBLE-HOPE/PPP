!*
SUBROUTINE corr_osb_mosb(CKF,OB,SAT,SIT,isit)
!!
!! Correction for Mult-frequency OSB, PRIDE or WUMFIN
!! Created by Shengyi Xu
!! 2023-11-29
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
  INTEGER(IT) :: isat,ifreq,iy,imon,id,ih,im,freq,ind,iepoch,ip,i
  REAL(RL) :: sec,osbInNs,cbia
  CHARACTER(LEN=4) :: flag

  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-------------------------------

  !! time tag
  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)
  WRITE(1008,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
  
  !code osb or cc2nocc
  IF (INDEX(TRIM(CKF.codebias),'osb').NE.0) THEN
    DO isat=1, CKF.nprn
      DO ifreq=1, MAXFREQ
        flag = ''
        IF (OB.obs(isat,ifreq+MAXFREQ) .NE. 0) THEN
          READ(OB.fob(isat,ifreq+MAXFREQ)(2:2),'(I1)')freq
          ind = INDEX(OBSTYPE,OB.fob(isat,ifreq+MAXFREQ)(3:3))
          ! 900s
          DO ip=1, 96
            IF (CKF%sod .GE. SAT(isat)%osbTimeSpan(1,freq,ind,ip,1) .AND. CKF%sod .LT. SAT(isat)%osbTimeSpan(1,freq,ind,ip,2)) THEN
              iepoch = ip
              EXIT
            END IF
          END DO

          !xsy:伪距osb没有对应通道用其他通道补全
          IF (SAT(isat)%osbBias(1,freq,ind,iepoch) .EQ. -999.d0) THEN
            DO i=1, 16
              IF (SAT(isat)%osbBias(1,freq,i,iepoch) .NE. -999.d0) THEN
                ind = i
                flag = 'NOT'
                EXIT
              END IF
            END DO
          END IF

          IF (SAT(isat)%osbBias(1,freq,ind,iepoch) .EQ. -999.d0) THEN
            flag = 'NONE'
            CKF%nofixsat = TRIM(CKF%nofixsat)//' '//CKF.cprn(isat)
            cbia = 0.d0
          ELSE
            cbia = SAT(isat)%osbBias(1,freq,ind,iepoch)*VEL_LIGHT*1.0d-9
            OB.obs(isat,ifreq+MAXFREQ) = OB.obs(isat,ifreq+MAXFREQ) - cbia
          END IF

          WRITE(1008,'(2X,A5,(A),2A4,1X,F14.4,A5,F14.3,A4,F14.3)')SIT.name,'  CODE: ',CKF.cprn(isat),OB.fob(isat,ifreq+MAXFREQ),SAT(isat)%osbBias(1,freq,ind,iepoch),&
                flag,OB.obs(isat,ifreq+MAXFREQ) + cbia,' => ',OB.obs(isat,ifreq+MAXFREQ)
        END IF
      END DO
    END DO
  END IF

  IF (CKF.liar .EQ. .TRUE.) THEN
    DO isat=1, CKF.nprn
      DO ifreq=1, MAXFREQ
        flag = ''
        IF (OB.obs(isat,ifreq) .NE. 0) THEN
          READ(OB.fob(isat,ifreq)(2:2),'(I1)')freq
          ind = INDEX(OBSTYPE,OB.fob(isat,ifreq)(3:3))
          ! 900s
          DO ip=1, 96
            IF (CKF%sod .GE. SAT(isat)%osbTimeSpan(2,freq,ind,ip,1) .AND. CKF%sod .LT. SAT(isat)%osbTimeSpan(2,freq,ind,ip,2)) THEN
              iepoch = ip
              EXIT
            END IF
          END DO

          !xsy:相位osb没有对应通道用其他通道补全
          IF (SAT(isat)%osbBias(2,freq,ind,iepoch) .EQ. -999.d0) THEN
            DO i=1, 16
              IF (SAT(isat)%osbBias(2,freq,i,iepoch) .NE. -999.d0) THEN
                ind = i
                flag = 'NOT'
                EXIT
              END IF
            END DO
          END IF

          IF (SAT(isat)%osbBias(2,freq,ind,iepoch) .EQ. -999.d0) THEN
            flag = 'NONE'
            CKF%nofixsat = TRIM(CKF%nofixsat)//' '//CKF.cprn(isat)
            cbia = 0.d0
          ELSE
            ! LOSB = b0 + slope(t-t0)
            osbInNs = SAT(isat)%osbBias(2,freq,ind,iepoch) + SAT(isat)%osbBiasSlope(freq,ind,iepoch)*(CKF%sod-SAT(isat)%osbTimeSpan(2,freq,ind,iepoch,1))
            ! IF (CKF%cprn(isat)(1:1) .EQ. 'G') THEN
            !   WRITE(*,'(A3,F10.4,F10.4,F16.10,F8.1)')CKF%cprn(isat),osbInNs,SAT(isat)%osbBias(2,freq,ind,iepoch),SAT(isat)%osbBiasSlope(freq,ind,iepoch),(CKF%sod-SAT(isat)%osbTimeSpan(2,freq,ind,iepoch,1))
            ! END IF  
            cbia = (osbInNs*VEL_LIGHT*1.0d-9)/(VEL_LIGHT/SAT(isat).freq(ifreq))
            OB.obs(isat,ifreq) = OB.obs(isat,ifreq) - cbia
          END IF
          
          WRITE(1008,'(2X,A5,(A),2A4,1X,F14.4,A5,F14.3,A4,F14.3)')SIT.name,' PHASE: ',CKF.cprn(isat),OB.fob(isat,ifreq),osbInNs,&
                flag,OB.obs(isat,ifreq)+cbia,' => ',OB.obs(isat,ifreq)
        END IF
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE