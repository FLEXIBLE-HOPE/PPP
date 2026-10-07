!*
SUBROUTINE corr_sdb_gxp(CKF,OB,SAT,SIT,HD)
!!
!! Correction SDB from https://www.researchgate.net/profile/Xiaopeng-Gong-2/research
!! Created by Shengyi Xu
!! 2024-01-17
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
TYPE(RNXHEAD) :: HD

    !*
    ! The local variables
    !!------------------------
    INTEGER(IT) :: isat,ifreq,irtype,freq,ichanel
    INTEGER(IT) :: iy,imon,id,ih,im
    REAL(RL) :: sec

    LOGICAL(LG) :: lfirst
    !!xsy: SDB correction for 50 rectype
    !!    10: frequency numbers
    !!    16: OBSTYPE length 
    CHARACTER(LEN=50) :: SDBrtype(MAXRECTYPE)
    REAL(RL) :: SDB(MAXSAT,10,16,MAXRECTYPE)

    DATA lfirst /.TRUE./
    SAVE lfirst,SDBrtype,SDB

    INTEGER(IT) :: pointer_string

    !*
    ! Start the exectuable code
    !!-------------------------------
    !! 读OSB文件
    IF (lfirst .EQ. .TRUE.) THEN
        CALL read_sdb_gxp(CKF,SDBrtype,SDB)
        lfirst = .FALSE.
    END IF

    irtype = pointer_string(MAXRECTYPE,SDBrtype,TRIM(HD%rectype))
    IF (irtype.EQ.0) RETURN

    DO isat=1,CKF%nprn
        DO ifreq=1,MAXFREQ
            IF (OB%obs(isat,ifreq+MAXFREQ).NE.0.d0) THEN
                READ(OB%fob(isat,ifreq+MAXFREQ)(2:2),'(I1)')freq
                ichanel = INDEX(OBSTYPE,OB%fob(isat,ifreq+MAXFREQ)(3:3))
                !WRITE(*,'(A4,2X,A3,2X,A3,2X,F4.2)')SIT%name,CKF%cprn(isat),OB%fob(isat,ifreq+MAXFREQ),SDB(isat,freq,ichanel,irtype)
                OB%obs(isat,ifreq+MAXFREQ) = OB%obs(isat,ifreq+MAXFREQ) + SDB(isat,freq,ichanel,irtype)*VEL_LIGHT*1.0d-9
            END IF
        END DO
    END DO

    RETURN

END SUBROUTINE
