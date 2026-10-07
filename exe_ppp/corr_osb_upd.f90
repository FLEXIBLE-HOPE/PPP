SUBROUTINE corr_osb_upd(CKF,OB,SAT,SIT,isit)
!!
!! XU SHENGYI: CREATED [2023-05-26]  
!! PURPOSE: READ OSB FILE IN THE MODE OF UPD [EPOCH, UPD => OSB]
!*
USE const
USE ckdctrl
USE observation
USE satellite
USE station
USE iso_fortran_env
IMPLICIT NONE
!*
! the argument input
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE)   :: SIT
TYPE(SATE)   :: SAT(MAXSAT)
INTEGER(IT)  :: isit

    !*
    ! the local variables
    !!------------------------
    INTEGER(IT) :: iy,im,id,ih,imi,ind,ifreq,isat,freq,i
    REAL(RL) :: sec
    CHARACTER(LEN=4) :: flag !YES/NOT/NONE

    !*
    ! Start the exectuable code
    !!-------------------------------
    !! time tag
    CALL mjd2date(CKF.mjd,CKF.sod,iy,im,id,ih,imi,sec)
    WRITE(1008,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,im,id,ih,imi,sec,CKF.mjd,CKF.sod
    
    !code osb or cc2nocc
    IF (INDEX(TRIM(CKF.codebias),'osb').NE.0) THEN
        DO isat=1, CKF.nprn
            DO ifreq=1, MAXFREQ
                flag = 'Y'
                IF (OB.obs(isat,ifreq+MAXFREQ) .NE. 0) THEN
                    READ(OB.fob(isat,ifreq+MAXFREQ)(2:2),'(I1)')freq
                    ind = INDEX(OBSTYPE,OB.fob(isat,ifreq+MAXFREQ)(3:3))
                    !xsy:伪距osb没有对应通道用其他通道补全
                    IF (SAT(isat).osbValue(1,freq,ind) .EQ. 0.d0) THEN
                        DO i=1, 16
                            IF (SAT(isat).osbValue(1,freq,i) .NE. 0.d0) THEN
                                ind = i
                                flag = 'NOT'
                                EXIT
                            END IF
                        END DO
                    END IF
                    IF (SAT(isat).osbValue(1,freq,ind) .EQ. 0.d0) THEN
                        flag = 'NONE'
                        CKF.nofixsat = TRIM(CKF.nofixsat)//' '//CKF.cprn(isat)
                    END IF

                    OB.obs(isat,ifreq+MAXFREQ) = OB.obs(isat,ifreq+MAXFREQ) - (SAT(isat).osbValue(1,freq,ind) *VEL_LIGHT*1.0d-9)
                    WRITE(1008,'(2X,A5,(A),A4,(A),A1,2X,A4,4X,A4,1X,2F14.4,F18.3)')SIT.name,'  CODE: ',OB.fob(isat,ifreq+MAXFREQ),' => ',OBSTYPE(ind:ind),flag,&
                                                    CKF.cprn(isat),SAT(isat).osbValue(1,freq,ind),SAT(isat).osbValue(1,freq,ind)*VEL_LIGHT*1.0d-9,OB.obs(isat,ifreq+MAXFREQ)
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
                    IF (SAT(isat).osbValue(2,freq,ind) .EQ. 0.d0) THEN
                        DO i=1, 16
                            IF (SAT(isat).osbValue(2,freq,i) .NE. 0.d0) THEN
                                ind = i
                                flag = 'NOT'
                                EXIT
                            END IF
                        END DO
                    ELSE
                        flag = 'Y'
                    END IF
                    IF (SAT(isat).osbValue(2,freq,ind) .EQ. 0.d0) THEN
                        flag = 'NONE'
                        CKF.nofixsat = TRIM(CKF.nofixsat)//' '//CKF.cprn(isat)
                    END IF

                    OB.obs(isat,ifreq) = OB.obs(isat,ifreq) - (SAT(isat).osbValue(2,freq,ind)*VEL_LIGHT*1.0d-9)/(VEL_LIGHT/SAT(isat).freq(ifreq))
                    WRITE(1008,'(2X,A5,(A),A4,(A),A1,2X,A4,4X,A4,1X,2F14.4,F18.3)')SIT.name,' PHASE: ',OB.fob(isat,ifreq),' => ',OBSTYPE(ind:ind),flag,&
                                                    CKF.cprn(isat),SAT(isat).osbValue(2,freq,ind),(SAT(isat).osbValue(2,freq,ind)*VEL_LIGHT*1.0d-9)/(VEL_LIGHT/SAT(isat).freq(ifreq)),OB.obs(isat,ifreq)
                END IF
            END DO

        END DO
    END IF
    
    RETURN

END SUBROUTINE
