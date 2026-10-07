!*
SUBROUTINE corr_hisi_1(CKF,OB,SAT,SIT,isit)
!!
!! Correction bias for HISI of 2023 ...
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
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT
INTEGER(IT) :: isit

    !*
    ! The local variables
    !!------------------------
    INTEGER(IT) :: i,isat,ifreq,nobs,nslip,isys
    LOGICAL(LG) :: lfirst(MAXSIT)
    REAL(RL) :: lcode(MAXSAT,MAXSIT),lphase(MAXSAT,MAXSIT),delta(MAXSAT),bias,slop

    DATA lfirst /MAXSIT*.TRUE./
    SAVE lfirst,lcode,lphase,delta

    INTEGER(IT) :: pointer_string
    !*
    ! Start the exectuable code
    !!-------------------------------

    ! --- 矫正伪距常偏和相位趋势项误差[P1-P2/L1-L2]
    DO isat=1,CKF%nprn
        IF (OB%obs(isat,1+MAXFREQ) .NE. 0) THEN

            IF (CKF%cprn(isat)(1:1).EQ.'C') THEN
                isys=INDEX(SYS,CKF%cprn(isat)(1:1))
                ! 2023299
                ! --- B1I:699m       B1C: 603m
                ! --- B1I:-3.6483d-4 B1C: 3.2858d-4
                IF (CKF%freq(1,isys).EQ.'L2') THEN
                    bias=699.d0
                    slop=-3.6483d-4
                ELSE IF (CKF%freq(1,isys).EQ.'L1') THEN 
                    bias=603.d0
                    slop=3.2858d-4
                END IF
                OB%obs(isat,MAXFREQ+1) = OB%obs(isat,MAXFREQ+1) + bias
                OB%obs(isat,1)         = OB%obs(isat,1)         + bias/SAT(isat)%lamda(1)                
                OB%obs(isat,1) = OB%obs(isat,1) + slop/SAT(isat)%lamda(1)*((CKF%mjd-CKF%mjd0)*86400.0+CKF%sod-CKF%sod0)        
            ELSE
                OB%obs(isat,MAXFREQ+1) = OB%obs(isat,MAXFREQ+1) + 603.0
                OB%obs(isat,1)         = OB%obs(isat,1)         + 603.0/SAT(isat)%lamda(1)                        
                OB%obs(isat,1) = OB%obs(isat,1) + ( 3.2858d-4)/SAT(isat)%lamda(1)*((CKF%mjd-CKF%mjd0)*86400.0+CKF%sod-CKF%sod0)
            END IF

        END IF
    END DO

    ! --- 修复钟跳的码相不一致[CMC STation Diff]
    ifreq = 1
    nobs = 0
    nslip = 0
    IF (lfirst(isit)) THEN
        lfirst(isit) = .FALSE.
        DO isat=1,CKF%nprn
            lcode(isat,isit) = 0.d0
            lphase(isat,isit) = 0.d0
        END DO
    ELSE
        DO isat=1,CKF%nprn
            IF (OB%obs(isat,ifreq) .NE. 0.d0) THEN
                nobs = nobs + 1
            ELSE
                CYCLE
            END IF
            IF (abs(OB%obs(isat,ifreq+MAXFREQ)-lcode(isat,isit)) .GE. 1.d-3*VEL_LIGHT-1000.d0) THEN
                nslip = nslip + 1
            END IF
        END DO
    END IF

    DO isat=1,CKF%nprn
        lcode(isat,isit)  = OB%obs(isat,ifreq+MAXFREQ)
        lphase(isat,isit) = OB%obs(isat,ifreq)
    END DO

    IF (nobs.EQ.nslip .AND. nobs.NE.0) THEN
        DO isat=1,CKF%nprn
            IF (CKF%cprn(isat)(1:1).EQ.'C') THEN
                delta(isat) = delta(isat) + 2.22d0/SAT(isat)%lamda(ifreq)
            ELSE
                delta(isat) = delta(isat) + 2.10d0/SAT(isat)%lamda(ifreq)
            END IF
        END DO
    END IF

    DO isat=1,CKF%nprn
        IF (OB%obs(isat,ifreq) .EQ. 0.d0) CYCLE
        OB%obs(isat,ifreq) = OB%obs(isat,ifreq) - delta(isat)
    END DO

    RETURN

END SUBROUTINE

!*
SUBROUTINE corr_hisi_2(CKF,OB,SAT,SIT,isit)
!!
!! Correction bias for HISI of 2025 ...
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
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT
INTEGER(IT) :: isit

    !*
    ! The local variables
    !!------------------------
    INTEGER(IT) :: i,isat,ifreq,nobs,nslip,isys
    LOGICAL(LG) :: lfirst(MAXSIT)
    REAL(RL) :: lcode(MAXSAT,MAXSIT),lphase(MAXSAT,MAXSIT),delta(MAXSAT,MAXFREQ),bias,slop

    DATA lfirst /MAXSIT*.TRUE./
    SAVE lfirst,lcode,lphase,delta

    INTEGER(IT) :: pointer_string
    !*
    ! Start the exectuable code
    !!-------------------------------

    ! --- 矫正伪距常偏和相位趋势项误差[P1-P2/L1-L2]
    DO isat=1,CKF%nprn
        IF (OB%obs(isat,1+MAXFREQ) .NE. 0) THEN

            isys=INDEX(SYS,CKF%cprn(isat)(1:1))
            DO ifreq=1,CKF%nfreq(isys)
                bias=0.d0
                slop=0.d0
                ! 2025113: L1-L2 [m]
                ! --- B1I:818m           B1C: 730m
                ! --- B1I:-0.003211568d0 B1C: 0.002893605d0
                IF (CKF%freq(ifreq,isys).EQ.'L2') THEN
                    bias=818.d0
                    slop=-0.003211568d0
                ELSE IF (CKF%freq(ifreq,isys).EQ.'L1') THEN 
                    bias=730.d0
                    slop= 0.002893605d0
                END IF            
                OB%obs(isat,MAXFREQ+ifreq) = OB%obs(isat,MAXFREQ+ifreq) + bias
                OB%obs(isat,ifreq)         = OB%obs(isat,ifreq)         + bias/SAT(isat)%lamda(ifreq)                
                OB%obs(isat,ifreq) = OB%obs(isat,ifreq) + slop/SAT(isat)%lamda(ifreq)*((CKF%mjd-CKF%mjd0)*86400.0+CKF%sod-CKF%sod0)   
            END DO     
        END IF
    END DO

    ! --- 修复钟跳的码相不一致[CMC STation Diff]
    nobs = 0
    nslip = 0
    IF (lfirst(isit)) THEN
        lfirst(isit) = .FALSE.
        DO isat=1,CKF%nprn
            lcode(isat,isit) = 0.d0
            lphase(isat,isit) = 0.d0
        END DO
    ELSE
        DO isat=1,CKF%nprn
            IF (OB%obs(isat,1) .NE. 0.d0 .AND. lphase(isat,isit).NE.0.d0) THEN
                nobs = nobs + 1
            ELSE
                CYCLE
            END IF
            IF (abs(OB%obs(isat,1+MAXFREQ)-lcode(isat,isit)) .GE. 1.d-3*VEL_LIGHT-1000.d0) THEN
                nslip = nslip + 1
            END IF
        END DO
    END IF

    DO isat=1,CKF%nprn
        lcode(isat,isit)  = OB%obs(isat,1+MAXFREQ)
        lphase(isat,isit) = OB%obs(isat,1)
    END DO

    IF (nobs.EQ.nslip .AND. nobs.NE.0) THEN
        WRITE(OUTPUT_UNIT,'(A,I6,2X,F8.1,2X,I6)') '... PLEASE BE MIND THE CLOCK SLIP FOR HISI DEVICE '//SIT%name//' ON ', CKF%mjd,CKF%sod,nobs
        DO isat=1,CKF%nprn
            isys=INDEX(SYS,CKF%cprn(isat)(1:1))
            DO ifreq=1,CKF%nfreq(isys)            
                IF (CKF%freq(ifreq,isys).EQ.'L2') THEN
                    delta(isat,ifreq) = delta(isat,ifreq) - 4.d0/SAT(isat)%lamda(ifreq)
                ELSE IF (CKF%freq(ifreq,isys).EQ.'L1') THEN 
                    delta(isat,ifreq) = delta(isat,ifreq) - 4.d0/SAT(isat)%lamda(ifreq)
                END IF
            END DO
        END DO
    END IF

    DO isat=1,CKF%nprn
        isys=INDEX(SYS,CKF%cprn(isat)(1:1))
        DO ifreq=1,CKF%nfreq(isys)           
            IF (OB%obs(isat,ifreq) .EQ. 0.d0) CYCLE
            OB%obs(isat,ifreq) = OB%obs(isat,ifreq) - delta(isat,ifreq)
        END DO
    END DO

    RETURN

END SUBROUTINE