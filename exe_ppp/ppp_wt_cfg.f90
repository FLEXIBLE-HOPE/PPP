!*
!! CREATED BY SHENGYI XU
!! Output the CFG to File Header
!*
SUBROUTINE ppp_wt_cfg(CKF,SIT,lfn)
!*
USE ckdctrl
USE station
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(SITE) :: SIT(1:*)
INTEGER(IT) :: lfn

    !*
    ! The local variables
    !!-----------------------
    INTEGER(IT) :: iy0,imon0,id0,ih0,im0,iy1,imon1,id1,ih1,im1,isys,ifreq,isit,ind,mjd0,mjd1
    REAL(RL) :: isec0,isec1,sod0,sod1
    CHARACTER(LEN_STRING) :: pppar,line
    INTEGER(IT) :: pointer_string
    !*
    ! Start the exectuable code
    !!-------------------------------
    
    CALL mjd2date(CKF%mjd0,CKF%sod0,iy0,imon0,id0,ih0,im0,isec0)
    CALL mjd2date(CKF%mjd1,CKF%sod1,iy1,imon1,id1,ih1,im1,isec1)

    WRITE(lfn,'(A)')'##-> OUTPUT by PANDA PPP software'
    WRITE(lfn,'(A30,I5,4I3,F5.1,I8,F10.1)')' %-> Start Time              :',iy0,imon0,id0,ih0,im0,isec0,CKF%mjd0,CKF%sod0
    WRITE(lfn,'(A30,I5,4I3,F5.1,I8,F10.1)')' %-> End   Time              :',iy1,imon1,id1,ih1,im1,isec1,CKF%mjd1,CKF%sod1
    WRITE(lfn,'(A30,A100)')                '#%-> Frequency Used          :',TRIM(CKF%freq_used)
    WRITE(lfn,'(A30,F20.1)')               ' %-> Sampling Rate [s]       :',CKF%dintv
    IF (CKF%lli_flag) THEN
        WRITE(lfn,'(A30,I8,F6.1,F6.1,A10)')' %-> QC Gap/LG/LW            :',CKF%gap,CKF%lg,CKF%lw,'LLI'
    ELSE
        WRITE(lfn,'(A30,I8,F6.1,F6.1)')    ' %-> QC Gap/LG/LW            :',CKF%gap,CKF%lg,CKF%lw
    END IF
    WRITE(lfn,'(A30,A10,A10)')             ' %-> Observation Combination :',TRIM(CKF%dobs),TRIM(CKF%cobs)
    WRITE(lfn,'(A30,A20)')                 ' %-> ION Model               :',TRIM(CKF%ionmod)
    WRITE(lfn,'(A30,A20)')                 ' %-> ZTD Model [min]         :',TRIM(CKF%ztdmod)
    WRITE(lfn,'(A30,A20)')                 ' %-> ZTD Gradient [min]      :',TRIM(CKF%grdmod)
    WRITE(lfn,'(A30,F10.3,F10.3)')         ' %-> PPPRTK for ION/ZTD [m]  :',CKF%ionConstrians,CKF%ztdConstrians
    IF (CKF%lsitcons.EQ..TRUE.) THEN
        WRITE(lfn,'(A30,3F10.3)')          ' %-> STA XYZ constrains [m]  :',CKF%sitConstrians(1:3)
    END IF
    IF (CKF%lisb) THEN
        WRITE(lfn,'(A30,A20)')             ' %-> Other System Rclock     :',' ISB CONSTANT'
    ELSE
        WRITE(lfn,'(A30,A20)')             ' %-> Other System Rclock     :',' ISB WHITE NOISE'
    END IF
    IF (CKF%liar.EQ..TRUE.)THEN
        pppar='TRUE'
    ELSE
        pppar='FALSE'
    END IF
    IF (CKF%lockout(1).NE.0.D0 .AND. CKF%lockout(2).NE.0.D0) THEN
        mjd0=INT(CKF%lockout(1))
        sod0=(CKF%lockout(1)-mjd0)*86400.d0
        mjd1=INT(CKF%lockout(2))
        sod1=(CKF%lockout(2)-mjd1)*86400.d0
        WRITE(lfn,'(A30,I10,F10.1,A4,I10,F10.1)') ' %-> LOCKOUT time span       :',mjd0,sod0,'-->',mjd1,sod1        
    END IF
    WRITE(lfn,'(A30,A20)')                 ' %-> PPP-AR                  :',TRIM(pppar)
    WRITE(lfn,'(A30,A20)')                 ' %-> PPP-AR FixMode          :',TRIM(CKF%ArMode)
    WRITE(lfn,'(A30,A20)')                 ' %-> PPP-AR Phase Bias Mode  :',TRIM(CKF%ArBiasMode)
    WRITE(lfn,'(A30,A20)')                 ' %-> PPP-AR OSB Mode         :',TRIM(CKF%osbMode)
    WRITE(lfn,'(A30,I10,I10)')             ' %-> PPP-PAR MAX Del and Sav :',CKF%nl_maxdel,INT(CKF%nl_minsav)
    WRITE(lfn,'(A30,F10.1,F10.1)')         ' %-> PPP-AR Ratio for NL/WL  :',CKF%nl_ratio,CKF%wl_ratio
    WRITE(lfn,'(A30,F10.1,F10.1)')         ' %-> PPP-AR CutElev for NL/WL:',CKF%cutoff,CKF%wlcutoff
    WRITE(lfn,'(A30,A20)')                 ' %-> Code Bias Mode          :',TRIM(CKF%codebias)
    WRITE(lfn,'(A30,I20)')                 ' %-> Mask SNR                :',CKF%SnrLimit
    WRITE(lfn,'(A30,F20.1)')               ' %-> DIA Thershold           :',CKF%DiaSig
    WRITE(lfn,'(A30,I20)')                 ' %-> Reconvergence [s]       :',CKF%ReConvTime
    DO isit=1,CKF%nsit
        WRITE(lfn,'(A30,A5,A5,7F7.1,2X,A20,7X,A20)')        ' %-> SIT X/Y/Z sig0/qx0 [m]  :',SIT(isit)%name,SIT(isit)%skd,SIT(isit)%cutoff*RAD2DEG,SIT(isit)%dx0(1:3),SIT(isit)%qx(1:3),&
                                            SIT(isit)%rectyp,SIT(isit)%anttyp
        line = ''
        DO isys=1, CKF%nsys
            line = TRIM(line)//' '//CKF%system(isys:isys)//':'
            ind = INDEX(SYS,CKF%system(isys:isys))
            
            DO ifreq=1,CKF%nfreq(ind)
                IF (TRIM(SIT(isit)%freqChannel(ind,ifreq)) .EQ. '') CYCLE
                line = TRIM(line)//SIT(isit)%freqChannel(ind,ifreq)
                IF (ifreq .NE. CKF%nfreq(ind)) THEN
                    line = TRIM(line)//'_'
                END IF
            END DO
        END DO
        WRITE(lfn,'(A30,A5,A5,A50)')        ' %-> SIT Signal Channel      :',SIT(isit)%name,SIT(isit)%skd,TRIM(line)
        WRITE(lfn,'(A31,A150)')             ' %-> Observation File        : ',ADJUSTL(SIT(isit)%obsfile)
    END DO
    WRITE(lfn,'(A)')'##-> END OF HEADER'
    RETURN
END SUBROUTINE
