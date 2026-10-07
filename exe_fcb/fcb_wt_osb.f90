!*
!! Created by shengyi xu
!! Date: 2023-07-09
!! Purpose: Map FCB to OSB in UPD mode and write OSB file in sinex
!*
SUBROUTINE fcb_wt_osb(CKF,UPD)
!*
!!
!*
USE ckdctrl
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(CKDCFG) :: CKF
TYPE(FCB) :: UPD

    !*
    ! The local variables
    !!----------------------
    LOGICAL(LG) :: lout,lfirst
    INTEGER(IT) :: lfn,i,k,ilast,iy,imon,id,ih,im,isys,isat,ifreq,year0,doy0,year1,doy1,mjd1,ind_freq,ind_channel
    CHARACTER(LEN=5) :: str_sod_0,str_sod_1
    CHARACTER(LEN=3) :: str_doy0,str_doy1,codetype
    CHARACTER(LEN=1) :: mark
    REAL(RL) :: sod1
    !单天OSB
    !   2 :1-C;2-L
    !   10:freq(1,2,5,6,7,8)
    !   16:channel(PWCIXSAQLDBYMZN)
    REAL(RL) :: OsbValue(MAXSAT,2,10,16)

    DATA lfn /0/ , lfirst /.TRUE./
    SAVE lfn,lfirst,OsbValue

    !*
    ! The function called
    !!----------------------
    INTEGER(IT) :: get_valid_unit

    IF (lfn .EQ. 0) THEN
        lfn=get_valid_unit(10)
        OPEN(lfn,file=CKF.flnosbupd)
    END IF

    IF (lfirst .EQ. .TRUE.) THEN
        OsbValue=0.d0
        lfirst=.FALSE.
        CALL fcb_wt_osbhead(lfn,CKF)
        CALL read_singledayosb(CKF,OsbValue)
    END IF

    CALL mjd2doy(CKF.mjd,year0,doy0)
    CALL timinc(CKF.mjd,CKF.sod,CKF.dintv,mjd1,sod1)
    CALL mjd2doy(mjd1,year1,doy1)

    IF (doy0 .LT. 10) THEN
        WRITE(str_doy0,'(A2,I1)')'00',doy0
    ELSE IF (doy0 .GE. 10 .AND. doy0 .LT. 100) THEN
        WRITE(str_doy0,'(A1,I2)')'0',doy0
    ELSE IF (doy0 .GE. 100) THEN
        WRITE(str_doy0,'(I3)')doy0
    ELSE
        WRITE(output_unit,'(A)') '***ERROR(fcb_wt_osbhead.f90): please check the start time'
        CALL exit(1)
    END IF
    
    IF (doy1 .LT. 10) THEN
        WRITE(str_doy1,'(A2,I1)')'00',doy1
    ELSE IF (doy1 .GE. 10 .AND. doy1 .LT. 100) THEN
        WRITE(str_doy1,'(A1,I2)')'0',doy1
    ELSE IF (doy1 .GE. 100) THEN
        WRITE(str_doy1,'(I3)')doy1
    ELSE
        WRITE(output_unit,'(A)') '***ERROR(fcb_wt_osbhead.f90): please check the end time'
        CALL exit(1)
    END IF

    str_sod_0 = '00000'
    IF (CKF.sod .LT. 10.d0) THEN
        WRITE(str_sod_0,'(A4,I1)')'0000',INT(CKF.sod)
    ELSE IF (CKF.sod .GE. 10.d0 .AND. CKF.sod .LT. 100.d0) THEN
        WRITE(str_sod_0,'(A3,I2)')'000',INT(CKF.sod)
    ELSE IF (CKF.sod .GE. 100.d0 .AND. CKF.sod .LT. 1000.d0) THEN
        WRITE(str_sod_0,'(A2,I3)')'00',INT(CKF.sod)
    ELSE IF (CKF.sod .GE. 1000.d0 .AND. CKF.sod .LT. 10000.d0) THEN
        WRITE(str_sod_0,'(A1,I4)')'0',INT(CKF.sod)
    ELSE IF (CKF.sod .GE. 10000.d0 .AND. CKF.sod .LT. 100000.d0) THEN
        WRITE(str_sod_0,'(I5)')INT(CKF.sod)
    END IF

    str_sod_1 = '00000'
    IF (sod1 .LT. 10.d0) THEN
        WRITE(str_sod_1,'(A4,I1)')'0000',INT(sod1)
    ELSE IF (sod1 .GE. 10.d0 .AND. sod1 .LT. 100.d0) THEN
        WRITE(str_sod_1,'(A3,I2)')'000',INT(sod1)
    ELSE IF (sod1 .GE. 100.d0 .AND. sod1 .LT. 1000.d0) THEN
        WRITE(str_sod_1,'(A2,I3)')'00',INT(sod1)
    ELSE IF (sod1 .GE. 1000.d0 .AND. sod1 .LT. 10000.d0) THEN
        WRITE(str_sod_1,'(A1,I4)')'0',INT(sod1)
    ELSE IF (sod1 .GE. 10000.d0 .AND. sod1 .LT. 100000.d0) THEN
        WRITE(str_sod_1,'(I5)')INT(sod1)
    END IF

    DO isat=1,CKF.nprn
        isys=INDEX(SYS,CKF.cprn(isat)(1:1))

        !写入伪距OSB
        DO ind_freq=1,10
            DO ind_channel=1,16
                IF (OsbValue(isat,1,ind_freq,ind_channel) .EQ. 0.d0) CYCLE
                WRITE(codetype,'(A1,I1,A1)')'C',ind_freq,OBSTYPE(ind_channel:ind_channel)
                WRITE(lfn,'(A11,A3,11X,A3,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                    CKF.cprn(isat),codetype,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',OsbValue(isat,1,ind_freq,ind_channel)
            END DO
        END DO

        DO ifreq=1,CKF.nfreq(isys)
            IF (CKF.cprn(isat)(1:1) .EQ. 'G') THEN
                IF (CKF.freq(ifreq,isys)(2:2) .EQ. '1') THEN
                    mark='C'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '2') THEN
                    mark='W'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '5') THEN
                    mark='I'
                END IF
            ELSE IF (CKF.cprn(isat)(1:1) .EQ. 'J') THEN
                IF (CKF.freq(ifreq,isys)(2:2) .EQ. '1') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '2') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '5') THEN
                    mark='X'
                END IF
            ELSE IF (CKF.cprn(isat)(1:1) .EQ. 'E') THEN
                IF (CKF.freq(ifreq,isys)(2:2) .EQ. '1') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '5') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '6') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '7') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '8') THEN
                    mark='X'
                END IF
            ELSE IF (CKF.cprn(isat)(1:1) .EQ. 'C') THEN
                IF (CKF.freq(ifreq,isys)(2:2) .EQ. '1') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '2') THEN
                    mark='I'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '5') THEN
                    mark='X'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '6') THEN
                    mark='I'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '7') THEN
                    mark='Z'
                ELSE IF (CKF.freq(ifreq,isys)(2:2) .EQ. '8') THEN
                    mark='X'
                END IF
            END IF
            IF (ifreq .EQ. 1) THEN
                IF (UPD.nfcb(isat) .EQ. 10.d0) CYCLE
                WRITE(lfn,'(A11,A3,11X,A2,A1,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                CKF.cprn(isat),CKF.freq(1,isys),mark,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',UPD.nfcb(isat)
            ELSE IF (ifreq .EQ. 2) THEN
                IF (UPD.wfcb(isat) .EQ. 10.d0) CYCLE
                WRITE(lfn,'(A11,A3,11X,A2,A1,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                CKF.cprn(isat),CKF.freq(2,isys),mark,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',UPD.wfcb(isat)
            ELSE IF (ifreq .EQ. 3) THEN
                IF (UPD.ewfcb(isat) .EQ. 10.d0) CYCLE
                WRITE(lfn,'(A11,A3,11X,A2,A1,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                CKF.cprn(isat),CKF.freq(3,isys),mark,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',UPD.ewfcb(isat)
            ELSE IF (ifreq .EQ. 4) THEN
                IF (UPD.eewfcb(isat) .EQ. 10.d0) CYCLE
                WRITE(lfn,'(A11,A3,11X,A2,A1,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                CKF.cprn(isat),CKF.freq(4,isys),mark,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',UPD.eewfcb(isat)
            ELSE IF (ifreq .EQ. 5) THEN
                IF (UPD.hewfcb(isat) .EQ. 10.d0) CYCLE
                WRITE(lfn,'(A11,A3,11X,A2,A1,7X,I4,A1,A3,A1,A5,1X,I4,A1,A3,A1,A5,A3,F24.4)')' OSB       ',&
                CKF.cprn(isat),CKF.freq(5,isys),mark,year0,':',str_doy0,':',str_sod_0,year1,':',str_doy1,':',str_sod_1,' ns',UPD.hewfcb(isat)
            END IF
        END DO
    END DO

    IF (mjd1+sod1/86400 .GT. CKF.mjd1+CKF.sod1/86400) THEN
        WRITE(lfn,'(A137)')'-BIAS/SOLUTION                                                                                                                           '
        WRITE(lfn,'(A137)')'*-------------------------------------------------------------------------------------------------------                                 '
        WRITE(lfn,'(A137)')'%=ENDBIA                                                                                                                                 '
    END IF
    RETURN

END SUBROUTINE

!*
SUBROUTINE fcb_wt_osbhead(lfn,CKF)
!*
!! WRITE OSB file Head
!*
USE ckdctrl
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(CKDCFG) :: CKF
INTEGER(IT) :: lfn

    !*
    ! The local variables
    !!----------------------
    INTEGER(IT) :: i,j,year0,year1,doy0,doy1
    REAL(RL) :: sec
    CHARACTER(LEN=25) :: run_tim
    CHARACTER(LEN=3) :: str_doy0,str_doy1

    !*
    ! Start the exectuable code
    !!-------------------------------
    !*
    CALL mjd2doy(CKF.mjd0,year0,doy0)
    CALL mjd2doy(CKF.mjd1,year1,doy1)

    IF (doy0 .LT. 10) THEN
        WRITE(str_doy0,'(A2,I1)')'00',doy0
    ELSE IF (doy0 .GE. 10 .AND. doy0 .LT. 100) THEN
        WRITE(str_doy0,'(A1,I2)')'0',doy0
    ELSE IF (doy0 .GE. 100) THEN
        WRITE(str_doy0,'(I3)')doy0
    ELSE
        WRITE(output_unit,'(A)') '***ERROR(fcb_wt_osbhead.f90): please check the start time'
        CALL exit(1)
    END IF
    
    IF (year0 .EQ. year1 .AND. doy0 .EQ. doy1) doy1=doy1+1
    IF (doy1 .LT. 10) THEN
        WRITE(str_doy1,'(A2,I1)')'00',doy1
    ELSE IF (doy1 .GE. 10 .AND. doy1 .LT. 100) THEN
        WRITE(str_doy1,'(A1,I2)')'0',doy1
    ELSE IF (doy1 .GE. 100) THEN
        WRITE(str_doy1,'(I3)')doy1
    ELSE
        WRITE(output_unit,'(A)') '***ERROR(fcb_wt_osbhead.f90): please check the end time'
        CALL exit(1)
    END IF
    WRITE(lfn,'(A15,A14,A5,I4,A1,A3,A7,I4,A1,A3,(A))')'%=BIA 1.00 WHU ',run_tim(),' WHU ',year0,':',str_doy0,':00000 ',year1,':',str_doy1,':00000 P 00000 0      '
    WRITE(lfn,'(A137)')'*-------------------------------------------------------------------------------------------------------                                 '
    WRITE(lfn,'(A137)')'* Solution INdependent EXchange Format (SINEX)                                                                                           '
    WRITE(lfn,'(A137)')'*-------------------------------------------------------------------------------------------------------                                 '
    WRITE(lfn,'(A137)')'+FILE/REFERENCE                                                                                                                          '
    WRITE(lfn,'(A137)')' DESCRIPTION        Multi-GNSS Code and Phase Biases for satellites                                                                      '
    WRITE(lfn,'(A137)')' INPUT              WHU Code OSB solutions in SINEX format                                                                               '
    WRITE(lfn,'(A137)')' OUTPUT             WHU Phase OSB solutions in SINEX format                                                                              '
    WRITE(lfn,'(A137)')' CONTACT            ShengYi Xu, shengyixu@whu.edu.cn                                                                                     '
    WRITE(lfn,'(A137)')'-FILE/REFERENCE                                                                                                                          '
    WRITE(lfn,'(A137)')'*-------------------------------------------------------------------------------------------------------                                 '
    WRITE(lfn,'(A137)')'+BIAS/DESCRIPTION                                                                                                                        '
    WRITE(lfn,'(A48,F6.1,A)')' PARAMETER_SPACING                            ',CKF.dintv,'                                                                         '
    WRITE(lfn,'(A137)')'-BIAS/DESCRIPTION                                                                                                                        '
    WRITE(lfn,'(A137)')'*-------------------------------------------------------------------------------------------------------                                 '
    WRITE(lfn,'(A137)')'+BIAS/SOLUTION                                                                                                                           '
    WRITE(lfn,'(A137)')'*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___'

    RETURN
END SUBROUTINE
