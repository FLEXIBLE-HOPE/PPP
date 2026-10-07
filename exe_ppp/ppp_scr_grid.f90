!! Identify biases using a grid
!! Jianghui Geng
!! Feb 26 2012
!!
!*
SUBROUTINE ppp_scr_grid(isit,CKF,SIT,OB,OBD,SAT)
!!
!*
    USE const
    USE ckdctrl
    USE station
    USE satellite
    USE observation
    USE ISO_FORTRAN_ENV
    USE info
    IMPLICIT NONE

    !*
    ! The arguments
    !!----------------------
    INTEGER(IT) :: isit
    TYPE(CKDCFG) :: CKF
    TYPE(SITE) :: SIT
    TYPE(RNXOBS) :: OB
    TYPE(RNXOBS) :: OBD
    TYPE(SATE) :: SAT(MAXSAT)

    !*
    ! The local variables
    !!---------------------------
    INTEGER(IT) :: pointer_string
    INTEGER(IT) :: i,k,kobs,isat,hdel,ndel,flg(MAXSAT),jst(MAXSAT),isys,ipar,ipt,num
    REAL(RL) :: avg,std,sig,omc(MAXSAT),phase(MAXSIT),wgt(MAXSAT)

    !*
    ! Start the exectuable code
    !!---------------------------

    !! fill in omc
    omc=0.d0
    DO isat=1,CKF.nprn
        IF (OBD.omc(isat,1).NE.0.d0 .AND. OBD.flag(isat,1).NE.1) THEN
            omc(isat)=SAT(isat).fac(1)*OBD.omc(isat,1)-SAT(isat).fac(2)*OBD.omc(isat,2)
        END IF
        !!the new ambiguity parameters
        IF (OBD.flag(isat,1) .NE. 0) THEN
            omc(isat)=0.d0
        END IF
    END DO
    IF (ALL(omc .EQ. 0.d0)) RETURN

    DO isys=1, CKF.nsys
        wgt=1.d0
        flg=2
        DO isat=1,CKF.nprn
            IF (omc(isat).NE.0.d0 .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(isys:isys)) THEN 
                ! print *,CKF%cprn(isat),omc(isat)
                flg(isat)=0
            END IF
        END DO
        IF (COUNT(flg(1:CKF.nprn).EQ.0) .GT. 0) THEN
            num=COUNT(flg(1:CKF.nprn).EQ.0)
            CALL get_wgt_mean(.TRUE.,omc(1),flg,wgt,CKF.nprn,ndel,avg,std,sig)
            DO WHILE(num-ndel.GT.2 .AND. std.GT.0.03d0)
                hdel=ndel
                CALL sign_robust(CKF.nprn,omc(1),flg,0.015d0,ndel)
                IF (ndel .EQ. hdel) EXIT
                CALL get_wgt_mean(.FALSE.,omc(1),flg,wgt,CKF.nprn,ndel,avg,std,sig)
            END DO
            IF (ndel .GT. 0) THEN
                DO isat=1,CKF.nprn
                    IF (omc(isat).NE.0.d0 .AND. CKF.cprn(isat)(1:1).EQ.CKF.system(isys:isys)) THEN
                        IF (flg(isat) .NE. 0) THEN
                            OB%flag(isat,:) = 1
                            IF (flg(isat) .EQ. 2) THEN
                                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4)') ' ... new amb (range_wgt)_grid ',&
                                    CKF%mjd,CKF%sod,SIT%name,CKF%cprn(isat),omc(isat)-avg,avg,std
                            ELSE
                                WRITE(OUTPUT_UNIT,'(A,I5,F9.1,1X,A4,1X,A3,3F14.4)') ' ... new amb (range_sign)_grid ',&
                                    CKF%mjd,CKF%sod,SIT%name,CKF%cprn(isat),omc(isat)-avg,avg,std
                            END IF                            
                        END IF
                    END IF
                END DO                
            END IF
        END IF        
    END DO

    RETURN

END SUBROUTINE
