!
!   CREATED BY SHENYIXU, 2023-07-25
!   PURPOSE: 
!       1. Considering TEC rate
!       2. [A new automated cycle slip detection and repair method for a single dual-frequency GPS receiver], Zhizhao Liu
!*
!
SUBROUTINE ppp_cmpt_lglw(CKF,OB,SIT,SAT,iepo)
!
!*
USE const
USE ckdctrl
USE observation
USE station
USE satellite
USE iso_fortran_env
IMPLICIT NONE

!*
! Start the exectuable
!!------------------------
TYPE(ckdcfg) :: CKF
TYPE(RNXOBS) :: OB(1:*)
TYPE(SITE) :: SIT(1:*)
TYPE(SATE) :: SAT(MAXSAT)
INTEGER(IT) :: iepo

    !*
    ! The local variables
    !!------------------------
    LOGICAL(LG) :: lfirst
    INTEGER(IT) :: isit,isat,isys,ljd(MAXSAT,MAXSIT),interrupt_flag(MAXSIT),iy,imon,id,ih,im
    REAL(RL) :: ltm(MAXSAT,MAXSIT),lsg(MAXSAT,MAXSIT),lsw(MAXSAT,MAXSIT),lsewl(MAXSAT,MAXSIT),lseewl(MAXSAT,MAXSIT),lshewl(MAXSAT,MAXSIT),isec
    REAL(RL) :: lcode(MAXSAT,MAXSIT),lphase(MAXSAT,MAXSIT),ldoppler(MAXSAT,MAXSIT)

    REAL(RL) :: Elw(MAXSAT,MAXSIT),SIGlw(MAXSAT,MAXSIT),detI
    REAL(RL) :: TEC_0(MAXSAT,MAXSIT),TECR_0(MAXSAT,MAXSIT)
    REAL(RL) :: tec_1(MAXSAT,MAXSIT),tec_2(MAXSAT,MAXSIT),tec_3(MAXSAT,MAXSIT)

    REAL(RL) :: dt,lgr,lw,lew,leewl,lhewl,wwl,lp,lc

    DATA lfirst /.TRUE./
    SAVE lfirst,ljd,ltm,lsg,lsw,lsewl,lseewl,lshewl,lcode,lphase,tec_1,tec_2,tec_3

    !*
    ! The function called
    !!------------------------
    REAL(RL) :: timdif

    !*
    ! Start the exectuable code
    !!----------------------------
    CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,isec)
    WRITE(1006,'((A),I5,4I3,F11.7,I7,F10.2,(A))')'TIM ',iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,'  SIT-SAT-LG-DeltI-LG_TECR-LW'

    IF (MOD(INT((CKF%mjd-CKF%mjd0)*86400.0+CKF%sod-CKF%sod0),CKF%ReConvTime) .EQ. 0) THEN
        lfirst=.TRUE.
    END IF

    IF (lfirst) THEN
        lfirst=.FALSE.
        DO isit=1,MAXSIT
            DO isat=1, CKF%nprn
                ljd(isat,isit)=0
                ltm(isat,isit)=0.d0

                lsg(isat,isit)=0.d0
                lsw(isat,isit)=0.d0
                lsewl(isat,isit)=0.d0
                lseewl(isat,isit)=0.d0
                lshewl(isat,isit)=0.d0

                lcode(isat,isit)=0.d0
                lphase(isat,isit)=0.d0
                ldoppler(isat,isit)=0.d0

                tec_1(isat,isit)=0.d0
                tec_2(isat,isit)=0.d0
                tec_3(isat,isit)=0.d0
            END DO            
        END DO
    END IF
    interrupt_flag=0
    DO isit=1, CKF%nsit 
        IF (COUNT(OB(isit)%obs(1:CKF%nprn,MAXFREQ+1).NE.0.d0) .LT. 4) THEN
        interrupt_flag(isit)=1
        END IF
    END DO

    !! real-time preprocessing
    DO isat=1,CKF%nprn
        isys=INDEX(SYS,CKF%cprn(isat)(1:1))

        DO isit=1,CKF%nsit
            IF (interrupt_flag(isit).EQ.1) CYCLE
            IF (OB(isit)%obs(isat,1).EQ.0.d0)CYCLE
            
            !! gap is so large that a new ambiguity has to be set
            IF (ljd(isat,isit).NE.0) THEN
                dt=timdif(CKF%mjd,CKF%sod,ljd(isat,isit),ltm(isat,isit))
                IF (dt.GT.CKF%gap) OB(isit)%flag(isat,1:MAXFREQ)=1
            END IF

            detI=0.d0
            !! compute the TEC, TECR and TECA
            IF (CKF%nfreq(isys).GE.2) THEN
                TEC_0(isat,isit)=SAT(isat)%freq(1)**2/(40.3d16*(SAT(isat)%g2-1))*(SAT(isat)%lamda(1)*OB(isit)%obs(isat,1)-SAT(isat)%lamda(2)*OB(isit)%obs(isat,2))
                IF (iepo.EQ.1) THEN
                    tec_1(isat,isit)=TEC_0(isat,isit)
                    tec_2(isat,isit)=0.d0
                    tec_3(isat,isit)=0.d0
                ELSE IF (iepo.EQ.2) THEN
                    tec_2(isat,isit)=tec_1(isat,isit)
                    tec_1(isat,isit)=TEC_0(isat,isit)
                    tec_3(isat,isit)=0.d0                   
                ELSE IF (iepo.EQ.3) THEN
                    tec_3(isat,isit)=tec_2(isat,isit)
                    tec_2(isat,isit)=tec_1(isat,isit)
                    tec_1(isat,isit)=TEC_0(isat,isit)
                ELSE
                    TECR_0(isat,isit)=(2*tec_1(isat,isit)-3*tec_2(isat,isit)+tec_3(isat,isit))/CKF%dintv
                    detI=(40.3d16*(SAT(isat)%g2-1))*CKF%dintv/SAT(isat)%freq(1)**2 * TECR_0(isat,isit)

                    tec_3(isat,isit)=tec_2(isat,isit)
                    tec_2(isat,isit)=tec_1(isat,isit)
                    tec_1(isat,isit)=TEC_0(isat,isit)                 
                END IF

                lgr=SAT(isat).lamda(1)*OB(isit).obs(isat,1)-SAT(isat).lamda(2)*OB(isit).obs(isat,2)
                lw=OB(isit).obs(isat,1)-OB(isit).obs(isat,2)-&
                    (SAT(isat).g*OB(isit).obs(isat,MAXFREQ+1)+OB(isit).obs(isat,MAXFREQ+2))/(1.d0+SAT(isat).g)/SAT(isat).lamdw
            END IF

            IF (CKF%nfreq(CKF%iref).GE.2) THEN
                IF (ljd(isat,isit).NE.0)THEN
                    !! only 3 epochs is used, maybe longer, and need to repair cycle slip
                    IF (iepo.LE.3) THEN
                        IF(DABS(lgr-lsg(isat,isit)).GT.CKF%lg*(SAT(isat).lamda(2)-SAT(isat).lamda(1)) .OR. DABS(lw-lsw(isat,isit)).GT.CKF.lw) THEN
                            OB(isit)%flag(isat,1:MAXFREQ)=1
                            WRITE(*,'((A),4F10.5)') ' |->... new amb (cycle slip): '//CKF%cprn(isat),-lgr+lsg(isat,isit),detI,detI-lgr+lsg(isat,isit),DABS(lw-lsw(isat,isit))                       
                            WRITE(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ->... new amb (cycle slip): '//CKF%cprn(isat)
                        END IF
                    ELSE
                        ! 0.2 TECU
                        IF(DABS(detI-lgr+lsg(isat,isit)).GT.0.20*0.16 .OR. DABS(lw-lsw(isat,isit)).GT.CKF.lw) THEN
                            OB(isit)%flag(isat,1:MAXFREQ)=1
                            WRITE(*,'((A),4F10.5)') ' |->... new amb (cycle slip): '//CKF%cprn(isat),-lgr+lsg(isat,isit),detI,detI-lgr+lsg(isat,isit),DABS(lw-lsw(isat,isit))
                            WRITE(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ->... new amb (cycle slip): '//CKF%cprn(isat)
                        END IF
                    END IF
                ELSE
                    OB(isit)%flag(isat,1:MAXFREQ)=1
                    WRITE(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' +>... new amb (new sat): '//CKF%cprn(isat)
                END IF
            END IF

            !! reverse obs
            ljd(isat,isit)=CKF%mjd
            ltm(isat,isit)=CKF%sod

            lsg(isat,isit)=lgr
            lsw(isat,isit)=lw
            lsewl(isat,isit)=lew
            lseewl(isat,isit)=leewl
            lshewl(isat,isit)=lhewl

            lphase(isat,isit)=OB(isit)%obs(isat,1)
            lcode(isat,isit)=OB(isit)%obs(isat,1+MAXFREQ)
            ldoppler(isat,isit)=OB(isit)%obs(isat,2*MAXFREQ+1)            
        END DO
    END DO

    RETURN
END SUBROUTINE