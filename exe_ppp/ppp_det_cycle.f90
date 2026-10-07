!
!   CREATED BY SHENYIXU, 2023-07-29
!   PURPOSE: 
!       1. Detect cycle slip for highrate
!       2. Considering the TECR change, suitable for continous cycle-slip
!*
! GF maybe change to double-difference for LEO
SUBROUTINE ppp_det_cycle(CKF,OB,SAT)
!!
!*
USE const
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! Start the exectuable
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB(1:*)
TYPE(SATE) :: SAT(MAXSAT)

    !*
    ! The local variables
    !!------------------------
    LOGICAL(LG) :: lfirst
    INTEGER(IT) :: isit,isat,isys,ljd(MAXSAT,MAXSIT), iy,imon,id,ih,im,interrupt_flag(MAXSIT),inlg,inlw
    REAL(RL) :: ltm(MAXSAT,MAXSIT),lg_2(MAXSAT,MAXSIT),lg_1(MAXSAT,MAXSIT),lsw(MAXSAT,MAXSIT),isec,teca

    REAL(RL) :: dt,lgr,lw

    DATA lfirst /.TRUE./
    SAVE lfirst,ljd,ltm,lg_1,lg_2,lsw

    !*
    ! The function called
    !!------------------------
    REAL(RL) :: timdif


    !*
    ! Start the exectuable code
    !!----------------------------
    
    CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,isec)
    WRITE(1006,'((A),I5,4I3,F11.7,I7,F10.2,(A))')'TIM ',iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,'  SIT-SAT-LG-LW-LWL-LEWL-LEEWL-LHEWL'
    IF (MOD(INT((CKF.mjd-CKF.mjd0)*86400.0+CKF.sod-CKF.sod0),CKF.ReConvTime) .EQ. 0) THEN
        lfirst=.TRUE.
    END IF

    IF (lfirst) THEN
        lfirst=.FALSE.
        DO isit=1, MAXSIT
        DO isat=1, CKF.nprn
            ljd(isat,isit)=0
            ltm(isat,isit)=0.d0
            lg_1(isat,isit)=0.d0
            lg_2(isat,isit)=0.d0
            lsw(isat,isit)=0.d0
        END DO
        END DO
    END IF

    interrupt_flag=0
    !!xsy: 如果当前历元有效卫星数小于4，那么继续保存上一历元的LG、LW、LEW,否则会对周条探测和修复造成误区
    DO isit=1, CKF.nsit 
        IF (COUNT(OB(isit).obs(1:CKF.nprn,MAXFREQ+1).NE.0.d0) .LT. 4) THEN
        interrupt_flag(isit)=1
        END IF
    END DO 

    !! real-time preprocessing
    DO isat=1, CKF.nprn

        isys=INDEX(SYS,CKF.cprn(isat)(1:1))

        DO isit=1, CKF.nsit
            IF (interrupt_flag(isit) .EQ. 1) CYCLE
            IF (OB(isit).obs(isat,1) .EQ. 0.d0) CYCLE
        
            !! gap is so large that a new ambiguity has to be set
            IF (ljd(isat,isit) .NE. 0) THEN
                dt=timdif(CKF.mjd,CKF.sod,ljd(isat,isit),ltm(isat,isit))
                !! for LEO, should be shorter
                IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
                    IF (dt .GT. 60) THEN
                        OB(isit).flag(isat,1:MAXFREQ)=1
                        write(*,'((A))') ' ... A big gap for : '//ckf.cprn(isat)
                        write(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ... A big gap for : '//ckf.cprn(isat)
                    END IF
                ELSE
                    IF (dt .GT. CKF.gap) THEN
                        OB(isit).flag(isat,1:MAXFREQ)=1
                        write(*,'((A))') ' ... A big gap for : '//ckf.cprn(isat)
                        write(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ... A big gap for : '//ckf.cprn(isat)
                    END IF
                END IF
            END IF

            IF (CKF.nfreq(isys) .EQ. 2) THEN
                !! LG ionosphere observation, widelane ambiguity
                lgr=(SAT(isat).lamda(1)*OB(isit).obs(isat,1)-SAT(isat).lamda(2)*OB(isit).obs(isat,2))/ &
                                    (SAT(isat).lamda(2)-SAT(isat).lamda(1))

                lw=OB(isit).obs(isat,1)-OB(isit).obs(isat,2)-&
                        (SAT(isat).g*OB(isit).obs(isat,MAXFREQ+1)+OB(isit).obs(isat,MAXFREQ+2))/(1.d0+SAT(isat).g)/SAT(isat).lamdw
            END IF


            IF(CKF.nfreq(CKF.iref) .GE. 2)THEN
                !! new ambiguity
                IF(ljd(isat,isit) .NE. 0) THEN
                    !2023-07-23 by xsy: 手动添加周跳
                    ! IF(CKF%sod.EQ.38180.d0 .AND. CKF%cprn(isat)(1:1).EQ.'L')THEN
                    !    WRITE(*,*)'L05 38180.00'
                    !    OB(isit).flag(isat,1:MAXFREQ)=1
                    ! END IF
                    IF(CKF%sod.EQ.37986.d0 .AND. CKF%cprn(isat)(1:1).EQ.'L')THEN
                        WRITE(*,*)' +> L05 37986.00'
                        OB(isit).flag(isat,1:MAXFREQ)=1
                    END IF
                    IF(CKF%sod.EQ.37992.d0 .AND. CKF%cprn(isat)(1:1).EQ.'L')THEN
                        WRITE(*,*)' +> L05 37992.00'
                        OB(isit).flag(isat,1:MAXFREQ)=1
                    END IF

                    WRITE(1006,'(I4,A5,6F9.3)')isit,CKF.cprn(isat),(lgr-lg_1(isat,isit)),(lw-lsw(isat,isit)),&
                                                (lgr-2*lg_1(isat,isit)+lg_2(isat,isit))*0.5d0,0.d0,0.d0,OB(isit)%elev(isat)*RAD2DEG                    
                    inlg=CKF%lg
                    IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
                      !TECR阈值大致为0.15TECU/s即2.4cm/s -> dintv=2 -> 4.8cm -> 1.2cycle,但是低轨卫星运动快,阈值或许需要适当提高
                      inlg=3.0
                    END IF
                    inlw=CKF%lw
                    ! IF (CKF%cprn(isat)(1:1).EQ.'L') THEN
                    !   !L05伪距精度优于0.5m,取4倍rms为2.0m => MW噪声大约为伪距的0.3倍 => 因此阈值可以给0.6,甚至可以稍微放大
                    !   inlw=0.6
                    ! END IF

                    !! The first detection
                    IF(DABS(lgr-lg_1(isat,isit)).GT.inlg .OR. DABS(lw-lsw(isat,isit)).GT.inlw) THEN
                        OB(isit).flag(isat,1:MAXFREQ)=1
                        write(*,'((A))') ' ... new amb (cycle slip): '//ckf.cprn(isat)
                        write(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ... new amb (cycle slip): '//ckf.cprn(isat)
                    END IF
                    !! The second detection
                    IF(lg_2(isat,isit).NE.0.d0 .AND. OB(isit)%flag(isat,1).NE.1) THEN
                        IF(CKF%cprn(isat)(1:1).EQ.'L') THEN
                            teca=0.6
                        ELSE
                            teca=inlg*0.5
                        END IF
                        IF(DABS((lgr-2*lg_1(isat,isit)+lg_2(isat,isit))*0.5d0) .GE.teca) THEN
                            OB(isit).flag(isat,1:MAXFREQ)=1
                            write(*,'((A))') ' ... new amb (cycle slip): '//ckf.cprn(isat)
                            write(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' ... new amb (cycle slip): '//ckf.cprn(isat)
                        END IF
                    END IF
                ELSE
                    !! first ambiguity for this satellite
                    OB(isit).flag(isat,1:MAXFREQ)=1
                    write(1001,'(I5,4I3,F11.7,I7,F10.2,(A))') iy,imon,id,ih,im,isec,ckf.mjd,ckf.sod,' +> new amb (new sat): '//ckf.cprn(isat)
                END IF
            END IF

            !! reserve obs
            ljd(isat,isit)=CKF.mjd
            ltm(isat,isit)=CKF.sod

            lg_2(isat,isit)=lg_1(isat,isit)
            lg_1(isat,isit)=lgr

            lsw(isat,isit)=lw
        
        END DO
    END DO
END SUBROUTINE
