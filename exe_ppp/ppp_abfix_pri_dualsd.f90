!*
! CREATED BY ShengYi XU
! DATE: 2023-11-06
! PURPOSE: Make SD ambiguity for IF and Dual-frequency including the WL(MW) and NL(IF map to N1)
!*
SUBROUTINE ppp_abfix_pri_dualsd(CKF,isit,OB,SIT,SAT,UPD,NM,PM,AM,QM,SL)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE ambiguity
USE observation
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
INTEGER(IT) :: isit
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT
TYPE(SATE) :: SAT(MAXSAT)
TYPE(FCB) :: UPD
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
TYPE(AMBT) :: AM(1:*)
TYPE(INVM) :: QM
TYPE(SOL) :: SL

    !*
    ! The local variables
    !!------------------------
    LOGICAL(LG) :: lfirst(MAXSIT)
    INTEGER(IT) :: kpt(MAXSYS),nxl(MAXSYS),refsat_AM(MAXSYS),num_wl
    INTEGER(IT) :: i,j,isys,isat,ifg(MAXSAT),ipt(MAXSAT),maxtim,ind,iamb,ndl,sd_num,k
    REAL(RL) :: maxele,dump,sd_wl(MAXSAT),sd_sigw(MAXSAT),fxl(MAXSYS,MAXSIT),wgt(MAXSAT)
    REAL(RL),ALLOCATABLE :: map(:,:),mapt(:,:),maptx(:,:)
    REAL(RL) :: rwl,wwl,rxl,pxl,vxl,sxl,alpha
    INTEGER(IT) :: irectype,nfix(MAXSYS)
    TYPE(AMBD), POINTER :: AB(:)

    !*
    ! The function called
    !!----------------------------
    INTEGER(IT) :: pointer_string

    DATA lfirst/MAXSIT*.TRUE./
    SAVE lfirst,fxl
    !*
    ! Start the exectuable code
    !!------------------------
    ! Map IFamb to N1
    ALLOCATE(map(QM%ntot,QM%ntot))
    ALLOCATE(mapt(QM%ntot,QM%ntot))
    ALLOCATE(maptx(QM%ntot,QM%ntot))

    irectype = pointer_string(MAXRECTYPE, UPD%rectype, TRIM(SIT%rectyp))
    num_wl = 0

    IF (lfirst(isit) .EQ. .TRUE.) then
        lfirst(isit)=.FALSE.
        fxl(:,isit)=10.d0
    END IF

    WRITE(1004,'((A),I7,F10.2,(A))') 'TIM',CKF%mjd,CKF%sod,'-------------------------------------------------------------------------------------------> UDIF'
    WRITE(1004,'(A)')'---PRN---------------ELE------------FLT------------SIG------------FIX------------UPD----------rBIAS-------------CMT'
    
    ! STEP 1. choose reference satellite based on elevation for ammbfix
    map = 0.d0
    mapt = 0.d0
    DO i=1, QM%ntot
        map(i,i)=1.d0
        mapt(i,i)=1.d0
    END DO
    DO isys=1, CKF%nsys
        nxl(isys) = 0
        kpt(isys) = 0
        ipt = 0
        DO i=1, NM%ns
            isat = AM(i)%psat
            IF(CKF%cprn(isat)(1:1) .NE. CKF%system(isys:isys)) CYCLE
            IF(OB%omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE
            IF(TRIM(AM(i)%pname) .NE. 'AMBL1') CYCLE
            ! For IF and Dual-frequency
            IF (CKF%nfq(CKF%iref) .EQ. 1) THEN
                IF (CKF%ArBiasMode .EQ. 'FCB') THEN
                    IF (AM(i)%iobs.EQ.0 .OR. UPD%wfcb(isat).EQ.10.d0 .OR. UPD%nfcb(isat).EQ.10.d0 .OR. AM(i)%elev/AM(i)%iobs.LE.CKF%cutoff) CYCLE
                ELSE IF (CKF%ArBiasMode .EQ. 'OSB') THEN
                    IF (AM(i)%iobs.EQ.0 .OR. INDEX(CKF%nofixsat,CKF%cprn(isat)).NE.0 .OR. AM(i)%elev/AM(i)%iobs.LE.CKF%cutoff) CYCLE
                END IF
            END IF
            nxl(isys) = nxl(isys)+1
            ifg(nxl(isys)) = 0
            ipt(nxl(isys)) = i
        END DO

        maxtim = 0
        maxele = 0.d0
        
        IF(nxl(isys) .GT. 0)THEN
            ! DO i=1, nxl(isys)
            !     IF(ifg(i) .NE. 0) CYCLE
            !     IF(AM(ipt(i))%iobs .GT. maxtim)THEN
            !         kpt(isys) = ipt(i)
            !         maxtim = AM(ipt(i))%iobs
            !         maxele = AM(ipt(i))%elev/AM(ipt(i))%iobs
            !     ELSE IF(AM(ipt(i))%iobs .EQ. maxtim .AND. AM(ipt(i))%elev/AM(ipt(i))%iobs .GT. maxele)THEN
            !         kpt(isys) = ipt(i)
            !         maxele = AM(ipt(i))%elev/AM(ipt(i))%iobs
            !     END IF
            ! END DO
            DO i=1, nxl(isys)
                IF(AM(ipt(i))%elev/AM(ipt(i))%iobs .GT. maxele)THEN
                    kpt(isys) = ipt(i)
                    maxele = AM(ipt(i))%elev/AM(ipt(i))%iobs
                END IF
            END DO            
            IF(kpt(isys) .NE. 0)THEN
                isat = AM(kpt(isys))%psat
                CKF%refnprn(isys) = isat
                CKF%refcprn(isys) = CKF%cprn(isat)
                ! WRITE(*,'((A),2I4,A6,I4,3X,(A),F5.2)')'REFSAT: ',isys,CKF%refnprn(isys),CKF%refcprn(isys),kpt(isys),AM(kpt(isys))%pname,maxele
            END IF
        ELSE
            CYCLE
        END IF

        ! STEP 2. MAP udIFamb to sdIFamb
        IF (CKF%nfq(CKF%iref) .EQ. 1) THEN
            DO i=QM%nxyz+1, QM%ntot
                isat = AM(i-QM%nxyz)%psat
                IF(CKF%cprn(isat)(1:1) .NE. CKF%system(isys:isys)) CYCLE
                IF(TRIM(AM(i-QM%nxyz)%pname) .NE. 'AMBL1') CYCLE
                IF(isat .EQ. CKF%refnprn(isys)) THEN
                    ! IF (CKF%llog .EQ. .FALSE.) THEN
                    !     IF (AM(i-QM%nxyz)%ptime(2).EQ.CKF%mjd+CKF%sod/86400.d0) THEN
                    !         map(i,i) = 0.d0
                    !         mapt(i,i)= 0.d0
                    !     END IF
                    ! ELSE
                    !     map(i,i) = 0.d0
                    !     mapt(i,i)= 0.d0
                    ! END IF
                    map(i,i) = 0.d0
                    mapt(i,i)= 0.d0                    
                ELSE
                    map(i                 ,i                ) =  1.d0                             !LC amb
                    map(i                 ,kpt(isys)+QM%nxyz) = -1.d0
                    mapt(i                ,i                ) =  1.d0
                    mapt(kpt(isys)+QM%nxyz,i                ) = -1.d0
                    AM(i-QM%nxyz)%famb = AM(i-QM%nxyz)%famb - AM(kpt(isys))%famb
                END IF
            END DO
            DO i=QM%nxyz+1, QM%ntot
                isat = AM(i-QM%nxyz)%psat
                IF (CKF%cprn(isat)(1:1) .NE. CKF%system(isys:isys)) CYCLE
                IF(TRIM(AM(i-QM%nxyz)%pname).NE.'AMBL1') CYCLE
                !参考星模糊度需要消掉
                IF(isat.EQ.CKF%refnprn(isys))THEN
                    AM(i  -QM.nxyz)%famb= 0.d0
                END IF
            END DO            
        END IF
    END DO
    DO i=1, QM%ntot
        DO j=i+1, QM%ntot
            QM%invx(i,j) = QM%invx(j,i)
        END DO
    END DO
    CALL matmpy(map,QM%invx(1:QM%ntot,1:QM%ntot),maptx, QM%ntot,QM%ntot,QM%ntot)
    CALL matmpy(maptx,mapt,QM%invx(1:QM%ntot,1:QM%ntot),QM%ntot,QM%ntot,QM%ntot)

    ! STEP 3. UD WL(MW)
    iamb = pointer_string(OB%npar,OB%pname,'AMBL1')
    DO isat=1, CKF%nprn
        IF (OB%omc(isat,1).NE.0.d0 .AND. OB%omc(isat,MAXFREQ+1).NE.0.d0) THEN
            ind = OB%ltog(iamb,isat) - NM%npc
            isys = INDEX(CKF%system,CKF%cprn(isat)(1:1))
            IF (isat .EQ. CKF%refnprn(isys)) THEN
                refsat_AM(isys) = ind
            END IF
            ! Correct PCO for MW
            IF (CKF%if_pco_corr) THEN
                ! xsy: only the Z PCO of satellite
                ! rwl=(OB%obs(isat,1) + SAT(isat)%xyz(3,1)*SAT(isat)%freq(1)/VEL_LIGHT)- &
                !     (OB%obs(isat,2) + SAT(isat)%xyz(3,2)*SAT(isat)%freq(2)/VEL_LIGHT)- &
                !     (SAT(isat)%g*(OB%obs(isat,MAXFREQ+1)+SAT(isat)%xyz(3,1)) + (OB%obs(isat,MAXFREQ+2)+SAT(isat)%xyz(3,2)))/(1.d0+SAT(isat)%g)/SAT(isat)%lamdw
                ! xsy: the Z PCO of satellite and SIT PCO
                ! rwl=(OB%obs(isat,1) + (SIT%enu(3,1,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,1))*SAT(isat)%freq(1)/VEL_LIGHT)- &
                !     (OB%obs(isat,2) + (SIT%enu(3,2,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,2))*SAT(isat)%freq(2)/VEL_LIGHT)- &
                !     (SAT(isat)%g*(OB%obs(isat,MAXFREQ+1)+SIT%enu(3,1,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,1)) + (OB%obs(isat,MAXFREQ+2)+SIT%enu(3,2,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,2)))/(1.d0+SAT(isat)%g)/SAT(isat)%lamdw
                rwl=(OB%obs(isat,1) + (SIT%enu(3,1,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,1)*dcos(OB%nadir(isat)))*SAT(isat)%freq(1)/VEL_LIGHT)- &
                    (OB%obs(isat,2) + (SIT%enu(3,2,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,2)*dcos(OB%nadir(isat)))*SAT(isat)%freq(2)/VEL_LIGHT)- &
                    (SAT(isat)%g*(OB%obs(isat,MAXFREQ+1)+SIT%enu(3,1,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,1)*dcos(OB.nadir(isat))) + (OB%obs(isat,MAXFREQ+2)+SIT%enu(3,2,isys)*dsin(OB%elev(isat)) + SAT(isat)%xyz(3,2)*dcos(OB.nadir(isat))))/(1.d0+SAT(isat)%g)/SAT(isat)%lamdw                    
            ELSE
                rwl=OB%obs(isat,1)-OB%obs(isat,2)- &
                    (SAT(isat)%g*OB%obs(isat,MAXFREQ+1)+OB%obs(isat,MAXFREQ+2))/(1.d0+SAT(isat)%g)/SAT(isat)%lamdw
            END IF

            ! Correct WL FCB
            IF (CKF%ArBiasMode .EQ. 'FCB') THEN
                IF (UPD%wfcb(isat) .NE. 10.d0) THEN
                    rwl = rwl - UPD%wfcb(isat)
                END IF
            END IF
            ! Correct the FPA(receiver-specific) bias for WL
            IF (irectype .NE. 0) THEN
                rwl = rwl - UPD%recfcb(isat,2,irectype)
            END IF
            ! Weight mean WL
            wwl = 1.d0
            IF (OB%elev(isat)*RAD2DEG .LE. 30.d0) wwl = wwl*2.d0*DSIN(OB%elev(isat))
            IF (OB%flag(isat,1) .NE. 0) THEN
                AM(ind)%abin = NINT(rwl)                ! a priori value, must be integer
            END IF
            rwl = rwl-AM(ind)%abin
            AM(ind)%xrwl = AM(ind)%xrwl+wwl*rwl       ! FPA mean
            AM(ind)%weig = AM(ind)%weig+wwl           ! FPA weight
            AM(ind)%xswl = AM(ind)%xswl+wwl*rwl**2    ! FPA sigma

            AM(ind)%abwl = AM(ind)%xrwl/AM(ind)%weig + AM(ind)%abin ! wl

            dump = AM(ind)%xrwl/AM(ind)%weig
            dump = AM(ind)%xswl - AM(ind)%weig*dump**2
            AM(ind)%sigw = DSQRT(ABS(dump)/AM(ind)%iobs/AM(ind)%weig)    ! wl_sigma
            ! WRITE(*,'(A4,2F10.3)')CKF%cprn(isat),AM(ind)%abwl,AM(ind)%sigw
        END IF
    END DO

    ! STEP 4. FIX the SD WL(MW)
    IF (CKF%liar) THEN
        SL%fixnum_wl = 0

        DO isys=1, CKF%nsys
            sd_num = 0
            ipt=-1
            DO ind=1, NM.ns
                isat=AM(ind).psat
                IF (CKF%cprn(isat)(1:1) .NE. CKF%system(isys:isys)) CYCLE
                IF (OB%omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE
                IF (CKF%llog .EQ. .FALSE.) THEN
                    IF (AM(ind)%ptime(2) .NE. CKF%mjd+CKF%sod/86400.d0) CYCLE                
                END IF
                IF (isat .EQ. CKF%refnprn(isys)) THEN
                    IF (irectype .NE. 0) THEN
                        WRITE(1004,'((A),A4,(A),F10.3,5F15.3,4X,(A))')'WL',CKF%cprn(isat),'--------',OB%elev(isat)*RAD2DEG,AM(refsat_AM(isys))%abwl,AM(refsat_AM(isys))%sigw,UPD%wfcb(isat),0.d0,UPD%recfcb(isat,2,irectype),'-------> REF_SAT'
                    ELSE
                        WRITE(1004,'((A),A4,(A),F10.3,5F15.3,4X,(A))')'WL',CKF%cprn(isat),'--------',OB%elev(isat)*RAD2DEG,AM(refsat_AM(isys))%abwl,AM(refsat_AM(isys))%sigw,UPD%wfcb(isat),0.d0,0.d0,'-------> REF_SAT'
                    END IF
                    CYCLE
                END IF
                ! IF ((AM(ind)%ptime(2)-AM(ind)%ptime(1))*86400.d0 .LT. CKF%minsec_common) CYCLE
                IF (CKF%ArBiasMode .EQ. 'FCB') THEN
                    IF (UPD%wfcb(isat) .EQ. 10.d0) THEN
                        WRITE(1004,'((A),A4,A)')'WL',CKF%cprn(isat),'---------NONE WL FCB'
                        CYCLE
                    END IF
                END IF

                sd_num = sd_num + 1
                sd_wl(sd_num) = AM(ind)%abwl - AM(refsat_AM(isys))%abwl
                ifg(sd_num) = 0
                wgt(sd_num) = 1.d0
                ipt(ind) = sd_num
            END DO
            IF (sd_num .EQ. 0) CYCLE

            CALL proc_xl_dirc(sd_num,sd_wl,wgt,ifg,ndl,fxl(isys,isit),vxl,sxl)

            ! IF ((sd_num-ndl)*1.d0/sd_num.GE.0.6d0 .AND. vxl.LE.0.25d0) THEN
            IF ((sd_num-ndl)*1.d0/sd_num.GE.0.6d0 .AND. vxl.LE.0.30d0) THEN
        
                DO ind=1, NM%ns
                    isat = AM(ind)%psat
                    IF (OB%omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE
                    IF (CKF%llog .EQ. .FALSE.) THEN
                        IF (AM(ind)%ptime(2) .NE. CKF%mjd+CKF%sod/86400.d0) CYCLE                
                    END IF                    
                    IF (ipt(ind) .EQ. -1) CYCLE

                    rwl = AM(ind)%abwl - AM(refsat_AM(isys))%abwl - fxl(isys,isit)
                    pxl = DSQRT(AM(ind)%sigw**2 + AM(refsat_AM(isys))%sigw**2)

                    CALL prob_resol(rwl,pxl,1,CKF%wl_maxdev,CKF%wl_maxsig,alpha)
                    IF (CKF%lamb .EQ. .TRUE.) THEN
                        IF (AM(ind)%ifab .LT. 2) THEN
                            IF (alpha .GT. CKF%wl_alpha) THEN
                                AM(ind)%ifab = 2
                                AM(ind)%abwl = NINT(rwl)
                                AM(ind)%fwl  = NINT(rwl)
                            ELSE
                                AM(ind)%abwl = rwl
                                AM(ind)%fwl  = 0.d0
                            END IF
                        END IF
                    ELSE
                        IF (alpha .GT. CKF%wl_alpha) THEN
                            AM(ind)%ifab = 2
                            AM(ind)%abwl = NINT(rwl)
                            AM(ind)%fwl  = NINT(rwl)
                            SL%fixnum_wl=SL%fixnum_wl+1
                        ELSE
                            AM(ind)%abwl = rwl
                            AM(ind)%fwl  = 0.d0
                        END IF
                    END IF
                    
                    IF (AM(ind)%ifab .GE. 2) THEN
                        num_wl = num_wl + 1
                        IF (irectype .NE. 0) THEN
                            WRITE(1004,'((A),A4,2I4,F10.3,5F15.3,5X,(A))')'WL',CKF%cprn(isat),ind,num_wl,OB%elev(isat)*RAD2DEG,rwl,pxl,AM(ind)%abwl,UPD%wfcb(isat),UPD%recfcb(isat,2,irectype),'FIX_WL'
                        ELSE
                            WRITE(1004,'((A),A4,2I4,F10.3,5F15.3,5X,(A))')'WL',CKF%cprn(isat),ind,num_wl,OB%elev(isat)*RAD2DEG,rwl,pxl,AM(ind)%abwl,UPD%wfcb(isat),0.d0,'FIX_WL'
                        END IF
                    ELSE
                        IF (irectype .NE. 0) THEN
                            WRITE(1004,'((A),A4,2I4,F10.3,5F15.3)')'WL',CKF%cprn(isat),ind,0.d0,OB%elev(isat)*RAD2DEG,rwl,pxl,AM(ind)%abwl,UPD%wfcb(isat),UPD%recfcb(isat,2,irectype)
                        ELSE
                            WRITE(1004,'((A),A4,2I4,F10.3,5F15.3)')'WL',CKF%cprn(isat),ind,0.d0,OB%elev(isat)*RAD2DEG,rwl,pxl,AM(ind)%abwl,UPD%wfcb(isat),0.d0
                        END IF
                    END IF      
                END DO
            ELSE
                WRITE(OUTPUT_UNIT,'(A,2F10.3)') '-> BAD SD WL FIXING: ',(sd_num-ndl)*1.d0/sd_num, vxl
            END IF
        END DO
    END IF

    SL.fix_wl = 0
    IF (SL.fixnum_wl .GT. 4) THEN
        SL.fix_wl = 1
    END IF

    ! STEP 5. FIX NL(N1)
    IF (NM%ns .GT. 0) THEN
        ALLOCATE(AB(NM%ns))
    ELSE
        RETURN
    END IF
    DO i=1, NM%ns
        IF (CKF%lamb.EQ..FALSE. .OR. (CKF%lamb.EQ..TRUE. .AND. AM(i)%ifab.LT.4)) AM(i)%abnl = 0.5d0
        AB(i)%pab=i
    END DO

    !! mapping IF LC to N1, is same as ppp_map_invx
    ! map=0.d0
    ! mapt=0.d0
    ! DO i=1, QM%ntot
    !     map(i,i) = 1.d0
    !     mapt(i,i)= 1.d0
    ! END DO
    ! DO i=QM%nxyz+1, QM%ntot
    !     IF (TRIM(AM(i-QM%nxyz)%pname) .NE. 'AMBL1') CYCLE
    !     ind = AB(i-QM%nxyz)%pab
    !     isat= AM(ind)%psat
    !     map(i,i)  = 1.d0/SAT(isat)%lamdn
    !     mapt(i,i) = 1.d0/SAT(isat)%lamdn
    ! END DO
    ! DO i=1, QM%ntot
    !     DO j=i+1, QM%ntot
    !         QM%invx(i,j)=QM%invx(j,i)
    !     END DO
    ! END DO
    ! CALL matmpy(map,QM%invx(1:QM%ntot,1:QM%ntot),maptx,QM%ntot,QM%ntot,QM%ntot)
    ! CALL matmpy(maptx,mapt,QM%invx(1:QM%ntot,1:QM%ntot),QM%ntot,QM%ntot,QM%ntot)

    DO isys=1, CKF%nsys
        nfix(isys) = 0
        DO i=1, NM%ns
            isat = AM(i)%psat
            IF (CKF%cprn(isat)(1:1) .NE. CKF%system(isys:isys)) CYCLE
            IF (CKF%cprn(isat) .EQ. CKF%refcprn(isys)) THEN
                AB(i)%abst = -1
                CYCLE
            END IF

            IF (CKF%ArBiasMode .EQ. 'FCB') THEN
                IF (AM(i)%iobs.EQ.0 .OR. AM(i)%ifab.LT.2 .OR. UPD%nfcb(isat).EQ.10.d0 .OR. ABS(AM(i)%abwl-int(AM(i)%abwl)).NE.0.d0 .OR. AM(i)%elev/AM(i)%iobs.LE.CKF%cutoff) THEN
                    AB(i)%abst = -1
                    AB(i)%abfr = AM(i)%famb
                ELSE
                    AB(i)%abst = 0
                    AB(i)%abfr = AM(i)%famb/SAT(isat)%lamdn - AM(i)%abwl/(SAT(isat)%g-1.d0) - (UPD%nfcb(isat)-UPD%nfcb(CKF%refnprn(isys)))
                    IF (AM(i)%ifab .EQ. 4) THEN
                        nfix(isys) = nfix(isys) + 1
                        AB(i)%abst = 1
                        AB(i)%abfx = AM(i)%abnl
                    END IF
                END IF
            ELSE IF (CKF%ArBiasMode .EQ. 'OSB') THEN
                IF (AM(i)%iobs.EQ.0 .OR. AM(i)%ifab.LT.2 .OR. INDEX(CKF%nofixsat,CKF%cprn(isat)).NE.0 .OR. ABS(AM(i)%abwl-int(AM(i)%abwl)).NE.0.d0 .OR. AM(i)%elev/AM(i)%iobs.LE.CKF%cutoff) THEN
                    AB(i)%abst = -1
                    AB(i)%abfr = AM(i)%famb
                ELSE
                    AB(i)%abst = 0
                    AB(i)%abfr = AM(i)%famb/SAT(isat)%lamdn - AM(i)%abwl/(SAT(isat)%g-1.d0)
                    IF (AM(i)%ifab .EQ. 4) THEN
                        nfix(isys) = nfix(isys) + 1
                        AB(i)%abst = 1
                        AB(i)%abfx = AM(i)%abnl
                    END IF
                END IF                
            END IF
        END DO
    END DO

    SL%ncad = 0
    SL%nfix = 0
    SL%ratio = 0.d0
    IF (SUM(nfix) .GE. 0) THEN
        !! The following code is to get the value and variance for the
        !! unfixed ambiguities if there are some ambiguities are fixed
        !! so place fixed ambiguities at the end
        CALL ppp_plc_fixed(AB,QM,QM%invx)

        !! transform invx from LC to L1
        CALL ppp_map_invx(AM,SAT,AB,QM,QM%invx)

        !! impose constraints from already fixed ambiguities
        CALL ppp_add_ambcon(PM,AB,QM,QM%invx)

        !! place precise ambiguities at the end
        CALL ppp_plc_float(AB,QM,QM%invx)

        !! ambiguity resolution using LAMBDA method
        SL%ncad = QM%ndam + QM%nfix
        IF (QM%ndam .GT. 0) THEN

            CALL ppp_abfix_lambda(AB,QM,QM%invx,CKF%nl_maxdel,CKF%nl_minsav,CKF%nl_ratio,SL%ratio)

            !! store fixed widelane
            IF (QM%ncad .GT. 0) THEN
                SL%nfix = QM%ncad + QM%nfix

                !! sort inverted normal matrix
                CALL ppp_plc_fixed(AB,QM,QM%invx)

                !! further apply newly fixed widelane
                CALL ppp_add_ambcon(PM,AB,QM,QM%invx)
            END IF

            !! record ambiguity information
            DO i=1, NM%ns
                IF (AB(i)%abst.EQ.1 .OR. AB(i)%abst.EQ.2 .AND. QM%ncad.GT.0) THEN
                    SIT%nfixnl = SIT%nfixnl + 1
                    IF (CKF%lamb.EQ..TRUE. .AND. AM(AB(i)%pab)%ifab.EQ.4) CYCLE
                    AM(AB(i)%pab)%ifab = 4
                    AM(AB(i)%pab)%abnl = AB(i)%abfx
                ELSE
                    AM(AB(i)%pab)%abnl = AB(i)%abfr
                END IF
            END DO
        ELSE
            IF (CKF%lamb .EQ. .TRUE.) THEN
                SL%nfix = QM%nfix
                DO i=1, NM%ns
                    IF (AB(i)%abst.EQ.1 .OR. AB(i)%abst.EQ.2) THEN
                        SIT%nfixnl = SIT%nfixnl + 1
                    END IF
                END DO
            END IF
        END IF       
    END IF

    SL%fix_nl = 0
    SL%fixnum_nl = SIT%nfixnl
    IF (SIT%nfixnl .GT. 4) THEN
        SL%fix_nl = 1
    END IF

    WRITE(1004,'(A)')'---PRN--------------------FLT------------FIX------------UPD----------rBIAS-------------CMT'
    k=0
    DO i=1, NM%ns
        IF (TRIM(AM(AB(i)%pab)%pname).NE.'AMBL1') CYCLE
        k = k+1
        ind = AB(i)%pab
        isat= AM(ind)%psat
        IF (OB%omc(isat,MAXFREQ+1).EQ.0.d0) CYCLE
        isys = INDEX(CKF%system,CKF%cprn(isat)(1:1))
        IF (irectype .NE. 0) THEN
            IF (CKF%cprn(isat) .EQ. CKF%refcprn(isys)) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3,5X,(A))') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),UPD%recfcb(isat,1,irectype),'  -------> REF_SAT'
            END IF
            IF (AB(i)%abst.EQ.1 .OR. AB(i)%abst.EQ.2) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3,5X,(A))') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),UPD%recfcb(isat,1,irectype),'  FIX_NL'
            END IF
            IF (CKF%cprn(isat) .NE. CKF%refcprn(isys) .AND. AB(i)%abst .LE. 0) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3)') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),UPD%recfcb(isat,1,irectype)
            END IF
        ELSE
            IF (CKF%cprn(isat) .EQ. CKF%refcprn(isys)) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3,5X,(A))') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),0.d0,'  -------> REF_SAT'
            END IF
            IF (AB(i)%abst.EQ.1 .OR. AB(i)%abst.EQ.2) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3,5X,(A))') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),0.d0,'  FIX_NL'
            END IF
            IF (CKF%cprn(isat) .NE. CKF%refcprn(isys) .AND. AB(i)%abst .LE. 0) THEN
                WRITE(1004,'((A),A4,2I4,4F15.3)') 'NL',CKF%cprn(isat),k,AB(i)%abst,AB(i)%abfr,AB(i)%abfx,UPD%nfcb(isat),0.d0
            END IF            
        END IF
    END DO

    !! clean memory
    IF (NM.ns .GT. 0) DEALLOCATE(AB)

    RETURN

END SUBROUTINE
