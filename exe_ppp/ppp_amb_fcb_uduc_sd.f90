!*
SUBROUTINE ppp_amb_fcb_uduc_sd(CKF,OB,QM,AM,NM,UPD)
! STEP1: MAP UNDIFFERENCE AMBIGUITY TO SATELLITE AMBIGUITY
!     Dy = K*Dx*kT
!     Dx = Qxx=(ATA)-1
!     K  = map
! STEP2: MAP L1,L2... TO WL.NL...
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE observation
IMPLICIT NONE

!*
! The argument
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(INVM) :: QM
TYPE(AMBT) :: AM(1:*)
TYPE(INFM) :: NM
TYPE(RNXOBS) :: OB
TYPE(FCB) :: UPD

    !*
    ! The local variables
    !!------------------------
    REAL(RL),ALLOCATABLE :: map(:,:),mapt(:,:),maptx(:,:)
    REAL(RL) :: maxele,a,b,c
    INTEGER(IT) :: i,j,isys,isat,ifreq,ifg(MAXSAT),ipt(MAXSAT),maxtim,ind
    INTEGER(IT) :: kpt(MAXSYS),nxl(MAXSYS)                  !kpt: reference satellite   nxl: satellite number

    !*
    ! Start the exectuable code
    !!------------------------
    
    ALLOCATE(map(QM%ntot,QM%ntot))
    ALLOCATE(mapt(QM%ntot,QM%ntot))
    ALLOCATE(maptx(QM%ntot,QM%ntot))

    !! @ CMT BY XSY: ============================================================-> STEP 1 <-============================================================
    ! map undiffereance ambiguity to satellite-difference ambiguity
    map=0.d0
    mapt=0.d0
    DO i=1, QM%ntot
        map(i,i)=1.d0
        mapt(i,i)=1.d0
    END DO

    ! choose reference satellite based on elevation for ammbfix
    DO isys=1,CKF%nsys
        nxl(isys)=0
        kpt(isys)=0
        ipt=0
        DO i=1,NM%ns
            isat=AM(i)%psat
            ifreq=AM(i)%ifreq
            IF (CKF%cprn(isat)(1:1).NE.CKF%system(isys:isys)) CYCLE
            IF (TRIM(AM(i)%pname).NE.'AMBL1') CYCLE
            IF (OB%omc(isat,ifreq).EQ.0.D0) CYCLE
            
            !@ CMT BY XSY: CHECK THE UPD WHETHER IS EXIST
            ind=INDEX(SYS,CKF%cprn(isat)(1:1))
            IF (CKF%nfreq(ind) .EQ. 2) THEN
                IF(UPD%nfcb(isat).EQ.10.d0 .OR. UPD%wfcb(isat).EQ.10.d0) CYCLE
            ELSE IF (CKF%nfreq(ind) .EQ. 3) THEN
                IF(UPD%nfcb(isat).EQ.10.d0 .OR. UPD%wfcb(isat).EQ.10.d0 .OR. UPD%ewfcb(isat).EQ.10.d0) CYCLE                
            ELSE IF (CKF%nfreq(ind) .EQ. 4) THEN
                IF(UPD%nfcb(isat).EQ.10.d0 .OR. UPD%wfcb(isat).EQ.10.d0 .OR. UPD%ewfcb(isat).EQ.10.d0 .OR. UPD%eewfcb(isat).EQ.10.d0) CYCLE
            ELSE IF (CKF%nfreq(ind) .EQ. 5) THEN
                IF(UPD%nfcb(isat).EQ.10.d0 .OR. UPD%wfcb(isat).EQ.10.d0 .OR. UPD%ewfcb(isat).EQ.10.d0 .OR. UPD%eewfcb(isat).EQ.10.d0 .OR. UPD%hewfcb(isat).EQ.10.d0) CYCLE
            END IF

            !@ CMT BY XSY: FOR GPS MULTI-FREQUENCY, THE DUAL-FREQUENCY SATELLITES IS NOT CHOOSED AS REFERENCE
            IF (CKF%cprn(isat)(1:1).EQ.'G' .AND. CKF%nfreq(INDEX(SYS,'G')).EQ.3) THEN
                IF (OB%omc(isat,3).EQ.0.D0) CYCLE
            END IF

            IF(AM(i)%iobs .EQ.0 .OR. AM(i)%elev/AM(i)%iobs .LE.CKF%cutoff) CYCLE
            nxl(isys)=nxl(isys)+1
            ifg(nxl(isys))=0
            ipt(nxl(isys))=i
        END DO

        maxtim=0
        maxele=0.d0
        
        IF(nxl(isys).GT.0)THEN
            ! DO i=1,nxl(isys)
            !     IF(ifg(i).NE.0) CYCLE
            !     IF(AM(ipt(i))%iobs .GT.maxtim)THEN
            !         kpt(isys)=ipt(i)
            !         maxtim=AM(ipt(i))%iobs
            !         maxele=AM(ipt(i))%elev/AM(ipt(i))%iobs
            !     ELSE IF(AM(ipt(i))%iobs .EQ.maxtim .AND. AM(ipt(i))%elev/AM(ipt(i))%iobs .GT.maxele)THEN
            !         kpt(isys)=ipt(i)
            !         maxele=AM(ipt(i))%elev/AM(ipt(i))%iobs
            !     END IF
            ! END DO
            DO i=1, nxl(isys)
                IF(AM(ipt(i))%elev/AM(ipt(i))%iobs .GT. maxele)THEN
                    kpt(isys) = ipt(i)
                    maxele = AM(ipt(i))%elev/AM(ipt(i))%iobs
                END IF
            END DO                  
            IF(kpt(isys).NE.0)THEN
                isat=AM(kpt(isys))%psat
                CKF%refnprn(isys)=isat
                CKF%refcprn(isys)=CKF%cprn(isat)
                ! WRITE(*,'((A),2I4,A6,I4,3X,(A),F5.2)')'REFSAT: ',isys,CKF%refnprn(isys),CKF%refcprn(isys),kpt(isys),AM(kpt(isys))%pname,maxele
            END IF
        ELSE
            CYCLE
        END IF

        !@ CMT BY XSY: SYS INDEX
        DO i=QM%nxyz+1,QM%ntot
            isat=AM(i-QM%nxyz)%psat
            ifreq=AM(i-QM%nxyz)%ifreq

            IF (CKF%cprn(isat)(1:1).NE.CKF%system(isys:isys)) CYCLE
            IF (OB%omc(isat,ifreq).EQ.0.D0) CYCLE

            IF (isat.EQ.CKF%refnprn(isys)) THEN
                map(i,i)=0.d0
                mapt(i,i)=0.d0
            ELSE
                map(i,i)=1.d0
                map(i,kpt(isys)+QM%nxyz+ifreq-1)=-1.d0
                mapt(i,i)=1.d0
                mapt(kpt(isys)+QM%nxyz+ifreq-1,i)=-1.d0
                !@ CMT BY XSY: SATELLITE-DIFFERENCE FLOAT AMBIGUITY
                AM(i-QM%nxyz)%famb=AM(i-QM%nxyz)%famb-AM(kpt(isys)+ifreq-1)%famb
            END IF
        END DO
        !@ CMT BY XSY: CLEAN THE REFERENCE AMBIGUITY
        DO i=QM%nxyz+1,QM%ntot
            isat=AM(i-QM%nxyz)%psat
            ifreq=AM(i-QM%nxyz)%ifreq
            IF (CKF%cprn(isat)(1:1).NE.CKF%system(isys:isys)) CYCLE
            IF (OB%omc(isat,ifreq).EQ.0.D0) CYCLE

            IF (isat.EQ.CKF%refnprn(isys)) THEN
                AM(i-QM%nxyz)%famb=0.d0
            END IF
        END DO      
    END DO
    
    DO i=1, QM%ntot
        DO j=i+1, QM%ntot
            QM%invx(i,j)=QM%invx(j,i)
        END DO
        ! WRITE(999,'(<QM.ntot>(F14.4,2X))')(map(i,j),j=1,QM.ntot)
    END DO

    CALL matmpy(map,QM%invx(1:QM%ntot,1:QM%ntot),maptx,QM%ntot,QM%ntot,QM%ntot)
    CALL matmpy(maptx,mapt,QM%invx(1:QM%ntot,1:QM%ntot),QM%ntot,QM%ntot,QM%ntot)
    
    !! @ CMT BY XSY: ============================================================-> STEP 2 <-============================================================
    ! map raw ambiguity to liner combination ambiguity
    map=0.d0
    mapt=0.d0
    DO i=1, QM%ntot
        map(i,i)=1.d0
        mapt(i,i)=1.d0
    END DO

    ! map QXX ambigtuiry for N1,N2,N5 to N1,Nwl,Newl
    DO i=QM%nxyz+1,QM%ntot
        isat=AM(i-QM%nxyz)%psat
        ifreq=AM(i-QM%nxyz)%ifreq
        isys=INDEX(SYS,CKF%cprn(isat)(1:1))    
        IF (TRIM(AM(i-QM%nxyz)%pname) .NE. 'AMBL1') CYCLE

        ! L1
        map(i,i)=1.d0
        mapt(i,i)=1.d0
        DO j=1,CKF%nfreq(isys)
            IF (OB%obs(isat,j).EQ.0.D0) CYCLE
            ! L1-L2
            IF (j.EQ.2) THEN
                map(i+1,i)=1.d0
                map(i+1,i+1)=-1.d0

                mapt(i,i+1)=1.d0
                mapt(i+1,i+1)=-1.d0            
            ! L2-Ln
            ELSE IF (j.GT.2) THEN
                map(i+j-1,i+1)= 1.d0
                map(i+j-1,i+j-1)=-1.d0 

                mapt(i+1,i+j-1)= 1.d0
                mapt(i+j-1,i+j-1)=-1.d0             
            END IF
        END DO
    END DO

    DO i=1, QM%ntot
        DO j=i+1, QM%ntot
            QM%invx(i,j)=QM%invx(j,i)
        END DO
    END DO

    ! mapping raw ambiguity to L1, WL, EWL combination
    CALL matmpy(map,QM%invx(1:QM%ntot,1:QM%ntot),maptx,QM%ntot,QM%ntot,QM%ntot)
    CALL matmpy(maptx,mapt,QM%invx(1:QM%ntot,1:QM%ntot),QM%ntot,QM%ntot,QM%ntot)

    DEALLOCATE(map)
    DEALLOCATE(mapt)
    DEALLOCATE(maptx)

    RETURN
  
END SUBROUTINE