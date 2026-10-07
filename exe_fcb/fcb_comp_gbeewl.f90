SUBROUTINE fcb_comp_gbeewl(CKF,SAT,UPD,NM,PM)
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! Start the exectuable code
!!-----------------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(INFM) :: NM(1:*)
TYPE(PRMT) :: PM(CKF.nsys+CKF.nprn+1+2+MAXFREQ*(CKF.nprn+CKF.nprn)+(MAXFREQ-1)*MAXSYS,1:*)
TYPE(FCB) :: UPD

  !*
  ! The local variables
  !!-------------------------
  LOGICAL(LG) :: lca(MAXPP),lfirst
  INTEGER(IT) :: i,j,k,l,m,isit,isat,jsat,tmp(2),itr
  INTEGER(IT) :: npp,ncd(MAXPP),ifg(MAXSIT*2),iss(2,MAXPP)
  REAL(RL) :: rwl,swl,rlc,alpha,sigm,resi,ptime(2,2)
  REAL(RL) :: wval(MAXSAT),sig(MAXSAT),rnl(MAXSIT*2),wgt(MAXSIT*2)
  REAL(RL) :: fnl(MAXPP),vnl(MAXPP),lval(MAXSAT),elev(MAXPP)
  INTEGER(IT) :: isys,nprn(MAXSYS),ipsat,jpsat
  CHARACTER(LEN_PRN) :: cprn(MAXSAT,MAXSYS)

  DATA lfirst /.TRUE./
  DATA wval /MAXSAT*10.d0/
  SAVE cprn,lfirst,nprn,wval

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    nprn=0
    cprn=''
    DO i=1, CKF.nprn
      j=INDEX(CKF.system,CKF.cprn(i)(1:1))
      IF (j .NE. 0) THEN
        nprn(j)=nprn(j)+1
        cprn(nprn(j),j)=CKF.cprn(i)
      END IF
    END DO
    lfirst=.FALSE.
  END IF

  DO isys=1, CKF.nsys

    IF (CKF.system(isys:isys) .EQ. 'R') CYCLE

    lval=10.d0
    DO i=1, nprn(isys)
      k=pointer_string(CKF.nprn,CKF.cprn,cprn(i,isys))
      IF (k .NE. 0) THEN
        lval(i)=wval(k)
      END IF
    END DO

    !!DO itr=1, 2

    lca=.FALSE.
    npp=0
    wgt=1.d0
    iss=0
    !! to find ambiguities for satellite pairs: ipsat, jpsat
    DO isat=1, nprn(isys)

      ipsat=pointer_string(CKF.nprn,CKF.cprn,cprn(isat,isys))

      DO jsat=isat+1,nprn(isys)

        jpsat=pointer_string(CKF.nprn,CKF.cprn,cprn(jsat,isys))

        k=0
        !! for each site
        DO isit=1, CKF.nsit
          !! xsy: Bad station
          IF (NM(isit).esig .GT. CKF.DiaSig) CYCLE

          DO i=NM(isit).np+1,NM(isit).imtx

            !! the frequency 2
            IF (PM(i,isit).pcode(3) .NE. 2) CYCLE
            !! high elevation
            IF (PM(i,isit).elev/PM(i,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
            !! ipsat or jpsat 
            IF (PM(i,isit).pcode(2).NE.ipsat .AND. PM(i,isit).pcode(2).NE.jpsat) CYCLE

            !! the second frequency
            DO l=NM(isit).np+1,NM(isit).imtx
              !! frequency 4
              IF (PM(l,isit).pcode(3) .NE. 4) CYCLE
              !! high elevation
              IF (PM(l,isit).elev/PM(l,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
              !! same satellite
              IF (PM(l,isit).pcode(2) .NE. PM(i,isit).pcode(2)) CYCLE
              !! same time period
              IF (PM(l,isit).ptime(2) .LT. PM(i,isit).ptime(1)) CYCLE
              IF (PM(l,isit).ptime(1) .GT. PM(i,isit).ptime(2)) CYCLE
              ptime(1,1)=MAX(PM(i,isit).ptime(1),PM(l,isit).ptime(1))
              !ptime(1,2)=MIN(PM(i,isit).ptime(2),PM(l,isit).ptime(2))
              ptime(1,2)=CKF.mjd+CKF.sod/86400.d0
              !! common time 
              IF ((ptime(1,2)-ptime(1,1))*86400.d0 .LE. CKF.minsec_common) CYCLE
              EXIT
            END DO

            ! find the frequency 2/4 ambiguity for ipsat/jpsate at isite

            IF (l.GT.NM(isit).imtx .OR. i.GT.NM(isit).imtx) CYCLE

            DO j=i+1, NM(isit).imtx

              IF (PM(j,isit).pcode(3) .NE. 2) CYCLE
              IF (PM(j,isit).elev/PM(j,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
              ! shoud be ipsat/jpsat
              IF (PM(j,isit).pcode(2).NE.ipsat.AND.PM(j,isit).pcode(2).NE.jpsat) CYCLE
              ! not the same satellite
              IF (PM(j,isit).pcode(2) .EQ. PM(i,isit).pcode(2)) CYCLE

              !! the second frequency
              DO m=NM(isit).np+1, NM(isit).imtx
                IF (PM(m,isit).pcode(3) .NE. 4) CYCLE
                IF (PM(m,isit).elev/PM(j,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
                ! same satellite
                IF (PM(m,isit).pcode(2) .NE. PM(j,isit).pcode(2)) CYCLE
                ! common perid
                IF (PM(m,isit).ptime(2) .LT. PM(j,isit).ptime(1)) CYCLE
                IF (PM(m,isit).ptime(1) .GT. PM(j,isit).ptime(2)) CYCLE
                ptime(2,1)=MAX(PM(j,isit).ptime(1),PM(m,isit).ptime(1))
                !ptime(2,2)=MIN(PM(j,isit).ptime(2),PM(m,isit).ptime(2))
                ptime(2,2)=CKF.mjd+CKF.sod/86400.d0
                IF ((ptime(2,2)-ptime(2,1))*86400.d0 .LE. CKF.minsec_common) CYCLE
                EXIT
              END DO

              IF (m.GT.NM(isit).imtx .OR. j.GT.NM(isit).imtx) CYCLE

              ! satellite differences
              ! the common time
              IF ((MIN(ptime(1,2),ptime(2,2))-MAX(ptime(1,1),ptime(2,1)))*86400.d0 .LE. CKF.minsec_common) CYCLE

              IF (PM(i,isit).pcode(2).EQ.ipsat .AND. PM(j,isit).pcode(2).EQ.jpsat) THEN
                ! N2-N4, isat-jsat
                rwl=(PM(i,isit).xest/SAT(ipsat).lamda(2)-PM(l,isit).xest/SAT(ipsat).lamda(4))-&
                    (PM(j,isit).xest/SAT(jpsat).lamda(2)-PM(m,isit).xest/SAT(jpsat).lamda(4))
                swl=DSQRT((PM(i,isit).xsig/SAT(ipsat).lamda(2))**2+(PM(l,isit).xsig/SAT(ipsat).lamda(4))**2+&
                    (PM(j,isit).xsig/SAT(jpsat).lamda(2))**2+(PM(m,isit).xsig/SAT(jpsat).lamda(4))**2)
              ELSE IF (PM(i,isit).pcode(2).EQ.jpsat .AND. PM(j,isit).pcode(2).EQ.ipsat) THEN
                rwl=(PM(j,isit).xest/SAT(ipsat).lamda(2)-PM(m,isit).xest/SAT(ipsat).lamda(4))-&
                    (PM(i,isit).xest/SAT(jpsat).lamda(2)-PM(l,isit).xest/SAT(jpsat).lamda(4))
                swl=DSQRT((PM(j,isit).xsig/SAT(ipsat).lamda(2))**2+(PM(m,isit).xsig/SAT(ipsat).lamda(4))**2+&
                    (PM(i,isit).xsig/SAT(jpsat).lamda(2))**2+(PM(l,isit).xsig/SAT(jpsat).lamda(4))**2)
              END IF

              k=k+1
              IF (k .GT. MAXSIT*2) THEN
                WRITE(ERROR_UNIT,'(A)') '***ERROR(fcb_comp_wide): too many ambiguities'
                CALL exit(1)
              END IF
              rnl(k)=rwl
              wgt(k)=1.d0
              ifg(k)=0
              ! NINT(k)+FCB=FLOAT
              rnl(k)=rnl(k)-NINT(rnl(k))
              elev(k)=(PM(i,isit).elev+PM(j,isit).elev)/2.d0

            END DO
          END DO

        END DO
        !! five stations at least
        ! IF (k .LT. 5) CYCLE
        IF (k .LT. CKF%commonsit) CYCLE

        ! rearrange with elevation
        DO i=1, k-1
          DO j=1, k
            IF (elev(i) .LT. elev(j)) THEN
              resi=elev(i)
              elev(i)=elev(j)
              elev(j)=resi

              resi=rnl(i)
              rnl(i)=rnl(j)
              rnl(j)=resi
            END IF
          END DO
        END DO

        !! average collected fractional parts
        npp=npp+1
        iss(1,npp)=isat
        iss(2,npp)=jsat
        lca(isat) =.TRUE.
        lca(jsat) =.TRUE.
        j=0
        fnl(npp)=10.d0
        CALL proc_xl_dirc(k,rnl,wgt,ifg,j,fnl(npp),vnl(npp),sigm)
        ncd(npp)=k-j
        !WRITE(OUTPUT_UNIT,'(2(A3,1X),f8.3,i4,i3,2f8.3)') CKF.cprn(ipsat),CKF.cprn(jpsat),fnl(npp),k,j,vnl(npp),sigm
        ! IF (ncd(npp) .LT. CKF%commonsit) npp=npp-1
      END DO
    END DO
    !IF (npp .EQ. 0) RETURN
    IF (npp .EQ. 0) CYCLE

    !! reset spval when necessary
    DO isat=1, nprn(isys)
      i=pointer_string(CKF.nprn,CKF.cprn,cprn(isat,isys))
      IF (lca(isat) .EQ. .FALSE.) THEN
        lval(isat)=10.d0
      END IF
    END DO

    !! sort fnl
    DO i=1,npp
      lca(i)=.FALSE.
      DO j=i+1,npp
        IF (ncd(i) .LT. ncd(j)) THEN
          k  =ncd(i); ncd(i)=ncd(j); ncd(j)=k
          sigm=fnl(i); fnl(i)=fnl(j); fnl(j)=sigm
          sigm=vnl(i); vnl(i)=vnl(j); vnl(j)=sigm
          tmp(1:2)=iss(1:2,i);iss(1:2,i)=iss(1:2,j);iss(1:2,j)=tmp(1:2)
        END IF
      END DO
    END DO

    !! coarse alignment
    DO WHILE(.TRUE.)
      k=0
      DO i=1,npp
        IF (lca(i)) CYCLE
        isat=iss(1,i)
        jsat=iss(2,i)
        IF (lval(isat).NE.10.d0 .AND. lval(jsat).NE.10.d0) THEN
          k=k+1
          lca(i)=.TRUE.
          fnl(i)=fnl(i)-NINT(fnl(i)-lval(isat)+lval(jsat))
        ELSE IF (lval(isat).NE.10.d0 .AND. lval(jsat).EQ.10.d0) THEN
          k=k+1
          lca(i)=.TRUE.
          lval(jsat)=lval(isat)-fnl(i)
          fnl(i)    =fnl(i)+NINT(lval(jsat))
          lval(jsat)=lval(jsat)-NINT(lval(jsat))
        ELSE IF (lval(isat).EQ.10.d0 .AND. lval(jsat).NE.10.d0) THEN
          k=k+1
          lca(i)=.TRUE.
          lval(isat)=lval(jsat)+fnl(i)
          fnl(i)     =fnl(i)-NINT(lval(isat))
          lval(isat)=lval(isat)-NINT(lval(isat))
        END IF
      END DO
      IF (k .EQ. 0) THEN
        IF (ALL(lca(1:npp) .EQ. .TRUE.)) THEN
          EXIT
        ELSE IF (ALL(lval(1:nprn(isys)) .EQ. 10.d0)) THEN
          lca(1)=.TRUE.
          lval(iss(1,1))=0.d0
          lval(iss(2,1))=-fnl(1)
        ELSE
          WRITE(OUTPUT_UNIT,'(A)') '###WARNING(fcb_comp_naro): over 1 reference '
          j=1
          DO WHILE(j .LE. npp)
            IF (lca(j) .EQ. .FALSE.) THEN
              WRITE(OUTPUT_UNIT,'(A,2(1X,A3))') ' Pair Removed: ',cprn(iss(1,j),isys),cprn(iss(2,j),isys)
              k=j
              DO WHILE(k.LE.npp-1)
                ncd(k)    =ncd(k+1)
                fnl(k)    =fnl(k+1)
                vnl(k)    =vnl(k+1)
                iss(1:2,k)=iss(1:2,k+1)
                lca(k)=lca(k+1)
                k=k+1
              END DO
              npp=npp-1
              CYCLE
            END IF
            j=j+1
          END DO
          EXIT
        END IF
      END IF
    END DO

    !! 2023-09-11 by xsy: only to constrans the reference satellite
    ! isat = pointer_string(CKF%nprn,CKF%cprn,cprn(iss(1,1),isys))
    ! IF (UPD%eewfcb(isat) .NE. 10.d0) THEN
    !   lval(iss(1,1))=UPD%eewfcb(isat)
    ! END IF    

    !! 2023-09-11 by xsy: the before epoch's FCB will be this epoch's initial value
    ! DO i=1,npp
    !   isat = pointer_string(CKF%nprn,CKF%cprn,cprn(iss(1,i),isys))
    !   IF (UPD%eewfcb(isat) .NE. 10.d0) THEN
    !     lval(iss(1,i))=UPD%eewfcb(isat)
    !   END IF
    !   isat = pointer_string(CKF%nprn,CKF%cprn,cprn(iss(2,i),isys))
    !   IF (UPD%eewfcb(isat) .NE. 10.d0) THEN
    !     lval(iss(2,i))=UPD%eewfcb(isat)
    !   END IF
    ! END DO

    CALL alinew(npp,iss,fnl,ncd,vnl,nprn(isys),cprn(1,isys),lval,sig)
    UPD.updrefsat(isys,4) = cprn(iss(1,1),isys)

    !!END DO

    DO i=1,nprn(isys)
      j=pointer_string(CKF.nprn,CKF.cprn,cprn(i,isys))
      wval(j)=lval(i)
      UPD.eewfcb(j)=lval(i)
      UPD.eewsl(j)=sig(i)
    END DO

    !! output residuals
    !WRITE(OUTPUT_UNIT,'(A)') 'EWL Residuals after precise alignment:'
    WRITE(4000,'((A),I7,F10.2,A3,I7)') 'EEWL Residuals ->',CKF.mjd,CKF.sod,CKF.system(isys:isys),npp
    DO i=1,npp
      isat=pointer_string(CKF.nprn,CKF.cprn,cprn(iss(1,i),isys))
      jsat=pointer_string(CKF.nprn,CKF.cprn,cprn(iss(2,i),isys))
      resi=wval(isat)-wval(jsat)-fnl(i)
      WRITE(4000,'(2(A3,1x),f8.3,i4,f8.3)') cprn(iss(1,i),isys),cprn(iss(2,i),isys),resi,ncd(i),vnl(i)
    END DO

  END DO

  RETURN

END SUBROUTINE
