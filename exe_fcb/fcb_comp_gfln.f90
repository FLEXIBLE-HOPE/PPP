! Compute narrow-lane FCBs
! Jianghui Geng
! March 1 2012

SUBROUTINE fcb_comp_gfln(CKF,SAT,UPD,NM,PM)
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
  REAL(RL) :: nval(MAXSAT),sig(MAXSAT),rnl(MAXSIT*2),wgt(MAXSIT*2)
  REAL(RL) :: fnl(MAXPP),vnl(MAXPP),lval(MAXSAT),elev(MAXPP)
  INTEGER(IT) :: isys,nprn(MAXSYS),ipsat,jpsat
  CHARACTER(LEN_PRN) :: cprn(MAXSAT,MAXSYS)

  DATA lfirst /.TRUE./
  DATA nval /MAXSAT*10.d0/
  SAVE cprn,lfirst,nprn,nval

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
        lval(i)=nval(k)
      END IF
    END DO

    !!DO itr=1, 2

    lca=.FALSE.
    npp=0
    wgt=1.d0
    iss=0
    DO isat=1, nprn(isys)

      ipsat=pointer_string(CKF.nprn,CKF.cprn,cprn(isat,isys))
      IF (UPD.wfcb(ipsat).EQ.10.d0) CYCLE

      DO jsat=isat+1,nprn(isys)

        jpsat=pointer_string(CKF.nprn,CKF.cprn,cprn(jsat,isys))
        IF (UPD.wfcb(jpsat).EQ.10.d0) CYCLE

        k=0
        !! for each site
        DO isit=1, CKF.nsit

          DO i=NM(isit).np+1,NM(isit).imtx

            !! the frist frequency
            IF (INDEX(PM(i,isit).pname,'AMBL1') .EQ. 0) CYCLE
            IF (PM(i,isit).elev/PM(i,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
            IF (PM(i,isit).pcode(2).NE.ipsat .AND. PM(i,isit).pcode(2).NE.jpsat) CYCLE

            !! the second frequency
            DO l=i+1,NM(isit).imtx
              IF (INDEX(PM(l,isit).pname,'AMBL2') .EQ. 0) CYCLE
              IF (PM(l,isit).elev/PM(l,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
              IF (PM(l,isit).pcode(2) .NE. PM(i,isit).pcode(2)) CYCLE
              IF (PM(l,isit).ptime(2) .LT. PM(i,isit).ptime(1)) CYCLE
              IF (PM(l,isit).ptime(1) .GT. PM(i,isit).ptime(2)) CYCLE
              ptime(1,1)=MAX(PM(i,isit).ptime(1),PM(l,isit).ptime(1))
              ptime(1,2)=MIN(PM(i,isit).ptime(2),PM(l,isit).ptime(2))
              IF ((ptime(1,2)-ptime(1,1))*86400.d0 .LE. CKF.minsec_common) CYCLE
              EXIT
            END DO

            IF (l .GT. NM(isit).imtx) THEN
              CYCLE
            END IF

            DO j=i+1, NM(isit).imtx

              IF (INDEX(PM(j,isit).pname,'AMBL1') .EQ. 0) CYCLE
              IF (PM(j,isit).elev/PM(j,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
              !! form single-difference wide-lane ambiguity
              IF (PM(j,isit).pcode(2).NE.ipsat.AND.PM(j,isit).pcode(2).NE.jpsat) CYCLE
              IF (PM(j,isit).pcode(2).EQ.PM(i,isit).pcode(2)) CYCLE

              !! the second frequency
              DO m=j+1, NM(isit).imtx
                IF (INDEX(PM(m,isit).pname,'AMBL2') .EQ. 0) CYCLE
                IF (PM(m,isit).elev/PM(j,isit).iobs*RAD2DEG .LT. CKF.cutoff) CYCLE
                IF (PM(m,isit).pcode(2) .NE. PM(j,isit).pcode(2)) CYCLE
                IF (PM(m,isit).ptime(2) .LT. PM(j,isit).ptime(1)) CYCLE
                IF (PM(m,isit).ptime(1) .GT. PM(j,isit).ptime(2)) CYCLE
                ptime(2,1)=MAX(PM(j,isit).ptime(1),PM(m,isit).ptime(1))
                ptime(2,2)=MIN(PM(j,isit).ptime(2),PM(m,isit).ptime(2))
                IF ((ptime(2,2)-ptime(2,1))*86400.d0 .LE. CKF.minsec_common) CYCLE
                EXIT
              END DO

              IF (m .GT. NM(isit).imtx) THEN
                CYCLE
              END IF

              IF ((MIN(ptime(1,2),ptime(2,2))-MAX(ptime(1,1),ptime(2,1)))*86400.d0 .LE. CKF.minsec_common) CYCLE

              IF (PM(i,isit).pcode(2).EQ.ipsat .AND. PM(j,isit).pcode(2).EQ.jpsat) THEN

                rwl=PM(l,isit).abwl-PM(m,isit).abwl-UPD.wfcb(ipsat)+UPD.wfcb(jpsat)

                !! LC ambiguities
                !rlc=(PM(i,isit).xest/SAT(ipsat).lamda(1)*SAT(ipsat).freq(1)/(SAT(ipsat).freq(1)-SAT(ipsat).freq(2))-&
                !     PM(l,isit).xest/SAT(ipsat).lamda(2)*SAT(ipsat).freq(2)/(SAT(ipsat).freq(1)-SAT(ipsat).freq(2)))-&
                !    (PM(j,isit).xest/SAT(jpsat).lamda(1)*SAT(jpsat).freq(1)/(SAT(jpsat).freq(1)-SAT(jpsat).freq(2))-&
                !     PM(m,isit).xest/SAT(jpsat).lamda(2)*SAT(jpsat).freq(2)/(SAT(jpsat).freq(1)-SAT(jpsat).freq(2)))


                rlc=(PM(i,isit).xest/SAT(ipsat).lamda(1)*SAT(ipsat).g**2/(SAT(ipsat).g**2-1)- &
                     PM(l,isit).xest/SAT(ipsat).lamda(2)*SAT(ipsat).g/(SAT(ipsat).g**2-1))-&
                    (PM(j,isit).xest/SAT(jpsat).lamda(1)*SAT(jpsat).g**2/(SAT(jpsat).g**2-1)- &
                     PM(m,isit).xest/SAT(jpsat).lamda(2)*SAT(jpsat).g/(SAT(jpsat).g**2-1))

              ELSE IF (PM(i,isit).pcode(2).EQ.jpsat .AND. PM(j,isit).pcode(2).EQ.ipsat) THEN
 
                rwl=PM(m,isit).abwl-PM(l,isit).abwl-UPD.wfcb(ipsat)+UPD.wfcb(jpsat)

                !! LC ambiguities
                !rlc=(PM(j,isit).xest/SAT(ipsat).lamda(1)*SAT(ipsat).freq(1)/(SAT(ipsat).freq(1)-SAT(ipsat).freq(2))-&
                !     PM(m,isit).xest/SAT(ipsat).lamda(2)*SAT(ipsat).freq(2)/(SAT(ipsat).freq(1)-SAT(ipsat).freq(2)))-&
                !    (PM(i,isit).xest/SAT(jpsat).lamda(1)*SAT(jpsat).freq(1)/(SAT(jpsat).freq(1)-SAT(jpsat).freq(2))-&
                !     PM(l,isit).xest/SAT(jpsat).lamda(2)*SAT(jpsat).freq(2)/(SAT(jpsat).freq(1)-SAT(jpsat).freq(2)))


                rlc=(PM(j,isit).xest/SAT(ipsat).lamda(1)*SAT(ipsat).g**2/(SAT(ipsat).g**2-1)- &
                     PM(m,isit).xest/SAT(ipsat).lamda(2)*SAT(ipsat).g/(SAT(ipsat).g**2-1))-&
                    (PM(i,isit).xest/SAT(jpsat).lamda(1)*SAT(jpsat).g**2/(SAT(jpsat).g**2-1)- &
                     PM(l,isit).xest/SAT(jpsat).lamda(2)*SAT(jpsat).g/(SAT(jpsat).g**2-1))

              END IF

              !! The sigma of widthline is essential for ambiguity resolution
              !swl=DSQRT(PM(l,isit).sigw**2+PM(m,isit).sigw**2+UPD.wsl(ipsat)**2+UPD.wsl(jpsat)**2)
              swl=DSQRT(UPD.wsl(ipsat)**2+UPD.wsl(jpsat)**2)

              !! attempt to fix this wide-lane ambiguity
              CALL prob_resol(rwl,swl,1,CKF.wl_maxdev,CKF.wl_maxsig,alpha)
              IF (alpha .GT. CKF.wl_alpha) THEN
                k=k+1          ! found one
                IF (k .GT. MAXSIT*2) THEN
                  WRITE(ERROR_UNIT,'(A)') '***ERROR(fcb_comp_naro): too many ambiguities'
                  CALL exit(1)
                END IF
                !rnl(k)=rlc-1.d0/(SAT(ipsat).g-1.D0)*NINT(rwl)
                rnl(k)=rlc*(SAT(ipsat).g+1)/SAT(ipsat).g-1/(SAT(ipsat).g-1)*NINT(rwl)
                wgt(k)=1.d0
                ifg(k)=0
                rnl(k)=rnl(k)-NINT(rnl(k))
                elev(k)=(PM(i,isit).elev+PM(j,isit).elev)/2.d0
              END IF

            END DO
          END DO

        END DO
        !! five stations at least
        ! IF (k .LT. 5) CYCLE
        IF (k .LT. CKF%commonsit) CYCLE

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
      END DO
    END DO
    !IF (npp .EQ. 0) RETURN
    IF (npp .EQ. 0) CYCLE

    !! reset spval when necessary
    DO isat=1, nprn(isys)
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

    CALL alinn(npp,iss,fnl,ncd,vnl,nprn(isys),cprn(1,isys),lval,sig)
    

    !!END DO

    DO i=1,nprn(isys)
      j=pointer_string(CKF.nprn,CKF.cprn,cprn(i,isys))
      nval(j)=lval(i)
      UPD.nfcb(j)=lval(i)
      UPD.nsl(j)=sig(i)
    END DO

    !! output residuals
    !WRITE(OUTPUT_UNIT,'(A)') 'NL Residuals after precise alignment:'
    DO i=1,npp
      isat=pointer_string(CKF.nprn,CKF.cprn,cprn(iss(1,i),isys))
      jsat=pointer_string(CKF.nprn,CKF.cprn,cprn(iss(2,i),isys))
      resi=nval(isat)-nval(jsat)-fnl(i)
      !!WRITE(OUTPUT_UNIT,'(2(A3,1x),f8.3,i4,f8.3)') cprn(isat,isys),cprn(jsat,isys),resi,ncd(i),vnl(i)
      WRITE(30000,'(2(A3,1x),f8.3,i4,f8.3)') cprn(isat,isys),cprn(jsat,isys),resi,ncd(i),vnl(i)
    END DO

  END DO

  RETURN

END SUBROUTINE
