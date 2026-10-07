!
!! purpose  : detection for biases and blunders
!! parameter:
!!    input : OB     -- observation
!!            NM,AM  -- information matrix & ambiguity parameters
!!            infs   -- one-dimensional array for information matrix
!!    output:
!! author   : Geng J
!! created  : Nov. 24, 2007
!!
!*
SUBROUTINE ppp_dat_screen(CKF,SAT,OB,OBH,AM,NM,infs)
!!
!*
USE info
USE ckdctrl
USE ambiguity
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------
TYPE(SATE) :: SAT(MAXSAT)
TYPE(RNXOBS) :: OB, OBH
TYPE(INFM) :: NM
TYPE(AMBT) :: AM(1:*)
TYPE(CKDCFG) :: CKF
REAL(RL) :: infs(1:*)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j,isat,iamb,nb,ifreq,Ins,InsSat,InsFreq
  INTEGER(IT) :: iset(MAXFREQ*MAXSAT),itg(MAXFREQ*MAXSAT),nbias,ibias(MAXFREQ*MAXSAT)
  REAL(RL) :: esig,slip(MAXFREQ*MAXSAT),xcor(MAXFREQ*MAXSAT),xest(MAXFREQ*MAXSAT),dump
  CHARACTER(LEN=5) :: camb

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!---------------------------

  !! whether cycle slips exist
  nb=-1   ! flag for no slips
  iset=0    ! corresponds to each carrier-phase obs
  esig=NM.esig
  !IF (NM.esig .GT. 4.0d0) THEN
  IF (NM.esig .GT. CKF.DiaSig) THEN

    !! step 1: compute sensitivity vectors
    j=0
    DO i=1,NM.nobs
      ! only phase
      IF (NM.ipob(i,2) .EQ. 1) THEN
        j=j+1
        itg(j)=i
        CALL ppp_scr_senvec(i,NM.weig(i),infs(NM.iptx(NM.imtx+1+j)+1),NM,infs)
      END IF
    END DO

    !! step 2: decide contaminated carrier phase
    IF (CKF.DiaSig .LT. 0.3d0) THEN
      CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
    ELSE
      CALL ppp_scr_identify(0.3d0,j,iset,NM,infs)
    END IF
    nb=COUNT(iset(1:NM.nobs/2) .NE. 0)

    IF (nb .EQ. 0) THEN
      IF (CKF.DiaSig .GT. 0.3d0) THEN
        IF (CKF.DiaSig .LT. 0.5d0) THEN
          CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
        ELSE
          CALL ppp_scr_identify(0.5d0,j,iset,NM,infs)
        END IF
        nb=COUNT(iset(1:NM.nobs/2) .NE. 0)

        IF (nb .EQ. 0) THEN
          IF (CKF.DiaSig .GT. 0.5d0) THEN
            IF (CKF.DiaSig .LT. 1.0d0) THEN
              CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
            ELSE
              CALL ppp_scr_identify(1.0d0,j,iset,NM,infs)
            END IF
            nb=COUNT(iset(1:NM.nobs/2) .NE. 0)
          END IF
        END IF
      END IF
    END IF

    !! step 3: try to repair cycle slips
    slip=0.d0
    IF (nb .GT. 0) THEN

      DO i=1,nb
        isat=NM.ipob(itg(iset(i)),1)
        ifreq=NM.ipob(itg(iset(i)),4)
        OB.flag(isat,ifreq)=1
      END DO

      IF (CKF.cobs(1:2) .EQ. 'IF') THEN
        CALL chsip(CKF,SAT,OB,OBH,AM,NM)
      ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
        !CALL cyslip(CKF,SAT,OB,OBH,AM,NM)
      END IF
      DO i=1,nb
        isat=NM.ipob(itg(iset(i)),1)
        ifreq=NM.ipob(itg(iset(i)),4)
        WRITE(camb,'(A4,I1)') 'AMBL',ifreq
        iamb=pointer_string(OB.npar,OB.pname,camb)
        IF (OB.flag(isat,ifreq) .EQ. 0)  THEN ! repaired
          iset(i)=-iset(i) ! flag to indicate repaired cycle slips
          j=OB.ltog(iamb,isat)-NM.npc
          IF (CKF.cobs(1:2) .EQ. 'IF') THEN
            slip(i)=SAT(isat).lamdw*AM(j).wcls/(1.d0+SAT(isat).g)+SAT(isat).lamdn*AM(j).ncls
          ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
            ! for raw, only ambiguity
            slip(i)=SAT(isat).lamda(ifreq)*AM(j).wcls
          END IF
        ! unrepaired
        ELSE
          j=OB.ltog(iamb,isat)-NM.npc
          IF (CKF.lamb .EQ. .TRUE.) AM(j).ifab=0
        END IF
      END DO
    END IF

    !! step 4: adapt information system
    CALL ppp_scr_adapt(iset,slip,NM,infs)
  END IF

  !! clean lower part of NM.infs storing u vectors
  DO j=1,NM.imtx
    DO i=j+1,NM.imtx+NM.nobs+2
      infs(NM.iptx(j)+i)=0.d0
    END DO
  END DO

  !! update OBH
  IF (NM.nobs .NE. 0) THEN
    OBH=OB
    OBH.jd=CKF.mjd
    OBH.tsec=CKF.sod
  END IF

  IF (nb .EQ. -1) RETURN

  !! adaptation of information matrix & ambiguity list
  IF (nb .GT. 0) THEN
    WRITE(OUTPUT_UNIT,'(A,I7,F9.2,2F7.1)') 'BIAS DEC',CKF.mjd,CKF.sod,esig,NM.esig
    !! set new ambiguity parameters
    DO i=1,nb
      isat =NM.ipob(itg(ABS(iset(i))),1)
      ifreq=NM.ipob(itg(ABS(iset(i))),4)
      !2023-07-23 by xsy: DIA cycle slip change the oringe AMB and Qxx
      CKF%nofixsat = TRIM(CKF.nofixsat)//' '//CKF.cprn(isat)
      
      AM(NM.ns+i).ifab =0
      !xsy
      AM(NM.ns+i).iobs =0
      AM(NM.ns+i).elev =OB.elev(isat)*RAD2DEG
      AM(NM.ns+i).psat =isat
      AM(NM.ns+i).ifreq =ifreq
      WRITE(AM(NM.ns+i).pname,'(A4,I1)') 'AMBL',ifreq
      iamb=pointer_string(OB.npar,OB.pname,AM(NM.ns+i).pname)
      AM(NM.ns+i).ptime=CKF.mjd+CKF.sod/86400.d0
      AM(NM.ns+i).xini =0.d0
      AM(NM.ns+i).xcor =0.d0
      AM(NM.ns+i).xsig =0.d0
      AM(NM.ns+i).abin =0.d0
      AM(NM.ns+i).weig =0.d0
      AM(NM.ns+i).heweig=0.d0    !xsy
      AM(NM.ns+i).eeweig=0.d0    !xsy
      AM(NM.ns+i).eweig =0.d0    !xsy
      AM(NM.ns+i).xrwl =0.d0
      AM(NM.ns+i).xswl =0.d0
      AM(NM.ns+i).xrewl =0.d0    !xsy
      AM(NM.ns+i).xsewl =0.d0    !xsy
      AM(NM.ns+i).xreewl =0.d0   !xsy
      AM(NM.ns+i).xseewl =0.d0   !xsy
      AM(NM.ns+i).xrhewl =0.d0   !xsy
      AM(NM.ns+i).xshewl =0.d0   !xsy
      AM(NM.ns+i).hewcls=0.5d0   !xsy
      AM(NM.ns+i).eewcls=0.5d0   !xsy
      AM(NM.ns+i).ewcls =0.5d0   !xsy
      AM(NM.ns+i).wcls =0.5d0
      AM(NM.ns+i).ncls =0.5d0

      IF (OB.flag(isat,ifreq) .EQ. 0) THEN
        !2023-07-22 by xsy: after repair cycle slip, iobs should not be 1 beacuase can affect REFsat selection
        AM(NM.ns+i).iobs=1
        DO Ins=1,NM%ns
          InsSat=AM(Ins)%psat
          InsFreq=AM(Ins)%ifreq
          IF(InsSat.EQ.isat .AND. InsFreq.EQ.ifreq) AM(NM.ns+i).iobs=AM(Ins).iobs
        END DO

        AM(OB.ltog(iamb,isat)-NM.npc).wcls=0.5d0
        AM(OB.ltog(iamb,isat)-NM.npc).ncls=0.5d0
      ELSE
        ! if unrepaired, original ambiguity will be eliminated
        !xsy:注意DIA引入的周跳，iobs设置负值，方便索引
        AM(NM.ns+i).iobs  =-(OB.ltog(iamb,isat)-NM.npc) !DIA new amb set as '-'
        OB.ltog(iamb,isat)=NM.imtx+i !DIA old ambiguity's ltog flag has changed,同一个模糊度的ltog值发生变化，新进来的匹配新的ltog，旧的同样被修改了，因此可以凭此删除旧模糊度
      END IF
      IF (CKF.liar) THEN
        AM(NM.ns+i).abhewl=0.5d0 !xsy
        AM(NM.ns+i).abeewl=0.5d0 !xsy
        AM(NM.ns+i).abewl=0.5d0  !xsy
        AM(NM.ns+i).abwl=0.5d0
        AM(NM.ns+i).abnl=0.5d0
      END IF
      WRITE(OUTPUT_UNIT,'(A,A3,2X,A,I1)') '---- SAT ',CKF.cprn(isat),'Carrier-phase bias on Freq. #',ifreq
    END DO
    !! ambiguity list
    NM.ns=NM.ns+nb
    NM.imtx=NM.imtx+nb
  ELSE
    WRITE(OUTPUT_UNIT,'(A,I7,F9.2,F7.1)') 'BIAS UNK', CKF.mjd, CKF.sod, NM.esig
  END IF

  !! estimate new ambiguities & remove out-of-date ambiguities
  IF (nb .GT. 0) THEN
    xcor=0.d0
    xest=0.d0
    DO i=NM.imtx,NM.npc+1,-1
      dump=infs(NM.iptx(NM.imtx+1)+i)
      DO j=i+1,NM.imtx
        dump=dump-infs(NM.iptx(j)+i)*xcor(j)
      END DO
      xcor(i)=dump/infs(NM.iptx(i)+i)
      xest(i)=xcor(i)+AM(i-NM.npc).xini
    END DO

    !! set initial value for new ambiguities
    DO i=1,NM.ns
      IF (AM(i).iobs .LT. 0) THEN
        AM(i).xini=xest(NM.npc-AM(i).iobs)
        IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
          isat=AM(i).psat
          READ(AM(i).pname,'(4X,I1)') ifreq
          AM(i).xini=xest(NM.npc-AM(i).iobs)*SAT(isat).lamda(ifreq)
        END IF
        !引入新的模糊度
        AM(i).iobs=1
      END IF
    END DO

    !! remove original ambiguities that came across slips
    nbias=0
    DO i=1, NM.ns
      isat=AM(i).psat
      iamb=pointer_string(OB.npar,OB.pname,AM(i).pname)
      IF (OB.ltog(iamb,isat) .NE. NM.npc+i) THEN
        nbias=nbias+1
        ibias(nbias)=i
      END IF
    END DO
    IF (nbias .GT. 0) THEN
      CALL ppp_del_ambi(nbias,ibias,AM,NM,infs)

      !! adjust remaining pointers to information matrix
      DO i=1,NM.ns
        iamb=pointer_string(OB.npar,OB.pname,AM(i).pname)
        OB.ltog(iamb,AM(i).psat)=NM.npc+i
      END DO

      !! to upper triangle for information matrix
      CALL ppp_msu_update(.FALSE.,0,NM,infs)
  
    END IF
  END IF

  RETURN

END SUBROUTINE
