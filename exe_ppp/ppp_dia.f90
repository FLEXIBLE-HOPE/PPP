!
! CREATED BY SHENGYI XU
! DIA detect cycle slips and go back to raw infs and creat former equation
!*
SUBROUTINE ppp_dia(CKF,SAT,OB,OBH,NM,infs,nb)
!!
!*
USE info
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------
TYPE(SATE) :: SAT(MAXSAT)
TYPE(RNXOBS) :: OB,OBH
TYPE(INFM) :: NM
TYPE(CKDCFG) :: CKF
REAL(RL) :: infs(1:*)
INTEGER(IT) :: nb

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j,isat,iamb,ifreq
  INTEGER(IT) :: iset(MAXFREQ*MAXSAT),itg(MAXFREQ*MAXSAT),nbias,ibias(MAXFREQ*MAXSAT)
  REAL(RL) :: esig,thre,slip(MAXFREQ*MAXSAT)

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!---------------------------

  !! whether cycle slips exist
  nb=-1     ! flag for no slips
  iset=0    ! corresponds to each carrier-phase obs
  esig=NM.esig
  thre=CKF%DiaSig
  !@CMT BY XSY: For some case,the esig suddenly gets bigger, maybe small cycle slip is exsit
  ! IF (NM%esig_before*3.LT.esig .AND. esig.GT.0.3D0 .AND. thre.GT.0.3D0) THEN
  !   thre=0.3D0
  ! END IF
  !@CMT BY XSY: This is not suitable for LUTAN BDS-3
  ! IF (esig.GT.0.12d0 .AND. esig.LT.thre .AND. NM%esig_before*3.LT.esig .AND. NM%esig_before.NE.0.d0) THEN
  !   thre=NM%esig_before*3
  ! END IF
  IF (NM.esig .GT. thre) THEN
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
    ! IF (CKF.DiaSig .LT. 0.3d0) THEN
    !   CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
    ! ELSE
    !   CALL ppp_scr_identify(0.3d0,j,iset,NM,infs)
    ! END IF
    CALL ppp_scr_identify(thre,j,iset,NM,infs)
    nb=COUNT(iset(1:NM.nobs/2) .NE. 0)

    ! IF (nb .EQ. 0) THEN
    !   IF (CKF.DiaSig .GT. 0.3d0) THEN
    !     IF (CKF.DiaSig .LT. 0.5d0) THEN
    !       CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
    !     ELSE
    !       CALL ppp_scr_identify(0.5d0,j,iset,NM,infs)
    !     END IF
    !     nb=COUNT(iset(1:NM.nobs/2) .NE. 0)

    !     IF (nb .EQ. 0) THEN
    !       IF (CKF.DiaSig .GT. 0.5d0) THEN
    !         IF (CKF.DiaSig .LT. 1.0d0) THEN
    !           CALL ppp_scr_identify(CKF.DiaSig,j,iset,NM,infs)
    !         ELSE
    !           CALL ppp_scr_identify(1.0d0,j,iset,NM,infs)
    !         END IF
    !         nb=COUNT(iset(1:NM.nobs/2) .NE. 0)
    !       END IF
    !     END IF
    !   END IF
    ! END IF

    IF (nb .GT. 0) THEN
      DO i=1,nb
        isat=NM.ipob(itg(iset(i)),1)
        ifreq=NM.ipob(itg(iset(i)),4)
        OB.flag(isat,:)=1
      END DO
    END IF
    slip=0.d0
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

  IF (nb .GT. 0) THEN
    WRITE(OUTPUT_UNIT,'(A,I7,F9.2,2F7.2)') 'BIAS DEC',CKF.mjd,CKF.sod,esig,NM.esig
    DO i=1,nb
      isat =NM.ipob(itg(ABS(iset(i))),1)
      ifreq=NM.ipob(itg(ABS(iset(i))),4)
      WRITE(OUTPUT_UNIT,'(A,A3,2X,A,I1)') '---- SAT ',CKF.cprn(isat),'Carrier-phase bias on Freq. #',ifreq
    END DO
  ELSE
    WRITE(OUTPUT_UNIT,'(A,I7,F9.2,F7.2)') 'BIAS UNK', CKF.mjd, CKF.sod, NM.esig
  END IF

  RETURN

END SUBROUTINE
