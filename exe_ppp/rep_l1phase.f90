!!
!! purpose  : repair cycle slips on L1 carrier-phase
!! parameter:
!!    input : nprn -- # of satellites
!!            irp  -- reference satellite
!!            pard -- partial derivatives of positions
!!            comc -- epoch-differenced observed minus computed values for carrier-phase
!!            wcls -- integer cycle slips for widelane
!!    output: ncls -- integer cycle slips for L1 carrier-phase
!! author   : Geng J
!! created  : Aug 29 2011
!!
!*
SUBROUTINE rep_l1phase(SAT,nprn,npar,irp,bsip,pard,pomc,comc,wcls,ncls)
!!
!*
USE par
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------------
TYPE(SATE) :: SAT(MAXSAT)
INTEGER(IT) :: nprn,npar,irp,bsip(MAXSAT)
REAL(RL) :: wcls(MAXSAT),ncls(MAXSAT),pard(MAXPARSIT,MAXSAT)
REAL(RL) :: pomc(2,MAXSAT),comc(2,MAXSAT)

  !*
  ! The local variables
  !!-----------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,j,k,nost,nsip
  REAL(RL) :: wwl,wl1,det,ltpl,tmp,chi,rto,inia(MAXSAT)
  REAL(RL), POINTER :: amat(:,:),nmat(:,:),lmat(:),esti(:),fixi(:)

  DATA lfirst/.TRUE./
  SAVE lfirst,wwl,wl1

  !*
  ! The function used
  !!----------------------
  INTEGER(IT) :: pointer_int

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! initialization
  IF (lfirst) THEN
    lfirst=.FALSE.
    wwl=1.d0
    wl1=40.d0
  END IF

  ! number of undifferenced cycle slips
  nsip=COUNT(bsip(1:nprn) .GT. 0)
  ! number of observed satellites
  nost=COUNT(comc(1,1:nprn) .NE. 0.d0)
  ! not enough measurements
  IF (npar .EQ. 0) THEN
    IF (nost .LT. nsip) RETURN
  ELSE
    IF ((nost-1)*2 .LT. npar+nsip) RETURN
  END IF

  !! array allocation
  ALLOCATE(amat(nost-1,npar+nsip+2))
  ALLOCATE(nmat(npar+nsip,npar+nsip))
  ALLOCATE(lmat(npar+nsip))
  ALLOCATE(esti(npar+nsip))  ! estimates
  ALLOCATE(fixi(nsip))    ! integers for cycle slips

  !! a priori integers of undifferenced cycle slips
  j=0
  inia=0.d0
  DO i=1,nprn
    IF (bsip(i) .GT. 0) THEN
      j=j+1
      bsip(i)=j
      inia(i)=NINT(comc(1,i)/SAT(i).lamda(1)-pomc(1,i)/SAT(i).lamda(1))
    END IF
  END DO
  j=pointer_int(nprn,bsip,-1)
  IF (j .NE. 0) inia(j)=ncls(j)

  !! form satellite pair
  amat(1:nost-1,1:npar+nsip+2)=0.d0
  nost=0
  DO i=1,nprn
    IF (i.EQ.irp .OR. comc(1,i).EQ.0.d0) CYCLE
    nost=nost+1
    !! partial derivatives
    DO j=1,npar
      amat(nost,j)=pard(j,i)-pard(j,irp)
    END DO
    tmp=0.d0
    !! WL fixed for i and irp
    IF (bsip(i) .GT. 0) THEN
      amat(nost,bsip(i)+npar)=SAT(i).lamda(1)
      tmp=wcls(i)
    END IF
    det=0.d0
    IF (bsip(irp) .GT. 0) THEN
      amat(nost,bsip(irp)+npar)=-SAT(irp).lamda(1)
      det=wcls(irp)
    END IF
    !! observed minus computed measurements
    ! unambiguous wide-lane
    amat(nost,npar+nsip+1)=SAT(i).lamdw*(comc(1,i)/SAT(i).lamda(1)-comc(2,i)/SAT(i).lamda(2)-tmp)- &
                         SAT(irp).lamdw*(comc(1,irp)/SAT(irp).lamda(1)-comc(2,irp)/SAT(irp).lamda(2)-det)
    ! L1 carrier-phase
    amat(nost,npar+nsip+2)=comc(1,i)-inia(i)*SAT(i).lamda(1)-comc(1,irp)+inia(irp)*SAT(irp).lamda(1)
  END DO

  !! 1. form normal equation
  nmat(1:npar+nsip,1:npar+nsip)=0.d0
  lmat(1:npar+nsip)=0.d0
  DO j=1,npar+nsip
    DO i=j,npar+nsip
      DO k=1,nost
        nmat(i,j)=nmat(i,j)+amat(k,i)*amat(k,j)*wl1
        IF (i.LE.npar .AND. j.LE.npar) nmat(i,j)=nmat(i,j)+amat(k,i)*amat(k,j)*wwl
      END DO
    END DO
    DO k=1,nost
      lmat(j)=lmat(j)+amat(k,j)*wl1*amat(k,npar+nsip+2)
      IF (j.LE.npar) lmat(j)=lmat(j)+amat(k,j)*wwl*amat(k,npar+nsip+1)
    END DO
  END DO
  ltpl=0.d0
  DO k=1,nost
    ltpl=ltpl+wl1*amat(k,npar+nsip+2)**2
    IF (npar .NE. 0) ltpl=ltpl+wwl*amat(k,npar+nsip+1)**2
  END DO

  !! 2. invert normal matrix
  DO j=1,npar+nsip  ! fill in upper part of nmat
    DO i=1,j-1
      nmat(i,j)=nmat(j,i)
    END DO
  END DO
  CALL matinv(nmat,npar+nsip,npar+nsip,det)
  IF (det .EQ. 0.d0)  THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(rep_l1phase): matrix inversion singularity'
    CALL exit(1)
  END IF

  !! 3. parameter estimates
  esti(1:npar+nsip)=0.d0
  DO i=1,npar+nsip
    DO j=1,npar+nsip
      esti(i)=esti(i)+nmat(i,j)*lmat(j)
    END DO
    ltpl=ltpl-esti(i)*lmat(i)
  END DO

  !! 4. search integer cycle slips
  fixi(1:nsip)=0.5d0      ! 0.5 indicates no fixing
  CALL rep_integer(npar,nost,nsip,esti,ltpl,nmat,fixi,chi,rto)
  DO i=1,nprn
    IF (bsip(i).GT.0 .AND. fixi(bsip(i)).NE.0.5d0) THEN
      ncls(i)=inia(i)+fixi(bsip(i))
    END IF
  END DO

  !! clean memory
  DEALLOCATE(amat)
  DEALLOCATE(nmat)
  DEALLOCATE(lmat)
  DEALLOCATE(esti)
  DEALLOCATE(fixi)

  RETURN

END SUBROUTINE
