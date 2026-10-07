!
!! purpose  : wide-lane cycle-slip resolution
!! parameter:
!!    input : nprn -- # of satellites
!!            irp  -- reference satellite
!!            bsip -- where are the slips
!!            elev -- satellite elevations in radians
!!            pard -- partial derivatives of positions
!!            pomc,comc -- epoch-differenced observed minus computed values for pseudorange & carrier-phase
!!    output: wcls -- integer wide-lane cycle slips
!! author   : Geng J
!! created  : Dec 31 2010
!
SUBROUTINE rep_widelane(SAT,nprn,npar,irp,bsip,pard,pomc,comc,wcls)
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
REAL(RL) :: wcls(MAXSAT),pard(MAXPARSIT,MAXSAT)
REAL(RL) :: pomc(2,MAXSAT),comc(2,MAXSAT)

  !*
  ! The local variables
  !!-----------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,j,k,nost,nsip
  REAL(RL) :: wrng,wphs,det,ltpl,chi,rto,inia(MAXSAT)
  REAL(RL), POINTER :: amat(:,:),nmat(:,:),lmat(:),esti(:),fixi(:)

  DATA lfirst/.TRUE./
  SAVE lfirst,wrng,wphs

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_int

  !*
  ! Start the exectuable code
  !!-------------------------------

  !! initialization
  IF (lfirst) THEN
    lfirst=.FALSE.
    wrng=1.d-3
    wphs=1.d3
  END IF

  ! number of undifferenced cycle slips
  nsip=COUNT(bsip(1:nprn) .GT. 0)
  ! number of observed satellites
  nost=COUNT(comc(1,1:nprn) .NE. 0.d0)
  ! not enough measurements
  IF (npar .EQ. 0) THEN
    IF (nost-1 .LT. nsip) RETURN
  ELSE
    IF ((nost-1)*2 .LT. npar+nsip) RETURN
  END IF

  !! array allocation
  ALLOCATE(amat(nost-1,npar+nsip+2))
  ALLOCATE(nmat(npar+nsip,npar+nsip))
  ALLOCATE(lmat(npar+nsip))
  ! estimates
  ALLOCATE(esti(npar+nsip))
  ! integers for cycle slips
  ALLOCATE(fixi(nsip))

  !! a priori integers of undifferenced cycle slips
  j=0
  inia=0.d0
  DO i=1, nprn
    IF (bsip(i) .GT. 0) THEN
      j=j+1
      bsip(i)=j
      inia(i)=NINT(comc(1,i)/SAT(i).lamda(1)-comc(2,i)/SAT(i).lamda(2)- &
                   pomc(1,i)/SAT(i).lamda(1)+pomc(2,i)/SAT(i).lamda(2))
    END IF
  END DO
  j=pointer_int(nprn,bsip,-1)
  IF (j .NE. 0) inia(j)=wcls(j)

  !! form satellite pair
  amat(1:nost-1,1:npar+nsip+2)=0.d0
  ! number of carrier-phase observations
  nost=0
  DO i=1,nprn
    IF (i.EQ.irp .OR. comc(1,i).EQ.0.d0) CYCLE
    nost=nost+1
    !! partial derivatives
    DO j=1,npar
      amat(nost,j)=pard(j,i)-pard(j,irp)
    END DO
    IF (bsip(i) .GT. 0) THEN
      amat(nost,bsip(i)+npar)=SAT(i).lamdw
    END IF
    IF (bsip(irp) .GT. 0) THEN
      amat(nost,bsip(irp)+npar)=-SAT(irp).lamdw
    END IF
    !! observed minus computed measurements
    !amat(nost,npar+nsip+1)=pomc(1,i)-pomc(1,irp) ! pseudorange
    amat(nost,npar+nsip+1)=pomc(1,i)*SAT(i).g/(SAT(i).g+1.d0)+pomc(2,i)/(1.d0+SAT(i).g)-&
                           pomc(1,irp)*SAT(irp).g/(SAT(irp).g+1.d0)-pomc(2,irp)/(1.d0+SAT(irp).g) ! pseudorange
    amat(nost,npar+nsip+2)=SAT(i).lamdw*(comc(1,i)/SAT(i).lamda(1)-comc(2,i)/SAT(i).lamda(2)-inia(i))-&
                         SAT(irp).lamdw*(comc(1,irp)/SAT(irp).lamda(1)-comc(2,irp)/SAT(irp).lamda(2)-inia(irp)) ! carrier-phase
  END DO

  !! 1. form normal equation
  nmat(1:npar+nsip,1:npar+nsip)=0.d0
  lmat(1:npar+nsip)=0.d0
  DO j=1,npar+nsip
    DO i=j,npar+nsip
      DO k=1,nost
        nmat(i,j)=nmat(i,j)+amat(k,i)*amat(k,j)*wphs
        IF (i.LE.npar .AND. j.LE.npar) nmat(i,j)=nmat(i,j)+amat(k,i)*amat(k,j)*wrng
      END DO
    END DO
    DO k=1,nost
      lmat(j)=lmat(j)+amat(k,j)*wphs*amat(k,npar+nsip+2)
      IF (j.LE.npar) lmat(j)=lmat(j)+amat(k,j)*wrng*amat(k,npar+nsip+1)
    END DO
  END DO
  ltpl=0.d0
  DO k=1,nost
    ltpl=ltpl+wphs*amat(k,npar+nsip+2)**2
    IF (npar .GT. 0) ltpl=ltpl+wrng*amat(k,npar+nsip+1)**2
  END DO

  !! 2. invert normal matrix
  DO j=1,npar+nsip  ! fill in upper part of nmat
    DO i=1,j-1
      nmat(i,j)=nmat(j,i)
    END DO
  END DO
  CALL matinv(nmat,npar+nsip,npar+nsip,det)
  IF (det .EQ. 0.d0) THEN
    WRITE(ERROR_UNIT,'(a)') '***ERROR(rep_widelane): matrix inversion singularity'
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
  ! 0.5 indicates no fixing
  fixi(1:nsip)=0.5d0
  CALL rep_integer(npar,nost,nsip,esti,ltpl,nmat,fixi,chi,rto)
  DO i=1,nprn
    IF (bsip(i).GT.0 .AND. fixi(bsip(i)).NE.0.5d0) THEN
      wcls(i)=inia(i)+fixi(bsip(i))
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
