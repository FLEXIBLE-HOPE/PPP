!*
SUBROUTINE oi_meg_flnorb(CKF,lfn,nepo)
!!
!*
USE orbit
USE satellite
IMPLICIT NONE

!*
! The arguments
!!--------------------------
INTEGER(IT) :: lfn
TYPE(ORBCFG) :: CKF
INTEGER(IT) :: nepo

  !*
  ! The local variables
  !!----------------------------
  TYPE(ORBHDR) :: OH
  TYPE(SATEPAR) :: ICS
  INTEGER(IT) :: k,j,lfnorb,lfnsat
  INTEGER(IT) :: isat,iepo,npar

  REAL(RL) :: t,funct(MAXSAT,MAXEQUS)
  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  lfnorb=get_valid_unit(10)
  OPEN(UNIT=lfnorb,FILE=CKF.flnorb,FORM='UNFORMATTED')

  OH.nprn=CKF.nprn; OH.cprn=CKF.cprn; OH.mjd0=CKF.mjd0
  OH.mjd1=CKF.mjd1; OH.rmjd=CKF.rmjd; OH.sod0=CKF.sod0;
  OH.sod1=CKF.sod1; OH.rsod=CKF.rsod; OH.dintv=CKF.dintv

  npar=1
  IF (CKF.lpart .EQ. .TRUE.) npar=CKF.nequ/6
  CKF.nequ=npar*3
  OH.nequ=CKF.nequ

  !! Control structure
  WRITE(lfnorb) OH

  DO isat=1, CKF.nprn
    lfnsat=lfn+isat
    READ(lfnsat) ICS
    WRITE(lfnorb) ICS
  END DO

  !! Pos&vel and partial records
  DO iepo=1, nepo
    DO isat=1, CKF.nprn
      lfnsat=lfn+isat
      READ(lfnsat) t,((funct(isat,(k-1)*6+j),j=1,6), k=1,npar)
    END DO
    !! form 6 to 3
    WRITE(lfnorb) (((funct(isat,(k-1)*6+j),j=1,3),k=1,npar),isat=1,CKF.nprn)
    
  END DO

  DO isat=1, CKF.nprn
    lfnsat=lfn+isat
    CLOSE(lfnsat)
  END DO
  CLOSE(lfnorb)

  RETURN

END SUBROUTINE
