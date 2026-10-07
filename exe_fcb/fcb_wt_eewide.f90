SUBROUTINE fcb_wt_eewide(CKF,UPD)
!!
!*
USE ckdctrl
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!-----------------------
TYPE(CKDCFG) :: CKF
TYPE(FCB) :: UPD

  !*
  ! The local variables
  !!----------------------
  LOGICAL(LG) :: lout
  INTEGER(IT) :: lfn,i,k,ilast,iy,imon,id,ih,im,isys
  REAL(RL) :: is

  DATA lfn /0/
  SAVE lfn

  !*
  ! The function called
  !!----------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfn .EQ. 0) THEN
    lfn=get_valid_unit(10)
    OPEN(lfn,file=CKF.flnfeewl)
  END IF

  !! write narrow-lane fractional parts for each satellite
  k=0
  DO i=1,CKF.nprn
    IF (UPD.eewfcb(i) .EQ. 10.d0) CYCLE
    k=k+1
    ilast=i
  END DO
  IF (k .EQ. 0) RETURN
  lout=.FALSE.
  i=1
  DO WHILE(i .LE. CKF.nprn)
    IF (.NOT.lout) THEN
      lout=.TRUE.
      CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,is)
      WRITE(lfn,'(A3,I5,4I3,F10.6,I7,F10.2,I4,2X,<MAXSYS>(A3,1X))') 'TIM',iy,imon,id,ih,im,is,CKF.mjd,CKF.sod,k,(UPD.updrefsat(isys,4),isys=1,MAXSYS)
    END IF
    DO WHILE(i .LE. CKF.nprn)
      IF (UPD.eewfcb(i) .NE. 10.d0) THEN
        IF (i .EQ. ilast) THEN
          WRITE(lfn,'(A3,2F8.3)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        ELSE
          WRITE(lfn,'(A3,2F8.3,$)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        END IF
        EXIT
      END IF
      i=i+1
    END DO
    i=i+1
    DO WHILE(i .LE. CKF.nprn)
      IF (UPD.eewfcb(i) .NE. 10.d0) THEN
        IF (i .EQ. ilast) THEN
          WRITE(lfn,'(3X,A3,2F8.3)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        ELSE
          WRITE(lfn,'(3X,A3,2F8.3,$)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        END IF
        EXIT
      END IF
      i=i+1
    END DO
    i=i+1
    DO WHILE(i .LE. CKF.nprn)
      IF (UPD.eewfcb(i) .NE. 10.d0) THEN
        IF (i .EQ. ilast) THEN
          WRITE(lfn,'(3X,A3,2F8.3)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        ELSE
          WRITE(lfn,'(3X,A3,2F8.3,$)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        END IF
        EXIT
      END IF
      i=i+1
    END DO
    i=i+1
    DO WHILE(i .LE. CKF.nprn)
      IF (UPD.eewfcb(i) .NE. 10.d0) THEN
        WRITE(lfn,'(3x,a3,2f8.3)') CKF.cprn(i),UPD.eewfcb(i),UPD.eewsl(i)
        EXIT
      END IF
      i=i+1
    END DO
    i=i+1
  END DO

  RETURN

END SUBROUTINE
