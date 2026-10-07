!
!! purpose  : remove non-active ambiguities
!! parameter:
!!    input : nbias,ibias -- bias to be removed
!!            NM,NW -- primary & auxiliary information matrix
!!            AM    -- undifferenced ambiguities
!! author   : Geng J
!! created  : Sep 8 2011
!
SUBROUTINE fcb_del_ambi(nbias,ibias,CKF,SIT,PM,NM,infs)
!*
USE info
USE const
USE ckdctrl
USE station
IMPLICIT NONE

!*
! The arguments
!!---------------------
INTEGER(IT) :: nbias, ibias(nbias)
TYPE(SITE) :: SIT
TYPE(CKDCFG) :: CKF
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
REAL(RL) :: infs(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,j,k,ic,ik
  REAL(RL) :: s,beta,u(MAXPARSIT),gama,dump

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!--------------------------

  !! output removed wide-lane ambiguities
  IF (CKF.imw .EQ. -1) THEN
    CKF.imw=get_valid_unit(10)
    OPEN(UNIT=CKF.imw,FILE=CKF.flnmw)
  END IF

  IF (CKF.cobs(1:2) .EQ. 'IF') THEN
    DO i=1,nbias
      ic=ibias(i)
      IF (PM(ic).iobs.GT.0 .AND. (PM(ic).ptime(2)-PM(ic).ptime(1))*86400.d0.GT.CKF.minsec_common) THEN
        WRITE(CKF.imw,'(A4,1X,A3,1X,2F22.6,2F18.10,2F9.4,F6.1)') SIT.name,CKF.cprn(PM(ic).pcode(2)),&
               PM(ic).abwl,PM(ic).abewl,PM(ic).ptime(1:2),PM(ic).sigw,PM(ic).sigew,PM(ic).elev/PM(ic).iobs*RAD2DEG
      END IF
    END DO
  END IF

  !! remove a column in the information matrix
  DO i=1,nbias
    ic=ibias(i)

    !! elemental transformation vector
    CALL fcb_ele_hht(NM.imtx,infs(NM.iptx(ic)+1),s,u,beta)

    !! apply to remaining columns and compact matrix
    DO j=1,NM.imtx+1
      IF (j.EQ.ic) CYCLE
      ik=j
      IF (j.GT.ic) THEN
        ik=j-1
        IF (j .NE. NM.imtx+1) PM(ik)=PM(j)
      END IF
      gama=0.d0
      DO k=1,NM.imtx
        IF (u(k).EQ.0.d0 .OR. infs(NM.iptx(j)+k).EQ.0.d0) CYCLE
        gama=gama+u(k)*infs(NM.iptx(j)+k)
      END DO
      gama=gama*beta
      DO k=1,NM.imtx-1
        infs(NM.iptx(ik)+k)=infs(NM.iptx(j)+k+1)+gama*u(k+1)
      END DO

      !! remove the last element of this column
      infs(NM.iptx(ik)+NM.imtx)=0.d0
    END DO

    !! remove the last column
    DO j=1,NM.imtx
      infs(NM.iptx(NM.imtx+1)+j)=0.d0
    END DO

    !! decrease number of ambiguities
    NM.ns  =NM.ns  -1
    NM.imtx=NM.imtx-1

    !! revise ibias accordingly
    DO j=i+1,nbias
      IF (ibias(j).GT.ic) ibias(j)=ibias(j)-1
    END DO
  END DO

  RETURN

END SUBROUTINE
