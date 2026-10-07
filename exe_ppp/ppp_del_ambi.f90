!
!! purpose  : remove non-active ambiguities
!! parameter:
!!    input : nbias,ibias -- bias to be removed
!!            NM,NW -- primary & auxiliary information matrix
!!            AM    -- undifferenced ambiguities
!! author   : Geng J
!! created  : Sep 8 2011
!
SUBROUTINE ppp_del_ambi(nbias,ibias,AM,NM,infs)
!!
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!-----------------
TYPE(INFM) :: NM
TYPE(AMBT) :: AM(1:*)
INTEGER(IT) :: nbias,ibias(nbias)
REAL(RL) infs(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,j,k,ic,ik
  REAL(RL) :: s,beta,u(MAXPARSIT+MAXSAT*2),gama

  !*
  ! Start the exectuable code
  !!---------------------------

  !! remove a column in the information matrix
  DO i=1,nbias
    ic=ibias(i)

    !! elemental transformation vector
    CALL ppp_ele_hht(NM.imtx,infs(NM.iptx(NM.npc+ic)+1),s,u,beta)

    !! apply to remaining columns and compact matrix
    DO j=1,NM.imtx+1
      IF (j .EQ. NM.npc+ic) CYCLE
      ik=j
      IF (j .GT. NM.npc+ic) THEN
        ik=j-1
        IF (j .NE. NM.imtx+1) AM(ik-NM.npc)=AM(j-NM.npc)
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
    NM.imtx=NM.imtx-1
    NM.ns=NM.ns-1

    !! revise ibias accordingly
    DO j=i+1,nbias
      IF (ibias(j) .GT. ic) ibias(j)=ibias(j)-1
    END DO
  END DO

  RETURN

END SUBROUTINE
