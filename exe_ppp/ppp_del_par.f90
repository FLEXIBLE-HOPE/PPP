!
!! purpose  : remove non-active ambiguities
!! parameter:
!!    input : nbias,ibias -- bias to be removed
!!            NM,NW -- primary & auxiliary information matrix
!!            AM    -- undifferenced ambiguities
!! author   : Geng J
!! created  : Sep 8 2011
!!
!*
SUBROUTINE ppp_del_par(nbias,ibias,AM,NM,infs)
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
  INTEGER(IT) :: i,j,k,ic,id
  REAL(RL) :: fac

  !*
  ! Start the exectuable code
  !!---------------------------

  !! remove a column in the information matrix
  DO i=1,nbias
    ic=ibias(i)
    id=NM.npc+ic

    ! remove a parameter
    DO j=1,id-1

      fac=infs(NM.iptx(id)+j)/infs(NM.iptx(id)+id)
      DO k=id+1, NM.imtx+1
        IF (infs(NM.iptx(id)+id) .NE. 0.D0) THEN
          infs(NM.iptx(k)+j)=infs(NM.iptx(k)+j)-fac*infs(NM.iptx(k)+id)
        END IF
      END DO

    END DO

    ! compact the information matrix and move the ambiguity array
    DO j=id,NM.imtx
      infs(NM.iptx(j)+1:NM.iptx(j)+id-1) = infs(NM.iptx(j+1)+1:NM.iptx(j)+id-1)
      infs(NM.iptx(j)+id:NM.iptx(j)+j) = infs(NM.iptx(j+1)+id+1:NM.iptx(j+1)+j+1)
    END DO
    DO j=1,NM.imtx+1
      infs(NM.iptx(j)+NM.imtx)=0.d0
    END DO

    DO j=ic, NM.ns-1
      AM(j)=AM(j+1)
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
