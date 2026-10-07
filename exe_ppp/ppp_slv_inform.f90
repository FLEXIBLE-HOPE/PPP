!
!! purpose   : solve information matrix in primary filter
!! parameters: SF -- filter
!!             isol  = 1, parameter only else parameter and inversion of the SRI 
!! author    : Geng J
!
SUBROUTINE ppp_slv_inform(liar,NM,PM,AM,infs,QM,invx,Qxyz)
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!---------------------
LOGICAL(LG) :: liar
TYPE(INFM) :: NM
TYPE(INVM) :: QM
TYPE(PRMT) :: PM(1:*)
TYPE(AMBT) :: AM(1:*)
REAL(RL) :: infs(1:*),invx(1:*)
REAL(RL) :: Qxyz(3,3)

  !*
  ! The local variables
  !!------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,j,k,isol,ierr
  REAL(RL) :: dummy,xcor
  INTEGER(IT), POINTER :: idx(:)
  REAL(RL), POINTER :: dx(:)

  !*
  ! Start the exectuable code
  !!------------------------------

  isol=2

  !! solve information matrix
  IF (isol .EQ. 1) THEN    ! parameter only, backward substitution

    DO i=NM.imtx,1,-1
      dummy=infs(NM.iptx(NM.imtx+1)+i)
      DO j=i+1,NM.imtx
        IF (j .LE. NM.npc) THEN
          xcor=PM(j).xcor       !npc参数
        ELSE
          xcor=AM(j-NM.npc).xcor!amb参数
        END IF
        dummy=dummy-infs(NM.iptx(j)+i)*xcor
      END DO
      IF (i .LE. NM.npc) THEN                ! parameters
        IF (PM(i).iobs .EQ. 0) THEN
          PM(i).xcor=0.d0
        ELSE
          PM(i).xcor=dummy/infs(NM.iptx(i)+i)
        END IF
        PM(i).xest=PM(i).xini+PM(i).xcor
        PM(i).xsig=0.d0
      ELSE                                   ! ambiguities
        IF (AM(i-NM.npc).iobs .EQ. 0) THEN
          AM(i-NM.npc).xcor=0.d0
        ELSE
          AM(i-NM.npc).xcor=dummy/infs(NM.iptx(i)+i)
        END IF
        AM(i-NM.npc).xsig=0.d0
      END IF
    END DO
  ELSE
    ! cofactor matrix, invert information matrix：R-1
    ALLOCATE(dx(NM.imtx*(NM.imtx+1)/2),stat=ierr)
    ALLOCATE(idx(NM.imtx),stat=ierr)
    DO i=1,NM.imtx
      idx(i)=(i-1)*i/2
    END DO
    dx=0.d0
    IF (PM(1).iobs .EQ. 0) THEN
      dx(idx(1)+1)=0.d0
    ELSE
      dx(idx(1)+1)=1.d0/infs(NM.iptx(1)+1)
    END IF
    DO j=2,NM.imtx
      !! xsy
      ! IF ((j.LE.NM.npc.AND.PM(j).iobs.EQ.0) .OR. (j.GT.NM.npc.AND.AM(j-NM.npc).iobs.EQ.0)) THEN
      !   dx(idx(j)+j)=0.d0
      ! ELSE
      !   dx(idx(j)+j)=1.d0/infs(NM.iptx(j)+j)
      ! END IF
      IF (j .le. NM.npc) THEN
        IF (PM(j).iobs .eq. 0) THEN
          dx(idx(j)+j)=0.d0
        ELSE
          dx(idx(j)+j)=1.d0/infs(NM.iptx(j)+j)
        END IF
      ELSE
        IF (AM(j-NM.npc).iobs .eq. 0) THEN
          dx(idx(j)+j)=0.d0
        ELSE
          dx(idx(j)+j)=1.d0/infs(NM.iptx(j)+j)
        END IF
      END IF
      
      DO k=1,j-1
        DO i=k,j-1
          dx(idx(j)+k)=dx(idx(j)+k)+dx(idx(i)+k)*infs(NM.iptx(j)+i)
        END DO
        dx(idx(j)+k)=-1.d0*dx(idx(j)+k)*dx(idx(j)+j)
      END DO
    END DO

    !! compute lower-trangle part of cofactor matrix: (R-1 * R-t)   xsy: 计算Qxx，协因数阵,用来固定模糊度、消模糊度
    IF (liar) THEN
      QM.invx(1:NM.imtx,1:NM.imtx)=0.d0
      DO j=1,NM.imtx
        DO i=j,NM.imtx
          DO k=i,NM.imtx
            invx(QM.idq(j)+i)=invx(QM.idq(j)+i)+dx(idx(k)+j)*dx(idx(k)+i)
          END DO
        END DO
      END DO
      QM.ntot=NM.imtx
      QM.nxyz=NM.npc
      QM.ndam=NM.ns
      QM.nfix=0
    END IF

    !! compute estimates & variance     xsy: 计算参数估值及其方差
    DO i=1,NM.imtx
      ! non-amb parameters
      IF (i .LE. NM.npc) THEN
        PM(i).xcor=0.d0
        PM(i).xsig=0.d0
        DO j=i,NM.imtx
          PM(i).xcor=PM(i).xcor+dx(idx(j)+i)*infs(NM.iptx(NM.imtx+1)+j)
          PM(i).xsig=PM(i).xsig+dx(idx(j)+i)**2   !Qxx = (Rt R)-1 = R-1 (R-1)t
        END DO
        PM(i).xest=PM(i).xini+PM(i).xcor
        PM(i).xsig=DSQRT(PM(i).xsig)*NM.esig
        ! xsy 2024-02-21: for Qxyz
        IF (INDEX(PM(i).pname,'STAPX').NE.0) THEN
           Qxyz(1,1)=invx(QM.idq(i)+i)
           Qxyz(2,1)=invx(QM.idq(i)+i+1)
           Qxyz(3,1)=invx(QM.idq(i)+i+2)
           Qxyz(2,2)=invx(QM.idq(i+1)+i+1)
           Qxyz(3,2)=invx(QM.idq(i+1)+i+2)
           Qxyz(3,3)=invx(QM.idq(i+2)+i+2)

           Qxyz(1,2)=Qxyz(2,1)
           Qxyz(1,3)=Qxyz(3,1)
           Qxyz(2,3)=Qxyz(3,2)
        END IF

        ! IF (INDEX(PM(i)%pname,'STAP').NE.0 .OR. INDEX(PM(i)%pname,'RECCLK').NE.0 .OR. INDEX(PM(i)%pname,'ZTD').NE.0 .OR. INDEX(PM(i)%pname,'ION').NE.0) THEN
        !     IF (INDEX(PM(i)%pname,'ION').NE.0 .AND. PM(i)%xini.EQ.0.D0) CYCLE
        !     WRITE(5000,'(A10,3(1X,F12.3),3X,F4.1,1X,F12.8)') '- '//PM(i)%pname,PM(i)%xini,PM(i)%xcor,PM(i)%xest,PM(i)%map,PM(i)%rw
        ! END IF
      ! LC ambiguities
      ELSE
        AM(i-NM.npc).xcor=0.d0
        AM(i-NM.npc).xsig=0.d0
        DO j=i,NM.imtx
          AM(i-NM.npc).xcor=AM(i-NM.npc).xcor+dx(idx(j)+i)*infs(NM.iptx(NM.imtx+1)+j)
          AM(i-NM.npc).xsig=AM(i-NM.npc).xsig+dx(idx(j)+i)**2
        END DO
        AM(i-NM.npc).famb=AM(i-NM.npc).xini+AM(i-NM.npc).xcor
        AM(i-NM.npc).xsig=DSQRT(AM(i-NM.npc).xsig)*NM.esig
        !write(*,*) AM(i-NM.npc).pname,AM(i-NM.npc).xcor,AM(i-NM.npc).xsig
      END IF
    END DO
    DEALLOCATE(dx)
    DEALLOCATE(idx)
  END IF

  RETURN

END SUBROUTINE
