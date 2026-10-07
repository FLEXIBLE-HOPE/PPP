!!
!! purpose   : prediction of SRIF for p and x parameters
!! parameters: NM,PM -- Information matrix & parameters
!!             infs  -- one-dimensional information matrix
!! auther    : Geng J
!! reference : Bierman GJ  1972. pp 142 (3.2) and Xiaolei Dai (2.4-16,2.4-17,2.4-18)
!*

SUBROUTINE ppp_tim_update_block(jd,sod,dintv,NM,PM,infs,SAT)
!*
USE info
USE satellite
USE f95_precision
USE lapack95
USE blas95
USE mkl_service
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------
TYPE(INFM) :: NM
TYPE(PRMT) :: PM(1:*)
INTEGER(IT) :: jd
REAL(RL) :: sod,infs(1:*),dintv
TYPE(SATE) :: SAT

  !*
  ! The local variables
  !!-----------------------
  LOGICAL(LG) :: lfirst
  REAL(RL), ALLOCATABLE :: infmx(:,:),rw(:,:),zw(:,:),w(:,:),phi(:,:)

  INTEGER(IT) :: i,ir,ic,j,k,isat,iw,minut,iupc,np
  REAL(RL) :: s,beta,gama,v(MAXPARSIT),u(MAXPARSIT)

  !*
  ! Start the exectuable code
  !!-----------------------------

  np=0
  DO i=1, NM.np
    iw=INDEX(PM(i).pname,':')
    IF (iw .NE. 0) THEN
      READ(PM(i).pname(iw+1:),*) minut
      IF ((jd-PM(i).ptime(2))*86400.d0+sod .LT. -MAXWND) CYCLE
    END IF
    np=np+1
  END DO

  ALLOCATE(w(np,1),STAT=i)
  IF (i .NE. 0) GOTO 100
  w=0.d0

  ALLOCATE(phi(np,np),STAT=i)
  IF (i .NE. 0) GOTO 100
  phi=0.d0

  ALLOCATE(rw(np,np),STAT=i)
  IF (i .NE. 0) GOTO 100
  rw=0.d0

  ALLOCATE(infmx(np+NM.imtx,np+NM.imtx+1),STAT=i)
  IF (i .NE. 0) GOTO 100
  infmx=0.d0

  !! orbit integeration to predict the state
  IF (SAT.npar.NE.0 .AND. PM(1).pname(1:5).EQ.'PXSAT') THEN
    !! estimated parameters at epoch i
    DO i=1, SAT.npar
      SAT.x(i)=PM(i).xest
      ! WRITE(*,'(A10,1X,3(E12.6,2X))')PM(i).pname(1:10),PM(i).xini,PM(i).xcor,PM(i).xest
    END DO

    !! predicted position and velocity
    CALL ppp_oi(jd,sod,dintv,SAT,SAT.x(1:6))
  END IF

  !! time update
  DO i=1, np
    PM(i).zw=0.d0
    !! piece-wise parameters
    iw=INDEX(PM(i).pname,':')
    IF (iw .NE. 0) THEN
      READ(PM(i).pname(iw+1:),*) minut
      IF ((jd-PM(i).ptime(2))*86400.d0+sod .LT. -MAXWND) CYCLE
    END IF

    !! for receiver clock, kinematic position, slant ionosphere
    !! SISRE parameter
    IF (PM(i).pname(4:6).EQ.'CLK' .OR. PM(i).pname(1:4).EQ.'STAP' .OR. &
        PM(i).pname(1:3).EQ.'SID') THEN
      PM(i).xini=PM(i).xest+PM(i).map*PM(i).xcor
    ELSE
      !! for ztd
      PM(i).xini=PM(i).map*PM(i).xest
    END IF
    PM(i).zw=-PM(i).rw*PM(i).map*PM(i).xcor

    rw(i,i)=PM(i).rw
    phi(i,i)=PM(i).map
    w(i,1)= PM(i).map*PM(i).xcor
  END DO


  IF (SAT.npar.NE.0 .AND. PM(1).pname(1:5).EQ.'PXSAT') THEN
    
    ! state transformation matrix
    DO i=1, SAT.npar
      w(i,1)=0.d0
      DO j=1, SAT.npar
        phi(i,j)=SAT.phi((j-1)*6+i)
        w(i,1)=w(i,1)+phi(i,j)*PM(j).xcor
      END DO
    END DO

    !! JG: check the xini for all parameters
    DO i=1, SAT.npar
      PM(i).xini=SAT.x(i)
      ! IF (i .GT. 6) w(i,1)=PM(i).xcor
      IF (i .GT. 6) w(i,1)=0.d0
      ! w(i,1)=0.d0
    END DO
  END IF

  !CALL matmpy(-rw,phi,infmx(1:np,1:np),np,np,np)
  !CALL matmpy( rw,w,infmx(1:np,NM.imtx+np+1),np,np,1)
  CALL GEMM(RW(1:np,1:np),PHI(1:np,1:np),infmx(1:np,1:np),'N','N',-1.d0,0.d0)
  CALL GEMM(RW(1:np,1:np),W(1:np,1:1),infmx(1:np,NM.imtx+np+1:NM.imtx+np+1),'N','N',-1.d0,0.d0)

  DO i=1, np
    DO k=np+1, np+np
      infmx(i,k)=rw(i,k-np)
      infmx(i+np,k-np)=NM.infs(i,k-np)
    END DO
  END DO

  DO i=np+1, NM.imtx+np
    DO j=np+np+1, NM.imtx+np
      infmx(i,j)=NM.infs(i-np,j-np)
    END DO
    infmx(i,NM.imtx+np+1)=NM.infs(i-np,NM.imtx+1)
  END DO
  
  CALL GEQRF(infmx)

  ! !! Householder transformation to pre-eliminate process parameters
  DO i=1, np
    iupc=0
    iw=INDEX(PM(i).pname,':')
    IF (iw .NE. 0) THEN
      READ(PM(i).pname(iw+1:),*) minut
      IF ((jd-PM(i).ptime(2))*86400.d0+sod .LT. -MAXWND) iupc=1
    END IF

    !! change time tags for ztd
    IF (iupc .EQ. 0) THEN
      PM(i).ptime(1)=PM(i).ptime(1)+minut/1440.d0
      PM(i).ptime(2)=PM(i).ptime(2)+minut/1440.d0
    END IF
      
  END DO

  DO i=1, NM.imtx
    NM.infs(i,1:NM.imtx+1)=infmx(i+np,np+1:np+NM.imtx+1)
  END DO

  IF (SAT.npar.NE.0 .AND. PM(1).pname(1:5).EQ.'PXSAT') THEN
    SAT.phi=0.d0
    DO i=1, 6
      SAT.phi(6*(i-1)+i)=1.d0
    END DO
  END IF

  DEALLOCATE(w)
  DEALLOCATE(rw)
  DEALLOCATE(phi)
  DEALLOCATE(infmx)

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_time_update): memory allocation'
  CALL exit(1)

END SUBROUTINE
