!
!! purpose  : screen data in real time
!! parameter:
!!    input : iset   -- problematic carrier phase
!!            NM     -- information matrix & parameters
!!            infs   -- one-dimensional array for information matrix
!! author   : Geng J
!! created  : Nov. 24, 2007
!
SUBROUTINE ppp_scr_adapt(iset,slip,NM,infs)
!*
USE info
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
INTEGER(IT) :: iset(1:*)
REAL(RL) :: infs(1:*),slip(1:*)
TYPE(INFM) :: NM

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: nb,i,j,k
  REAL(RL) :: s,beta,gama,u(MAXSAT*2)
  REAL(RL), POINTER :: ures(:)

  !*
  ! Start the exectuable code
  !!-------------------------------

  !! number of problmatic observations
  nb=COUNT(iset(1:NM.nobs/2) .NE. 0)
  IF (NM.imtx+1+NM.nobs+nb .GT. NM.nmtx) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_scr_adapt): information matrix too small'
    CALL exit(1)
  END IF

  !! save upper part of vector (imtx+1)
  IF (nb .GT. 0) THEN
    ALLOCATE(ures(NM.imtx))
    DO i=1,NM.imtx
      ures(i)=infs(NM.iptx(NM.imtx+1)+i)
    END DO
  END IF

  !! split and move sensitivity vectors of problematic measurements
  DO j=1,nb
    DO i=1,NM.imtx+nb+NM.nobs
      infs(NM.iptx(NM.imtx+j)+i)=0.d0
    END DO
    DO i=1,NM.imtx+NM.nobs
      IF (i .LE. NM.imtx) THEN
        infs(NM.iptx(NM.imtx+j)+i)   =infs(NM.iptx(NM.imtx+1+ABS(iset(j)))+i)
      ELSE
        infs(NM.iptx(NM.imtx+j)+i+nb)=infs(NM.iptx(NM.imtx+1+ABS(iset(j)))+i)
      END IF
    END DO
    IF (iset(j) .LT. 0) THEN
      infs(NM.iptx(NM.imtx+j)+NM.imtx+j)=1.d6   ! repaired cycle slips
    ELSE
      infs(NM.iptx(NM.imtx+j)+NM.imtx+j)=1.d-4
    END IF
  END DO

  !! restore vector (imtx+1)
  IF (nb .GT. 0) THEN
    DO i=1,NM.imtx+nb+NM.nobs
      infs(NM.iptx(NM.imtx+nb+1)+i)=0.d0
      IF (i.GT.NM.imtx .AND. i.LE.NM.imtx+nb) THEN
        IF (iset(i-NM.imtx).LT.0) THEN
          infs(NM.iptx(NM.imtx+nb+1)+i)=1.d6*slip(i-NM.imtx) ! repaired cycle slips
        END IF
      END IF
    END DO
    DO i=1,NM.imtx+NM.nobs
      IF (i .LE. NM.imtx) THEN
        infs(NM.iptx(NM.imtx+nb+1)+i)   =ures(i)
      ELSE
        infs(NM.iptx(NM.imtx+nb+1)+i+nb)=NM.resi(i-NM.imtx)
      END IF
    END DO
    DEALLOCATE(ures)
  END IF

  !! clean remaining sensitivity vectors
  DO j=NM.imtx+nb+1+1,NM.imtx+1+NM.nobs/2
    DO i=1,NM.imtx+NM.nobs
      infs(NM.iptx(j)+i)=0.d0
    END DO
  END DO

  !! adaptation, transform to upper triangle
  DO j=1,nb

    !! get elemental transformation vector
    CALL ppp_ele_hht(NM.nobs+nb-j+1,infs(NM.iptx(NM.imtx+j)+NM.imtx+j),s,u,beta)
    infs(NM.iptx(NM.imtx+j)+NM.imtx+j)=s

    !! transform the next columns
    DO k=j+1,nb+1
      gama=0.d0
      DO i=NM.imtx+j,NM.imtx+nb+NM.nobs
        IF (u(i-NM.imtx-j+1).EQ.0.d0 .OR. infs(NM.iptx(NM.imtx+k)+i).EQ.0.d0) CYCLE
        gama=gama+u(i-NM.imtx-j+1)*infs(NM.iptx(NM.imtx+k)+i)
      END DO
      IF (gama .EQ. 0.d0) CYCLE
      gama=beta*gama
      DO i=NM.imtx+j,NM.imtx+nb+NM.nobs
        IF (u(i-NM.imtx-j+1).EQ.0.d0) CYCLE
        infs(NM.iptx(NM.imtx+k)+i)=infs(NM.iptx(NM.imtx+k)+i)+gama*u(i-NM.imtx-j+1)
      END DO
    END DO

    !! clean lower part
    DO i=NM.imtx+j+1,NM.imtx+nb+NM.nobs
      infs(NM.iptx(NM.imtx+j)+i)=0.d0
    END DO
  END DO

  !! store new residuals after mitigating slips
  IF (nb .GT. 0) THEN
    NM.esig=0.d0
    DO i=NM.imtx+nb+1,NM.imtx+nb+NM.nobs
      NM.resi(i-NM.imtx-nb)=infs(NM.iptx(NM.imtx+nb+1)+i)
      infs(NM.iptx(NM.imtx+nb+1)+i)=0.d0
      NM.esig=NM.esig+NM.resi(i-NM.imtx-nb)**2
    END DO
    IF (NM.nobs .GT. 0) NM.esig=DSQRT(NM.esig/NM.nobs)
  END IF

  RETURN

END SUBROUTINE
