!
!! purpose  : measurement processing of the square root information filter
!!              inforamtion     R  Z             R  Z
!!              observations    A  Z     ==>        V
!! paramters: lsavu -- save u vectors for data screening
!!            lsavr -- save residuals
!!            nobs  -- number of incorporated observations
!!            NM   -- information matrix
!!            infs -- one-dimentional normal equation
!! author   : Geng J
!! reference: Bierman G. J., 1977
!
SUBROUTINE ppp_msu_update(lsavu,nobs,NM,infs)
!*
USE info
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------
TYPE(INFM) :: NM
LOGICAL(LG) :: lsavu
INTEGER(IT) :: nobs
REAL(RL) :: infs(1:*)

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: s,u(MAXPARSIT+MAXSAT*2),beta,gama

  !*
  ! Start the exectuable code
  !!--------------------------

  !! Householder transformation
  DO j=1,NM.imtx

    !! whether necessary for this column
    IF (ALL(infs(NM.iptx(j)+j+1:NM.iptx(j)+NM.imtx+nobs) .EQ. 0.d0)) THEN
      IF (lsavu) THEN
        infs(NM.iptx(j)+NM.imtx+nobs+1)=0.d0
        infs(NM.iptx(j)+NM.imtx+nobs+2)=0.d0
      END IF
    ELSE

      !! get elemental transformation vector
      CALL ppp_ele_hht(NM.imtx+nobs-j+1,infs(NM.iptx(j)+j),s,u,beta)
      infs(NM.iptx(j)+j)=s

      !! transform the next columns
      DO k=j+1,NM.imtx+1
        gama=0.d0
        DO i=j,NM.imtx+nobs
          IF (u(i-j+1).EQ.0.d0 .OR. infs(NM.iptx(k)+i).EQ.0.d0) CYCLE
          gama=gama+u(i-j+1)*infs(NM.iptx(k)+i)
        END DO
        IF (gama .EQ. 0.d0) CYCLE
        gama=beta*gama
        DO i=j,NM.imtx+nobs
          IF (u(i-j+1) .EQ. 0.d0) CYCLE
          infs(NM.iptx(k)+i)=infs(NM.iptx(k)+i)+gama*u(i-j+1)
        END DO
      END DO

      !! store u vectors for quality control
      IF (lsavu) THEN
        IF (NM.imtx+nobs+2 .GT. NM.nmtx) THEN
          WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_msu_update): information matrix'
          CALL exit(1)
        END IF

        !! store elemental transformation vector
        DO i=j+1, NM.imtx+nobs+1
          infs(NM.iptx(j)+i)=u(i-j)
        END DO
        infs(NM.iptx(j)+NM.imtx+nobs+2)=beta
      ELSE

        !! clean current column
        DO i=j+1,NM.imtx+nobs
          infs(NM.iptx(j)+i)=0.d0
        END DO
      END IF
    END IF

  !! next column
  END DO

  ! residual include phase and code
  !! fill in postfit residuals
  NM.esig=(NM.nobs-nobs)*NM.esig**2
  DO i=1,nobs
    NM.resi(NM.nobs-nobs+i)=infs(NM.iptx(NM.imtx+1)+NM.imtx+i)
    infs(NM.iptx(NM.imtx+1)+NM.imtx+i)=0.d0
    NM.esig=NM.esig+NM.resi(NM.nobs-nobs+i)**2
  END DO
  IF (NM.nobs .GT. 0) NM.esig=DSQRT(NM.esig/NM.nobs)
  
  ! DO i=1,NM.imtx+1
  !   WRITE(6000,'(<NM.imtx>(F8.3,2X))')(infs(NM.iptx(i)+j),j=1,NM.imtx)
  ! END DO
  ! WRITE(6000,*)'  '

  RETURN

END SUBROUTINE
