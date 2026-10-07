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
!! 2026-10-07 by xsy: 参考 ckdrt_irc_cr_orb_sdb 的 fcb_msu_update, 用 LAPACK
!!    dgeqrf + dormqr (BLAS3 分块) 替换手写标量三重循环, 数学等价; dgeqrf 的
!!    LAPACK 反射子存储 (v1=1, tau, 次对角存 v_k) 再还原成 ppp_scr_senvec
!!    读取的 Geng 约定 (u(1)=-s*tau, u(k)=-s*tau*v_{k-1}, beta=-1/(tau*s^2))
!!    infs 的第 j 列起于 (j-1)*nmtx+1, 即 lda=nmtx 的列主序, 可直接喂 LAPACK
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
  INTEGER(IT) :: i,j,m,n,lda,lwork,ierr
  REAL(RL) :: s,wq(1)
  REAL(RL), ALLOCATABLE :: tau(:),work(:)

  !*
  ! Start the exectuable code
  !!--------------------------

  m=NM.imtx+nobs
  n=NM.imtx
  lda=NM.nmtx

  ALLOCATE(tau(n))
  CALL dgeqrf(m,n,infs,lda,tau,wq,-1_IT,ierr)   ! workspace query
  lwork=MAX(1,INT(wq(1)))
  ALLOCATE(work(lwork))

  !! Householder transformation of the parameter columns only
  CALL dgeqrf(m,n,infs,lda,tau,work,lwork,ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A,I0)') '***ERROR(ppp_msu_update): dgeqrf failed, info=',ierr
    CALL exit(1)
  END IF

  !! apply Q^T to the right hand side column (NM.imtx+1) : dgeqrf 未包含该列
  CALL dormqr('L','T',m,1,n,infs,lda,tau,infs(n*lda+1),lda,wq,-1_IT,ierr)
  i=MAX(1,INT(wq(1)))
  IF (i .GT. lwork) THEN
    DEALLOCATE(work)
    ALLOCATE(work(i))
    lwork=i
  END IF
  CALL dormqr('L','T',m,1,n,infs,lda,tau,infs(n*lda+1),lda,work,lwork,ierr)
  IF (ierr .NE. 0) THEN
    WRITE(ERROR_UNIT,'(A,I0)') '***ERROR(ppp_msu_update): dormqr failed, info=',ierr
    CALL exit(1)
  END IF

  IF (lsavu) THEN
    IF (NM.imtx+nobs+2 .GT. NM.nmtx) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_msu_update): information matrix'
      CALL exit(1)
    END IF

    !! restore elemental transformation vector (Geng convention, for ppp_scr_senvec)
    DO j=1,n
      s=infs((j-1)*lda+j)
      IF (tau(j) .EQ. 0.d0) THEN
        infs((j-1)*lda+j+1:(j-1)*lda+NM.imtx+nobs+2)=0.d0
      ELSE
        infs((j-1)*lda+j+2:(j-1)*lda+NM.imtx+nobs+1)=-infs((j-1)*lda+j+1:(j-1)*lda+NM.imtx+nobs)*s*tau(j)
        infs((j-1)*lda+j+1)=-s*tau(j)
        infs((j-1)*lda+NM.imtx+nobs+2)=-1.d0/(tau(j)*s*s)
      END IF
    END DO
  ELSE

    !! clean current column below the diagonal
    DO j=1,n
      infs((j-1)*lda+j+1:(j-1)*lda+NM.imtx+nobs+2)=0.d0
    END DO
  END IF

  DEALLOCATE(tau,work)

  ! residual include phase and code
  !! fill in postfit residuals
  NM.esig=(NM.nobs-nobs)*NM.esig**2
  DO i=1,nobs
    NM.resi(NM.nobs-nobs+i)=infs(n*lda+NM.imtx+i)
    infs(n*lda+NM.imtx+i)=0.d0
    NM.esig=NM.esig+NM.resi(NM.nobs-nobs+i)**2
  END DO
  IF (NM.nobs .GT. 0) NM.esig=DSQRT(NM.esig/NM.nobs)

  RETURN

END SUBROUTINE