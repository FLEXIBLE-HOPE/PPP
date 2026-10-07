!
!! purpose  : identify problematic carrier phase
!! parameter:
!!    input : sigm -- threshold sigma
!!            maxcph -- maximum number of problematic phase
!!            NM,infs -- information matrix
!!    output: iset -- carrier phase to be selected
!! author   : Jianghui Geng
!! created  : Dec 12 2011
!
SUBROUTINE fcb_scr_identify(sigm,maxcph,iset,NM,infs)
!!
!*
USE par
USE info
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: maxcph,iset(1:*)
REAL(RL) :: sigm,infs(1:*)
TYPE(INFM) :: NM

  !*
  ! The local variables
  !!-------------------------
  LOGICAL(LG) :: lfnd
  INTEGER(IT) :: i,j,imax,nb,nobs,sset(NM.nobs),num
  REAL(RL) :: newsig,ssig

  !*
  ! The function called
  !!--------------------------
  LOGICAL(LG) :: chos,chitst
  REAL(RL) :: psig

  !*
  ! Start the exectuable code
  !!---------------------------


  !! loop over number of observations
  iset(1:NM.nobs/2)=0
  imax=min(NM.nobs/2,maxcph)

  IF (imax .GT. 5) THEN
   num=imax-5
   IF (imax .GT. 15) num=5
   IF (imax .GT. 20) num=4
   IF (imax .GT. 30) num=3
   IF (imax .GT. 60) num=2
  ELSE
   num=imax
  END IF

  DO i=1,num
  !DO i=1,imax
    ssig=1.d10
    DO WHILE(chos(i,NM.nobs/2,iset))
      newsig=psig(i,iset,0,NM,infs)
      IF (newsig .LT. ssig) THEN
        ssig=newsig
        sset(1:NM.nobs/2)=iset(1:NM.nobs/2)
      END IF
    END DO
    IF (ssig .LT. sigm) THEN
      iset(1:NM.nobs/2)=sset(1:NM.nobs/2)
      EXIT
    ELSE
      iset(1:NM.nobs/2)=0
    END IF
  END DO

  !! identify fake biases
  nb=COUNT(iset(1:NM.nobs/2).NE.0)
  IF (nb .GT. 0) THEN
    lfnd=.TRUE.
    DO WHILE(lfnd .AND. nb.GT.0)
      lfnd=.FALSE.
      DO i=1,nb
        newsig=psig(nb,iset,i,NM,infs)
        IF (chitst(1,NM.nobs-NM.np,ssig,newsig,0.9d0)) THEN
          lfnd=.TRUE.
          ssig=newsig
          DO j=i,nb-1
            iset(j)=iset(j+1)
          END DO
          iset(nb) =0
          nb       =nb-1
          EXIT
        END IF
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE

!
!! purpose  : least squares adjustment to generate new sigma
!! parameter:
!!    input : nb -- number of cycle slips
!!            iset -- problematic observations identified
!!            NM,infs -- information matrix
!! author   : Jianghui Geng
!! created  : Dec 12 2011
!

REAL(RL) FUNCTION psig(nb,iset,idel,NM,infs)
!!
!*
USE par
USE info
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: nb,iset(1:*),idel
REAL(RL) :: infs(1:*)
TYPE(INFM) NM

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j,k
  REAL(RL) :: dump
  REAL(RL), POINTER :: c(:,:),w(:)

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! compute LtPL
  psig=0.d0
  DO k=1,NM.nobs
    psig=psig+NM.resi(k)**2
  END DO

  !! compute new sigma0
  ALLOCATE(c(nb,nb))
  ALLOCATE(w(nb))
  c=0.d0
  w=0.d0
  DO i=1,nb
    DO j=1,nb
      IF (j.GE.i) THEN
        DO k=1,NM.nobs
          c(i,j)=c(i,j)+infs(NM.iptx(NM.imtx+1+iset(i))+NM.imtx+k)* &
                        infs(NM.iptx(NM.imtx+1+iset(j))+NM.imtx+k)
        END DO
      ELSE
        c(i,j)=c(j,i)
      END IF
    END DO

    !! right hand side
    DO k=1,NM.nobs
      w(i)=w(i)+infs(NM.iptx(NM.imtx+1+iset(i))+NM.imtx+k)*NM.resi(k)
    END DO
  END DO

  !! deleted elements
  IF (idel .NE. 0) THEN
    c(idel,idel)=c(idel,idel)+1.d8
  END IF

  !! normal matrix inversion
  CALL matinv(c,nb,nb,dump)
  IF (dump .EQ. 0.d0) THEN
    WRITE(ERROR_UNIT,'(a)') '***ERROR(ckd_scr_identify): matrix inversion'
    CALL exit(1)
  END IF
  DO i=1,nb
    dump=0.d0
    DO j=1,nb
      dump=dump+c(i,j)*w(j)
    END DO
    psig=psig-dump*w(i)
  END DO
  psig=dsqrt(psig/NM.nobs)
  DEALLOCATE(c)
  DEALLOCATE(w)

  RETURN

END FUNCTION
