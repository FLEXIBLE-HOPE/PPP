!
!! purpose  : integer resolution of cycle slips (partial search is allowed)
!! parameter:
!!    input : nsip -- # of cycle slips
!!            esti -- estimates
!!            nmat -- variance-covariance
!!    output: fixi -- integer cycle slips
!! author   : Geng J
!! created  : Jan 1 2011
!
SUBROUTINE rep_integer(npar,nost,nsip,esti,ltpl,nmat,fixi,chisq,ratio)
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: npar,nsip,nost
REAL(RL) :: ltpl,ratio,chisq
REAL(RL) :: esti(npar+nsip),nmat(npar+nsip,npar+nsip),fixi(nsip)
INTEGER(IT), PARAMETER :: max_del=2,min_sav=5
REAL(RL), PARAMETER :: crit=3.d0

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i,j,k,l,m,nfix
  REAL(RL) :: temp,disall(2)
  INTEGER(IT), POINTER :: idel(:),sdel(:)
  REAL(RL), POINTER :: q22(:),bias(:),sbias(:)

  !*
  ! The function used
  !!--------------------------
  LOGICAL(LG) :: selec
  INTEGER(IT) :: pointer_int

  !*
  ! Start the exectuable code
  !!--------------------------

  !! initialization
  ! number of good carrier-phase measurements
  nfix=nost-nsip
  ALLOCATE(q22(nsip*(nsip+1)/2))
  ALLOCATE(bias(nsip))
  k=1
  DO j=npar+1,npar+nsip
    DO i=j,npar+nsip
      q22(k)=nmat(i,j)
      k=k+1
    END DO
    bias(j-npar)=esti(j)
  END DO

  !! try firstly integer search for all cycle slips
  CALL ambslv(nsip,q22,bias,disall)
  IF (npar .EQ. 0) THEN
    chisq=(disall(1)+ltpl)/nost/ltpl*(nost-nsip)
  ELSE
    chisq=(disall(1)+ltpl)/(nost*2-npar)/ltpl*(nost*2-npar-nsip)
  END IF
  ratio=disall(2)/disall(1)

  !! partial search
  IF (ratio .GE. crit) THEN
    DO i=1, nsip
      fixi(i)=NINT(bias(i))
    END DO
  ELSE
    ALLOCATE(idel(max_del))
    ALLOCATE(sdel(max_del))
    ALLOCATE(sbias(nsip))
    DO i=1, max_del
      IF (nfix+nsip-i.LT.min_sav .OR. i.GE.nsip) EXIT ! at least ...
      ratio=0.d0
      idel(1:max_del)=0
      sdel(1:max_del)=0
      DO WHILE(selec(nsip,i,idel))
        k=0
        m=0
        DO j=npar+1,npar+nsip
          IF (pointer_int(i,idel,j-npar) .NE. 0) CYCLE
          m=m+1
          bias(m)=esti(j)
          DO l=j,npar+nsip
            IF (pointer_int(i,idel,l-npar) .NE. 0) CYCLE
            k=k+1
            q22(k)=nmat(l,j)
          END DO
        END DO
        CALL ambslv(nsip-i,q22,bias,disall)
        temp=disall(2)/disall(1)

        !! save most possible solutions in terms of ratio values
        IF (temp .GT. ratio) THEN
          ratio=temp
          IF (npar .EQ. 0) THEN
            chisq=(disall(1)+ltpl)/nost/ltpl*(nost-nsip+i)
          ELSE
            chisq=(disall(1)+ltpl)/(nost*2-npar)/ltpl*(nost*2-npar-nsip+i)
          END IF
          sdel(1:i)=idel(1:i)
          sbias(1:nsip-i)=bias(1:nsip-i)
        END IF
      END DO
      IF (ratio .GE. crit) THEN
        k=0
        DO j=1,nsip
          fixi(j)  =0.5d0
          IF (pointer_int(i,sdel,j) .EQ. 0) THEN
            k=k+1
            fixi(j)=NINT(sbias(k))
          END IF
        END DO
        EXIT
      END IF
    END DO
    DEALLOCATE(idel)
    DEALLOCATE(sdel)
    DEALLOCATE(sbias)
  END IF

  !! clean memory
  DEALLOCATE(q22)
  DEALLOCATE(bias)

  RETURN

END SUBROUTINE


!! Try different ambiguity combination to be removed
LOGICAL(LG) FUNCTION selec(nsip,ndel,idel)
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: nsip,ndel,idel(1:*)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i, ic

  !*
  ! Start the exectuable code
  !!---------------------------

  !! initialize ambiguity array
  IF (idel(1) .EQ. 0) THEN
    DO i=1,ndel
      idel(i)=i
    END DO
    selec=.TRUE.
    RETURN
  END IF

  !! change deleted ambiguities
  ic=ndel
  DO WHILE(ic .GT. 0)
    idel(ic)=idel(ic)+1
    IF (idel(ic) .GT. nsip) THEN
      ic=ic-1
      CYCLE
    END IF
    i=ic+1
    DO WHILE(i .LE. ndel)
      idel(i)=idel(i-1)+1
      IF (idel(i) .GT. nsip) THEN
        ic=ic-1
        EXIT
      END IF
      i=i+1
    END DO
    IF (i .GT. ndel) EXIT
  END DO
  selec=.TRUE.
  IF (ic .LE. 0) selec=.FALSE.

  RETURN

END FUNCTION
