!*
SUBROUTINE oi_thrust_acc(lpart,npar,pname,xics,nman,tman,acr,mjd,cmat,acc0)
!!
!*
USE const
IMPLICIT NONE


!*
! The arguments
!!-----------------------------
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar,nman
CHARACTER(LEN=*) :: pname(1:*)
REAL(RL) :: xics(1:*),acr(3,3),tman(2,MAXMAN)
REAL(RL) :: cmat(1:*),acc0(1:*),mjd

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i,j,k,iterm,icomp,iman
  REAL(RL) :: acc(3),term(6),fact,dt

  LOGICAL(LG) :: lin

  INTEGER(IT) :: MAXPARLOC
  PARAMETER(MAXPARLOC=18)
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(15) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname &
  /'TH_a1   ','TH_a2   ','TH_a3   ','TH_a4   ','TH_a5   ','TH_a6   ', &
   'TH_c1   ','TH_c2   ','TH_c3   ','TH_c4   ','TH_c5   ','TH_c6   ', &
   'TH_r1   ','TH_r2   ','TH_r3   ','TH_r4   ','TH_r5   ','TH_r6   '/

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the executable code
  !!----------------
  
  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO


  acc=0.d0
  DO i=1, MAXPARLOC
    param(i)=0.d0
    IF (ltog(i) .EQ. 0) CYCLE
    IF (lpart .EQ. .TRUE.) THEN
      j=(ltog(i)-6-1)*3
      DO k=1, 3
        cmat(j+k)=0.d0
      END DO
    END IF
  END DO
 

  lin=.FALSE.
  DO i=1, nman
    IF (mjd.LT.tman(1,i) .OR. mjd.GT.tman(2,i)) CYCLE
    lin=.TRUE.
    iman=i
    dt=(mjd-tman(1,i))*86400.d0
    EXIT
  END DO
  
  IF (lin .EQ. .FALSE.) RETURN
  
  term(1:6)=1.d-12
  DO i=1, MAXPARLOC
    param(i)=0.d0
    IF (ltog(i) .EQ. 0) CYCLE
    param(i)=param(i)+xics(ltog(i))

    icomp=(i-1)/6+1
    iterm=i-(icomp-1)*6
    IF (iterm .NE. iman) CYCLE
    IF (lpart .EQ. .TRUE.) THEN
      j=(ltog(i)-6-1)*3 
      DO k=1, 3
        cmat(j+k)=term(iterm)*acr(k,icomp)
      END DO
    END IF
    DO k=1,3 
      acc(k)=acc(k)+term(iterm)*param(i)*acr(k,icomp)
    END DO
  END DO

  DO i=1, 3
    acc0(i)=acc0(i)+acc(i)
  END DO

  RETURN
      
END SUBROUTINE
