!
! purpose  : generally undifferenced partial ambiguity fixing for EWL and WL
!! parameter:
!!    input : namb-- number of ambiguty to be fixed
!!            bias-- float ambiguity
!!            qxx -- Qxx
!!    output: 
!!            bias -- fixed ambiguity
!!            ratio -- resulted test statistics
!! author   : shengyi xu
!! created  : 2022-11-2
!
SUBROUTINE ppp_abfix_wlpar(namb,bias,qxx,maxdel,minrto,ratio)
!!
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: namb,maxdel
REAL(RL) :: qxx(1:*),bias(1:*),minrto,ratio

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i,j,k,m,l
  REAL(RL) :: disall(2),tratio
  INTEGER(IT) :: idel(1:namb),sdel(1:namb)
  REAL(RL) :: q22(namb*(namb+1)/2),bbi(1:namb),cobias(1:namb),note(1:namb),svbias(1:namb),conote(1:namb)

  !*
  ! The function called
  !!----------------------------
  LOGICAL(LG) :: chos
  INTEGER(IT) :: pointer_int

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! try search for all integer ambiguities
  svbias(1:namb)=bias(1:namb)
  CALL ambslv(namb,qxx,bias,disall)
  ratio=disall(2)/disall(1)

  !! partial ambiguity resolution
  IF (ratio .LE. minrto) THEN

    !! partial ambiguity fixing
    IF (maxdel .GT. 0) THEN

      !! remove possibly biased ambiguities
      DO i=1,maxdel
        IF (namb-i.LE.0) EXIT
        ratio=0.d0
        idel(1:maxdel)=0
        sdel(1:maxdel)=0
        ! selected according the variance and co-variance
        DO WHILE(chos(i,namb,idel))
          k=0
          m=0
          note=0
          DO j=1,namb
            ! removed
            IF (pointer_int(maxdel,idel,j) .NE. 0) CYCLE
            m=m+1
            bbi(m)=svbias(j)
            note(m)=j
            DO l=j,namb
              ! removed
              IF (pointer_int(maxdel,idel,l) .NE. 0) CYCLE
              k=k+1
              q22(k)=qxx((j-1)*(2*namb-j+2)/2+l-j+1)
            END DO
          END DO
          CALL ambslv(m,q22,bbi,disall)
          tratio=disall(2)/disall(1)

          !! save most possible solutions in terms of ratio values
          IF (tratio .GT. ratio) THEN
            ratio=tratio
            sdel(1:i)=idel(1:i)
            bias(1:namb-i)=bbi(1:namb-i)
            conote(1:namb-i)=note(1:namb-i)
          END IF
        END DO
        ! i the number of removed
        IF (ratio .GT. minrto) THEN
          cobias(1:namb-i)=bias(1:namb-i)
          DO j=1,i
            !ambset(sdel(j))=1.d0
            bias(sdel(j))=svbias(sdel(j))
          END DO
          DO j=1,namb-i
            bias(conote(j))=cobias(j)
          END DO

          RETURN
        END IF
      END DO
    END IF
  ELSE
    RETURN
  END IF

  RETURN

END SUBROUTINE
