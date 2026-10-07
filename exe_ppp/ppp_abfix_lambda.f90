!
!! purpose  : generally undifferenced ambiguity fixing
!! parameter:
!!    input : QM,invx -- inverted normal matrix
!!    output: AB -- ambiguity struct
!!            bias -- ambiguity increment
!!            ratio -- resulted test statistics
!! author   : Geng J
!! created  : Sep 6 2011
!
SUBROUTINE ppp_abfix_lambda(AB,QM,invx,maxdel,minsav,minrto,ratio)
!!
!*
USE info
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(AMBD) :: AB(1:*)
TYPE(INVM) :: QM
INTEGER(IT) :: maxdel,minsav
REAL(RL) :: invx(1:*),minrto,ratio

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i,j,k,m,l,imax,nxyz,itim,nnl
  REAL(RL) :: disall(2),tratio,adop,det
  INTEGER(IT), POINTER :: idel(:),sdel(:)
  REAL(RL), POINTER :: bias(:),q22(:),bbi(:),qamb(:,:)

  !*
  ! The function called
  !!----------------------------
  LOGICAL(LG) :: chos
  INTEGER(IT) :: pointer_int

  !*
  ! Start the exectuable code
  !!-----------------------------

  !! new nxyz (may include redundant ambiguities)
  nxyz=QM%ntot-QM%ndam                              !xsy:  QM.ndam:待固定的窄巷模糊度，在Qxx(invx)最后

  !! memory allocation
  ALLOCATE(bias(QM%ntot))
  ALLOCATE(q22(QM%ndam*(QM%ndam+1)/2))              !xsy: 模糊度Qxx(方差协方差，上三角)

  !! copy cofactors to an array
  k=0
  nnl=0
  adop=1.d0
  DO j=nxyz+1,nxyz+QM%ndam
    AB(j-QM%nxyz)%abst=1  ! presumed fixed: 假设固定,要固定的模糊度朝后放
    bias(j-nxyz)=AB(j-QM%nxyz)%abfr
    DO i=j,nxyz+QM%ndam
      k=k+1
      q22(k)=invx(QM%idq(j)+i)
      IF (j.EQ.i) THEN
        adop=adop*invx(QM%idq(j)+i)
        nnl=nnl+1
      END IF
    END DO
  END DO
  QM%adop=SQRT(adop)**(1.d0/nnl)

  ! ALLOCATE(qamb(nnl,nnl))
  ! DO j=nxyz+1,nxyz+QM%ndam
  !   DO i=j,nxyz+QM%ndam
  !     qamb(i,j)=invx(QM%idq(j)+i)
  !     qamb(j,i)=qamb(i,j)
  !   END DO
  ! END DO
  ! WRITE(100,*)'QAMB ->'
  ! DO j=1,nnl
  !   WRITE(100,'(<nnl>(F14.10,2X))')(qamb(i,j),i=1,nnl)
  ! END DO
  ! CALL Ddet(qamb,nnl,det)
  ! QM%adop=det

  !! try search for all integer ambiguities
  QM%ncad=QM%ndam
  CALL ambslv(QM%ncad,q22,bias,disall)
  ratio=disall(2)/disall(1)

  !! partial ambiguity resolution
  IF (ratio .LE. minrto) THEN

    !! the last removable ambiguity
    QM.ncad=0
    imax=QM.ndam

    !! partial ambiguity fixing
    IF (maxdel .GT. 0) THEN
      ALLOCATE(idel(maxdel))
      ALLOCATE(sdel(maxdel))
      ALLOCATE(bbi (QM.ndam))
      itim = 0
      !! remove possibly biased ambiguities
      DO i=1,maxdel
        IF (QM.ndam-i.EQ.0 .OR. QM.nfix+QM.ndam-i.LT.minsav) EXIT
        ratio=0.d0
        idel(1:maxdel)=0
        sdel(1:maxdel)=0
        ! selected according the variance and co-variance
        DO WHILE(chos(i,imax,idel))
          itim = itim+1
          k=0
          m=0
          DO j=nxyz+1,nxyz+QM.ndam
            ! removed
            IF (pointer_int(maxdel,idel,j-nxyz) .NE. 0) CYCLE
            m=m+1
            bbi(m)=AB(j-QM.nxyz).abfr
            DO l=j,nxyz+QM.ndam
              ! removed
              IF (pointer_int(maxdel,idel,l-nxyz) .NE. 0) CYCLE
              k=k+1
              q22(k)=invx(QM.idq(j)+l)
            END DO
          END DO
          CALL ambslv(QM.ndam-i,q22,bbi,disall)
          tratio=disall(2)/disall(1)
          ! WRITE(*,*)i,maxdel,QM.ndam,itim,tratio
          !! save most possible solutions in terms of ratio values
          IF (tratio .GT. ratio) THEN
            ratio=tratio
            sdel(1:i)=idel(1:i)
            bias(1:QM.ndam-i)=bbi(1:QM.ndam-i)
          END IF
        END DO
        ! i the number of removed
        IF (ratio .GT. minrto) THEN
          DO j=1,i
            AB(nxyz-QM.nxyz+sdel(j)).abst=0
          END DO
          QM.ncad=QM.ndam-i                         !xsy: 待固定模糊度
          EXIT
        END IF
      END DO
      DEALLOCATE(idel)
      DEALLOCATE(sdel)
      DEALLOCATE(bbi)
    END IF
  END IF
  !! collect results: adapt ambiguity estimates
  IF (QM.ncad .GT. 0) THEN
    k=0
    DO i=nxyz+1,nxyz+QM.ndam
      IF (AB(i-QM.nxyz).abst .EQ. 1) THEN
        k=k+1
        AB(i-QM.nxyz).abfx=NINT(bias(k))
      END IF
    END DO
  ELSE
    DO i=nxyz+1,nxyz+QM.ndam
      AB(i-QM.nxyz).abst=0 ! reset
    END DO
  END IF
  DEALLOCATE(bias)
  DEALLOCATE(q22)
  ! DEALLOCATE(qamb)

  RETURN

END SUBROUTINE
