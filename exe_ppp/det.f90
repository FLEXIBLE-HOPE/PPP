!*
SUBROUTINE Ddet(A,N,D)
!*
!! PARAMETERS
!!      A: matric
!!      N: nxn
!!      D: return determination
!! PURPOSE
!!      CACULATE THE ADOP
!! CREATED BY SHENGYI XU
!!      2024-03-02
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: N
REAL(RL) :: A(N,N),D

    !*
    ! The local variables
    !!------------------------
    INTEGER(IT) :: i,j,k,nind
    REAL(RL) :: s,u(MAXPARSIT+MAXSAT*2),beta,gama
    
    !*
    ! Start the exectuable code
    !!--------------------------

    !! Householder transformation
    DO j=1,N
        !! whether necessary for this column
        IF (ALL(A(j+1:N,j).EQ.0.D0)) THEN
            CYCLE
        ELSE
            !! get elemental transformation vector
            CALL ppp_ele_hht(N-j+1,A(j:N,j),s,u,beta)
            A(j,j)=s

            !! transform the next columns
            DO k=j+1,N
                gama=0.d0
                DO i=j,N
                    IF (u(i-j+1).EQ.0.D0 .OR. A(i,k).EQ.0.D0) CYCLE
                    gama=gama+u(i-j+1)*A(i,k)
                END DO
                IF (gama.EQ.0.D0) CYCLE
                gama=beta*gama
                DO i=j,N
                    IF (u(i-j+1).EQ.0.D0) CYCLE
                    A(i,k)=A(i,k)+gama*u(i-j+1)
                END DO
            END DO
            !! clean current column
            DO i=j+1,N
                A(i,j)=0.D0
            END DO
        END IF
    END DO
    
    D=-1.D0
    nind=0
    DO i=1,N
        IF (A(i,i).EQ.0.D0) CYCLE
        D=D*A(i,i)
        nind=nind+1
    END DO
    WRITE(*,*)D,N,nind
    ! ADOP
    D=SQRT(D)**(1.d0/nind)

END SUBROUTINE

! subroutine determinant(A,N,d)
!     !----------------------------------------------------------------------
!     !--------输入-------
!     !A——方阵A(N,N)
!     !N——A的维度
!     !--------输出--------
!     !d——行列式
!     !-----------------------------------------------------------------------
!     implicit real*8(a-h,o-z), integer*4 (i-n)
!     real*8 A(N,N),A_New(N,N)
!     real*8 Paixu(N,1)

!     xuhao=0.d0
!     A_New=A
!     d=1.d0
!     do i=1,N
!         Paixu(i,1)=0.d0
!     enddo
!     do j=1,N
!         do i=1,N
!         if((A_New(i,j)/=0.d0).and.(Paixu(i,1)==0.d0))then
!             xuhao=xuhao+1.d0
!             Paixu(i,1)=xuhao
!             d=d*A_New(i,j)
!             do m=1,N
!                 if(Paixu(m,1)==0.d0)then
!                     A_New(m,:)=A_New(m,:)-A_New(i,:)*A_New(m,j)/A_New(i,j)
!                 endif
!             enddo
!             exit
!         endif
!         enddo
!     enddo
    
!     !write(*,*) Paixu
!     !write(*,*)
    
!     !逆序数
!     Nixu_num=0.d0
!     do i=1,N-1
!         do j=i+1,N
!             if(Paixu(i,1)>Paixu(j,1))then
!                 Nixu_num=Nixu_num+1.d0
!             endif
!         enddo
!     enddo
!     !write(*,*) Nixu_num
!     !write(*,*)
    
!     !d——行列式
!     d=(-1.d0)**Nixu_num*d
    
!  end