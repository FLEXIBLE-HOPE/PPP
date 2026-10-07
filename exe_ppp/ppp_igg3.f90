!*
!! 添加IGG3调权/方差
!! xsy-2022-07-27
!
SUBROUTINE ppp_igg3(nobs,CKF,OB,NM)
!!
!*
USE info
USE const
USE ckdctrl
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(INFM) :: NM
INTEGER(IT) :: nobs
   
   !*
   ! The local variables
   !!-----------------------
   INTEGER(IT) :: i,isat,ifreq,lctype
   REAL(RL) :: obsres,wfac

   i=1
   DO WHILE(i .LE. nobs)
      isat =NM.ipob(i,1)
      ifreq=NM.ipob(i,4)
      lctype=NM.ipob(i,2) !xsy: 1:phase 2:range 3:ion const 4:ztd const
      !2023-07-22 by xsy: IGG3 change
      IF (lctype.EQ.2) THEN !Range
         obsres=NM.resi(i)/NM.weig(i)
         wfac=CKF.Igg3_k0/abs(obsres)*((CKF.Igg3_k1-abs(obsres))/(CKF.Igg3_k1-CKF.Igg3_k0))**2
         wfac=1.d0/wfac

         IF (abs(obsres) .GT. CKF.Igg3_k0 .AND. abs(obsres) .LE. CKF.Igg3_k1) THEN
            OB.var(isat,MAXFREQ*(lctype-1)+ifreq)=OB.var(isat,MAXFREQ*(lctype-1)+ifreq)*wfac
         ELSE IF (abs(obsres) .GT. CKF.Igg3_k1) THEN
            ! OB.var(isat,MAXFREQ*(lctype-1)+ifreq)=OB.var(isat,MAXFREQ*(lctype-1)+ifreq)*10000000000.d0
            ! OB.var(isat,ifreq)=OB.var(isat,ifreq)*10000000000.d0
            ! OB%omc(isat,MAXFREQ*(lctype-1)+ifreq)=0.d0
            ! OB%omc(isat,                   ifreq)=0.d0
            OB%omc(isat,:)=0.d0
            WRITE(OUTPUT_UNIT,'(A,1X,F6.2,1X,A,1X,F5.1)') ' ... [IGG3] CODE RESIDUAL FOR '//CKF%cprn(isat)//':',abs(obsres),'K1: ',CKF%Igg3_k1
         END IF
      ! ELSE IF (lctype.EQ.1) THEN !Phase
      !    IF (abs(obsres) .GT. 0.1d0) THEN
      !       OB.var(isat,MAXFREQ*(lctype-1)+ifreq)=OB.var(isat,MAXFREQ*(lctype-1)+ifreq)*10000000000.d0
      !    END IF
      END IF
      i=i+1
   END DO

   NM.nobs=0
   RETURN
  
END SUBROUTINE

