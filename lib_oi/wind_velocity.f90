! !*
! SUBROUTINE windVelocity(mjd,gsc,vel)
! !!
! !*
! USE const                  
! IMPLICIT NONE

! !*
! ! The arguments
! !!-----------------------------
! REAL(RL)  :: mjd,gsc(1:*),vel(1:*)

  
!   !*
!   ! The local variables
!   !!-----------------------------
!   INTEGER(IT) :: i, year, doy, iyd
!   REAL(RL) :: sec,ap(7),LNOF2TRF(3,3)
!   REAL(RL) :: w(3)

!   !*
!   ! Start the exectable code
!   !!----------------

!   ap=0d0
!   CALL tableLinearInterpolate('geomap', .FALSE., mjd, ap)
  
!   CALL MjdToDoy(INT(mjd),year,doy)
!   sec = (mjd-INT(mjd))*86400d0
!   iyd = MOD(year,100)*1000+doy
  
!   w=0d0
!   CALL hwm07(iyd,REAL(sec),REAL(GSc(3)/1d3),REAL(GSc(1)*rad2deg),REAL(GSc(2)*rad2deg),0.0,0.0,0.0,REAL(ap(1:2)),w(1:2))
  
!   ! m/s -> km/s
!   DO i=1, 3
!     Vel(i) = REAL(w(i),RL)*1d-3
!   ENDDO
  
!   CALL Rot_enu2xyz(GSc(1),GSc(2),LNOF2TRF)
!   CALL MatMpy(LNOF2TRF,Vel,Vel,3,3,1)
  
!   !*
!   ! RETURN
!   !!------------
!   RETURN

! END SUBROUTINE windVelocity