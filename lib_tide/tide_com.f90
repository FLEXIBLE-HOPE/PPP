!*
SUBROUTINE tide_com(mjd,dx)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!------------------
REAL(RL) :: mjd  ! in UT1
REAL(RL) :: dx(1:*) ! in meter


   !*
   ! The local variables
   !!-------------------------
   REAL(RL) :: f,arg(11)
   INTEGER(IT) :: iy,id,i,k

   REAL(RL) :: coef(6,11)

   ! http://holt.oso.chalmers.se/loading/CMC/
   DATA coef &  ! FES2004
   /-1.2661d-03,-1.4298d-03,-1.3724d-03, 8.2077d-04, 1.1479d-03, 2.3005d-04, &
    -1.7763d-04,-5.7273d-04,-5.3350d-04,-3.1591d-04,-5.1370d-05, 2.8184d-04, &
    -3.2372d-04,-2.8986d-04,-2.7121d-04, 1.9849d-04, 2.6018d-04,-1.4302d-04, &
    -1.1814d-04,-1.5250d-04,-1.1223d-04,-1.0889d-05,-1.5751d-05, 1.2367d-04, &
    -1.1370d-03, 4.4839d-03,-1.8539d-03,-8.6426d-04,-9.1022d-04,-1.7823d-03, &
    -1.6802d-04, 2.9702d-03,-1.3985d-03,-2.2975d-04,-8.8858d-04,-6.4989d-04, &
    -3.6495d-04, 1.4941d-03,-6.1436d-04,-2.9129d-04,-2.9261d-04,-5.7461d-04, &
     3.0709d-05, 4.5472d-04,-2.7831d-04,-2.9313d-05,-2.1734d-04,-4.1637d-05, &
    -5.0643d-04,-7.3040d-05,-2.2065d-04, 4.1472d-04,-1.0212d-04, 8.2276d-05, &
    -2.7885d-04, 2.0596d-05, 4.6882d-05, 1.8399d-04,-7.4897d-06, 1.3209d-05, &
    -1.4899d-04, 2.6146d-06, 1.3687d-04, 3.5475d-05,-2.4093d-05, 3.1666d-07 /

   !*
   ! Start the exectuable code
   !!--------------------------


   !! S1-S2 atmospheric pressure loading
   f=2*PI*(mjd-INT(mjd))
   dx(1)= 2.1188d-4*DCOS(f)-7.6861d-4*DSIN(f)+1.4472d-4*DCOS(2*f)-1.7844d-4*DSIN(2*f)
   dx(2)=-7.2766d-4*DCOS(f)-2.3582d-4*DSIN(f)-3.2691d-4*DCOS(2*f)-1.5878d-4*DSIN(2*f)
   dx(3)=-1.2176d-5*DCOS(f)+3.2243d-5*DSIN(f)-9.6271d-5*DCOS(2*f)+1.6978d-5*DSIN(2*f)

   CALL mjd2doy(INT(mjd),iy,id)
   CALL ARG2(iy,id+mjd-INT(mjd),arg)

   DO i=1, 3
     DO k=1, 11
       dx(i)=dx(i)+coef((i-1)*2+1,k)*DCOS(arg(k))+coef(i*2,k)*DSIN(arg(k))
     END DO
   END DO

   RETURN

END SUBROUTINE
