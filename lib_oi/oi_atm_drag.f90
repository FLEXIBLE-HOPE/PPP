!*
SUBROUTINE oi_atm_drag(model,lpart,mjdutc,PAN,mass,npar,pname,xics,rot,esat,esun,acc,cmat)
!!
!*
USE const
USE tables
USE satellite
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
TYPE(SATEPAN) :: PAN
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar
CHARACTER(LEN=*) :: model,pname(1:*)
REAL(RL) :: mjdutc,mass,xics(1:*),acc(1:*)
REAL(RL) :: cmat(1:*),rot(3,3),esat(1:*),esun(1:*)

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i,j,k

  REAL(RL) :: gsat(3),gsun(3)
  REAL(RL) :: kp,f107(2),den,avv(3)
  REAL(RL) :: v_unit(3),vel,f(3),cosa,vwind(3)

  INTEGER(IT), PARAMETER :: MAXPARLOC=1
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(LEN_ORBPAR) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname / 'DRAG_c     '/

  REAL(RL) :: fpt(2),fbar(2),akp(4)
  REAL(RL) :: tz,tinf,tp120,ro,d(6),wmm,hl
  CHARACTER(LEN_STRING) :: line
  LOGICAL(LG) :: lfirst,lexist
  DATA lfirst /.TRUE./
  SAVE lfirst

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!-----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    IF (INDEX(model,'DTM13') .NE. 0) THEN
      line=f_tableFileName('dtm13')

      lexist=.TRUE.
      INQUIRE(FILE=line,EXIST=lexist)
      IF (lexist .EQ. .FALSE.) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_atm_drag): '//TRIM(line)//' is not exist.'
        CALL exit(1)
      END IF
      i=get_valid_unit(10)
      CALL P_ReadDTM12(i,line)
    END IF
  END IF

  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO

  DO i=1, MAXPARLOC
    param(i)=1.d0
    IF (ltog(i) .NE. 0) param(i)=param(i)+xics(ltog(i))
  END DO

  CALL xyzblh(esat(1:3)*1.d3,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,gsat)
  CALL xyzblh(esun(1:3)*1.d3,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,gsun)

  IF (INDEX(model,'DTM94') .NE. 0) THEN
    cosa=(6+(3.d0-6.d0)*(DABS(gsat(1))*RAD2DEG-0.d0)/90.d0)/24.d0
    CALL table_linear_interpolate('geomag_kp',.FALSE.,mjdutc-cosa,kp)
    CALL table_linear_interpolate('solar_flux',.FALSE.,mjdutc-1.d0,f107)
    ! kg/m3
    CALL dtm94(mjdutc,gsat,gsun,f107(1),f107(2),kp,den)
  ELSE IF (INDEX(model,'DTM13') .NE. 0) THEN
    fpt=0.d0
    fbar=0.d0
    akp=0.d0
    CALL table_linear_interpolate('geomag_kp',.FALSE.,mjdutc-0.125d0,kp)
    akp(1)=kp
    DO i=0, 8
      CALL table_linear_interpolate('geomag_kp',.FALSE.,mjdutc-i*-0.125d0,kp)
      akp(3)=akp(3)+kp
    END DO
    akp(3)=akp(3)/9.d0
    CALL table_linear_interpolate('solar_flux',.FALSE.,mjdutc-1.d0,f107)
    fpt(1)=f107(1)
    CALL table_linear_interpolate('solar_flux',.FALSE.,mjdutc,f107)
    fbar(1)=f107(1)

    hl=gsat(2)*rad2deg
    IF (hl .LT. 0.d0) hl=hl+360d0
    CALL mjd2doy(INT(mjdutc),i,k)
    wmm=(mjdutc-INT(mjdutc))*86400d0
    hl=DMOD(wmm/3600d0+hl/15d0,24d0)*15.0*DEG2RAD

    CALL dtm2012(DBLE(k),fpt,fbar,akp,gsat(3)/1.d3,hl,gsat(1),gsat(2),tz,tinf,tp120,ro,d,wmm)
    den=ro*1.d3
  ELSE
    den = 2.125d-11
  END IF

  ! Wind velocity = omega x r_satellite
  ! assuming that the atmosphere co-ratates with the Earth
  vwind(1)=-esat(2)*E_ROTATE
  vwind(2)= esat(1)*E_ROTATE
  vwind(3)= 0.d0
  !IF (INDEX(model,'HWM07') .NE. 0 ) THEN
  !  CALL wind_velocity(mjdutc,gsat,avv)
  !  vwind=vwind+avv
  !END IF

  ! Velocity in earth fixed system       
  DO i=1, 3
    avv(i)=esat(i+3)-vwind(i)
  END DO

  ! Rotate it to inertial system
  CALL matmpy(rot,avv,v_unit, 3, 3, 1)
  CALL unit_vector(3,v_unit,v_unit,vel)

  f=0.d0
  DO i=1, PAN.npan
    cosa=dot(3,v_unit,PAN.normj(1,i))
    IF (cosa .LE. 0.d0) CYCLE
    DO j=1,3
      f(j)=f(j)+PAN.area(i)*cosa*v_unit(j)*PAN.drag(i)
    END DO
  END DO

  DO i=1,3
    f(i)=-0.5/mass*den*f(i)*vel*vel*1.d3
    acc(i)=acc(i)+param(1)*f(i)
  END DO

  IF (.NOT. lpart) RETURN

  DO i=1, MAXPARLOC
    IF (ltog(i) .EQ. 0) CYCLE
    j=(ltog(i)-6-1)*3
    DO k=1, 3
      SELECT CASE(i)
        CASE(1)
          cmat(j+k)=f(k)
        CASE DEFAULT
      END SELECT
    END DO
  END DO

  RETURN

END SUBROUTINE
