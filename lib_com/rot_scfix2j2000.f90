!*
SUBROUTINE rot_scfix2j2000(mjd,sod,cprn,csvn,blk,xsat,xsun,xscf,yscf,zscf)
!!
!! purpose   : rotation matrix from spacecraft(GNSS)-fixed system to J2000
!!             The yaw-error correction should be implemented later.
!!
!! parameters:
!!        xsat -- satellite coordinates J2000
!!        xsun -- sun coordinates in J2000
!!        xscf,yscf,zscf -- unit vector of sc-fixed x-axes in J2000
!!                 They define the transformation from x(scf) to j2000 as following
!!                 x(j2000) =  [ xscf,yscf,zscf ] x(scf)
!! 
!! created by: Ge Maorong
!!
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
INTEGER(IT) :: mjd
CHARACTER(LEN=*) :: cprn,blk,csvn
REAL(RL) :: sod,xsat(1:*),xsun(1:*)
REAL(RL) :: xscf(1:*),yscf(1:*),zscf(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i
  REAL(RL) :: beta,u,yaw
  REAL(RL) :: alpha,u_sc2sun(3)
  LOGICAL(LG) :: lfix

  !*
  ! The functions used
  !!----------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuablde code
  !!--------------------------------

  lfix=.FALSE.

  !! sc-fixed z-axis, from sc to earth
  DO i=1 ,3
    zscf(i)=-xsat(i)
  END DO
  CALL unit_vector(3,zscf,zscf,yaw)

  !! unint vector from sc to sun
  DO i=1,3
    u_sc2sun(i)=xsun(i)-xsat(i)
  END DO
  CALL unit_vector(3,u_sc2sun,u_sc2sun,yaw)

  SELECT CASE(TRIM(blk))
    CASE('BEIDOU-2G','BEIDOU-3G-SECM','BEIDOU-3G-CAST')
      lfix=.TRUE.
    CASE('BEIDOU-2I','BEIDOU-2M','QZSS')

      !CALL betau(xsat,xsun,beta,u)
      !yaw=DATAN2(-DTAN(beta),DSIN(u))

      !! The epoch has been added by eclipse, so we remove the code
      IF (TRIM(csvn).EQ.'017' .OR. TRIM(csvn).EQ.'015' .OR. TRIM(csvn).EQ.'005') THEN
        lfix=.FALSE.
      ELSE
        CALL bds_s2f_epoch(mjd,sod,cprn,lfix)
      END IF

      ! if the start day is alreday in yaw-fixed mode, this code has some problems
      !IF (lfix .EQ. .FALSE.) THEN
      !  IF (DABS(beta).LE.4.d0*DEG2RAD .AND. DABS(yaw).LE.5.d0*DEG2RAD) THEN
      !    CALL bds_s2f_out(.TRUE.,mjd,sod,cprn,beta*RAD2DEG,u*RAD2DEG)
      !    lfix=.TRUE.
      !  END IF
      !END IF

      !IF (lfix .EQ. .TRUE.) THEN
      !  IF (DABS(beta).GE.4.d0*DEG2RAD .AND. DABS(yaw).GE.5.d0*DEG2RAD) THEN
      !    CALL bds_s2f_out(.FALSE.,mjd,sod,cprn,beta*RAD2DEG,u*RAD2DEG)
      !    lfix=.FALSE.
      !  END IF
      !END IF
  END SELECT

  IF (lfix .EQ. .TRUE.) THEN

     DO i=1,3
       xscf(i)=xsat(i+3)
     END DO
     CALL unit_vector(3,xscf,xscf,yaw)

     !! the reverse orbit normal direction
     CALL cross(zscf,xscf,yscf)
     CALL unit_vector(3,yscf,yscf,yaw)

     CALL cross(yscf,zscf,xscf)
     CALL unit_vector(3,xscf,xscf,yaw)

  ELSE

    !! angle between u_sc2sun and zscf, if they near parallel yscf as cross product
    !! of them can not be defined. Values from the last epoch should be used
    alpha=dot(3,zscf,u_sc2sun)
    IF (DABS(alpha).GT.1.d0) alpha=NINT(alpha)
    alpha=DACOS(alpha)*RAD2DEG
    IF (alpha.GE.1.d-6 .AND. alpha.LE.180.d0-1.d-6) THEN

      !! the sc-fixed y-axis.
      CALL cross(zscf,u_sc2sun,yscf)
      CALL unit_vector(3,yscf,yscf,yaw)

      !! the sc-fixed z-axis
      CALL cross(yscf,zscf,xscf)
      CALL unit_vector(3,xscf,xscf,yaw)

    ELSE
      WRITE(OUTPUT_UNIT,'(A)') '###WARNING(rot_scfix2j200): scx/y no definition, previous epoch used '
    END IF

    SELECT CASE(TRIM(blk))
      CASE('BLOCK IIA','BLOCK IIR-A','BLOCK IIR-B','BLOCK IIR-M','BLOCK IIF','BLOCK IIIA', &
!      CASE('BLOCK IIA','BLOCK IIR-A','BLOCK IIR-B','BLOCK IIR-M','BLOCK IIF', &
           'GLONASS-M', 'GLONASS-K1', &
           'GALILEO-1','GALILEO-2', &
           'BEIDOU-3I-SECM','BEIDOU-3I-CAST','BEIDOU-3M-SECM','BEIDOU-3M-CAST', &
           'BEIDOU-3IS-SECM','BEIDOU-3IS-CAST','BEIDOU-3MS-SECM','BEIDOU-3MS-CAST', &
           'BEIDOU-3SI-SECM','BEIDOU-3SI-CAST','BEIDOU-3SM-SECM','BEIDOU-3SM-CAST')
        CALL yawatt(mjd+sod/86400.d0,cprn,csvn,blk,xsat,xsun,xscf,yscf,zscf)
      CASE('BEIDOU-2I','BEIDOU-2M')
        IF (TRIM(csvn).EQ.'017' .OR. TRIM(csvn).EQ.'015' .OR. TRIM(csvn).EQ.'005') THEN
          CALL yawatt(mjd+sod/86400.d0,cprn,csvn,blk,xsat,xsun,xscf,yscf,zscf)
        END IF
    END SELECT

  END IF

  RETURN

END SUBROUTINE
