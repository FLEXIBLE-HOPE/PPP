!* 
SUBROUTINE get_freq(SAT,type,nfreq,cfreq)
!!
!! TO DETERMINATE THE FREQUENCIES USED
!!
!*
USE const
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
TYPE(SATE) :: SAT
CHARACTER(LEN=*) :: type
CHARACTER(LEN_FREQ) :: cfreq(MAXFREQ)
INTEGER(IT) :: nfreq
 
  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  INTEGER(IT) :: i=0,j=0
  REAL(RL) :: temp=0.d0
  CHARACTER(LEN_FREQ) :: tfreq

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  IF (nfreq .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(get_freq): The frequency used for '//SAT.cprn//' is not selected.'
    CALL exit(1)
  END IF

  SELECT CASE(TRIM(type))
    ! IONOSPHERE-FREE COMBINATION
    CASE('IF')
      !IF (nfreq .NE. 2) THEN
      !  WRITE(ERROR_UNIT,'(A)') '***ERROR(get_freq): The number of frequency used for '//SAT.cprn
      !  WRITE(ERROR_UNIT,'(A)') '                    is greater than 2 for ionosphere-free observation.'
      !  CALL exit(1)
      !END IF

      DO i=1, nfreq
        CALL getFreq(SAT.cprn(1:1),SAT.type,cfreq(i),SAT.ifreq,SAT.freq(i))
      END DO

      IF (SAT.freq(1) .LT. SAT.freq(2)) THEN
        temp=SAT.freq(2)
        SAT.freq(2)=SAT.freq(1)
        SAT.freq(1)=temp

        tfreq=cfreq(2)
        cfreq(2)=cfreq(1)
        cfreq(1)=tfreq
      END IF

      SAT.g=SAT.freq(1)/SAT.freq(2)
      SAT.g2=SAT.g*SAT.g
      SAT.lamdw=VEL_LIGHT/(SAT.freq(1)-SAT.freq(2))
      SAT.lamdn=VEL_LIGHT/(SAT.freq(1)+SAT.freq(2))
      SAT.fac(1)=SAT.g2/(SAT.g2-1.D0)
      SAT.fac(2)=1.D0/(SAT.g2-1.D0)

      DO i=1, nfreq
        SAT.lamda(i)=VEL_LIGHT/SAT.freq(i)
      END DO

    CASE('RAW','GRAPHIC')

      DO i=1, nfreq
        CALL getFreq(SAT.cprn(1:1),SAT.type,cfreq(i),SAT.ifreq,SAT.freq(i))
      END DO

      DO i=1, nfreq
        SAT.lamda(i)=VEL_LIGHT/SAT.freq(i)
      END DO

      IF (nfreq .GE. 2) THEN
        IF (SAT.freq(1) .LT. SAT.freq(2)) THEN
          temp=SAT.freq(2)
          SAT.freq(2)=SAT.freq(1)
          SAT.freq(1)=temp

          tfreq=cfreq(2)
          cfreq(2)=cfreq(1)
          cfreq(1)=tfreq
        END IF

        SAT.g=SAT.freq(1)/SAT.freq(2)
        SAT.g2=SAT.g*SAT.g
        SAT.lamdw=VEL_LIGHT/(SAT.freq(1)-SAT.freq(2))
        SAT.lamdn=VEL_LIGHT/(SAT.freq(1)+SAT.freq(2))
        SAT.fac(1)=SAT.g2/(SAT.g2-1.d0)
        SAT.fac(2)=1.d0/(SAT.g2-1.d0)
      END IF

    ! GEOMETRY-FREE AND IONOSPHERE-FREE COMBINATION
    CASE('GFIF')
      IF (nfreq .NE. 3) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(get_freq): The number of frequency used for '//SAT.cprn
        WRITE(ERROR_UNIT,'(A)') '                    do not equal 3 for geometry-free and ionosphere-free observation.'
        CALL exit(1)
      END IF

      DO i=1, nfreq
        CALL getFreq(SAT.cprn(1:1),SAT.type,cfreq(i),SAT.ifreq,SAT.freq(i))
      END DO

      ! do not array the frequency, but use the input as the default
!      DO i=1, 3
!        DO j=1, 3-i
!          IF (SAT.freq(j) .GT. SAT.freq(j+1)) THEN
!            temp=SAT.freq(j)
!            SAT.freq(j)=SAT.freq(j+1)
!            SAT.freq(j+1)=temp
!
!            tfreq=cfreq(j)
!            cfreq(j)=cfreq(j+1)
!            cfreq(j+1)=tfreq
!          END IF
!        END DO
!      END DO

      SAT.g=SAT.freq(1)/SAT.freq(2)
      SAT.g2=SAT.g*SAT.g
      SAT.lamdw=VEL_LIGHT/(SAT.freq(1)-SAT.freq(2))
      SAT.lamdn=VEL_LIGHT/(SAT.freq(1)+SAT.freq(2))
      SAT.fac(1)=SAT.g2/(SAT.g2-1.D0)
      SAT.fac(2)=1.D0/(SAT.g2-1.D0)

      DO i=1, nfreq
        SAT.lamda(i)=VEL_LIGHT/SAT.freq(i)
      END DO
    ! GEOMETRY-FREE COMBINATION
    CASE('GF')
    ! WIDELANE COMBINATION
    CASE('WL')
    CASE DEFAULT
      WRITE(ERROR_UNIT,'(A)') '***ERROR(get_freq): uknown observation type '//TRIM(type)
      CALL exit(1)
  END SELECT

  RETURN

END SUBROUTINE



!*
SUBROUTINE getfreq(csys,satp,freq,ifreq,val)
!!
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
CHARACTER :: csys
CHARACTER(LEN=*) :: satp,freq
INTEGER(IT) :: ifreq
REAL(RL) :: val

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  SELECT CASE(csys)
    CASE('G')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=GPS_L1
        CASE('L2')
          val=GPS_L2
        CASE('L5')
          val=GPS_L5
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for GPS '//TRIM(freq)
      END SELECT
    CASE('R')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=GLS_L1+ifreq*GLS_dL1
        CASE('L2')
          val=GLS_L2+ifreq*GLS_dL2
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for GLONASS '//TRIM(freq)
      END SELECT
    CASE('E')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=GAL_E1
        CASE('L8')
          val=GAL_E5
        CASE('L6')
          val=GAL_E6
        CASE('L5')
          val=GAL_E5a
        CASE('L7')
          val=GAL_E5b
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for GALILEO '//TRIM(freq)
      END SELECT
    CASE('C')
      IF (INDEX(TRIM(satp),'BEIDOU-2') .NE. 0) THEN
        SELECT CASE(TRIM(freq))
          CASE('L1')
            val=BDS_B1
          CASE('L2')
            val=BDS_B1
          CASE('L7')
            val=BDS_B2
          CASE('L6')
            val=BDS_B3
          CASE DEFAULT
            WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for BEIDOU '//TRIM(freq)
        END SELECT
      ELSE
        SELECT CASE(TRIM(freq))
          CASE('L1')
            val=BDS_B1C
          CASE('L2')
            val=BDS_B1I
          CASE('L5')
            val=BDS_B2A
          CASE('L7')
            val=BDS_B2B
          CASE('L8')
            val=BDS_B2C
          CASE('L6')
            val=BDS_B3I
          CASE DEFAULT
            WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for BEIDOU '//TRIM(freq)
        END SELECT
      END IF
    CASE('J')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=QZS_L1
        CASE('L2')
          val=QZS_L2
        CASE('L5')
          val=QZS_L5
        CASE('L6')
          val=QZS_LEX
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for QZSS '//TRIM(freq)
      END SELECT      
    CASE('I')
      SELECT CASE(TRIM(freq))
        CASE('L5')
          val=IRS_L5
        CASE('L9')
          val=IRS_S
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for IRNSS '//TRIM(freq)
      END SELECT      
    CASE('L')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=LEO_L1
        CASE('L2')
          val=LEO_L2          
        CASE('L5')
          val=LEO_L5
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for LEO '//TRIM(freq)
      END SELECT    
    CASE('S')
      SELECT CASE(TRIM(freq))
        CASE('L1')
          val=LEO_L1
        CASE('L2')
          val=LEO_L2          
        CASE('L5')
          val=LEO_L5
        CASE DEFAULT
          WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency for SEO '//TRIM(freq)
      END SELECT         
    CASE DEFAULT
      WRITE(OUTPUT_UNIT,'(A)') '***ERROR(getFreq): unknown frequency '//TRIM(freq)
  END SELECT

  RETURN

END SUBROUTINE
