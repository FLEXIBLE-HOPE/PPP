C     Last change:  OO   20 Nov 2000    1:01 pm
      PROGRAM TEST_KLOB
C     Driver program to demonstrate the klobuchar.for subroutine
C     using a numerical example.
C     Note: In this example leap seconds are ignored, i.e., it  
C     is assumed that UTC time equals GPS Time.  If 13.0 leap
C     seconds were added to the 20:45 time in the numerical example,
C     it would change the ion delay by less than one millimeter.

      IMPLICIT NONE
      REAL*8 alfa(4),beta(4),sec,week,tow,elev,fi,lambda,azimuth,dIon1
      INTEGER*4 i,year,month,day,hour,MIN
      CHARACTER*80 string80
      CHARACTER*12 navfile,outfile

      OPEN(10,FILE='test_klo.inp')
      READ(10,'(a12)')  navfile               ! File with alphas and
                                              ! betas (RINEX NAV format)
      READ(10,'(a12)')  outfile               ! Out-file
      READ(10,*) year,month,day,hour,MIN,sec  ! Time (UTC)
      READ(10,*) elev                         ! Elevation (deg)
      READ(10,*) azimuth                      ! Azimuth   (deg)
      READ(10,*) fi                           ! Latitude (deg min sec)
      READ(10,*) lambda                       ! Longtitude (deg min sec)
      CLOSE (10)

      OPEN(20,FILE=navfile)
      OPEN(30,FILE=outfile)

                                              ! Read alphas
      string80(61:69)='         '
      DO WHILE(string80(61:69) .NE. 'ION ALPHA')
        READ(20,'(a80)') string80
      END DO
      READ(string80,'(2x,4D12.4)') (alfa(i), i=1,4)

      REWIND 20
                                              ! Read betas
      string80(61:69)='         '
      DO WHILE(string80(61:68) .NE. 'ION BETA')
        READ(20,'(a80)') string80
      END DO
      READ(string80,'(2x,4D12.4)') (beta(i), i=1,4)

      CLOSE (20)

                                              ! Compute Time-of-Week in
                                              ! seconds (UTC)
					      ! (leap secs are ignored)
      CALL GPSTIME(year,month,day,hour,MIN,sec,week,tow)

                                              ! Ionospheric L1 corr.
      CALL KLOBUCHAR(fi,lambda,elev,azimuth,tow,alfa,beta,dIon1)

      WRITE(30,'(a6,f9.0,tr5,a6,f8.3)') ' Tow :',tow, 'dION: ',dIon1
      CLOSE (30)

      STOP
      END
