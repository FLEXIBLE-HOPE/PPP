C     Last change:  OO    7 Feb 1999    1:56 pm
C     Subroutine for computing time of week (in seconds) and GPS-week
C     from Gregorian date, hour, min and sec.
C     (NOTE: Leap seconds are disregarded here, GPS = UTC + leap secs)
C     GPS-week=0 , time of week=0  at Sunday 06.01.1980 00:00 o'clock
C
C     Author  : Ola Ovstedal
C     Created : 06.08.1991                Last modified: 07.02.1999
C     ------------------------------------------------------------------
C     Changes:
C     07.02.1999: Changed INTEGER*2 --> INTEGER, AND INT2 --> INT
C                 Compiled With Lahey F90 4.5
C                 Use 4 digits for year
C     ------------------------------------------------------------------
      SUBROUTINE GPSTIME(year,month,day,hour,min,sec,week,tow)
C     INTEGER*2 year,month,day,hour,min,mday(12),i,nday  ! 4/2.99
      INTEGER   year,month,day,hour,min,mday(0:12),i,nday
      REAL*8 sec,tow,doweek,week
      DATA (mday(i),i=0,12) /0,31,59,90,120,151,181,212,243,273,304,334,
     +                       365/
C     ***Tot. nr of days since 05.01.1980*******************************
      nday=(year-1980)*365+mday(month-1)+day-5                  ! 7/2.99
C     **Correction for leap-years***************************************
      nday=nday+(year-1980)/4                                   ! 7/2.99
      IF (MOD((year-1980),4).EQ.0 .AND. month.LT.3) nday=nday-1 ! 7/2.99
C     **GPS-week, day of week and tow (sec)*****************************
C     week=INT2(nday/7)             ! 4/2.96
      week = INT(nday/7)
      doweek=((nday/7.0D0)-week)*7.0D0+1.0D0
      tow=(doweek-1.0D0)*86400.0D0+hour*3600.0D0+min*60.0D0+sec
      RETURN
      END 
