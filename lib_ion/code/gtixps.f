C*
      SUBROUTINE GTIXPS(INXSP1,IGP,INPARG,OUTARG,IPOS)
CC
CC NAME       :  GTIXPS
CC
CC PURPOSE    :  GET POSITION OF (1) NEAREST OR (2) LOWER LEFT HAND
CC               GRID POINT - OR ARGUMENTS BASED ON POSITION
CC
CC PARAMETERS :
CC         IN :  INXSP1 : SPECIFICATIONS 1                    R*8(*)
CC                        ( 1): EPOCH OF FIRST MAP (IN MJD)
CC                        ( 2): EPOCH OF LAST MAP (IN MJD)
CC                        ( 3): INTERVAL (IN SEC)
CC                        ( 4): FROM LATITUDE
CC                        ( 5): TO LATITUDE
CC                        ( 6): WITH INCREMENT (IN DEG)
CC                        ( 7): FROM LONGITUDE
CC                        ( 8): TO LONGITUDE
CC                        ( 9): WITH INCREMENT (IN DEG)
CC                        (10): ELEVATION CUTOFF (IN DEG)
CC                        (11): BASE RADIUS (IN KM)
CC                        (12): FROM HEIGHT
CC                        (13): TO HEIGHT
CC                        (14): WITH INCREMENT (IN KM)
CC               IGP    : GRID POINT TO BE RETURNED           I*4
CC                        =0: RETURN ARGUMENTS "OUTARG" BASED
CC                            ON POSITION "IPOS"
CC                        =1: NEAREST
CC                        =2: LOWER LEFT HAND
CC               INPARG : INPUT ARGUMENTS                     R*8(*)
CC                        (1): EPOCH (IN MJD)
CC                        (2): LATITUDE (IN DEG)
CC                        (3): LONGITUDE (IN DEG)
CC                        (4): HEIGHT (IN KM)
CC        OUT :  OUTARG : OUTPUT ARGUMENTS                    R*8(*)
CC                        (1): EPOCH (IN MJD)
CC                        (2): LATITUDE (IN DEG)
CC                        (3): LONGITUDE (IN DEG)
CC                        (4): HEIGHT (IN KM)
CC     IN/OUT :  IPOS   : POSITION                            I*4
CC                        =0: UNDEFINED
CC
CC SR CALLED  :  ---
CC
CC REMARKS    :  IONEX VERSION 1.0
CC
CC AUTHOR     :  S.SCHAER
CC
CC VERSION    :  4.1
CC
CC CREATED    :  01-MAR-98             LAST MODIFIED :  08-JUN-98
CC
CC CHANGES    :  08-APR-98 : SS: CONSIDER "IGP"
CC               08-JUN-98 : SS: RETURN ARGUMENTS, IF "IGP=0"
CC
CC COPYRIGHT  :  ASTRONOMICAL INSTITUTE
CC      1998      UNIVERSITY OF BERNE
CC                    SWITZERLAND
CC
C*
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8        INXSP1(*),INPARG(*),OUTARG(*)
C
      INCLUDE 'I:COMLFNUM'
C
C SET DEFAULT RETURN CODE
C -----------------------
      IRC=1
C
C GET NUMBER OF MAPS AND GRID DIMENSIONS
C --------------------------------------
      XINT=INXSP1(3)/86400.D0
      NMAP=IDNINT((INXSP1(2)-INXSP1(1))/XINT)+1
C
      NLAT=IDNINT((INXSP1(5)-INXSP1(4))/INXSP1(6))+1
      NLON=IDNINT((INXSP1(8)-INXSP1(7))/INXSP1(9))+1
C
      IF (INXSP1(14).EQ.0.D0) THEN
        NHGT=1
      ELSE
        NHGT=IDNINT((INXSP1(13)-INXSP1(12))/INXSP1(14))+1
      ENDIF
C
C GET EPOCH, LATITUDE, LONGITUDE, AND HEIGHT BASED ON POSITION
C ------------------------------------------------------------
      IF (IGP.EQ.0) THEN
        IVAL=IPOS
        IMAP=(IVAL-1)/(NHGT*NLAT*NLON)+1
        IVAL=IVAL-(IMAP-1)*NHGT*NLAT*NLON
        IHGT=(IVAL-1)/(NLAT*NLON)+1
        IVAL=IVAL-(IHGT-1)*NLAT*NLON
        ILAT=(IVAL-1)/NLON+1
        IVAL=IVAL-(ILAT-1)*NLON
        ILON=IVAL
C
        OUTARG(1)=INXSP1( 1)+(IMAP-1)*XINT
        OUTARG(2)=INXSP1( 4)+(ILAT-1)*INXSP1( 6)
        OUTARG(3)=INXSP1( 7)+(ILON-1)*INXSP1( 9)
        OUTARG(4)=INXSP1(12)+(IHGT-1)*INXSP1(14)
C
        RETURN
      ENDIF
C
C GET POSITION OF NEAREST OR LOWER LEFT HAND GRID POINT
C -----------------------------------------------------
C
C - EPOCH
      XMAP=(INPARG(1)-INXSP1(1))/XINT+1.D0
      IF (IGP.EQ.1) THEN
        IMAP=IDNINT(XMAP)
      ELSE
        IMAP=IDINT(XMAP)
      ENDIF
      IF (IMAP.GE.1 .AND. IMAP.LE.NMAP) THEN
        OUTARG(1)=INXSP1(1)+(IMAP-1)*XINT
      ELSE
        WRITE(LFNERR,1901) INPARG(1)
1901    FORMAT(/,' ### SR GTIXPS: IRREGULAR EPOCH',
     1    /,16X,'EPOCH: ',F12.5,' MJD',/)
        IRC=0
      ENDIF
C
C - LATITUDE
      XLAT=(INPARG(2)-INXSP1(4))/INXSP1(6)+1.D0
      IF (IGP.EQ.1) THEN
        ILAT=IDNINT(XLAT)
      ELSE
        ILAT=IDINT(XLAT)
      ENDIF
      IF (ILAT.GE.1 .AND. ILAT.LE.NLAT) THEN
        OUTARG(2)=INXSP1(4)+(ILAT-1)*INXSP1(6)
      ELSE
        WRITE(LFNERR,1902) INPARG(2)
1902    FORMAT(/,' ### SR GTIXPS: IRREGULAR LATITUDE',
     1    /,16X,'LATITUDE: ',F8.3,' DEG',/)
        IRC=0
      ENDIF
C
C - LONGITUDE
      XLON=(INPARG(3)-INXSP1(7))/INXSP1(9)+1.D0
      IF (IGP.EQ.1) THEN
        ILON=IDNINT(XLON)
      ELSE
        ILON=IDINT(XLON)
      ENDIF
C
      NCYC=IDNINT(360.D0/DABS(INXSP1(9)))
      IF (ILON.LT.1) THEN
        ILON=ILON+NCYC
      ELSEIF (ILON.GT.NLON) THEN
        ILON=ILON-NCYC
      ENDIF
C
      IF (ILON.GE.1 .AND. ILON.LE.NLON) THEN
        OUTARG(3)=INXSP1(7)+(ILON-1)*INXSP1(9)
      ELSE
        WRITE(LFNERR,1903) INPARG(3)
1903    FORMAT(/,' ### SR GTIXPS: IRREGULAR LONGITUDE',
     1    /,16X,'LONGITUDE: ',F8.3,' DEG',/)
        IRC=0
      ENDIF
C
C - HEIGHT
      IF (INXSP1(14).EQ.0.D0) THEN
        IHGT=1
      ELSE
        XHGT=(INPARG(4)-INXSP1(12))/INXSP1(14)+1.D0
        IF (IGP.EQ.1) THEN
          IHGT=IDNINT(XHGT)
        ELSE
          IHGT=IDINT(XHGT)
        ENDIF
        IF (IHGT.GE.1 .AND. IHGT.LE.NHGT) THEN
          OUTARG(4)=INXSP1(12)+(IHGT-1)*INXSP1(14)
        ELSE
          WRITE(LFNERR,1904) INPARG(4)
1904      FORMAT(/,' ### SR GTIXPS: IRREGULAR HEIGHT',
     1      /,16X,'HEIGHT: ',F8.3,' KM',/)
          IRC=0
        ENDIF
      ENDIF
C
C RETURN POSITION OF GRID POINT
C -----------------------------
      IF (IRC.EQ.1) THEN
        IPOS=(IMAP-1)*NHGT*NLAT*NLON+
     1    (IHGT-1)*NLAT*NLON+
     2    (ILAT-1)*NLON+ILON
      ELSE
        IPOS=0
      ENDIF
C
      RETURN
      END
