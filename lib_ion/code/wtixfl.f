C*
      SUBROUTINE WTIXFL(INXFIL,MAXVAL,MAXLIN,MAXADT,SATSTR,PGMSTR,
     1                  AGESTR,DATSTR,TITTXT,DESTXT,MAPSTR,OBSTXT,
     2                  COMTXT,INXSP1,INXSP2,ADTLBL,ADTLST,INXINF,
     3                  TECMAP,IEXP0 )
CC
CC NAME       :  WTIXFL
CC
CC PURPOSE    :  WRITE IONEX FILE
CC
CC PARAMETERS :
CC         IN :  INXFIL : EXTERNAL IONEX OUTPUT FILE NAME     CH*32
CC               MAXVAL : MAXIMUM NUMBER OF TEC/RMS VALUES    I*4
CC               MAXLIN : MAXIMUM NUMBER OF DES/COM LINES     I*4
CC               MAXADT : MAXIMUM NUMBER OF AUX RECORDS       I*4
CC               SATSTR : SATELLITE SYSTEM                    CH*20
CC               PGMSTR : PROGRAM                             CH*20
CC               AGESTR : AGENCY                              CH*20
CC               DATSTR : DATE AND TIME                       CH*20
CC               TITTXT : TITLE LINE                          CH*60
CC                        =' ': UNDEFINED
CC               DESTXT(I),I=1,..,MAXLIN: MULTI-LINE          CH*60(*)
CC                        DESCRIPTION
CC                        =' ': UNDEFINED
CC               MAPSTR : MAPPING FUNCTION                    CH*4
CC               OBSTXT : OBSERVABLES USED                    CH*60
CC               COMTXT(I),I=1,..,MAXLIN: MULTI-LINE COMMENT  CH*60(*)
CC                        =' ': UNDEFINED
CC               INXSP1 : SPECIFICATIONS 1                    R*8(*)
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
CC               INXSP2 : SPECIFICATIONS 2                    I*4(*)
CC                        =0: UNDEFINED
CC                        (1): NUMBER OF STATIONS
CC                        (2): NUMBER OF SATELLITES
CC               ADTLBL : LABEL OF AUX DATA                   CH*60
CC                        =' ': NO AUX DATA
CC               ADTLST(I),I=1,..,MAXADT: LIST OF AUX DATA    CH*80(*)
CC                        =' ': UNDEFINED
CC               INXINF : INFORMATION TO BE SAVED             I*4(*)
CC                        =0/1: NO/YES
CC                        (1): TEC MAPS
CC                        (2): RMS MAPS
CC               TECMAP(I,J),I=1,..,MAXVAL,J=1,2: TEC/RMS MAP R*8(*,*)
CC                        (IN TECU)
CC                        =999.9: UNDEFINED
CC               IEXP0  : DEFAULT EXPONENT                    I*4
CC                        =-1: RECOMMENDED VALUE
CC                        =99: EXPONENT AUTOMATICALLY ADAPTED
CC
CC SR CALLED  :  OPNERR, OPNFIL, JMT   , RADGMS, WTIXDT, WTIXHD
CC
CC REMARKS    :  IONEX VERSION 1.0
CC
CC AUTHOR     :  S.SCHAER
CC
CC VERSION    :  4.1
CC
CC CREATED    :  16-SEP-97             LAST MODIFIED :  23-AUG-99
CC
CC CHANGES    :  01-OCT-97 : SS: "IEXP0=99" AS SPECIAL OPTION
CC               23-AUG-99 : SS: DECLARE "VORZ" AS CH*1
CC
CC COPYRIGHT  :  ASTRONOMICAL INSTITUTE
CC      1997      UNIVERSITY OF BERNE
CC                    SWITZERLAND
CC
C*
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*80  ADTLST(*)
      CHARACTER*60  TITTXT,DESTXT(*),OBSTXT,COMTXT(*),ADTLBL
      CHARACTER*32  INXFIL
      CHARACTER*20  SATSTR,PGMSTR,AGESTR,DATSTR,INXREC(6),AUXREC
      CHARACTER*4   MAPSTR
      CHARACTER*3   MAPTYP(2)
      CHARACTER*1   VORZ
C
      REAL*8        INXSP1(*),TECMAP(MAXVAL,*)
C
      INTEGER*4     INXSP2(*),INXINF(*)
C
      INCLUDE 'I:COMLFNUM'
C
      DATA MAPTYP/'TEC','RMS'/
C
C SET IONEX RECORDS
C -----------------
      DATA INXREC/'START OF     MAP    ','EPOCH OF CURRENT MAP',
     1            'EXPONENT            ','LAT/LON1/LON2/DLON/H',
     2            'END OF     MAP      ','END OF FILE         '/
C
C OPEN IONEX OUTPUT FILE
C ----------------------
      CALL OPNFIL(LFNLOC,INXFIL,'UNKNOWN','FORMATTED',
     1  ' ',' ',IOSTAT)
      CALL OPNERR(LFNERR,LFNLOC,IOSTAT,INXFIL,'WTIXFL')
C
C WRITE IONEX HEADER
C ------------------
      CALL WTIXHD(MAXLIN,MAXADT,SATSTR,PGMSTR,AGESTR,DATSTR,
     1            TITTXT,DESTXT,MAPSTR,OBSTXT,COMTXT,INXSP1,
     2            INXSP2,ADTLBL,ADTLST,IEXP0 )
C
C WRITE TEX/RMS MAPS
C ------------------
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
      DO 110 ITYP=1,2
        IF (INXINF(ITYP).EQ.0) GOTO 110
        IVAL=1
C
        DO 120 IMAP=1,NMAP
          XEPO=INXSP1(1)+(IMAP-1)*XINT
C
          CALL JMT(XEPO,IYEAR,IMONTH,XDAY)
          IDAY=IDINT(XDAY)
          CALL RADGMS(3,XDAY,VORZ,IHOUR,IMIN,XSEC)
          ISEC=IDNINT(XSEC)
C
C WRITE RECORD OPENING TEC/RMS MAP
          AUXREC=INXREC(1)
          WRITE(AUXREC(10:12),'(A3)') MAPTYP(ITYP)
C
          WRITE(LFNLOC,910) IMAP,AUXREC
910       FORMAT(I6,54X,A20)
C
C WRITE EPOCH OF CURRENT TEC/RMS MAP
          WRITE(LFNLOC,920) IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC,
     1      INXREC(2)
920       FORMAT(6I6,24X,A20)
C
          DO 130 IHGT=1,NHGT
            XHGT=INXSP1(12)+(IHGT-1)*INXSP1(14)
C
C AUTOMATICALLY ADAPT EXPONENT
            IF (IEXP0.EQ.99) THEN
              IEXP=-99
              DO 135 ITEC=IVAL,IVAL+NLAT*NLON-1
                XTEC=TECMAP(ITEC,ITYP)
                IF (XTEC.NE.0.D0 .AND. XTEC.NE.999.9D0)
     1            IEXP=MAX0(IDINT(DLOG10(DABS(XTEC)/9998.D0)-
     2            100)+100,IEXP)
135           CONTINUE
              IF (IEXP.EQ.-99) IEXP=-1
C
C WRITE ACTIVE EXPONENT
              WRITE(LFNLOC,910) IEXP,INXREC(3)
            ELSE
              IEXP=IEXP0
            ENDIF
C
            DO 140 ILAT=1,NLAT
              XLAT=INXSP1(4)+(ILAT-1)*INXSP1(6)
C
C WRITE TEC/RMS DATA SEQUENCE
              WRITE(LFNLOC,930) XLAT,(INXSP1(I),I=7,9),XHGT,INXREC(4)
930           FORMAT(2X,5F6.1,28X,A20)
C
              CALL WTIXDT(NLON,IEXP,TECMAP(IVAL,ITYP))
              IVAL=IVAL+NLON
140         CONTINUE
130       CONTINUE
C
C WRITE RECORD CLOSING TEC/RMS MAP
          AUXREC=INXREC(5)
          WRITE(AUXREC(8:10),'(A3)') MAPTYP(ITYP)
C
          WRITE(LFNLOC,940) IMAP,AUXREC
940       FORMAT(I6,54X,A20)
120     CONTINUE
110   CONTINUE
C
C WRITE LAST RECORD CLOSING IONEX OUTPUT FILE
C -------------------------------------------
      WRITE(LFNLOC,900) INXREC(6)
900   FORMAT(60X,A20)
C
C CLOSE IONEX OUTPUT FILE
C -----------------------
      CLOSE (UNIT=LFNLOC)
C
      RETURN
      END
