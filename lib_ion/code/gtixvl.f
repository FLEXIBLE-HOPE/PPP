C*
      SUBROUTINE GTIXVL(MAXVAL,INXSP1,INXINF,TECMAP,TECARG,IIS   ,
     1                  TECVAL,IRC   )
CC
CC NAME       :  GTIXVL
CC
CC PURPOSE    :  GET IONEX TEC AND RMS VALUES AS A FUNCTION OF
CC               EPOCH, LATITUDE, LONGITUDE, AND NOMINAL HEIGHT.
CC               FOUR INTERPOLATION STRATEGIES ARE SUPPORTED (SEE
CC               ALSO IONEX MANUAL):
CC                 (1) TAKE NEAREST MAP,
CC                 (2) INTERPOLATE BETWEEN TWO CONSECUTIVE MAPS,
CC                 (3) INTERPOLATE BETWEEN TWO CONSECUTIVE ROTATED
CC                     MAPS, OR
CC                 (4) TAKE NEAREST ROTATED MAP.
CC               A SIMPLE 4-POINT FORMULA IS APPLIED TO INTERPOLATE
CC               BETWEEN THE GRID POINTS.
CC
CC PARAMETERS :
CC         IN :  MAXVAL : MAXIMUM NUMBER OF TEC/RMS VALUES    I*4
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
CC               INXINF : INFORMATION AVAILABLE               I*4(*)
CC                        =0/1: NO/YES
CC                        (1): TEC MAPS
CC                        (2): RMS MAPS
CC               TECMAP(I,J),I=1,..,MAXVAL,J=1,2: TEC/RMS MAP R*8(*,*)
CC                        (IN TECU)
CC                        =999.9: UNDEFINED
CC               TECARG : ARGUMENTS                           R*8(*)
CC                        (1): EPOCH (IN MJD)
CC                        (2): LATITUDE (IN DEG)
CC                        (3): LONGITUDE (IN DEG)
CC                        (4): HEIGHT (IN KM)
CC               IIS    : INTERPOLATION STRATEGY              I*4
CC                        =1: TAKE NEAREST MAP
CC                        =2: INTERPOLATE BETWEEN CONSECUTIVE
CC                            MAPS
CC                        =3: INTERPOLATE BETWEEN CONSECUTIVE
CC                            ROTATED MAPS
CC                        =4: TAKE NEAREST ROTATED MAP
CC        OUT :  TECVAL : TEC AND RMS VALUES (IN TECU)        R*8(*)
CC                        (1): TEC
CC                        (2): RMS
CC               IRC    : RETURN CODE                         I*4
CC                        =0: ERROR OCCURRED
CC                        =1: OK
CC
CC SR CALLED  :  GTIXPS
CC
CC REMARKS    :  IONEX VERSION 1.0
CC
CC AUTHOR     :  S.SCHAER
CC
CC VERSION    :  4.1
CC
CC CREATED    :  20-NOV-97             LAST MODIFIED :  15-APR-98
CC
CC CHANGES    :  ---
CC
CC COPYRIGHT  :  ASTRONOMICAL INSTITUTE
CC      1997      UNIVERSITY OF BERNE
CC                    SWITZERLAND
CC
C*
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8        INXSP1(*),TECMAP(MAXVAL,*)
      REAL*8        TECARG(*),TECVAL(*)
      REAL*8        INPARG(4),OUTARG(4),ORIARG(4)
      REAL*8        TWOEPO(2),TWOFAC(2),PNTVAL(4,2)
C
      INTEGER*4     INXINF(*)
      INTEGER*4     PNTPOS(4)
C
      INCLUDE 'I:COMLFNUM'
C
C SET DEFAULT RETURN CODE
C -----------------------
      IRC=1
C
C DEFINE NUMBER OF MAPS TO BE CONSIDERED
C --------------------------------------
      IF (IIS.EQ.1) THEN
        NMAP=1
      ELSEIF (IIS.EQ.2) THEN
        NMAP=2
      ELSEIF (IIS.EQ.3) THEN
        NMAP=2
      ELSEIF (IIS.EQ.4) THEN
        NMAP=1
      ELSE
        WRITE(LFNERR,1901)
1901    FORMAT(/,' ### SR GTIXVL: INVALID INTERPOLATION STRATEGY',/)
        GOTO 900
      ENDIF
C
C LOOK FOR VALID IONEX MAPS
C -------------------------
      XINT=INXSP1(3)/86400.D0
      IF (NMAP.EQ.1) THEN
        CALL GTIXPS(INXSP1,1,TECARG,OUTARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
C
        TWOEPO(1)=OUTARG(1)
        TWOFAC(1)=1.D0
      ELSE
        CALL GTIXPS(INXSP1,2,TECARG,OUTARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
C
        TWOEPO(1)=OUTARG(1)
        TWOEPO(2)=OUTARG(1)+XINT
        TWOFAC(1)=(TWOEPO(2)-TECARG(1))/XINT
        TWOFAC(2)=(TECARG(1)-TWOEPO(1))/XINT
      ENDIF
C
      TECVAL(1)=0.D0
      TECVAL(2)=0.D0
      DO IMAP=1,NMAP
C
C GET POSITION OF LOWER LEFT HAND GRID POINT E00
C ----------------------------------------------
        INPARG(1)=TWOEPO(IMAP)+XINT/2.D0
        INPARG(2)=TECARG(2)
        IF (IIS.EQ.1 .OR. IIS.EQ.2) THEN
          INPARG(3)=TECARG(3)
          INPARG(4)=TECARG(4)
        ELSE
          INPARG(3)=TECARG(3)+360.D0*(TECARG(1)-TWOEPO(IMAP))
          INPARG(4)=TECARG(4)
        ENDIF
        CALL GTIXPS(INXSP1,2,INPARG,ORIARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
        PNTPOS(1)=IPOS
C
C CHECK HEIGHT ARGUMENT
C ---------------------
        IF (INXSP1(14).NE.0.D0) THEN
          XDIF=TECARG(4)-ORIARG(4)
          IF (XDIF.NE.0.D0) THEN
            WRITE(LFNERR,1902) TECARG(4)
1902        FORMAT(/,' ### SR GTIXVL: IRREGULAR HEIGHT',
     1        /,16X,'HEIGHT: ',F8.3,' KM',/)
            GOTO 900
          ENDIF
        ENDIF
C
C COMPUTE FACTORS P AND Q
C -----------------------
        XP=(INPARG(3)-ORIARG(3))/INXSP1(9)
        XP=DMOD(XP+360.D0/DABS(INXSP1(9)),1.D0)
        XQ=(INPARG(2)-ORIARG(2))/INXSP1(6)
C
        INPARG(1)=TWOEPO(IMAP)
        INPARG(4)=ORIARG(4)
C
C GET E01'S POSITION
        INPARG(2)=ORIARG(2)
        INPARG(3)=ORIARG(3)+INXSP1(9)
        CALL GTIXPS(INXSP1,1,INPARG,OUTARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
        PNTPOS(2)=IPOS
C
C GET E10'S POSITION
        INPARG(2)=ORIARG(2)+INXSP1(6)
        INPARG(3)=ORIARG(3)
        CALL GTIXPS(INXSP1,1,INPARG,OUTARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
        PNTPOS(3)=IPOS
C
C GET E11'S POSITION
        INPARG(2)=ORIARG(2)+INXSP1(6)
        INPARG(3)=ORIARG(3)+INXSP1(9)
        CALL GTIXPS(INXSP1,1,INPARG,OUTARG,IPOS)
        IF (IPOS.EQ.0) GOTO 900
        PNTPOS(4)=IPOS
C
        DO ITYP=1,2
          IF (INXINF(ITYP).EQ.1) THEN
            DO IPNT=1,4
              XVAL=TECMAP(PNTPOS(IPNT),ITYP)
              IF (XVAL.NE.999.9D0) THEN
                PNTVAL(IPNT,ITYP)=XVAL
              ELSE
                WRITE(LFNERR,1903)
1903            FORMAT(/,' ### SR GTIXVL: UNDEFINED TEC/RMS VALUE(S)',/)
                GOTO 900
              ENDIF
            ENDDO
C
            XSUM=(1.D0-XP)*(1.D0-XQ)*PNTVAL(1,ITYP)+
     1        XP*(1.D0-XQ)*PNTVAL(2,ITYP)+
     2        XQ*(1.D0-XP)*PNTVAL(3,ITYP)+
     3        XP*XQ*PNTVAL(4,ITYP)
            TECVAL(ITYP)=TECVAL(ITYP)+TWOFAC(IMAP)*XSUM
          ENDIF
        ENDDO
      ENDDO
C
C OK
C --
      RETURN
C
C ERROR OCCURRED
C --------------
900   CONTINUE
      IRC=0
C
      RETURN
      END
