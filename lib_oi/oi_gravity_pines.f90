!*
SUBROUTINE oi_gravity_pines(lpart,mjd,xsat,mate2j,ndegree,dce,dse,n_ocean_degree,dco,dso,n_opt_degree,dcop,dsop, &
              n_aod_degree,dca,dsa,lgfm,mindeg,maxdeg,ltog,npar,acc,amat,cmat)
!!
!!
!*
USE const
USE tables
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) :: ndegree,n_ocean_degree
INTEGER(IT) :: n_opt_degree,n_aod_degree
INTEGER(IT) :: mindeg,maxdeg,npar
LOGICAL(LG) :: lpart,lgfm
INTEGER(IT) :: ltog(0:MAXGRADEG,0:MAXGRADEG,2)
REAL(RL) :: mjd,mate2j(3,3),xsat(1:*)
REAL(RL) :: dce(4,0:4),dse(4,0:4)
REAL(RL) :: dco(MAXOCNDEG,0:MAXOCNDEG)
REAL(RL) :: dso(MAXOCNDEG,0:MAXOCNDEG)
REAL(RL) :: dcop(MAXOPTDEG,0:MAXOPTDEG)
REAL(RL) :: dsop(MAXOPTDEG,0:MAXOPTDEG)
REAL(RL) :: dca(MAXAODDEG,0:MAXAODDEG)
REAL(RL) :: dsa(MAXAODDEG,0:MAXAODDEG)
REAL(RL) :: acc(1:*),amat(3,3),cmat(1:*)

  !*
  ! The local variables
  !!-----------------------

  INTEGER(IT) :: i,j,k,l,n,m,ierr
  INTEGER(IT) :: lfn,n_gravity_degree

  REAL(RL) ::  cp(4),sp(4)
  REAL(RL), ALLOCATABLE :: cpn(:,:),spn(:,:),pcnm(:,:,:),psnm(:,:,:)

  REAL(RL) :: r,s,t,u,rat1,rat2
  REAL(RL) :: rmjd,gm,radius,xpm,ypm
  REAL(RL) :: ratio01,ratio11,ratio12,ratio22
  REAL(RL) :: dummy,Dnm,Enm,Fnm,Gnm,Hnm

  REAL(RL) :: Cnm(0:MAXGRADEG,0:MAXGRADEG)
  REAL(RL) :: Snm(0:MAXGRADEG,0:MAXGRADEG)
  REAL(RL) :: Cnms(0:MAXGRADEG,0:MAXGRADEG)
  REAL(RL) :: Snms(0:MAXGRADEG,0:MAXGRADEG)
  REAL(RL) :: A(0:MAXGRADEG+2,0:MAXGRADEG+2)
  REAL(RL) :: Rm(-2:MAXGRADEG),Im(-2:MAXGRADEG)
  REAL(RL) :: Ro(0:MAXGRADEG),SUM_F(4),SUM_Fn(4)
  REAL(RL) :: SUM_P(4,4), SUM_Pn(4,4)

  REAL(RL) :: Ctrnd(0:50,0:50),Strnd(0:50,0:50)
  REAL(RL) :: CCpe1(0:50,0:50),SCpe1(0:50,0:50)
  REAL(RL) :: CSpe1(0:50,0:50),SSpe1(0:50,0:50)
  REAL(RL) :: CCpe2(0:50,0:50),SCpe2(0:50,0:50)
  REAL(RL) :: CSpe2(0:50,0:50),SSpe2(0:50,0:50)
  REAL(RL) :: CCpe3,CSpe3
  REAL(RL) :: cos1,sin1,cos2,sin2,cos3,sin3

  CHARACTER(LEN_STRING) :: line,flnegm

  LOGICAL(LG) :: lfirst,lexist
  DATA lfirst /.TRUE./

  SAVE lfirst,n_gravity_degree,gm,radius,rmjd,Cnms,Snms,Cnm,Snm
  SAVE Ctrnd,Strnd,CCpe1,SCpe1,CSpe1,SSpe1,CCpe2,SCpe2,CSpe2,SSpe2,CCpe3,CSpe3

  !*
  ! The function called
  !!------------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!------------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    Cnm=0.d0
    Snm=0.d0

    ! For test case for tide model is on and gravity is off
    flnegm=f_tableFileName('EGM')
    lexist=.TRUE.
    INQUIRE(FILE=flnegm,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_gravity_pines): '//TRIM(flnegm)//' is not exist.'
      CALL exit(1)
    END IF
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnegm)

    n_gravity_degree=0
    rmjd=0.d0
    gm=0.d0
    radius=0.d0
    line=' '
    DO WHILE (INDEX(line,'+coefficient') .EQ. 0)
      READ(lfn,'(A)') line
      IF (line(1:5) .EQ. 'EARTH') THEN
        READ(line(6:), *, IOSTAT=ierr) gm,radius,i
        IF (ierr .NE. 0 ) THEN
          WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_gravity_pines): read '//TRIM(line)
          CALL exit(1)
        END IF
        gm=gm*1d-9
        radius=radius*1d-3
        j=i/10000
        m=(i-j*10000)/100
        n=i-j*10000-m*100
        rmjd=DBLE(modified_julday(n,m,j))
      END IF
    END DO
    IF (rmjd.EQ.0.d0 .OR. gm.EQ.0.d0 .OR. radius.EQ.0.d0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_gravity_pines): there is no GM, EARTH RADIUS, or TIME for gravity model'
      CALL exit(1)
    END IF

    IF (INDEX(flnegm,'EIGEN_6C') .NE. 0) THEN
      Ctrnd=0.d0
      Strnd=0.d0
      CCpe1=0.d0
      SCpe1=0.d0
      CSpe1=0.d0
      SSpe1=0.d0
      CCpe2=0.d0
      SCpe2=0.d0
      CSpe2=0.d0
      SSpe2=0.d0
      CCpe3=0.d0
      CSpe3=0.d0

      DO WHILE (INDEX(line,'-coefficient') .EQ. 0 )
        IF (INDEX(line,'+coefficient') .EQ. 0) THEN
          READ(line(5:),*,IOSTAT=ierr) n,m
          IF (n.GT.MAXGRADEG .OR. m.GT.MAXGRADEG) THEN
            GOTO 10
          ELSE
            IF (n_gravity_degree .LT. n) n_gravity_degree=n
          END IF
        END IF
        IF (line(1:3) .EQ. 'gfc') THEN
          READ(line(5:),*,IOSTAT=ierr) n,m,Cnm(n,m),Snm(n,m)
        ELSE IF (line(1:4) .EQ. 'trnd') THEN
          READ(line(5:),*,IOSTAT=ierr) n,m,Ctrnd(n,m),Strnd(n,m)
        ELSE IF (line(1:4) .EQ. 'acos') THEN
          IF (line(75:77) .EQ. '1.0') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CCpe1(n,m),SCpe1(n,m)
          ELSE IF (line(75:77) .EQ. '0.5') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CCpe2(n,m),SCpe2(n,m)
          ELSE IF (line(75:81) .EQ. '18.6129') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CCpe3
          END IF
        ELSE IF (line(1:4) .EQ. 'asin') THEN
          IF (line(75:77) .EQ. '1.0') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CSpe1(n,m),SSpe1(n,m)
          ELSE IF (line(75:77) .EQ. '0.5') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CSpe2(n,m),SSpe2(n,m)
          ELSE IF (line(75:81) .EQ. '18.6129') THEN
            READ(line(5:),*,IOSTAT=ierr) n,m,CSpe3
          END IF
        END IF
10      READ(lfn,'(A)',IOSTAT=ierr) line
      END DO
      CLOSE(lfn)
    ELSE
      !! EGM2008
      ierr=0
      DO WHILE (INDEX(line,'-coefficient') .EQ. 0)
        IF (line(1:1) .EQ. ' ') THEN
          READ(line, *, IOSTAT=ierr) n,m
          IF (n.GT.MAXGRADEG .OR. m.GT.MAXGRADEG) CYCLE
          READ(line, *, IOSTAT=ierr) n,m,Cnm(n,m),Snm(n,m)
          IF (n_gravity_degree .LT. n) n_gravity_degree=n
        END IF
20      READ(lfn,'(A)') line
      ENDDO
      CLOSE(lfn)
    END IF

    !!! IERS Conventions 2010 (zero-tide)
    !! Cnm(2,0)=-0.48416948d-3
    !! Cnm(3,0)=0.9571612d-6
    !! Cnm(4,0)=0.5399659d-6

    A=0.d0
    Cnms=Cnm
    Snms=Snm

    IF (lgfm .EQ. .TRUE.) THEN
      IF (maxdeg .gt. ndegree) then
        WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_gravity_pines): the partial gravity is greater then the background gravity max. degree'
        CALL exit(1)
      END IF
    END IF
  END IF

  IF (ndegree .GT. n_gravity_degree) THEN
    WRITE(OUTPUT_UNIT,'(A,I5)') '#WARNING(oi_gravity_pines): the maxmum degree is ',n_gravity_degree
    ndegree=n_gravity_degree
  END IF

  ! Restore the original geopotential coefficients for that link to tidal correction.
  DO n=0, MAX(4,ndegree,n_ocean_degree,n_opt_degree,n_aod_degree)
    DO m=0, n
      Cnm(n,m)=Cnms(n,m)
      Snm(n,m)=Snms(n,m)
    END DO
  END DO

  ! Initial the parameter for gravity partical
  IF (lgfm .EQ. .TRUE.) THEN
    ALLOCATE(cpn(4,0:maxdeg+2),STAT=ierr)
    IF (ierr .NE. 0) GOTO 100
    
    ALLOCATE(spn(4,0:maxdeg+2),STAT=ierr)
    IF (ierr .NE. 0) GOTO 100

    ALLOCATE(pcnm(6,-2:maxdeg+2,-2:maxdeg+2),STAT=ierr)
    IF (ierr .NE. 0) GOTO 100

    ALLOCATE(psnm(6,-2:maxdeg+2,-2:maxdeg+2),STAT=ierr)
    IF (ierr .NE. 0) GOTO 100

    cp=0.d0
    sp=0.d0
    cpn=0.d0
    spn=0.d0
    pcnm=0.d0
    psnm=0.d0
  END IF

  IF (INDEX(flnegm,'EIGEN_6C') .NE. 0) THEN
    cos1=DCOS((mjd-rmjd)/365.25d0*2.d0*PI)
    sin1=DSIN((mjd-rmjd)/365.25d0*2.d0*PI)
    cos2=DCOS((mjd-rmjd)/365.25d0*4.d0*PI)
    sin2=DSIN((mjd-rmjd)/365.25d0*4.d0*PI)
    cos3=DCOS((mjd-rmjd)/365.25d0*2.d0*PI/18.6129)
    sin3=DSIN((mjd-rmjd)/365.25d0*2.d0*PI/18.6129)

    k=MAX(4,ndegree,n_ocean_degree,n_opt_degree,n_aod_degree)
    IF (k .GT. 50) k=50
    DO n=0, k
      DO m=0, n
        Cnm(n,m)=Cnms(n,m)+Ctrnd(n,m)*(mjd-rmjd)/365.25d0+CCpe1(n,m)*cos1+CSpe1(n,m)*sin1+CCpe2(n,m)*cos2+CSpe2(n,m)*sin2
        Snm(n,m)=Snms(n,m)+Strnd(n,m)*(mjd-rmjd)/365.25d0+SCpe1(n,m)*cos1+SSpe1(n,m)*sin1+SCpe2(n,m)*cos2+SSpe2(n,m)*sin2
      END DO
    END DO
    Cnm(2,0)=Cnm(2,0)+CCpe3*cos3+CSpe3*sin3
  END IF

  !! Get the max degree for gravity, solid earth tide, ocean tide
  ndegree=MAX(4,ndegree,n_ocean_degree,n_opt_degree,n_aod_degree)

  ! JG: The following correction is only for EGM2008 as description of IERS2010
  IF ((ndegree.GT.2) .AND. (INDEX(flnegm,'EIGEN_6C').EQ.0)) THEN
    ! IERS Conventions 2010, P.80, Table 6.2
    CALL mean_pole('IERS2010',mjd,xpm,ypm)
    !CALL mean_pole('IERS2010_C21S21',mjd,xpm,ypm)
    !Cnm(2,1)= (DSQRT(3.d0)*xpm*Cnms(2,0)-xpm*Cnms(2,2)+ypm*Snms(2,2))*ARCSEC2RAD
    !Snm(2,1)=(-DSQRT(3.d0)*ypm*Cnms(2,0)-ypm*Cnms(2,2)-xpm*Snms(2,2))*ARCSEC2RAD
    !Cnm(2,1)= (DSQRT(3.d0)*xpm*-0.48416948d-3-xpm*2.4393836d-6+ypm*-1.4002737d-6)*ARCSEC2RAD
    !Snm(2,1)=(-DSQRT(3.d0)*ypm*-0.48416948d-3-ypm*2.4393836d-6-xpm*-1.4002737d-6)*ARCSEC2RAD
    Cnm(2,1)=(DSQRT(3.d0)*xpm*-0.48416948d-3-xpm*Cnms(2,2)+ypm*Snms(2,2))*ARCSEC2RAD
    Snm(2,1)=(-DSQRT(3.d0)*ypm*-0.48416948d-3-ypm*Cnms(2,2)-xpm*Snms(2,2))*ARCSEC2RAD

    ! for using updated mean pole, no correction for C21, S21

    !Cnm(2,1)=Cnm(2,1)-3.37d-12*(mjd-rmjd)/365.25d0
    !Snm(2,1)=Snm(2,1)+16.06d-12*(mjd-rmjd)/365.25d0

    Cnm(2,0)=Cnms(2,0)+11.6d-12*(mjd-rmjd)/365.25d0
  END IF
  IF ((ndegree.GT.3) .AND. (INDEX(flnegm,'EIGEN_6C').EQ.0)) Cnm(3,0)=Cnms(3,0)+4.9d-12*(mjd-rmjd)/365.25d0
  IF ((ndegree.GT.4) .AND. (INDEX(flnegm,'EIGEN_6C').EQ.0)) Cnm(4,0)=Cnms(4,0)+4.7d-12*(mjd-rmjd)/365.25d0
  
  !! correction due to solid earth tide
  !Cnm(2,1)=Cnms(2,1)
  !Snm(2,1)=Snms(2,1)
  !Cnm(2,0)=Cnms(2,0)
  !Cnm(3,0)=Cnms(3,0)
  !Cnm(4,0)=Cnms(4,0)
  DO n=2, 4
    DO m=0, n
      Cnm(n,m)=Cnm(n,m)+dce(n,m)
      Snm(n,m)=Snm(n,m)+dse(n,m)
    END DO
  END DO

  !! correction due to ocean tide
  DO n=2, n_ocean_degree
    DO m=0, n
      Cnm(n,m)=Cnm(n,m)+dco(n,m)
      Snm(n,m)=Snm(n,m)+dso(n,m)
    END DO
  END DO

  !! correction due to ocean pole tide
  DO n=2, n_opt_degree
    DO m=0, n
      Cnm(n,m)=Cnm(n,m)+dcop(n,m)
      Snm(n,m)=Snm(n,m)+dsop(n,m)
    END DO
  END DO

  !! correction due to atmosphere and ocean dealising
  DO n=2, n_aod_degree
    DO m=0, n
      Cnm(n,m)=Cnm(n,m)+dca(n,m)
      Snm(n,m)=Snm(n,m)+dsa(n,m)
    END DO
  END DO

  r=DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)
  s=xsat(1)/r
  t=xsat(2)/r
  u=xsat(3)/r


  A(0,0)=1.d0
  A(1,1)=SQRT(3.d0)
  DO n=2, ndegree+2
    A(n,n)=A(n-1,n-1)*SQRT(1.d0+0.5d0/n)
  END DO
  DO j=1, ndegree+2
    rat1=SQRT(DBLE((j+j+1)*(j+1))/DBLE(j*(4*j+2)))
    rat2=SQRT(DBLE((j+j+1)*(j-1))/DBLE(j*(4*j-2)))
    A(j,0)=rat1*u*A(j,1)-rat2*A(j-1,1)
    IF (j+1 .LE. ndegree+2) THEN
      rat1=SQRT(DBLE((j+j+3)*j)/DBLE((j+j+1)*(j+2)))
      rat2=SQRT(DBLE((4*j+6)*(j+1))/DBLE((j+j+1)*(j+2)))
      A(j+1,1)=rat2*A(j,0)+rat1*u*A(j,1)
    END IF
    DO n=j+2, ndegree+2
      m=n-j
      rat1=SQRT(DBLE((n+n+1)*(n-m))/DBLE((n+n-1)*(n+m)))
      rat2=rat1*SQRT(DBLE(n+m-1)/DBLE(n-m))
      A(n,m)=rat2*A(n-1,m-1)+rat1*u*A(n-1,m)
    END DO
  END DO

  Rm(-2)=1.d0
  Rm(-1)=1.d0
  Rm( 0)=1.d0
  Im(-2)=0.d0
  Im(-1)=0.d0
  Im( 0)=0.d0
  DO m=1, ndegree
    Rm(m)=s*Rm(m-1)-t*Im(m-1)
    Im(m)=s*Im(m-1)+t*Rm(m-1)
  END DO

  dummy=radius/r
  Ro(0)=gm/r/radius
  DO n=1, ndegree+2
    Ro(n)=dummy*Ro(n-1)
  END DO

  DO i=1, 4
    SUM_F(i)=0.d0
    IF (lpart .EQ. .TRUE.)  THEN
      DO j=1, 4
        SUM_P(i,j)=0.d0
      END DO
    END IF
  END DO

  ! Loop over n
  DO n=2, ndegree
    DO i=1, 4
      SUM_Fn(i)=0.d0
      IF (lpart .EQ. .TRUE.) THEN
        DO j=1, 4
          SUM_Pn(i,j)=0.d0
        END DO
      END IF
      IF (lgfm .EQ. .TRUE.) THEN
        DO m=0, maxdeg+2
          cpn(i,m)=0.d0
          spn(i,m)=0.d0
        END DO
      END IF
    END DO

    DO m=0, n

      ! (De)normalization
      ratio01=SQRT(DBLE((n+m+1)*(n-m)))                   ! N(n,m)/N(n,m+1)
      ratio11=SQRT(DBLE((2*n+1)*(n+m+1)*(n+m+2))/(2*n+3)) ! N(n,m)/N(n+1,m+1)
      IF (m .EQ. 0) THEN
        ratio01=ratio01*SQRT(0.5d0)
        ratio11=ratio11*SQRT(0.5d0)
      END IF
      ratio12=ratio11*SQRT(DBLE((n-m)*(n+m+3)))           ! N(n,m)/N(n+1,m+2)
      ratio22=ratio11*SQRT(DBLE((2*n+3)*(n+m+3)*(n+m+4))/(2*n+5)) !N/N(n+2,m+2)

      Dnm=Cnm(n,m)*Rm(m)  +Snm(n,m)*Im(m)
      Enm=Cnm(n,m)*Rm(m-1)+Snm(n,m)*Im(m-1)
      Fnm=Snm(n,m)*Rm(m-1)-Cnm(n,m)*Im(m-1)

      ! inite the partial of gravity field coefficients. 
      IF (lgfm.EQ..TRUE. .AND. n.LE.maxdeg) THEN
        cpn(3,m)=Rm(m)
        spn(3,m)=Im(m)
        cpn(1,m)=Rm(m-1)
        spn(1,m)=Im(m-1)
        cpn(2,m)=-Im(m-1)
        spn(2,m)=Rm(m-1)
      END IF

      IF (lpart .EQ. .TRUE.) THEN
        Gnm=Cnm(n,m)*Rm(m-2)+Snm(n,m)*Im(m-2)
        Hnm=Snm(n,m)*Rm(m-2)-Cnm(n,m)*Im(m-2)
      END IF

      SUM_Fn(1)=SUM_Fn(1)+A(n,m)*m*Enm
      SUM_Fn(2)=SUM_Fn(2)+A(n,m)*m*Fnm
      SUM_Fn(3)=SUM_Fn(3)+A(n,m+1)*Dnm*ratio01
      SUM_Fn(4)=SUM_Fn(4)-A(n+1,m+1)*Dnm*ratio11

      IF (lgfm.EQ..TRUE. .AND. n.LE.maxdeg) THEN
        CPn(1,m)=A(n,m)*m*CPn(1,m)
        SPn(1,m)=A(n,m)*m*SPn(1,m)
        CPn(2,m)=A(n,m)*m*CPn(2,m)
        SPn(2,m)=A(n,m)*m*SPn(2,m)
        CPn(4,m)=-A(n+1,m+1)*ratio11*CPn(3,m)
        SPn(4,m)=-A(n+1,m+1)*ratio11*SPn(3,m)
        CPn(3,m)=A(n,m+1)*ratio01*CPn(3,m)
        SPn(3,m)=A(n,m+1)*ratio01*SPn(3,m)
      END IF

      IF (lpart .EQ. .TRUE.) THEN
        SUM_Pn(1,1)=SUM_Pn(1,1)+A(n,m)*m*(m-1.d0)*Gnm
        SUM_Pn(1,2)=SUM_Pn(1,2)+A(n,m)*m*(m-1.d0)*Hnm
        SUM_Pn(1,3)=SUM_Pn(1,3)+A(n,m+1)*m*Enm*ratio01
        SUM_Pn(1,4)=SUM_Pn(1,4)-A(n+1,m+1)*m*Enm*ratio11
        SUM_Pn(2,3)=SUM_Pn(2,3)+A(n,m+1)*m*Fnm*ratio01
        SUM_Pn(2,4)=SUM_Pn(2,4)-A(n+1,m+1)*m*Fnm*ratio11
        SUM_Pn(3,4)=SUM_Pn(3,4)-A(n+1,m+2)*Dnm*ratio12
        SUM_Pn(4,4)=SUM_Pn(4,4)+A(n+2,m+2)*Dnm*ratio22
      END IF
    END DO

    DO i=1, 4
      SUM_F(i)=SUM_F(i)+Ro(n+1)*SUM_Fn(i)
    END DO

    IF (lgfm.EQ..TRUE. .AND. n.LE.maxdeg) THEN
      DO m=0, n
        PCnm(1,m,n)=ro(n+1)*(CPn(1,m)+CPn(4,m)*s)
        PCnm(2,m,n)=ro(n+1)*(CPn(2,m)+CPn(4,m)*t)
        PCnm(3,m,n)=ro(n+1)*(CPn(3,m)+CPn(4,m)*u)
        PSnm(1,m,n)=ro(n+1)*(SPn(1,m)+SPn(4,m)*s)
        PSnm(2,m,n)=ro(n+1)*(SPn(2,m)+SPn(4,m)*t)
        PSnm(3,m,n)=ro(n+1)*(SPn(3,m)+SPn(4,m)*u)
      END DO
    END IF

    IF (lpart .EQ. .TRUE.) THEN
      DO i=1, 4
        DO j=i, 4
          SUM_P(i,j)=SUM_P(i,j)+Ro(n+2)/radius*SUM_Pn(i,j)
        END DO
      END DO
    END IF
  END DO

  SUM_F(1)=SUM_F(1)+SUM_F(4)*s
  SUM_F(2)=SUM_F(2)+SUM_F(4)*t
  SUM_F(3)=SUM_F(3)+SUM_F(4)*u

  ! Acceleration to inertial system
  DO i=1, 3
    dummy=mate2j(i,1)*SUM_F(1)+mate2j(i,2)*SUM_F(2)+mate2j(i,3)*SUM_F(3)
    acc(i)=acc(i)+dummy
  END DO

  ! To inertial system
  IF (lgfm .EQ. .TRUE.) THEN
    DO n=2, maxdeg
      DO m=0, n
        DO i=1, 3
          PCnm(i+3,m,n)=(mate2j(i,1)*PCnm(1,m,n)+mate2j(i,2)*PCnm(2,m,n)+mate2j(i,3)*PCnm(3,m,n))
          PSnm(i+3,m,n)=(mate2j(i,1)*PSnm(1,m,n)+mate2j(i,2)*PSnm(2,m,n)+mate2j(i,3)*PSnm(3,m,n))
        END DO
      END DO
    END DO 

    DO n=mindeg,maxdeg
      DO m=0, n
        j=(ltog(n,m,1)+npar-6-1)*3
        DO k=1, 3
          cmat(j+k)=PCnm(k+3,m,n)
        END DO
      END DO

      DO m=1, n
        j=(ltog(n,m,2)+npar-6 -1)*3
        DO k=1, 3
          cmat(j+k)=PSnm(k+3,m,n)
        END DO
      END DO
    END DO

    DEALLOCATE(cpn)
    DEALLOCATE(spn)
    DEALLOCATE(pcnm)
    DEALLOCATE(psnm)
  END IF

  IF (lpart .EQ. .TRUE.) THEN
    SUM_Pn(1,1) = SUM_P(1,1)+s*s*SUM_P(4,4)+SUM_F(4)/r+2.d0*s*SUM_P(1,4)
    SUM_Pn(1,2) = SUM_P(1,2)+s*t*SUM_P(4,4)+s*SUM_P(2,4)+t*SUM_P(1,4)
    SUM_Pn(1,3) = SUM_P(1,3)+s*u*SUM_P(4,4)+s*SUM_P(3,4)+u*SUM_P(1,4)
    SUM_Pn(2,2) =-SUM_P(1,1)+t*t*SUM_P(4,4)+SUM_F(4)/r+2.d0*t*SUM_P(2,4)
    SUM_Pn(2,3) = SUM_P(2,3)+t*u*SUM_P(4,4)+u*SUM_P(2,4)+t*SUM_P(3,4)
    SUM_Pn(3,3) = -(SUM_Pn(1,1)+SUM_Pn(2,2))
    SUM_Pn(2,1) = SUM_Pn(1,2)
    SUM_Pn(3,1) = SUM_Pn(1,3)
    SUM_Pn(3,2) = SUM_Pn(2,3)

    ! Partial rotate to inertial
    DO i=1, 3
      DO j=1, 3
        dummy=0.d0
        DO k=1, 3
          DO l=1, 3
            dummy=dummy+mate2j(i,k)*SUM_Pn(k,l)*mate2j(j,l)
          END DO
        END DO
        amat(i,j)=amat(i,j)+dummy
      END DO
    END DO
  END IF

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_gravity_pines): memory allocation'
  CALL exit(1)


END SUBROUTINE
