program test_gim
use ion
implicit none

  character*256 fgim,msg
  real(rl) lat_i, lon_i
  real(rl) vtec,timdif
  real(rl) intv,sod,is
  type(ionex) iod

  logical lexist
  integer*4 i,narg,iargc,mjd
  integer*4 iy,im,id,ih,imin
  
  !! Read arguements
  narg=iargc()
  IF (narg .LT. 1) THEN
    WRITE(*,'(A)') ' USAGE: test_gim -gim fgim -lat lat -lon lon -intv intv'
    CALL exit(1)
  ENDIF

  i=1
  DO WHILE(i .LE. narg)
    CALL getarg(i,msg)
    SELECT CASE(TRIM(msg))
      CASE('-gim')
        i=i+1
        CALL getarg(i,msg)
        fgim=trim(msg)
      CASE('-lat')
        i=i+1
        CALL getarg(i,msg)
        read(msg,*) lat_i
      CASE('-lon')
        i=i+1
        CALL getarg(i,msg)
        read(msg,*) lon_i
      CASE('-intv')
        i=i+1
        CALL getarg(i,msg)
        read(msg,*) intv
      CASE DEFAULT
        WRITE(*,'(A)') '%%%WARNING(get_lsq_args): unknown arguments: '//TRIM(msg)
    END SELECT
    i=i+1
  END DO
  write(*,*) trim(fgim),lat_i,lon_i,intv

  INQUIRE(FILE=fgim, EXIST=lexist)
  IF (lexist .EQ. .FALSE.) THEN
    WRITE(*,'(A)') '%%%WARNING(get_lsq_args): '//TRIM(fgim)//' is not exist.'
  END IF
  CALL read_ionex(fgim,IOD)


  mjd=iod.hd.mjd0
  sod=iod.hd.sod0
  write(*,*) iod.hd.mjd0,iod.hd.sod0,iod.hd.mjd1,iod.hd.sod1
  do while(timdif(mjd,sod,iod.hd.mjd1,iod.hd.sod1) .le. 0.d0)
  write(*,*) mjd,sod
    CALL gim2(mjd,sod,lat_i*deg2rad,lon_i*deg2rad,IOD,vtec)
    
    call mjd2date(mjd,sod,iy,im,id,ih,imin,is)
    WRITE(1000,'(A3,A3,I6,4I3,F10.6,I3,E22.12)') 'AS ','G01',&
                iy,im,id,ih,imin,is,1,vtec

    call timinc(mjd,sod,intv,mjd,sod)
  end do


end program
