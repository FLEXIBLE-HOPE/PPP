
subroutine uclsrp(xhat,yhat,zhat,phat,satpos_eci,sunpos_eci,shadow_factor,acc_eci)

! these variables need to be global in other program
USE tables
implicit none

double precision shadow_factor
    
!the unit vector of BFS axis in ECI, convertion between ECI and BFS
double precision, dimension(3) :: phat,xhat,yhat,zhat,sunpos_eci,satpos_eci,acc_eci

character*1024 x_filename, y_filename, z_filename
double precision, dimension(181,361) :: xdata , ydata, zdata
integer*4 i

logical first
data first /.true./
save first,xdata,ydata,zdata
   
  if (first .eq. .true.) then
     first = .false.
    
    x_filename=f_tablefilename('uclgdx')
    y_filename=f_tablefilename('uclgdy')
    z_filename=f_tablefilename('uclgdz')

    ! call reading grid files
    call ucl_grid_reader(xdata, x_filename)
    call ucl_grid_reader(ydata, y_filename)
    call ucl_grid_reader(zdata, z_filename)

  end if

  ! from igs convension to actual definitation
  !print *,'zdata(2,3)', zdata(2,3)

  call getAcc(satpos_eci,sunpos_eci,xdata,ydata,zdata,phat,xhat,yhat,zhat,shadow_factor,acc_eci)

  !print *, acc_eci

  return

end subroutine
