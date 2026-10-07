!*
SUBROUTINE getAcc(satpos_eci,sunpos_eci,gridX,gridY,gridZ,phat,xhat,yhat,zhat,shadow_factor,acc_eci)
!declaration
implicit none

double precision shadow_factor
double precision, dimension(3) :: satpos_eci, sunpos_eci,phat,xhat,yhat,zhat,acc_eci
double precision, dimension(181,361) :: gridX,gridY,gridZ    
    
double precision AU,TSI,PI,D2R,R2D,x_acc_bus,y_acc_bus,z_acc_bus
double precision tmp,tsi_real,dis_factor,mass,area_panel,reflectivity_panel,specularity_panel,lon, lat
double precision, dimension(3) :: sunpos_bfs,flux_dir, force_panel_eci,force_bus_bfs,force_bus_eci


  TSI =1368 ! w/m^2
  AU = 149597870.7 !km	
  PI = 3.14159265357
  R2D = 180.0/PI

  D2R = PI/180.0
  flux_dir(1) = satpos_eci(1) - sunpos_eci(1)
  flux_dir(2) = satpos_eci(2) - sunpos_eci(2)
  flux_dir(3) = satpos_eci(3) - sunpos_eci(3)

  tmp = DSQRT((flux_dir(1))**2 + (flux_dir(2))**2 +(flux_dir(3))**2)
  ! here dis_factor should be set according to the distance from sat to sun 
  dis_factor = (AU/tmp)**2
  flux_dir(1) = flux_dir(1)/tmp
  flux_dir(2) = flux_dir(2)/tmp
  flux_dir(3) = flux_dir(3)/tmp

  !print *, 'flux_dir',flux_dir

  tsi_real = TSI*dis_factor*shadow_factor

  !print *, 'dis_factor',dis_factor

  !optical property of solar panel
  !for galileo IOV
  mass = 696.815
  area_panel = 10.82
  reflectivity_panel = 0.245
  specularity_panel = 1.0

  !calculate the lon and lat of sun in bfs
  sunpos_bfs(1) = -(xhat(1)*flux_dir(1) + xhat(2)*flux_dir(2) + xhat(3)*flux_dir(3))
  sunpos_bfs(2) = -(yhat(1)*flux_dir(1) + yhat(2)*flux_dir(2) + yhat(3)*flux_dir(3))
  sunpos_bfs(3) = -(zhat(1)*flux_dir(1) + zhat(2)*flux_dir(2) + zhat(3)*flux_dir(3))
  tmp = DSQRT((sunpos_bfs(1))**2 + (sunpos_bfs(2))**2 +(sunpos_bfs(3))**2)
  sunpos_bfs(1) = sunpos_bfs(1)/tmp
  sunpos_bfs(2) = sunpos_bfs(2)/tmp
  sunpos_bfs(3) = sunpos_bfs(3)/tmp

  lon = ATAN2(sunpos_bfs(2), sunpos_bfs(1) )*R2D
  lat = ASIN(sunpos_bfs(3))*R2D 

  ! get force for the bus, here x_acc_bus, y_acc_bus,z_acc_bus are both in BFS
  call myinterpolation(lon,lat,gridX,x_acc_bus)
  call myinterpolation(lon,lat,gridY,y_acc_bus)
  call myinterpolation(lon,lat,gridZ,z_acc_bus)

  !print *, 'lon',lon,'lat',lat
  !print *, 'x_acc_bus',x_acc_bus
  !print *, 'y_acc_bus',y_acc_bus
  !print *, 'z_acc_bus',z_acc_bus

  !here force_bus is in BFS
  force_bus_bfs(1) = (x_acc_bus)/1386.0*tsi_real
  force_bus_bfs(2) = (y_acc_bus)/1386.0*tsi_real
  force_bus_bfs(3) = (z_acc_bus)/1386.0*tsi_real
 
  !print *, 'force_bus_bfs',force_bus_bfs
  !convert the force_bus_bfs from BFS to ECI
  force_bus_eci(1) = force_bus_bfs(1)*xhat(1) + force_bus_bfs(2)*yhat(1) + force_bus_bfs(3)*zhat(1)
  force_bus_eci(2) = force_bus_bfs(1)*xhat(2) + force_bus_bfs(2)*yhat(2) + force_bus_bfs(3)*zhat(2)
  force_bus_eci(3) = force_bus_bfs(1)*xhat(3) + force_bus_bfs(2)*yhat(3) + force_bus_bfs(3)*zhat(3)

  !calculate the force of solar panel, because phat is in ECI, force_panel is in ECI
  call radiationForce(phat,flux_dir,tsi_real, area_panel, specularity_panel, reflectivity_panel, force_panel_eci)

  !print *, 'force_panel_eci',force_panel_eci
  !print *, 'force_bus_eci',force_bus_eci

  ! get the final result
  acc_eci(1) = (force_bus_eci(1) + force_panel_eci(1))/mass
  acc_eci(2) = (force_bus_eci(2) + force_panel_eci(2))/mass
  acc_eci(3) = (force_bus_eci(3) + force_panel_eci(3))/mass

  !print *, 'acc_eci',acc_eci

  return

end subroutine
