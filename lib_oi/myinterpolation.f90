!$myinterpolation
SUBROUTINE myinterpolation(lon,lat,grid_data,result)
double precision lon, lat, result
double precision, dimension(181,361) :: grid_data
double precision min_longitude,max_longitude,min_latitude,max_latitude,u_raw,v_raw,tol,del_u,del_v,row1,row2
integer i1,i2,j1,j2,u_floor,v_floor

min_longitude = -180.0
max_longitude = 180.0
min_latitude = -90.0
max_latitude = 90.0
tol=1.0D-15

IF (lon < min_longitude) THEN
	lon = min_longitude
END IF
IF (lon > max_longitude) THEN	
	lon = max_longitude
END IF

IF (lat < min_latitude) THEN
	lat = min_latitude
END IF

IF (lon > max_longitude)  THEN	
	lat = max_latitude
END IF

u_raw = lon - min_longitude
v_raw = lat - min_latitude

u_floor = FLOOR(u_raw)
v_floor = FLOOR(v_raw)

del_u = u_raw - u_floor
del_v = v_raw - v_floor

i1 = v_floor+1
i2 = i1 +1
j1 = u_floor+1
j2 = j1 + 1



IF (del_u > tol)  THEN

    IF (del_v > tol) THEN
    		
    		
         !// Interpolation point falls between grid points in both directions:
         row1 =0.0
         row2 =0.0
         
          !// Horizontal interpolation along upper border.
          row1 = (grid_data(i1,j2) - grid_data(i1,j1) ) * del_u + grid_data(i1,j1);
          
          !// Horizontal interpolation along lower border.
          row2 = (grid_data(i2,j2) - grid_data(i2,j1) ) * del_u + grid_data(i2,j1);
          !// Vertical interpolation (downward) between row1 and row2.
          result = (row2 - row1) * del_v + row1;
            
             
     ELSE 
              
					!// Interpolation point (pretty much) falls on horizontal grid line:
					!// 1D rightward interpolation (increasing longitude).
					result = (grid_data(i1,j2) - grid_data(i1,j1)) * del_u + grid_data(i1,j1);
					 
     END IF
         
ELSE
      IF  (del_v  > tol)  THEN
                !// Interpolation point (pretty much) falls on vertical grid line:
                !// 1D downward interpolation (increasing latitude).
                result = (grid_data(i2,j1) - grid_data(i1,j1)) * del_v + grid_data(i1,j1)
                 
         ELSE
                !// Interpolation point falls on grid node exactly (ish):
                result = grid_data(i1,j1)
                 
        END IF
END IF

	  return 
END      