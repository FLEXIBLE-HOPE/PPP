! I/O PARAMETERS:                                                       
!                                                                       
!   NAME    I/O  A/S   DESCRIPTION OF PARAMETERS                        
!   ------  ---  ---   ------------------------------------------------                 
!   n     	 O    A    the normal of this surface
!		flux_dir I				 the unit direction of radiation flux
!   flux     I    S    the radiation flux, must contain shadow factor, i.e. flux*shadow_factor
!   area   	 I    S    the area of the surface
!		u				 I				 the specularity of the material
!		v				 I				 the reflectivity of the material
        
subroutine radiationForce(n,flux_dir,flux, area, u, v, force)

implicit none


double precision c,area,flux, u,v,tmp
double precision, dimension(3) :: flux_dir
double precision, dimension(3) :: n
double precision, dimension(3) :: force
    
double precision dp,cos_theta,W
double precision, dimension(3) :: h,reflection_direction,f1,f2,f3
    
    
  c=299792458.0
  force(1) = 0
  force(2) = 0
  force(3) = 0

  ! dot product
  cos_theta = -( flux_dir(1)*n(1) + flux_dir(2)*n(2) + flux_dir(3)*n(3))
  W = flux*cos_theta
  h(1) = cos_theta*n(1)
  h(2) = cos_theta*n(2)
  h(3) = cos_theta*n(3)

  reflection_direction(1) = 2.0*h(1) + flux_dir(1)
  reflection_direction(2) = 2.0*h(2) + flux_dir(2)
  reflection_direction(3) = 2.0*h(3) + flux_dir(3)
  tmp = dsqrt(reflection_direction(1)**2 + reflection_direction(2)**2 + reflection_direction(3)**2)
  reflection_direction(1) = reflection_direction(1)/tmp
  reflection_direction(2) = reflection_direction(2)/tmp
  reflection_direction(3) = reflection_direction(3)/tmp

  f1(1) = area*W/c*flux_dir(1)
  f1(2) = area*W/c*flux_dir(2)
  f1(3) = area*W/c*flux_dir(3)

  f2(1) = -area*W/c*u*v*reflection_direction(1)
  f2(2) = -area*W/c*u*v*reflection_direction(2)
  f2(3) = -area*W/c*u*v*reflection_direction(3)

  f3(1) = -2.0/3.0*area*W*cos_theta*v*(1.0-u)/c*n(1)
  f3(2) = -2.0/3.0*area*W*cos_theta*v*(1.0-u)/c*n(2)
  f3(3) = -2.0/3.0*area*W*cos_theta*v*(1.0-u)/c*n(3)

  force(1) = f1(1) + f2(1) + f3(1)
  force(2) = f1(2) + f2(2) + f3(2)
  force(3) = f1(3) + f2(3) + f3(3)

  return
end
