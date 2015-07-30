  SUBROUTINE cyl_to_tilt

    USE constants
    USE fields
    USE tracers
    USE diagnostic
    USE injection
    !
    IMPLICIT none
    !
    !
    ! Rotates from a nontilted cylindrical coordinate system to a 
    ! tilted magnetic plane.
    !	

    DOUBLE PRECISION :: Rcurv_inv
  
    IF (rot_type=='pos') THEN
      ! could be some mistakes with Rcurv factors here: KBG 17 Feb 2010
      pos_tilt(1) = 	   pos_cyl(1)
      ! must use r\phi instead of just \phi - using the r = R_c approximation for consistency with Bratio not depdnt on r
      pos_tilt(2) =     Rcurv*pos_cyl(2)*cosineoftheta + pos_cyl(3)*sineoftheta 
      pos_tilt(3) = -1.*Rcurv*pos_cyl(2)*sineoftheta   + pos_cyl(3)*cosineoftheta
      
      vel_tilt(1) =      vel_cyl(1)
      vel_tilt(2) = 	   vel_cyl(2)*cosineoftheta + vel_cyl(3)*sineoftheta 
      vel_tilt(3) =  -1.*vel_cyl(2)*sineoftheta   + vel_cyl(3)*cosineoftheta
    ELSE IF (rot_type=='Efield') THEN
      stop 'E_rotate cylinder to tilted cylinder is not set up...stopping'
      E_tilt(1) = 0.
      E_tilt(2) = 0.
      E_tilt(3) = 0.
    END IF

  END SUBROUTINE cyl_to_tilt
