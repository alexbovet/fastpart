  SUBROUTINE tilt_to_cyl
    !
    USE fields
    USE tracers
    USE diagnostic
    USE injection
    IMPLICIT none
    !
    ! Rotates from a tilted (magnetic plane) cylindrical basis to a cylindrical basis.
    !	

    IF (rot_type=='pos') THEN
       stop 'position rotate from tilt to cyl is not set up yet...stopping'
       pos_cyl(1) = 0.
       pos_cyl(2) = 0.
       pos_cyl(3) = 0.
    ELSE IF (rot_type=='Efield') THEN
       E_cyl(1) = E_tilt(1)
       ! Rcurv/r is consistent with the theta = theta(Rcurv) approximation
       E_cyl(2) = -1.*Rcurv*E_tilt(3)*sineoftheta/pos_tilt(1) ! terms missing because the parallel component in...
       E_cyl(3) = E_tilt(3)*cosineoftheta ! ...the magnetic plane is zero by definition
    END IF
    
  END SUBROUTINE tilt_to_cyl
