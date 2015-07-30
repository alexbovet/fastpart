  SUBROUTINE cart_to_cyl

    USE constants
    USE fields
    USE tracers
    USE diagnostic
    USE injection
    IMPLICIT none
    !
    !
    ! Rotates from a Cartesian basis to a cylindrical basis.
    !	
    DOUBLE PRECISION :: rpos=0., rinv=0., Rcurv_inv = 0., phidiff=0., phidiffcheck=0., rawphi=0.
    
    phidiffcheck = pi

    IF (rot_type=='pos') THEN
      rpos = sqrt(pos_cart(1)**2 + pos_cart(2)**2)
      IF (rpos/=0) THEN
	rinv = 1./rpos
      ELSE
        write(*,*) 'rpos = 0 LIKELY ERROR!...continuing'
        rinv = 100000000.  ! this is a bit flippant but shouldn't ever be necessary
      END IF
      pos_cyl(1) = rpos
      rawphi = atan2(pos_cart(2),pos_cart(1))
      IF (rawphi<0) THEN
        rawphi = rawphi + 2*pi
      END IF
      phidiff = phicheck(tracerid) - rawphi
      IF (phidiff>phidiffcheck) THEN
	nturnphi(tracerid) =nturnphi(tracerid) + 1
      ELSE IF (phidiff<-1.*phidiffcheck) THEN
	nturnphi(tracerid) =nturnphi(tracerid) - 1
      ELSE
        nturnphi(tracerid) =nturnphi(tracerid)
      END IF
      phicheck(tracerid) = rawphi
      pos_cyl(2) = rawphi + 2*pi*nturnphi(tracerid)
      pos_cyl(3) = pos_cart(3)
      vel_cyl(1) = pos_cart(1)*rinv*vel_cart(1) + pos_cart(2)*rinv*vel_cart(2)
      vel_cyl(2) = -1.0*pos_cart(2)*rinv*vel_cart(1) + pos_cart(1)*rinv*vel_cart(2)
      !vel_cyl(2) = vel_cyl(2)*Rcurv_inv ! finding the linear, rather than angular velocity...is this appropriate?
      vel_cyl(3) = vel_cart(3)
    ELSE IF (rot_type=='Efield') THEN
      stop 'E_rotate cartesian to cylinder not set up'
      E_cyl(1) = 0.
      E_cyl(2) = 0.
      E_cyl(3) = 0.
    END IF
      
  END SUBROUTINE cart_to_cyl
