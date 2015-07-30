  SUBROUTINE cyl_to_cart
    !
    USE fields
    USE tracers
    USE diagnostic
    USE injection
    IMPLICIT none
    !
    !
    ! Rotates from a nontilted cylindrical coordinate system to a 
    ! Cartesian basis.
    !	
    
    DOUBLE PRECISION :: rpos=0., rinv=0.

    IF (rot_type=='pos') THEN
       write(*,*) 'pos_rotate cyl_to_cart is not in service...stopping'
       ! 	pos_cart(1) = pos_cyl(1)*cos(pos_cyl(2))
       ! 	pos_cart(2) = pos_cyl(1)*sin(pos_cyl(2))
       ! 	pos_cart(3) = pos_cyl(3)
       ! 	vel_cart(1) = vel_cyl(1)*pos_cart(1)*rinv - vel_cyl(2)*pos_cart(2)*rinv
       ! 	vel_cart(2) = vel_cyl(1)*pos_cart(2)*rinv + vel_cyl(2)*pos_cart(1)*rinv
       ! 	vel_cart(3) = vel_cyl(3)
    ELSE IF (rot_type=='Efield') THEN
       rpos = sqrt(pos_cart(1)**2 + pos_cart(2)**2)
       IF (rpos/=0) THEN 
          rinv = 1./rpos
       ELSE
          write(*,*) 'rpos = 0 LIKELY ERROR!...continuing'
          rinv = 100000000.  ! this is a bit flippant but shouldn't ever be necessary
       END IF
       E_cart(1) = E_cyl(1)*pos_cart(1)*rinv - E_cyl(2)*pos_cart(2)*rinv
       E_cart(2) = E_cyl(1)*pos_cart(2)*rinv + E_cyl(2)*pos_cart(1)*rinv
       E_cart(3) = E_cyl(3)
    END IF
    
  END SUBROUTINE cyl_to_cart
