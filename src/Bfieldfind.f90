  SUBROUTINE Bfieldfind
    !    
    USE basic
    USE numerics
    USE tracers
    USE fields
    !
    IMPLICIT none
    !
    ! Find the magnetic field for several case:
    ! 1) 'torpex'
    ! 2) 'torus' which is 'torpex' without spatial dependence
    ! 3) 'slab' has the B field in the x direction with magnitude Bnaught
    !
    DOUBLE PRECISION                                      :: A=0.
    DOUBLE PRECISION, DIMENSION(ndim_tracer)              :: B

    B_cart = 0.
    B = 0.

    IF (Bfield_type=='torpex') THEN
      A    = Rcurv/(pos_cart(1)**2+pos_cart(2)**2) ! find ratio of curvature radius to r^2
      B(1) = -1.0*A*pos_cart(2)*B_t
      B(2) =      A*pos_cart(1)*B_t
      B(3) = Bratio
    ELSE IF (Bfield_type=='torus') THEN
      B(1) = -1.0*B_t
      B(2) =  1.0*B_t
      B(3) = Bratio
    ELSE IF (Bfield_type=='slab') THEN
      B(1) = 0.
      B(2) = By
      B(3) = 0.
    ELSE
      B = 0.
    END IF

    B_cart = B*red_mass

  END SUBROUTINE Bfieldfind
