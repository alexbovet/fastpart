MODULE injection
!
!   Module for injection parameters
!
!
  DOUBLE PRECISION    :: RCM_injection=0., ZCM_injection=0., Rwidth_injection=0., Zwidth_injection=0., &
			 XCM_injection=0., YCM_injection=0., Xwidth_injection=0., Ywidth_injection=0., &
		         alpha0 = 0., beta0 = 0., V0 = 0., devalpha = 0., devbeta = 0., devV = 0. 
  CHARACTER(LEN=10)   :: inject_type_pos='pointlike', inject_type_vel='stationary'
  LOGICAL	      :: inject_coord_cart=.FALSE., inject_coord_cyl=.FALSE., inject_coord_tilt=.FALSE., rand_seed = .TRUE.
  
!  
END MODULE injection

