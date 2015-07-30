MODULE tracers
  !
  ! Declare variables for the tracer pushing.
  !
  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE  :: tracer_pos_cart, tracer_pos_cyl, tracer_pos_tilt, tracer_vel_cart, tracer_vel_cyl, tracer_vel_tilt
  INTEGER                                         :: ntracer=1,stagpart=1,stagbin=1,tracerid=1
  INTEGER, DIMENSION(:), ALLOCATABLE		  :: nturnphi ! number of turns in phi coordinate
  DOUBLE PRECISION          			  :: tracerdtfactor=1.0,tracerdeltat,tracertime=0.,red_mass=1.
  DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE  :: phicheck, x_in , y_in, z_in ! initial position
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: pos_cart,pos_cyl,pos_tilt,E_cart,E_cyl,E_tilt,B_cart,vel_cart,vel_cyl,vel_tilt

END MODULE tracers
