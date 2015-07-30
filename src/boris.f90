  SUBROUTINE boris
    !    
    USE basic
    USE numerics
    USE tracers
    USE fields
    !
    IMPLICIT none
    !
    ! Moves a set of tracers in the Cartesian plane.  The tracers obey the entire Lorentz
    ! force law, which is implemented through the Boris algorithm (Birdsall/Langdon book). 
    ! The technique gives the tracer coordinates in phase space, including all three
    ! velocities and positions.
    !
    ! When the new value of the tracer attributes are 
    ! returned to the main program, the new values must be set up to be passed into the
    ! Boris integrator at the next step.

    DOUBLE PRECISION                                 :: A=0., t2=0., deltat_t2=0., tdt05=0.
    DOUBLE PRECISION, DIMENSION(ndim_tracer)         :: vminus,vprime,vplus,t,s,B

    B_cart=0.
    B=0.
    vminus=0.
    vprime=0.
    vplus=0.
    t=0.
    s=0.

    tdt05 = tracerdeltat*0.5

    pos_cart = tracer_pos_cart(:,tracerid)
    vel_cart = tracer_vel_cart(:,tracerid)
    
    CALL Efieldfind
    
    vminus = vel_cart + E_cart*tdt05 ! find v-minus, equation 7 in B&L book, pg 58 (Boris 1970)
    ! when the B field is parallel to the z axis...motion due to B is just rotation
    
    CALL Bfieldfind
    B = B_cart
    
    t(1) = B(1)*tdt05  ! definition of the t vector from Boris
    t(2) = B(2)*tdt05  ! factor of A comes from cylindrial to cartesian
    t(3) = B(3)*tdt05     
    t2 = t(1)**2. + t(2)**2. + t(3)**2.
    deltat_t2 = tracerdeltat/(1.+t2) ! for quicker simplification of s
    s(1) = B(1)*deltat_t2 
    s(2) = B(2)*deltat_t2
    s(3) = B(3)*deltat_t2
    
    vprime(1) = vminus(1) + t(3)*vminus(2) - t(2)*vminus(3) ! the next six lines are applicable when
    vprime(2) = vminus(2) + t(1)*vminus(3) - t(3)*vminus(1) ! B and v are in arbitrary directions ...
    vprime(3) = vminus(3) + t(2)*vminus(1) - t(1)*vminus(2) ! ... which is necessary because of the toroidal field
    
    vplus(1) = vminus(1) + s(3)*vprime(2) - s(2)*vprime(3)
    vplus(2) = vminus(2) + s(1)*vprime(3) - s(3)*vprime(1)
    vplus(3) = vminus(3) + s(2)*vprime(1) - s(1)*vprime(2)
    
    ! find the updated velocity, which will be accessible to the main program
    tracer_vel_cart(:,tracerid) = vplus + E_cart*tdt05  ! now this is the velocity to output after the rotation
    
    ! find the updated position, which will be accessible to the main program
    tracer_pos_cart(:,tracerid) = pos_cart + tracerdeltat*tracer_vel_cart(:,tracerid) 
    
  END SUBROUTINE boris
