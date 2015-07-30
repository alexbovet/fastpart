SUBROUTINE rk4_3d
  !
  USE basic
  USE numerics
  USE tracers
  USE fields
  !
  IMPLICIT NONE
  !
  ! Runge-Kutta 4th order method for position and velocity in 
  ! 3 dimensions.
  !

  DOUBLE PRECISION :: tt,es_pot_B
  DOUBLE PRECISION, DIMENSION(ndim_tracer) :: accel, pos0,pos1,pos2,pos3, vel0,vel1,vel2,vel3, &
                                              dpos_1,dpos_2,dpos_3,dpos_4, dvel_1,dvel_2,dvel_3,dvel_4
  tt=tracerdeltat*0.5

  B_cart=0.
  E_cart=0.
  accel=0.

  pos0 = tracer_pos_cart(:,tracerid)
  vel0 = tracer_vel_cart(:,tracerid)

  pos_cart = pos0
  vel_cart = vel0

  rk_iter = 1

  ! find the first Runge-Kutta term
  CALL Bfieldfind ! uses pos_cart
  CALL Efieldfind ! uses pos_cart
  accel(1) = E_cart(1) + vel_cart(2)*B_cart(3) - vel_cart(3)*B_cart(2)
  accel(2) = E_cart(2) + vel_cart(3)*B_cart(1) - vel_cart(1)*B_cart(3)
  accel(3) = E_cart(3) + vel_cart(1)*B_cart(2) - vel_cart(2)*B_cart(1)
  
  dpos_1 = vel0
  dvel_1 = accel
  
  ! find the position for second RK4 term
  pos1 = pos0 + dpos_1*tt

  ! find the velocity for second RK4 term  
  vel1 = vel0 + dvel_1*tt

  ! reset for next RK4 evaluation
  pos_cart = pos1  
  B_cart=0.  
  E_cart=0.
  accel=0.

  rk_iter = 2

  ! note: do not repeat the time interpolation...
  ! ...assume static field during a tracer time step...
  ! ...this is good to less than 1% for a GBS field KBG

  ! find the second Runge-Kutta term
  CALL Bfieldfind ! uses pos_cart
  CALL Efieldfind ! uses pos_cart 
  accel(1) = E_cart(1) + vel1(2)*B_cart(3) - vel1(3)*B_cart(2)
  accel(2) = E_cart(2) + vel1(3)*B_cart(1) - vel1(1)*B_cart(3)
  accel(3) = E_cart(3) + vel1(1)*B_cart(2) - vel1(2)*B_cart(1)
  
  dpos_2 = vel1
  dvel_2 = accel
  
  ! find the position for the third Runge-Kutta term
  pos2 = pos0 + dpos_2*tt

  ! find the velocity for the third RK4 term
  vel2 = vel0 + dvel_2*tt

  ! reset for next field evaluation
  pos_cart = pos2
  B_cart=0.
  E_cart=0.
  accel=0.

  rk_iter = 3

  ! find the third Runge-Kutta term
  CALL Bfieldfind ! uses pos_cart
  CALL Efieldfind ! uses pos_cart
  accel(1) = E_cart(1) + vel2(2)*B_cart(3) - vel2(3)*B_cart(2)
  accel(2) = E_cart(2) + vel2(3)*B_cart(1) - vel2(1)*B_cart(3)
  accel(3) = E_cart(3) + vel2(1)*B_cart(2) - vel2(2)*B_cart(1)

  dpos_3 = vel2
  dvel_3 = accel
 
  ! find position for evaluating field for fourth RK4 accel term  
  pos3 = pos0 + dpos_3*tracerdeltat

  ! find velocity for fourth RK4 position term
  vel3 = vel0 + dvel_3*tracerdeltat

  ! again, assume that the field is stationary across a 
  ! tracer time step - this will not be good for all fields

  ! reset for next field evaluation
  pos_cart = pos3
  B_cart=0.
  E_cart=0.
  accel=0.

  rk_iter = 4

  ! find the fourth Runge-Kutta term
  CALL Bfieldfind ! uses pos_cart
  CALL Efieldfind ! uses pos_cart
  accel(1) = E_cart(1) + vel3(2)*B_cart(3) - vel3(3)*B_cart(2)
  accel(2) = E_cart(2) + vel3(3)*B_cart(1) - vel3(1)*B_cart(3)
  accel(3) = E_cart(3) + vel3(1)*B_cart(2) - vel3(2)*B_cart(1)
  
  dpos_4 = vel3
  dvel_4 = accel

  ! combine the guesses into the final answer for this particle, for this time step

  tracer_pos_cart(:,tracerid) = &
       pos0 + tracerdeltat*(dpos_1/6. + dpos_2/3. + dpos_3/3. + dpos_4/6.)
  
  tracer_vel_cart(:,tracerid) = &
       vel0 + tracerdeltat*(dvel_1/6. + dvel_2/3. + dvel_3/3. + dvel_4/6.)
  
!  tracer_pos_cart(2,tracerid) = &
!       pos0(2) + tracerdeltat*(dpos_1(2)/6. + dpos_2(2)/3. + dpos_3(2)/3. + dpos_4(2)/6.)
!  
!  tracer_vel_cart(2,tracerid) = &
!       vel0(2) + tracerdeltat*(dvel_1(2)/6. + dvel_2(2)/3. + dvel_3(2)/3. + dvel_4(2)/6.)
!  
!  tracer_pos_cart(3,tracerid) = &
!       pos0(3) + tracerdeltat*(dpos_1(3)/6. + dpos_2(3)/3. + dpos_3(3)/3. + dpos_4(3)/6.)
!  
!  tracer_vel_cart(3,tracerid) = &
!       vel0(3) + tracerdeltat*(dvel_1(3)/6. + dvel_2(3)/3. + dvel_3(3)/3. + dvel_4(3)/6.)

END SUBROUTINE rk4_3d
