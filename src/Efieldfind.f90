  SUBROUTINE Efieldfind

    USE basic
    USE numerics
    USE fields
    USE tracers
    USE constants
    USE diagnostic
    USE injection
    USE interpolate
    USE random
    IMPLICIT none
    !
    ! Contains several options for the electric field used for pushing the
    ! particles, including (use parameter Efield_type): 
    ! 1) 'constant' in E_dir direction either 'r' or 'z'
    ! 2) 'sinusoid' for a nongridded sine wave in E_dir direction
    ! 3) 'sinegrid' for a gridded sine wave using Lr and Lz, numrpts and numzpts
    !    which is specified in sinegrid.f90 
    ! 4) 'GBS' which uses GBSload.f90 and has a 'static_field' T/F option
    ! 5) 'kicks'
    ! The gridded fields use lineart and bcuint to find the field at particle
    ! positions.
    ! The electric field is automatically returned in Cartesian coordinates,
    ! using the rotation routines.
    
    DOUBLE PRECISION :: ErA, ErB, EzA, EzB, cart_x, cart_z, tilt_r, tilt_z, es_pot_A, es_pot_B, &
                        rkick1=0.,rkick2=0.
    DOUBLE PRECISION, SAVE :: kick1=0.,kick2=0
    INTEGER :: k1,foundit=0

    IF (Efield_type=='kicks') THEN
       IF (rk_iter==1 .OR. ODE_solver=='boris') THEN
          IF (E_basis=='tilt') THEN
             old_Er_tilt = -1.*tracer_E_tilt(1,tracerid)
             old_Ez_tilt = -1.*tracer_E_tilt(3,tracerid)
          ELSE IF (E_basis=='cart') THEN
             old_Ex_cart  = -1.*tracer_E_cart(1,tracerid)
             old_Ez_cart  = -1.*tracer_E_cart(3,tracerid)
          END IF
       END IF
    END IF

    IF (E_basis=='tilt') THEN
       tilt_r = 0.
       tilt_z = 0.
       pos_cyl = 0.
       pos_tilt = 0.
       rot_type='pos' ! just do all the rotations now
       CALL cart_to_cyl ! it would be faster, when there are many particles, to avoid this rotation when possible
       rot_type='pos'
       CALL cyl_to_tilt
    ELSE IF (E_basis=='cart') THEN
       cart_x = 0.
       cart_z = 0.
    ELSE
       STOP 'invalid E_basis...stopping'
    END IF

    IF (periodic_field) THEN
      IF (E_basis=='tilt') THEN
         tilt_r = modulo(pos_tilt(1) - deltarmod0, deltarmod)
         tilt_z = modulo(pos_tilt(3) - deltazmod0, deltazmod)
      ELSE IF (E_basis=='cart') THEN
         cart_x = modulo(pos_cart(1) - deltarmod0, deltarmod)
         cart_z = modulo(pos_cart(3) - deltazmod0, deltazmod)
      ELSE 
         STOP 'invalid E_basis...stopping'
      END IF
   ELSE
      IF (E_basis=='tilt') THEN
         tilt_r = pos_tilt(1)
         tilt_z = pos_tilt(3)
      ELSE IF (E_basis=='cart') THEN
         cart_x = pos_cart(1)
         cart_z = pos_cart(3)
      ELSE
         STOP 'invalid E_basis...stopping'
      END IF
   END IF
   
   IF (E_basis=='tilt') THEN
      IF (Efield_type=='constant') THEN
         IF (E_dir=='r') THEN
            es_pot = E_mag1*tilt_r ! analytic form based on the constant E routine
            E_tilt(1) = E_mag1 ! include (-) later
            E_tilt(2) = 0.
            E_tilt(3) = 0.
         ELSE IF (E_dir=='z') THEN
            es_pot = E_mag1*tilt_z ! inverts E = -\nabla\phi analytically
            E_tilt(1) = 0.
            E_tilt(2) = 0.
            E_tilt(3) = E_mag1 ! include (-) later
         ELSE IF (E_dir=='rz') THEN
            es_pot = E_mag1*(tilt_r+tilt_z) ! inverts E = -\nabla\phi analytically
            E_tilt(1) = E_mag1
            E_tilt(2) = 0.
            E_tilt(3) = E_mag1 ! include (-) later
         ELSE 
            stop 'invalid E_dir...stopping'
         END IF
      ELSE IF (Efield_type=='sinusoid') THEN
         es_pot   = E_max*exp(-1.*tilt_r**4.*rscale)*cos(kr*tilt_r*sinmoder + sinmodet3*tracertime)* &
              (E_mag1*cos(kz*tilt_z*sinmodez1 + sinmodet1*tracertime) + &
              E_mag2*cos(kz*tilt_z*sinmodez2 + sinmodet2*tracertime))
         IF (energy_diag .AND. MOD(cstep, nsave_1d) == 0) THEN
            es_pot_A = E_max*exp(-1.*tilt_r**4.*rscale)*cos(kr*tilt_r*sinmoder + sinmodet3*(tracertime-nsave_1d*tracerdeltat))* &
                 (E_mag1*cos(kz*tilt_z*sinmodez1 + sinmodet1*(tracertime-nsave_1d*tracerdeltat)) + &
                 E_mag2*cos(kz*tilt_z*sinmodez2 + sinmodet2*(tracertime-nsave_1d*tracerdeltat)))
            dphidt(tracerid)=(es_pot - es_pot_A)/(nsave_1d*tracerdeltat)
         END IF
         E_tilt(1) = -1.0*E_max*exp(-1.*tilt_r**4.*rscale)*kr*sinmoder*sin(kr*tilt_r*sinmoder + sinmodet3*tracertime)* &
              (E_mag1*cos(kz*tilt_z*sinmodez1 + sinmodet1*tracertime) + &
              E_mag2*cos(kz*tilt_z*sinmodez2 + sinmodet2*tracertime))
         E_tilt(2) = 0.
         E_tilt(3) = -1.0*E_max*exp(-1.*tilt_r**4.*rscale)*kz*cos(kr*tilt_r*sinmoder + sinmodet3*tracertime)* &
              (E_mag1*sinmodez1*sin(kz*tilt_z*sinmodez1 + sinmodet1*tracertime) + &
              E_mag2*sinmodez2*sin(kz*tilt_z*sinmodez2 + sinmodet2*tracertime))
      ELSE IF (Efield_type=='kicks') THEN
         IF (cstep==1) THEN
            E_tilt(1) = E_max
            E_tilt(2) = 0.
            E_tilt(3) = E_max2
         ELSE IF (MOD(cstep,kickstep1)==0 .OR. MOD(cstep,kickstep2)==0) THEN
            IF (MOD(cstep,kickstep1)==0) THEN
               kick1 = 0.
               kick2 = 0.
               CALL gasdev_s(kick1)
               CALL gasdev_s(kick2)
               E_tilt(1) = E_max  + E_kick1a*kick1
               E_tilt(2) = 0.
               E_tilt(3) = E_max2 + E_kick1b*kick2
            ELSE IF (MOD(cstep,kickstep2)==0) THEN
               kick1 = 0.
               kick2 = 0.
               CALL gasdev_s(kick1)
               CALL gasdev_s(kick2)
               E_tilt(1) = E_max  + E_kick2a*kick1
               E_tilt(2) = 0.
               E_tilt(3) = E_max2 + E_kick2b*kick2
            END IF
         ELSE 
            E_tilt(1) = old_Er_tilt 
            E_tilt(2) = 0.
            E_tilt(3) = old_Ez_tilt 
         END IF
!!!!!!!!!!!          es_pot = (E_max + old_Er_tilt)*(tilt_r + tilt_z)+ kick1*tilt_r + kick2*tilt_z
      ELSE IF (Efield_type=='GBS' .OR. Efield_type=='sinegrid') THEN
         IF (energy_diag .AND. MOD(cstep, nsave_1d) == 0) THEN
            CALL bcuint(tilt_r,tilt_z,es_potgrid_t(:,:,tA),dpotdrtiltgrid_t(:,:,tA),dpotdztiltgrid_t(:,:,tA),&
                 d2potdrdztiltgrid_t(:,:,tA),es_pot_A,ErA,EzA)
            CALL bcuint(tilt_r,tilt_z,es_potgrid_t(:,:,tB),dpotdrtiltgrid_t(:,:,tB),dpotdztiltgrid_t(:,:,tB),&
                 d2potdrdztiltgrid_t(:,:,tB),es_pot_B,ErB,EzB)
            dphidt(tracerid)=(es_pot_B - es_pot_A)/fieldtimestep
         END IF
         CALL bcuint(tilt_r,tilt_z,es_potgrid,dpotdrtiltgrid,dpotdztiltgrid,d2potdrdztiltgrid,es_pot,E_tilt(1),E_tilt(3))
         E_tilt(2) = 0.     
      ELSE IF (Efield_type=='none') THEN
         E_tilt = 0.
      ELSE
         STOP 'invalid Efield_type...stopping'
      END IF
   ELSE IF (E_basis=='cart') THEN
      IF (Efield_type=='GBS' .OR. Efield_type=='sinegrid') THEN
         IF (energy_diag .AND. MOD(cstep, nsave_1d) == 0) THEN
            CALL bcuint(cart_x,cart_z,es_potgrid_t(:,:,tA),dpotdrtiltgrid_t(:,:,tA),dpotdztiltgrid_t(:,:,tA),&
                 d2potdrdztiltgrid_t(:,:,tA),es_pot_A,ErA,EzA)
            CALL bcuint(cart_x,cart_z,es_potgrid_t(:,:,tB),dpotdrtiltgrid_t(:,:,tB),dpotdztiltgrid_t(:,:,tB),&
                 d2potdrdztiltgrid_t(:,:,tB),es_pot_B,ErB,EzB)
            dphidt(tracerid)=(es_pot_B - es_pot_A)/fieldtimestep
          END IF
          CALL bcuint(cart_x,cart_z,es_potgrid,dpotdrtiltgrid,dpotdztiltgrid,d2potdrdztiltgrid,es_pot,E_cart(1),E_cart(3))
          E_cart(2) = 0.
       ELSE IF (Efield_type=='none') THEN
          E_cart = 0.
       ELSE IF (Efield_type=='kicks') THEN
          IF (cstep==1) THEN
             E_cart(1) = E_max
             E_cart(2) = 0.
             E_cart(3) = E_max2
          ELSE IF (MOD(cstep,kickstep1)==0 .OR. MOD(cstep,kickstep2)==0) THEN
             IF (MOD(cstep,kickstep1)==0) THEN
                kick1 = 0.
                kick2 = 0.
                CALL gasdev_s(kick1)
                CALL gasdev_s(kick2)
                E_cart(1) = E_max  + E_kick1a*kick1
                E_cart(2) = 0.
                E_cart(3) = E_max2 + E_kick1b*kick2
             ELSE IF (MOD(cstep,kickstep2)==0) THEN
                IF (readkicks) THEN
                   foundit = 0
                   DO k1 = 1,323
                      IF (kicksteps(tracerid,k1)==cstep) THEN
                         foundit = 1
                         IF (rk_iter==1) THEN
                            rkick1 = 0.
                            rkick2 = 0.
                            CALL gasdev_s(rkick1)
                            CALL gasdev_s(rkick2)
                            kick1 = rkick1
                            kick2 = rkick2
                         END IF
                         E_cart(1) = E_max  + E_kick2a*kick1
                         E_cart(2) = 0.
                         E_cart(3) = E_max2 + E_kick2b*kick2
                      ELSE IF (foundit == 0) THEN
                         E_cart(1) = old_Ex_cart
                         E_cart(2) = 0.
                         E_cart(3) = old_Ez_cart
                      END IF
                   END DO
                ELSE 
                   kick1 = 0.
                   kick2 = 0.
                   CALL gasdev_s(kick1)
                   CALL gasdev_s(kick2)
                   E_cart(1) = E_max  + E_kick2a*kick1
                   E_cart(2) = 0.
                   E_cart(3) = E_max2 + E_kick2b*kick2
                END IF
             END IF
          ELSE
             E_cart(1) = old_Ex_cart
             E_cart(2) = 0.
             E_cart(3) = old_Ez_cart
          END IF
       ELSE
          STOP 'invalid Efield_type in the cartesian basis...stopping'
       END IF
    ELSE
       STOP 'invalid E_basis...stopping'
    END IF
    
    tracer_es_pot(tracerid) = es_pot*red_mass
    
    IF (E_basis=='tilt') THEN
       E_tilt = -1.*E_tilt*red_mass ! must do after since bicubic wants y', not (-y') as input
       ! yes, this will happen each time Efieldfind is called in rk4_3d, but it is cheap and it will be correct the last time
       tracer_E_tilt(:,tracerid) = E_tilt
       E_cyl = 0.
       E_cart = 0.
       rot_type='Efield' ! flipping this switch is required before calling each rotate routine
       CALL tilt_to_cyl
       rot_type='Efield'
       CALL cyl_to_cart
    ELSE IF (E_basis=='cyl') THEN
       STOP 'no cylindrical E fields...stopping'
       E_cart = 0.
       rot_type='Efield'
       CALL cyl_to_cart
    ELSE IF (E_basis=='cart') THEN
       E_cart = -1.*E_cart*red_mass
       tracer_E_cart(:,tracerid) = E_cart
    ELSE
       STOP 'invalid E_basis...stopping'
    END IF


! next three lines must be moved

!      driftzB_tilt(tracerid) = ((tracer_vel_tilt(1,tracerid)**2+tracer_vel_tilt(3,tracerid)**2)*0.5 + tracer_vel_tilt(2,tracerid)**2)*sqrt(1.+Bratio**2)/Rcurv
!      driftzExB_tilt(tracerid) = E_tilt(1,tracerid)*tracer_pos_tilt(1,tracerid)/(Rcurv*sqrt(1.+Bratio**2))
!      driftrExB_tilt(tracerid) =  -1.*E_tilt(3,tracerid)*tracer_pos_tilt(1,tracerid)/(Rcurv*sqrt(1.+Bratio**2))

  END SUBROUTINE Efieldfind
