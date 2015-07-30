SUBROUTINE initial
  !
  !   Set initial conditions
  !
  !
  USE constants
  USE basic
  USE injection
  USE tracers
  USE fields
  USE random
  USE numerics
  IMPLICIT none
  INCLUDE 'mpif.h'
  !
  !
  ! Specifies the initial positions and velocity vectors for tracer particles.  The user
  ! should give init_type and kxinit (if necessary), ngp, mu, sigma, xstart and ystart.
  !
  
  DOUBLE PRECISION :: harvest1=0.,harvest1a=0.,harvest1b=0.,harvest2=0.,harvest3=0., &
                      deltax,deltay,deltaz,alpha,beta,Vmag,Vmaga,Vmagb

  DOUBLE PRECISION, DIMENSION(:,:),  ALLOCATABLE  :: tracer_pos_cart_tmp, tracer_pos_cyl_tmp, tracer_vel_cart_tmp
  DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE  :: phicheck_tmp

  INTEGER :: myid, master, numprocs, processor,ierr,stat(MPI_STATUS_SIZE)
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  master = 0
  !________________________________________________________________________________
  WRITE(*,'(a/)') '=== Set initial conditions ==='
  !________________________________________________________________________________
  !

  tracer_es_pot = 0.
  tracer_E_cart = 0.
  tracer_E_cyl  = 0.
  tracer_E_tilt = 0.
  
  phicheck = 0.
  nturnphi = 0

  deltax = Xwidth_injection/ntracer
  deltay = Ywidth_injection/ntracer
  deltaz = Zwidth_injection/ntracer
  
  IF (rand_seed) THEN
      CALL init_random_seed()
  ELSE
      CALL init_default_seed()
  END IF

  
  IF (inject_coord_cart) THEN
    IF (inject_type_pos=='pointlike') THEN
	DO tracerid = 1,ntracer
	tracer_pos_cart(1,tracerid)= XCM_injection
	tracer_pos_cart(2,tracerid)= YCM_injection
	! for comparison with the next value of phi to see if a turn has occurred\
	phicheck(tracerid) = atan2(tracer_pos_cart(2,tracerid),tracer_pos_cart(1,tracerid))
	IF (phicheck(tracerid)<0) THEN
	  phicheck(tracerid) = phicheck(tracerid) + 2*pi
	END IF
	tracer_pos_cart(3,tracerid)= ZCM_injection
	END DO
    ELSE IF (inject_type_pos=='gaussian') THEN

	IF (myid == master) THEN

          ALLOCATE(tracer_pos_cart_tmp(ndim_tracer,ntracer)) 
          ALLOCATE(phicheck_tmp(ntracer))
          DO tracerid = 1,ntracer
	  harvest1 = 0.
	  harvest2 = 0.
	  harvest3 = 0.
	  CALL gasdev_s(harvest1)
	  tracer_pos_cart_tmp(1,tracerid) = XCM_injection + Xwidth_injection*harvest1
	  CALL gasdev_s(harvest2)
	  tracer_pos_cart_tmp(2,tracerid) = YCM_injection + Ywidth_injection*harvest2
	  phicheck_tmp(tracerid) = atan2(tracer_pos_cart_tmp(2,tracerid),tracer_pos_cart_tmp(1,tracerid))
	  IF (phicheck_tmp(tracerid)<0) THEN
	      phicheck_tmp(tracerid) = phicheck_tmp(tracerid) + 2*pi
	  END IF
	  CALL gasdev_s(harvest3)
	  tracer_pos_cart_tmp(3,tracerid) = ZCM_injection + Zwidth_injection*harvest3
	  END DO


	DO processor = master+1,numprocs-1
	  DO tracerid = 1,ntracer
	  harvest1 = 0.
	  harvest2 = 0.
	  harvest3 = 0.
	  CALL gasdev_s(harvest1)
	  tracer_pos_cart(1,tracerid) = XCM_injection + Xwidth_injection*harvest1
	  CALL gasdev_s(harvest2)
	  tracer_pos_cart(2,tracerid) = YCM_injection + Ywidth_injection*harvest2
	  phicheck(tracerid) = atan2(tracer_pos_cart(2,tracerid),tracer_pos_cart(1,tracerid))
	  IF (phicheck(tracerid)<0) THEN
	      phicheck(tracerid) = phicheck(tracerid) + 2*pi
	  END IF
	  CALL gasdev_s(harvest3)
	  tracer_pos_cart(3,tracerid) = ZCM_injection + Zwidth_injection*harvest3
	  END DO
	  IF (processor .NE. master) THEN
		CALL MPI_SEND(tracer_pos_cart,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
		CALL MPI_SEND(phicheck,ntracer,MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
	  END IF
	END DO
        tracer_pos_cart = tracer_pos_cart_tmp
        phicheck = phicheck_tmp
        DEALLOCATE(tracer_pos_cart_tmp) 
        DEALLOCATE(phicheck_tmp)


	END IF
        IF (myid .NE. master) THEN
	   CALL MPI_RECV(tracer_pos_cart,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
	   CALL MPI_RECV(phicheck,ntracer,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
        END IF

    ELSE IF (inject_type_pos=='uniform') THEN
	DO tracerid = 1,ntracer
	tracer_pos_cart(1,tracerid) = XCM_injection - Xwidth_injection*0.5 + deltax*tracerid
	tracer_pos_cart(2,tracerid) = YCM_injection - Ywidth_injection*0.5 + deltay*tracerid
	tracer_pos_cart(3,tracerid) = ZCM_injection - Zwidth_injection*0.5 + deltaz*tracerid
	END DO
    ELSE
	STOP 'position injection type not available'
    END IF


    IF (E_basis=='tilt') THEN
	DO tracerid = 1,ntracer
	pos_cart=tracer_pos_cart(:,tracerid)
	rot_type='pos'    ! initial rotation position for finding the field
	CALL cart_to_cyl
	tracer_pos_cyl(:,tracerid) = pos_cyl
	rot_type='pos'
	CALL cyl_to_tilt
	tracer_pos_tilt(:,tracerid) = pos_tilt
	END DO
    END IF

  ELSE IF (inject_coord_cyl) THEN
    STOP 'cylindrical coordinate injection in not well-defined...stopping'

    ! Options for injection position of tracers
    IF (inject_type_pos=='pointlike') THEN
	DO tracerid = 1,ntracer
	tracer_pos_cyl(1,tracerid)= RCM_injection
	tracer_pos_cyl(2,tracerid)= 0.
	phicheck(tracerid) = tracer_pos_cyl(2,tracerid)
	tracer_pos_cyl(3,tracerid)= ZCM_injection
	END DO
    ELSE IF (inject_type_pos=='gaussian') THEN

	IF (myid == master) THEN

        ALLOCATE(phicheck_tmp(ntracer)) 
        ALLOCATE(tracer_pos_cyl_tmp(ndim_tracer,ntracer))

	DO tracerid = 1,ntracer
	  harvest1 = 0.
	  harvest2 = 0.
	  CALL gasdev_s(harvest1)
	  tracer_pos_cyl_tmp(1,tracerid) = RCM_injection + Rwidth_injection*harvest1
	  CALL gasdev_s(harvest2)
	  tracer_pos_cyl_tmp(2,tracerid) = 0.
	  phicheck_tmp(tracerid) = tracer_pos_cyl_tmp(2,tracerid)
	  tracer_pos_cyl_tmp(3,tracerid) = ZCM_injection + Zwidth_injection*harvest2
	END DO


	DO processor = master+1,numprocs-1
	  DO tracerid = 1,ntracer
	  harvest1 = 0.
	  harvest2 = 0.
	  CALL gasdev_s(harvest1)
	  tracer_pos_cyl(1,tracerid) = RCM_injection + Rwidth_injection*harvest1
	  CALL gasdev_s(harvest2)
	  tracer_pos_cyl(2,tracerid) = 0.
	  phicheck(tracerid) = tracer_pos_cyl(2,tracerid)
	  tracer_pos_cyl(3,tracerid) = ZCM_injection + Zwidth_injection*harvest2
          END DO
	  IF (processor .NE. master) THEN
		CALL MPI_SEND(tracer_pos_cyl,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
		CALL MPI_SEND(phicheck,ntracer,MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
	  END IF
	END DO

        tracer_pos_cyl = tracer_pos_cyl_tmp
        phicheck = phicheck_tmp
        DEALLOCATE(phicheck_tmp) 
        DEALLOCATE(tracer_pos_cyl_tmp)

	END IF
        IF (myid .NE. master) THEN
	   CALL MPI_RECV(tracer_pos_cyl,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
	   CALL MPI_RECV(phicheck,ntracer,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
        END IF


    ELSE IF (inject_type_pos=='sinweight') THEN ! sin curve weighting for phase-mixing study
	
	STOP 'weights for phase-mixing are not available'
	!           tracer_pos_cyl(1,tracerid) = 
	!           tracer_pos_cyl(2,tracerid) = 
	! 9 Feb 2010 KBG
	! wgt(tracerid) = sin(kxinit*pi*x0(tracerid)/Lx)
	
    ELSE IF (inject_type_pos=='fromfile') THEN ! initial positions from file
	STOP 'initial positions from file is not available now...stopping' ! 9 Feb 2010 KBG
    ELSE
	STOP 'invalid initial tracer condition...stopping'
    END IF


    IF (E_basis=='tilt') THEN
	DO tracerid = 1,ntracer
	pos_cyl = tracer_pos_cyl(:,tracerid)
	rot_type='pos'
	CALL cyl_to_cart
	tracer_pos_cart(:,tracerid) = pos_cart
	rot_type='pos'
	CALL cyl_to_tilt
	tracer_pos_tilt(:,tracerid) = pos_tilt
	END DO
    END IF
  ELSE IF (inject_coord_tilt) THEN
    STOP 'tilted coordinate injection position not available now...stopping'
  ELSE
    STOP 'choice of injection coordinate basis is invalid...stopping'
  END IF

  x_in = tracer_pos_cart(1,:)
  y_in = tracer_pos_cart(2,:)
  z_in = tracer_pos_cart(3,:)
  
  ! Options for injection velocity vectors of tracers
  IF (inject_type_vel=='stationary') THEN
     DO tracerid = 1,ntracer 
        tracer_vel_cart(:,tracerid)= 0.
     END DO
  ELSE IF (inject_type_vel=='alphabeta') THEN
     IF (V0 == 0 .AND. devV == 0.) THEN
        STOP 'velocity cone cannot be vanishing'
     END IF
     tracerid = 0

     IF (myid == master) THEN
	ALLOCATE(tracer_vel_cart_tmp(ndim_tracer,ntracer)) 
	tracerid = 0
	DO WHILE (tracerid < ntracer)
	    harvest1a = 0.
	    harvest1b = 0.
	    harvest2 = 0.
	    harvest3 = 0.
	    Vmag = 0.
	    Vmaga = 0.
	    Vmagb = 0.
	    alpha = 0.
	    beta = 0.
	    CALL gasdev_s(harvest1a)
	    Vmaga = V0 + devV*harvest1a
	    !CALL gasdev_s(harvest1b)
	    !Vmagb = V0 + devV*harvest1b
	    !Vmag = sqrt((Vmaga**2 + Vmagb**2)*0.5)
	    IF (Vmaga > 0.) THEN
	      tracerid = tracerid + 1
	      CALL gasdev_s(harvest2)
	      alpha = alpha0 + devalpha*harvest2
	      CALL gasdev_s(harvest3)
	      beta = beta0 + devbeta*harvest3
	      tracer_vel_cart_tmp(1,tracerid) = Vmaga*cos(alpha)*sin(beta)
	      tracer_vel_cart_tmp(2,tracerid) = Vmaga*cos(alpha)*cos(beta)
	      tracer_vel_cart_tmp(3,tracerid) = Vmaga*sin(alpha)
	    END IF
	END DO


	DO processor = master+1,numprocs-1
	tracerid = 0
	DO WHILE (tracerid < ntracer)
	    harvest1a = 0.
	    harvest1b = 0.
	    harvest2 = 0.
	    harvest3 = 0.
	    Vmag = 0.
	    Vmaga = 0.
	    Vmagb = 0.
	    alpha = 0.
	    beta = 0.
	    CALL gasdev_s(harvest1a)
	    Vmaga = V0 + devV*harvest1a
	    !CALL gasdev_s(harvest1b)
	    !Vmagb = V0 + devV*harvest1b
	    !Vmag = sqrt((Vmaga**2 + Vmagb**2)*0.5)
	    IF (Vmaga > 0.) THEN
	      tracerid = tracerid + 1
	      CALL gasdev_s(harvest2)
	      alpha = alpha0 + devalpha*harvest2
	      CALL gasdev_s(harvest3)
	      beta = beta0 + devbeta*harvest3
	      tracer_vel_cart(1,tracerid) = Vmaga*cos(alpha)*sin(beta)
	      tracer_vel_cart(2,tracerid) = Vmaga*cos(alpha)*cos(beta)
	      tracer_vel_cart(3,tracerid) = Vmaga*sin(alpha)
	    END IF
	END DO
	IF (processor .NE. master) THEN
	    CALL MPI_SEND(tracer_vel_cart,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
	END IF
        END DO

        tracer_vel_cart = tracer_vel_cart_tmp
        DEALLOCATE(tracer_vel_cart_tmp) 


     END IF
     IF (myid .NE. master) THEN
	CALL MPI_RECV(tracer_vel_cart,ndim_tracer*ntracer,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
     END IF

  ELSE
     STOP 'invalid initial tracer condition...stopping'
  END IF
  
  !     IF (ODE_solver=='gyrav' .AND. ngp .gt. 0) THEN
  !        harvest1 = 0.
  !        harvest2 = 0.
  !        DO tracerid = 1,ntracer
  !           call gasdev_s(harvest1)
  !           harvest1 = mu + sigma*harvest1
  !           call gasdev_s(harvest2)
  !           harvest2 = mu + sigma*harvest2
  !           gyrorad(tracerid) = sqrt(harvest1**2 + harvest2**2) ! gyroradius: 2D maxwellian
  !        END DO
  !     END IF
  
END SUBROUTINE initial
  
 

