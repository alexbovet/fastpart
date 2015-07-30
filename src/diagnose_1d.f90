SUBROUTINE diagnose_1d

  USE basic
  USE futils
  USE diagnostic
  USE tracers
  USE injection
  USE fields
  USE numerics
  IMPLICIT none
  
  INCLUDE 'mpif.h'

  INTEGER                                        :: iframe
  CHARACTER(64)                                  :: dset_name

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tmp_part_KE , tmp_part_PE, tmp_dphidt, tmp_driftrExB_tilt, tmp_driftzB_tilt, tmp_driftzExB_tilt
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tmp_tracer_pos_cart, tmp_tracer_pos_cyl, tmp_tracer_pos_tilt, tmp_tracer_vel_cart, tmp_tracer_vel_cyl, tmp_tracer_vel_tilt

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: tmp_tracer_E_cyl, tmp_tracer_E_cart, tmp_tracer_E_tilt

  INTEGER :: myid, master, numprocs, ierr
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  master = 0

!________________________________________________________________________________
!
  IF ( myid == master ) THEN

       WRITE(*,'(a,1x,i10.10,a1,i10.10,20x,a,1pe10.2)') &
       &    '*** Timestep (this run/total) =', step, '/', cstep, 'Time =', time

  END IF
!________________________________________________________________________________
!

  IF (radialbc == 'recycle') THEN
    DO tracerid=1,ntracer
      IF (tracer_pos_tilt(1,tracerid) < rmin .OR. tracer_pos_tilt(1,tracerid) > rmax) THEN
        tracer_pos_cart(1,tracerid) = x_in(tracerid) ! reset to the initial position to keep providing statistics
        tracer_pos_cart(2,tracerid) = y_in(tracerid)
        tracer_pos_cart(3,tracerid) = z_in(tracerid)
      END IF
    END DO
  ELSE
  	STOP 'invalid radialbc...stopping'
  END IF

  !doing the rotations here is redundant but not very expensive when nsave_1d is large
  IF (E_basis=='tilt') THEN
     DO tracerid=1,ntracer
        pos_cart = tracer_pos_cart(:,tracerid)
        vel_cart = tracer_vel_cart(:,tracerid)
        rot_type='pos' ! the rotation is also done in Efieldfind, and this one simply updates the position for the diagnostic
        CALL cart_to_cyl  
        tracer_pos_cyl(:,tracerid) = pos_cyl
        tracer_vel_cyl(:,tracerid) = vel_cyl
        rot_type='pos'
        CALL cyl_to_tilt
        tracer_pos_tilt(:,tracerid) = pos_tilt
        tracer_vel_tilt(:,tracerid) = vel_tilt
     END DO
  END IF

  IF (myid == master) THEN

	CALL append(fidres, "/data/var1d/time", time)
	     
	IF (cstep==0) THEN 
	    iframe=1
	ELSE
	    CALL getatt(fidres,"/data/var1d/" , "frames", iframe) 
	    iframe=iframe+1

	END IF 
	  
	CALL attach(fidres,"/data/var1d/" , "frames", iframe) 

  END IF

  IF (energy_diag) THEN

      
    part_KE = 0.
    part_PE = 0.
    
    DO tracerid=1,ntracer
      part_KE(tracerid) = (tracer_vel_cart(1,tracerid)**2 + tracer_vel_cart(2,tracerid)**2 + tracer_vel_cart(3,tracerid)**2)*0.5
      part_PE(tracerid) = tracer_es_pot(tracerid)
    END DO



    ALLOCATE(tmp_part_KE(ntracer*numprocs))
    ALLOCATE(tmp_part_PE(ntracer*numprocs))
    ALLOCATE(tmp_dphidt(ntracer*numprocs))

    CALL MPI_GATHER(part_KE, ntracer, MPI_DOUBLE_PRECISION, tmp_part_KE, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(part_PE, ntracer, MPI_DOUBLE_PRECISION, tmp_part_PE, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(dphidt, ntracer, MPI_DOUBLE_PRECISION, tmp_dphidt, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    IF (myid == master) THEN

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "part_KE", iframe
	CALL putarr(fidres,dset_name,  tmp_part_KE)
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "part_PE", iframe
	CALL putarr(fidres,dset_name,  tmp_part_PE)
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "dphidt", iframe
	CALL putarr(fidres,dset_name,  tmp_dphidt)
	CALL attach(fidres,dset_name, "time", time)

    END IF

    DEALLOCATE(tmp_part_KE)
    DEALLOCATE(tmp_part_PE)
    DEALLOCATE(tmp_dphidt)




  END IF

  IF (save_cyl) THEN


    ALLOCATE(tmp_tracer_pos_cyl(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_vel_cyl(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_E_cyl(ndim_tracer,ntracer*numprocs))

    CALL MPI_GATHER(tracer_pos_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_vel_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_E_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cyl(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cyl(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cyl(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    IF (myid == master) THEN

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "rcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cyl(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "phicoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cyl(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "zcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cyl(3,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velr_cyl", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cyl(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velphi_cyl", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cyl(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velz_cyl", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cyl(3,:))
	CALL attach(fidres,dset_name, "time", time)
	
	IF (save_Efield) THEN

	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Er_cyl", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_cyl(1,:))
	    CALL attach(fidres,dset_name, "time", time)

	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Ez_cyl", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_cyl(3,:))
	    CALL attach(fidres,dset_name, "time", time)

	END IF

    END IF

    DEALLOCATE(tmp_tracer_pos_cyl)
    DEALLOCATE(tmp_tracer_vel_cyl)
    DEALLOCATE(tmp_tracer_E_cyl)

  END IF



  IF (save_tilt) THEN

    ALLOCATE(tmp_tracer_pos_tilt(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_vel_tilt(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_E_tilt(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_driftrExB_tilt(ntracer*numprocs))
    ALLOCATE(tmp_driftzB_tilt(ntracer*numprocs))
    ALLOCATE(tmp_driftzExB_tilt(ntracer*numprocs))

    CALL MPI_GATHER(tracer_pos_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_vel_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_E_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_tilt(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_tilt(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_tilt(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(driftrExB_tilt, ntracer, MPI_DOUBLE_PRECISION, tmp_driftrExB_tilt, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(driftzB_tilt, ntracer, MPI_DOUBLE_PRECISION, tmp_driftzB_tilt, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(driftzExB_tilt, ntracer, MPI_DOUBLE_PRECISION, tmp_driftzExB_tilt, ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    IF (myid == master) THEN

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "rtiltcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_tilt(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "phitiltcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_tilt(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "ztiltcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_tilt(3,:))
	CALL attach(fidres,dset_name, "time", time)
	 
	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velr_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_tilt(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velphi_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_tilt(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velz_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_tilt(3,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "driftrExB_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_driftrExB_tilt)
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "driftzB_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_driftzB_tilt)
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "driftzExB_tilt", iframe
	CALL putarr(fidres,dset_name,  tmp_driftzExB_tilt)
	CALL attach(fidres,dset_name, "time", time)

	IF (save_Efield) THEN
	     
	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Er_tilt", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_tilt(1,:))
	    CALL attach(fidres,dset_name, "time", time)
       
	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Ez_tilt", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_tilt(3,:))
	    CALL attach(fidres,dset_name, "time", time)
       
	 END IF
    
    END IF

    DEALLOCATE(tmp_tracer_pos_tilt)
    DEALLOCATE(tmp_tracer_vel_tilt)
    DEALLOCATE(tmp_tracer_E_tilt)
    DEALLOCATE(tmp_driftrExB_tilt)
    DEALLOCATE(tmp_driftzB_tilt)
    DEALLOCATE(tmp_driftzExB_tilt)
  END IF



  IF (save_cart) THEN

    ALLOCATE(tmp_tracer_pos_cart(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_vel_cart(ndim_tracer,ntracer*numprocs))
    ALLOCATE(tmp_tracer_E_cart(ndim_tracer,ntracer*numprocs))

    CALL MPI_GATHER(tracer_pos_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_pos_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_pos_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_vel_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_vel_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_vel_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    CALL MPI_GATHER(tracer_E_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cart(1,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cart(2,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    CALL MPI_GATHER(tracer_E_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, tmp_tracer_E_cart(3,:), ntracer, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

    IF (myid == master) THEN



	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "xcartcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cart(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "ycartcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cart(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "zcartcoord", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_pos_cart(3,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velx_cart", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cart(1,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "vely_cart", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cart(2,:))
	CALL attach(fidres,dset_name, "time", time)

	WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "velz_cart", iframe
	CALL putarr(fidres,dset_name,  tmp_tracer_vel_cart(3,:))
	CALL attach(fidres,dset_name, "time", time)

	IF (save_Efield) THEN

	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Ex_cart", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_cart(1,:))
	    CALL attach(fidres,dset_name, "time", time)

	    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var1d", "Ez_cart", iframe
	    CALL putarr(fidres,dset_name,  tmp_tracer_E_cart(3,:))
	    CALL attach(fidres,dset_name, "time", time)

	END IF

    END IF

    DEALLOCATE(tmp_tracer_pos_cart)
    DEALLOCATE(tmp_tracer_vel_cart)
    DEALLOCATE(tmp_tracer_E_cart)

   END IF

   

END SUBROUTINE diagnose_1d