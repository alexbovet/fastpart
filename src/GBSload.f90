
  SUBROUTINE GBSload

    USE basic
    USE fields
    USE constants
    USE futils
    INCLUDE 'mpif.h'
    !
    ! Loads a GBS streamfunction and computes the derivatives which 
    ! will be used to find the electric field.
    INTEGER :: k,l,i,tau,tauplus
    DOUBLE PRECISION, DIMENSION(numzpts,numrpts)         :: a
    DOUBLE PRECISION, DIMENSION(numrpts,numzpts)         :: b
    DOUBLE PRECISION, DIMENSION(numrpts)                 :: avg_pot_r
    CHARACTER(LEN=64)			                 :: arraystrmf,tail

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: es_potgrid_t_tmp

    INTEGER :: myid, master, numprocs, processor,ierr,stat(MPI_STATUS_SIZE),GBSmaxtimesteps
    CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
    master = 0
    IF (myid == master) THEN  ! Loading on one processor only, that will then send portion of the perturbation to other processors

!________________________________________________________________________________
  WRITE(*,'(a/)') '=== Load GBS phi(r,z_tilt,t) grid ==='
!________________________________________________________________________________
! 
    es_potgrid_t = 0.
    a=0.
    b=0.

    IF (strmfstart == 0) THEN
      write(*,*) 'user must specify a starting record number (strmfstart) from the GBS strmf output -'
      STOP 'leading zeros are built in...stopping'
    ELSE
      CALL openf(GBSfile, fidGBS, 'r')
      CALL getsize(fidGBS, 'data/var2d/time', GBSmaxtimesteps)
    END IF

    IF ((all_bins_same_pot == .FALSE. .AND. GBSmaxtimesteps .LE. strmfstart + numprocs*(numtpts-stalength+stagdist) - stagdist + staglength) .OR. (all_bins_same_pot == .TRUE. .AND. GBSmaxtimesteps .LE. strmfstart + staglength)) THEN
	CALL closef(fidGBS)
	STOP 'GBS simulation is not longer enough for requested simulation...stopping'
    END IF





	    IF (ALLOCATED(es_potgrid_fluct)) THEN
		DEALLOCATE(es_potgrid_fluct)
		ALLOCATE(es_potgrid_fluct(0:numrpts,0:numzpts,0:numtpts))
	    END IF
		ALLOCATE(es_potgrid_t_tmp(0:numrpts,0:numzpts,0:numtpts))

	    tauplus = 0

	    DO tau = 0,numtpts
	      a = 0.0
	      b = 0.0
	      tauplus = strmfstart + tau !+ stagdist*(j1 - 1)
	      write(tail,7001) tauplus
	      7001 FORMAT(i6.6)
	      arraystrmf = 'data/var2d/strmf/'//tail	
	      write(*,*) 'arraystrmf ',arraystrmf
	      CALL getarr(fidGBS, arraystrmf, a)
	      ! transpose the array due to the different GBS convention
	      b = transpose(a)
	      ! here, shift the array into the 0:,0: position
	      DO k = 0,numrpts-1
		  DO l = 0,numzpts-1 
		    ! multiply by this normalization factor
		    es_potgrid_t_tmp(k,l,tau) = b(k+1,l+1)*phibeamnorm
		  END DO
	      END DO
	    END DO

	    ! adjust the GBS potential to increase the fluctuations
	    ! by a factor of 'GBS_factor' over the background, using 
	    ! a file 'avgpotefile' that contains a profile of phi(r) averaged over
	    ! z and time

	    ! adjust the GBS potential to increase the background
	    ! by a factor of 'GBS_bg' while holding fluctuations constant
	    ! using a file 'avgpotefile' that contains phi(r) averaged over
	    ! z and time

	    IF (GBS_factor /= 1 .OR. GBS_bg /= 1) THEN
		open(unit=3,file=avgpotefile, status='old', readonly)
		DO k = 0,numrpts
		  es_potgrid_fluct(k,:,:) = es_potgrid_t_tmp(k,:,:) - phibeamnorm*avg_pot_r(k+1)
		END DO
		es_potgrid_t_tmp = GBS_bg*(es_potgrid_t_tmp - es_potgrid_fluct) + GBS_factor*es_potgrid_fluct 
	    END IF

	    ! multiply the e.s. potential by a factor 'GBS_mag'

	    IF (GBS_mag /= 1) THEN
		es_potgrid_t_tmp = GBS_mag*es_potgrid_t_tmp
	    END IF

	    IF (vert_reverse) THEN
		DO l = 0,numzpts-1
		  es_potgrid_t_tmp(:,l,:) = es_potgrid_t_tmp(:,numzpts-1-l,:)
		END DO
	    END IF

	    IF (all_bins_same_pot == .FALSE.) THEN
		strmfstart = strmfstart + numtpts - staglength + stagdist ! Each processor recieve a different portion of the perturbation
	    END IF

            



    DO processor = master+1,numprocs-1

	    IF ( ALLOCATED(es_potgrid_t) .AND. ALLOCATED(es_potgrid_fluct)) THEN
		DEALLOCATE(es_potgrid_t)
		ALLOCATE(es_potgrid_t(0:numrpts,0:numzpts,0:numtpts))
		DEALLOCATE(es_potgrid_fluct)
		ALLOCATE(es_potgrid_fluct(0:numrpts,0:numzpts,0:numtpts))
	    END IF

	    tauplus = 0

	    DO tau = 0,numtpts
	      a = 0.0
	      b = 0.0
	      tauplus = strmfstart + tau !+ stagdist*(j1 - 1)
	      write(tail,7001) tauplus
	      arraystrmf = 'data/var2d/strmf/'//tail	
	      write(*,*) 'arraystrmf ',arraystrmf
	      CALL getarr(fidGBS, arraystrmf, a)
	      ! transpose the array due to the different GBS convention
	      b = transpose(a)
	      ! here, shift the array into the 0:,0: position
	      DO k = 0,numrpts-1
		  DO l = 0,numzpts-1 
		    ! multiply by this normalization factor
		    es_potgrid_t(k,l,tau) = b(k+1,l+1)*phibeamnorm
		  END DO
	      END DO
	    END DO

	    ! adjust the GBS potential to increase the fluctuations
	    ! by a factor of 'GBS_factor' over the background, using 
	    ! a file 'avgpotefile' that contains a profile of phi(r) averaged over
	    ! z and time

	    ! adjust the GBS potential to increase the background
	    ! by a factor of 'GBS_bg' while holding fluctuations constant
	    ! using a file 'avgpotefile' that contains phi(r) averaged over
	    ! z and time

	    IF (GBS_factor /= 1 .OR. GBS_bg /= 1) THEN
		open(unit=3,file=avgpotefile, status='old', readonly)
		DO k = 0,numrpts
		  es_potgrid_fluct(k,:,:) = es_potgrid_t(k,:,:) - phibeamnorm*avg_pot_r(k+1)
		END DO
		es_potgrid_t = GBS_bg*(es_potgrid_t - es_potgrid_fluct) + GBS_factor*es_potgrid_fluct 
	    END IF

	    ! multiply the e.s. potential by a factor 'GBS_mag'

	    IF (GBS_mag /= 1) THEN
		es_potgrid_t = GBS_mag*es_potgrid_t
	    END IF

	    IF (vert_reverse) THEN
		DO l = 0,numzpts-1
		  es_potgrid_t(:,l,:) = es_potgrid_t(:,numzpts-1-l,:)
		END DO
	    END IF

	    IF (all_bins_same_pot == .FALSE.) THEN
		strmfstart = strmfstart + numtpts - staglength + stagdist ! Each processor recieve a different portion of the perturbation
	    END IF

	    ! sending turbulence to other processors
	    
	    IF (processor .NE. master) THEN
		CALL MPI_SEND(es_potgrid_t,(numrpts+1)*(numzpts+1)*(numtpts+1),MPI_DOUBLE_PRECISION,processor,0,MPI_COMM_WORLD,ierr)
	    END IF
    END DO

    CALL closef(fidGBS)
    END IF


    IF (myid .NE. master) THEN
	CALL MPI_RECV(es_potgrid_t,(numrpts+1)*(numzpts+1)*(numtpts+1),MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
    END IF

    IF (myid == master) THEN
        es_potgrid_t = es_potgrid_t_tmp
        DEALLOCATE(es_potgrid_t_tmp)
    END IF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  END SUBROUTINE GBSload
