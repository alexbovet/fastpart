SUBROUTINE auxval
!
!   Set auxiliary values
!
  USE basic
  USE fields 
  USE tracers
  USE numerics
  IMPLICIT none
!
  INCLUDE 'mpif.h'
!   Local vars and arrays
  INTEGER :: k1,j2
  DOUBLE PRECISION :: length
!
  INTEGER :: ierr, myid, master, numprocs
!    Preparing for parallelization
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  master = 0
!________________________________________________________________________________
  IF (myid == master) THEN
      WRITE(*,'(a/)') '=== Set auxiliary values ==='
  END IF
!________________________________________________________________________________
! 
! Pay attention that stagbin is a multiple of numprocs (number of processors)
   stagbin = stagbin/numprocs
   ntracer = stagpart*stagbin
   IF (all_bins_same_pot == .FALSE.) THEN
	numtpts = staglength+(stagbin-1)*stagdist+1
   ELSE
	numtpts = staglength
   END IF


   length = fieldtimestep*staglength/tracerdtfactor

   IF (length > nrun) THEN
      STOP 'length is greater than nrun...stopping'
   END IF

   IF (Bfield_type=='slab') THEN
      Rcurv = 0.
      E_basis = 'cart'
      IF (Bx > 0 .OR. Bz > 0) THEN
         STOP 'in slab field, only abs(By) can be larger than zero...stopping'
      END IF
   ELSE
      E_basis = 'tilt'
   END IF

   IF (Bratio > 0) THEN
      theta=atan(Bratio)
      sineoftheta = sin(theta)
      cosineoftheta = cos(theta)
   ELSE 
      theta = 0.
      sineoftheta = 0.
      cosineoftheta = 1.
   END IF

   tracerdeltat=dt*tracerdtfactor

   IF (Lr .NE. 0.) THEN
     kr = 2*pi/Lr
   ELSE 
     kz = 0.
   END IF
   
   IF (Lz .NE. 0.) THEN
     kz = 2*pi/Lz
   ELSE 
     kz = 0.
   END IF

   dr = Lr/numrpts
   dz = Lz/numzpts

   IF (ODE_solver=='boris' .OR. ODE_solver=='rk4_3d') THEN
        ndim_tracer=3
   ELSE IF (ODE_solver=='sympl') THEN
        ndim_tracer=2
   ELSE 
       write(*,*)'Error: Invalid ODE solver'
   END IF
!
! Allocate memory for global arrays
!
  CALL memory('i')
!
  
  IF (Efield_type=='kicks' .AND. readkicks) THEN
     IF (myid == master)THEN
	OPEN(UNIT = 7,FILE = "./kickrph.Gauss")
	DO j2 = 1,1000
	   DO k1 = 1,323
	      READ(7,*) kicksteps(j2,k1)
	   END DO
	END DO
     END IF
     CALL MPI_BCAST(kicksteps,1000*323,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  END IF
  
END SUBROUTINE auxval
