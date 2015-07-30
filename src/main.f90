PROGRAM main
!
!   Skeleton for a time dependent program
!   Note: Even in this sequential version, MPI is required
!         because of FUTILS (more specifcally because
!         of the HASTABLE module)!
!
  USE basic
  USE fields, ONLY: Efield_type, sinecentdiff
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
  INTEGER :: ierr, myid
  INTEGER :: master
!    Preparing for parallelization
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  master = 0
!--------------------------------------------------------------------------------
!                         1.   Prologue
  IF (myid == master) THEN
                                                       CALL timera(0, 'Prologue')
	CALL daytim('Start at ')
  END IF
!
!   Define data specific to run
!
  CALL basic_data
!
  IF( .NOT. nlres ) THEN
     CALL newrun
  ELSE
     CALL restart
  END IF
!
!   Compute auxilliary values
!
  CALL auxval
!
! Load GBS data
  IF (Efield_type=='GBS') THEN
     CALL GBSload
     CALL centdiff
  END IF
!
!   Load sinegrid data
  IF (Efield_type=='sinegrid') THEN
     CALL sinegrid
     IF (sinecentdiff) THEN
        CALL centdiff
     END IF
  END IF
!   Initial conditions
!
  IF( .NOT. nlres ) THEN
     CALL initial
  ELSE
     CALL resume
  END IF
!
!   Start or restart the run
!
  CALL start
!
!   Initial diagnostics
!
  CALL diagnose(0)

  IF (myid == master) THEN

                                                       CALL timera(1, 'Prologue')
  END IF
!--------------------------------------------------------------------------------
!                         2.   Time stepping
  IF (myid == master) THEN
                                                       CALL timera(0, 'Main loop')
  END IF
!
  step = 0
  DO
     step = step+1
     cstep = cstep+1
     time = time+dt
     
     CALL stepon
     CALL diagnose(step)
     CALL tesend ! allows an "emergency" stop file with a number of timesteps before end
     IF( nlend ) EXIT
  END DO

  IF (myid == master) THEN
                                                       CALL timera(1, 'Main loop')
  END IF
!--------------------------------------------------------------------------------
!                         9.   Epilogue
  IF (myid == master) THEN
                                                       CALL timera(0, 'Epilogue')

  CALL diagnose(-1)

                                                       CALL timera(1, 'Epilogue')
                                                       CALL timera(9, '')
                                                       CALL timera(-1, '')
  
  END IF

  CALL endrun


  IF (myid == master) THEN

  CALL daytim('Done at ')
  
  END IF
  CALL mpi_finalize(ierr)
END PROGRAM main
