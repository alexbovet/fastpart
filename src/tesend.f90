 SUBROUTINE tesend
!
!   Test for run completion
!
  USE basic
!
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
!   Local vars and arrays
  CHARACTER(len=*), PARAMETER :: stop_file = 'mystop'
  LOGICAL                     :: mlexist
  REAL(kind=db)                :: eltime, step_time, tremain
  INTEGER                     :: rem_steps


  INTEGER :: myid, master, numprocs, ierr
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  master = 0

  IF(myid == master) THEN
!________________________________________________________________________________
!!$  WRITE(*,'(a/)') '=== Test for run completion ==='
!________________________________________________________________________________
!                   1.  Some processors had set nlend
!
  IF( nlend ) THEN
     WRITE(*,'(/a)') 'NLEND set to .TRUE.!'  
     RETURN
  END IF
!________________________________________________________________________________
!                   2.  NRUN modified throught "stop file"
!
  INQUIRE(file=stop_file, exist=mlexist)
  IF( mlexist ) THEN
     OPEN(lu_stop, file=stop_file)
     READ(lu_stop,*) rem_steps     ! Modify remaining steps "on the fly"
     nrun = step + rem_steps
     CLOSE(lu_stop, status='delete')
     WRITE(*,'(/a,i4,a)') 'Stop file found: will exit in', rem_steps, ' steps'
  END IF
  END IF
!________________________________________________________________________________
!                   3.  Test on NRUN
!
  nlend = step .GE. nrun
  IF ( nlend ) THEN 
     WRITE(*,'(/a)') 'NRUN steps done', myid
     RETURN
  END IF
!________________________________________________________________________________
!                   4.  Test on TMAX
!
  nlend = time .GE. tmax
  IF ( nlend ) THEN 
     WRITE(*,'(/a)') 'TMAX reached', myid
     RETURN
  END IF
!________________________________________________________________________________
!                   5.  Test on time allocated to job
!
  IF(myid == master) THEN
  CALL timera(-1, '', eltime)     ! Current elapsed time
  tremain = job_time - eltime
!
  CALL timera(1, 'Main loop', step_time)
  step_time = 1.2 * step_time / step     ! Averaged time per step + 20%
!
  nlend = tremain .LT. (step_time+extra_time)
  IF( nlend ) THEN
     WRITE(*,'(/a,f8.3)') &
          &          'Allocated Job time exhausted:, remaining time =', tremain
  END IF
  

  END IF

  CALL MPI_BCAST(nlend,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)

  RETURN
!
  nlend = .FALSE.

!________________________________________________________________________________
END SUBROUTINE tesend
