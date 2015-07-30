SUBROUTINE endrun
!
!   Terminate the run
!
  USE basic
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
  INTEGER :: ierr, myid
  INTEGER :: master

  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  master = 0
!
!   Local vars and arrays
!
!________________________________________________________________________________
!
  IF( nlend ) THEN
!
!----------------------------------------------------------------------
!              1.   Normal end of run
!
     WRITE(*,*) '   Terminate the run', myid
!
!----------------------------------------------------------------------
!              2.   Abnormal exit
!
  ELSE
     WRITE(*,*) '   Abnormal exit', myid
  END IF
!
!----------------------------------------------------------------------
!              9.   Epilogue
!
!   Create restart file
  CALL chkrst(1)
!
!   Closing all files
!  IF ( myid == master ) THEN
!	CLOSE(lu_in, status='delete')
!  END IF

  CALL memory('f')

END SUBROUTINE endrun
