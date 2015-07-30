SUBROUTINE chkrst(flag)
!
!   Process checkpoint/restart file
!
  USE basic
  USE futils
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER, INTENT(in) :: flag

  INTEGER :: myid, master, numprocs, ierr
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  master = 0

  IF(myid == master) THEN

!
!  Local vars and arrays
!________________________________________________________________________________
!
  SELECT CASE(flag)
!________________________________________________________________________________
!                   1.  Open and read restart file
!
  CASE(0)
     CALL openf(rstfile, fidrst)
     CALL getatt(fidrst, '/Basic', 'cstep', cstep)
     CALL getatt(fidrst, '/Basic', 'time', time)
     CALL closef(fidrst)
     WRITE(*,'(3x,a)') &
          &     "Reading from restart file "//TRIM(rstfile)//" completed!"
!________________________________________________________________________________
!                   2.  Create and write to restart file (DP reals)
!
  CASE(1)
     IF( .NOT. nlsave ) RETURN
     CALL mv2bk(rstfile)
     CALL creatf(rstfile, fidrst, real_prec='d', desc='Restart file')
     CALL creatg(fidrst, '/Basic', 'Basic data')
     CALL attach(fidrst, '/Basic', 'cstep', cstep)
     CALL attach(fidrst, '/Basic', 'time', time)
     CALL attach(fidrst, '/Basic', 'jobnum', jobnum)
     CALL closef(fidrst)
     WRITE(*,'(3x,a)') &
          &     "Writing to restart file "//TRIM(rstfile)//" completed!"
!
  END SELECT

  END IF
!
 END SUBROUTINE chkrst