MODULE basic
!
  USE hashtable
  USE constants
  IMPLICIT NONE
!
!   Basic module for time dependent problems
!
  CHARACTER(len=80) :: label1, label2, label3, label4
!
!   BASIC Namelist
!
  LOGICAL          :: nlres=.FALSE.    ! Restart flag
  LOGICAL          :: nlsave=.TRUE.    ! Checkpoint (save) flag
  LOGICAL          :: newres=.FALSE.   ! New result HDF5 file
  INTEGER          :: nrun=1           ! Number of time steps to run
  DOUBLE PRECISION    :: job_time=3600.0  ! Time allocated to this job in seconds
  DOUBLE PRECISION    :: tmax=100000.0    ! Maximum simulation time (in timesteps)
  DOUBLE PRECISION    :: extra_time=60.0  ! Extra time allocated
  DOUBLE PRECISION    :: dt=1.0           ! Time step
  DOUBLE PRECISION    :: time=0           ! Current simulation time (Init from restart file)
 
!
!   Other basic global vars and arrays
!
  INTEGER :: jobnum                    ! Job number
  INTEGER :: step                      ! Calculation step of this run
  INTEGER :: cstep=0                   ! Current step number (Init from restart file)
  LOGICAL :: nlend                     ! Signal end of run
!
!  List of logical file units
  INTEGER :: lu_in   = 90              ! File duplicated from STDIN
  INTEGER :: lu_stop = 91              ! stop file, see subroutine TESEND
!
!  HDF5 file
  CHARACTER(len=64) :: resfile = "results.h5"   ! Main result file
  CHARACTER(len=64) :: rstfile = "restart.h5"   ! Restart file
  INTEGER           :: fidres                   ! FID for resfile
  INTEGER           :: fidrst                   ! FID for restart file
  TYPE(BUFFER_TYPE) :: hbuf0                    ! Hashtable for 0d var
!
CONTAINS
!
!================================================================================
  SUBROUTINE basic_data
!
!   Define basic data
!
    IMPLICIT NONE
    INCLUDE 'mpif.h'
!
!   Local vars and arrays
    CHARACTER(len=256) :: line
!
    NAMELIST /BASIC/ job_time, extra_time, nrun, tmax, dt, nlres, nlsave, newres

    INTEGER :: myid, master, numprocs, ierr
    CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
    master = 0

!________________________________________________________________________________
!                   1.   Process Standard Input File
!
    IF (myid == master) THEN

    DO
       READ(*,'(a)', END=110) line
       WRITE(lu_in, '(a)') TRIM(line)
    END DO
110 CONTINUE
    REWIND(lu_in)
!________________________________________________________________________________
!                   1.   Label the run
!
    READ(lu_in,'(a)') label1
    READ(lu_in,'(a)') label2
    READ(lu_in,'(a)') label3
    READ(lu_in,'(a)') label4
!
    WRITE(*,'(12x,a/)') label1(1:len_trim(label1))
    WRITE(*,'(12x,a/)') label2(1:len_trim(label2))
    WRITE(*,'(12x,a/)') label3(1:len_trim(label3))
    WRITE(*,'(12x,a/)') label4(1:len_trim(label4))
!________________________________________________________________________________
!                   2.   Read in basic data specific to run
!
    READ(lu_in,basic)
    WRITE(*,basic)

    END IF
    
    CALL MPI_BCAST(job_time,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(extra_time,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nrun,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(tmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nlres,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nlsave,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(newres,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)


!
  END SUBROUTINE basic_data
!================================================================================
  SUBROUTINE daytim(str)
!
!   Print date and time
!
    IMPLICIT NONE
!
    CHARACTER(len=*), INTENT(in) :: str
!
!   Local vars and arrays
    CHARACTER(len=16) :: d, t, dat, clock
!________________________________________________________________________________
!
    CALL DATE_AND_TIME(d,t)
    dat=d(7:8) // '/' // d(5:6) // '/' // d(1:4)
    clock=t(1:2) // ':' // t(3:4) // ':' // t(5:10)
    WRITE(*,'(a,1x,a,1x,a)') str, dat(1:10), clock(1:12)
!
  END SUBROUTINE daytim
!================================================================================
  SUBROUTINE timera(cntrl, str, eltime)
!
!   Timers (cntrl=0/1 to Init/Update)
!
    IMPLICIT NONE
    INTEGER, INTENT(in)                  :: cntrl
    CHARACTER(len=*), INTENT(in)         :: str
    DOUBLE PRECISION, OPTIONAL, INTENT(out) :: eltime
!
    INTEGER, PARAMETER                    :: ncmax=128
    INTEGER, SAVE                         :: icall=0, nc=0
    DOUBLE PRECISION, SAVE                   :: startt0=0.0
    CHARACTER(len=16), SAVE               :: which(ncmax)
    INTEGER                               :: lstr, found, i
    DOUBLE PRECISION                         :: seconds
    DOUBLE PRECISION, DIMENSION(ncmax), SAVE :: startt = 0.0, endt = 0.0
!________________________________________________________________________________
    IF( icall .EQ. 0 ) THEN
       icall = icall+1
       startt0 = seconds()
    END IF

    lstr = LEN_TRIM(str)
    IF( lstr .GT. 0 ) found = loc(str)
!________________________________________________________________________________
!
    SELECT CASE (cntrl)
!
    CASE(-1)    !  Current wall time
       IF( PRESENT(eltime) ) THEN
          eltime = seconds() - startt0
       ELSE
          WRITE(*,'(/a,a,1pe10.3/)') "++ ", ' Wall time used so far = ', seconds() - startt0
       END IF
!
    CASE(0)    !  Init Timer
       IF( found .EQ. 0 ) THEN  !  Called for the 1st time for 'str'
          nc = nc+1
          which(nc) = str(1:lstr)
          found = nc
       END IF
       startt(found) = seconds()
!
    CASE(1)   !  Update timer
       endt(found) = seconds() - startt(found)
       IF( PRESENT(eltime) ) THEN
          eltime = endt(found)
       ELSE
          WRITE(*,'(/a,a,1pe10.2/)') "++ "//str, ' wall clock time = ', endt(found)
       END IF
!
    CASE(2)   !  Update and reset timer
       endt(found) = endt(found) + seconds() - startt(found)
       startt(found) = seconds()
       IF( PRESENT(eltime) ) THEN
          eltime = endt(found)
       END IF
!
    CASE(9)   !  Display all timers
       IF( nc .GT. 0 ) THEN
          WRITE(*,'(a)') "Timer Summary"
          WRITE(*,'(a)') "============="
          DO i=1,nc
             WRITE(*,'(a20,2x,2(1pe12.3))') TRIM(which(i))//":", endt(i)
          END DO
       END IF
!
    END SELECT
!
  CONTAINS
    INTEGER FUNCTION loc(str)
      CHARACTER(len=*), INTENT(in) :: str
      INTEGER :: i, ind
      loc = 0
      DO i=1,nc
          ind = INDEX(which(i), str(1:lstr))
          IF( ind .GT. 0 .AND. LEN_TRIM(which(i)) .EQ. lstr ) THEN
             loc = i
             EXIT
          END IF
       END DO
    END FUNCTION loc
  END SUBROUTINE timera
!================================================================================
END MODULE basic
