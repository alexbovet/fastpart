SUBROUTINE resume
!
!   Resume from previous run
!
  IMPLICIT NONE
!
!   Local vars and arrays
!________________________________________________________________________________
  WRITE(*,'(a)') '   Resume from previous run'
!________________________________________________________________________________
! 
!   Open and read initial conditions from restart file
  CALL chkrst(0)
!
END SUBROUTINE resume
