SUBROUTINE stepon
!
!   Advance one time step
!
  USE basic
  USE numerics
  USE tracers
  USE fields
  USE random
  USE futils
  USE interpolate
!
  IMPLICIT none
  INCLUDE 'mpif.h'
!
!   Local vars and arrays
  INTEGER :: j1,j2
!________________________________________________________________________________
! 

  tracerid = 0
  tracertime = step*tracerdeltat
  
  DO j1 = 1,stagbin
! I have a desire to load GBS in pieces
! so that the array is not so big, but
! the overlap makes this harder to deal
! with, so just load the whole thing upfront.
!     IF (Efield_type=='GBS') THEN
!        CALL GBSload(j1)
!        CALL centdiff
!     END IF
     
     IF (static_field) THEN
        tA = (j1-1)*stagdist
        tB = tA
     ELSE
        tA = FLOOR(tracertime/fieldtimestep) + (j1-1)*stagdist
        tB = tA + 1
     END IF

     IF (tA .gt. size(es_potgrid_t,3) .or. tB .gt. size(es_potgrid_t,3))THEN
	write(*,*) 'Too far in time', tA, tB
     END IF

     IF (Efield_type=='GBS' .OR. Efield_type=='sinegrid') THEN
        CALL lineart
     END IF

     DO j2 = 1,stagpart
        tracerid = j2 + (j1 - 1)*stagpart

	IF (Efield_type=='kicks') THEN
	    kickstep1 = FLOOR(kickfreq1*(ran1() + 1))
	    IF (kickrph) THEN
		kickstep2 = FLOOR(kickfreq2*(ran1() + 1))
	    ELSE
		kickstep2 = kickfreq2 +  FLOOR(FLOAT(tracerid)/(FLOAT(ntracer)/200))
	    END IF
	END IF

        IF (ODE_solver=='boris') THEN
           CALL boris
        ELSE IF (ODE_solver=='rk4_3d') THEN
           CALL rk4_3d
        ELSE
           STOP 'invalid trajectory ODE_solver...stopping'
        END IF
     END DO
  END DO

END SUBROUTINE stepon
