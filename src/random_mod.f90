MODULE random

!  
!  Pay attention when using random numbers in parallel simulation (For now only master process uses them, for initialization)
!  

CONTAINS

  SUBROUTINE gasdev_s(harvest)
    !
    ! Returns in harvest a normally distributed deviate with zero mean and unit standard deviation
    ! using intrinsic random_number as the source of uniform deviates
    !
    IMPLICIT none

    DOUBLE PRECISION, intent(out) :: harvest

    DOUBLE PRECISION        :: rsq, v1,v2
    DOUBLE PRECISION, SAVE  :: g
    LOGICAL, SAVE           :: gaus_stored = .false.

    IF (gaus_stored) THEN       ! We have an extra deviate,
       harvest = g              ! so return it
       gaus_stored = .false.    ! and unset the flag.
    ELSE
       DO
          v1 = 2.0*ran1() - 1.0     ! move the values to the box -1 to +1
          v2 = 2.0*ran1() - 1.0
          rsq = v1**2 + v2**2   ! see if they are in the unit circle
          if (rsq > 0.0 .and. rsq < 1.0) exit
       END DO                   ! otherwise, try again
       rsq = sqrt(-2.0*log(rsq)/rsq) ! Now make the Box-Muller transformation to
       harvest = v1*rsq         ! get two normal deviates. Return one and
       g = v2*rsq               ! save the other for next time
       gaus_stored = .true.                    ! Set flag.
    END IF

  END SUBROUTINE gasdev_s

  FUNCTION ran1()
    IMPLICIT none
    DOUBLE PRECISION :: ran1, x
    CALL RANDOM_NUMBER(x)
    ran1 = x
  END FUNCTION ran1

  SUBROUTINE init_random_seed()
    INTEGER :: i,n,clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37*(/ (i-1, i=1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)

  END SUBROUTINE init_random_seed

  SUBROUTINE init_default_seed()
    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    seed = 1245

    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)

  END SUBROUTINE init_default_seed

END MODULE random
