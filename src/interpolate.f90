MODULE interpolate
!
  IMPLICIT none
!
CONTAINS
!
  SUBROUTINE bilinear (x1,x2,ygrid,y1grid,y2grid,ansy,ansy1,ansy2)
    !
    USE nrutil, ONLY : nrerror
    USE fields, ONLY: numrpts,numzpts,Lr,Lz
    USE basic
    IMPLICIT none
    DOUBLE PRECISION, INTENT(IN) :: x1,x2
    DOUBLE PRECISION, DIMENSION(0:numrpts,0:numzpts), INTENT(IN)  :: ygrid,y1grid,y2grid
    DOUBLE PRECISION, INTENT(OUT) :: ansy,ansy1,ansy2
    ! Linear interpolation method for two dimensional function...
    ! Given an m by n array of
    ! function values ya(0:m, 0:n), tabulated at the grid points defined elsewhere; and
    ! given values x1 and x2 of the independent variables; this routine returns an interpolated
    ! function value y.
    INTEGER :: j,k,m,n,ndum,x1l,x2l,x1u,x2u
    DOUBLE PRECISION :: t,u,d1,d2
    DOUBLE PRECISION, DIMENSION(4) :: y,y1,y2

    m = numrpts ! these require fields_mod and are input parameters
    n = numzpts

    d1 = Lr/m
    d2 = Lz/n

    IF (x1 .lt. 0 .or. x2 .lt. 0) STOP 'killed: trying to interpolate for value lt zero'

    ! determine the interval in which the point of interest lies
    ! the following method finds the index
    x1l = int(x1/d1)
    x1u = x1l + 1
    x2l = int(x2/d2)
    x2u = x2l + 1
    ! find the tabulated function values for performing interpolation
    IF (x1u == x1l .or. x2u == x2l) CALL &
      nrerror('bilinear: problem with input values - boundary pair equal?')

    y(1) = ygrid(x1l,x2l)      ! numbered counterclockwise from lower left, as required by bcucof
    y(2) = ygrid(x1u,x2l)
    y(3) = ygrid(x1u,x2u)
    y(4) = ygrid(x1l,x2u)

    y1(1) = y1grid(x1l,x2l)     ! numbered counterclockwise from lower left, as required by bcucof
    y1(2) = y1grid(x1u,x2l)
    y1(3) = y1grid(x1u,x2u)
    y1(4) = y1grid(x1l,x2u)

    y2(1) = y2grid(x1l,x2l)     ! numbered counterclockwise from lower left, as required by bcucof
    y2(2) = y2grid(x1u,x2l)
    y2(3) = y2grid(x1u,x2u)
    y2(4) = y2grid(x1l,x2u)

    ! compute the bilinear coefficients - the normalized distance from the
    ! interpolant point to the grid point
    t = x1/d1 - x1l
    u = x2/d2 - x2l 
    ! find the interpolated function value
    ansy  = (1. - t)*(1. - u)*y (1) + t*(1. - u)*y (2) + t*u*y (3) + (1. - t)*u*y (4)
    ansy1 = (1. - t)*(1. - u)*y1(1) + t*(1. - u)*y1(2) + t*u*y1(3) + (1. - t)*u*y1(4)
    ansy2 = (1. - t)*(1. - u)*y2(1) + t*(1. - u)*y2(2) + t*u*y2(3) + (1. - t)*u*y2(4)

  END SUBROUTINE bilinear

  SUBROUTINE bcuint(x1,x2,ygrid,y1grid,y2grid,y12grid,ansy,ansy1,ansy2)
    USE fields, ONLY: numrpts,numzpts,Lr,Lz
    USE basic
    USE tracers
    USE nrutil, ONLY : nrerror, bcucof
    IMPLICIT none
    DOUBLE PRECISION, INTENT(IN)  :: x1,x2
    DOUBLE PRECISION, DIMENSION(0:numrpts,0:numzpts), INTENT(IN)  :: ygrid,y1grid,y2grid,y12grid
    DOUBLE PRECISION, INTENT(OUT) :: ansy,ansy1,ansy2
    ! This routine has been modified from the NR original to work more smoothly with the 
    ! framework of this code (KBG@CRPP, 5/3/2010).
    !Bicubic interpolation within a grid square. Needed quantities are y,y1,y2,y12 (as described
    !in bcucof); x1l and x1u, the lower and upper coordinates of the grid square in the 1-
    !direction; x2l and x2u likewise for the 2-direction; and x1,x2, the coordinates of the
    !desired point for the interpolation. The interpolated function value is returned as ansy,
    !and the interpolated gradient values as ansy1 and ansy2. This routine calls bcucof.
    INTEGER :: i,m,n,x1l,x2l,x1u,x2u
    DOUBLE PRECISION :: t,u,d1,d2
    DOUBLE PRECISION, DIMENSION(4) :: y,y1,y2,y12
    DOUBLE PRECISION, DIMENSION(4,4) :: c

    m = numrpts ! these require fields_mod and are input parameters
    n = numzpts

    d1 = Lr/m
    d2 = Lz/n
    ! determine the interval in which the point of interest lies
    ! the following method finds the index
    x1l = int(x1/d1) 
    x1u = x1l + 1.
    x2l = int(x2/d2)
    x2u = x2l + 1.

    IF (x1 .lt. 0 .or. x2 .lt. 0) STOP 'killed: trying to interpolate for value lt zero'

    IF (x1u .gt. m .or. x2u .gt. n .or. x1l .gt. m .or. x2l .gt. n) THEN
	write(*,*) 'Index problem: x1l=', x1l, ' x1u=', x1u, ' x2l=', x2l, ' x2u=', x2u, ' numrpts=', numrpts, ' numzpts=', numzpts
	STOP 'Killed'
    END IF

    y(1) = ygrid(x1l,x2l)      ! numbered counterclockwise from lower left, as required by bcucof
    y(2) = ygrid(x1u,x2l)
    y(3) = ygrid(x1u,x2u)
    y(4) = ygrid(x1l,x2u)

    y1(1) = y1grid(x1l,x2l)     ! numbered counterclockwise from lower left, as required by bcucof
    y1(2) = y1grid(x1u,x2l)
    y1(3) = y1grid(x1u,x2u)
    y1(4) = y1grid(x1l,x2u)
   
    y2(1) = y2grid(x1l,x2l)     ! numbered counterclockwise from lower left, as required by bcucof
    y2(2) = y2grid(x1u,x2l)
    y2(3) = y2grid(x1u,x2u)
    y2(4) = y2grid(x1l,x2u)

    y12(1) = y12grid(x1l,x2l) ! numbered counterclockwise from lower left, as required by bcucof
    y12(2) = y12grid(x1u,x2l)
    y12(3) = y12grid(x1u,x2u)
    y12(4) = y12grid(x1l,x2u)

    CALL bcucof(y,y1,y2,y12,d1,d2,c) ! get the c values from nrutil.f90

    IF (x1u == x1l .or. x2u == x2l) CALL &
      nrerror('bcuint: problem with input values - boundary pair equal?')

    t= x1/d1 - x1l ! same as for bilinear
    u= x2/d2 - x2l 

    ansy=0.0
    ansy1=0.0
    ansy2=0.0
 
    DO i=4,1,-1 !Equation (3.6.6) in NR Fortran90 - linear combinations of the c values
      ansy =  t*ansy + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i,1)
      ansy1 = u*ansy1 + (3.*c(4,i)*t + 2.*c(3,i))*t + c(2,i)
      ansy2 = t*ansy2 + (3.*c(i,4)*u + 2.*c(i,3))*u + c(i,2)
    END DO
    ansy1=ansy1/d1 ! renormalize the derivatives to the grid spacing
    ansy2=ansy2/d2

  END SUBROUTINE bcuint

  SUBROUTINE lineart
    !
    USE basic
    USE fields
    USE tracers
    IMPLICIT none
!  
    ! One dimension linear interpolation in time for electrostatic potential.
    ! Uses the fieldtimestep and numtpts input parameters and gives
    ! values to es_potgrid using the 3 dimensional (r,z',t)
    ! array es_potgrid_t. 
    DOUBLE PRECISION    :: intert, sintert
    !INTEGER             :: i,GBSidx                                                                         
    !CHARACTER(LEN=64)   :: tail 

    ! this definition means that tracerdeltat should be a multiple of fieldtimestep
    intert = tracertime/fieldtimestep - FLOOR(tracertime/fieldtimestep)
    sintert = 1. - intert
!    i = 1
!    DO WHILE (i < stagbin) 
!      IF (strmfstart < tA .AND. tA < strmfstart + i*stagdist) THEN
!        GBSidx = i
!      ELSE
!        i = i + 1
!      END IF
!    END DO 
!
!    write(tail,7001) i                                                                                                    
!    7001 FORMAT(i2.2)                                                                                                           
!    arraystrmf = 'es_pot'//tail

    es_potgrid        = sintert*es_potgrid_t(:,:,tA) + intert*es_potgrid_t(:,:,tB)
    dpotdrtiltgrid    = sintert*dpotdrtiltgrid_t(:,:,tA) + intert*dpotdrtiltgrid_t(:,:,tB)
    dpotdztiltgrid    = sintert*dpotdztiltgrid_t(:,:,tA) + intert*dpotdztiltgrid_t(:,:,tB)
    d2potdrdztiltgrid = sintert*d2potdrdztiltgrid_t(:,:,tA) + intert*d2potdrdztiltgrid_t(:,:,tB)

  END SUBROUTINE lineart

END MODULE interpolate
