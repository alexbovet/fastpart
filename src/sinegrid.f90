  SUBROUTINE sinegrid

    USE basic
    USE fields
    USE constants
    USE interpolate
    !
    ! Generates a sinusoidal electric field on a grid.

    INTEGER :: k,l,tau

    rtiltgrid           = 0.
    ztiltgrid           = 0.
    es_potgrid_t        = 0.
    es_potgrid          = 0.
    dpotdrtiltgrid      = 0.
    dpotdztiltgrid      = 0.
    d2potdrdztiltgrid   = 0.
    dpotdrtiltgrid_t    = 0.
    dpotdztiltgrid_t    = 0.
    d2potdrdztiltgrid_t = 0.


    DO k = 0,numrpts - 1
      rtiltgrid(k) = Lr*k/numrpts
      DO l = 0,numzpts - 1
        ztiltgrid(l) = Lz*l/numzpts
	DO tau = 0,numtpts
	  es_potgrid_t(k,l,tau)      =       E_mag1*sin(kr*rtiltgrid(k)*sinmoder)*sin(kz*ztiltgrid(l)*sinmodez)*cos(sinmodet*tau)
          IF (sinecentdiff==.FALSE.) THEN
             dpotdrtiltgrid_t(k,l,tau)    =             E_mag1*kr*sinmoder*cos(kr*rtiltgrid(k)*sinmoder)*sin(kz*ztiltgrid(l)*sinmodez)*cos(sinmodet*tau)
             dpotdztiltgrid_t(k,l,tau)    =             E_mag1*kz*sinmodez*sin(kr*rtiltgrid(k)*sinmoder)*cos(kz*ztiltgrid(l)*sinmodez)*cos(sinmodet*tau)
             d2potdrdztiltgrid_t(k,l,tau) = E_mag1*kz*kr*sinmoder*sinmodez*cos(kr*rtiltgrid(k)*sinmoder)*cos(kz*ztiltgrid(l)*sinmodez)*cos(sinmodet*tau)
          END IF
	END DO
      END DO
    END DO

    rtiltgrid(numrpts) = rtiltgrid(0)
    ztiltgrid(numzpts) = ztiltgrid(0)

    IF (sinecentdiff==.FALSE.) THEN
       es_potgrid_t(numrpts,:,:) = es_potgrid_t(0,:,:)
       es_potgrid_t(:,numzpts,:) = es_potgrid_t(:,0,:)

       dpotdrtiltgrid_t(numrpts,:,:) = dpotdrtiltgrid_t(0,:,:)
       dpotdrtiltgrid_t(:,numzpts,:) = dpotdrtiltgrid_t(:,0,:)
       
       dpotdztiltgrid_t(numrpts,:,:) = dpotdztiltgrid_t(0,:,:)
       dpotdztiltgrid_t(:,numzpts,:) = dpotdztiltgrid_t(:,0,:)
       
       d2potdrdztiltgrid_t(numrpts,:,:) = d2potdrdztiltgrid_t(0,:,:)
       d2potdrdztiltgrid_t(:,numzpts,:) = d2potdrdztiltgrid_t(:,0,:)
    END IF

  END SUBROUTINE sinegrid
