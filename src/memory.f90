SUBROUTINE memory(choice)  
  
  USE injection
  USE tracers
  USE numerics
  USE fields
  USE diagnostic

  IMPLICIT NONE

  CHARACTER (1), INTENT (IN) :: choice

  SELECT CASE (choice)

  CASE('i')
     ALLOCATE(tracer_pos_cart(ndim_tracer,ntracer))
     ALLOCATE(tracer_pos_cyl(ndim_tracer,ntracer))
     ALLOCATE(phicheck(ntracer))
     ALLOCATE(nturnphi(ntracer))
     ALLOCATE(x_in(ntracer))
     ALLOCATE(y_in(ntracer))
     ALLOCATE(z_in(ntracer))
     ALLOCATE(pos_cart(ndim_tracer))
     ALLOCATE(pos_cyl(ndim_tracer))
     ALLOCATE(pos_tilt(ndim_tracer))
     ALLOCATE(vel_cart(ndim_tracer))
     ALLOCATE(vel_cyl(ndim_tracer))
     ALLOCATE(vel_tilt(ndim_tracer))
     ALLOCATE(tracer_pos_tilt(ndim_tracer,ntracer))
     ALLOCATE(tracer_vel_cart(ndim_tracer,ntracer)) 
     ALLOCATE(tracer_vel_cyl(ndim_tracer,ntracer))  
     ALLOCATE(tracer_vel_tilt(ndim_tracer,ntracer)) 
     ALLOCATE(tracer_es_pot(ntracer))
     ALLOCATE(es_pot_old(ntracer))
     ALLOCATE(dphidt(ntracer))
     ALLOCATE(rtiltgrid(0:numrpts))
     ALLOCATE(ztiltgrid(0:numzpts))
     ALLOCATE(es_potgrid_t(0:numrpts,0:numzpts,0:numtpts))
     ALLOCATE(es_potgrid_fluct(0:numrpts,0:numzpts,0:numtpts))
     ALLOCATE(dpotdrtiltgrid_t(0:numrpts,0:numzpts,0:numtpts))
     ALLOCATE(dpotdztiltgrid_t(0:numrpts,0:numzpts,0:numtpts))
     ALLOCATE(d2potdrdztiltgrid_t(0:numrpts,0:numzpts,0:numtpts))
     ALLOCATE(es_potgrid(0:numrpts,0:numzpts))
     ALLOCATE(dpotdrtiltgrid(0:numrpts,0:numzpts))
     ALLOCATE(dpotdztiltgrid(0:numrpts,0:numzpts))
     ALLOCATE(d2potdrdztiltgrid(0:numrpts,0:numzpts))
     ALLOCATE(tracer_E_cart(ndim_tracer,ntracer))  
     ALLOCATE(tracer_E_cyl(ndim_tracer,ntracer))  
     ALLOCATE(tracer_E_tilt(ndim_tracer,ntracer))  
     ALLOCATE(B_cart(ndim_tracer))
     ALLOCATE(E_cart(ndim_tracer))
     ALLOCATE(E_cyl(ndim_tracer))
     ALLOCATE(E_tilt(ndim_tracer))
     ALLOCATE(part_KE(ntracer))
     ALLOCATE(part_PE(ntracer))
     ALLOCATE(driftrExB_tilt(ntracer))
     ALLOCATE(driftzB_tilt(ntracer))
     ALLOCATE(driftzExB_tilt(ntracer))

  CASE('f')
     DEALLOCATE(tracer_pos_cart)
     DEALLOCATE(tracer_pos_cyl)
     DEALLOCATE(phicheck)
     DEALLOCATE(nturnphi)
     DEALLOCATE(x_in)
     DEALLOCATE(y_in)
     DEALLOCATE(z_in)
     DEALLOCATE(pos_cart)
     DEALLOCATE(pos_cyl)
     DEALLOCATE(pos_tilt)
     DEALLOCATE(vel_cart)
     DEALLOCATE(vel_cyl)
     DEALLOCATE(vel_tilt)
     DEALLOCATE(tracer_pos_tilt)
     DEALLOCATE(tracer_vel_cart)
     DEALLOCATE(tracer_vel_cyl)
     DEALLOCATE(tracer_vel_tilt)
     DEALLOCATE(tracer_es_pot)
     DEALLOCATE(es_pot_old)
     DEALLOCATE(dphidt)
     DEALLOCATE(rtiltgrid)
     DEALLOCATE(ztiltgrid)
     DEALLOCATE(es_potgrid_t)
     DEALLOCATE(es_potgrid_fluct)
     DEALLOCATE(dpotdrtiltgrid_t)
     DEALLOCATE(dpotdztiltgrid_t)
     DEALLOCATE(d2potdrdztiltgrid_t)
     DEALLOCATE(es_potgrid)
     DEALLOCATE(dpotdrtiltgrid)
     DEALLOCATE(dpotdztiltgrid)
     DEALLOCATE(d2potdrdztiltgrid)
     DEALLOCATE(B_cart)
     DEALLOCATE(E_cart)
     DEALLOCATE(E_cyl)
     DEALLOCATE(E_tilt)
     DEALLOCATE(tracer_E_cart)
     DEALLOCATE(tracer_E_cyl)
     DEALLOCATE(tracer_E_tilt)
     DEALLOCATE(part_KE)
     DEALLOCATE(part_PE)
     DEALLOCATE(driftrExB_tilt)
     DEALLOCATE(driftzB_tilt)
     DEALLOCATE(driftzExB_tilt)
     
  END SELECT


END SUBROUTINE memory
