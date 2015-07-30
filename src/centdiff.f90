 SUBROUTINE centdiff

    USE basic
    USE fields
    USE constants
    USE futils
    !
    ! Computes the derivatives which 
    ! will be used to find the electric field.

    INTEGER           :: k,l,i,j
    DOUBLE PRECISION :: d1,d2,d1inv,d2inv,d12inv
    INTEGER :: phi_unit=45,Er_unit=46,Ez_unit=47,Erz_unit=48
    character (305) :: file_phi,file_Er,file_Ez,file_Erz
    file_phi = "phi1"
    file_Er = "Er1"
    file_Ez = "Ez1"
    file_Erz = "Erz1"

    d1 = dr
    d2 = dz

    d1inv = 0.5/d1
    d2inv = 0.5/d2
    d12inv = d1inv*d2inv

    ! next, fill in the last row/column to match the zeroths
    es_potgrid_t(numrpts,:,:) = es_potgrid_t(0,:,:)
    es_potgrid_t(:,numzpts,:) = es_potgrid_t(:,0,:)

    ! find the specialty values for the zeroth and last rows
    dpotdrtiltgrid_t(0,:,:) = d1inv*(es_potgrid_t(1,:,:) - es_potgrid_t(numrpts - 1,:,:))
    dpotdrtiltgrid_t(numrpts,:,:) = dpotdrtiltgrid_t(0,:,:)
    ! do the usual centered differencing
    DO i=1,numrpts-1
      dpotdrtiltgrid_t(i,:,:) = d1inv*(es_potgrid_t(i + 1,:,:) - es_potgrid_t(i - 1,:,:))
    END DO
    
    ! find the specialty values for the zeroth and last columns
    dpotdztiltgrid_t(:,0,:) = d2inv*(es_potgrid_t(:,1,:)  - es_potgrid_t(:,numzpts - 1,:))
    dpotdztiltgrid_t(:,numzpts,:) = dpotdztiltgrid_t(:,0,:)
    ! do the usual centered differencing
    DO j=1,numzpts-1
      dpotdztiltgrid_t(:,j,:) = d2inv*(es_potgrid_t(:,j + 1,:) - es_potgrid_t(:,j - 1,:))
    END DO

    ! find the specialty values for the zeroth and last columns and rows
    DO j=1,numrpts-1
      d2potdrdztiltgrid_t(j,0,:) = &
      d12inv*(es_potgrid_t(j+1,1,:) - es_potgrid_t(j+1,numzpts-1,:) - es_potgrid_t(j-1,1,:) + & 
      es_potgrid_t(j-1,numzpts-1,:))
    END DO
    DO k=1,numzpts-1
      d2potdrdztiltgrid_t(0,k,:) = &
      d12inv*(es_potgrid_t(1,k+1,:) - es_potgrid_t(numrpts-1,k+1,:) - es_potgrid_t(1,k-1,:) + &
      es_potgrid_t(numrpts-1,k-1,:))
    END DO
    ! the 0,0 position requires extra care
    d2potdrdztiltgrid_t(0,0,:) = &
    d12inv*(es_potgrid_t(1,1,:) - es_potgrid_t(numrpts-1,1,:) - es_potgrid_t(1,numzpts-1,:) + &
    es_potgrid_t(numrpts-1,numzpts-1,:))
    ! fill in the last row/column to match the zeroth
    d2potdrdztiltgrid_t(:,numzpts,:) = d2potdrdztiltgrid_t(:,0,:)
    d2potdrdztiltgrid_t(numrpts,:,:) = d2potdrdztiltgrid_t(0,:,:)
    ! do the usual centered differencing for the bidirectional derivative
    DO j=1,numrpts-1
      DO k=1,numzpts-1
        d2potdrdztiltgrid_t(j,k,:) = &
           d12inv*(es_potgrid_t(j+1,k+1,:) - es_potgrid_t(j+1,k-1,:) - es_potgrid_t(j-1,k+1,:) + &
           es_potgrid_t(j-1,k-1,:))
      END DO
    END DO
    
    IF (firstframeout) THEN
      open(phi_unit, file=file_phi, STATUS='REPLACE') ! the choice of tracer_unit here is completely arbitrary
      open(Er_unit, file=file_Er, STATUS='REPLACE') ! the choice of tracer_unit here is completely arbitrary
      open(Ez_unit, file=file_Ez, STATUS='REPLACE') ! the choice of tracer_unit here is completely arbitrary
      open(Erz_unit, file=file_Erz, STATUS='REPLACE') ! the choice of tracer_unit here is completely arbitrary
      do j=0,numrpts
      	do k=0,numzpts
	        write(phi_unit,*) j, k, es_potgrid_t(j,k,0) 
	        write(Er_unit,*)  j, k, dpotdrtiltgrid_t(j,k,0) 
	        write(Ez_unit,*)  j, k, dpotdztiltgrid_t(j,k,0) 
	        write(Erz_unit,*) j, k, d2potdrdztiltgrid_t(j,k,0) 
     	end do
      end do
      close(phi_unit)
      close(Er_unit)
      close(Ez_unit)
      close(Erz_unit)
   END IF

  END SUBROUTINE centdiff
