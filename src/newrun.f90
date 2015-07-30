SUBROUTINE newrun
!
!   Additional data specific for a new run
!
 USE basic
 USE fields
 USE injection
 USE diagnostic
 USE numerics
 USE tracers
 USE constants
!
  IMPLICIT none
  INCLUDE 'mpif.h'

  INTEGER :: myid, master, numprocs, ierr
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  master = 0
  IF (myid == master) THEN ! Only master processor read variables and send to MPI_COMM_WORLD only those needed

!
!   Local vars and arrays
!________________________________________________________________________________
  WRITE(*,'(a/)') '=== Define additional data for a new run ==='
!________________________________________________________________________________

!
!   Load field (GBS, marymima, etc...) parameters from input file
!
    NAMELIST /fields_in/ Bratio, B_dir, Bx,By,Bz,B_t,Rcurv, Bfield_type, Lr, Lz, E_basis, Efield_type, static_field, vert_reverse, &
                         GBS_mag, GBS_factor, GBS_bg, E_dir, E_mag1, E_mag2, deltarmod0, E_max, E_max2, rscale, kickrph, kickfreq1, kickfreq2, &
			 E_kick1a, E_kick1b, E_kick2a, E_kick2b, kickrphfile, readkicks, &
                         deltarmod, deltazmod0, deltazmod, rmin, rmax, radialbc, periodic_field, interpolate_type, &
			 numrpts, numzpts, staglength, stagdist, strmfstart,taumax, avgpotefile, GBSfile, fieldtimestep, phibeamnorm, &
			 firstframeout, sinecentdiff, sinmoder, sinmodez1, sinmodet1,sinmodez2, sinmodet2, sinmodet3, sinecentdiff
    READ(lu_in,fields_in)
    WRITE(*,fields_in)
!
!
!   Load injection parameters from input file
!
    NAMELIST /injection_in/ stagpart, stagbin,tracerdtfactor, inject_type_pos, inject_type_vel, inject_coord_cart, &
                            inject_coord_cyl, inject_coord_tilt, red_mass, &
			    XCM_injection, YCM_injection, RCM_injection, ZCM_injection, Xwidth_injection, Ywidth_injection, &
			    Rwidth_injection, Zwidth_injection, alpha0, beta0, V0, devalpha, devbeta, devV, all_bins_same_pot, rand_seed
    READ(lu_in,injection_in)
    WRITE(*,injection_in)
!
!
!   Load output parameters from input file
!
    NAMELIST /diagnose_in/ dt_out_0d, nsave_0d, dt_out_1d, nsave_1d, write_result_dp, save_cart, save_cyl, &
		           save_Efield, save_tilt, energy_diag
    READ(lu_in,diagnose_in)
    WRITE(*,diagnose_in)
!
!
!
!   Load numerics parameters from input file
!
    NAMELIST /numerics_in/ ODE_solver, ndim_tracer
    READ(lu_in,numerics_in)
    WRITE(*,numerics_in)
!

  END IF

! Send field parameters to other processor
  CALL MPI_BCAST(Bratio,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(B_dir,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Bx,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(By,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Bz,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(B_t,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rcurv,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Bfield_type,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Lr,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Lz,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_basis,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Efield_type,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(static_field,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(vert_reverse,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(GBS_mag,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(GBS_factor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(GBS_bg,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_dir,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_mag1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_mag2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(deltarmod0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_max,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_max2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rscale,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(kickrph,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(kickfreq1,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(kickfreq2,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_kick1a,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_kick1b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_kick2a,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(E_kick2b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(readkicks,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(deltarmod,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(deltazmod0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(deltazmod,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(radialbc,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(periodic_field,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(interpolate_type,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(numrpts,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(numzpts,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(staglength,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(stagdist,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(strmfstart,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(taumax,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(fieldtimestep,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(phibeamnorm,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(firstframeout,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinecentdiff,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmoder,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmodez1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmodet1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmodez2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmodet2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinmodet3,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sinecentdiff,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)

! Send injection parameters to other processors
  CALL MPI_BCAST(stagpart,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(stagbin,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(tracerdtfactor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inject_type_pos,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inject_type_vel,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inject_coord_cart,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inject_coord_cyl,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inject_coord_tilt,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(red_mass,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(XCM_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(YCM_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(RCM_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ZCM_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Xwidth_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Ywidth_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rwidth_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Zwidth_injection,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(alpha0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(beta0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(V0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(devalpha,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(devbeta,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(devV,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(all_bins_same_pot,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rand_seed,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)

! Send output parameters to other processors
  CALL MPI_BCAST(write_result_dp,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt_out_0d,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nsave_0d,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt_out_1d,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nsave_1d,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(energy_diag,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(save_cart,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(save_cyl,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(save_Efield,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(save_tilt,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)

! Send numeric parameters to other processors
  CALL MPI_BCAST(ODE_solver,10,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ndim_tracer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

END SUBROUTINE newrun
