SUBROUTINE diagnose(kstep)
!
!   Diagnostics
!
  USE basic
  USE futils
  USE hashtable
  USE tracers
  USE injection
  USE fields
  USE diagnostic
  USE numerics
  IMPLICIT none
  INCLUDE 'mpif.h'

!
  INTEGER, INTENT(in) :: kstep
!
!   Local vars and arrays
  INTEGER, PARAMETER :: BUFSIZE = 20
  CHARACTER(len=128) :: str, fname
  INTEGER            :: rank = 0, dims(2), ind
  

  INTEGER :: myid, master, numprocs, ierr
  CALL mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
  master = 0


 ! WRITE(*,*) 'Diagnose ',myid,time
!____________________________________________________________________________
!                   1.   Initial diagnostics
 IF (kstep .EQ. 0 .AND. myid == master) THEN
!
     WRITE(*,'(a)') '   Initial diagnostics'
!
!                       1.1   Initial run or when NEWRES set to .TRUE.
!
     IF( .NOT. nlres .OR. newres) THEN
	IF (write_result_dp) THEN
	  CALL creatf(resfile, fidres, real_prec='d')
	ELSE 
	  CALL creatf(resfile, fidres)
	END IF
        WRITE(*,'(3x,a,a)') TRIM(resfile), ' created'
!
!  Label the run
        IF( LEN_TRIM(label1).GT.0 ) CALL attach(fidres, "/", "label1", TRIM(label1))
        IF( LEN_TRIM(label2).GT.0 ) CALL attach(fidres, "/", "label2", TRIM(label2))
        IF( LEN_TRIM(label3).GT.0 ) CALL attach(fidres, "/", "label3", TRIM(label3))
        IF( LEN_TRIM(label4).GT.0 ) CALL attach(fidres, "/", "label4", TRIM(label4))
!
!  Job number
        jobnum = 0
!
!  Data group
        CALL creatg(fidres, "/data", "data")
        CALL creatg(fidres, "/data/var0d", "0d history arrays")
        CALL creatg(fidres, "/data/var1d", "1d profiles")
        CALL creatg(fidres, "/data/var2d", "2d profiles")
        CALL creatg(fidres, "/data/var3d", "3d profiles")
!
!  File group
        CALL creatg(fidres, "/files", "files")
        CALL attach(fidres, "/files",  "jobnum", jobnum)
!
! 
!  Create groups for saving 1d variables 

        CALL creatd(fidres, rank, dims, "/data/var1d/time", "1/Omega_{ci}")

	IF (energy_diag) THEN

	  CALL creatg(fidres, "/data/var1d/part_KE", "part_KE")
	  CALL attach(fidres, "/data/var1d/part_KE", "PlotOrder", 1)
	  CALL creatg(fidres, "/data/var1d/part_PE", "part_PE")
	  CALL attach(fidres, "/data/var1d/part_PE", "PlotOrder", 1)	  
	  CALL creatg(fidres, "/data/var1d/dphidt", "dphidt")
	  CALL attach(fidres, "/data/var1d/dphidt", "PlotOrder", 1)

	END IF

	IF (save_cart) THEN

	  CALL creatg(fidres, "/data/var1d/xcartcoord", "xcartcoord")
	  CALL attach(fidres, "/data/var1d/xcartcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/ycartcoord", "ycartcoord")
	  CALL attach(fidres, "/data/var1d/ycartcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/zcartcoord", "zcartcoord")
	  CALL attach(fidres, "/data/var1d/zcartcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velx_cart", "velx_cart")
	  CALL attach(fidres, "/data/var1d/velx_cart", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/vely_cart", "vely_cart")
	  CALL attach(fidres, "/data/var1d/vely_cart", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velz_cart", "velz_cart")
	  CALL attach(fidres, "/data/var1d/velz_cart", "PlotOrder", 1)

          IF (save_Efield) THEN
             CALL creatg(fidres, "/data/var1d/Ex_cart", "Ex_cart")
             CALL attach(fidres, "/data/var1d/Ex_cart", "PlotOrder", 1)
             
             CALL creatg(fidres, "/data/var1d/Ez_cart", "Ez_cart")
             CALL attach(fidres, "/data/var1d/Ez_cart", "PlotOrder", 1)
          END IF

	END IF
	
	IF (save_cyl) THEN

	  CALL creatg(fidres, "/data/var1d/rcoord", "rcoord")
	  CALL attach(fidres, "/data/var1d/rcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/phicoord", "phicoord")
	  CALL attach(fidres, "/data/var1d/phicoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/zcoord", "zcoord")
	  CALL attach(fidres, "/data/var1d/zcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velr_cyl", "velr_cyl")
	  CALL attach(fidres, "/data/var1d/velr_cyl", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velphi_cyl", "velphi_cyl")
	  CALL attach(fidres, "/data/var1d/velphi_cyl", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velz_cyl", "velz_cyl")
	  CALL attach(fidres, "/data/var1d/velz_cyl", "PlotOrder", 1)

	END IF

	IF (save_tilt) THEN

	  CALL creatg(fidres, "/data/var1d/rtiltcoord", "rtiltcoord")
	  CALL attach(fidres, "/data/var1d/rtiltcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/phitiltcoord", "phitiltcoord")
	  CALL attach(fidres, "/data/var1d/phitiltcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/ztiltcoord", "ztiltcoord")
	  CALL attach(fidres, "/data/var1d/ztiltcoord", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velr_tilt", "velr_tilt")
	  CALL attach(fidres, "/data/var1d/velr_tilt", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velphi_tilt", "velphi_tilt")
	  CALL attach(fidres, "/data/var1d/velphi_tilt", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/velz_tilt", "velz_tilt")
	  CALL attach(fidres, "/data/var1d/velz_tilt", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/driftrExB_tilt", "driftrExB_tilt")
	  CALL attach(fidres, "/data/var1d/driftrExB_tilt", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/driftzB_tilt", "driftzB_tilt")
	  CALL attach(fidres, "/data/var1d/driftzB_tilt", "PlotOrder", 1)

	  CALL creatg(fidres, "/data/var1d/driftzExB_tilt", "driftzExB_tilt")
	  CALL attach(fidres, "/data/var1d/driftzExB_tilt", "PlotOrder", 1)

          IF (save_Efield) THEN
             CALL creatg(fidres, "/data/var1d/Er_tilt", "Er_tilt")
             CALL attach(fidres, "/data/var1d/Er_tilt", "PlotOrder", 1)
             
             CALL creatg(fidres, "/data/var1d/Ez_tilt", "Ez_tilt")
             CALL attach(fidres, "/data/var1d/Ez_tilt", "PlotOrder", 1)
          END IF

	END IF

        
!                       1.2   Restart run
!
     ELSE
        CALL cp2bk(resfile)    ! backup previous result file
        CALL openf(resfile, fidres)
        WRITE(*,'(3x,a,a)') TRIM(resfile), ' open'
        CALL getatt(fidres, "/files",  "jobnum", jobnum)
        jobnum = jobnum+1
        WRITE(*,'(3x,a,i3)') "Current Job Number =", jobnum
        CALL attach(fidres, "/files",  "jobnum", jobnum)
     END IF
!
!  Add input namelist variables as attributes of /data/input
     WRITE(str,'(a,i2.2)') "/data/input.",jobnum
     CALL creatg(fidres, TRIM(str))
     CALL attach(fidres, TRIM(str), "job_time", job_time)
     CALL attach(fidres, TRIM(str), "extra_time", extra_time)
     CALL attach(fidres, TRIM(str), "dt", dt)
     CALL attach(fidres, TRIM(str), "tmax", tmax)
     CALL attach(fidres, TRIM(str), "nrun", nrun)
     CALL attach(fidres, TRIM(str), "nlres", nlres)
     CALL attach(fidres, TRIM(str), "nlsave", nlsave)
     CALL attach(fidres, TRIM(str), "newres", newres)

     CALL attach(fidres, TRIM(str), "Bratio", Bratio)
     CALL attach(fidres, TRIM(str), "Rcurv", Rcurv)
     CALL attach(fidres, TRIM(str), "Bfield_type", Bfield_type)
     CALL attach(fidres, TRIM(str), "Efield_type", Efield_type)
     CALL attach(fidres, TRIM(str), "Lr", Lr)
     CALL attach(fidres, TRIM(str), "Lz", Lz)
     CALL attach(fidres, TRIM(str), "kr", kr)
     CALL attach(fidres, TRIM(str), "kz", kz)
     CALL attach(fidres, TRIM(str), "periodic_field", periodic_field)
     CALL attach(fidres, TRIM(str), "radialbc", radialbc)
     CALL attach(fidres, TRIM(str), "deltarmod0", deltarmod0)
     CALL attach(fidres, TRIM(str), "deltazmod0", deltazmod0)
     CALL attach(fidres, TRIM(str), "deltarmod", deltarmod)
     CALL attach(fidres, TRIM(str), "deltazmod", deltazmod)
     CALL attach(fidres, TRIM(str), "rmin", rmin)
     CALL attach(fidres, TRIM(str), "rmax", rmax)
     CALL attach(fidres, TRIM(str), "numrpts", numrpts)
     CALL attach(fidres, TRIM(str), "numzpts", numzpts)
     CALL attach(fidres, TRIM(str), "numtpts", numtpts)
     CALL attach(fidres, TRIM(str), "GBSfile", GBSfile)
     CALL attach(fidres, TRIM(str), "avgpotefile", avgpotefile)
     CALL attach(fidres, TRIM(str), "strmfstart", strmfstart)
     CALL attach(fidres, TRIM(str), "taumax", taumax)
     CALL attach(fidres, TRIM(str), "fieldtimestep", fieldtimestep)

     CALL attach(fidres, TRIM(str), "ntracer", ntracer)
     CALL attach(fidres, TRIM(str), "inject_type_pos", inject_type_pos)     
     CALL attach(fidres, TRIM(str), "inject_type_vel", inject_type_vel)
     CALL attach(fidres, TRIM(str), "inject_coord_cart", inject_coord_cart)
     CALL attach(fidres, TRIM(str), "inject_coord_cyl", inject_coord_cyl)
     CALL attach(fidres, TRIM(str), "inject_coord_tilt", inject_coord_tilt)
     CALL attach(fidres, TRIM(str), "tracerdtfactor", tracerdtfactor)
     CALL attach(fidres, TRIM(str), "XCM_injection", XCM_injection)
     CALL attach(fidres, TRIM(str), "YCM_injection", YCM_injection)
     CALL attach(fidres, TRIM(str), "ZCM_injection", ZCM_injection)
     CALL attach(fidres, TRIM(str), "RCM_injection", YCM_injection)
     CALL attach(fidres, TRIM(str), "Xwidth_injection", Xwidth_injection)
     CALL attach(fidres, TRIM(str), "Ywidth_injection", Ywidth_injection)
     CALL attach(fidres, TRIM(str), "Zwidth_injection", Zwidth_injection)
     CALL attach(fidres, TRIM(str), "Rwidth_injection", Rwidth_injection)
     CALL attach(fidres, TRIM(str), "V0", V0)
     CALL attach(fidres, TRIM(str), "alpha0", alpha0)
     CALL attach(fidres, TRIM(str), "beta0", beta0)
     CALL attach(fidres, TRIM(str), "devV", devV)
     CALL attach(fidres, TRIM(str), "devalpha", devbeta)
     CALL attach(fidres, TRIM(str), "devbeta", devbeta)
     CALL attach(fidres, TRIM(str), "beta0", beta0)

     CALL attach(fidres, TRIM(str), "dt_out_0d", dt_out_0d)
     CALL attach(fidres, TRIM(str), "nsave_0d", nsave_0d)
     CALL attach(fidres, TRIM(str), "save_cart", save_cart)
     CALL attach(fidres, TRIM(str), "save_cyl", save_cyl)
     CALL attach(fidres, TRIM(str), "save_tilt", save_tilt)
     CALL attach(fidres, TRIM(str), "energy_diag", energy_diag)

     CALL attach(fidres, TRIM(str), "ODE_solver", ODE_solver)
     CALL attach(fidres, TRIM(str), "ndim_tracer", ndim_tracer)

!
!  Save STDIN of this run
     WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
     INQUIRE(unit=lu_in, name=fname)
     CALL putfile(fidres, TRIM(str), TRIM(fname))
!
!  Initialize buffers for 0d history arrays
     CALL htable_init(hbuf0, BUFSIZE)
     CALL set_htable_fileid(hbuf0, fidres, "/data/var0d")
! 

END IF
!___________________________end of initial diagnostics___________________________
!________________________________________________________________________________
!                   2.   Periodic diagnostics
 IF ( ((kstep .EQ. 0) .AND. (.NOT. nlres)) .OR. (kstep .GT. 0) ) THEN


!
!                       2.1   0d history arrays
!
     IF (myid == master) THEN
	CALL add_record(hbuf0, "time", "simulation time", time)
     END IF
!
!                       2.2   1d profiles
! 

    
     IF (write_out_nstep .AND. nsave_1d .NE. 0 ) THEN
        IF (MOD(cstep, nsave_1d) == 0) THEN
		CALL diagnose_1d
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        END IF
     ELSE IF ( ((.NOT.write_out_nstep).AND.((time.GT.(last_timeout_1d + dt_out_1d)) ) ) ) THEN      
            last_timeout_1d = time - MOD(time,dt_out_1d) 
            CALL diagnose_1d
     END IF
     

!                       2.3   2d profiles
!
!                       2.4   3d profiles
!
     IF (myid == master) THEN
	CALL htable_endstep(hbuf0)
     END IF
!________________________________________________________________________________
!                   3.   Final diagnostics
 ELSEIF (kstep .EQ. -1) THEN
!
!   Flush 0d history array buffers
     IF (myid == master) THEN
	CALL htable_hdf5_flush(hbuf0)
     END IF


!
!   Close all diagnostic files
     IF (myid == master) THEN
	CALL closef(fidres)
     END IF

!________________________________________________________________________________
  END IF
!
END SUBROUTINE diagnose
