MODULE diagnostic
!
!   Module for diagnostic parameters
!
!
  LOGICAL             				:: write_out_nstep=.TRUE., write_result_dp=.TRUE.
  DOUBLE PRECISION    				:: dt_out_0d=1., last_timeout_0d=0.
  DOUBLE PRECISION    				:: dt_out_1d=1., last_timeout_1d=0.
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: part_PE, part_KE, driftrExB_tilt, driftzB_tilt, driftzExB_tilt

  INTEGER             				:: nsave_0d=1., nsave_1d=1.

  LOGICAL            				:: save_cyl=.FALSE., save_cart=.FALSE., save_tilt=.TRUE.
  LOGICAL					:: save_Efield = .TRUE., energy_diag=.FALSE.
!  
END MODULE diagnostic

