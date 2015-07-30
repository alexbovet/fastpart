MODULE fields
!
!   Module for fields parameters
!
  IMPLICIT none
!
  ! ratio Bvertical/Btoroidal, radius of curvature
  DOUBLE PRECISION  :: Bratio=0., Rcurv=1., Bx=0.,By=0.,Bz=0.,B_t=1.     
  DOUBLE PRECISION  :: theta=0.,sineoftheta=0.,cosineoftheta=1.! Angle between the B field and the floor
  DOUBLE PRECISION  :: Lr=0.,Lz=0.,kr,kz,dr,dz, sinmoder,sinmodez1, sinmodet1,sinmodez2, sinmodet2, sinmodet3
  INTEGER, DIMENSION(1000,323) :: kicksteps
  INTEGER           :: numrpts=128, numzpts=128
  INTEGER           :: numtpts,staglength,stagdist,tA,tB
  LOGICAL           :: periodic_field=.FALSE., vert_reverse=.FALSE., static_field=.FALSE., sinecentdiff=.TRUE.
  CHARACTER(LEN=10) :: E_basis='tilt', Efield_type='sinusoid', E_dir='r', Bfield_type='torpex', rot_type, B_dir
  CHARACTER(LEN=20) :: interpolate_type='bilinear', radialbc='loss'
  DOUBLE PRECISION  :: rscale=0., E_max=0., E_max2=0., E_mag1=0., E_mag2=0., es_pot=0., E_kick1a=0., E_kick1b=0., &
                       E_kick2a=0., E_kick2b=0., old_Ex_cart=0., old_Ez_cart=0.,   old_Er_tilt=0., old_Ez_tilt=0.
  ! The sign is reversed on the electric field magnitude E_mag
  DOUBLE PRECISION  :: deltarmod0=0., deltarmod=0., deltazmod0=0., deltazmod=0., rmin=0., rmax=100.
  DOUBLE PRECISION  :: fieldtimestep, phibeamnorm=0.00235, GBS_factor=1., GBS_mag=1., GBS_bg=1.
  CHARACTER(len=64) :: avgpotefile, kickrphfile, GBSfile ! results file from GBS - this should be user specified
  INTEGER           :: fidavgpote,fidkickrph, fidGBS, strmfstart=0, taumax=7500, kickstep1=1, kickfreq1=0, kickstep2=1, kickfreq2=0
  LOGICAL           :: firstframeout=.FALSE., kickrph=.TRUE., readkicks=.FALSE.
  !fidGBS is a placeholder, strmfstart is the last 3 digits of the starting strmf file	       
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: tracer_es_pot, es_pot_old, dphidt
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: rtiltgrid, ztiltgrid
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: es_potgrid, dpotdrtiltgrid, dpotdztiltgrid, d2potdrdztiltgrid
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: es_potgrid_t, es_potgrid_fluct, & 
                                                     dpotdrtiltgrid_t, dpotdztiltgrid_t, d2potdrdztiltgrid_t
  ! These are the fields at the tracer positions in various coordinate systems
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: tracer_E_cyl, tracer_E_cart, tracer_E_tilt 
  ! Use same perturbation for all bins or not
  LOGICAL :: all_bins_same_pot = .FALSE.
!
END MODULE fields

! 
