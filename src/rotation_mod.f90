MODULE rotation
!
!   Module for rotations of coordinates including
!   1) cylindrical to tilted (magnetic plane) cylindrical
!   2) tilted cylindrical to cartesian
!   3) cartesian to cylindrical  
!
!
! !  DOUBLE PRECISION ::  


!  r = sqrt(x*x+y*y)
!       rinv = 1.0/r
!       Er = Ertilt
!       Ez = Eztilt*sintheta
!       Ex = Er*x*rinv + Ez*y*rinv
!       Ey = Er*y*rinv - Ez*x*rinv
!  
END MODULE rotation

