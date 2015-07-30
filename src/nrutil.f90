!module reals
!  INTEGER, parameter, public :: dp = kind(0.0D0), sp = kind(0.0E0)
!  INTEGER, parameter, public :: pr = dp
!END module reals
!
MODULE nrutil ! this module is not equivalent to the module of the same name in the NR book

!  use reals

  IMPLICIT none

CONTAINS
  
  FUNCTION reallocate_rv(p,n)
    ! Reallocate a pointer to a new size, preserving its previous contents.
    DOUBLE PRECISION, DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold, ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv
  
  FUNCTION iminloc(arr) ! index of minloc on an array
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin = minloc(arr(:))
    iminloc = imin(1)
  END function iminloc

  FUNCTION assert_eq2(n1,n2,string)  ! report and die if integers not all equal
    CHARACTER(LEN=*), INTENT(in) :: string
    INTEGER, INTENT(in) :: n1,n2
    INTEGER :: assert_eq2
      
    IF (n1 == n2) THEN
     assert_eq2=n1
    ELSE
      write(*,*) 'nrerror: an assert_eq failed with this tag:', &
        string
      STOP 'program terminated by assert_eq2'
    END IF
  END function assert_eq2

  SUBROUTINE nrerror(string) ! report a message, then die
    CHARACTER(LEN=*), INTENT(in) :: string
    write(*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: d1,d2
    DOUBLE PRECISION, DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
    DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: c
      !Given arrays y, y1, y2, and y12, each of length 4, containing the function, gradients, and
      !cross derivative at the four grid points of a rectangular grid cell (numbered counterclockwise
      !from the lower left), and given d1 and d2, the length of the grid cell in the 1- and 2-
      !directions, this routine returns the 4×4 table c that is used by routine bcuint for bicubic
      !interpolation.
    DOUBLE PRECISION, DIMENSION(16) :: x,c1
    DOUBLE PRECISION :: xx
!     DOUBLE PRECISION, DIMENSION(4,4) :: c1
    DOUBLE PRECISION, DIMENSION(16,16) :: wt
    INTEGER :: i,j,k,l
    DATA wt /1., 0.,-3., 2., 0., 0., 0., 0.,-3., 0., 9.,-6., 2., 0.,-6., 4., &
	     0., 0., 0., 0., 0., 0., 0., 0., 3., 0.,-9., 6.,-2., 0., 6.,-4., &
	     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 9.,-6., 0., 0.,-6., 4., &
             0., 0., 3.,-2., 0., 0., 0., 0., 0., 0.,-9., 6., 0., 0., 6.,-4., &
             0., 0., 0., 0., 1., 0.,-3., 2.,-2., 0., 6.,-4., 1., 0.,-3., 2., &
             0., 0., 0., 0., 0., 0., 0., 0.,-1., 0., 3.,-2., 1., 0.,-3., 2., &
             0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-3., 2., 0., 0., 3.,-2., &
             0., 0., 0., 0., 0., 0., 3.,-2., 0., 0.,-6., 4., 0., 0., 3.,-2., &
             0., 1.,-2., 1., 0., 0., 0., 0., 0.,-3., 6.,-3., 0., 2.,-4., 2., &
             0., 0., 0., 0., 0., 0., 0., 0., 0., 3.,-6., 3., 0.,-2., 4.,-2., &
             0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-3., 3., 0., 0., 2.,-2., &
             0., 0.,-1., 1., 0., 0., 0., 0., 0., 0., 3.,-3., 0., 0.,-2., 2., &
             0., 0., 0., 0., 0., 1.,-2., 1., 0.,-2., 4.,-2., 0., 1.,-2., 1., &
             0., 0., 0., 0., 0., 0., 0., 0., 0.,-1., 2.,-1., 0., 1.,-2., 1., &
             0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,-1., 0., 0.,-1., 1., &
             0., 0., 0., 0., 0., 0.,-1., 1., 0., 0., 2.,-2., 0., 0.,-1., 1./

    x(1:4)=y !Pack a temporary vector x.
    x(5:8)=y1*d1
    x(9:12)=y2*d2
    x(13:16)=y12*d1*d2
    x=matmul(wt,x) !Matrix multiply by the stored table.
    c=reshape(x,(/4,4/),order=(/2,1/)) !Unpack the result into the output table.

  END SUBROUTINE bcucof

END MODULE nrutil
