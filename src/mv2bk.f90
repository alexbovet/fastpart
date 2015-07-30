SUBROUTINE mv2bk(file)
!
!   Move file to a backup file if it exists
!
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(in) :: file
!
  LOGICAL            :: ex
  CHARACTER(len=32)  :: bkfile
  CHARACTER(len=128) :: cmd
!
  INQUIRE(file=file, exist=ex)
  IF( ex ) THEN
     bkfile = TRIM(file) // '_old'
     CALL move_file(file, LEN_TRIM(file), bkfile, LEN_TRIM(bkfile))
     WRITE(*,'(a)')  "   Moving existing " // TRIM(file) // " to " &
          &     // TRIM(bkfile)
  END IF
END SUBROUTINE mv2bk
!
SUBROUTINE cp2bk(file)
!
!   Copy file to a backup file if it exists
!
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(in) :: file
  LOGICAL :: ex
  CHARACTER(len=32) :: bkfile
  CHARACTER(len=128) :: cmd
!
  INQUIRE(file=file, exist=ex)
  IF( ex ) THEN
     bkfile = TRIM(file) // '_old'
     CALL copy_file(file, LEN_TRIM(file), bkfile, LEN_TRIM(bkfile))
     WRITE(*,'(a)')  "   Copying existing " // TRIM(file) // " to " &
          &     // TRIM(bkfile)
  END IF
END SUBROUTINE cp2bk
