MODULE constants
!
!   Define some constants
!
  INTEGER, PARAMETER       :: db     = SELECTED_REAL_KIND(8)
!
  DOUBLE PRECISION, PARAMETER :: vlight = 299792458.0_db         ! c
  DOUBLE PRECISION, PARAMETER :: vacimp = 376.73031346177066_db  ! \mu_0*c
  DOUBLE PRECISION, PARAMETER :: eev    = 510998.89613320108_db  ! m_e*c^2/e
  DOUBLE PRECISION, PARAMETER :: pi     = 3.1415926535897931_db
END MODULE constants
