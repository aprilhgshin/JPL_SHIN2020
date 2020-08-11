C
C     *===============================================================*
C     | USER_INPUT.h Contains user input for open boundary extraction.
C     *===============================================================*
C
C     nOB_mask  :: Number of open boundaries.
C     nOB_fld2D  :: Number of 2D fields to output.
C     nOB_fld3D  :: Number of 3D fields to output.
C     nOB_fld  :: Total number of fields to output.

C     State Logical Variables: Indicate True/False to/to not output field
C     etaN_ob  :: free-surface r-anomaly (r unit) at current time level
C     uVel_ob  :: zonal velocity (m/s, i=1 held at western face)
C     vVel_ob  :: meridional velocity (m/s, j=1 held at southern face)
C     theta_ob :: potential temperature (oC, held at pressure/tracer point)
C     salt_ob  :: salinity (ppt, held at pressure/tracer point)
C     gU_ob, gV_ob :: Time tendencies at current and previous time levels.
C     etaH_ob  :: surface r-anomaly, advanced in time consistently
C              with 2.D flow divergence (Exact-Conservation):
C                etaH^n+1 = etaH^n - delta_t*Div.(H^n U^n+1)

C---------------------------------------
C     nFldOpt :: Number of field options (i.e. state logical variables)
C     n3D_fld :: Number of 3D fields
C     n2D_fld :: Number of 2D fields
C     fld3D_depth :: int array of depth of only 3D fields

      INTEGER, PARAMETER :: nFldOpt = 15
      INTEGER, PARAMETER :: n3D_fld = 6
      INTEGER, PARAMETER :: n2D_fld = nFldOpt - n3D_fld

C     Please do not change unless adding on additional fields

C---------------------------------------

      INTEGER, PARAMETER :: nOB_mask =  1
      INTEGER, PARAMETER :: nOB_fld2D = 0
      INTEGER, PARAMETER :: nOB_fld3D = 2
      INTEGER, PARAMETER :: nOB_fld = nOB_fld2D + nOB_fld3D

C ----2D Fields (x,y):
      LOGICAL area_ob
      LOGICAL heff_ob
      LOGICAL hsnow_ob
      LOGICAL hsalt_ob
      LOGICAL uice_ob
      LOGICAL vice_ob
      LOGICAL etaN_ob
      LOGICAL etaH_ob
C ----2D Fields (y,z)
      LOGICAL obnw_ob

C ----3D Fields (x,y,z):
      LOGICAL uVel_ob
      LOGICAL vVel_ob
      LOGICAL theta_ob
      LOGICAL salt_ob
      LOGICAL gU_ob
      LOGICAL gV_ob


      PARAMETER (
     &            area_ob = .FALSE.,
     &            heff_ob = .FALSE.,
     &           hsnow_ob = .FALSE.,
     &           hsalt_ob = .FALSE.,
     &            uice_ob = .FALSE.,
     &            vice_ob = .FALSE.,
     &            etaN_ob = .FALSE.,
     &            etaH_ob = .FALSE.,
     &            obnw_ob = .FALSE.,

     &            uVel_ob = .FALSE.,
     &            vVel_ob = .FALSE.,
     &           theta_ob = .TRUE.,
     &            salt_ob = .TRUE.,
     &              gU_ob = .FALSE.,
     &              gV_ob = .FALSE.)
