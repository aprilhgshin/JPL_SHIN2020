C
C     *===============================================================*
C     | USER_INPUT.h Contains user input for open boundary extraction.
C     *===============================================================*
C
C     nOB_mask  :: Number of open boundary masks.
C     nOB_fld  :: Number of fields.

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

      INTEGER, PARAMETER :: nFldOpt = 8
C     Please do not change unless adding on state logical variables
C---------------------------------------

      INTEGER, PARAMETER :: nOB_mask =  1
      INTEGER, PARAMETER :: nOB_fld = 1

      LOGICAL etaN_ob
      LOGICAL uVel_ob
      LOGICAL vVel_ob
      LOGICAL theta_ob
      LOGICAL salt_ob
      LOGICAL gU_ob
      LOGICAL gV_ob
      LOGICAL etaH_ob

      PARAMETER (
     &            etaN_ob = .FALSE.,
     &            uVel_ob = .FALSE.,
     &            vVel_ob = .FALSE.,
     &           theta_ob = .FALSE.,
     &            salt_ob = .TRUE.,
     &              gU_ob = .FALSE.,
     &              gV_ob = .FALSE.,
     &            etaH_ob = .FALSE.)
