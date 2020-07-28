#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_INIT_VARIA

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize DIAGNOSTICS_OB variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS_OB

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER i,j,k,bi,bj

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
        DO k=1,Nr
         DO j=1-Oly,sNy+OLy
          DO i=1-Olx,sNx+Olx
#ifdef DIAGNOSTICS_OB_3D_STATE
           diagOB_StatScal1(i,j,k,bi,bj) = 0. _d 0
           diagOB_StatScal2(i,j,k,bi,bj) = 0. _d 0
           diagOB_StatVelU (i,j,k,bi,bj) = 0. _d 0
           diagOB_StatVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* DIAGNOSTICS_OB_3D_STATE */
#ifdef DIAGNOSTICS_OB_TENDENCY
           diagOB_TendScal1(i,j,k,bi,bj) = 0. _d 0
           diagOB_TendScal2(i,j,k,bi,bj) = 0. _d 0
           diagOB_TendVelU (i,j,k,bi,bj) = 0. _d 0
           diagOB_TendVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* DIAGNOSTICS_OB_TENDENCY */
          ENDDO
         ENDDO
        ENDDO

C     2-dim.
#ifdef DIAGNOSTICS_OB_2D_STATE
        DO j=1-Oly,sNy+OLy
         DO i=1-Olx,sNx+Olx
          diagOB_Surf1(i,j,bi,bj) = 0. _d 0
          diagOB_Surf2(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* DIAGNOSTICS_OB_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef DIAGNOSTICS_OB_3D_STATE
       IF ( diagOB_Scal1File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_Scal1File, ' ',
     &                         diagOB_StatScal1, 0, myThid )
         CALL EXCH_3D_RL( diagOB_StatScal1, Nr, myThid )
       ENDIF
       IF ( diagOB_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_Scal2File, ' ',
     &                         diagOB_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( diagOB_StatScal2, Nr, myThid )
       ENDIF

       IF ( diagOB_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_VelUFile, ' ',
     &                         diagOB_StatVelU, 0, myThid )
       ENDIF
       IF ( diagOB_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_VelVFile, ' ',
     &                         diagOB_StatVelV, 0, myThid )
       ENDIF
       IF ( diagOB_VelUFile .NE. ' ' .OR.
     &      diagOB_VelVFile .NE. ' ' ) THEN
        IF ( diagOB_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( diagOB_StatVelU, diagOB_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( diagOB_StatVelU, diagOB_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* DIAGNOSTICS_OB_3D_STATE */

#ifdef DIAGNOSTICS_OB_2D_STATE
       IF ( diagOB_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( diagOB_Surf1File,' ',
     &                        diagOB_Surf1, 0, myThid )
         CALL EXCH_XY_RL( diagOB_Surf1, myThid )
       ENDIF
       IF ( diagOB_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( diagOB_Surf2File,' ',
     &                        diagOB_Surf2, 0, myThid )
         CALL EXCH_XY_RL( diagOB_Surf2, myThid )
       ENDIF
#endif /* DIAGNOSTICS_OB_2D_STATE */

      ELSE
C--   restart from a pickup:
         CALL DIAGNOSTICS_OB_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
