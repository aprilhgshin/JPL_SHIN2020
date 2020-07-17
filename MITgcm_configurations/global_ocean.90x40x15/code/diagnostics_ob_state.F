#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_DIAGNOSTICS_STATE

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_STATE(
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Fill-in diagnostics array for DIAGNOSTICS_OB state variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C  put following include within #ifdef in case this file is copied to local
C  code dir but pkg/DIAGNOSTICS_OB is not compiled (missing file DIAGNOSTICS_OB.h).
#ifdef ALLOW_DIAGNOSTICS_OB
#include "DIAGNOSTICS_OB.h"
#endif

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS_OB
C     !LOCAL VARIABLES:

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

#ifdef DIAGNOSTICS_OB_3D_STATE
         CALL DIAGNOSTICS_FILL( diagOB_StatScal1,'diagOBSta1',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( diagOB_StatScal2,'diagOBSta2',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( diagOB_StatVelU, 'diagOBStaU',
     &                          0,Nr, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( diagOB_StatVelV, 'diagOBStaV',
     &                          0,Nr, 0, 1, 1, myThid )
#endif /* DIAGNOSTICS_OB_3D_STATE */

#ifdef DIAGNOSTICS_OB_2D_STATE
         CALL DIAGNOSTICS_FILL( diagOB_Surf1 ,   'diagOBSur1',
     &                          0, 1, 0, 1, 1, myThid )
         CALL DIAGNOSTICS_FILL( diagOB_Surf2 ,   'diagOBSur2',
     &                          0, 1, 0, 1, 1, myThid )
#endif /* DIAGNOSTICS_OB_2D_STATE */

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
