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
C  code dir but pkg/diagnostics_ob is not compiled (missing file DIAGNOSTICS_OB.h).
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

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
