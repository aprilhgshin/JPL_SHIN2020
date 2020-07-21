#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_INIT_FIXED

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize DIAGNOSTICS_OB variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===

#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL DIAGNOSTICS_OB_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_OB_INIT( myThid )
      ENDIF
#endif

C      CALL OB_EXTRACT( myThid )
      CALL OB_EXTRACT_LOOKUP( myThid )

      RETURN
      END
