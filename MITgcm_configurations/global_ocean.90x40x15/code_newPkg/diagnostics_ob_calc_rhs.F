#include "DIAGNOSTICS_OB_OPTIONS.h"

CBOP
C !ROUTINE: DIAGNOSTICS_OB_CALC_RHS

C !INTERFACE: ==========================================================
      SUBROUTINE DIAGNOSTICS_OB_CALC_RHS(
     I                     myTime, myIter, myThid )

C !DESCRIPTION:
C Calculate custom tendency terms outside k-loop in DO_OCEANIC_PHYS

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
c#include "PARAMS.h"
c#include "DIAGNOSTICS_OB.h"

C !INPUT PARAMETERS: ===================================================
C myTime    :: Current time in simulation
C myIter    :: Current time-step number
C myThid    :: my Thread Id number
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================

C !LOCAL VARIABLES: ====================================================
C bi,bj     :: Tile indices
      INTEGER bi, bj
CEOP

C--   Loops on tile indices bi,bj
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C--   Calculate custom tendency terms ...

C--   end bi,bj loops.
       ENDDO
      ENDDO

      RETURN
      END
