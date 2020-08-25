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
#include "EESUPPORT.h"
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
      INTEGER i,j,k,l, bi,bj


      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C      ELSE
C--   restart from a pickup:
C         CALL DIAGNOSTICS_OB_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

#endif /* ALLOW_DIAGNOSTICS_OB */


      RETURN
      END
