#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: DIAGNOSTICS_OB_MNC_INIT

C     !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_MNC_INIT( myThid )

C     !DESCRIPTION:
C     Initialize MNC lookup tables for DIAGNOSTICS_OB.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS_OB
#ifdef ALLOW_MNC

      IF (useMNC) THEN

        CALL MNC_CW_ADD_VNAME('diagOB_StatScal1',
     &      'Cen_xy_Hn__L__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatScal1',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatScal1','description',
     &      'scalar state variable diagOB_StatScal1',myThid)

        CALL MNC_CW_ADD_VNAME('diagOB_StatScal2',
     &      'Cen_xy_Hn__L__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatScal2',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatScal2','description',
     &      'scalar state variable diagOB_StatScal2',myThid)

        CALL MNC_CW_ADD_VNAME('diagOB_StatVelU',
     &      'U_xy_Hn__C__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatVelU',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatVelU','description',
     &       'uVel-type state variable diagOB_StatVelU',myThid)

        CALL MNC_CW_ADD_VNAME('diagOB_StatVelV',
     &      'V_xy_Hn__C__t', 4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatVelV',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_StatVelV','description',
     &       'vVel-type state variable diagOB_StatVelV',myThid)

        CALL MNC_CW_ADD_VNAME('diagOB_Surf1',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_Surf1',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_Surf1','description',
     &       'Surface scalar field diagOB_Surf1',myThid)

        CALL MNC_CW_ADD_VNAME('diagOB_Surf2',
     &      'Cen_xy_Hn__-__t', 3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_Surf2',
     &      'units','[unit]', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('diagOB_Surf2','description',
     &       'Surface scalar field diagOB_Surf2',myThid)

      ENDIF

#endif /* ALLOW_MNC */
#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
