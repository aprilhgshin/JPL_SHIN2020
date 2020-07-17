#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_CHECK

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_CHECK( myThid )

C     !DESCRIPTION:
C     Check dependances with other packages

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT PARAMETERS:
C     myThid   :: my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS_OB
C     !LOCAL VARIABLES:
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'DIAGNOSTICS_OB_CHECK: #define DIAGNOSTICS_OB'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , 1)

C--  Print out some kee parameters :
       CALL WRITE_0D_L( diagOB_StaV_Cgrid, INDEX_NONE,
     &  'diagOB_StaV_Cgrid =', ' /* state vector on C-grid */')
       CALL WRITE_0D_L( diagOB_Tend_Cgrid, INDEX_NONE,
     &  'diagOB_Tend_Cgrid =', ' /* vector tendency on C-grid */')
       CALL WRITE_0D_L( diagOB_applyTendT, INDEX_NONE,
     &  'diagOB_applyTendT =',
     &  ' /* apply diagOB tendency to Temperature Eq.*/')
       CALL WRITE_0D_L( diagOB_applyTendS, INDEX_NONE,
     &  'diagOB_applyTendS =',
     &  ' /* apply diagOB tendency to Salinity Eq.*/')
       CALL WRITE_0D_L( diagOB_applyTendU, INDEX_NONE,
     &  'diagOB_applyTendU =',
     &  ' /* apply diagOB tendency to U momentum Eq.*/')
       CALL WRITE_0D_L( diagOB_applyTendV, INDEX_NONE,
     &  'diagOB_applyTendV =',
     &  ' /* apply diagOB tendency to V momentum Eq.*/')
C-   add some other as example:
       CALL WRITE_0D_I( diagOB_index1, INDEX_NONE, ' diagOB_index1 =',
     &  ' /* user defined parameter */')
       CALL WRITE_0D_RL( diagOB_param1, INDEX_NONE, ' diagOB_param1 =',
     &  ' /* user defined parameter */')
       CALL WRITE_0D_C( diagOB_string1,-1,INDEX_NONE,'diagOB_string1 =',
     &  ' /* user defined parameter */')

C--  Check parameters:

C-     Put stops here if some flags are wrongly chosen.
C-     For example, require this package to use tempStepping
       IF ( .NOT.tempStepping .AND. diagOB_applyTendT ) THEN
        WRITE(msgBuf,'(A)') 'An example requirement is:'
        CALL PRINT_ERROR( msgBuf , 1)
        WRITE(msgBuf,'(A)') 'diagOB_applyTendT needs tempStepping=.true.'
        CALL PRINT_ERROR( msgBuf , 1)
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R DIAGNOSTICS_OB_CHECK'
       ENDIF

      _END_MASTER(myThid)

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
