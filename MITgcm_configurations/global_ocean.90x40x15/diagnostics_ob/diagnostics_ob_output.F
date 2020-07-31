#include "DIAG_OPTIONS.h"

      SUBROUTINE DIAGNOSTICS_OB_OUTPUT( myTime, myIter, myThid )

        IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myTime :: my time in simulation ( s )
C     myIter :: my Iteration number
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     == Local variables ==
      CHARACTER*10 suff
      INTEGER j, k, bi, bj, fld, sz
      LOGICAL DIFFERENT_MULTIPLE

C----------------------------------------------------------------
C Diagnostics_ob time averaging
C----------------------------------------------------------------
C     Set field values to subField array
      print *, "DIAGNOSTICS_OB_OUTPUT CALLED"
      print *, "OB_CONTROL: mpiMYId:", mpiMYId
      print *, "DIAGNOSTICS_OB_OUTPUT time", myTime

      sz = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy
      CALL PLOT_FIELD_XYRL( salt,
     &                   'DIAGNOSTICS_OB_OUTPUT salt', -1, myThid )
      print *, "DIAGNOSTICS_OB_OUTPUT: Printing salt"
      CALL PRINT_FLOAT_ARR(salt(:,:,1,:,:),sz)

       DO fld=1, nOB_fld
         DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
             DO j=1-OLy,sNy+OLy
               DO k=1-OLx,sNx+OLx
                 subField(fld,k,j,bi,bj) = salt(k,j,1,bi,bj)
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       print *, "DIAGNOSTICS_OB_OUTPUT: Printing subField"
       CALL PRINT_FLOAT_ARR(subField(1,:,:,:,:),sz)

C       _RL fld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      END
