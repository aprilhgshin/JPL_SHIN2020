#include "DIAGNOSTICS_OB_OPTIONS.h"
C--  File DIAGNOSTICS_OB_tendency_apply.F: Routines to apply DIAGNOSTICS_OB tendencies
C--   Contents
C--   o DIAGNOSTICS_OB_TENDENCY_APPLY_U
C--   o DIAGNOSTICS_OB_TENDENCY_APPLY_V
C--   o DIAGNOSTICS_OB_TENDENCY_APPLY_T
C--   o DIAGNOSTICS_OB_TENDENCY_APPLY_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_TENDENCY_APPLY_U

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_TENDENCY_APPLY_U(
     U                     gU_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add DIAGNOSTICS_OB tendency terms to U tendency.  Routine works for one
C     level at a time. Assumes that U and V tendencies are on C-Grid or
C     on A-grid.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gU_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gU_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef DIAGNOSTICS_OB_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( diagOB_applyTendU .AND. diagOB_Tend_Cgrid ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gU_arr(i,j) = gU_arr(i,j)
     &     + maskW(i,j,k,bi,bj) * diagOB_TendVelU(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ELSEIF ( diagOB_applyTendU ) THEN
C--   Assume Agrid position:
        DO j=jMin,jMax
         DO i=iMin,iMax
          gU_arr(i,j) = gU_arr(i,j)
     &     + maskW(i,j,k,bi,bj)*( diagOB_TendVelU(i-1,j,k,bi,bj)
     &                           +diagOB_TendVelU( i ,j,k,bi,bj)
     &                          )*0.5 _d 0
         ENDDO
        ENDDO
      ENDIF

#endif /* DIAGNOSTICS_OB_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_TENDENCY_APPLY_V

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_TENDENCY_APPLY_V(
     U                     gV_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add DIAGNOSTICS_OB tendency terms to V tendency.  Routine works for one
C     level at a time. Assumes that U and V tendencies are on C-Grid or
C     on A-grid.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gV_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gV_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef DIAGNOSTICS_OB_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( diagOB_applyTendV .AND. diagOB_Tend_Cgrid ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gV_arr(i,j) = gV_arr(i,j)
     &     + maskS(i,j,k,bi,bj) * diagOB_TendVelV(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ELSEIF ( diagOB_applyTendV ) THEN
C--   Assume Agrid position:
        DO j=jMin,jMax
         DO i=iMin,iMax
          gV_arr(i,j) = gV_arr(i,j)
     &     + maskS(i,j,k,bi,bj)*( diagOB_TendVelV(i,j-1,k,bi,bj)
     &                           +diagOB_TendVelV(i, j ,k,bi,bj)
     &                          )*0.5 _d 0
         ENDDO
        ENDDO
      ENDIF

#endif /* DIAGNOSTICS_OB_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_TENDENCY_APPLY_T(
     U                     gT_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add DIAGNOSTICS_OB tendency terms to T (theta) tendency.  Routine works
C     for one level at a time

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef DIAGNOSTICS_OB_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( diagOB_applyTendT ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gT_arr(i,j) = gT_arr(i,j)
     &     + maskC(i,j,k,bi,bj) * diagOB_TendScal1(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ENDIF

#endif /* DIAGNOSTICS_OB_TENDENCY */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_TENDENCY_APPLY_S(
     U                     gS_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add DIAGNOSTICS_OB tendency terms to S tendency.  Routine works for one
C     level at a time

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
c#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef DIAGNOSTICS_OB_TENDENCY
C     !LOCAL VARIABLES:
      INTEGER i, j

      IF ( diagOB_applyTendS ) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          gS_arr(i,j) = gS_arr(i,j)
     &     + maskC(i,j,k,bi,bj) * diagOB_TendScal2(i,j,k,bi,bj)
         ENDDO
        ENDDO
      ENDIF

#endif /* DIAGNOSTICS_OB_TENDENCY */

      RETURN
      END
