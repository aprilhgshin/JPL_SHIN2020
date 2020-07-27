#ifdef ALLOW_DIAGNOSTICS_OB
#include "USER_INPUT.h"

C     Package flag
      LOGICAL diagOB_MNC
      LOGICAL diagOB_MDSIO
      COMMON /DIAGOB_PACKAGE/
     &                     diagOB_MNC, diagOB_MDSIO

C     DIAGOB parameters
      LOGICAL diagOB_StaV_Cgrid
      LOGICAL diagOB_Tend_Cgrid
      LOGICAL diagOB_applyTendT
      LOGICAL diagOB_applyTendS
      LOGICAL diagOB_applyTendU
      LOGICAL diagOB_applyTendV

C-    additional parameters:
      LOGICAL diagOB_doSwitch1
      LOGICAL diagOB_doSwitch2
      INTEGER diagOB_index1
      INTEGER diagOB_index2
      _RL diagOB_param1
      _RL diagOB_param2
      CHARACTER*(MAX_LEN_FNAM) diagOB_string1
      CHARACTER*(MAX_LEN_FNAM) diagOB_string2

C-    file names for initial conditions:
      CHARACTER*(MAX_LEN_FNAM) diagOB_Scal1File
      CHARACTER*(MAX_LEN_FNAM) diagOB_Scal2File
      CHARACTER*(MAX_LEN_FNAM) diagOB_VelUFile
      CHARACTER*(MAX_LEN_FNAM) diagOB_VelVFile
      CHARACTER*(MAX_LEN_FNAM) diagOB_Surf1File
      CHARACTER*(MAX_LEN_FNAM) diagOB_Surf2File

      COMMON /DIAGOB_PARAMS_L/
     &       diagOB_StaV_Cgrid, diagOB_Tend_Cgrid,
     &       diagOB_applyTendT, diagOB_applyTendS,
     &       diagOB_applyTendU, diagOB_applyTendV,
     &       diagOB_doSwitch1, diagOB_doSwitch2
      COMMON /diagOB_PARAMS_I/ diagOB_index1, diagOB_index2
      COMMON /diagOB_PARAMS_R/ diagOB_param1, diagOB_param2
      COMMON /diagOB_PARAMS_C/ diagOB_string1, diagOB_string2,
     &       diagOB_Scal1File, diagOB_Scal2File,
     &       diagOB_VelUFile,  diagOB_VelVFile,
     &       diagOB_Surf1File, diagOB_Surf2File

C     ================== Global Variables for open boundary ====================
      CHARACTER*(40) ob_fileName
      CHARACTER*(40) t_fileName
      CHARACTER*8 ob_fldNames(nOB_fld)
      CHARACTER*80 ob_fnames(nOB_mask + nOB_fld)

      _RL ob_subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL bath_subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL temp_subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)

      _RL subBathOnMask(sNx + sNy + 1)
      _RL subTempOnMask(sNx + sNy + 1)
      _RL sub_global_indices(sNx + sNy + 1)
      INTEGER lookup_table(nOB_mask, Ny*Nx)
      _RL global_ob((sNy+sNx)*(nPx*nPy))

      INTEGER avgPeriod_ob
      INTEGER deltaT_ob
      INTEGER totPhase_ob
C     ==========================================================================

      COMMON /DIAG_OB_EXTRACT/
     &     ob_fldNames, ob_fnames,
     &     ob_subMask,lookup_table, ob_fileName,
     &     global_ob, sub_global_indices,t_fileName,
     &     temp_subMask, bath_subMask,
     &     subBathOnMask, subTempOnMask,
     &     avgPeriod_ob, deltaT_ob, totPhase_ob

C      COMMON / DIAG_OB_EXTRACT_R /
C     &     ob_subMask,
C     &     global_ob, sub_global_indices,
C     &     temp_subMask, bath_subMask,
C     &     subBathOnMask, subTempOnMask
C      COMMON / DIAG_OB_EXTRACT_I /
C     &     lookup_table
C      COMMON / DIAG_OB_EXTRACT_C /
C     &     ob_fileName, t_fileName

#ifdef DIAGNOSTICS_OB_3D_STATE
C     DIAGOB 3-dim. fields
      _RL diagOB_StatScal1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_StatScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_StatVelU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_StatVelV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /diagOB_STATE_3D/
     &    diagOB_StatScal1, diagOB_StatScal2,
     &    diagOB_StatVelU,  diagOB_StatVelV
#endif /* DIAGNOSTICS_OB_3D_STATE */
#ifdef DIAGNOSTICS_OB_2D_STATE
C     DIAGOB 2-dim. fields
      _RL diagOB_Surf1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL diagOB_Surf2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /diagOB_STATE_2D/
     &    diagOB_Surf1, diagOB_Surf2
#endif /* DIAGNOSTICS_OB_2D_STATE */

#ifdef DIAGNOSTICS_OB_TENDENCY
      _RL diagOB_TendScal1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_TendScal2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_TendVelU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL diagOB_TendVelV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /diagOB_TENDENCY/
     &    diagOB_TendScal1, diagOB_TendScal2,
     &    diagOB_TendVelU,  diagOB_TendVelV
#endif /* DIAGNOSTICS_OB_TENDENCY */

#endif /* ALLOW_DIAGNOSTICS_OB */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
