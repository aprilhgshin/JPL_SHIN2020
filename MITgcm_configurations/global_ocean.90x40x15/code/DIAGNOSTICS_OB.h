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
C     These common block variables are initialized in diagnostics_ob_init_varia.F

      CHARACTER*10 ob_fldNames(nOB_fld)
C      CHARACTER*98 ob_filenames(nOB_mask)
      CHARACTER*16 ob_filenames(nOB_mask)

      LOGICAL fld_choice(nFldOpt)


      _RL ob_subMask(nOB_mask,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
C      _RL subField_avg(nOB_fld,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
C      _RL subField(nOB_fld,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)

C     First row contains local i's. Second row contains local j's.
      INTEGER sub_local_ij_ob(nOB_mask, 2, sNx + sNy)

      INTEGER sub_glo_indices_allproc(nOB_mask, nPx*nPy, sNx + sNy)
      INTEGER numOBPnts_allproc(nOB_mask, nPx*nPy)

      _RL subFieldOnMask(nOB_mask,nOB_fld, sNx + sNy)
      _RL subFieldOnMask_avg(nOB_mask,nOB_fld, sNx + sNy)

      INTEGER lookup_table(nOB_mask, Ny*Nx)
      _RL global_ob((sNy+sNx)*(nPx*nPy))

      _RL global_ob_mask(nOB_mask,Nx, Ny,nSx,nSy)

      INTEGER num_ob_points(nOB_mask)

      _RL avgPeriod_ob
      _RL deltaT_ob
      _RL startTime_ob
      _RL endTime_ob
      _RL nTimeSteps_ob
      _RL time_passed
C     ==========================================================================

C      COMMON / DIAG_OB_EXTRACT_R /
C     &     ob_subMask, subField_avg, subField,
C     &     global_ob, subFieldOnMask_avg, global_ob_mask
C     &     subFieldOnMask, nTimeSteps_ob, time_passed
C      COMMON / DIAG_OB_EXTRACT_I /
C     &     lookup_table, num_ob_points
C     &     avgPeriod_ob, deltaT_ob
C      COMMON / DIAG_OB_EXTRACT_C /
C     &     ob_fldNames, ob_filenames
C      COMMON / DIAG_OB_EXTRACT_L /
C     &     fld_choice

      COMMON / DIAG_OB_EXTRACT_R /
     &     ob_subMask, global_ob,
     &     global_ob_mask, subFieldOnMask_avg,
     &     subFieldOnMask, nTimeSteps_ob, time_passed,
     &     startTime_ob, endTime_ob, avgPeriod_ob, deltaT_ob
      COMMON / DIAG_OB_EXTRACT_I /
     &     lookup_table, sub_local_ij_ob, sub_glo_indices_allproc,
     &     numOBPnts_allproc, num_ob_points
      COMMON / DIAG_OB_EXTRACT_C /
     &     ob_fldNames, ob_filenames
      COMMON / DIAG_OB_EXTRACT_L /
     &     fld_choice

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
