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


C     MAX_NMASKS :: Int value for assumed maximum number of open boundary masks.
C     MAX_NFLDS  :: Int value for assumed maximum number of fields per open boundary mask.
C     ob_allFlds :: Char array of names of all fields.

C     Defined by user input in data.diagnostics_ob:
C     ob_flds2D    :: Char array of names of 2D fields for each open boundary mask.
C     ob_flds3D    :: Char array of names of 3D fields for each open boundary mask.
C     ob_fnames    :: Char array of filenames for open boundary mask files.
C     ob_levels3D  :: Int array of depths of 3D fields for each open boundary mask.
C     ob_nFlds2D   :: Int array of number of 2D fields in each open boundary mask
C     ob_nFlds3D   :: Int array of number of 3D fields in each open boundary mask
C     ob_tags      :: Int array of unique tags assigned to every field in every mask


      INTEGER, PARAMETER :: MAX_NMASKS = 12
      INTEGER, PARAMETER :: MAX_NFLDS = 20

      INTEGER ob_filePrec
      CHARACTER*8 ob_flds2D(MAX_NFLDS, nOB_mask)
      CHARACTER*8 ob_flds3D(MAX_NFLDS, nOB_mask)

      CHARACTER*30 ob_fnames(nOB_mask)
      INTEGER ob_levels3D(MAX_NFLDS, nOB_mask)
      INTEGER ob_nFlds2D(nOB_mask)
      INTEGER ob_nFlds3D(nOB_mask)

      INTEGER ob_tags(nOB_mask, 2, MAX_NFLDS)

      _RL ob_subMask(nOB_mask,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)

C     First row contains local i's. Second row contains local j's.
      INTEGER sub_local_ij_ob(nOB_mask, 2, sNx + sNy)

      INTEGER sub_glo_indices_allproc(nOB_mask, nPx*nPy, sNx + sNy)
      INTEGER numOBPnts_allproc(nOB_mask, nPx*nPy)

      _RL subFieldOnMask_2D(nOB_mask,MAX_NFLDS, sNx + sNy)
C      _RL subFieldOnMask_3D(nOB_mask,MAX_NFLDS, Nr, sNx + sNy)
      _RL subFieldOnMask_3D(nOB_mask,MAX_NFLDS, sNx + sNy, Nr)
      _RL subFieldOnMask_2Davg(nOB_mask,MAX_NFLDS, sNx + sNy)
C      _RL subFieldOnMask_3Davg(nOB_mask,MAX_NFLDS, Nr, sNx + sNy)
      _RL subFieldOnMask_3Davg(nOB_mask,MAX_NFLDS, sNx + sNy, Nr)

      INTEGER lookup_table(nOB_mask, Ny*Nx)
      REAL*8 global_ob2D((sNy+sNx)*(nPx*nPy))
      REAL*8 global_ob3D((sNy+sNx)*(nPx*nPy), Nr)


      _RL global_ob_mask(nOB_mask,Nx, Ny,nSx,nSy)

      INTEGER num_ob_points(nOB_mask)

      _RL avgPeriod_ob
      _RL deltaT_ob
      _RL startTime_ob
      _RL endTime_ob
      _RL nTimeSteps_ob
      _RL time_passed
C     ==========================================================================

      COMMON / DIAG_OB_EXTRACT_R /
     &     ob_subMask,
     &     global_ob_mask,
     &     nTimeSteps_ob, time_passed,
     &     startTime_ob, endTime_ob, avgPeriod_ob, deltaT_ob,
     &     global_ob2D, global_ob3D,
     &     subFieldOnMask_2D, subFieldOnMask_2Davg,
     &     subFieldOnMask_3D, subFieldOnMask_3Davg

      COMMON / DIAG_OB_EXTRACT_I /
     &     lookup_table, sub_local_ij_ob, sub_glo_indices_allproc,
     &     numOBPnts_allproc, num_ob_points, ob_filePrec,
     &     ob_levels3D, ob_nFlds2D, ob_nFlds3D, ob_tags
      COMMON / DIAG_OB_EXTRACT_C /
     &     ob_flds2D, ob_flds3D, ob_fnames


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
