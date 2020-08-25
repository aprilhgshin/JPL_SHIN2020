#ifdef ALLOW_DIAGNOSTICS_OB
#include "USER_INPUT.h"


C     ================== Global Variables for open boundary ====================

C     Defined by user input in data.diagnostics_ob:
C     ob_flds2D    :: Char array of names of 2D fields for each open boundary mask.
C     ob_flds3D    :: Char array of names of 3D fields for each open boundary mask.
C     ob_fnames    :: Char array of filenames for open boundary mask files.
C     ob_levels3D  :: Int array of depths of 3D fields for each open boundary mask.
C     ob_nFlds2D   :: Int array of number of 2D fields in each open boundary mask
C     ob_nFlds3D   :: Int array of number of 3D fields in each open boundary mask
C     ob_tags      :: Int array of unique tags assigned to every field in every mask

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

      INTEGER time_level
      LOGICAL combineMaskTimeLevels
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
     &     ob_levels3D, ob_nFlds2D, ob_nFlds3D, ob_tags, time_level
      COMMON / DIAG_OB_EXTRACT_C /
     &     ob_flds2D, ob_flds3D, ob_fnames
      COMMON / DIAG_OB_EXTRACT_C /
     &     combineMaskTimeLevels

#endif /* ALLOW_DIAGNOSTICS_OB */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
