#ifdef ALLOW_DIAGNOSTICS_OB
#include "USER_INPUT.h"


C     ================== Global Variables for open boundary ====================

C     Defined by or derived from user input in data.diagnostics_ob:
C     ob_flds2D    :: Char array of names of 2D fields for each open boundary mask.
C     ob_flds3D    :: Char array of names of 3D fields for each open boundary mask.
C     ob_fnames    :: Char array of filenames for open boundary mask files.
C     ob_levels3D  :: Int array of depths of 3D fields for each open boundary mask.
C     ob_nFlds2D   :: Int array of number of 2D fields in each open boundary mask
C     ob_nFlds3D   :: Int array of number of 3D fields in each open boundary mask
C     ob_tags      :: Int array of unique tags assigned to every field in every mask
C     combineMaskTimeLevels :: Logical for whether all time levels should be written
C                              to a single bin file or to separate bin files
C     ob_filePrec :: file precision for binary file output
C     avgPeriod_ob :: averaging period
C     startTime_ob :: start time for writing output
C     endTime_ob :: end time for writing output

C     Everything else:

C     global_ob2D :: _RL array for final output of combined time-averaged 2D field values on OB points
C     global_ob2D :: _RL array for final output of combined time-averaged 3D field values on OB points
C     nTimeSteps_ob :: Integer value for number of time steps taken within averaging period
C     time_passed :: total time passed
C     time_level :: number of time levels passed i.e. number of averaging period passed

C     Containing information for all open boundary masks:
C     ob_subMask :: _RL array for portion of open boundary global mask assigned to process
C     sub_local_ij_ob :: _RL array for i,j indices of open boundary points wrt ob_subMask domain
C     sub_glo_indices_allproc :: _RL array for global indices of open boundary points
C                                for each process
C     numOBPnts_allproc :: _RL array for number of open boundary points in each process
C     subFieldOnMask_2D :: _RL array for 2D field values on open boundary points
C     subFieldOnMask_3D :: _RL array for 3D field values on open boundary points
C     subFieldOnMask_2Davg :: _RL array for time-averaged 2D field values on open boundary points
C     subFieldOnMask_3Davg :: _RL array for time-averaged 3D field values on open boundary points
C     lookup_table :: _RL array containing all open boundary global masks
C     global_ob_mask :: _RL array containing all open boundary masks
C-------------------------------------------------------------------------------

      CHARACTER*8 ob_flds2D(MAX_NFLDS, nOB_mask)
      CHARACTER*8 ob_flds3D(MAX_NFLDS, nOB_mask)
      CHARACTER*30 ob_fnames(nOB_mask)
      INTEGER ob_levels3D(MAX_NFLDS, nOB_mask)
      INTEGER ob_nFlds2D(nOB_mask)
      INTEGER ob_nFlds3D(nOB_mask)
      INTEGER ob_tags(nOB_mask, 2, MAX_NFLDS)
      LOGICAL combineMaskTimeLevels
      INTEGER ob_filePrec
      _RL avgPeriod_ob
      _RL startTime_ob
      _RL endTime_ob


      REAL*8 global_ob2D((sNy+sNx)*(nPx*nPy))
      REAL*8 global_ob3D((sNy+sNx)*(nPx*nPy), Nr)
      _RL nTimeSteps_ob
      _RL time_passed
      INTEGER time_level


      _RL ob_subMask(nOB_mask,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)

C     First row contains local i's. Second row contains local j's.
      INTEGER sub_local_ij_ob(nOB_mask, 2, sNx + sNy)

      INTEGER sub_glo_indices_allproc(nOB_mask, nPx*nPy, sNx + sNy)
      INTEGER numOBPnts_allproc(nOB_mask, nPx*nPy)

      _RL subFieldOnMask_2D(nOB_mask,MAX_NFLDS, sNx + sNy)
      _RL subFieldOnMask_3D(nOB_mask,MAX_NFLDS, sNx + sNy, Nr)
      _RL subFieldOnMask_2Davg(nOB_mask,MAX_NFLDS, sNx + sNy)
      _RL subFieldOnMask_3Davg(nOB_mask,MAX_NFLDS, sNx + sNy, Nr)

      INTEGER lookup_table(nOB_mask, Ny*Nx)
      _RL global_ob_mask(nOB_mask,Nx, Ny,nSx,nSy)
C     ==========================================================================

      COMMON / DIAG_OB_EXTRACT_R /
     &     ob_subMask,
     &     global_ob_mask,
     &     nTimeSteps_ob, time_passed,
     &     startTime_ob, endTime_ob, avgPeriod_ob,
     &     global_ob2D, global_ob3D,
     &     subFieldOnMask_2D, subFieldOnMask_2Davg,
     &     subFieldOnMask_3D, subFieldOnMask_3Davg

      COMMON / DIAG_OB_EXTRACT_I /
     &     lookup_table, sub_local_ij_ob, sub_glo_indices_allproc,
     &     numOBPnts_allproc, ob_filePrec,
     &     ob_levels3D, ob_nFlds2D, ob_nFlds3D, ob_tags, time_level
      COMMON / DIAG_OB_EXTRACT_C /
     &     ob_flds2D, ob_flds3D, ob_fnames
      COMMON / DIAG_OB_EXTRACT_C /
     &     combineMaskTimeLevels

#endif /* ALLOW_DIAGNOSTICS_OB */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
