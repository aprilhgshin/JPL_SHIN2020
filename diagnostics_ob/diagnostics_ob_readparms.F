#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_READPARMS

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_READPARMS( myThid )

C     !DESCRIPTION:
C     Initialize DIAGNOSTICS_OB variables and constants.

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

C     !LOCAL VARIABLES:
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
      INTEGER bi,bj,i,j,k,l, ob, fld, counter3D

      NAMELIST /DIAGNOSTICS_OB_PARM01/
     &       diagOB_MNC,
     &       diagOB_StaV_Cgrid, diagOB_Tend_Cgrid,
     &       diagOB_applyTendT, diagOB_applyTendS,
     &       diagOB_applyTendU, diagOB_applyTendV,
     &       diagOB_doSwitch1,  diagOB_doSwitch2,
     &       diagOB_index1,     diagOB_index2,
     &       diagOB_param1,     diagOB_param2,
     &       diagOB_string1,    diagOB_string2,
     &       diagOB_Scal1File,  diagOB_Scal2File,
     &       diagOB_VelUFile,   diagOB_VelVFile,
     &       diagOB_Surf1File,  diagOB_Surf2File

C ----------------------------
C     Open Boundary Extraction Namelist:
C     nml_fldNames : field names for ocean state and flux variables

C -----------------------------

      CHARACTER*8 nml_fields2D(MAX_NFLDS, nOB_mask)
      CHARACTER*8 nml_fields3D(MAX_NFLDS, nOB_mask)

      CHARACTER*30 nml_fileNames(nOB_mask)
      INTEGER nml_levels3D(MAX_NFLDS, nOB_mask)

      INTEGER nml_filePrec
      _RL avgPeriod_nml
      _RL startTime_nml
      _RL endTime_nml
      _RL undefRL

      NAMELIST / DIAG_OB_EXTRACT /
     &     avgPeriod_nml, startTime_nml, endTime_nml,
     &     nml_filePrec,
     &     nml_fields2D, nml_fields3D,
     &     nml_fileNames, nml_levels3D

      IF ( .NOT.useDIAGNOSTICS_OB ) THEN
C-    pkg DIAGNOSTICS_OB is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.DIAGNOSTICS_OB is found
         CALL PACKAGES_UNUSED_MSG( 'useDIAGNOSTICS_OB', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values for DIAGNOSTICS_OB
      diagOB_MNC         = useMNC
      diagOB_StaV_Cgrid  = .TRUE.
      diagOB_Tend_Cgrid  = .TRUE.
      diagOB_applyTendT  = .FALSE.
      diagOB_applyTendS  = .FALSE.
      diagOB_applyTendU  = .FALSE.
      diagOB_applyTendV  = .FALSE.
C-    additional parameters:
      diagOB_doSwitch1   = .FALSE.
      diagOB_doSwitch2   = .FALSE.
      diagOB_index1      = 0
      diagOB_index2      = 0
      diagOB_param1      = 0. _d 0
      diagOB_param2      = 0. _d 0
      diagOB_string1     = ' '
      diagOB_string2     = ' '
C-    file names for initial conditions:
      diagOB_Scal1File   = ' '
      diagOB_Scal2File   = ' '
      diagOB_VelUFile    = ' '
      diagOB_VelVFile    = ' '
      diagOB_Surf1File   = ' '
      diagOB_Surf2File   = ' '

C-------------------------------------------------------------------------------
C       Initializing local variables and variables from DIAG_OB_EXTRACT COMMON block:
C-------------------------------------------------------------------------------

      avgPeriod_nml = undefRL
      startTime_nml = undefRL
      endTime_nml = undefRL
      time_passed = 0.0
      nml_filePrec = 64

      DO j=1, nOB_mask
        DO i=1, MAX_NFLDS
          nml_fields2D(i,j) = ' '
          nml_fields3D(i,j) = ' '
          nml_levels3D(i,j) = 0
        ENDDO
        nml_fileNames(j) = ' '
      ENDDO

      DO j=1, nOB_mask
        DO i=1, MAX_NFLDS
          ob_flds2D(i,j) = ' '
          ob_flds3D(i,j) = ' '
          ob_levels3D(i,j) = 0
        ENDDO
        ob_fnames(j) = ' '
        ob_nFlds2D(j) = 0
        ob_nFlds3D(j) = 0
      ENDDO


C End of Initializing-----------------------------------------------------------

      WRITE(msgBuf,'(A)')
     &      'DIAGNOSTICS_OB_READPARMS: opening data.diagnostics_ob'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I             'data.diagnostics_ob', 'DIAGNOSTICS_OB_READPARMS',
     O              iUnit,
     I              myThid )

C     Read parameters from open data file
      WRITE(msgBuf,'(2A)') 'S/R DIAGNOSTICS_READPARMS,',
     &     ' read namelist "DIAG_OB_EXTRACT": start'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      READ  (UNIT=iUnit,NML=DIAG_OB_EXTRACT)
      WRITE(msgBuf,'(2A)') 'S/R DIAGNOSTICS_READPARMS,',
     &     ' read namelist "DIAG_OB_EXTRACT": OK'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Make sure that we locally honor the global MNC on/off flag
      diagOB_MNC = diagOB_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      diagOB_MNC = .FALSE.
#endif
      diagOB_MDSIO = (.NOT. diagOB_MNC) .OR. outputTypesInclusive


C     Fill DIAGNOSTICS_OB.h Common Block with Namelist Info

      DO ob=1, nOB_mask
        ob_fnames(ob) = nml_fileNames(ob)

        DO i=1, MAX_NFLDS
          ob_flds2D(i,ob) = nml_fields2D(i,ob)
          if (nml_fields2D(i,ob) .ne. ' ') then
            ob_nFlds2D(ob) = ob_nFlds2D(ob) + 1
          endif
          ob_flds3D(i,ob) = nml_fields3D(i,ob)
          if (nml_fields3D(i,ob) .ne. ' ') then
            ob_nFlds3D(ob) = ob_nFlds3D(ob) + 1
          endif
          ob_levels3D(i,ob) = nml_levels3D(i,ob)
        ENDDO
      ENDDO


      IF (avgPeriod_nml .eq. undefRL) then
        avgPeriod_ob = deltaTClock
      ELSE
        avgPeriod_ob = avgPeriod_nml
      ENDIF

      IF (startTime_nml .eq. undefRL) then
        startTime_ob = startTime
      ELSE
        startTime_ob = startTime_nml
      ENDIF

      IF (endTime_nml .eq. undefRL) then
        endTime_ob = endTime
      ELSE
        endTime_ob = endTime_nml
      ENDIF

      ob_filePrec = nml_filePrec




      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
