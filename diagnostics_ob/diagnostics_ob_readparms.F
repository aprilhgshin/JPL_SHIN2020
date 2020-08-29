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

C ----------------------------
C     Open Boundary Extraction Namelist:
C
C     nml_fields2D :: 2D field names for ocean state and flux variables
C     nml_fields3D :: 3D field names for ocean state and flux variables
C     nml_fileNames :: OB mask filenames
C     nml_levels3D :: depths starting from surface for each 3D field chosen by user
C     nml_combMaskTimeLevels :: Logical for either writing all time levels to binary file
C                               or write in separate binary files
C     avgPeriod_nml :: averaging period (frequency) of output
C     startTime_nml :: start time for output
C     endTime_nml :: end time for output
C
C -----------------------------

      CHARACTER*8 nml_fields2D(MAX_NFLDS, nOB_mask)
      CHARACTER*8 nml_fields3D(MAX_NFLDS, nOB_mask)

      CHARACTER*30 nml_fileNames(nOB_mask)
      INTEGER nml_levels3D(MAX_NFLDS, nOB_mask)
      LOGICAL nml_combMaskTimeLevels

      INTEGER nml_filePrec
      _RL avgPeriod_nml
      _RL startTime_nml
      _RL endTime_nml
      _RL undefRL

      NAMELIST / DIAG_OB_EXTRACT /
     &     avgPeriod_nml, startTime_nml, endTime_nml,
     &     nml_filePrec,
     &     nml_fields2D, nml_fields3D,
     &     nml_fileNames, nml_levels3D,
     &     nml_combMaskTimeLevels

      IF ( .NOT.useDiagnostics_ob ) THEN
C-    pkg DIAGNOSTICS_OB is not used
C-    Track pkg activation status:
C     print a (weak) warning if data.diagnostics_ob is found
         CALL PACKAGES_UNUSED_MSG( 'useDiagnostics_ob', ' ', ' ' )
        RETURN
      ENDIF


      IF (.NOT. useSEAICE ) then
C     if pkg SEAICE is not used, print warning
        CALL PACKAGES_UNUSED_MSG( 'useSEAICE', ' ', ' ' )
C        CALL ALL_PROC_DIE( myThid )
        STOP 'ABNORMAL END: S/R DIAG_OB: SEAICE PACKAGE NOT ACTIVATED'
      ENDIF

      _BEGIN_MASTER(myThid)


C-------------------------------------------------------------------------------
C       Initializing local variables and variables from DIAG_OB_EXTRACT COMMON block:
C-------------------------------------------------------------------------------

      startTime_nml = endTime
      endTime_nml = startTime
      nml_filePrec = 64
      nml_combMaskTimeLevels = .TRUE.


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

C------------------------------------------------------------------------------
C     Fill DIAGNOSTICS_OB.h Common Block with Namelist Info
C------------------------------------------------------------------------------

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

      avgPeriod_ob = avgPeriod_nml
      startTime_ob = startTime_nml
      endTime_ob = endTime_nml

C     If startTime, endTime, and avgPeriod are not specified by user then:
C      IF (avgPeriod_nml .eq. undefRL) then
C        avgPeriod_ob = deltaTClock
C      ELSE
C        avgPeriod_ob = avgPeriod_nml
C      ENDIF

C      IF (startTime_nml .eq. undefRL) then
C        startTime_ob = startTime
C      ELSE
C        startTime_ob = startTime_nml
C      ENDIF

C      IF (endTime_nml .eq. undefRL) then
C        endTime_ob = endTime
C      ELSE
C        endTime_ob = endTime_nml
C      ENDIF

      ob_filePrec = nml_filePrec
      combineMaskTimeLevels = nml_combMaskTimeLevels




      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
