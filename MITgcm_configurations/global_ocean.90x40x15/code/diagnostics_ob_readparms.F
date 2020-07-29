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
      INTEGER bi,bj,i,j,k,l

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
C     nml_fldFiles : filenames for field files used in OB_EXTRACT_LOOKUP
C     avgPeriod_nml : averaging period
C     deltaT_nml : size of time step
C     totPhase_nml : total time phase
C -----------------------------
      CHARACTER*10 nml_fldNames(3)
      CHARACTER*40 nml_fldFiles(3)
      INTEGER avgPeriod_nml
      INTEGER deltaT_nml
      INTEGER totPhase_nml
      LOGICAL fld_choice(3)

C      NAMELIST / DIAG_OB_EXTRACT /
C     &       ob_fldNames, ob_fnames,
C     &       ob_subMask, lookup_table,
C     &       global_ob, sub_global_indices,
C     &       temp_subMask, bath_subMask,
C     &       subBathOnMask, subTempOnMask,
C     &       avgPeriod_ob, deltaT_ob, totPhase_ob

      NAMELIST / DIAG_OB_EXTRACT /
     &     nml_fldFiles, nml_fldNames,
     &     avgPeriod_nml, deltaT_nml, totPhase_nml

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


C     Initializing local array:
      DO j=1,3
        fld_choice(j) = .FALSE.
      ENDDO

C     Initializing open boundary sub mask
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
            DO k=1-OLx,sNx+OLx
              ob_subMask(k,j,bi,bj)   = 0.0
              DO i=1, nOB_fld
                subField_avg(i,k,j,bi,bj) = 0.0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO l = 1, sNy+sNx+1
        sub_global_indices(l) = 0
        DO k=1, nOB_fld
          subFieldOnMask(k,l) = 0.0
        ENDDO
      ENDDO

C     k limited to total number of open boundaries
      DO j = 1, nOB_mask
        DO l = 1, Ny*Nx
        lookup_table(j,l) = 0.0
        ENDDO
      ENDDO

      DO j=1, (sNy+sNx)*(nPx*nPy)
        global_ob(j) = 0.0
      ENDDO

C     Initializing to empty string
      DO k=1, nOB_fld
        ob_fldNames(k) = ' '
        ob_fldFiles(k) = ' '
      ENDDO

      DO k=1, nOB_mask
        ob_filenames(k) = ' '
      ENDDO

      DO k=1, 3
        nml_fldNames(k) = ' '
      ENDDO
      DO k=1, 3
        nml_fldFiles(k) = ' '
      ENDDO

C     Initializing averaging period (avgPeriod_ob) and time counter (deltaT_ob)
      avgPeriod_nml = 86400
      deltaT_nml = 86400
      totPhase_nml = 86400



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

C     Organize user input booleans into an array:
      fld_choice(1) = ob_hydrogSalt
      fld_choice(2) = ob_hydrogTheta
      fld_choice(3) = ob_bathymetry

C     Fill DIAGNOSTICS_OB.h Common Block with Namelist Info
      k = 1
      DO j=1, 3
        if ((fld_choice(j) .eqv. .TRUE.) .AND.
     &                (k .le. nOB_fld)) then
          ob_fldFiles(k) = nml_fldFiles(j)
          ob_fldNames(k) = nml_fldNames(j)
          k = k + 1
        end if
      ENDDO

      avgPeriod_ob = avgPeriod_nml
      deltaT_ob = deltaT_nml
      totPhase_ob = totPhase_nml

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END
