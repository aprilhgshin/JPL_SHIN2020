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
      CHARACTER*8 nml_fldNames(2,10)
      INTEGER nml_maskTypes(2,10)
      INTEGER nml_3Dfld_depths(n3D_fld)

      _RL avgPeriod_nml
      _RL startTime_nml
      _RL endTime_nml

      NAMELIST / DIAG_OB_EXTRACT /
     &     avgPeriod_nml, startTime_nml, endTime_nml,
     &     nml_fldNames, nml_maskTypes, nml_3Dfld_depths

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

      DO i=1, 2
        DO j=1, 10
          nml_fldNames(i,j) = ' '
          nml_maskTypes(i,j) = 0
          fld_choice(i,j) = .FALSE.
        ENDDO
      ENDDO

      DO i=1, nOB_fld
        ob_fldNames(i) = ' '
      ENDDO

      DO i=1, n3D_fld
        nml_3Dfld_depths(i) = 0
      ENDDO

      DO fld=1, nOB_fld
        fld_to_output(fld) = 0
        fld_maskType(fld) = 0
        fld_nDim(fld) = 0
        fld_depth(fld) = 0
      ENDDO

C     Initializing averaging period (avgPeriod_ob) and time counter (deltaT_ob)
      avgPeriod_nml = 0.0
      startTime_nml = 0.0
      endTime_nml = 0.0
C     CHANGE THESE VALUES TO THE UNDEFINED CONSTANT


      time_passed = 0


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
      fld_choice(1,1) = area_ob
      fld_choice(1,2) = heff_ob
      fld_choice(1,3) = hsnow_ob
      fld_choice(1,4) = hsalt_ob
      fld_choice(1,5) = uice_ob
      fld_choice(1,6) = vice_ob
      fld_choice(1,7) = etaN_ob
      fld_choice(1,8) = etaH_ob

      fld_choice(1,9) = obnw_ob
      fld_choice(2,1) = uVel_ob
      fld_choice(2,2) = vVel_ob
      fld_choice(2,3) = theta_ob
      fld_choice(2,4) = salt_ob
      fld_choice(2,5) = gU_ob
      fld_choice(2,6) = gV_ob

C      INTEGER fld_to_output(nOB_fld)
C      INTEGER fld_maskType(nOB_fld)
C      INTEGER fld_nDim(nOB_fld)
C      INTEGER fld_depth(nOB_fld)

C     Fill DIAGNOSTICS_OB.h Common Block with Namelist Info
      k = 1
      counter3D = 1

      DO i=1, nFldOpt
        if (i .le. n2D_fld) then
          if ((fld_choice(1,i) .eqv. .TRUE.) .AND.
     &                (k .le. nOB_fld)) then
            ob_fldNames(k) = nml_fldNames(1,i)
            fld_to_output(k) = i
            fld_maskType(k) = nml_maskTypes(1,i)
            fld_nDim(k) = 2
            fld_depth(k) = 1
            k = k + 1
          endif
        else
          if ((fld_choice(2,i-n2D_fld) .eqv. .TRUE.) .AND.
     &                (k .le. nOB_fld)) then
            ob_fldNames(k) = nml_fldNames(2,i-n2D_fld)
            fld_to_output(k) = i
            fld_nDim(k) = 3
            fld_depth(k) = nml_3Dfld_depths(counter3D)
            k = k + 1
            counter3D = counter3D + 1
          endif
        endif
      ENDDO

      print *,"OB readparms: fld_to_output",fld_to_output
      print *,"OB readparms: fld_nDim",fld_nDim
      print *,"OB readparms: fld_depth",fld_depth

      print *,"OB readparms: calling CREATE_MASK_DICT"
      call CREATE_MASK_DICT()


      avgPeriod_ob = avgPeriod_nml
      startTime_ob = startTime_nml
      endTime_ob = endTime_nml



      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_DIAGNOSTICS_OB */

      RETURN
      END


      SUBROUTINE CREATE_MASK_DICT()
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER i, j, mask_id
        INTEGER mask_count(nMaskTypes)
        INTEGER mask_flds(nMaskTypes,nOB_fld)
C       Initializing arrays
        DO i=1, nMaskTypes
          mask_count(i) = 0
          DO j=1, nOB_fld
            mask_flds(i,j) = 0
          ENDDO
        ENDDO

        DO i=1, nMask_perOB
          maskId_used(i) = 0
        ENDDO

        print *,"CREATE_MASK_DICT: Num of fields to output:", nOB_fld
        print *,"CREATE_MASK_DICT: Mask IDs for each field:",
     &           fld_maskType

        DO i=1, nOB_fld
          if (fld_maskType(i) .eq. 1) then
            mask_count(1) = mask_count(1) + 1
            mask_flds(1, mask_count(1)) = i
          else if (fld_maskType(i) .eq. 2) then
            mask_count(2) = mask_count(2) + 1
            mask_flds(2, mask_count(2)) = i
          else if (fld_maskType(i) .eq. 3) then
            mask_count(3) = mask_count(3) + 1
            mask_flds(3, mask_count(3)) = i
          end if
        ENDDO

        print *,"CREATE_MASK_DICT: fld count for each mask",mask_count
        print *,"CREATE_MASK_DICT: fields for mask 1:",mask_flds(1,:)
        print *,"CREATE_MASK_DICT: fields for mask 2:",mask_flds(2,:)
        print *,"CREATE_MASK_DICT: fields for mask 3:",mask_flds(3,:)


        mask_id = 1
        DO i=1, nMaskTypes
          if (mask_count(i) .gt. 0) then
             maskId_used(mask_id) = i
             mask_id = mask_id + 1
             DO j=1, mask_count(i)
               mask_fields(mask_id, j) = mask_flds(i,j)
             ENDDO
             mask_id = mask_id + 1
           endif
        ENDDO

        DO i=1, nMask_perOB
          print *,"CREATE_MASK_DICT: fields for mask",i,":"
          print *, mask_fields(i,:)
        ENDDO

        print *,"CREATE_MASK_DICT: maskId_used",maskId_used

      END
