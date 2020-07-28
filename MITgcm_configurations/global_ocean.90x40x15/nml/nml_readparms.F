#include "DIAG_OPTIONS.h"



C Example namelist routine

C     !INTERFACE:
      SUBROUTINE NML_READPARMS( myThid)

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "NML.h"

C     !LOCAL VARIABLES:
      INTEGER j, k, ku, myThid
      INTEGER dimx, dimy, dimz
      CHARACTER*10 field_names(4)
      CHARACTER*80 msgBuf

      NAMELIST / NML_PARMS /
     & dimx, dimy, dimz, field_names

C     Initializing variables
      dimx = 0
      dimy = 0
      dimz = 0
      DO j=1, 4
        field_names(j) = ' '
        nml_fields(j) = ' '
      ENDDO

C     Initialzing array in NML.h common block (Note: values not set with data.nml)
      DO j=1, 10
        DO k=1, 10
          nml_map(j,k) = 0.0
        ENDDO
      ENDDO

C     Printing before setting namelist values to data.nml values:
      print *,"NML_READPARMS: Printing NAMELIST variables"
      print *,"NML_READPARMS:dimx before reading data.nml:", dimx
      print *,"NML_READPARMS:dimy before reading data.nml:", dimy
      print *,"NML_READPARMS:dimz before reading data.nml:", dimz
      print *,"NML_READPARMS:field_names before reading data.nml:",
     &         field_names

C     Reading from data.namelist
      WRITE(msgBuf,'(2A)')
     &     ' DIAGNOSTICS_READPARMS: opening data.nml'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,SQUEEZE_RIGHT,1)

      CALL OPEN_COPY_DATA_FILE('data.nml',
     &     'NML_READPARMS', ku, myThid )

      WRITE(msgBuf,'(2A)') 'S/R NML_READPARMS,',
     &     ' read namelist "NML_PARMS": start'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      READ  (ku,NML=NML_PARMS)
      WRITE(msgBuf,'(2A)') 'S/R NML_READPARMS,',
     &     ' read namelist "NML_PARMS": OK'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

C     Printing after setting namelist values to data.nml values:
      print *,"NML_READPARMS:dimx after reading data.nml:", dimx
      print *,"NML_READPARMS:dimy after reading data.nml:", dimy
      print *,"NML_READPARMS:dimz after reading data.nml:", dimz
      print *,"NML_READPARMS:field_names after reading data.nml:",
     &         field_names

#ifdef SINGLE_DISK_IO
      CLOSE(ku)
#else
      CLOSE(ku,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Printing before COMMON block variable values set to Namelist values:
      print *,"NML_READPARMS:Now printing NML.h COMMON block variables:"
      print *,"NML_READPARMS:dimension_x before namelist:",dimension_x
      print *,"NML_READPARMS:dimension_y before namelist:",dimension_y
      print *,"NML_READPARMS:dimension_z before namelist:",dimension_z
      print *,"NML_READPARMS:nml_fields before namelist:",nml_fields

C     Fill NML_PARMS Common Block with Namelist Info
      dimension_x = dimx
      dimension_y = dimy
      dimension_z = dimz
      DO j=1, 4
        nml_fields(j) = field_names(j)
      ENDDO

C     Printing after COMMON block variable values set to Namelist values:
      print *,"NML_READPARMS:dimension_x after namelist:",dimension_x
      print *,"NML_READPARMS:dimension_y after namelist:",dimension_y
      print *,"NML_READPARMS:dimension_z after namelist:",dimension_z
      print *,"NML_READPARMS:nml_fields after namelist:",nml_fields

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER
      RETURN
      END
