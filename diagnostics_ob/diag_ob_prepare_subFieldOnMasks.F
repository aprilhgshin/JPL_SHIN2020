#include "DIAG_OPTIONS.h"


      SUBROUTINE SET_TO_INI_STATE(myThid)
C     Initializes field value average  to zero
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER nFlds2D, nFlds3D
        INTEGER i, k, ob_id, fld2D, fld3D

        DO ob_id=1, nOB_mask
C        Total number of 2D fields
         nFlds2D = ob_nFlds2D(ob_id)
         DO fld2D=1, nFlds2D
          DO i=1, sNx+sNy
           subFieldOnMask_2Davg(ob_id,fld2D,i) = 0.0
          ENDDO
         ENDDO

c        Total number of 3D fields
         nFlds3D = ob_nFlds3D(ob_id)
         DO fld3D=1, nFlds3D
          DO i=1, sNx+sNy
           DO k=1, ob_levels3D(fld3D, ob_id)
            subFieldOnMask_3Davg(ob_id,fld3D,i,k) = 0.0
           ENDDO
          ENDDO
         ENDDO
        ENDDO

      END


      SUBROUTINE SET_SUBFIELDS(myThid, myIter)

C     Stores desired fields chosen by user into array subFields
C     Assign unique tags to all 2D and 3D fields for MPI_SEND
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#endif
#include "DIAGNOSTICS_OB.h"


      INTEGER ob_id, myThid, loc_count, myIter
      INTEGER i, j, k, bi, bj, numValid, tag
      INTEGER nFlds2D, nFlds3D, fld2D, fld3D
      INTEGER, PARAMETER :: print_debug = 0

      IF (print_debug .eq. 1) then
       print *,"Inside SET_SUBFIELDS"
       print*,"SET_SUBFIELDS: ob_levels3D ",ob_levels3D
       DO ob_id=1, nOB_mask
         print*,"SET_SUBFIELDS: mask:",ob_id
         print *,"sub_local_ij_ob j" ,sub_local_ij_ob(ob_id,2,:)
         print *,"sub_local_ij_ob i" ,sub_local_ij_ob(ob_id,1,:)
       ENDDO
      ENDIF

C     Initialize tag number to 1. Tag number will be incrememted after each field is extracted.
      tag = 1
      DO ob_id=1, nOB_mask

C       Number of open boundary points for mask ob_id and process mpiMyId (added one because of fortran's indices starting from 1 )
        numValid = numOBPnts_allproc(ob_id,mpiMyId+1)

C      Looping through 2D fields first, if any.
       nFlds2D = ob_nFlds2D(ob_id)
       DO fld2D=1, nFlds2D

C       Save unique tag for field and incremement tag for next field to be extracted
        ob_tags(ob_id, 1, fld2D) =  tag
        tag = tag + 1

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO loc_count=1, numValid
C          i,j indices with respect to process' sub mask for open boundary points (Note: all processes have their own declared sub_local_ij_ob)
           j = sub_local_ij_ob(ob_id,2,loc_count)
           i = sub_local_ij_ob(ob_id,1,loc_count)

C          Comparing strings provided by user to fields in diagnostics_ob
           IF (ob_flds2D(fld2D,ob_id) .eq. 'AREA') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = AREA(i,j,bi,bj)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'HEFF') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = HEFF(i,j,bi,bj)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'HSNOW') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = HSNOW(i,j,bi,bj)
#ifdef SEAICE_VARIABLE_SALINITY
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'HSALT') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = HSALT(i,j,bi,bj)
#endif /* SEAICE_VARIABLE_SALINITY */
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'UICE') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = UICE(i,j,bi,bj)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'VICE') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = VICE(i,j,bi,bj)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'ETAN') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = etaN(i,j,bi,bj)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'ETAH') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = etaH(i,j,bi,bj)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
C      !end of loop for nFlds2D
       ENDDO

C      Looping through 3D fields, if any.
       nFlds3D = ob_nFlds3D(ob_id)
       DO fld3D=1, nFlds3D

C       Save unique tag for field and incremement tag for next field to be extracted
        ob_tags(ob_id, 2, fld3D) = tag
        tag = tag + 1

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO loc_count=1, numValid

C          i,j indices with respect to process' sub mask for open boundary points
           j = sub_local_ij_ob(ob_id,2,loc_count)
           i = sub_local_ij_ob(ob_id,1,loc_count)
           DO k=1, ob_levels3D(fld3D,ob_id)

C          Comparing strings provided by user to fields in diagnostics_ob
            IF (ob_flds3D(fld3D,ob_id) .eq. 'UVEL') then
             subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        uVel(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'VVEL') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        vVel(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'WVEL') then
             subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        wVel(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'THETA') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        theta(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'SALT') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k)=
     &                                       salt(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'GU') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        gU(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'GV') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        gV(i,j,k,bi,bj)
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO

        IF (print_debug .eq. 1) then
          print*,"DOUBLE CHECKING:"
          print*,"ob_flds3D(fld3D,ob_id)",ob_flds3D(fld3D,ob_id)
          print*,"fld3D", fld3D
          print*,"ob_flds3D: ", ob_flds3D
          print*,"theta at time level", time_level
          print*,"numvalid:",numValid
          print*,"mask",ob_id
          print*,"i", sub_local_ij_ob(ob_id,1,:)
          print*,"j",  sub_local_ij_ob(ob_id,2,:)
        ENDIF

C       !End of loop for nFlds3D
       ENDDO
C     !End of loop for nOB_mask
      ENDDO

      IF (print_debug .eq. 1) then
       print *, "FINAL set 2D:"
       call PRINT_FLOAT_ARR(subFieldOnMask_2D(1,1,:), sNx+sNy)
       print *,"FINAL set 3D first level:"
       call PRINT_FLOAT_ARR(subFieldOnMask_3D(1,1,:,1), (sNx+sNy))
      ENDIF

      END


      SUBROUTINE CUMULATE_FLD_TIME(myThid)
C     Accumulate extracted field values on OB points
C       - Add extracted field values on OB points to all extracted field values on OB points
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

      INTEGER myThid
      INTEGER i, k, ob_id, fld, fld2D, fld3D

      DO ob_id=1, nOB_mask
       DO i=1, sNx+sNy
        DO fld2D=1, ob_nFlds2D(ob_id)
         subFieldOnMask_2Davg(ob_id,fld2D,i)
     &                      = subFieldOnMask_2Davg(ob_id,fld2D,i)
     &                        + subFieldOnMask_2D(ob_id,fld2D,i)
        ENDDO
        DO fld3D=1, ob_nFlds3D(ob_id)
         DO k=1, ob_levels3D(fld3D, ob_id)
          subFieldOnMask_3Davg(ob_id,fld3D,i,k)
     &                    = subFieldOnMask_3Davg(ob_id,fld3D,i,k)
     &                      + subFieldOnMask_3D(ob_id,fld3D,i,k)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      END


      SUBROUTINE TIMEAVE_OB_FLD(myThid, myTime)
C       Computes average of Diagnostics_ob fields using the number of time steps taken in current period
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER i, k,fld, ob_id, fld2D, fld3D
        _RL myTime

        DO ob_id=1, nOB_mask
          DO i=1, sNx+sNy
           DO fld2D=1, ob_nFlds2D(ob_id)
            subFieldOnMask_2Davg(ob_id,fld2D,i) =
     &              subFieldOnMask_2Davg(ob_id,fld2D,i)/nTimeSteps_ob
           ENDDO
           DO fld3D=1, ob_nFlds3D(ob_id)
            DO k=1, ob_levels3D(fld3D,ob_id)
             subFieldOnMask_3Davg(ob_id,fld3D,i,k) =
     &              subFieldOnMask_3Davg(ob_id,fld3D,i,k)/nTimeSteps_ob
            ENDDO
           ENDDO
          ENDDO

       ENDDO
C       Set number of time steps within averaging period back to zero for the next averaging period
       nTimeSteps_ob = 0
      END
