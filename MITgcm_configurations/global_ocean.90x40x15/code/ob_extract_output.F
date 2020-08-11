C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678
#include "DIAG_OPTIONS.h"

C Changes to implement:
C question: does everyting get run one process at a time - init_vaira, do the model io

C DONE Make global_mask a coommon block variable
C DONE Initialize global_mask
C DONE Lookup_table
C DONE count total_valid for each open boundary - save these values in another array

C DONE MPI_send and recv only arrays with valid elements
C DONE Don't send global indices
C DONE So message tag will correspond to fld number

C make this use 3D fields
C y,z
C changing declaratoin in diagnostics_ob.h
C order of fields correspond to depth too. like depth 1 is eta. deptha 2 - 51 contains theta, etc.
C So subFldsOnMask (nOB_fld, sNx+sNy+1, nZ )

C fields to input :
C Everything i have but gu and gv
C Include HEFF, AREA, HSNON, HSALT, uice, vice from seaice.h
C OBNw from OBCS.h

C DONE dont use the word phase
C DONE change totPhase to having start time and end time instead

C Have user specify what depth it wants
C In common block of diagnostics_ob.h, create an array of depth dimensions of each field

C change k to i --- fix all do loops with this done like this

C change data file comment line 42, and change 55 to have 7 days + 10 hours

C change fld names to be lenght 8

C DONE Make do loop only go through 1 through sNx etc. not -OLx

C change data.diagnostics to output fields 2D and 3D
C Change frequency for both to be 1day  i.e. 86400 seconds

C apply program sturcutre on diagnostics output and compare with my outputs

C Very last addition: passive tracers





       SUBROUTINE OB_EXTRACT_OUTPUT( myTime, myIter, myThid )
C-------------------------------------------------------------------------------
C      OB_EXTRACT_OUTPUT:
C      Extracts 2D virtual transects, slices of ECCO 3D global model fields
C      in depth vs. distance. The tiles over which the virtual transects span
C      contain lateral open boundary conditions of ocean state variables and
C      fluxes for regional ocean models.
C      To reduce the runtime and memory usage, we divide each tile (open
C      boundary mask input and fields) into multiple sub-tiles of equal dimension,
C      use the Message Passage Interface (MPI) standard to assign each sub-mask
C      and sub-fields to a unique process. Each process will identify open
C      boundary points in their assigned sub-mask and sends an array of field
C      values and indices with respect to the global domain corresponding to
C      each of the open boundary points back to the parent process. Parent
C      process combines these arrays and writes result to binary files
C
C      Program outputs at the end of each averaginng period.
C      input parameters: myTime, myIter, myThid
C
C-------------------------------------------------------------------------------

#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#ifdef ALLOW_OBCS
C#include "OBCS_FIELDS.h"
#endif
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#endif
#include "DIAGNOSTICS_OB.h"

       _RL myTime
       integer myThid
       integer myIter

C     From PARAMS.h:
C     nIter0              :: Start time-step number of for this run
C     startTime/ endTime
C     deltaTClock

C     From EESUPPORT.h : all initialized by subroutine INI_PROCS:
C     mpiNprocs - No. of MPI processes.
C     mpiMyId   - MPI process id of me.
C     mpiComm   - MPI communicator to use.
C     mpiPx     - My MPI proc. grid X coord
C     mpiPy     - My MPI proc. grid Y coord
C     mpiXGlobalLo - My bottom-left (south-west) x-coordinate in
C                    global domain.
C     mpiYGlobalLo - My bottom-left (south-west) y-coordinate in
C                    global domain.
C     mpi_myXGlobalLo :: List of all processors bottom-left X-index in global domain
C     mpi_myYGlobalLo :: List of all processors bottom-left Y-index in global domain

C-------------------------------------------------------------------------------
C      !LOCAL VARIABLES:

C      global_mask      - 4 dimensional _RL array that combines all open boundarys
C      recv_arr         - buffer array to store incoming arrays from MPI_RECV
C      g_indices        - array that holds global indices of sub fields from each child process

C      msgBuf           - Used to write statements into the output file
C      global_ob_index  - Integer indicating bottom-left global element number
C      ob_index         - Integer obtained from lookup table in corresponding global index
C      new_goi          - Integer value holding computed global index
C      ob_subMask_size  - Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      ob_id            - Integer used for do loop through each ob mask
C      total_valid      - Integer value of total count of valid elements in each open boundary mask/ Count done when looping through arrays sent from processes + local field on mask array
C      num_valid        - Integer value of count of valid elements in each open boundary mask. Count done when looping through mask
C      tag              - integer value for looping through all message tags

C      pid              - Integer value for process ID
C      ierror           - Error code for MPI functions
C      status           - Integer array with values characterizing MPI function i.e. message tag, recipient, sender, etc.

C      debug            - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C      print_debug      - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore


C      Although we read REAL*4 values, most built in functions have _RS input which by default is _RL.
C      So global_mask, ob_subMask, and bath_subMask are set as _RL or _RS.
C       _RL global_mask(Nx,Ny,nSx,nSy)
       _RL recv_arr2D(sNx + sNy)
       _RL recv_arr3D(Nr, sNx + sNy)

       _RL g_indices(sNx + sNy)

       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, k, fld, iG, jG, bi, bj
       integer global_ob_index, ob_index, new_goi
       integer ob_subMask_size, counter, ob_id
       integer total_valid, num_valid, tag

       integer, parameter :: debug = 1
       integer, parameter :: print_debug = 0

       logical DIFFERENT_MULTIPLE

#ifdef ALLOW_USE_MPI
       integer pid, ierror, status(MPI_STATUS_SIZE)
#endif

C       common / local_var / global_mask, total_valid
       common / local_var / total_valid

C-------------------------------------------------------------------------------

        DO i=1, sNx + sNy
          recv_arr2D(i) = 0.0
        ENDDO

        DO k=1, Nr
         DO i=1, sNx + sNy
           recv_arr3D(k,i) = 0.0
         ENDDO
        ENDDO

        if (debug .eq. 1) then
          print *, "OB_EXTRACT_OUTPUT timepassed:", time_passed
          print *, "OB_EXTRACT_OUTPUT nTimeSteps_ob:",nTimeSteps_ob
          print *, "OB_EXTRACT_OUTPUT deltaTClock:", deltaTClock
          print *, "OB_EXTRACT_OUTPUT avgPeriod_ob:", avgPeriod_ob
          print *, "OB_EXTRACT_OUTPUT myTime:", myTime
          print *, "OB_EXTRACT_OUTPUT startTime:",startTime
          print *, "OB_EXTRACT_OUTPUT model endTIme", endTime
          print *, "OB_EXTRACT_OUTPUT user start time:",startTime_ob
          print *, "OB_EXTRACT_OUTPUT user end time:",endTime_ob
          print *, "OB_EXTRACT_OUTPUT DIFFERENT_MULTIPLE",
     &          DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock)
        end if

C-------------------------------------------------------------------------------
c      All Processes Tasks:
C-------------------------------------------------------------------------------

c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
       ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy

       print *,"OB_EXTRACT_OUTPUT: fld_nDim",fld_nDim
       print *,"OB_EXTRACT_OUTPUT: fld_to_output",fld_to_output
       fld_to_output(1) = 1

       CALL SET_SUBFIELDS(myThid)
       if (debug .eq. 1) then
         print *,"OB_EXTRACT_OUTPUT: num valid elements:",
     &              numOBPnts_allproc(1,mpiMyId+1)
         print *,"OB_EXTRACT_OUTPUT: subFieldOnMask set"
         call PRINT_FLOAT_ARR(subFieldOnMask_2D(1,1,:), sNx+sNy)
         call PRINT_FLOAT_ARR(subFieldOnMask_3D(1,1,:,:), 2*(sNx+sNy))
       end if


       IF (nTimeSteps_ob .eq. 0) then
C        Initialize averages to zero
         CALL SET_TO_INI_STATE(myThid)
         if (debug .eq. 1) then
           print *,"OB_EXTRACT_OUTPUT: subFieldOnMask_avg ini"
C           call PRINT_FLOAT_ARR(subFieldOnMask_avg(1,1,:), sNx+sNy)

         endif
       ELSE
C          Cumulate Diagnostics_ob fields and number of time steps (for Time Average)
         CALL CUMULATE_FLD_TIME(myThid)
         if (debug .eq. 1) then
           print *,"OB_EXTRACT_OUTPUT: subFieldOnMask_avg cumulate"
C           call PRINT_FLOAT_ARR(subFieldOnMask_avg(1,1,:), sNx+sNy)
         endif
       ENDIF

       print *, "nTimeSteps_ob",nTimeSteps_ob
       print*, "D",DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock)
       print*, "time_passed",time_passed
C      IF (((myIter.NE.nIter0)


       IF (((nTimeSteps_ob .NE. 0)
     &       .AND.DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock))
     &       .AND. (myTime .le. endTime_ob)) THEN
         print *, "HERE "
C        Compute time average of all sub fields
         CALL TIMEAVE_OB_FLD(myThid, myTime)
         if (debug .eq. 1) then
           print *,"OB_EXTRACT_OUTPUT: subFieldOnMask_avg avg"
C           call PRINT_FLOAT_ARR(subFieldOnMask_avg(1,1,:), sNx+sNy)
         endif


C          DO ob_id = 1, nOB_mask

C          Initializing counting integers
C           num_valid = 0

C           iG = mpi_myXGlobalLo(mpiMyId+1)
C           jG = mpi_myYGlobalLo(mpiMyId+1)
C           global_ob_index = Nx*(jG -1)  + iG

C           counter = 2
C           new_goi = global_ob_index
C          Computing global indices, storing values on open boundary points
C           DO bj = myByLo(myThid), myByHi(myThid)
C            DO bi = myBxLo(myThid), myBxHi(myThid)
C              DO j=1,sNy
C                DO i=1,sNx
C                  new_goi = (j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
C     &                       j*(iG-1)+i
C                  if (ob_subMask(ob_id,i,j,bi,bj) .ge. 1.0 ) then
C                    DO fld=1, nOB_fld
C                     subFieldOnMask(fld,counter) =
C     &                                     subField_avg(fld,i,j,bi,bj)
C                    ENDDO
C                    sub_global_indices(counter) = new_goi
C                    counter = counter + 1
C                    num_valid = num_valid + 1
C                  end if
C                 ENDDO
C               ENDDO
C             ENDDO
C           ENDDO

C          Inserting number of valid elements as first element of each array
C           DO fld=1,nOB_fld
C             subFieldOnMask(fld,1) = num_valid
C           ENDDO
C           sub_global_indices(1) = num_valid
c1234567
        DO ob_id=1, nOB_mask
         DO pid=1,nPx*nPy
          print *, "numOBPnts_allproc process", pid-1
          print *,"numOBPnts_allproc:",numOBPnts_allproc(ob_id,pid)
         ENDDO
        ENDDO
        print *, "numOBPnts_allproc at 1st index",
     &    numOBPnts_allproc(1,1)

        print *, "numOBPnts_allproc at 2nd index",
     &    numOBPnts_allproc(1,2)


        DO ob_id=1, nOB_mask

#ifdef ALLOW_USE_MPI
         if (mpiMyId .gt. 0) then
C        Processes 1 and onward send sub fields to Parent process
          do fld=1, nOB_fld
           print *, "num valid:",numOBPnts_allproc(ob_id,mpiMyId+1)
           if (numOBPnts_allproc(ob_id,mpiMyId+1) .ge. 1) then
            print *,"Yes going to send"
            if (fld .le. nOB_fld2D) then
             call MPI_SEND(subFieldOnMask_2Davg(ob_id,fld,:),sNx+sNy
     &              ,MPI_REAL,0, fld, MPI_COMM_MODEL, ierror)
            else
              call MPI_SEND(
     &                subFieldOnMask_3Davg(ob_id,fld-nOB_fld2D,:,:),
     &                fld_depth(fld)*(sNx+sNy),MPI_REAL,0, fld*10,
     &                MPI_COMM_MODEL, ierror)
            end if
           end if
          end do
         end if
C123456789
#endif

C-------------------------------------------------------------------------------
C           Master Process Tasks:
C-------------------------------------------------------------------------------

#ifdef ALLOW_USE_MPI
C123456789
         IF (mpiMyId .eq. 0) then
          ob_index = 1

C           Appending subFields to final global_ob array
C           Field numbers, fld, are used as message tags.
C           DO tag=2, nOB_fld + 1
          DO fld=1, nOB_fld
           total_valid = 0
C          First append process 0 local arrays to final global array if number of valid elements in those local arrays > 0:
C          Indices referring to process IDs begin with 1
           if (numOBPnts_allproc(ob_id,1).ge.1) then
            if (fld .le. nOB_fld2D) then
             call APPEND_OB_VALUE2D(ob_id, 0,
     &                              subFieldOnMask_2Davg(ob_id,fld,:))
            else
C            append_ob_values3D params: ob_id, pid, fld, ob_values
             call APPEND_OB_VALUE3D(ob_id, 0, fld,
     &                    subFieldOnMask_3Davg(ob_id,fld-nOB_fld2D,:,:))
            endif
            total_valid = total_valid +
     &                               numOBPnts_allproc(ob_id, 1)
            if (debug .eq. 1) then
             print *,"APPEND_OB_VALUE: appending local "
     &                    ,"field to global_ob with ",
     &                numOBPnts_allproc(ob_id, 1),"valid elements:"
C            Indices referring to process IDs begin with 1
             print *,"num_ob_points(ob_id)",num_ob_points(ob_id)
             call PRINT_FLOAT_ARR(global_ob2D, total_valid)
            end if
           end if

C             Receive all sub arrays from processes with pid > 0 and valid elements > 0.
C             Append receive sub arrays to final array.
           DO pid=1, mpiNprocs-1
C              Indices referring to process 0 begin with 1
C              If number of valid elements is greater than or equal to 1
            IF (numOBPnts_allproc(ob_id,pid+1) .ge. 1) then
             if (fld .le. nOB_fld2D) then
              call MPI_RECV(recv_arr2D, sNx+sNy, MPI_REAL,
     &                   pid, fld ,MPI_COMM_MODEL, status, ierror)

C                 if (debug .eq. 1) then
C                   print *, "Received from process", pid
C                   print *, "Received recv_arr for fields"
C                   print *, "shape:", shape(recv_arr)
C                   print *, "Displaying FULL array, beyond number of "
C     &                       ,"valid elements"
C                   call PRINT_FLOAT_ARR(recv_arr, sNx+sNy)
C                 end if

              call APPEND_OB_VALUE2D(ob_id,pid, recv_arr2D)
            else
              call MPI_RECV(recv_arr3D, fld_depth(fld)*(sNx+sNy),
     &               MPI_REAL, pid, fld ,MPI_COMM_MODEL, status, ierror)
C             append_ob_values3D params: ob_id, pid, fld, ob_values
              call APPEND_OB_VALUE3D(ob_id,pid, fld, recv_arr3D)
              total_valid = total_valid +
     &                       numOBPnts_allproc(ob_id, pid+1)
C                 if (debug .eq. 1) then
C                   print *,"APPEND_OB_VALUE: appending recv_arr "
C     &                    ,"field to global_ob with ",
C     &                numOBPnts_allproc(ob_id, pid+1),"valid elements:"
C                   call PRINT_FLOAT_ARR(global_ob, total_valid)
C                 end if
              endif

C              !ENDIF for number of valid elements
              ENDIF
C             !ENDDO for pid
             ENDDO

C           !ENDDO for fld
           ENDDO
C         !ENDIF for mpiMyId .eq. 0
         ENDIF
C       !ENDDO for ob_id
        ENDDO
#endif

C        Move on to next time step
      ENDIF
         nTimeSteps_ob = nTimeSteps_ob + 1
         time_passed = time_passed + deltaTClock

      RETURN
      END


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      SUBROUTINE SET_TO_INI_STATE(myThid)
C     Initializes average fields to zero
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER i, ob_id, fld

        DO ob_id=1, nOB_mask
          DO fld=1, nOB_fld
            DO i=1, sNx+sNy
              subFieldOnMask_2Davg(ob_id,fld,i) = 0.0
            ENDDO
          ENDDO
        ENDDO

      END

      SUBROUTINE SET_SUBFIELDS(myThid)
C       Stores desired fields chosen by user into array subFields
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


        INTEGER ob_id, myThid, counter, loc_count
        INTEGER i, j, k, bi, bj, fld, opt, numValid, fld_out
        print *,"Inside SET_SUBFIELDS"
C        print *,"SET_SUBFIELDS: shape fld_choice",shape(fld_choice)
        print *,"SET_SUBFIELDS: fld_nDim",fld_nDim
        print *,"SET_SUBFIELDS: fld_to_output",fld_to_output
C        print *,"SET_SUBFIELDS: fld_choice(nFldOpt)",fld_choice(nFldOpt)

        DO ob_id=1, nOB_mask
         numValid = numOBPnts_allproc(ob_id,mpiMyId+1)
         DO fld=1, nOB_fld
          counter = 1
          fld_out = fld_to_output(fld)
          print *,"SET_SUBFIELDS: fld_out",fld_out
          print *,"SET_SUBFIELDS: numValid",numValid
          DO loc_count=1, numValid
           DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
             j = sub_local_ij_ob(ob_id,2,loc_count)
             i = sub_local_ij_ob(ob_id,1,loc_count)
C             IF (fld_nDim(fld) .eq. 2) then
             IF (fld_nDim(fld) .eq. 2) then
              IF (fld_out .eq. 1) then
               subFieldOnMask_2D(ob_id,fld,counter) = AREA(i,j,bi,bj)
              ELSE IF (fld_out .eq. 2) then
               subFieldOnMask_2D(ob_id,fld,counter) = HEFF(i,j,bi,bj)
              ELSE IF (fld_out .eq. 3) then
               subFieldOnMask_2D(ob_id,fld,counter) = HSNOW(i,j,bi,bj)
#ifdef SEAICE_VARIABLE_SALINITY
              ELSE IF (fld_out .eq. 4) then
               subFieldOnMask_2D(ob_id,fld,counter) = HSALT(i,j,bi,bj)
#endif /* SEAICE_VARIABLE_SALINITY */
              ELSE IF (fld_out .eq. 5) then
               subFieldOnMask_2D(ob_id,fld,counter) = UICE(i,j,bi,bj)
              ELSE IF (fld_out .eq. 6) then
               subFieldOnMask_2D(ob_id,fld,counter) = VICE(i,j,bi,bj)
              ELSE IF (fld_out .eq. 7) then
               subFieldOnMask_2D(ob_id,fld,counter) = etaH(i,j,bi,bj)
              ELSE IF (fld_out .eq. 8) then
               subFieldOnMask_2D(ob_id,fld,counter) = etaN(i,j,bi,bj)
C              ELSE IF (fld_out .eq. 9) then
C               subFieldOnMask_2D(ob_id,fld,counter) = obnw(i,j,bi,bj)
              ENDIF

            ELSE IF (fld_nDim(fld) .eq. 3) then
              DO k=1, fld_depth(fld)
               IF (fld_out .eq. 10) then
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                        uVel(i,j,k,bi,bj)
               ELSE IF (fld_out .eq. 11) then
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                        vVel(i,j,k,bi,bj)
               ELSE IF (fld_out .eq. 12) then
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                        theta(i,j,k,bi,bj)
               ELSE IF (fld_out .eq. 13) then
                print *, "TEST ij:",i,j
                print *, "TEST salt before: ",salt(i,j,k,bi,bj)
                print *, "TEST sFOM before: ",
     &                     subFieldOnMask_3D(ob_id,fld,k,counter)
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                       salt(i,j,k,bi,bj)
                print *, "TEST salt after: ",salt(i,j,k,bi,bj)
                print *, "TEST sFOM after: ",
     &                     subFieldOnMask_3D(ob_id,fld,k,counter)

               ELSE IF (fld_out .eq. 14) then
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                        gU(i,j,k,bi,bj)
               ELSE IF (fld_out .eq. 15) then
                subFieldOnMask_3D(ob_id,fld,k,counter) =
     &                                        gV(i,j,k,bi,bj)
               ENDIF
              ENDDO
             ENDIF
             counter = counter + 1
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C        print *, "FINAL sFOM:", subFieldOnMask_3D
        print *,"sFOM final:"
        call PRINT_FLOAT_ARR(subFieldOnMask_3D, 2*(sNx+sNy))


C      DO ob_id=1, nOB_mask
C          numValid = numOBPnts_allproc(ob_id,mpiMyId+1)
C         DO fld=1, nOB_fld
C          DO opt=1, nFldOpt
C           IF (fld_choice(opt)) then
C            counter = 1
C            DO bj = myByLo(myThid), myByHi(myThid)
C             DO bi = myBxLo(myThid), myBxHi(myThid)
C              DO loc_count=1, numValid
C                j = sub_local_ij_ob(ob_id,2,loc_count)
C                i = sub_local_ij_ob(ob_id,1,loc_count)

C                IF (opt .eq. 1) then
C                 subFieldOnMask(ob_id,fld,counter) = etaN(i,j,bi,bj)
C                ELSE IF (opt .eq. 2) then
C                 subFieldOnMask(ob_id,fld,counter) = uVel(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 3) then
C                 subFieldOnMask(ob_id,fld,counter) = vVel(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 4) then
C                 subFieldOnMask(ob_id,fld,counter) = theta(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 5) then
C                 subFieldOnMask(ob_id,fld,counter) = salt(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 6) then
C                 subFieldOnMask(ob_id,fld,counter) = gU(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 7) then
C                 subFieldOnMask(ob_id,fld,counter) = gV(i,j,1,bi,bj)
C                ELSE IF (opt .eq. 8) then
C                 subFieldOnMask(ob_id,fld,counter) = etaH(i,j,bi,bj)
C                ENDIF
C                counter = counter + 1
C              ENDDO
C            ENDDO
C           ENDDO
C          ENDIF
C         ENDDO
C        ENDDO
C       ENDDO
      END

      SUBROUTINE CUMULATE_FLD_TIME(myThid)
C       Cumulate Diagnostics_ob fields
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

      INTEGER myThid
      INTEGER i, k, ob_id, fld, fld2D, fld3D

      DO ob_id=1, nOB_mask
       DO i=1, sNx+sNy
        DO fld2D=1, nOB_fld2D
         subFieldOnMask_2Davg(ob_id,fld2D,i)
     &                      = subFieldOnMask_2Davg(ob_id,fld2D,i)
     &                        + subFieldOnMask_2D(ob_id,fld2D,i)
        ENDDO
        DO fld3D=1, nOB_fld3D
         DO k=1, fld_depth(nOB_fld2D + fld3D)
          subFieldOnMask_3Davg(ob_id,fld3D,k,i)
     &                    = subFieldOnMask_3Davg(ob_id,fld3D,k,i)
     &                      + subFieldOnMask_3D(ob_id,fld3D,k,i)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      END

      SUBROUTINE TIMEAVE_OB_FLD(myThid, myTime)
C       Computes average of Diagnostics_ob fields using the number of time steps taken in current period
C       Write averaged fields to binary file
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
           DO fld2D=1, nOB_fld2D
            subFieldOnMask_2Davg(ob_id,fld2D,i) =
     &              subFieldOnMask_2Davg(ob_id,fld2D,i)/nTimeSteps_ob
           ENDDO
           DO fld3D=1, nOB_fld3D
            DO k=1, fld_depth(nOB_fld2D + fld3D)
             subFieldOnMask_3Davg(ob_id,fld3D,k,i) =
     &              subFieldOnMask_3Davg(ob_id,fld3D,k,i)/nTimeSteps_ob
            ENDDO
           ENDDO
          ENDDO
C         Write subField_avg to binary file
C         CALL WRITE_SUB_BIN(subFieldOnMask_avg(ob_id,fld,i), fld,
C     &                      myTime)
C        Set number of time steps wihtin averaging period back to zero
        print *, "setting nTimeSteps_ob back to zero at avg"
        nTimeSteps_ob = 0

       ENDDO
      END


      SUBROUTINE APPEND_OB_VALUE2D(ob_id, pid, ob_values)
C       Appends field value corresponding to open boundary points to final global array global_ob.
C       Use values of global_indices to extract value in lookup_table at position of value of global_indices.
C       Value in lookup_table used to identify position of current point on open boundary.
C       Store field value at position in final global array global_ob.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

         _RL ob_values(sNx+sNy)
         integer num_valid, i, ob_index, ob_id, pid
         integer, parameter :: print_debug = 1

         if (print_debug .eq. 1) then
          print *, "APPEND_OB_VALUE: sub_glo_indices_allproc:"
          call PRINT_INT_ARR(sub_glo_indices_allproc(ob_id,pid+1,:)
     &                         ,sNx+sNy )
         end if

         num_valid = numOBPnts_allproc(ob_id,pid+1)

         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: num_valid", num_valid
         end if

         do i=1, num_valid
           ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))
C           if (ob_index .le. (sNy+sNx)*(nPx*nPy)) then
C             if (print_debug .eq. 1) then
C               print *, "APPEND_OB_VALUE: ob_index", ob_index
C             print *, "APPEND_OB_VALUE: counter i",i
C               print *, "APPEND_OB_VALUE: ob_value", ob_values(i)
C             end if
             global_ob2D(int(ob_index)) = ob_values(i)
C           end if
       end do
         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: global_ob2D:"
           call PRINT_FLOAT_ARR(global_ob2D, (sNy+sNx)*(nPx*nPy))
         end if

       end

       SUBROUTINE APPEND_OB_VALUE3D(ob_id, pid, fld, ob_values)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

        _RL ob_values(Nr, sNx+sNy)
        integer num_valid, i, k, ob_index, ob_id, pid, fld
        integer, parameter :: print_debug = 1

        num_valid = numOBPnts_allproc(ob_id,pid+1)

        do k=1, fld_depth(fld)
         do i=1, num_valid
          ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))
          global_ob3D(k,int(ob_index)) = ob_values(k,i)
         enddo
       enddo
      end


       SUBROUTINE CREATE_FILENAME(dir, fldNm, ob_id, filename,
     &                            is_child, myTime)
c       Creating filenames for all binary files.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

       integer is_child, ob_id
       _RL myTime
       character*40 filename
       character*4 fldNm
       character*2 dir
       character*3 str_px, str_py
       character*11 str_time

       character str_ob

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       001 for index 1 in x axis
c       002 for index 2 in x axis
       write(str_ob, "(I1)") ob_id
       write(str_time, "(F11.0)") myTime

       if (is_child .eq. 1) then
         write(str_px, "(I3.3)") mpiPx
         write(str_py, "(I3.3)") mpiPy
         write(filename, "(11A)") dir, fldNm,"_avg.",
     &                      str_time,str_px, ".", str_py, ".bin"
       else
c         Creating a filename for global arrays
         write(filename, "(7A)") dir, fldNm, "_on_mask", str_ob,
     &                           "_global.",str_time,"bin"
       end if

       end


       SUBROUTINE WRITE_SUB_BIN(subMask, fld, myTime)
C      Writing all sub arrays read from binary fie containing global array
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

       integer fld
       _RL myTime
       character*2 dir
       character*40 filename
       _RL subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
       integer, parameter :: print_write = 1

       dir = "./"
       call CREATE_FILENAME(dir, ob_fldNames(fld), 0, filename, 1,
     &                      myTime)
       open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(1) subMask
       close(1)
       if (print_write .eq. 1) then
         print *, "WRITE_SUB_BIN: ",ob_fldNames(fld),"process",
     &             mpiMyId, "is written to ",filename
       end if

       end


       SUBROUTINE WRITE_GLOBAL_BIN(total_valid,fld,ob_id,myTime)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

       integer total_valid, ob_id, fld

       character*2 dir
       character*40 filename
       _RL myTime

c      Writing global array to binary file
       dir = "./"
       call CREATE_FILENAME(dir, ob_fldNames(fld), ob_id, filename, 0,
     &                      myTime)
       open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(1) global_ob(1:total_valid)
       close(1)

       end


       SUBROUTINE PLOT_TEST_GLOBAL_OB(myThid, ob_id, fld)
c      Unravels global_ob and plots onto an array with dimensions Nx by Ny
C       include "SIZE.h"
C       include "USER_INPUT.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"
C       _RL global_mask(Nx,Ny,nSx,nSy)
       _RL global_correct((sNy+sNx)*(nPx*nPy))
       _RL global_NxNy(Nx,Ny,nSx,nSy)

       integer i,j, bi, bj, counter, ob_counter, correct
       integer fld, total_valid, ob_id, myThid
       integer, parameter :: plot_on_NxNy = 1
       integer, parameter :: print_correctness = 1
       character*4 fldNm
       character*30 plotTitle

C       common / local_var / global_mask, total_valid
       common / local_var / total_valid

       fldNm = ob_fldNames(fld)

C      Initializing global_correct, used to indicate if elements of global_ob is correct or not
       DO i=1,(sNy+sNx)*(nPx*nPy)
         global_correct(i) = 1.0
       ENDDO

C       Initialize NxNy array:

        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1, Ny
              DO i=1, Nx
                global_NxNy(i,j,bi,bj) = 0.0
              ENDDO
            ENDDO
          ENDDO
        ENDDO

C      If lookup table contains ob mask point location (1-90) then extract corresponding point in global_ob at location
       correct = 1
       counter = 1
       ob_counter = 1
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1, Ny
             DO i=1, Nx
               if (ob_counter .le. total_valid) then
                 if (lookup_table(ob_id,counter).eq.ob_counter) then
C                Inserting field points into NxNy array
                   global_NxNy(i,j,bi,bj) = global_ob(ob_counter)
C                   if (global_mask(i,j,bi,bj) .eq. 0) then
                   if (global_ob_mask(ob_id,i,j,bi,bj) .eq. 0) then
                     correct = 0
                     global_correct(ob_counter) = 0.0
                   end if
                   ob_counter = ob_counter + 1
                 end if
               end if
               counter = counter + 1
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       if (correct .eq. 0) then
         if (print_correctness .eq. 1) then
           print *, "PLOT_TEST_GLOBAL_OB: global_ob final array is ",
     &            "incorrect for ", fldNm, "for ob", ob_id
           call PRINT_FLOAT_ARR(global_correct, (sNy+sNx)*(nPx*nPy))
         end if
       else
         if (print_correctness .eq. 1) then
           print *, "PLOT_TEST_GLOBAL_OB: global_ob final array is ",
     &         "correct for the ", fldNm, " field for ob", ob_id
         end if
         if (plot_on_NxNy .eq. 1) then
           print *, "PLOT_TEST_GLOBAL_OB: plotting global_ob ",
     &              "on NxNy dimension array:"
           write(plotTitle, "(2A)") "PLOT_TEST_GLOBAL_OB: ", fldNm
           call PLOT_GLO_FIELD_XYRL(global_NxNy,plotTitle,myThid)
         end if
       end if
       end


       SUBROUTINE PRINT_INT_ARR(array, arr_length)
c      Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
        integer arr_length, i
        INTEGER array(arr_length)
        do, i=1, arr_length
          write(*,"(I5,$)") array(i)
            if (i .EQ. arr_length) then
               write(*, '(A,/)') ''
            endif
        enddo
        print *, " "
        end


       SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
c      Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
        integer arr_length, i
        _RL array(arr_length)
        do, i=1, arr_length
          write(*,"(F9.3,$)") array(i)
            if (i .EQ. arr_length) then
               write(*, '(A,/)') ''
            endif
        enddo
        print *, " "
        end



       SUBROUTINE PLOT_GLO_FIELD_XYRL(fld, fldNam, myThid)

C      Building on top of subroutine PLOT_FIELD_XYRS:
C     /==========================================================\
C     | SUBROUTINE PLOT_FIELD_XYRS                               |
C     | Print out an XY _RS field using text map.                |
C     |==========================================================|
C     | This routine references "numerical model" parameters like|
C     | like the integration time. It uses these to create a     |
C     | title for the field before calling a generic execution   |
C     | environment support routine.                             |
C     | This routine can also be edited to cause only some region|
C     | of a field to be printed by default, or every other      |
C     | point etc..                                              |
C     | Other plot formats can also be substituted here.         |
C     | _RS is usually REAL*4                                    |
C     \==========================================================/
       IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     == Routine arguments ==
C     fld - Field to plot
C     fldNam - Name of field
C     myIter - Iteration number for plot
C     myThid - Thread id of thread instance calling plot_field
       _RL fld(Nx, Ny, nSx, nSy)

C       _RL fld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       CHARACTER*(*) fldNam
       INTEGER myThid

C     == Local variables ==
       CHARACTER*(MAX_LEN_MBUF) fldTitle
       INTEGER iStart, iEnd, iStride
       INTEGER jStart, jEnd, jStride
       INTEGER kStart, kEnd, kStride
       INTEGER biStart, biEnd, biStride
       INTEGER bjStart, bjEnd, bjStride

C--   To get around synchronisation and multi-threaded I/O issues
C--   thread 1 will do all the writes.
       _BARRIER
C       IF ( myThid .EQ. 1 ) THEN
C--    Form name for identifying "plot"
        WRITE(fldTitle,'(A,A)') '// Field ', fldNam
C       ENDIF
C--    Do "plot" using textual contour map "execution environment" routine
C      Substitute other plotting utilities here!
        iStart   =  1
        iEnd     =  Nx
        iStride  =  1
        jStart   =  Ny
        jEnd     =  1
        jStride  = -1
        kStart   =  1
        kEnd     =  1
        kStride  =  1
        biStart  =  1
        biEnd    =  nSx
        biStride =  1
        bjStart  =  nSy
        bjEnd    =  1
        bjStride = -1
        CALL PRINT_MAPRL(
     &        fld, fldTitle, PRINT_MAP_XY,
     &         1,Nx,1,Ny,1,1,  nSx,  nSy,
     &         iStart,   iEnd,  iStride,
     &         jStart,   jEnd,  jStride,
     &         kStart,   kEnd,  kStride,
     &        biStart,  biEnd, biStride,
     &        bjStart,  bjEnd, bjStride )

       _BARRIER

       RETURN
       END
