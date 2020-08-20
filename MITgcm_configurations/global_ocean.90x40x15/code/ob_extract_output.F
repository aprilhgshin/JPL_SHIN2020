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

C DONE make this use 3D fields
C DONE changing declaratoin in diagnostics_ob.h
C NOT DOING order of fields correspond to depth too. like depth 1 is eta. deptha 2 - 51 contains theta, etc.
C DONE So subFldsOnMask (nOB_fld, sNx+sNy, nZ )

C fields to input :
C DONE Include HEFF, AREA, HSNON, HSALT, uice, vice from seaice.h
C OBNw from OBCS.h

C DONE dont use the word phase
C DONE change totPhase to having start time and end time instead

C DONE Have user specify what depth it wants
C DONE In common block of diagnostics_ob.h, create an array of depth dimensions of each field

C DONE change k to i --- fix all do loops with this done like this

C DONEchange data file comment line 42, and change 55 to have 7 days + 10 hours

C DONE change fld names to be lenght 8

C DONE Make do loop only go through 1 through sNx etc. not -OLx

C DONE change data.diagnostics to output fields 2D and 3D
C DONE Change frequency for both to be 1day  i.e. 86400 seconds

C DONE fields should not be dependent on any order what so ever. When looping for
C DONE cumulating, compare strings (names of fields) to match.

C DONE A mask is an open boundary

C DONE Try to follow diagnostics package as closely as possible. -- assume biggest size we need for arrays

C DONE in data.diagnostics_ob array of field names chosen by user with one of the
C DONE dimensions corresponding to the mask number

C DONE dimensions in data.diagnostics
C DONE delete obnw

C TODO:
C 1. Write to binary files using MITGCM's existing subroutine
C 2. apply program sturcutre on diagnostics output and compare with my outputs
C 3. Very last addition: passive tracers
C 4. Delete all diagnostics_ob files and variables not used





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
C      msgBuf           - Used to write statements into the output file
C      ob_subMask_size  - Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      ob_id            - Integer used for do loop through each ob mask
C      total_valid      - Integer value of total count of valid elements in each open boundary mask/ Count done when looping through arrays sent from processes + local field on mask array
C      pid              - Integer value for process ID
C      ierror           - Error code for MPI functions
C      status           - Integer array with values characterizing MPI function i.e. message tag, recipient, sender, etc.
C      debug            - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C      print_debug      - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore


C      Although we read REAL*4 values, most built in functions have _RS input which by default is _RL.


       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, k, iG, jG, bi, bj, fld2D, fld3D, fld
       integer ob_subMask_size, ob_id, tag
       integer total_valid
       real*8 buffer2D(sNx+sNy)
       real*8 buffer3D(sNx+sNy, Nr)

       integer, parameter :: debug = 1
       integer, parameter :: print_debug = 0
       integer, parameter :: send_debug = 0


       logical DIFFERENT_MULTIPLE

#ifdef ALLOW_USE_MPI
       integer pid, ierror, status(MPI_STATUS_SIZE)
       common / mpi_var / status, ierror
C  All MPI objects (e.g., MPI_Datatype, MPI_Comm) are of type INTEGER in Fortran.
       integer full_dim(2), sub_dim(2), start_ind(2)
#endif

       common / local_var / total_valid



       IF (myTime .le. endTime_ob) then

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
        print *,"OB_EXTRACT_OUTPUT: ob_nFlds2D",ob_nFlds2D
        print *,"OB_EXTRACT_OUTPUT: ob_nFlds3D",ob_nFlds3D
        print *,"OB_EXTRACT_OUTPUT: ob_levels3D",ob_levels3D

        DO i=1,sNx+sNy
          buffer2D(i) = 0.0
          DO k=1,Nr
            buffer3D(i,k) = 0.0
          ENDDO
        ENDDO

        DO ob_id=1, nOB_mask
         DO i=1, ob_nFlds3D(ob_id)
          IF (ob_levels3D(i,ob_id) .gt. 90) then
           ob_levels3D(i,ob_id) = ob_levels3D(i,ob_id) - 90
          ENDIF
         ENDDO
        ENDDO
c       ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c       field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
        ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy


        CALL SET_SUBFIELDS(myThid, myIter)
        if (print_debug .eq. 1) then
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

C        print *, "nTimeSteps_ob",nTimeSteps_ob
C        print*, "D",DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock)
C        print*, "time_passed",time_passed

C      !ENDIF for IF myTime .le. endTIme_ob
       ENDIF


C------If at end of averaging period:

       IF (((nTimeSteps_ob .NE. 0)
     &       .AND.DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock))
     &       .AND. (myTime .le. endTime_ob)) THEN

C       Compute time average of all sub fields
        CALL TIMEAVE_OB_FLD(myThid, myTime)


        IF (print_debug .eq. 1) then
         DO ob_id=1, nOB_mask
           print*,"SFOM3Davg fld1:",subFieldOnMask_3Davg(ob_id,1,:,1),
     &          "at ob_id", ob_id
         ENDDO
        ENDIF




C       Now looping through each mask and communicating fields between master process and child processes
        DO ob_id=1, nOB_mask


#ifdef ALLOW_USE_MPI
         if (mpiMyId .gt. 0) then

C          print *,"OB_EXTRACT_OUTPUT: Process",mpiMyId,"at mask",ob_id
C          print *,"OB_EXTRACT_OUTPUTnumOBPnts_allproc(ob_id,mpiMyId+1)"
C     &             ,numOBPnts_allproc(ob_id,mpiMyId+1)

C        Processes 1 and onward send sub fields to Parent process
C        Field numbers, fld2D and fld3D, are used as message tags.

          if (numOBPnts_allproc(ob_id,mpiMyId+1) .ge. 1) then

C           print *, "num valid:",numOBPnts_allproc(ob_id,mpiMyId+1)
           print *,"Yes going to send for mask", ob_id

           do fld2D=1, ob_nFlds2D(ob_id)
            call OB_PASS_RL_to_R8(buffer2D,buffer3D,ob_id,2,fld2D)
            tag = ob_tags(ob_id, 1, fld2D)
            if (send_debug .eq. 1) then
             print *,"Sending field", ob_flds2D(fld2D,ob_id)
             print*,"Sending with tag:",tag,"for ob_id",ob_id,
     &             "fld2D",fld2D, "at myIter", myIter, "proc",mpiMyId
            endif
C            call MPI_SEND(subFieldOnMask_2Davg(ob_id,fld2D,:),sNx+sNy
C     &      ,MPI_REAL8,0, tag, MPI_COMM_MODEL, ierror)
            call MPI_SEND(buffer2D, sNx+sNy, MPI_REAL8, 0, tag,
     &                    MPI_COMM_MODEL, ierror)

            if (send_debug .eq. 1) then
             print *,"Sending array:",ob_flds2D(fld2D,ob_id)
             print *,buffer2D
            endif

           end do

           do fld3D=1, ob_nFlds3D(ob_id)
            tag = ob_tags(ob_id, 2, fld3D)
            call OB_PASS_RL_to_R8(buffer2D,buffer3D,ob_id,3,fld3D)
            if (send_debug .eq. 1) then
             print *,"Sending field", ob_flds3D(fld3D,ob_id), "for mask"
     &              , ob_id
             print*,"Sending with tag:",tag,"for ob_id",ob_id,
     &             "fld3D",fld3D, "at myIter", myIter, "proc",mpiMyId
             print *,"Sending array shape: ", shape(buffer3D)
            endif
C            call MPI_SEND(subFieldOnMask_3Davg(ob_id,fld3D,:,:),
C     &   ob_levels3D(fld3D,ob_id)*(sNx+sNy), MPI_REAL8, 0, tag,
C     &   MPI_COMM_MODEL, ierror)
            call MPI_SEND(buffer3D, ob_levels3D(fld3D,ob_id)*(sNx+sNy),
     &                    MPI_REAL8, 0, tag, MPI_COMM_MODEL, ierror)

            if (send_debug .eq. 1) then
             DO k=1, ob_levels3D(fld3D,ob_id)
              print *,"printing at level k",k, "with tag", tag
              print*,buffer3D(:,k)
             ENDDO
            endif

            print *,"FINISHED SENDING."

           end do
          end if
         end if

#endif


C-------------------------------------------------------------------------------
C           Master Process Tasks:
C-------------------------------------------------------------------------------
         IF (mpiMyId .eq. 0) then

C           print *,"HELLO HERE"
C           DO pid=1, mpiNprocs-1
C            DO fld=1, ob_nFlds3D(ob_id)
C             print *,"HELLO HERE IN fld",fld


C         ENDDO
C        ENDDO

C          print *,"OB_EXTRACT_OUTPUT: MASK", ob_id,
C     &              "Before calling on master tasks"
C          print *,"numOBPnts_allproc(ob_id,mpiMyId+1)",
C     &             numOBPnts_allproc(ob_id,mpiMyId+1)
          DO fld2D=1, ob_nFlds2D(ob_id)
            call MASTER_PROC_TASKS(fld2D, ob_id, 2, myThid, myTime,
     &                             myIter)
          ENDDO

          DO fld3D=1, ob_nFlds3D(ob_id)
            call MASTER_PROC_TASKS(fld3D, ob_id, 3, myThid, myTime,
     &                             myIter)
          ENDDO

C         !ENDIF for mpiMyId .eq. 0
         ENDIF
C       !ENDDO for ob_id
        ENDDO
C     !ENDIF for end of averaging period

C      Set averages to back zero
       CALL SET_TO_INI_STATE(myThid)
C       if (debug .eq. 1) then
C         print *,"OB_EXTRACT_OUTPUT: subFieldOnMask_avg ini"
C           call PRINT_FLOAT_ARR(subFieldOnMask_avg(1,1,:), sNx+sNy)
C       endif

      ENDIF
C      Move on to next time step
       nTimeSteps_ob = nTimeSteps_ob + 1
       time_passed = time_passed + deltaTClock

      RETURN
      END




C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      SUBROUTINE OB_PASS_RL_to_R8(buffer2D, buffer3D, ob_id, nDim,fld)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_OB.h"

        REAL*8 buffer2D(sNx + sNy)
        REAL*8 buffer3D(sNx + sNy, Nr)
        INTEGER ob_id, nDim, fld
        INTEGER i,k
        INTEGER, PARAMETER :: debug = 0

        IF (nDim .eq. 2) then
          DO i=1, sNx+sNy
            buffer2D(i) = subFieldOnMask_2Davg(ob_id,fld,i)
          ENDDO

        ELSE IF (nDim .eq. 3) then
          DO i=1, sNx+sNy
            DO k=1,ob_levels3D(fld,ob_id)
              buffer3D(i,k) = subFieldOnMask_3Davg(ob_id,fld,i,k)
            ENDDO
          ENDDO

          if (debug .eq. 1) then
           DO k=1, ob_levels3D(fld,ob_id)
            print *,"At level k:",k
            print*, "buffer3D:",buffer3D(:,k)
            print*,"subFieldOnMask_3Davg:",
     &                   subFieldOnMask_3Davg(ob_id,fld,:,k)
           ENDDO
          endif




        ENDIF

      END


      SUBROUTINE MASTER_PROC_TASKS(fld, ob_id, nDim, myThid, myTime,
     &                             myIter)
C     Appending subFields to final global_ob array
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

C     Input Parameters:
C     fld :: current field number
C     ob_id :: current mask number
C     nDim :: number of dimensions of field to output
C     myThid :: thread number
C     myTime :: current time in model
      INTEGER fld, ob_id, nDim, myThid, myIter
      _RL myTime
      REAL*8 buffer2D(sNx+sNy)
      REAL*8 buffer3D(sNx+sNy, Nr)

      _RS dummyRS(1)
      INTEGER pid, status(MPI_STATUS_SIZE), ierror, tag
      INTEGER total_valid, i, k, count
      INTEGER, PARAMETER :: debug = 0
      INTEGER, PARAMETER :: recv_debug = 0
      integer, parameter :: write_to_binary = 0


      common / local_var / total_valid
      common / mpi_var / status, ierror

      DO i=1,sNx+sNy
        buffer2D(i) = 0.0
        DO k=1,Nr
          buffer3D(i,k) = 0.0
        ENDDO
      ENDDO

      total_valid = 0

C-----Append Local Array--------------------------------------------------------
C     - First append process 0 local arrays to final global array if number of
C          valid elements in those local arrays > 0
C     - Indices referring to process IDs begin with 1 in numOBPnts_allproc and
C          subFieldOnMask_2Davg
C-------------------------------------------------------------------------------
      if (debug .eq. 1) then
       print *,"MASTER: numOBPnts_allproc ",
     &            numOBPnts_allproc(ob_id,1), "at ob_id",ob_id
      endif

      IF (numOBPnts_allproc(ob_id,1).ge.1) then
       if (nDim .eq. 2) then

        print *,"MASTER2D: 2 dimensional field appending local field",
     &           ob_flds2D(fld,ob_id)
        call OB_PASS_RL_to_R8(buffer2D, buffer3D, ob_id, 2,fld)
        call APPEND_OB_VALUE2D(ob_id, 0, buffer2D)
        total_valid = total_valid + numOBPnts_allproc(ob_id, 1)

        if (debug .eq. 1) then
         print *,"MASTER 2D: mask number:", ob_id
         print *,"MASTER 2D: appended local field to global array"
     &          ," with ",numOBPnts_allproc(ob_id, 1),"valid elements:"
C         print *,"MASTER 2D: num_ob_points(ob_id)",num_ob_points(ob_id)
         print *,"MASTER 2D: global_ob2D: with total_valid",total_valid
         call PRINT_FLOAT_ARR(global_ob2D, total_valid)
        end if

       else if (nDim .eq. 3) then
        print *,"MASTER3D: 3 dimensional field appending local field",
     &           ob_flds3D(fld,ob_id)
        call OB_PASS_RL_to_R8(buffer2D, buffer3D, ob_id, 3,fld)

        call APPEND_OB_VALUE3D(ob_id, 0, fld, buffer3D)
        total_valid = total_valid + numOBPnts_allproc(ob_id, 1)

        if (debug .eq. 1) then
         print *,"MASTER 3D: mask number:", ob_id
         print *,"MASTER 3D: appended local field to global array"
     &          ," with ",numOBPnts_allproc(ob_id, 1),"valid elements:"
C        print *,"MASTER 3D: num_ob_points(ob_id)",num_ob_points(ob_id)
         do k=1, ob_levels3D(fld,ob_id)
          print *,"MASTER 3D: global_ob3D at depth",k,"with total_valid"
     &            ,total_valid
          call PRINT_FLOAT_ARR(global_ob3D(:,k), total_valid)
         enddo
        end if
       end if

      ENDIF

C-----End of appending local array----------------------------------------------


C------Receiving arrays through MPI and appending-------------------------------
C      - Receive all sub arrays from processes with pid>0 and valid elements>0.
C      - Append receive sub arrays to final array.
C-------------------------------------------------------------------------------
#ifdef ALLOW_USE_MPI

       DO pid=1, mpiNprocs-1
C       If number of valid elements is greater than or equal to 1
        IF (numOBPnts_allproc(ob_id,pid+1) .ge. 1) then

         if (debug .eq. 1) then
          print *,"MASTER: DO pid loop to receive with MPI"
          print*,"MASTER: ",numOBPnts_allproc(ob_id,pid+1),
     &           "elements for mask", ob_id, "at pid",pid
          print *,"MASTER: time to receive",
     &                    numOBPnts_allproc(ob_id,pid+1), "elements"
         endif

         IF (nDim .eq. 2) then
          tag = ob_tags(ob_id, 1, fld)
          call MPI_RECV(buffer2D, sNx+sNy, MPI_REAL8,
     &                   pid, tag ,MPI_COMM_MODEL, status, ierror)

          if (recv_debug .eq. 1) then
           print *, "MASTER 2D Received from process", pid
           print *, "MASTER Received buffer2D for field ",
     &                ob_flds2D(fld, ob_id)
           print *, "MASTER buffer2D shape:", shape(buffer2D)
           print *, "MASTER Displaying FULL 2D array, beyond number of"
     &                       ," valid elements"
           print *,buffer2D
          end if

          call APPEND_OB_VALUE2D(ob_id, pid, buffer2D)
C          print *,"MASTER: global_ob2D",global_ob2D

         ELSE IF (nDim .eq. 3) then

          tag = ob_tags(ob_id, 2, fld)

          call MPI_RECV(buffer3D(:,1:ob_levels3D(fld,ob_id)),
     &           ob_levels3D(fld,ob_id)*(sNx+sNy), MPI_REAL8, pid,
     &           tag ,MPI_COMM_MODEL, status, ierror)

          if (recv_debug .eq. 1) then
            print *, "MASTER 3D Received from process",pid,"tag",tag
            print *, "MASTER Received buffer3D for field ",
     &                ob_flds3D(fld, ob_id)
            print *, "MASTER buffer3D shape:", shape(buffer3D)
            print *, "MASTER Displaying FULL 3D array, beyond number of"
     &                       ," valid elements"
            DO k=1, ob_levels3D(fld,ob_id)
             print *,"printing at level k", k
             print*,buffer3D(:,k)
           ENDDO
          end if

C         !append_ob_values3D params: ob_id, pid, fld, ob_values
          call APPEND_OB_VALUE3D(ob_id,pid, fld, buffer3D)
          print *,"MASTER: global_ob3D",global_ob3D

         ENDIF
         total_valid = total_valid + numOBPnts_allproc(ob_id, pid+1)

C       !ENDIF for number of valid elements
        ENDIF

        if (debug .eq. 1) then
         print *,"MASTER: total_valid", total_valid
        end if

C      !ENDDO for pid
       ENDDO
#endif
C------End of receiving arrays through MPI and appending------------------------

       if (debug .eq. 1) then
        print *,"MASTER: For mask", ob_id
        print *,"Plotting final array at time", myTime
       end if
CC       call PLOT_TEST_GLOBAL_OB(myThid, ob_id, fld, nDim, myIter)

       if (debug .eq. 1) then
        print*,"MASTER: myIter:",myIter
       end if
       call WRITE_GLOBAL_BIN(total_valid,fld,ob_id,myIter,nDim)


       IF (write_to_binary .eq. 1) then
        IF (nDim .eq. 2) then
         call MDS_WRITE_FIELD(ob_flds2D(fld,ob_id), 32, .FALSE.,
     &    .TRUE., 'RL',1, 1, 1,
     &    global_ob2D(1:total_valid), dummyRS,
     &    1, myIter, myThid)
        ELSE IF (nDim .eq. 3) then
         call MDS_WRITE_FIELD(ob_flds3D(fld,ob_id), 32, .FALSE.,.TRUE.,
     &    'RL',ob_levels3D(fld,ob_id), 1, ob_levels3D(fld,ob_id),
     &    global_ob3D(1:total_valid, 1:ob_levels3D(fld,ob_id)), dummyRS,
     &    1, myIter, myThid)
        ENDIF
       ENDIF


C     I   fName,
C     I   filePrec,
C     I   globalFile,
C     I   useCurrentDir,
C     I   arrType,
C     I   kSize,kLo,kHi,
C     I   fldRL, fldRS,
C     I   jrecord,
C     I   myIter,
C     I   myThid )
C       WRITE_REC_XY_RL( fName,field,iRec,myIter,myThid )
C       WRITE_REC_XYZ_RL( fName,field,iRec,myIter,myThid )

      END


      SUBROUTINE SET_TO_INI_STATE(myThid)
C     Initializes average fields to zero
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER nFlds2D, nFlds3D
        INTEGER i, k, ob_id, fld2D, fld3D

        DO ob_id=1, nOB_mask
         nFlds2D = ob_nFlds2D(ob_id)
         DO fld2D=1, nFlds2D
          DO i=1, sNx+sNy
           subFieldOnMask_2Davg(ob_id,fld2D,i) = 0.0
          ENDDO
         ENDDO

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
      ENDIF

      tag = 1
      DO ob_id=1, nOB_mask
        numValid = numOBPnts_allproc(ob_id,mpiMyId+1)

C      Looping through 2D fields first, if any.
       nFlds2D = ob_nFlds2D(ob_id)
       DO fld2D=1, nFlds2D
        ob_tags(ob_id, 1, fld2D) =  tag
        IF (print_debug .eq. 1) then
         print *,"SET_SUBFIELDS: tag",tag, "in myIter",myIter
        ENDIF
        tag = tag + 1
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO loc_count=1, numValid
           j = sub_local_ij_ob(ob_id,2,loc_count)
           i = sub_local_ij_ob(ob_id,1,loc_count)
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
C            print *,"SET_SUBFIELDS: ij",i,j,"etan:",etaN(i,j,bi,bj)
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = etaN(i,j,bi,bj)
C            print *,"SET_SUBFIELDS:subFieldOnMask_2D",
C     &                subFieldOnMask_2D(ob_id,fld2D,loc_count)
           ELSE IF (ob_flds2D(fld2D,ob_id) .eq. 'ETAH') then
            subFieldOnMask_2D(ob_id,fld2D,loc_count) = etaH(i,j,bi,bj)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C      Looping through 3D fields, if any.
       nFlds3D = ob_nFlds3D(ob_id)
       DO fld3D=1, nFlds3D
        ob_tags(ob_id, 2, fld3D) = tag
        IF (print_debug .eq. 1) then
         print *,"SET_SUBFIELDS: tag",tag, "in myIter",myIter
        ENDIF
        tag = tag + 1
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO loc_count=1, numValid
           j = sub_local_ij_ob(ob_id,2,loc_count)
           i = sub_local_ij_ob(ob_id,1,loc_count)
           DO k=1, ob_levels3D(fld3D,ob_id)
C            print *,"SET_SUBFIELDS: k",k,"loc_count",loc_count
            IF (ob_flds3D(fld3D,ob_id) .eq. 'UVEL') then
             subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        uVel(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'VVEL') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        vVel(i,j,k,bi,bj)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'THETA') then
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k) =
     &                                        theta(i,j,k,bi,bj)
            print*,"DOUBLE CHECKING THETA:"
            print*,"orig then copy",theta(i,j,k,bi,bj),
     &            subFieldOnMask_3D(ob_id,fld3D,loc_count,k)
            ELSE IF (ob_flds3D(fld3D,ob_id) .eq. 'SALT') then
C              print *,"SET_SUBFIELDS: ij",i,j,"salt:",salt(i,j,k,bi,bj)
              subFieldOnMask_3D(ob_id,fld3D,loc_count,k)=
     &                                       salt(i,j,k,bi,bj)
C              print *,"SET_SUBFIELDS:subFieldOnMask_3D",
C     &                subFieldOnMask_3D(ob_id,fld3D,loc_count,k)
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
       ENDDO
      ENDDO

      IF (print_debug .eq. 1) then
       print *, "FINAL set 2D:"
       call PRINT_FLOAT_ARR(subFieldOnMask_2D(1,1,:), sNx+sNy)
       print *,"FINAL set 3D first level:"
       call PRINT_FLOAT_ARR(subFieldOnMask_3D(1,1,:,1), (sNx+sNy))
      ENDIF

      END

      SUBROUTINE CUMULATE_FLD_TIME(myThid)
C       Cumulate Diagnostics_ob fields
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

      INTEGER myThid
      INTEGER i, k, ob_id, fld, fld2D, fld3D
      INTEGER, PARAMETER :: print_debug = 0

      IF (print_debug .eq. 1) then
       print *,"Cumulating at ",nTimeSteps_ob
      ENDIF

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
        INTEGER, PARAMETER :: debug = 0

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
C       Set number of time steps wihtin averaging period back to zero
       IF (debug .eq. 1) then
        print *, "setting nTimeSteps_ob back to zero at avg"
       ENDIF
       nTimeSteps_ob = 0
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

         REAL*8 ob_values(sNx+sNy)
         integer num_valid, i, ob_index, ob_id, pid
         integer, parameter :: print_debug = 0

         ob_index = 1
         num_valid = numOBPnts_allproc(ob_id,pid+1)

         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE2D: num_valid", num_valid
         end if

         do i=1, num_valid
           ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))
C           if (ob_index .le. (sNy+sNx)*(nPx*nPy)) then
             if (print_debug .eq. 1) then
               print *, "APPEND_OB_VALUE2D: ob_index", ob_index
               print *, "APPEND_OB_VALUE2D: counter i",i
               print *, "APPEND_OB_VALUE2D: ob_value", ob_values(i)
             end if
             global_ob2D(ob_index) = ob_values(i)
C           end if
       end do
         if (print_debug .eq. 1) then
           print *,"global_ob2D shape",shape(global_ob2D)
           print *, "APPEND_OB_VALUE2D: global_ob2D:"
           call PRINT_FLOAT_ARR(global_ob2D, (sNy+sNx)*(nPx*nPy))
         end if

       end

       SUBROUTINE APPEND_OB_VALUE3D(ob_id, pid, fld, ob_values)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

        REAL*8 ob_values(sNx+sNy, Nr)
        integer num_valid, i, k, ob_index, ob_id, pid, fld
        integer, parameter :: print_debug = 0

        ob_index = 1
        num_valid = numOBPnts_allproc(ob_id,pid+1)

        if (print_debug .eq. 1) then
          print *,"APPEND_OB_VALUE3D: pid",pid
        endif

        do i=1, num_valid
         ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))
C         if (print_debug .eq. 1) then
C           print*,"APPEND_OB_VALUE3D: ob_index",ob_index
C         endif
         do k=1, ob_levels3D(fld,ob_id)
          global_ob3D(ob_index,k) = ob_values(i,k)
C          print*,"ob_values",ob_values(i,k),"at k",k,"i",i
         enddo
       enddo

       if (print_debug .eq. 1) then
         do k=1, ob_levels3D(fld,ob_id)
          print *, "APPEND_OB_VALUE: global_ob3D with total depth:",
     &       ob_levels3D(fld, ob_id), "printing for depth #",k
          print*,"APPEND_OB_VALUE ob_values:",ob_values(1:num_valid,k)
          print *,"APPEND_OB_VALUE global_ob3D", global_ob3D(:,k)
         enddo
       end if

      end


       SUBROUTINE CREATE_FILENAME(dir, fldNm, ob_id, filename,
     &                            is_child, myIter)
c       Creating filenames for all binary files.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

       integer is_child, ob_id, myIter
       _RL myTime
       character*50 filename
       character*8 fldNm
       character*2 dir
       character*3 str_px, str_py
       character*11 str_time

       character*2 str_ob
       character*8 str_iter

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       001 for index 1 in x axis
c       002 for index 2 in x axis
       write(str_ob, "(I2.2)") ob_id
       write(str_iter, "(I8.8)") myIter
       write(str_time, "(F11.0)") myTime

       if (is_child .eq. 1) then
         write(str_px, "(I3.3)") mpiPx
         write(str_py, "(I3.3)") mpiPy
         write(filename, "(11A)") dir, fldNm,"_avg.",
     &                      str_time,str_px, ".", str_py, ".bin"
       else
c         Creating a filename for global arrays
C         write(filename, "(7A)") dir, fldNm, "_on_mask", str_ob,
C     &                           "_global.",str_time,".bin"

         write(filename, "(7A)") "MASK_", str_ob, "_", fldNm,
     &          "_", str_iter, ".bin"
       end if

       end


       SUBROUTINE WRITE_GLOBAL_BIN(total_valid,fld,ob_id,myIter,nDim)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"
#include "DIAGNOSTICS_OB.h"

       integer total_valid, ob_id, fld, nDim, myIter

       character*2 dir
       character*50 filename
       real*4 bufferR4_2D((sNy+sNx)*(nPx*nPy))
       real*4 bufferR4_3D((sNy+sNx)*(nPx*nPy), Nr)
       integer, parameter :: debug = 1

c      Writing global array to binary file
       dir = "./"

       IF (ob_filePrec .eq. 32) then

        IF (nDim .eq. 2) then

         call OB_R8_to_R4(bufferR4_2D, bufferR4_3D,2,fld,ob_id)
         call CREATE_FILENAME(dir, ob_flds2D(fld,ob_id), ob_id,filename
     &                      ,0, myIter)
         open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
         write(1) bufferR4_2D(1:total_valid)
         if (debug .eq. 1) then
          print *,"WRITE_GLOBAL_BIN: 2D field",ob_flds2D(fld,ob_id),
     &           "written to file:",filename
          print *, "WRITE_GLOBAL_BIN: bufferR4_2D",
     &     bufferR4_2D(1:total_valid)

         endif

        ELSE IF (nDim .eq. 3) then
         call OB_R8_to_R4(bufferR4_2D, bufferR4_3D,3,fld,ob_id)
         call CREATE_FILENAME(dir, ob_flds3D(fld,ob_id), ob_id, filename
     &                      ,0, myIter)
         open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
         write(1) bufferR4_3D(1:total_valid,1:ob_levels3D(fld,ob_id))
         if (debug .eq. 1) then
          print *,"WRITE_GLOBAL_BIN: 3D field ",ob_flds3D(fld,ob_id),
     &           "written to file:",filename
          print *, "WRITE_GLOBAL_BIN: bufferR4_3D",
     &     bufferR4_3D(1:total_valid,1:ob_levels3D(fld,ob_id))
         endif

        ENDIF

       ELSE IF (ob_filePrec .eq. 64) then

        IF (nDim .eq. 2) then

         call CREATE_FILENAME(dir, ob_flds2D(fld,ob_id), ob_id,filename
     &                      ,0, myIter)
         open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
         write(1) global_ob2D(1:total_valid)
         if (debug .eq. 1) then
          print *,"WRITE_GLOBAL_BIN: 2D field",ob_flds2D(fld,ob_id),
     &           "written to file:",filename
          print *, "WRITE_GLOBAL_BIN: global_ob2D",global_ob2D
         endif
        ELSE IF (nDim .eq. 3) then
         call CREATE_FILENAME(dir, ob_flds3D(fld,ob_id), ob_id, filename
     &                      ,0, myIter)
         open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
         write(1) global_ob3D(1:total_valid,1:ob_levels3D(fld,ob_id))
         if (debug .eq. 1) then
          print *,"WRITE_GLOBAL_BIN: 3D field ",ob_flds3D(fld,ob_id),
     &           "written to file:",filename
          print *, "WRITE_GLOBAL_BIN: global_ob3D",global_ob3D

         endif
        ENDIF

       ENDIF

       close(1)

       end


       SUBROUTINE OB_R8_to_R4(bufferR4_2D, bufferR4_3D,nDim,fld,ob_id)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

         real*4 bufferR4_2D((sNy+sNx)*(nPx*nPy))
         real*4 bufferR4_3D((sNy+sNx)*(nPx*nPy), Nr)
         integer nDim, fld, ob_id
         integer, parameter :: debug = 0
         integer i, k

C        Initialize buffer arrays:
        do i=1, (sNy+sNx)*(nPx*nPy)
          bufferR4_2D(i) = 0.0
          do k=1, Nr
            bufferR4_3D(i,k) = 0.0
          enddo
        enddo

C        Set real*8 values to real*4

         if (nDim .eq. 2) then
           do i=1, (sNy+sNx)*(nPx*nPy)
             bufferR4_2D(i) = global_ob2D(i)
           enddo

         else if (nDim .eq. 3) then
           do i=1, (sNy+sNx)*(nPx*nPy)
             do k=1, ob_levels3D(fld,ob_id)
               bufferR4_3D(i,k) = global_ob3D(i,k)
             enddo
           enddo

           if (debug .eq. 1) then
            do k=1,ob_levels3D(fld,ob_id)
             print *,"OB: bufferR4_3D:",bufferR4_3D(:,k)
             print *,"OB: global_ob3D:",global_ob3D(:,k)
            enddo
           endif

         endif

       END


       SUBROUTINE PLOT_TEST_GLOBAL_OB(myThid, ob_id, fld, nDim, myIter)
c      Unravels global_ob and plots onto an array with dimensions Nx by Ny

#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

       REAL*8 global_correct((sNy+sNx)*(nPx*nPy))
       REAL*8 global_NxNy(Nx,Ny,nSx,nSy)
       _RS dummyRS(1)
       REAL*8 global_theta(Nx,Ny)

       integer i,j, bi, bj, k, counter, ob_counter, correct
       integer fld, total_valid, ob_id, myThid, nDim, myIter
       integer, parameter :: plot_on_NxNy = 1
       integer, parameter :: print_correctness = 1
       integer, parameter :: write_to_binary = 0

       character*8 fldNm
       character*30 plotTitle

C       common / local_var / global_mask, total_valid
       common / local_var / total_valid

       fldNm = ' '

       if (nDim .eq. 2) then
         fldNm = ob_flds2D(fld, ob_id)
       else if (nDim .eq. 3) then
         fldNm = ob_flds3D(fld, ob_id)
       endif

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

        print *, "PLOT_TEST_GLOBAL_OB: total_valid", total_valid

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
C                 !Inserting field points into NxNy array
                  if (nDim .eq. 2) then
                   global_NxNy(i,j,bi,bj) = global_ob2D(ob_counter)
                  else if  (nDim .eq. 3) then
                   global_NxNy(i,j,bi,bj) =
     &             global_ob3D(ob_counter,(ob_levels3D(fld,ob_id)/2)+1)
                  endif
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


       IF (write_to_binary .eq. 1) then
        IF (nDim .eq. 2) then

         print *,"PLOT_TEST: shape2D:",shape(global_NxNy)
         call MDS_WRITE_FIELD(ob_flds2D(fld,ob_id), 32, .TRUE., .TRUE.,
     &    'RL',1, 1, 1,
     &    global_NxNy(:,:,1,1), dummyRS,
     &    1, myIter, myThid)

C        print *,"PLOT_TEST: reading back from file "
C        call MDS_READ_FIELD(ob_flds2D(fld,ob_id), 32, .TRUE.,
C     &    'RL',1, 1, 1,
C     &    global_NxNy(:,:,1,1), dummyRS,
C     &    1, myIter, myThid)


        ELSE IF (nDim .eq. 3) then
         print *,"PLOT_TEST: shape3D L1:",shape(global_NxNy)
         call MDS_WRITE_FIELD(ob_flds3D(fld,ob_id), 32, .TRUE., .TRUE.,
     &    'RL',1, 1, 1,
     &    global_NxNy(:,:,1,1), dummyRS,
     &    1, myIter, myThid)
C        call MDS_WRITE_FIELD(ob_flds3D(fld,ob_id), 32, .TRUE., .TRUE.,
C     &    'RL',1, 1, ob_levels3D(fld,ob_id),
C     &    global_ob3D(1:ob_levels3D(fld,ob_id),1:total_valid), dummyRS,
C     &    1, myIter, myThid)

C     I   fName,
C     I   filePrec,
C     I   globalFile,
C     I   useCurrentDir,
C     I   arrType,
C     I   kSize,kLo,kHi,
C     I   fldRL, fldRS,
C     I   jrecord -- starting record within the output file that will contain the array
C     I   myIter,
C     I   myThid )

        ENDIF
       ENDIF
       END


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
