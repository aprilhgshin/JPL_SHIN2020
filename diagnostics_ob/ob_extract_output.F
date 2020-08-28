C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C

#include "DIAG_OPTIONS.h"

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
C      boundary points in their assigned sub-mask and extract field
C      values at open boundary points and sends to the master process. Master
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
C     startTime/ endTime
C     deltaTClock

C     From EESUPPORT.h : all initialized by subroutine INI_PROCS:
C     mpiNprocs - No. of MPI processes.
C     mpiMyId   - MPI process id of me.
C     mpiPx     - My MPI proc. grid X coord
C     mpiPy     - My MPI proc. grid Y coord
C     mpiXGlobalLo - My bottom-left (south-west) x-coordinate in
C                    global domain.
C     mpiYGlobalLo - My bottom-left (south-west) y-coordinate in
C                    global domain.

C-------------------------------------------------------------------------------
C      !LOCAL VARIABLES:

C      msgBuf           - Used to write statements into the output file
C      ob_subMask_size  - Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      ob_id            - Integer used for do loop through each ob mask
C      total_valid      - Integer value of total count of valid elements in each open boundary mask/ Count done when looping through arrays sent from processes + local field on mask array
C      pid              - Integer value for process ID
C      ierror           - Error code for MPI functions
C      status           - Integer array with values characterizing MPI function i.e. message tag, recipient, sender, etc.
C      debug            - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C      print_debug      - Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C      buffer2D         - Real*8 array storing subFieldOnMask_2Davg values converted to Real*8
C                         Also acts as buffer for MPI_RECV
C      buffer3D         - Real*8 array storing subFieldOnMask_3Davg values converted to Real*8
C                         Also acts as buffer for MPI_RECV

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
C         At each time step:
C          1. Extract field values on open boundary points
C          2. Accumulate extracted field values
C         At end of averaging period:
C          1. Generate time-averaged field values
C          2. Processes > 0 send time-averaged field values to master process
C             (Convert _RL array to R8 before sending. Send R8 array)
C-------------------------------------------------------------------------------

C       Initializing local arrays
        DO i=1,sNx+sNy
          buffer2D(i) = 0.0
          DO k=1,Nr
            buffer3D(i,k) = 0.0
          ENDDO
        ENDDO

c       ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c       field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
        ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy

C       Extracting field values on open boundary points at each time step
        CALL SET_SUBFIELDS(myThid, myIter)


        IF (nTimeSteps_ob .eq. 0) then
C        Initialize averages to zero
         CALL SET_TO_INI_STATE(myThid)
        ELSE
C        Cumulate Diagnostics_ob fields and number of time steps (for Time Average)
         CALL CUMULATE_FLD_TIME(myThid)
        ENDIF

C      !ENDIF for IF myTime .le. endTIme_ob
       ENDIF


C------If at end of averaging period:
       IF (((nTimeSteps_ob .NE. 0)
     &       .AND.DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock))
     &       .AND. (myTime .le. endTime_ob)) THEN

C       Compute time average of all sub fields
        CALL TIMEAVE_OB_FLD(myThid, myTime)

C       Now looping through each mask and communicating fields between master process and child processes
        DO ob_id=1, nOB_mask

#ifdef ALLOW_USE_MPI
         if (mpiMyId .gt. 0) then

C        Processes 1 and onward send sub fields to Parent process
          if (numOBPnts_allproc(ob_id,mpiMyId+1) .ge. 1) then

           print *,"Yes going to send for mask", ob_id

           DO fld2D=1, ob_nFlds2D(ob_id)
C           Convert _RL array to R8 array
            call OB_PASS_RL_to_R8(buffer2D,buffer3D,ob_id,2,fld2D)

            tag = ob_tags(ob_id, 1, fld2D)

            if (send_debug .eq. 1) then
             print *,"Sending field", ob_flds2D(fld2D,ob_id)
             print*,"Sending with tag:",tag,"for ob_id",ob_id,
     &             "fld2D",fld2D, "at myIter", myIter, "proc",mpiMyId
             print *,"Sending array"
             print *,buffer2D
            endif

            call MPI_SEND(buffer2D, sNx+sNy, MPI_REAL8, 0, tag,
     &                    MPI_COMM_MODEL, ierror)

           ENDDO


           DO fld3D=1, ob_nFlds3D(ob_id)
            tag = ob_tags(ob_id, 2, fld3D)
            call OB_PASS_RL_to_R8(buffer2D,buffer3D,ob_id,3,fld3D)

            if (send_debug .eq. 1) then
             print *,"Sending field", ob_flds3D(fld3D,ob_id), "for mask"
     &              , ob_id
             print*,"Sending with tag:",tag,"for ob_id",ob_id,
     &             "fld3D",fld3D, "at myIter", myIter, "proc",mpiMyId
             print *,"Sending array shape: ", shape(buffer3D)
            endif

            call MPI_SEND(buffer3D, ob_levels3D(fld3D,ob_id)*(sNx+sNy),
     &                    MPI_REAL8, 0, tag, MPI_COMM_MODEL, ierror)

            if (send_debug .eq. 1) then
             DO k=1, ob_levels3D(fld3D,ob_id)
              print *,"printing at level k",k, "with tag", tag
              print*,buffer3D(:,k)
             ENDDO
            endif

           ENDDO
          end if
         end if

#endif


C-------------------------------------------------------------------------------
C           Master Process Tasks:
C           1. Append Local Arrays in process 0
C           2. Receive time-averaged field values on open boundary points from
C              child processes > 0 through MPI and append to final array.
C           3. Write final output to binary file
C-------------------------------------------------------------------------------
         IF (mpiMyId .eq. 0) then


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

C      Set averages back to zero
       CALL SET_TO_INI_STATE(myThid)

C      Incrementing number of time levels passed
       time_level = time_level + 1

      ENDIF
C      Move on to next time step
       nTimeSteps_ob = nTimeSteps_ob + 1
       time_passed = time_passed + deltaTClock

      RETURN
      END




C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



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
      INTEGER, PARAMETER :: debug = 1
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

      print *,"ENTERED MASTER PROCESS"
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

        call OB_PASS_RL_to_R8(buffer2D, buffer3D, ob_id, 2,fld)
        call APPEND_OB_VALUE2D(ob_id, 0, buffer2D)
        total_valid = total_valid + numOBPnts_allproc(ob_id, 1)

        if (debug .eq. 1) then
         print *,"MASTER 2D: mask number:", ob_id, ",field",
     &            ob_flds2D(fld,ob_id)
         print *,"MASTER 2D: appended local field to global array"
     &          ," with ",numOBPnts_allproc(ob_id, 1),"valid elements:"
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
          print *,"MASTER: global_ob3D"
          print *,"shape(global_ob3D)",shape(global_ob3D)
          DO k=1, ob_levels3D(fld,ob_id)
           print *,"printing at level k", k
           print*,global_ob3D(:,k)
         ENDDO

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
       call PLOT_TEST_GLOBAL_OB(myThid, ob_id, fld, nDim, myIter)

       if (debug .eq. 1) then
        print*,"MASTER: myIter:",myIter
       end if

       print *,"WRITE TO FILE FOR MASK",ob_id
       call WRITE_GLOBAL_BIN(total_valid,fld,ob_id,myIter,nDim)


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
         integer num_valid, i, ob_index, ob_id, pid, glo
         integer, parameter :: print_debug = 1

         ob_index = 1
         num_valid = numOBPnts_allproc(ob_id,pid+1)

         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE2D: num_valid", num_valid
         end if

         do i=1, num_valid
           ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))
C           if (ob_index .le. (sNy+sNx)*(nPx*nPy)) then
            glo = sub_glo_indices_allproc(ob_id,pid+1,i)
             if (print_debug .eq. 1) then
               print *,"APPEND_OB_VALUE2D: glo", glo
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
C       Appends field value corresponding to open boundary points to final global array global_ob.
C       Use values of global_indices to extract value in lookup_table at position of value of global_indices.
C       Value in lookup_table used to identify position of current point on open boundary.
C       Store field value at position in final global array global_ob.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

        REAL*8 ob_values(sNx+sNy, Nr)
        integer num_valid, i, k, ob_index, ob_id, pid, fld, glo
        integer, parameter :: print_debug = 1

        ob_index = 1
        num_valid = numOBPnts_allproc(ob_id,pid+1)

        if (print_debug .eq. 1) then
          print *,"APPEND_OB_VALUE3D: mask", ob_id
          print *,"APPEND_OB_VALUE3D: field", ob_flds3D(fld,ob_id)
          print *,"APPEND_OB_VALUE3D: pid",pid
          print *,"APPEND_OB_VALUE3D: subglo",
     &      sub_glo_indices_allproc(ob_id,pid+1,:)

        endif

        do i=1, num_valid
         ob_index = lookup_table(ob_id,
     &                 sub_glo_indices_allproc(ob_id,pid+1,i))

         glo = sub_glo_indices_allproc(ob_id,pid+1,i)

         if (print_debug .eq. 1) then
           print *,"APPEND_OB_VALUE3D: glo", glo
           print *, "APPEND_OB_VALUE3D: ob_index", ob_index
           print *, "APPEND_OB_VALUE3D: counter i",i
           print*,"ob_values at k=1",ob_values(i,1)
         endif

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
          print*,"For fld",fld," and mask",ob_id
          print *, "APPEND_OB_VALUE: global_ob3D with total depth:",
     &       ob_levels3D(fld, ob_id), "printing for depth #",k
          print*,"APPEND_OB_VALUE ob_values:",ob_values(1:num_valid,k)
          print *,"APPEND_OB_VALUE global_ob3D",
     &      global_ob3D(1:num_valid,k)
         enddo
         print*,"FINAL: pid",pid,"obid:",ob_id
         print*,"FINAL 3D glo arr",
     &           sub_glo_indices_allproc(ob_id,pid+1,:)
       end if

       END


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

       REAL*8 global_NxNy(Nx,Ny,nSx,nSy)
       _RS dummyRS(1)
       REAL*8 global_theta(Nx,Ny)

       integer i,j, bi, bj, k, counter, ob_counter, correct
       integer fld, total_valid, ob_id, myThid, nDim, myIter
       integer, parameter :: plot_on_NxNy = 1

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
C                 !Inserting field points into NxNy array
                  if (nDim .eq. 2) then
                   global_NxNy(i,j,bi,bj) = global_ob2D(ob_counter)
                  else if  (nDim .eq. 3) then
                   global_NxNy(i,j,bi,bj) =
     &             global_ob3D(ob_counter,(ob_levels3D(fld,ob_id)/2)+1)
                  endif
                  ob_counter = ob_counter + 1
                 end if
               end if
               counter = counter + 1
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       if (plot_on_NxNy .eq. 1) then
         print *, "PLOT_TEST_GLOBAL_OB: plotting global_ob ",
     &              "on NxNy dimension array:"
         write(plotTitle, "(2A)") "PLOT_TEST_GLOBAL_OB: ", fldNm
         call PLOT_GLO_FIELD_XYRL(global_NxNy,plotTitle,myThid)

         write(plotTitle, "(A,I1)") "PLOT_TEST_GLOBAL_OB mask ",ob_id
         CALL PLOT_GLO_FIELD_XYRL(global_ob_mask(ob_id,:,:,:,:),
     &                    plotTitle, myThid)
       end if


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
