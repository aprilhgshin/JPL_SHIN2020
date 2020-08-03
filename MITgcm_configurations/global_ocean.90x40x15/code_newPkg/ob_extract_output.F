C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678
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
C      boundary points in their assigned sub-mask and sends an array of field
C      values and indices with respect to the global domain corresponding to
C      each of the open boundary points back to the parent process. Parent
C      process combines thessse arrays and writes result to binary files. 
C
C      input parameters: myTime, myIter, myThid
C
C-------------------------------------------------------------------------------

#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

       _RL myTime
C     myThid  ::  my Thread Id number
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

C      global_mask: 4 dimensional _RL array that combines all open boundarys
C      recv_arr: buffer array to store incoming arrays from MPI_RECV
C      g_indices: array that holds global indices of sub fields from each child process

C      msgBuf: Used to write statements into the output file
C      global_ob_index: Integer indicating bottom-left global element number
C      ob_index: Integer obtained from lookup table in corresponding global index
C      new_goi: Integer value holding computed global index
C      ob_subMask_size: Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      ob_id: Integer used for do loop through each ob mask
C      total_valid: Integer value of total count of valid elements in each open boundary mask/ Count done when looping through arrays sent from processes + local field on mask array
C      num_valid: Integer value of count of valid elements in each open boundary mask. Count done when looping through mask
C      tag: integer value for looping through all message tags

C      debug: Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C      print_debug: Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore
C-------------------------------------------------------------------------------


C      Although we read REAL*4 values, most built in functions have _RS input which by default is _RL.
C      So global_mask, ob_subMask, and bath_subMask are set as _RL or _RS.
       _RL global_mask(Nx,Ny,nSx,nSy)
       _RL recv_arr(sNx + sNy + 1)
       _RL g_indices(nPx*nPy, sNx + sNy + 1)

       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, k, fld, iG, jG, bi, bj, step_ob, incr
       integer global_ob_index, ob_index, new_goi
       integer ob_subMask_size, counter, ob_id
       integer total_valid, num_valid, tag
       integer myByLo_myThid,myByHi_myThid,myBxLo_myThid,myBxHi_myThid

       integer, parameter :: debug = 1
       integer, parameter :: print_debug = 0

       character*10 fldNm
       logical DIFFERENT_MULTIPLE


#ifdef ALLOW_USE_MPI
       integer pid, ierror, status(MPI_STATUS_SIZE)
#endif

       common / local_var / global_mask, total_valid
       common / myBxy_LoHi / myByLo_myThid, myByHi_myThid,
     &                       myBxLo_myThid, myBxHi_myThid

       myByLo_myThid = myByLo(myThid)
       myByHi_myThid = myByHi(myThid)
       myBxLo_myThid = myBxLo(myThid)
       myBxHi_myThid = myBxHi(myThid)


C       Initialize local arrays:
        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1, Ny
              DO i=1, Nx
                global_mask(i,j,bi,bj) = 0.0
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO i=1, sNx + sNy + 1
          recv_arr(i) = 0.0
        ENDDO


        DO pid=1, mpiNprocs
          DO i=1, sNx + sNy + 1
            g_indices(pid, i)=0.0
          ENDDO
        ENDDO

C-------------------------------------------------------------------------------
c      All Processes Tasks:
C-------------------------------------------------------------------------------

       IF (myTime .eq. startTime) then
         call INI_OB_FILENAMES()
C        Reading from binary file to retrieve open boundary sub mask
         DO i=1, nOB_mask
           CALL READ_REC_XY_RL( ob_filenames(i),
     &               ob_subMask(i,:,:,:,:),1, 0, myThid )
           _EXCH_XY_RL(ob_subMask(i,:,:,:,:), myThid)
           IF (debug .eq. 1) then
C            Plotting for debugging
             CALL PLOT_FIELD_XYRL( ob_subMask(i,:,:,:,:),
     &                        'Open Boundary Submask', -1, myThid )
           ENDIF
         ENDDO
       ENDIF

c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
       ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy

        CALL SET_SUBFIELDS(myThid)

        IF (print_debug .eq. 1) then
         CALL PRINT_FLOAT_ARR(subField(1,:,:,:,:),ob_subMask_size)
         CALL PRINT_FLOAT_ARR(salt(:,:,1,:,:),ob_subMask_size)
        ENDIF

        IF (debug .eq. 1) then
         CALL  PLOT_FIELD_XYRL(subField(1,:,:,:,:),
     &       'OB_EXTRACT_OUTPUT: salt copy in subField', -1, myThid )
         CALL  PLOT_FIELD_XYRL(salt(:,:,1,:,:),
     &            'OB_EXTRACT_OUTPUT: originalsalt', -1, myThid )
       ENDIF

C       deltaTClock = 43200.0
C       endTime = endTime + totPhase_ob
C        startTime = 0.0


       IF (myIter .eq. nIter0) then
C          Initialize averages to zero and number of time steps to one
         CALL SET_TO_INI_STATE(myThid)
       ELSE
C       Cumulate Diagnostics_ob fields and number of time steps (for Time Average)
         CALL CUMULATE_FLD_TIME(myThid)
       ENDIF

       if (debug .eq. 1) then
         print *, "OB_EXTRACT_OUTPUT timepassed:", time_passed
         print *, "OB_EXTRACT_OUTPUT deltaTClock:", deltaTClock
         print *, "OB_EXTRACT_OUTPUT avgPeriod_ob:", avgPeriod_ob
         print *, "OB_EXTRACT_OUTPUT myTime:", myTime
         print *, "OB_EXTRACT_OUTPUT startTime:",startTime
         print *, "OB_EXTRACT_OUTPUT model endTIme", endTime
         print *, "OB_EXTRACT_OUTPUT user end time:",
     & startTime + totPhase_ob
         print *, "OB_EXTRACT_OUTPUT DIFFERENT_MULTIPLE",
     &          DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock)
      end if



       IF (((myIter.NE.nIter0)
     &       .AND.DIFFERENT_MULTIPLE(avgPeriod_ob,myTime,deltaTClock))
     &       .AND. (time_passed .le. totPhase_ob)) THEN

C          Compute time average of all sub fields
         CALL TIMEAVE_OB_FLD(myThid, myTime)

         DO ob_id = 1, nOB_mask

C          Initializing counting integers
           num_valid = 0
           step_ob = 0

           iG = mpi_myXGlobalLo(mpiMyId+1)
           jG = mpi_myYGlobalLo(mpiMyId+1)
           global_ob_index = Nx*(jG -1)  + iG

           counter = 2
           new_goi = global_ob_index
C          Computing global indices, storing values on open boundary points
           DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO i=1-OLx,sNx+OLx
                  if ((j .ge. 1) .and. (i .ge. 1)) then
                      new_goi = (j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
     &                              j*(iG-1)+i
                  end if
                  if (ob_subMask(ob_id,i,j,bi,bj) .ge. 1.0 ) then

                    DO fld=1, nOB_fld
                     subFieldOnMask(fld,counter) =
     &                                     subField_avg(fld,i,j,bi,bj)
                    ENDDO
                    sub_global_indices(counter) = new_goi
                    counter = counter + 1
                    if (((i .ge. 1) .and. (j .ge. 1))
     &                  .and. ((i .le. sNx) .and. (j .le. sNy))) then
                      num_valid = num_valid + 1
                    end if
                  end if

                ENDDO
              ENDDO
            ENDDO
           ENDDO

           if (print_debug .eq. 1) then
             print *, "OB_EXTRACT_OUTPUT: num_valid:",num_valid
           end if

C          Inserting number of valid elements as first element of each array
           DO fld=1,nOB_fld
             subFieldOnMask(fld,1) = num_valid
           ENDDO
           sub_global_indices(1) = num_valid


#ifdef ALLOW_USE_MPI
            if (usingMPI) then
              if (mpiMyId .gt. 0) then

C               Processes 1 and onward send sub arrays to Parent process
                call MPI_SEND(sub_global_indices,(sNx+sNy+1),
     &               MPI_REAL, 0, 1, MPI_COMM_MODEL, ierror)
                do fld=1, nOB_fld
                  call MPI_SEND(subFieldOnMask(fld,:), sNx+sNy+1,
     &               MPI_REAL,0, fld+1, MPI_COMM_MODEL, ierror)
                end do

              else if (mpiMYId .eq. 0) then
C               Process 0 does not need to send arrays to itself.
                if (debug .eq. 1) then
                  print *,"OB_EXTRACT_LOOKUP: l sub_global_indices ",
     &                    "count:", int(sub_global_indices(1))
                  print *,"OB_EXTRACT_LOOKUP: local sub_global_indices"
                  call PRINT_FLOAT_ARR( sub_global_indices,num_valid)
                end if
              end if

            end if
#endif


C      Using GATHER_2D_R8 to correctly combine ob_subMasks into a single global array

C     !INPUT/OUTPUT PARAMETERS FOR GATHER_2D_R4:
C gloBuff   ( _R4 ) :: full-domain 2D IO-buffer array             (Output)
C myField   ( _R4 ) :: tiled, local (i.e. my Proc. tiles) 2D array (Input)
C xSize    (integer):: global buffer 1rst dim (x)
C ySize    (integer):: global buffer 2nd  dim (y)
C useExch2GlobLayOut:: =T: Use Exch2 global-map layout (only with EXCH2)
C zeroBuff (logical):: =T: initialise the buffer to zero before copy
C myThid   (integer):: my Thread Id number

C           Combining all sub masks to create one global mask
            call GATHER_2D_R8(global_mask(1:Nx,1:Ny,1,1),
     &                     ob_subMask(ob_id,1:sNx,1:sNy,1,1),
     &                     Nx, Ny, .true., .true., myThid)

            if (mpiMYId .eq. 0) then
              call PLOT_GLO_FIELD_XYRL(global_mask(1:,1:,1,1),
     &                    "OB_EXTRACT_LOOKUP: global mask", myThid)
            end if


C-------------------------------------------------------------------------------
C           Parent Process Tasks:
C-------------------------------------------------------------------------------
#ifdef ALLOW_USE_MPI
            if (usingMPI) then
              if (mpiMyId .eq. 0) then

C               Store global_ob_masks into lookup table
                counter = 1
                do j=1, Ny, 1
                  do i=1, Nx, 1
                    lookup_table(ob_id,counter) =
     &                               int(global_mask(i,j,1,1))
                    counter = counter + 1
                  end do
                end do

                ob_index = 1
C               Retrieving sub_global_indices first from all processes > 0
                do pid=1, mpiNprocs-1
                  call MPI_RECV(recv_arr, sNx+sNy+1, MPI_REAL,
     &                   pid, 1, MPI_COMM_MODEL, status, ierror)
C                 Save all global indices into array for future reference
                  do i=1, sNx+sNy+1
                    g_indices(pid,i) = recv_arr(i)
                  end do
                end do

C               Appending subFields to final global_ob array
C               Field numbers, fld, are used as message tags.
C               DO tag=2, nOB_fld + 1
                DO fld=1, nOB_fld

                  total_valid = 0
                  DO pid=1, mpiNprocs-1

C                 First append process 0 local arrays to final global array if number of valid elements in those local arrays > 0:

                    if (int(sub_global_indices(1)) .ge. 1) then

                      if (print_debug .eq. 1) then
                       print *,"OB_EXTRACT_LOOKUP: Appending local"
                       call PRINT_FLOAT_ARR(subFieldOnMask(1,:),
     &                                 int(sub_global_indices(1) + 1))
                      end if
                      total_valid = total_valid +
     &                                  int(sub_global_indices(1))
                      call APPEND_OB_VALUE(sub_global_indices(:),
     &                                subFieldOnMask(fld,:), ob_id)

                      if (debug .eq. 1) then
                        print *,"APPEND_OB_VALUE: appending local "
     &                    ,"field to global_ob with ",
     &                    int(sub_global_indices(1)),"valid elements:"
                        call PRINT_FLOAT_ARR(global_ob, total_valid)
                      end if

                    end if

C                   Receive all sub arrays from processes with pid > 0.
                    call MPI_RECV(recv_arr, (sNx+sNy+1),MPI_REAL,
     &                   pid, fld+1 ,MPI_COMM_MODEL, status, ierror)
                    if (debug .eq. 1) then
                      print *, "Received recv_arr for fields"
                      print *, "shape:", shape(recv_arr)
                      print *, "Displaying FULL array, beyond number"
     &                       ," of valid elements"
                      call PRINT_FLOAT_ARR(recv_arr, sNx+sNy+1)
                    end if

C                   Appending all other processes' local arrays to final array
C                                   if number of valid elements in arrays > 0:
                    if (int(g_indices(pid, 1)) .ge. 1) then
                      print *,"OB: recv fields:"
                      print *,"OB: recv shape:", shape(recv_arr)
                      call PRINT_FLOAT_ARR(recv_arr,sNx+sNy+1)
                     call APPEND_OB_VALUE(g_indices(pid,:),
     &                 recv_arr, ob_id)

                     total_valid = total_valid + g_indices(pid,1)

                    end if

                    print *,"OB_EXTRACT_LOOKUP: updated global ob by "
     &                  ,"process",pid,"with",int(g_indices(pid,1)),
     &                   "valid elements"
                    call PRINT_FLOAT_ARR(global_ob,total_valid)
                  ENDDO

                  print *, "OB_EXTRACT_LOOKUP: final global ob by ",
     &     "for tag", fld+1, "with total", total_valid, "valid elements"
                  call PRINT_FLOAT_ARR(global_ob,total_valid)

                  call WRITE_GLOBAL_BIN(total_valid, fld,
     &                           ob_id, myTime)

                  call PLOT_TEST_GLOBAL_OB(global_ob, lookup_table,
     &                                myThid, ob_id, ob_fldNames(fld))

                ENDDO
              end if
            end if
#endif
        ENDDO
       ENDIF

       time_passed = time_passed + deltaTClock

      END


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



      SUBROUTINE INI_OB_FILENAMES()
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"
        INTEGER k
C       Fill in array with filenames for each ob mask
        DO k=1, nOB_mask
          write(ob_filenames(k), "(A17,I1,A4)") "domain_flt32_mask",
     &          k, '.bin'
        ENDDO
      END


      SUBROUTINE SET_TO_INI_STATE(myThid)
C     Initializes average fields and number of time steps to zero
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER j, k, bi, bj, fld

        nTimeSteps_ob = 0
        DO fld=1, nOB_fld
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO k=1-OLx,sNx+OLx
                  subField_avg(fld,k,j,bi,bj) = 0.0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      END

      SUBROUTINE SET_SUBFIELDS(myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER j, k, bi, bj, fld, opt
        print *,"Inside SET_SUBFIELDS"
        print *,"SET_SUBFIELDS: shape fld_choice",shape(fld_choice)
        print *,"SET_SUBFIELDS: nFldOpt",nFldOpt
        print *,"SET_SUBFIELDS: fld_choice(nFldOpt)",fld_choice(nFldOpt)

        DO fld=1, nFldOpt
          print*, "SET_SUBFIELDS: fld_choice(opt)",fld_choice(fld)
        ENDDO

        DO fld=1, nOB_fld
         DO opt=1, nFldOpt
          IF (fld_choice(opt)) then
            DO bj = myByLo(myThid), myByHi(myThid)
              DO bi = myBxLo(myThid), myBxHi(myThid)
                DO j=1-OLy,sNy+OLy
                  DO k=1-OLx,sNx+OLx
                    IF (opt .eq. 1) then
                      subField(fld,k,j,bi,bj) = etaN(k,j,bi,bj)
                    ELSE IF (opt .eq. 2) then
                      subField(fld,k,j,bi,bj) = uVel(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 3) then
                      subField(fld,k,j,bi,bj) = vVel(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 4) then
                      subField(fld,k,j,bi,bj) = theta(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 5) then
                      subField(fld,k,j,bi,bj) = salt(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 6) then
                      subField(fld,k,j,bi,bj) = gU(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 7) then
                      subField(fld,k,j,bi,bj) = gV(k,j,1,bi,bj)
                    ELSE IF (opt .eq. 8) then
                      subField(fld,k,j,bi,bj) = etaH(k,j,bi,bj)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
       ENDDO
      END

      SUBROUTINE CUMULATE_FLD_TIME(myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER j, k, bi, bj, fld

        nTimeSteps_ob = nTimeSteps_ob + 1
        DO fld=1, nOB_fld
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO k=1-OLx,sNx+OLx
                  subField_avg(fld,k,j,bi,bj) =
     &           subField_avg(fld,k,j,bi,bj)+subField(fld,k,j,bi,bj)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END

      SUBROUTINE TIMEAVE_OB_FLD(myThid, myTime)
C       Cumulate Diagnostics_ob fields and number of time steps (for Time Average)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER myThid
        INTEGER j, k, bi, bj, fld
        _RL myTime

        DO fld=1, nOB_fld
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO k=1-OLx,sNx+OLx
                  subField_avg(fld,k,j,bi,bj) =
     &           subField_avg(fld,k,j,bi,bj)/nTimeSteps_ob
                ENDDO
              ENDDO
            ENDDO
          ENDDO
C         Write subField_avg to binary file
         CALL WRITE_SUB_BIN(subField_avg(fld,:,:,:,:), fld,
     &                      myTime)
        ENDDO


      END




      SUBROUTINE APPEND_OB_VALUE(global_indices,ob_values, ob_id)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

         _RL global_indices(sNx + sNy + 1)
         _RL ob_values(sNx + sNy + 1)
         integer num_valid, i, ob_index, ob_id
         integer, parameter :: print_debug = 0

         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: global_indices:"
           call PRINT_FLOAT_ARR(global_indices, sNx + sNy + 1)
         end if

         num_valid = int(ob_values(1))
         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: num_valid", num_valid
           print *, "APPEND_OB_VALUE: shape values:", shape(ob_values)
           print *, "APPEND_OB_VALUE: shape global_ob", shape(global_ob)
         end if

         do i=2, num_valid + 1
           ob_index = lookup_table(ob_id, int(global_indices(i)))
           if (ob_index .le. (sNy+sNx)*(nPx*nPy)) then
             if (print_debug .eq. 1) then
               print *, "APPEND_OB_VALUE: ob_index", ob_index
               print *, "APPEND_OB_VALUE: counter i",i
               print *, "APPEND_OB_VALUE: ob_value", ob_values(i)
             end if
             global_ob(ob_index) = ob_values(i)
           end if
         end do
         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: global_ob:"
           call PRINT_FLOAT_ARR(global_ob, (sNy+sNx)*(nPx*nPy))
         end if

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

c       From EESUPPORT.h
C       mpiPx     - My MPI proc. grid X coord
C       mpiPy     - My MPI proc. grid Y coord

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


       SUBROUTINE PLOT_TEST_GLOBAL_OB(global_ob, lookup_table, myThid,
     &                                ob_id, fldNm)
c      Unravels global_ob and plots onto an array with dimensions Nx by Ny
       include "SIZE.h"
       include "USER_INPUT.h"
       _RL global_ob((sNy+sNx)*(nPx*nPy))
       _RL global_mask(Nx,Ny,nSx,nSy)
       INTEGER lookup_table(nOB_mask, Ny*Nx)
       _RL global_correct((sNy+sNx)*(nPx*nPy))
       _RL global_NxNy(Nx,Ny,nSx,nSy)

       integer i,j, bi, bj, counter, ob_counter, correct, myThid
       integer myByLo_myThid, myByHi_myThid, total_valid, ob_id
       integer myBxLo_myThid, myBxHi_myThid
       integer, parameter :: plot_on_NxNy = 1
       character*4 fldNm
       character*30 plotTitle

       common / local_var / global_mask, total_valid
       common / myBxy_LoHi / myByLo_myThid, myByHi_myThid,
     &                       myBxLo_myThid, myBxHi_myThid

C      Initializing global_correct, used to indicate if elements of global_ob is correct or not
       DO i=1,(sNy+sNx)*(nPx*nPy)
         global_correct(i) = 1.0
       ENDDO

C       Initialize NxNy array:

        DO bj = myByLo_myThid, myByHi_myThid
          DO bi = myBxLo_myThid, myBxHi_myThid
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
       DO bj = myByLo_myThid, myByHi_myThid
         DO bi = myBxLo_myThid, myBxHi_myThid
           DO j=1, Ny
             DO i=1, Nx
               if (ob_counter .le. total_valid) then
                 if (lookup_table(ob_id,counter).eq.ob_counter) then
C                Inserting field points into NxNy array
                   global_NxNy(i,j,bi,bj) = global_ob(ob_counter)
                   if (global_mask(i,j,bi,bj) .eq. 0) then
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
         print *, "PLOT_TEST_GLOBAL_OB: global_ob final array is ",
     &            "incorrect for ", fldNm, "for ob", ob_id
         call PRINT_FLOAT_ARR(global_correct, (sNy+sNx)*(nPx*nPy))
       else
         print *, "PLOT_TEST_GLOBAL_OB: global_ob final array is ",
     &         "correct for the ", fldNm, " field for ob", ob_id
         if (plot_on_NxNy .eq. 1) then
           print *, "PLOT_TEST_GLOBAL_OB: plotting global_ob ",
     &              "on NxNy dimension array:"
           write(plotTitle, "(2A)") "PLOT_TEST_GLOBAL_OB: ", fldNm
           call PLOT_GLO_FIELD_XYRL(global_NxNy,plotTitle,myThid)
         end if
       end if
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
