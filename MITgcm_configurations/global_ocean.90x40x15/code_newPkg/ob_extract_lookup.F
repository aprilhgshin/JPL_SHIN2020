C All variables in header file are global to header files included below
#include "DIAG_OPTIONS.h"





       SUBROUTINE OB_EXTRACT_LOOKUP( myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
C#ifdef ALLOW_FIZHI
#include "PARAMS.h"
C#endif
#include "DIAGNOSTICS_SIZE.h"
C#include "DIAGNOSTICS.h"
#include "DIAGNOSTICS_OB.h"
C No need to include "mpif.h" becaues EESUPPORT.h already includes it.

C Putting everything under do loop ob_id (EVERYTHING)


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678

       integer myThid

C      !LOCAL VARIABLES:
C      global_mask: 4 dimensional _RL array that combines all open boundarys
C      msgBuf: Used to write statements into the output file
C      ob_subMask_size: Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      global_ob_index: Integer indicating bottom-left global element number
C      stop_program: Integer used to hold either 1 or 0, 1 to stop program or 0 to ignore
C      debug: Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore

C      Although we read REAL*4 values, most built in functions have _RS input which by default is _RL.
C      So global_mask, ob_subMask, and bath_subMask are set as _RL or _RS.
       _RL global_mask(Nx,Ny,nSx,nSy)
       _RL recv_arr(sNx + sNy + 1)
       _RL g_indices(nPx*nPy, sNx + sNy + 1)
C       _RL timeave_cum(Nr, nSx, nSy)
       _RL subField(nOB_fld,1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)

       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, k, fld, iG, jG, bi, bj, step_ob, incr
       integer ob, global_ob_index, ob_index, new_goi
       integer ob_subMask_size, counter, ob_id
       integer total_valid, tag, time, num_valid, myIter
       integer myByLo_myThid,myByHi_myThid,myBxLo_myThid,myBxHi_myThid
       real avg
       integer, parameter :: stop_program = 0
       integer, parameter :: debug = 1
       integer, parameter :: print_debug = 0
       character*10 fldNm
       logical DIFFERENT_MULTIPLE

C     == Routine arguments Function DIFFERENT_MULTIPLE ==
C     freq       :: Frequency by which time is divided.
C     val1       :: time that is checked
C     step       :: length of time interval (around val1) that is checked
C      _RL  freq, val1, step
      _RL freq(1)
      _RL val1(1)
      _RL step(1)


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

        DO i=1, nOB_fld
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO k=1-OLx,sNx+OLx
                  subField(i,k,j,bi,bj) = 0.0
                ENDDO
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

C        DO i=1,Nr
C          DO j=1,nSx
C            DO k=1,nSy
C              timeave_cum(Nr, nSx, nSy) = deltaT_ob
C            ENDDO
C          ENDDO
C        ENDDO

c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy

C       myThid  ::  my Thread Id number

C     From PARAMS.h:
C     nIter0              :: Start time-step number of for this run
C     startTime


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

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      All Processes :



       DO ob_id = 1, nOB_mask

         step_ob = 0
         myIter = nIter0
         print *, "OB: deltaTMom, ", deltaTMom
         DO incr=0, totPhase_ob, int(deltaTMom)
C     freq       :: Frequency by which time is divided.
C     val1       :: time that is checked
C     step       :: length of time interval (around val1) that is checked
           freq(1) = avgPeriod_ob
           val1(1) = step_ob
           step(1) = deltaT_ob

           CALL INI_OB_FILENAMES(ob_filenames)

C          ob_filenames(1) = 'domain_flt32_mask1.bin'
C          ob_fildFiles(1) = 'bathymetry.bin'
C          ob_fildFiles(2) = 'lev_s.bin'
C          ob_fldNames(1) = 'bath'
C          ob_fldNames(2) = 'salt'


C          Loop through chosen fields to read from corresponding binary files
           DO fld=1, nOB_fld
C             CALL READ_PICKUP(incr, myThid)

C            Read field file. Passing in time step # as argument.
             print *,"time step #:",myIter
             CALL READ_REC_XY_RL( ob_fldFiles(fld),
     &           subField(fld,:,:,:,:),1, myIter, myThid)

C            Fill in halo regions
             _EXCH_XY_RL(subField(fld,:,:,:,:), myThid)

C            Write subField to binary file
             CALL WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy,
     &                  subField(fld,:,:,:,:), ob_fldNames(fld), ob_id)

             if (debug .eq. 1) then
               CALL PLOT_FIELD_XYRL(subField(fld,:,:,:,:),
     &                            ob_fldNames(fld), -1, myThid )
             end if

             DO bj = myByLo(myThid), myByHi(myThid)
               DO bi = myBxLo(myThid), myBxHi(myThid)
                 DO j=1-OLy,sNy+OLy
                   DO k=1-OLx,sNx+OLx
                     subField_avg(fld,k,j,bi,bj) =
     &           subField_avg(fld,k,j,bi,bj) + subField(fld,k,j,bi,bj)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

C          If we incremented to/near the end of an averaging period:
           if ((step_ob .eq. avgPeriod_ob) .or.
     &          (DIFFERENT_MULTIPLE(freq, val1, step) .and.
     &          (step_ob .le. avgPeriod_ob))) then


C             CALL DIAGNOSTICS_OB_WRITE_PICKUP(.FALSE.,
C     &         'flds',startTime+incr, incr/deltaTMom, myThid)


C            Initializing counting integers
             total_valid = 0
             num_valid = 0
             step_ob = 0

C            Computing field averages
             DO fld=1, nOB_fld
               DO bj = myByLo(myThid), myByHi(myThid)
                 DO bi = myBxLo(myThid), myBxHi(myThid)
                   DO j=1-OLy,sNy+OLy
                     DO k=1-OLx,sNx+OLx
                       subField_avg(fld,k,j,bi,bj) =
     &                    subField_avg(fld,k,j,bi,bj) / avgPeriod_ob
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

C            Reading from binary file to retrieve open boundary sub mask
             CALL READ_REC_XY_RL( ob_filenames(ob_id),
     &               ob_subMask,1, 0, myThid )
             print *, "OB: READ SUCCESSFUL"
             
             _EXCH_XY_RL(ob_subMask, myThid)

             CALL WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, ob_subMask,
     &                          "mask", ob_id)

C             DO bj = myByLo(myThid), myByHi(myThid)
C              DO bi = myBxLo(myThid), myBxHi(myThid)
C               CALL TIMEAVE_NORMALIZ(ob_subMask, timeave_cum,1,bi,bj,
C     &                                myThid)
C              ENDDO
C            ENDDO

             if (debug .eq. 1) then
C              Plotting for debugging
               CALL PLOT_FIELD_XYRL( ob_subMask,
     &                        'Open Boundary Submask', -1, myThid )
             end if


             ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy

             iG = mpi_myXGlobalLo(mpiMyId+1)
             jG = mpi_myYGlobalLo(mpiMyId+1)
             global_ob_index = Nx*(jG -1)  + iG

             counter = 2
             new_goi = global_ob_index

C            Computing global indices, storing values on open boundary points
             DO bj = myByLo(myThid), myByHi(myThid)
               DO bi = myBxLo(myThid), myBxHi(myThid)
                 DO j=1-OLy,sNy+OLy
                   DO i=1-OLx,sNx+OLx
                     if ((j .ge. 1) .and. (i .ge. 1)) then
                         new_goi = (j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
     &                              j*(iG-1)+i
                     end if
                     if (ob_subMask(i,j,bi,bj) .ge. 1.0 ) then
                       if (print_debug .eq. 1) then
                        print *, "OB_EXTRACT_LOOKUP: ob_subMask gt 1:",
     &                    ob_subMask(i,j,bi,bj)
                        print *,"OB_EXTRACT_LOOKUP: i j of ob_subMask",
     &                    " gt 1", i, j
                       end if

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

C           Inserting number of valid elements as first element of each array
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
                if (print_debug .eq. 1) then
                  print *,"OB_EXTRACT_LOOKUP: l sub_global_indices ",
     &                    "count:", int(sub_global_indices(1))
                  print *,"OB_EXTRACT_LOOKUP: local sub_global_indices"
                  call PRINT_FLOAT_ARR( sub_global_indices(:),num_valid)
                end if
              end if

            end if
#endif

C           Combining all sub masks to create one global mask
            call GATHER_2D_R8(global_mask(1:Nx,1:Ny,1,1),
     &                     ob_subMask(1:sNx,1:sNy,1,1),
     &                     Nx, Ny, .true., .true., myThid)

            if (mpiMYId .eq. 0) then
              call PLOT_GLO_FIELD_XYRL(global_mask(1:,1:,1,1),
     &                    "OB_EXTRACT_LOOKUP: global mask", myThid)
            end if


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
     &                                subFieldOnMask(fld,:), global_ob,
     &                                lookup_table, ob_id)

                      if (debug .eq. 1) then
                        print *,"APPEND_OB_VALUE: appending local "
     &                    ,"bath to global_ob with ",
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
     &                 recv_arr,global_ob, lookup_table, ob_id)

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

                  call WRITE_GLOBAL_BIN(global_ob,total_valid,
     &                          ob_fldNames(fld), ob_id, mpiPx, mpiPy)

                  call PLOT_TEST_GLOBAL_OB(global_ob, lookup_table,
     &                                myThid, ob_id, ob_fldNames(fld))

                ENDDO
              end if
            end if
#endif
          endif
          step_ob = step_ob + int(deltaTMom)
          myIter = myIter + 1
        ENDDO

      ENDDO


      END



C      Writing ob_subMask, bath_subMask, temp_subMask to bin file
C         call WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy,ob_subMask,
C     &                          bath_subMask, temp_subMask)


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


C      Using GATHER_2D_R8 to correctly combine ob_subMasks into a single global array

C     !INPUT/OUTPUT PARAMETERS FOR GATHER_2D_R4:
C gloBuff   ( _R4 ) :: full-domain 2D IO-buffer array             (Output)
C myField   ( _R4 ) :: tiled, local (i.e. my Proc. tiles) 2D array (Input)
C xSize    (integer):: global buffer 1rst dim (x)
C ySize    (integer):: global buffer 2nd  dim (y)
C useExch2GlobLayOut:: =T: Use Exch2 global-map layout (only with EXCH2)
C zeroBuff (logical):: =T: initialise the buffer to zero before copy
C myThid   (integer):: my Thread Id number

C      SUBROUTINE AVERAGE (avg, arr)
C      include "USER_INPUT.h"

C      _RL arr(nfld_time)
C      integer i, sum
C      real avg
C      sum = 0
C      do i=1, nfld_time
C        sum = sum + arr(i)
C      end do
C      avg = real(sum/nfld_time)

C      END


      SUBROUTINE INI_OB_FILENAMES(ob_filenames)
        include "USER_INPUT.h"
        CHARACTER*40 ob_filenames(nOB_mask)
        INTEGER k
C       Fill in array with filenames for each ob mask

        DO k=1, nOB_mask
          write(ob_filenames(k), "(A17,I1,A4)") "domain_flt32_mask",
     &          k, '.bin'
        ENDDO
      END


      SUBROUTINE APPEND_OB_VALUE(global_indices,ob_values,global_ob,
     &                           lookup_table, ob_id)
         include "SIZE.h"
         include "USER_INPUT.h"
         _RL global_indices(sNx + sNy + 1)
         _RL ob_values(sNx + sNy + 1)
         _RL global_ob((sNy+sNx)*(nPx*nPy))
         INTEGER lookup_table(nOB_mask, Ny*Nx)
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
     &                            is_child, mpiPx, mpiPy)
c       Creating filenames for all binary files.

       integer is_child, mpiPx, mpiPy, ob_id
       character*30 filename
       character*4 fldNm
       character*2 dir
       character*3 str_px, str_py
       character str_ob

c       From EESUPPORT.h
C       mpiPx     - My MPI proc. grid X coord
C       mpiPy     - My MPI proc. grid Y coord

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       001 for index 1 in x axis
c       002 for index 2 in x axis
       write(str_ob, "(I1)") ob_id
       if (is_child .eq. 1) then
         write(str_px, "(I3.3)") mpiPx
         write(str_py, "(I3.3)") mpiPy
         write(filename, "(10A)") dir, "ob",str_ob,".",fldNm,".",
     &                            str_px, ".", str_py, ".bin"
       else
c         Creating a filename for global arrays
         write(filename, "(5A)") dir, fldNm, "_on_mask", str_ob,
     &                           "_global.bin"
       end if

       end


       SUBROUTINE WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, subMask,
     &                          fldNm, ob_id)
C      Writing all sub arrays read from binary fie containing global array
       include "SIZE.h"

       integer mpiMyId, mpiPx, mpiPy, ob_id
       character*4 fldNm
       character*2 dir
       character*30 filename
       _RL subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
       integer, parameter :: print_write = 1

       dir = "./"
       call CREATE_FILENAME(dir, fldNm, ob_id, filename, 1,
     &                      mpiPx, mpiPy)
       open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(1) subMask
       close(1)
       if (print_write .eq. 1) then
         print *, fldNm, " for ob",ob_id,"process", mpiMyId,
     &            "is written to ",filename
       end if

       end


       SUBROUTINE WRITE_GLOBAL_BIN(global_ob,total_valid,fldNm,ob_id,
     &                             mpiPx, mpiPy)
       include "SIZE.h"

       _RL global_ob((sNy+sNx)*(nPx*nPy))
       integer total_valid, ob_id, mpiPx, mpiPy

       character*2 dir
       character*4 fldNm
       character*30 filename

c      Writing global array to binary file
       dir = "./"
       call CREATE_FILENAME(dir, fldNm, ob_id, filename, 0,
     &                      mpiPx, mpiPy)
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
C                Inserting bathymetry points into NxNy array
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
          write(*,"(F6.3,$)") array(i)
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
