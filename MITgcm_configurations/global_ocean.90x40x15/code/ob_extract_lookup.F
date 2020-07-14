C All variables in header file are global to header files included below
#include "DIAG_OPTIONS.h"


       SUBROUTINE OB_EXTRACT_LOOKUP( myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#ifdef ALLOW_FIZHI
#include "PARAMS.h"
#endif
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"
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
       _RL global_mask(nOB_mask, Nx,Ny,nSx,nSy)
       _RL recv_arr(nOB_mask, sNx + sNy + 1)
       _RL g_indices(nPx*nPy, nOB_mask, sNx + sNy + 1)
       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, iG, jG, bi, bj
       integer ob, global_ob_index, ob_index, new_goi
       integer ob_subMask_size, counter, ob_id
       integer total_valid, tag
       integer myByLo_myThid,myByHi_myThid,myBxLo_myThid,myBxHi_myThid
       integer, parameter :: stop_program = 0
       integer, parameter :: debug = 1
       integer, parameter :: print_debug = 0
       character*10 fldNm
       character*40 b_fileName
#ifdef ALLOW_USE_MPI
       integer pid, ierror, status(MPI_STATUS_SIZE)
#endif

       common / local_var / global_mask, total_valid, g_indices
       common / myBxy_LoHi / myByLo_myThid, myByHi_myThid,
     &                       myBxLo_myThid, myBxHi_myThid

       myByLo_myThid = myByLo(myThid)
       myByHi_myThid = myByHi(myThid)
       myBxLo_myThid = myBxLo(myThid)
       myBxHi_myThid = myBxHi(myThid)

C      Initializing before adding onto total_valid
       total_valid = 0

C       Initialize local arrays:
        DO ob_id=1, nOB_mask
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1, Ny
                DO i=1, Nx
                  global_mask(ob_id,i,j,bi,bj) = 0.0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO


        DO ob_id=1, nOB_mask
          DO i=1, sNx + sNy + 1
            recv_arr(ob_id, i) = 0.0
          ENDDO
        ENDDO


        DO pid=1, mpiNprocs
          DO ob_id=1, nOB_mask
            DO i=1, sNx + sNy + 1
              g_indices(pid,ob_id,i)=0.0
            ENDDO
          ENDDO
        ENDDO


c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy

C       myThid  ::  my Thread Id number

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

C        Differing filenames yet to be figured out.
         b_fileName = 'bathymetry.bin'
         t_fileName = 'domain_temp.bin'
         ob_fileName = 'domain_flt32_mask.bin'

         CALL READ_REC_XY_RL( t_fileName,
     &               temp_subMask(ob_id,:,:,:,:),1, 0, myThid )
         CALL READ_REC_XY_RL( b_fileName,
     &               bath_subMask(ob_id,:,:,:,:),1, 0, myThid )
         CALL READ_REC_XY_RL( ob_fileName,
     &               ob_subMask(ob_id,:,:,:,:),1, 0, myThid )

         _EXCH_XY_RL(temp_subMask(ob_id,:,:,:,:), myThid)
         _EXCH_XY_RL(bath_subMask(ob_id,:,:,:,:), myThid)
         _EXCH_XY_RL(ob_subMask(ob_id,:,:,:,:), myThid)

         CALL WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, temp_subMask,
     &                      "temp", ob_id)

         CALL WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, bath_subMask,
     &                      "bath", ob_id)

         CALL WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, temp_subMask,
     &                      "mask", ob_id)


         ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*(abs(1-OLy)+(sNy+OLy))
     &                    *nSx*nSy

         iG = mpi_myXGlobalLo(mpiMyId+1)
         jG = mpi_myYGlobalLo(mpiMyId+1)
         global_ob_index = Nx*(jG -1)  + iG

         if (print_debug .eq. 1) then
           print *, "process id: ", mpiMyId
           print *, "iG:", iG, "jG:", jG
           print *, "OB_EXTRACT_LOOKUP:global_ob_index",global_ob_index
         end if

         counter = 2
         new_goi = global_ob_index

        if (debug .eq. 1) then
C         Plotting for debugging
          CALL PLOT_FIELD_XYRL( ob_subMask(ob_id,:,:,:,:),
     &                      'Open Boundary Submask', -1, myThid )

          CALL PLOT_FIELD_XYRL( bath_subMask(ob_id,:,:,:,:),
     &                      'Bathymetry Submask', -1, myThid )

          CALL PLOT_FIELD_XYRL( temp_subMask(ob_id,:,:,:,:),
     &                      'Temperature Submask', -1, myThid )

        end if

         if (print_debug .eq .1) then
           print *,"OB_EXTRACT_LOOKUP: MPI total num procs:",mpiNprocs
           print *,"OB_EXTRACT_LOOKUP: myThid: ", myThid
         end if


         if (print_debug .eq .1) then
           print *,"OB_EXTRACT_LOOKUP: MPI total num procs:",mpiNprocs
           print *,"OB_EXTRACT_LOOKUP: myThid: ", myThid
         end if


C      In this program, we only have one open boundary. Therefore, first dimension
C      of subBathOnMask and sub_global_indices is 1, referring to open boundary # 1
         DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
             DO j=1-OLy,sNy+OLy
               DO i=1-OLx,sNx+OLx
                 if ((j .ge. 1) .and. (i .ge. 1)) then
                     new_goi = (j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
     &                        j*(iG-1)+i
                   if ((i .eq. 1) .and. (j .eq. 1)) then
                     print *, "at 1,1: ", new_goi
                   end if
                 end if
                 if (ob_subMask(ob_id,i,j,bi,bj) .gt. 0.0 ) then
                   if (print_debug .eq. 1) then
                     print *, "OB_EXTRACT_LOOKUP: ob_subMask gt 1:",
     &                    ob_subMask(ob_id,i,j,bi,bj)
                     print *, "OB_EXTRACT_LOOKUP: i j of ob_subMask",
     &                    " gt 1", i, j
                     print *, "OB_EXTRACT_LOOKUP: bath_subMask gt 1",
     &                    bath_subMask(ob_id,i,j,bi,bj)
                   end if
                   subBathOnMask(ob_id, counter) =
     &                                bath_subMask(ob_id,i,j,bi,bj)
                   subTempOnMask(ob_id, counter) =
     &                                temp_subMask(ob_id,i,j,bi,bj)
                   sub_global_indices(ob_id, counter) = new_goi
                   if (print_debug .eq. 1) then
                     print *, "OB_EXTRACT_LOOKUP: new_goi:", new_goi
                   end if
                   counter = counter + 1
                 end if

               ENDDO
             ENDDO
           ENDDO
         ENDDO

C       Inserting number of valid elements as first element of each array
        subBathOnMask(ob_id, 1) = counter - 2
        subTempOnMask(ob_id, 1) = counter - 2
        sub_global_indices(ob_id, 1) = counter - 2



       ENDDO


#ifdef ALLOW_USE_MPI
          if (usingMPI) then
            if (mpiMyId .gt. 0) then

C             Processes 1 and onward, send sub arrays to Parent process
              call MPI_SEND(sub_global_indices,nOB_mask*(sNx+sNy+1),
     &                MPI_REAL, 0, 1, MPI_COMM_WORLD, ierror)
              call MPI_SEND(subBathOnMask, nOB_mask*(sNx+sNy+1),
     &               MPI_REAL,0, 2, MPI_COMM_MODEL, ierror)
              call MPI_SEND(subTempOnMask, nOB_mask*(sNx+sNy+1),
     &               MPI_REAL,0, 3, MPI_COMM_MODEL, ierror)

            else if (mpiMYId .eq. 0) then
              print *,"OB_EXTRACT_LOOKUP: l sub_global_indices count:"
     &               ,  sub_global_indices(1,1)
              print *, "OB_EXTRACT_LOOKUP: local sub_global_indices"
              call PRINT_FLOAT_ARR( sub_global_indices(1,:),sNx+sNy+1)
            end if

          end if
#endif



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

       DO ob_id=1, nOB_mask
         call GATHER_2D_R8(global_mask(ob_id,1:Nx,1:Ny,1,1),
     &                     ob_subMask(ob_id,1:sNx,1:sNy,1,1),
     &                     Nx, Ny, .true., .true., myThid)

        if (mpiMYId .eq. 0) then
          call PLOT_GLO_FIELD_XYRL(global_mask(ob_id,1:,1:,1,1),
     &                      "OB_EXTRACT_LOOKUP: global mask", myThid)
        end if
       ENDDO



#ifdef ALLOW_USE_MPI
       if (usingMPI) then
         if (mpiMyId .eq. 0) then

C          Store global_ob_masks into lookup table
            DO ob_id=1, nOB_mask, 1
              counter = 1
              do j=1, Ny, 1
                do i=1, Nx, 1
                  lookup_table(ob_id, counter)=
     &                            int(global_mask(ob_id,i,j,1,1))
                  if (print_debug .eq. 1) then
                    if (lookup_table(ob_id, counter) .ge. 1) then
                      print *,"OB_EXTRACT_LOOKUP:global index positive:"
     &                       ,lookup_table(ob_id, counter),"at",counter
                    end if
                  end if
                  counter = counter + 1
                end do
              end do
            ENDDO


C          Retrive sub_global_indicies and sub fields from every process > 0.
            ob_index = 1


C           Retrieving sub_global_indices first from all processes > 0
            do pid=1, mpiNprocs-1
              call MPI_RECV(recv_arr, nOB_mask*(sNx+sNy+1)
     &                  ,MPI_INTEGER, pid, 1, MPI_COMM_WORLD, status,
     &                  ierror)
C             Save all global indices into array for future reference
              do ob_id=1, nOB_mask
                do i=1, sNx+sNy+1
                  g_indices(pid,ob_id,i) = recv_arr(ob_id,i)
                end do
              end do
            end do

            do ob_id=1, nOB_mask

C             Retrieving sub fields and appending to final global_ob array
              do tag=2, nOB_fld
                do pid=1, mpiNprocs-1

C                 First append process 0 local arrays to final global array if number of valid elements in those local arrays > 0:
                  if (print_debug .eq. 1) then
                    print *, "OB_EXTRACT_LOOKUP: l sub_global_indices ",
     &                     "after count: ",sub_global_indices(1,1)
                  end if
                  if (int(sub_global_indices(ob_id,1)) .ge. 1) then
                    if (print_debug .eq. 1) then
                      print *,"OB_EXTRACT_LOOKUP: Appending local"
                      print *, "OB_EXTRACT_LOOKUP: subBathOnMask local:"
                      call PRINT_FLOAT_ARR(subBathOnMask(ob_id,:),
     &                                   sNx + sNy + 1)
                      print *, "OB_EXTRACT_LOOKUP: subTempOnMask local:"
                      call PRINT_FLOAT_ARR(subTempOnMask(ob_id,:),
     &                                   sNx + sNy + 1)
                    end if
                    if (tag .eq. 2) then
                      call APPEND_OB_VALUE(sub_global_indices(ob_id,:),
     &                                   subBathOnMask(ob_id,:),
     &                                   global_ob, ob_id, lookup_table)
                      print *, "OB_EXTRACT_LOOKUP:local global_ob:",
     &                       global_ob
                    else if (tag .eq. 3) then
                      call APPEND_OB_VALUE(sub_global_indices(ob_id,:),
     &                                   subTempOnMask(ob_id,:),
     &                                   global_ob, ob_id, lookup_table)
                    end if
                  end if

C                 Receive all sub arrays from processes with pid > 0. Only needs to be done once.
                    call MPI_RECV(recv_arr, nOB_mask*(sNx+sNy+1)
     &                  ,MPI_INTEGER, pid, tag, MPI_COMM_WORLD, status,
     &                  ierror)
                   if (debug .eq. 1) then
                     print *, "Received recv_arr with tag",
     &                        status(MPI_TAG)
                     call PRINT_FLOAT_ARR(recv_arr, nOB_mask*(sNx+sNy+1))
                   end if


C                 Appending all other processes' local arrays to final array if number of valid elements in arrays > 0:
                  if (int(g_indices(pid, ob_id, 1)) .ge. 1) then
                    call APPEND_OB_VALUE(g_indices(pid,ob_id,:),
     &                                  recv_arr(ob_id,:),
     &                                  global_ob, ob_id, lookup_table)
                    total_valid = total_valid + g_indices(pid,ob_id,1)

C                    do i=2, int(g_indices(pid, ob_id, 1)) + 1
C                      ob_index = lookup_table(ob_id,
C     &                           int(g_indices(pid, ob_id, 1)))
C                    end do
                  end if

               end do
               if (tag .eq. 2) then
                 fldNm = "bath"
               else if (tag .eq. 3) then
                 fldNm = "temp"
               end if
C                 call WRITE_GLOBAL_BIN(global_ob, total_valid, fldNm,
C     &                               ob_id, tag)
                print *, "OB_EXTRACT_LOOKUP: updated global ob by ",
     &                   "process", pid
               call PRINT_FLOAT_ARR(global_ob,(sNy+sNx)*(nPx*nPy))
             end do
             print *, "OB_EXTRACT_LOOKUP: final global ob by ",
     &                   "for tag", tag
             call PRINT_FLOAT_ARR(global_ob,(sNy+sNx)*(nPx*nPy))
         end do
         end if
       end if
#endif


       end

       SUBROUTINE APPEND_OB_VALUE(global_indices,ob_values,global_ob,
     &                            ob_id, lookup_table)
         include "SIZE.h"
         include "USER_INPUT.h"
         _RL global_indices(sNx + sNy + 1)
         _RL ob_values(sNx + sNy + 1)
         _RL global_ob((sNy+sNx)*(nPx*nPy))
         _RL lookup_table(nOB_mask, Ny*Nx)
         integer ob_id
         integer num_valid, i, ob_index
         integer, parameter :: print_debug = 0

         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: global_indices:"
           call PRINT_FLOAT_ARR(global_indices, sNx + sNy + 1)
         end if

         num_valid = int(global_indices(1))
         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: num_valid", num_valid
         end if

         do i=2, num_valid + 1
           ob_index = lookup_table(ob_id, int(global_indices(i)))
           if (print_debug .eq. 1) then
             print *, "APPEND_OB_VALUE: ob_index", ob_index
             print *, "APPEND_OB_VALUE: ob_value", ob_values(i)
           end if
           global_ob(ob_index) = ob_values(i)
         end do
         if (print_debug .eq. 1) then
           print *, "APPEND_OB_VALUE: global_ob:"
           call PRINT_FLOAT_ARR(global_ob, (sNy+sNx)*(nPx*nPy))
         end if
       end


       SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
c     Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
       integer arr_length, i
       _RL array(arr_length)
       do, i=1, arr_length
         write(*,"(F6.0,$)") array(i)
           if (i .EQ. arr_length) then
              write(*, '(A,/)') ''
           endif
       enddo
       print *, " "
       end



       SUBROUTINE CREATE_FILENAME(dir, fldNm, ob_id, filename,
     &                            is_child, mpiPx, mpiPy)
c       Creating filenames for all binary files.

       integer is_child, mpiPx, mpiPy, ob_id
       character*40 filename
       character*4 fldNm
       character*5 dir
       character*3 str_px, str_py

c       From EESUPPORT.h
C       mpiPx     - My MPI proc. grid X coord
C       mpiPy     - My MPI proc. grid Y coord

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       001 for index 1 in x axis
c       002 for index 2 in x axis
       if (is_child .eq. 1) then
         write(str_px, "(I3.3)") mpiPx
         write(str_py, "(I3.3)") mpiPy
         write(filename, "(10A)") dir, "ob",ob_id,".",fldNm,".",
     &                            str_px(1:3), ".", str_py(1:3), ".bin"
       else
c         Creating a filename for global arrays
         write(filename, "(6A)") dir, "ob", ob_id,"_",fldNm,
     &                           "_global.bin"
       end if

       end

       SUBROUTINE WRITE_SUB_BIN(mpiMyId, mpiPx, mpiPy, subMask,
     &                          fldNm, ob_id)
C      Writing all sub arrays read from binary fie containing global array
       include "SIZE.h"

       integer mpiMyId, mpiPx, mpiPy, ob_id
       character*4 fldNm
       character*5 dir
       character*40 filename
       _RL subMask(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
       integer, parameter :: print_write = 1

       dir = "./"
       call CREATE_FILENAME(dir, fldNm, ob_id, filename, 1,
     &                      mpiPx, mpiPy)
       open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(1) subMask
       close(1)
       if (print_write .eq. 1) then
         print *, fldNm, "for ob",ob_id,"process", mpiMyId,
     &            "is written to ",filename
       end if

       end


       SUBROUTINE WRITE_GLOBAL_BIN(global_ob, total_valid, prefix)
       include "SIZE.h"

       _RL global_ob(2, (sNy+sNx)*(nPx*nPy))
       integer total_valid

       character*15 prefix
       character*40 filename

c       Writing global tempOnMask to binary file
       prefix = "./T_on_mask1"
       call CREATE_FILENAME(prefix, 0, filename, 0, 0)
       open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(1) global_ob(1, 1:total_valid)
       close(1)

c       Writing global tempOnMask to binary file
       prefix = "./bath_on_mask1"
       call CREATE_FILENAME(prefix, 0, filename, 0, 0)
       open(2, FILE=filename, FORM="unformatted", ACCESS="stream")
       write(2) global_ob(2, 1:total_valid)
       close(2)

       end

       SUBROUTINE PLOT_TEST_GLOBAL_OB(global_ob, myThid,
     &     loc_global_indices)
c      Unravels global_ob and plots onto an array with dimensions Nx by Ny
       include "SIZE.h"
       _RL global_ob((sNy+sNx)*(nPx*nPy))
       _RL global_mask(Nx,Ny,nSx,nSy)
       _RL loc_global_indices(sNx + sNy + 1)
       _RL g_indices(mpiNprocs, nOB_mask, sNx + sNy + 1)

       _RL global_correct((sNy+sNx)*(nPx*nPy))
       _RL global_NxNy(Nx,Ny,nSx,nSy)

       integer i,j, bi, bj, counter, ob_counter, correct, myThid
       integer myByLo_myThid, myByHi_myThid, total_valid
       integer myBxLo_myThid, myBxHi_myThid
       integer, parameter :: plot_on_NxNy = 1

       common / local_var / global_mask, total_valid, g_indices
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


       correct = 1
       counter = 1
       ob_counter = 1
       DO bj = myByLo_myThid, myByHi_myThid
         DO bi = myBxLo_myThid, myBxHi_myThid
           DO j=1, Ny
             DO i=1, Nx
               if (ob_counter .le. total_valid) then
                 if (counter .eq. global_ob(1,ob_counter)) then
C                Inserting bathymetry points into NxNy array
                   global_NxNy(i,j,bi,bj) = global_ob(2,ob_counter)
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
     &            "incorrect."
         call PRINT_FLOAT_ARR(global_correct, (sNy+sNx)*(nPx*nPy))
       else
         print *, "PLOT_TEST_GLOBAL_OB: global_ob final array is ",
     &            "correct."
         if (plot_on_NxNy .eq. 1) then
           print *, "PLOT_TEST_GLOBAL_OB: plotting global_ob ",
     &              "on NxNy dimension array:"
           call PLOT_GLO_FIELD_XYRL(global_NxNy,"PLOT_TEST_GLOBAL_OB",
     &                              myThid)
         end if
       end if
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
