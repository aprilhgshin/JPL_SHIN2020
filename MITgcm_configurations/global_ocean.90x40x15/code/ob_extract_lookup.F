C All variables in header file are global to header files included below

C     == Global variables ===
C     ORDER SHOWN IN OTHER .F FILES:

C#include "SIZE.h"
C#include "EEPARAMS.h"
C#include "PARAMS.h"
C#include "DIAGNOSTICS_SIZE.h"
C#include "DIAGNOSTICS.h"

C     == Global variables ===
C#include "EEPARAMS.h"
C#include "SIZE.h"
C#include "DIAGNOSTICS_SIZE.h"
C#include "DIAGNOSTICS.h"

C#include "SIZE.h"
C#include "EEPARAMS.h"
C#include "PARAMS.h"
C#include "GRID.h"
C#include "DIAGNOSTICS_SIZE.h"
C#include "DIAGNOSTICS.h"

C#include "EEPARAMS.h"
C#include "SIZE.h"
C#include "DIAGNOSTICS_SIZE.h"
C#include "PARAMS.h"
C#include "DIAGNOSTICS.h"

C#include "SIZE.h"
C#include "EEPARAMS.h"
C#include "EESUPPORT.h"

#include "DIAG_OPTIONS.h"


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678

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


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678

       integer myThid


C      !LOCAL VARIABLES:
C      global_mask: 4 dimensional REAL*8 array that combines all open boundarys
C      msgBuf: Used to write statements into the output file
C      ob_subMask_size: Integer size of ob_subMask (ob_subMask from DIAGNOSTICS.h)
C      global_ob_index: Integer indicating bottom-left global element number
C      stop_program: Integer used to hold either 1 or 0, 1 to stop program or 0 to ignore
C      debug: Integer used to hold either 1 or 0, 1 to print debug statements or 0 to ignore

C      Although we read REAL*4 values, most built in functions have _RS input which by default is real*8.
C      So global_mask, ob_subMask, and bath_subMask are set as real*8 or _RS.
       real*8 global_mask(Nx,Ny,nSx,nSy)

       CHARACTER*(MAX_LEN_MBUF) msgBuf
       integer i,j, iG, jG, bi, bj
       integer ob, global_ob_index, ob_index, new_goi
       integer ob_subMask_size, counter
       integer total_valid1
       integer, parameter :: num_openB = 1
       integer, parameter :: stop_program = 0
       integer, parameter :: debug = 1
       CHARACTER*40 b_fileName
#ifdef ALLOW_USE_MPI
       integer pid, ierror, status(MPI_STATUS_SIZE)
C     &    , np, itemp(2)
#endif



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



c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy

       if (debug .eq. 1) then
         print *, "April Shin: ob_subMask indices:",
     &            1-Olx, sNx+Olx, 1-Oly, sNy+Oly
         print *, "April Shin ob_subMask shape: ", shape(ob_subMask)
       end if

c       All Processes read:

       b_fileName = 'bathymetry.bin'
       t_fileName = 'domain_temp.bin'
       ob_fileName = 'domain_flt32_mask.bin'
C       CALL READ_REC_XY_RS( t_fileName, temp_subMask, 1, 0, myThid )
       CALL READ_REC_XY_RS( b_fileName, bath_subMask, 1, 0, myThid)
       CALL READ_REC_XY_RS( ob_fileName, ob_subMask, 1, 0, myThid)


       _EXCH_XY_RS(bath_subMask, myThid)
       _EXCH_XY_RS(ob_subMask, myThid)

       iG = mpi_myXGlobalLo(mpiMyId+1)
       jG = mpi_myYGlobalLo(mpiMyId+1)
       global_ob_index = Nx*(jG -1)  + iG

       if (debug .eq. 1) then
         print *, "process id: ", mpiMyId
         print *, "iG:", iG, "jG:", jG
         print *, "OB_EXTRACT_LOOKUP:global_ob_index", global_ob_index
       end if

       counter = 2
       new_goi = global_ob_index



      if (debug .eq. 1) then
C        Plotting for debugging
        print *, "OB_EXTRACT_LOOKUP: Number of valid elements:",
     &            counter-2
        CALL PLOT_FIELD_XYRS( ob_subMask, 'Open Boundary Submask',
     &                         -1, myThid )
C         call PRINT_FLOAT_ARR(ob_subMask, ob_subMask_size)

        CALL PLOT_FIELD_XYRS( bath_subMask, 'Bathymetry Submask',
     &                         -1, myThid )
C         call PRINT_FLOAT_ARR(bath_subMask, ob_subMask_size)

C        print *, "April subBathOnMask in process", mpiMyId,
C     &            "before sending to Parent:"
C         call PRINT_FLOAT_ARR(subBathOnMask, num_openB*(sNx+sNy+1))
      end if

       if (debug .eq .1) then
         print *, "OB_EXTRACT_LOOKUP: MPI total num procs: ",mpiNprocs
         print *, "OB_EXTRACT_LOOKUP: myThid: ", myThid
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
               if (ob_subMask(i,j,bi,bj) .gt. 0.0 ) then
                 if (debug .eq. 1) then
                   print *, "OB_EXTRACT_LOOKUP: ob_subMask gt 1:",
     &                    ob_subMask(i,j,bi,bj)
                   print *, "OB_EXTRACT_LOOKUP: i j of ob_subMask gt 1",
     &                     i, j
                   print *, "OB_EXTRACT_LOOKUP: bath_subMask gt 1",
     &                    bath_subMask(i,j,bi,bj)
                 end if
                 subBathOnMask(1, counter) = bath_subMask(i,j,bi,bj)
                 sub_global_indices(1, counter) = new_goi
                 if (debug .eq. 1) then
                   print *, "OB_EXTRACT_LOOKUP: new_goi:", new_goi
                 end if
                 counter = counter + 1
               end if

             ENDDO
           ENDDO
         ENDDO
       ENDDO


      ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*(abs(1-OLy)+(sNy+OLy))
     &                    *nSx*nSy

      subBathOnMask(1, 1) = counter - 2
      sub_global_indices(1, 1) = counter - 2



#ifdef ALLOW_USE_MPI
       if (usingMPI) then
         if (mpiMyId .gt. 0) then

C#ifndef DISABLE_MPI_READY_TO_RECEIVE
           call MPI_SEND(sub_global_indices,num_openB*(sNx+sNy+1),
     &                MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
           call MPI_SEND(subBathOnMask, num_openB*(sNx+sNy+1), MPI_REAL
     &               ,0, 2, MPI_COMM_MODEL, ierror)

         end if
       end if
C#endif
#endif



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

       call GATHER_2D_R8(global_mask(1:Nx,1:Ny,1,1),
     &                   ob_subMask(1:sNx, 1:sNy, 1, 1),
     &                   Nx, Ny, .true., .true., myThid)




#ifdef ALLOW_USE_MPI
      if (usingMPI) then
        if (mpiMyId .eq. 0) then
C         print *, "April Shin global_mask shape:", shape(global_mask)
          call PLOT_GLO_FIELD_XYRS(global_mask,
     &                      "OB_EXTRACT_LOOKUP: global mask", myThid)
C          Store global_ob_masks into lookup table
           do ob=1, num_openB, 1

             counter = 1
             do j=1, Ny, 1
               do i=1, Nx, 1
                 lookup_table(ob, counter) = global_mask(i,j,1,1)
                 if (lookup_table(ob, counter) .ge. 1) then
                   print *, "OB_EXTRACT_LOOKUP:global index positive:",
     &              lookup_table(ob, counter), "at", counter
                 end if
                 counter = counter + 1
               end do
             end do

C            Append subBathOnMask values to global_ob1 if valid elements exist.
             if (sub_global_indices(ob, 1) .gt. 0) then
               print *, "OB_EXTRACT_LOOKUP: local subBathOnMask has",
     &                     "valid elements:"
               call PRINT_FLOAT_ARR(subBathOnMask, sNx+sNy+1)
               do i=2, sub_global_indices(ob, 1) + 1
                 ob_index = lookup_table(ob,
     &                              sub_global_indices(ob, i))
                 if (debug .eq. 1) then
                   print *, "OB_EXTRACT_LOOKUP: sub_global_index:",
     &                     sub_global_indices(ob, i)
                   print *, "OB_EXTRACT_LOOKUP:process 0 ob_index:",
     &                     ob_index
                   print *, "OB_EXTRACT_LOOKUP:subBathOnMask at ob ",
     &                     "points:" ,subBathOnMask(ob, i)
                 end if
                 global_ob1(1, ob_index) = subBathOnMask(ob, i)
               end do
             end if

           end do

           if (debug .eq. 1) then
             print *, "OB_EXTRACT_LOOKUP:Locally added: global_ob1:"
             call PRINT_FLOAT_ARR(global_ob1(1,1:(sNy+sNx)*(nPx*nPy)),
     &                          1*(sNy+sNx)*(nPx*nPy))
           end if

C          Retrive sub_global_indicies and subBathOnMask from every process > 0.
           ob_index = 1

           do ob=1, num_openB
             do pid=1, mpiNprocs-1

              call MPI_RECV(sub_global_indices, num_openB*(sNx+sNy+1)
     &                  ,MPI_INTEGER, pid, 1, MPI_COMM_WORLD, status,
     &                  ierror)
              call MPI_RECV(subBathOnMask, num_openB*(sNx+sNy+1),
     &             MPI_REAL,pid, 2 ,MPI_COMM_MODEL, status, ierror)

               if (sub_global_indices(ob, 1) .ge. 1) then
                 do i=2, sub_global_indices(ob, 1) + 1
                   ob_index = lookup_table(ob,
     &                              sub_global_indices(ob, i))
                   if (debug .eq. 1) then
                     print *, "OB_EXTRACT_LOOKUP: ob_index:",ob_index
                   end if
                   global_ob1(1, ob_index) = subBathOnMask(ob, i)
                 end do
               end if
               if (debug .eq. 1) then
                 print *, "OB_EXTRACT_LOOKUP:Parent received from",
     &                  " Process",status(MPI_SOURCE)
                 print *, "OB_EXTRACT_LOOKUP:subBathOnMask from ",
     &                  "Process:",status(MPI_SOURCE)
               end if

               call PRINT_FLOAT_ARR(subBathOnMask, sNx+sNy+1)
               print *, "OB_EXTRACT_LOOKUP:Process",pid,"sub array",
     &                 " appended to global_ob1:"
               call PRINT_FLOAT_ARR(global_ob1, 1*(sNy+sNx)*(nPx*nPy))

            end do
          end do
          print *, "OB_EXTRACT_LOOKUP: FINAL global_ob1:"
          call PRINT_FLOAT_ARR(global_ob1, 1*(sNy+sNx)*(nPx*nPy))
        end if
      end if
#endif
C     !DESCRIPTION:
C     *=====================================================================*
C     | SUBROUTINE BAR2
C     | o Global barrier routine.
C     *=====================================================================*
C     | Implements a simple true shared memory barrier that uses a global
C     | heap array that all threads can access to synchronise. Each thread
C     | writes to a predefined location. One thread polls the locations. Other
C     | threads poll an all clear assertion location. Once the polling
C     | thread that is looping over locations sees writes for each thread is
C     | writes the all clear assertion location and everyone proceeds. A
C     | cyclic series of locations is used to ensure that race conditions do
C     | not occur. A few simple statistics are recorded giving number of
C     | barrier calls, max, min and aggregate polling loop counts.
C     *=====================================================================*
C         call BAR2(myThid)

       end



       SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
c     Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
       integer arr_length, i
       real*8 array(arr_length)
       do, i=1, arr_length
         write(*,"(F6.0,$)") array(i)
           if (i .EQ. arr_length) then
              write(*, '(A,/)') ''
           endif
       enddo
       print *, " "
       end




       SUBROUTINE PLOT_GLO_FIELD_XYRS(fld, fldNam, myThid)

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
       real*8 fld(Nx, Ny, nSx, nSy)

C       _RS fld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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
        CALL PRINT_MAPRS(
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
