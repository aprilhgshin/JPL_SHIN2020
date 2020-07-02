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

       SUBROUTINE OB_EXTRACT( myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#ifdef ALLOW_FIZHI
#include "PARAMS.h"
#endif
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"


c       IMPLICIT NONE
C#include "EEPARAMS.h"
C#include "EESUPPORT.h"
C#ifdef ALLOW_FIZHI
C#include "PARAMS.h"
C#endif
C#include "GRID.h"
C#include "DIAGNOSTICS_SIZE.h"
C#include "DIAGNOSTICS.h"

c No need to include "mpif.h" becaues EESUPPORT.h already includes it.

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678


       integer myThid


C      local variables:
#ifdef ALLOW_USE_MPI
       integer pid, ierror, np, itemp(2)
#endif

       integer i,j, iG, jG, ob_subMask_size
       integer, parameter :: stop_program = 0
       integer, parameter :: debug = 0
C      I thought global_mask should be real*4 as we are reading in float32 values
C      however, it causes errors when GATHER_2D_R4 is called because its input is _RS. So it is set as real*8 or _RS.
       real*8 global_mask(Nx,Ny,nSx,nSy)

       CHARACTER*(MAX_LEN_MBUF) msgBuf

C      Testing print:
       if (debug .eq. 1) then
         print *, "Printing from diagnostics_init_fixed: April Shin"
         write(msgBuf, "(A)") "April Shin Open Boundary Test"
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                              SQUEEZE_RIGHT, myThid)
C        Testing to see if code is even running:
         open(10, FILE="aprilshin.txt", FORM="formatted")
         write(10, "(A)") "Written from diagnostics_init_fixed.F"
         close(10)
       end if

C#ifdef ALLOW_USE_MPI
C      Initialize MPI environment
C       call MPI_INIT(ierror)
C#endif

C       myThid  ::  my Thread Id number

c       From EEPARAMS.h :
C       numberOfProcs :: Number of processes computing in parallel
C       myProcId      :: My own "process" id.
C       myPx          :: My X coord on the proc. grid.
C       myPy          :: My Y coord on the proc. grid.
C       myXGlobalLo   :: My bottom-left (south-west) x-index global domain.
C                      The x-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C       myYGlobalLo   :: My bottom-left (south-west) y-index in global domain.
C                      The y-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C
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
C    Initialize global_mask (define an empty domain)
c     GATHER_2D_R4 initializes for us:
CC       DO j=1,Ny
CC         DO i=1,Nx
CC           global_mask(j,i) = 0.0
CC         ENDDO
CC       ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy

      print *, "April Shin: ob_subMask indices:",
     & 1-Olx, sNx+Olx, 1-Oly, sNy+Oly
      print *, "April Shin ob_subMask shape: ", shape(ob_subMask)
c       All Processes read:

       ob_fileName = 'bathymetry.bin'
       print *, "April Shin file name: ", ob_fileName
       CALL READ_REC_XY_RS( ob_fileName, ob_subMask, 1, 0, myThid )

       _EXCH_XY_RS(ob_subMask, myThid)



      ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*(abs(1-OLy)+(sNy+OLy))
     &                    *nSx*nSy
      print *, "obsubMask at -2: ", ob_subMask(-2,-2,1,1)
C     Plotting for debugging
      CALL PLOT_FIELD_XYRS( ob_subMask, 'April Shin Submask',
     &                         -1, myThid )


       !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
c       if (mpiMyId .ne. 0) then
c         call MPI_SEND(ob_subMask, ob_subMask_size, MPI_REAL, 0, 1,
c     &                 MPI_COMM_WORLD, ierror)
c      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C#ifdef ALLOW_USE_MPI
CC       print*,"ob_subMask for proc", myProcId
CC       call PRINT_FLOAT_ARR(ob_subMask,ob_subMask_size)
C12345678
       print *, "April Shin MPI Process ID: ", mpiMyId
       print *, "April Shin MPI total num procs: ", mpiNprocs
       print *, "April Shin myThid: ", myThid




C      call PLOT_GLO_FIELD_XYRS(global_mask, "April global mask",
C     &                            PRINT_MAP_XY)

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

         call GATHER_2D_R8(global_mask(1:Nx,1:Ny,1,1),
     &                     ob_subMask(1:sNx, 1:sNy, nSx, nSy),
     &                     Nx, Ny, .true., .true., myThid)
         if (mpiMyId .eq. 0) then
           print *, "April Shin global_mask shape:", shape(global_mask)
           call PLOT_GLO_FIELD_XYRS(global_mask, "April global mask",
     &                            myThid)
         end if




C      This plotting routine only writes fields with dimensions equal to ob_subMask (myIterL -1)
C       CALL PLOT_FIELD_XYRS( ob_subMask, 'April Shin Submask',
C     &                         2, myThid )

CC       if (myProcId .eq. 0) then
C     To speed-up mpi gather and scatter routines, myXGlobalLo
C     and myYGlobalLo from each process are transferred to
C     a common block array.  This allows process 0 to know
C     the location of the domains controlled by each process.
C         DO np = 1, nPx*nPy-1
C            itemp(1) = myXGlobalLo
C            itemp(2) = myYGlobalLo
C            pid = np - 1
C            CALL MPI_BCAST(itemp, 2, MPI_INTEGER, pid,
C     &                     MPI_COMM_MODEL, ierror)
C            mpi_myXGlobalLo(np) = itemp(1)
C            mpi_myYGlobalLo(np) = itemp(2)
C            print *,"mpi_myXGlobalLo",mpi_myXGlobalLo(np)
C            print *,"mpi_myYGlobalLo",mpi_myYGlobalLo(np)

C         ENDDO
C      end if
C#endif

C        Inputting elements from ob_subMask read by process 0:
C         DO bj = myByLo(myThid), myByHi(myThid)
C           DO bi = myBxLo(myThid), myBxHi(myThid)
C             DO j=1-OLy,sNy+OLy
C               DO i=1-OLx,sNx+OLx
C                 iG = mpi_myXGlobalLo(0)-1+(bi-1)*sNx+i
C                 jG = mpi_myYGlobalLo(0)-1+(bj-1)*sNy+j
C                 global_mask(jG,iG) = ob_subMask(j,i,bj,bi)
C               ENDDO
C             ENDDO
C           ENDDO
C         ENDDO

C     !INPUT/OUTPUT PARAMETERS FOR GATHER_2D_R4:
C gloBuff   ( _R4 ) :: full-domain 2D IO-buffer array             (Output)
C myField   ( _R4 ) :: tiled, local (i.e. my Proc. tiles) 2D array (Input)
C xSize    (integer):: global buffer 1rst dim (x)
C ySize    (integer):: global buffer 2nd  dim (y)
C useExch2GlobLayOut:: =T: Use Exch2 global-map layout (only with EXCH2)
C zeroBuff (logical):: =T: initialise the buffer to zero before copy
C myThid   (integer):: my Thread Id number
C12345678




CC         print *, "April Shin global_mask: "
CC         call PRINT_FLOAT_ARR(global_mask, Ny*Nx)


C    DO pid=1, numberOfProcs-1, 1
C      call MPI_RECV(ob_subMask, ob_subMask_size, MPI_REAL, pid, 1,
C&                    MPI_COMM_WORLD, status, ierror)

C      DO bj = myByLo(myThid), myByHi(myThid)
C        DO bi = myBxLo(myThid), myBxHi(myThid)
C          DO j=1-OLy,sNy+OLy
C            DO i=1-OLx,sNx+OLx
C              iG = myXGlobalLo-1+(bi-1)*sNx+i
C              jG = myYGlobalLo-1+(bj-1)*sNy+j
C              global_mask(jG,iG) = ob_subMask(j,i,bj,bi)
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO

C    ENDDO




CC      end if
C#ifdef ALLOW_USE_MPI
C       call MPI_FINALIZE(ierror)
C#endif
       end


       SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
c     Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
       integer arr_length, i
       real*4 array(arr_length)
       do, i=1, arr_length
         write(*,"(F4.0,$)") array(i)
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
