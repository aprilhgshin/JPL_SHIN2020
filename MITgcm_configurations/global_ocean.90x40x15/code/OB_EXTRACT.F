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



C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678

       SUBROUTINE OB_EXTRACT( myThid )

c       IMPLICIT NONE
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#ifdef ALLOW_FIZHI
#include "PARAMS.h"
#endif
#include "GRID.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"

#ifdef ALLOW_USE_MPI
#include "mpif.h"
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C12345678


       integer myThid

C      local variables:
       integer i,j, iG, jG, ob_subMask_size, pid, ierror, np
       integer, parameter :: stop_program = 1
       integer, parameter :: debug = 1
       real global_mask(Ny,Nx)
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

#ifdef ALLOW_USE_MPI
c      Initialize MPI environment
       call MPI_INIT(ierror)
#endif

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
       DO j=Ny
         DO i=Nx
           global_mask(j,i) = 0.0
         ENDDO
       ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c      ob_subMask initialized in diagnostics_readparms, dimensions match READ_REC_XY_RS field:
c      field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy


c       All Processes read:

       CALL READ_REC_XY_RS( ob_fileName, ob_subMask, 1, 0, myThid )

       ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*(abs(1-OLy)+(sNy+OLy))*nSx*nSy

c      To test if code is running to this point, print stop message if stop_program is set to 1
       if (stop_program .eq. 1) then
         stop "init_f"
       end if

C12345678
       !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
c       if (mpiMyId .ne. 0) then
c         call MPI_SEND(ob_subMask, ob_subMask_size, MPI_REAL, 0, 1,
c     &                 MPI_COMM_WORLD, ierror)
c      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       if (myProcId .eq. 0) then
C     To speed-up mpi gather and scatter routines, myXGlobalLo
C     and myYGlobalLo from each process are transferred to
C     a common block array.  This allows process 0 to know
C     the location of the domains controlled by each process.
C         DO np = 1, nPx*nPy
C            itemp(1) = myXGlobalLo
C            itemp(2) = myYGlobalLo
C            pid = np - 1
C            CALL MPI_BCAST(itemp, 2, MPI_INTEGER, pid,
C     &                     MPI_COMM_MODEL, ierror)
C            mpi_myXGlobalLo(np) = itemp(1)
C            mpi_myYGlobalLo(np) = itemp(2)
C         ENDDO

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
         call GATHER_2D_R4(global_mask, ob_subMask, Nx, Ny, .true.,
     &                     .true., myThid)


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
         CALL PLOT_FIELD_XYRS( global_mask, 'Values on Open Boundary ',
     &                         -1, myThid )

       end if

#ifdef ALLOW_USE_MPI
       call MPI_FINALIZE(ierror)
#endif
       end
