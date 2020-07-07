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
       integer, parameter :: run_program = 0
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


C       myThid  ::  my Thread Id number

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

      if (run_program .eq. 1) then

        print *, "April Shin: ob_subMask indices:",
     &            1-Olx, sNx+Olx, 1-Oly, sNy+Oly
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

         print *, "April Shin MPI Process ID: ", mpiMyId
         print *, "April Shin MPI total num procs: ", mpiNprocs
         print *, "April Shin myThid: ", myThid

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

      end if

       end
