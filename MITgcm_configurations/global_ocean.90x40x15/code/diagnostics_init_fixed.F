#include "DIAG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: DIAGNOSTICS_INIT_FIXED

C     !INTERFACE:
      SUBROUTINE DIAGNOSTICS_INIT_FIXED(myThid)

C     !DESCRIPTION:
C     finish setting up the list of available diagnostics and
C     prepare for storing selected diagnostics and statistics-diags.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"
#include "GRID.h"
#ifdef ALLOW_FIZHI
#include "PARAMS.h"
#endif

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:

C--   Set number of levels for all available diagnostics
C     (cannot add diags to list anymore after this call)
      CALL DIAGNOSTICS_SET_LEVELS( myThid )

C--   Calculate pointers for diagnostics set to non-zero frequency
      CALL DIAGNOSTICS_SET_POINTERS( myThid )

C--   Set parameters and variables used in post-processing diagnostics
      CALL DIAGNOSTICS_SET_CALC( myThid )

C--   Define region-mask for regional statistics diagnostics
      CALL DIAGSTATS_SET_REGIONS( myThid )

C--   Calculate pointers for statistics-diags set to non-zero frequency
      CALL DIAGSTATS_SET_POINTERS( myThid )

      CALL DIAGSTATS_INI_IO( myThid )

      CALL aprilshin( myThid )

      SUBROUTINE aprilshin( myThid )
        !include "../code_mpi/SIZE.h"
        include "/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/"
     &         ,"global_ocean.90x40x15/code_mpi/SIZE.h"
        include "mpif.h"
        !include "/home/mitgcm/Work/MITgcm/verification/tutorial_barotropic_gyre/build_mpi/mpif.h"

        integer myThid

c       local variables:
        integer i,j, iG, jG, ob_subMask_size, pid, ierror
        real global_mask(Ny,Nx)

        !open(10, FILE="aprilshin.txt", FORM="formatted")
        !write(10, "(A)") "hello this is april shin"
        !close(10)


c       Initialize MPI environment
        call MPI_INIT(ierror)


C       myThid  ::  my Thread Id number

c       From EEPARAMS.h :
C       numberOfProcs :: Number of processes computing in parallel
C       pidIO         :: Id of process to use for I/O.
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

C------
C   0) Initialize ob_subMask and global_mask (define an empty domain)
C------


        DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            iG = myXGlobalLo-1+i
            jG = myYGlobalLo-1+j
            global_mask(jG,iG) = 0.0
          ENDDO
        ENDDO


c       _RS field(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
c           iG=bi+(myXGlobalLo-1)/sNx
c           jG=bj+(myYGlobalLo-1)/sNy
c       All Processes read:
        CALL READ_REC_XY_RS( ob_fileName, ob_subMask, 1, 0, myThid )

        ob_subMask_size = (abs(1-OLx)+(sNx+OLx))*(abs(1-OLy)+(sNy+OLy))*nSx*nSy

        !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
        if (myProcId .ne. 0) then
          call MPI_SEND(ob_subMask, ob_subMask_size, MPI_REAL, 0, 1,
     &                  MPI_COMM_WORLD, ierror)
        else
          DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
              DO j=1-OLy,sNy+OLy
                DO i=1-OLx,sNx+OLx
                  iG = myXGlobalLo-1+(bi-1)*sNx+i
                  jG = myYGlobalLo-1+(bj-1)*sNy+j
                  global_mask(jG,iG) = ob_subMask(j,i,bj,bi)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

        end if


        if (myProcId .eq. 0) then

          DO pid=1, numberOfProcs-1, 1
            call MPI_RECV(ob_subMask, ob_subMask_size, MPI_REAL, pid, 1,
     &                    MPI_COMM_WORLD, status, ierror)

            DO bj = myByLo(myThid), myByHi(myThid)
              DO bi = myBxLo(myThid), myBxHi(myThid)
                DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLx
                    iG = myXGlobalLo-1+(bi-1)*sNx+i
                    jG = myYGlobalLo-1+(bj-1)*sNy+j
                    global_mask(jG,iG) = ob_subMask(j,i,bj,bi)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

          CALL PLOT_FIELD_XYRS( global_mask, 'Values on Open Boundary ',
     &                         -1, myThid )

        end if

      call MPI_FINALIZE(ierror)

      end


#ifdef ALLOW_FIZHI
      if( useFIZHI) then
      call fizhi_diagalarms(myThid)
      endif
#endif

      RETURN
      END
