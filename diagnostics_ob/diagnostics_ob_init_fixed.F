#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_INIT_FIXED

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize DIAGNOSTICS_OB variables that are kept fixed during the run.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER i,j,k,l, bi,bj, fld, ob_id, counter, pid
      INTEGER, PARAMETER :: debug = 1
CEOP


#ifdef ALLOW_MNC
      IF (useMNC) THEN
        CALL DIAGNOSTICS_OB_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_OB_INIT( myThid )
      ENDIF
#endif

C-------------------------------------------------------------------------------
C     Initialize diagnostics_ob COMMON block variables:
C-------------------------------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS_OB

      DO i=1, (sNy+sNx)*(nPx*nPy)
        global_ob2D(i) = 0.0
        DO k=1, Nr
          global_ob3D(i,k) = 0.0
        ENDDO
      ENDDO

      DO ob_id=1, nOB_mask

        DO l = 1, Ny*Nx
          lookup_table(ob_id,l) = 0
        ENDDO

        DO pid=1, nPx*nPy
          numOBPnts_allproc(ob_id,pid) = 0
        ENDDO

        DO i=1,2
          DO j=1, sNx+sNy
            sub_local_ij_ob(ob_id, i, j) = 0
          ENDDO
        ENDDO

        DO i=1, sNx+sNy
          DO fld=1, MAX_NFLDS
            subFieldOnMask_2D(ob_id, fld, i) = 0.0
            subFieldOnMask_2Davg(ob_id, fld, i) = 0.0
            DO k=1, Nr
              subFieldOnMask_3D(ob_id, fld, i, k) = 0.0
              subFieldOnMask_3Davg(ob_id, fld, i, k) = 40.0
            ENDDO
          ENDDO
        ENDDO

        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO i=1, Nx
            DO j=1, Ny
              global_ob_mask(ob_id,i,j,bi,bj) = 0.0
            ENDDO
          ENDDO

          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
                ob_subMask(ob_id,i,j,bi,bj)   = 0.0
            ENDDO
          ENDDO
         ENDDO
        ENDDO

        DO i=1, nPx*nPy
          numOBPnts_allproc(ob_id, i) = 0
          DO j=1, sNx+sNy
            sub_glo_indices_allproc(ob_id,i,j) = 0
          ENDDO
        ENDDO

C     !End of loop for nOB_mask
      ENDDO

C     Initializing number of time level passed to 1
      time_level = 1
      time_passed = 0.0

C-----End of COMMON block initialization----------------------------------------



C     Looping through each open boundary mask
      DO ob_id=1, nOB_mask

C       Read open boundary mask files and assign portion of open boundary mask to each process
        CALL READ_REC_XY_RL( ob_fnames(ob_id),
     &               ob_subMask(ob_id,:,:,:,:),1, 0, myThid )
C       Filling in halo regions of sub masks
        _EXCH_XY_RL(ob_subMask(ob_id,:,:,:,:), myThid)

C       Plotting for debugging
        print *,"INIT_FIXED: mask number", ob_id
        CALL PLOT_FIELD_XYRL( ob_subMask(ob_id,:,:,:,:),
     &                'INIT_FIXED:Open Boundary Submask', -1, myThid )

#ifdef ALLOW_USE_MPI
        CALL IDENTIFY_OB_POINTS(ob_id, myThid)
#endif
C     !End of loop for nOB_mask
      ENDDO

#ifdef ALLOW_USE_MPI
      CALL SEND_CHILD_TO_MASTER()
#endif



C -------------------------LOOKUP_TABLE INITIALIZATION--------------------------

C MAIN IDEA: The Master process stores all global open boundary masks into the
C lookup_table.
C
C METHOD:
C 1. Uses GATHER_2D_R8 to correctly combine sub masks in to a single global mask
C    because when the open boundary mask is read it is only saved as sub masks
C    for each child process.
C
C 2. Stores each global mask into each ROW of lookup table
C    - Doing so essentially rehapes the global mask into a one dimensional array
C      in the lookup table.
C    - The indices of the elements of the global mask in each row of the lookup
C      table starts from 1 and ends at Nx*Ny.
C----------------------------

#ifdef ALLOW_USE_MPI
      DO ob_id=1, nOB_mask
C       Using GATHER_2D_R8 to correctly combine sub masks into a single global mask
        CALL GATHER_2D_R8(global_ob_mask(ob_id,:,:,1,1),
     &                  ob_subMask(ob_id,1:sNx,1:sNy,1,1),
     &                  Nx, Ny, .true., .true., myThid)
C       Plotting for debugging
        CALL PLOT_GLO_FIELD_XYRL(global_ob_mask(ob_id,:,:,:,:),
     &                    "INIT_FIXED: global mask", myThid)

      IF (mpiMyId .eq. 0) then
C       Store global_ob_masks into lookup table
        counter = 1
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          do j=1, Ny, 1
            do i=1, Nx, 1
C             Storing elements of each mask in global_ob_mask into each row of lookup_table
              lookup_table(ob_id,counter) =
     &         int(global_ob_mask(ob_id,i,j,bi,bj))
              counter = counter + 1
            end do
           end do
          end do
        end do
C     !ENDIF for mpiMyId .eq. 0
      ENDIF

C     !End of loop for nOB_mask
      ENDDO
#endif
C-----End LOOKUP_TABLE initialization-------------------------------------------

#endif /* ALLOW_DIAGNOSTICS_OB */


      RETURN
      END



      SUBROUTINE IDENTIFY_OB_POINTS(ob_id, myThid)
C     Located open boundary points
C     Saved local i j indices and global indices
C     Saved count of open boundary points in process' piece of global mask
C
C     Global Indices:
C     - Global indices are the positions of elements in a Nx by Ny array.
C       Position number starts from 1 and ends at Nx*Ny
C
C       [ 1,     2,     3,     4     ...     ,Nx]
C       [ Nx+1,  Nx+2,  Nx+3,  Nx+4, ...  ,Nx+Nx]
C       [ 2Nx+1, 2Nx+2, 2Nx+3, 2Nx+4,... ,2Nx+Nx]
C       [ 3Nx+1, 3Nx+2, 3Nx+3, 3Nx+4,... ,3Nx+Nx]
C       [ :                                     ]
C       [ :                                     ]
C       [ :                                     ]
C       [ (Ny-1)*Nx+1,        ...         ,Nx*Ny]
C
C     - We compute the global indices of open boundary points we find in each
C       child process in their assigned sub mask.

C       To do so, we essentially count the number of elements starting from the
C       first row and column of the Nx by Ny array until the element right before
C       the current element we are computing the global index for.
C
C       Example:
C       [ 1,     2,     3,               4,    ...    ,Nx]
C       [ Nx+1,  Nx+2,  Nx+3,            Nx+4, ... ,Nx+Nx]
C       [ 2Nx+1, 2Nx+2, current element                  ]
C       [                                                ]
C       [                                                ]
C       [                                                ]
C       [                                                ]
C
C       Total count of elements before current element = 2Nx+2.
C       Therefore, global index of current element is 2Nx+3.
C
C       However, we must compute the global indices of elements only located in
C       a submask without iterating through the entire Nx by Ny array (need to
C       reduce runtime).
C       Therefore, we first identify the i,j indices of the row 1 and column 1
C       of the sub mask in the Nx by Ny global mask with MITgcm's mpi_myXGlobalLo
C       and mpi_myYGlobalLo.
C       With the aid of these i,j indices in addition to knowing the dimensions
C       of the sub mask, we are able to know the number of rows before the row
C       containing point P of example with respect to the Nx by Ny global mask
C       as well as the number of elements before P in the same row.
C
C       Example: * are elements of sub mask inside Nx by Ny global mask:
C                P is an open boundary point for which we compute the global index
C       [               |                        ]
C       [            r rows                      ]
C       [               |                        ]
C       [           * * | * * * * *              ]
C       [---c points--- P * * * * *              ]
C       [           * * * * * * * *              ]
C       [           * * * * * * * *              ]
C       [           * * * * * * * *              ]
C       [           * * * * * * * *              ]
C       [                                        ]
C       [                                        ] :Nx by Ny array
C
C       Then, as we iterate through the points in sNx by sNy sub mask to locate
C       the open boundary points, we can apply the following formula to compute
C       the corresponding global indices:
C       where i denotes the point along the x axis of sNx by sNy submask
C             j denotes the point along the y axis of sNx by sNy submask
C             iG denotes the global i index of row 1 col 1 of submask
C             jG denotes the global j index of row 1 col 1 of submask
C
C       global_index = (r * Nx) + (c + 1)
C
C       Total number of elements from all rows above: r * Nx = (jG + (j-2))*Nx
C       Number of points including P in th same row: c + 1 = iG + (i-1)
C
C       Therefore,
C               global_index = (jG + (j-2))*Nx + iG + (i-1)
C
C     - These global indices will be used to extract values from the lookup table.

#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

      INTEGER ob_id, myThid, global_index
      INTEGER counter, num_valid, iG, jG, i, j, bi, bj, ob

      iG = mpi_myXGlobalLo(mpiMyId+1)
      jG = mpi_myYGlobalLo(mpiMyId+1)

      counter = 1
      num_valid = 0
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
           DO i=1,sNx

             if (ob_subMask(ob_id,i,j,bi,bj) .ge. 1.0 ) then
C              Computing global index
C              old formula: works correctly, just too complicated
C              global_index = int(((j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
C     &                       j*(iG-1)+i)+((jG-1)*Nx))
C              new formula: works correctly, simple
               global_index = int((jG + (j-2))*Nx + iG + (i-1))

C              First row contains local i's. Second row contains local j's.
               sub_local_ij_ob(ob_id, 1, counter) = i
               sub_local_ij_ob(ob_id, 2, counter) = j

C              Indices referring to process IDs begin with 1
               sub_glo_indices_allproc(ob_id,mpiMyId+1,counter) =
     &           global_index

               counter = counter + 1
               num_valid = num_valid + 1
             end if
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END


      SUBROUTINE SEND_CHILD_TO_MASTER()
C     Sending global indices and number of open boundary points from processes > 0 to master process
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

        INTEGER pid, ob_id, i, status(MPI_STATUS_SIZE), ierror
        INTEGER g_indices(sNx+sNy)
        INTEGER num_pts

        IF (mpiMyId .gt. 0) then
          DO ob_id=1, nOB_mask
            print *,"INIT_FIXED: numOBPnts_allproc from proc",
     &              mpiMyId,":", numOBPnts_allproc
            print *,"INIT_FIXED: BEFORE sub_glo_indices_allproc",
     &              sub_glo_indices_allproc(ob_id,mpiMyId+1,:)
C           !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
            call MPI_SEND(sub_glo_indices_allproc(ob_id,mpiMyId+1,:),
     &             sNx+sNy, MPI_INT,0, ob_id, MPI_COMM_MODEL, ierror)
            call MPI_SEND(numOBPnts_allproc(ob_id,mpiMyId+1),
     &             1, MPI_INT,0, ob_id*100, MPI_COMM_MODEL, ierror)
          ENDDO
        ENDIF

        IF (mpiMyId .eq. 0) then
          DO pid=1, nPx*nPy-1
            DO ob_id=1, nOB_mask
C             !MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierror
              call MPI_RECV(g_indices, sNx+sNy, MPI_INT,
     &                   pid, ob_id ,MPI_COMM_MODEL, status, ierror)
              call MPI_RECV(num_pts, 1, MPI_INT,
     &                   pid, ob_id*100 ,MPI_COMM_MODEL, status, ierror)

              DO i=1, sNx+sNy
               sub_glo_indices_allproc(ob_id,pid+1,i) = g_indices(i)
              ENDDO
              numOBPnts_allproc(ob_id,pid+1) = num_pts
            ENDDO
          ENDDO
         ENDIF


       END
