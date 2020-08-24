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

      INTEGER i,j,k,l, bi,bj, fld, ob_id, counter, pid
      character*6 write_fname
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===

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

C#ifdef ALLOW_DIAGNOSTICS_OB
C      IF ( useDiagnostics_ob ) THEN
C        CALL OB_EXTRACT_LOOKUP( myThid )
C      ENDIF
C#endif


C---------------------------------------------------------
C     Initialize diagnostics_ob COMMON block variables:
C---------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS_OB
      DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              DO ob_id=1, nOB_mask
                ob_subMask(ob_id,i,j,bi,bj)   = 0.0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, nOB_mask
        DO l = 1, Ny*Nx
        lookup_table(j,l) = 0
        ENDDO
      ENDDO

      DO i=1, (sNy+sNx)*(nPx*nPy)
        global_ob2D(i) = 0.0
        DO k=1, Nr
          global_ob3D(i,k) = 0.0
        ENDDO
      ENDDO

      DO ob_id=1, nOB_mask
        num_ob_points(ob_id) = 0
        DO pid=1, nPx*nPy
          numOBPnts_allproc(ob_id,pid) = 0
        ENDDO
      ENDDO


      DO ob_id=1, nOB_mask
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
      ENDDO

      time_level = 1


      DO ob_id=1, nOB_mask
        print *,"INIT_VARIA: nOB_mask",nOB_mask
        print*,"ob_id,",ob_id
        print*,"ob_fnames", ob_fnames(ob_id)
        CALL READ_REC_XY_RL( ob_fnames(ob_id),
     &               ob_subMask(ob_id,:,:,:,:),1, 0, myThid )
        _EXCH_XY_RL(ob_subMask(ob_id,:,:,:,:), myThid)
C       Plotting for debugging
        print *,"INIT_VARIA: ob_id", ob_id
        CALL PLOT_FIELD_XYRL( ob_subMask(ob_id,:,:,:,:),
     &                'init_varia:Open Boundary Submask', -1, myThid )
        print *,"IDENTIFY OB POINTS BEFORE1: mask", ob_id

#ifdef ALLOW_USE_MPI
        print *,"IDENTIFY OB POINTS BEFORE2: mask", ob_id
        CALL IDENTIFY_OB_POINTS(ob_id, myThid)
#endif
        print *,"Time to move on to next mask from", ob_id,"to",ob_id+1
      ENDDO

#ifdef ALLOW_USE_MPI
      CALL SEND_CHILD_TO_MASTER()
#endif

C     Using GATHER_2D_R8 to correctly combine sub masks into a single global mask
      DO ob_id=1, nOB_mask
        print *,"INIT_VARIA: ob_id in loop", ob_id
        CALL GATHER_2D_R8(global_ob_mask(ob_id,:,:,1,1),
     &                  ob_subMask(ob_id,1:sNx,1:sNy,1,1),
     &                  Nx, Ny, .true., .true., myThid)
      CALL PLOT_GLO_FIELD_XYRL(global_ob_mask(ob_id,:,:,:,:),
     &                    "INIT_VARIA: global mask", myThid)

       write(write_fname, "(A1,I1,A4)") "m",ob_id,".bin"


#ifdef ALLOW_USE_MPI

      if (mpiMyId .eq. 0) then
C       Store global_ob_masks into lookup table
        counter = 1
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          do j=1, Ny, 1
            do i=1, Nx, 1
              lookup_table(ob_id,counter) =
     &         int(global_ob_mask(ob_id,i,j,bi,bj))
              if ((ob_id .eq. 2).AND.
     &            (global_ob_mask(ob_id,i,j,bi,bj) .gt. 0)) then
               print *,"INIT_V g_m:",
     &            int(global_ob_mask(ob_id,i,j,bi,bj))
              endif
              counter = counter + 1
              if (global_ob_mask(ob_id,i,j,bi,bj) .ge. 1) then
                num_ob_points(ob_id) = num_ob_points(ob_id) + 1
              end if
            end do
           end do
          end do
        end do
      end if
      ENDDO
#endif

C      num_ob_points = 0
C     DO ob_id=1, nOB_mask
C        counter = 1
C        open(ob_id, FILE=ob_filenames(ob_id),FORM="unformatted",
C     &      ACCESS="stream")
C        print *, "global_ob_mask shape", shape(global_ob_mask)
C        print *, "global_ob_mask Nx", Nx, "Ny", Ny
C        read(ob_id) global_ob_mask2(ob_id,:,:)

C        read(ob_id) ob_subMask(ob_id,1:sNx,1:sNy,1,1)

C        print *, "ob_subMask read successful"
C        close(ob_id)

C        print *, "HERE IN VARIA"
C        call PLOT_FIELD_XYRL(ob_subMask(1,:,:,:,:),
C     &                 "OB_EXTRACT_OUTPUT: ob_subMask", myThid)


C      ENDDO


#endif /* ALLOW_DIAGNOSTICS_OB */



      RETURN
      END



      SUBROUTINE IDENTIFY_OB_POINTS(ob_id, myThid)
C     Located open boundary points
C     Saved local i j indices and global indices
C     Saved count of open boundary points in process' piece of global mask
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"

      INTEGER ob_id, myThid, global_index
      INTEGER counter, num_valid, iG, jG, i, j, bi, bj, ob

      iG = mpi_myXGlobalLo(mpiMyId+1)
      jG = mpi_myYGlobalLo(mpiMyId+1)

      print *,"IDENTIFY_OB_POINTS, iG", iG
      print *,"IDENTIFY_OB_POINTS, jG", jG

      counter = 1
      num_valid = 0
      print *,"IDENTIFY_OB_POINTS: ob_id", ob_id
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
           DO i=1,sNx
             global_index = int((j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
     &                       j*(iG-1)+i)
C             print *,"IDENTIFY_OB_POINTS: global_index:", global_index
C             IF (ob_id .eq. 2) then
C              print *,"INIT_VARIA: ob_subMask:",
C     &                (ob_subMask(ob_id,i,j,bi,bj)),
C     &              "num_valid:",num_valid
C             ENDIF
             if (ob_subMask(ob_id,i,j,bi,bj) .ge. 1.0 ) then
               print *,"IDENTIFY_OB_POINTS YES mask",ob_id
C               print *,"IDENTIFY_OB_POINTS: YES",
C     &               ob_subMask(ob_id,:,:,bi,bj)

C     First row contains local i's. Second row contains local j's.
               sub_local_ij_ob(ob_id, 1, counter) = i
               sub_local_ij_ob(ob_id, 2, counter) = j

               print *,"IDENTIFY_OB_POINTS: ij:",i,j, "mask:",ob_id
C              Indices referring to process IDs begin with 1
               sub_glo_indices_allproc(ob_id,mpiMyId+1,counter) =
     &           global_index

C                 print *,"IDENTIFY_OB_POINTS: sub_glo_indices_allproc",
C     &            sub_glo_indices_allproc(ob_id,mpiMyId+1,counter)
               counter = counter + 1
               num_valid = num_valid + 1
             end if
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      print *,"IDENTIFY_OB_POINTS: sub_glo_indices_allproc",
     &    sub_glo_indices_allproc(ob_id,mpiMyId+1,:)
      print *,"INIT_VARIA: num_valid for mask",ob_id,":", num_valid
      numOBPnts_allproc(ob_id,mpiMyId+1) = num_valid

      print*,"INIT_VARIA: mask:",ob_id
      print *,"sub_local_ij_ob j" ,sub_local_ij_ob(ob_id,2,:)
      print *,"sub_local_ij_ob i" ,sub_local_ij_ob(ob_id,1,:)
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
            print *,"init_varia: numOBPnts_allproc from proc",
     &              mpiMyId,":", numOBPnts_allproc
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
