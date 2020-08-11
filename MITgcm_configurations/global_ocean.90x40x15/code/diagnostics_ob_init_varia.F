#include "DIAGNOSTICS_OB_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: DIAGNOSTICS_OB_INIT_VARIA

C !INTERFACE:
      SUBROUTINE DIAGNOSTICS_OB_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize DIAGNOSTICS_OB variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIAGNOSTICS_OB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS_OB

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER i,j,k,l, bi,bj, fld, ob_id, counter, pid

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Initialize arrays in common blocks :
C     3-dim.
        DO k=1,Nr
         DO j=1-Oly,sNy+OLy
          DO i=1-Olx,sNx+Olx
#ifdef DIAGNOSTICS_OB_3D_STATE
           diagOB_StatScal1(i,j,k,bi,bj) = 0. _d 0
           diagOB_StatScal2(i,j,k,bi,bj) = 0. _d 0
           diagOB_StatVelU (i,j,k,bi,bj) = 0. _d 0
           diagOB_StatVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* DIAGNOSTICS_OB_3D_STATE */
#ifdef DIAGNOSTICS_OB_TENDENCY
           diagOB_TendScal1(i,j,k,bi,bj) = 0. _d 0
           diagOB_TendScal2(i,j,k,bi,bj) = 0. _d 0
           diagOB_TendVelU (i,j,k,bi,bj) = 0. _d 0
           diagOB_TendVelV (i,j,k,bi,bj) = 0. _d 0
#endif /* DIAGNOSTICS_OB_TENDENCY */
          ENDDO
         ENDDO
        ENDDO

C     2-dim.
#ifdef DIAGNOSTICS_OB_2D_STATE
        DO j=1-Oly,sNy+OLy
         DO i=1-Olx,sNx+Olx
          diagOB_Surf1(i,j,bi,bj) = 0. _d 0
          diagOB_Surf2(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* DIAGNOSTICS_OB_2D_STATE */

C- end bi,bj loops
       ENDDO
      ENDDO


C ----END OF DIAG_OB_EXTRACT COMMON block initializing ---------------

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

C     Prevent loading of new data before everyone has finished initialisation
C     this is done within IO routines => no longer needed
c      _BARRIER

#ifdef DIAGNOSTICS_OB_3D_STATE
       IF ( diagOB_Scal1File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_Scal1File, ' ',
     &                         diagOB_StatScal1, 0, myThid )
         CALL EXCH_3D_RL( diagOB_StatScal1, Nr, myThid )
       ENDIF
       IF ( diagOB_Scal2File .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_Scal2File, ' ',
     &                         diagOB_StatScal2, 0, myThid )
         CALL EXCH_3D_RL( diagOB_StatScal2, Nr, myThid )
       ENDIF

       IF ( diagOB_VelUFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_VelUFile, ' ',
     &                         diagOB_StatVelU, 0, myThid )
       ENDIF
       IF ( diagOB_VelVFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL( diagOB_VelVFile, ' ',
     &                         diagOB_StatVelV, 0, myThid )
       ENDIF
       IF ( diagOB_VelUFile .NE. ' ' .OR.
     &      diagOB_VelVFile .NE. ' ' ) THEN
        IF ( diagOB_StaV_Cgrid ) THEN
         CALL EXCH_UV_3D_RL( diagOB_StatVelU, diagOB_StatVelV,
     &                       .TRUE., Nr, myThid )
        ELSE
C-      Assume Agrid position:
         CALL EXCH_UV_AGRID_3D_RL( diagOB_StatVelU, diagOB_StatVelV,
     &                       .TRUE., Nr, myThid )
        ENDIF
       ENDIF
#endif /* DIAGNOSTICS_OB_3D_STATE */

#ifdef DIAGNOSTICS_OB_2D_STATE
       IF ( diagOB_Surf1File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( diagOB_Surf1File,' ',
     &                        diagOB_Surf1, 0, myThid )
         CALL EXCH_XY_RL( diagOB_Surf1, myThid )
       ENDIF
       IF ( diagOB_Surf2File .NE. ' ' ) THEN
         CALL READ_FLD_XY_RL( diagOB_Surf2File,' ',
     &                        diagOB_Surf2, 0, myThid )
         CALL EXCH_XY_RL( diagOB_Surf2, myThid )
       ENDIF
#endif /* DIAGNOSTICS_OB_2D_STATE */

C      ELSE
C--   restart from a pickup:
C         CALL DIAGNOSTICS_OB_READ_PICKUP( nIter0, myThid )

C--   end start-from-iter-zero if/else block
      ENDIF

C---------------------------------------------------------
C     Initialize diagnostics_ob COMMON block variables:
C---------------------------------------------------------

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


      DO ob_id=1, nOB_mask
        num_ob_points(ob_id) = 0
        DO pid=1, nPx*nPy
          numOBPnts_allproc(ob_id,pid) = 0
        ENDDO
      ENDDO


      DO ob_id=1, nOB_mask
        DO i=1, sNx+sNy
          DO fld=1, nOB_fld2D
            subFieldOnMask_2D(ob_id, fld, i) = 0.0
            subFieldOnMask_2Davg(ob_id, fld, i) = 0.0
          ENDDO
          DO fld=1, nOB_fld3D
            DO k=1, Nr
              subFieldOnMask_3D(ob_id, fld, k, i) = 0.0
              subFieldOnMask_3Davg(ob_id, fld, k, i) = 0.0
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      print *,"init_varia:",subFieldOnMask_3D

      print *,"init_varia: fld_nDim",fld_nDim
      print *,"init_varia: fld_to_output",fld_to_output

      print *,"TESTTEST init_varia: loc(subFieldOnMask_3D)",
     &          loc(subFieldOnMask_3D)


      call INI_OB_FILENAMES2()


      DO ob_id=1, nOB_mask
        CALL READ_REC_XY_RL( ob_filenames(ob_id),
     &               ob_subMask(1,:,:,:,:),1, 0, myThid )
        _EXCH_XY_RL(ob_subMask(1,:,:,:,:), myThid)
C       Plotting for debugging
        CALL PLOT_FIELD_XYRL( ob_subMask(1,:,:,:,:),
     &                'init_varia:Open Boundary Submask', -1, myThid )
#ifdef ALLOW_USE_MPI
        CALL IDENTIFY_OB_POINTS(ob_id, myThid)
#endif
      ENDDO

#ifdef ALLOW_USE_MPI
      CALL SEND_CHILD_TO_MASTER()
#endif

C     Using GATHER_2D_R8 to correctly combine sub masks into a single global mask
      DO ob_id=1, nOB_mask
        CALL GATHER_2D_R8(global_ob_mask(ob_id,:,:,1,1),
     &                  ob_subMask(ob_id,1:sNx,1:sNy,1,1),
     &                  Nx, Ny, .true., .true., myThid)
C      CALL PLOT_GLO_FIELD_XYRL(global_ob_mask(ob_id,:,:,:,:),
C     &                    "OB_EXTRACT_OUTPUT: global mask", myThid)
      ENDDO

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
              counter = counter + 1
              if (global_ob_mask(ob_id,i,j,bi,bj) .ge. 1) then
                num_ob_points(ob_id) = num_ob_points(ob_id) + 1
              end if
            end do
           end do
          end do
        end do
      end if

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


      SUBROUTINE INI_OB_FILENAMES2()
C       Fill in array with filenames for each ob mask
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_OB.h"
        INTEGER ob_id
C        character*13 dir
C        character*1 ob_num
C        write(dir, "(A11)") "domain_mask"
C../../../../
C        DO ob_id=1, nOB_mask
C          write(ob_num, "(I1)") ob_id
C          print *,"ob_num",ob_num
C          print *, "ob_filenames(ob_id)",ob_filenames(ob_id)

C        write(ob_filenames(ob_id), "(3A)") dir, ob_num,".bin"
C        print *, "ob_filenames(ob_id)",ob_filenames(ob_id)
C        ENDDO
        DO ob_id=1, nOB_mask
         write(ob_filenames(ob_id),"(A11,I1,A4)") "domain_mask",
     &          ob_id, '.bin'
        ENDDO
C        ob_filenames(1) = "domain_flt32_mask1.bin"
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
      INTEGER counter, num_valid, iG, jG, i, j, bi, bj

      iG = mpi_myXGlobalLo(mpiMyId+1)
      jG = mpi_myYGlobalLo(mpiMyId+1)

      counter = 1
      num_valid = 0
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1,sNy
           DO i=1,sNx
             global_index = int((j-1)*sNx+(j-1)*(Nx-(iG+sNx-1))+
     &                       j*(iG-1)+i)
             if (ob_subMask(ob_id,i,j,bi,bj) .ge. 1.0 ) then
C               print *,"IDENTIFY_OB_POINTS:",
C     &               ob_subMask(ob_id,i,j,bi,bj)
C     First row contains local i's. Second row contains local j's.
               sub_local_ij_ob(ob_id, 1, counter) = i
               sub_local_ij_ob(ob_id, 2, counter) = j
               print *,"IDENTIFY_OB_POINTS: ij:",i,j
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
      print *,"INIT_VARIA: num_valid:", num_valid
      numOBPnts_allproc(ob_id,mpiMyId+1) = num_valid
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
