#include "DIAG_OPTIONS.h"


       SUBROUTINE WRITE_GLOBAL_BIN(total_valid,fld,ob_id,myIter,nDim)
C-------------------------------------------------------------------------------
c      This program writes the final output (global final array) to a binary file.

C      Writes either a 32 bit precision or 64 bit precision file.
C      Two ways to write output:
C         1. Write all times levels into a single binary files
C                 - Each time level is written into its own record in the binary file
C         2. Write a separates binary file for each time level
C                 - File names consist of iter number: MASK_NN_FIELD_ITER.bin
C
C      input params:
C      total_valid - total number of open boundary points to write to file
C      fld - field number
C      ob_id - open boundary number
C      myIter - current iter number
C      nDim - number of dimensions in field (either 2D or 3D)
C-------------------------------------------------------------------------------


#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"
#include "DIAGNOSTICS_OB.h"

       integer total_valid, ob_id, fld, nDim, myIter

       character*2 dir
       character*50 filename
       real*4 bufferR4_2D((sNy+sNx)*(nPx*nPy))
       real*4 bufferR4_3D((sNy+sNx)*(nPx*nPy), Nr)
       real*8 temp((sNy+sNx)*(nPx*nPy), Nr)
       integer, parameter :: debug = 1
       logical exists
       integer LEN, i

       LEN = 1
       dir = "./"

       IF (ob_filePrec .eq. 32) then

        IF (nDim .eq. 2) then
C       --------------------------------
C       File precision: 32
C       Field dimension: 2
C       --------------------------------
         call OB_R8_to_R4(bufferR4_2D, bufferR4_3D,2,fld,ob_id)
         call CREATE_FILENAME(dir, ob_flds2D(fld,ob_id), ob_id,filename
     &                      ,0, myIter)
         IF (combineMaskTimeLevels) then
C         Writing all time levels into a single binary file
          INQUIRE (IOLENGTH = LEN) bufferR4_2D(1:total_valid)
          INQUIRE (FILE=filename, EXIST=exists)
C         If file already exists, set status as 'old'. Otherwise, 'replace'
          if (exists) then
           open(1, FILE=filename, STATUS='OLD', ACCESS='DIRECT',
     &          RECL = LEN )
          else
           open(1, FILE=filename, STATUS='REPLACE', ACCESS='DIRECT',
     &          RECL = LEN )
          endif
          write(1, REC=time_level) bufferR4_2D(1:total_valid)
         ELSE
C         Write each time level into its own binary file
          open(1, FILE=filename, FORM="unformatted",ACCESS='STREAM')
          write(1) bufferR4_2D(1:total_valid)
         ENDIF
         close(1)


        ELSE IF (nDim .eq. 3) then
C       --------------------------------
C       File precision: 32
C       Field dimension: 3
C       --------------------------------
         call OB_R8_to_R4(bufferR4_2D, bufferR4_3D,3,fld,ob_id)
         call CREATE_FILENAME(dir, ob_flds3D(fld,ob_id), ob_id, filename
     &                      ,0, myIter)
         IF (combineMaskTimeLevels) then
C        Writing all time levels into a single binary file
          INQUIRE (IOLENGTH = LEN) bufferR4_3D(1:total_valid,
     &          1:ob_levels3D(fld,ob_id))
          INQUIRE (FILE=filename, EXIST=exists)
C         If file already exists, set status as 'old'. Otherwise, 'replace'
          if (exists) then
           open( 1, FILE=filename, STATUS='OLD', ACCESS='DIRECT',
     &          RECL = LEN)
          else
           open( 1, FILE=filename, STATUS='REPLACE', ACCESS='DIRECT',
     &          RECL =  LEN)
          endif
          write(1, REC=time_level) bufferR4_3D(1:total_valid,
     &          1:ob_levels3D(fld,ob_id))
         ELSE
C         Write each time level into its own binary file
          open(1, FILE=filename, FORM="unformatted",ACCESS='STREAM')
          write(1) bufferR4_3D(1:total_valid,1:ob_levels3D(fld,ob_id))
         ENDIF
         close(1)

        ENDIF


       ELSE IF (ob_filePrec .eq. 64) then

        IF (nDim .eq. 2) then
C       --------------------------------
C       File precision: 64
C       Field dimension: 2
C       --------------------------------
         call CREATE_FILENAME(dir, ob_flds2D(fld,ob_id), ob_id,filename
     &                      ,0, myIter)
         IF (combineMaskTimeLevels) then
          INQUIRE (IOLENGTH = LEN) global_ob2D(1:total_valid)
          INQUIRE (FILE=filename, EXIST=exists)
C         If file already exists, set status as 'old'. Otherwise, 'replace'
          if (exists) then
           open( 1, FILE=filename, STATUS='OLD', ACCESS='DIRECT',
     &          RECL = LEN )
          else
           open( 1, FILE=filename, STATUS='REPLACE', ACCESS='DIRECT',
     &          RECL = LEN )
          endif

          write(1, REC=time_level) global_ob2D(1:total_valid)
         ELSE
          open(1, FILE=filename, FORM="unformatted",ACCESS='STREAM')
          write(1) global_ob2D(1:total_valid)
         ENDIF
         close(1)


        ELSE IF (nDim .eq. 3) then
C       --------------------------------
C       File precision: 64
C       Field dimension: 3
C       --------------------------------
         call CREATE_FILENAME(dir, ob_flds3D(fld,ob_id), ob_id, filename
     &                      ,0, myIter)

         IF (combineMaskTimeLevels) then
          INQUIRE (IOLENGTH = LEN) global_ob3D(1:total_valid,
     &                                   1:ob_levels3D(fld,ob_id))
          INQUIRE (FILE=filename, EXIST=exists)
C         If file already exists, set status as 'old'. Otherwise, 'replace'
          if (exists) then
           open( 1, FILE=filename, STATUS='OLD', ACCESS='DIRECT',
     &          RECL = LEN)
          else
           open( 1, FILE=filename, STATUS='REPLACE', ACCESS='DIRECT',
     &          RECL = LEN)
          endif
          write(1, REC=time_level) global_ob3D(1:total_valid,
     &                                      1:ob_levels3D(fld,ob_id))
          if (debug .eq. 1) then
C          Reading back from binary file for debugging
           do i=1, time_level
            read(1, REC=i) temp(1:total_valid,
     &                                    1:ob_levels3D(fld,ob_id))
            print*,"time level:",i
            print *,"TEMP: ",temp(1:total_valid,
     &                                    1:ob_levels3D(fld,ob_id))
           enddo
          endif
         ELSE
          open(1, FILE=filename, FORM="unformatted",ACCESS='STREAM')
          write(1) global_ob3D(1:total_valid,1:ob_levels3D(fld,ob_id))
         ENDIF

         close(1)

        ENDIF

       ENDIF

      END



       SUBROUTINE CREATE_FILENAME(dir, fldNm, ob_id, filename,
     &                            is_child, myIter)
c       Creating filenames for all binary files.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

       integer is_child, ob_id, myIter, i
       _RL myTime
       character*50 filename
       character*8 fldNm
       character*2 dir
       character*3 str_px, str_py
       character*11 str_time

       character*2 str_ob
       character*8 str_iter

       integer len, count


c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       001 for index 1 in x axis
c       002 for index 2 in x axis
       write(str_ob, "(I2.2)") ob_id
       write(str_iter, "(I8.8)") myIter
       write(str_time, "(F11.0)") myTime

C       Counting number of valid characters to trim trailing whitespaces
        count = 0
        DO i=1, len(fldNm)
          IF (fldNm(i:i) .ne. ' ') then
            count = count + 1
          ENDIF
        ENDDO

       if (is_child .eq. 1) then
         write(str_px, "(I3.3)") mpiPx
         write(str_py, "(I3.3)") mpiPy
         write(filename, "(8A)") dir, fldNm(1:count),
     &                   "_avg.",str_iter,str_px, ".", str_py, ".bin"
       else
c         Creating a filename for global arrays
         IF (combineMaskTimeLevels) then
          write(filename, "(5A)") "MASK_", str_ob, "_",
     &                            fldNm(1:count),".bin"
         ELSE
           write(filename, "(7A)") "MASK_", str_ob, "_",
     &                   fldNm(1:count),"_",str_iter,".bin"
         ENDIF

       end if

       END




      SUBROUTINE OB_R8_to_R4(bufferR4_2D, bufferR4_3D,nDim,fld,ob_id)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "DIAGNOSTICS_OB.h"

        real*4 bufferR4_2D((sNy+sNx)*(nPx*nPy))
        real*4 bufferR4_3D((sNy+sNx)*(nPx*nPy), Nr)
        integer nDim, fld, ob_id
        integer, parameter :: debug = 0
        integer i, k

C        Initialize buffer arrays:
       do i=1, (sNy+sNx)*(nPx*nPy)
         bufferR4_2D(i) = 0.0
         do k=1, Nr
           bufferR4_3D(i,k) = 0.0
         enddo
       enddo

C        Set real*8 values to real*4

        if (nDim .eq. 2) then
          do i=1, (sNy+sNx)*(nPx*nPy)
            bufferR4_2D(i) = global_ob2D(i)
          enddo

        else if (nDim .eq. 3) then
          do i=1, (sNy+sNx)*(nPx*nPy)
            do k=1, ob_levels3D(fld,ob_id)
              bufferR4_3D(i,k) = global_ob3D(i,k)
            enddo
          enddo

          if (debug .eq. 1) then
           do k=1,ob_levels3D(fld,ob_id)
            print *,"OB: bufferR4_3D:",bufferR4_3D(:,k)
            print *,"OB: global_ob3D:",global_ob3D(:,k)
           enddo
          endif

        endif

      END
