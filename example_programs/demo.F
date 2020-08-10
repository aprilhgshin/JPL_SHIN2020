
      program DEMO
#include "DEMO.h"
        character*4 suffix
        REAL arr(4,4)
        INTEGER arr100(5)
        INTEGER sh(3), D(5)
        character*103 filename(2)
        character*98 dir
        character*1 ob_num
        write(dir, "(3A)") "/home/mitgcm/Work/JPL_SHIN2020",
     &    "/MITgcm_configurations/",
     &    "global_ocean.90x40x15/input/domain_flt32_mask"
        write(suffix, "(A4)") ".bin"
        write(ob_num, "(I1)") 1

        write(filename(2), "(3A)") dir, ob_num, suffix
        print *, filename(2)
        print *, sh(2)
C        print *, "header_arr before",header_arr
C        call DD()
C        print *, "header_arr after",header_arr

C        print *, "function D(100):",D(100)
C        print *, "function D(1):",D(1)

C        print *, "arr",arr100


      END

      FUNCTION D()
#include "DEMO.h"
C        POINTER (D, arr2)
        INTEGER fld, D(5)
C        arr2(1) = 0
C        arr2(2) = 0
C        if (fld .eq. 100) then
C          arr(1) = 1
C          arr(2) = 2
C          D = loc(arr)
          DO i=1,5
            header_arr(i) = i
          ENDDO
C        endif

      RETURN
      END

      SUBROUTINE DD()
#include "DEMO.h"
C        INTEGER local_arr(5), i
C        POINTER (ptr, header_arr)

C        DO i=1,5
C          local_arr(i) = i
C        ENDDO
C        ptr = loc(local_arr)
C        print *,"subroutine DD:",header_arr
      END
