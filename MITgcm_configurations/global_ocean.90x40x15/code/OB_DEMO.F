#include "DIAG_OPTIONS.h"


       SUBROUTINE OB_DEMO(myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "DIAGNOSTICS_OB.h"

      integer myThid
      integer sz

      sz = (abs(1-OLx)+(sNx+OLx))*
     &                   (abs(1-OLy)+(sNy+OLy))*nSx*nSy

      print *, "OB_DEMO: Printing subField"
      CALL PRINT_FLOAT_ARR(subField(1,:,:,:,:),sz)

      END
