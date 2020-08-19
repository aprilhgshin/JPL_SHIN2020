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

C      print *, "OB_DEMO: Printing subField"
C      CALL PRINT_FLOAT_ARR(subField(1,:,:,:,:),sz)

      CALL PLOT_FIELD_XYRL( ob_subMask(1,:,:,:,:),
     &                'OB_DEMO:Open Boundary Submask', -1, myThid )



      END
