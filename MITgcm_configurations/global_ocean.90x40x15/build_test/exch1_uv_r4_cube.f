












CBOP
C     !ROUTINE: CPP_EEOPTIONS.h
C     !INTERFACE:
C     include "CPP_EEOPTIONS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP\_EEOPTIONS.h                                         |
C     *==========================================================*
C     | C preprocessor "execution environment" supporting        |
C     | flags. Use this file to set flags controlling the        |
C     | execution environment in which a model runs - as opposed |
C     | to the dynamical problem the model solves.               |
C     | Note: Many options are implemented with both compile time|
C     |       and run-time switches. This allows options to be   |
C     |       removed altogether, made optional at run-time or   |
C     |       to be permanently enabled. This convention helps   |
C     |       with the data-dependence analysis performed by the |
C     |       adjoint model compiler. This data dependency       |
C     |       analysis can be upset by runtime switches that it  |
C     |       is unable to recoginise as being fixed for the     |
C     |       duration of an integration.                        |
C     |       A reasonable way to use these flags is to          |
C     |       set all options as selectable at runtime but then  |
C     |       once an experimental configuration has been        |
C     |       identified, rebuild the code with the appropriate  |
C     |       options set at compile time.                       |
C     *==========================================================*
CEOP


C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C=== Macro related options ===
C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working set size.
C     However, on vector CRAY systems this degrades performance.
C     Enable to switch REAL4_IS_SLOW from genmake2 (with LET_RS_BE_REAL4):

C--   Control use of "double" precision constants.
C     Use D0 where it means REAL*8 but not where it means REAL*16

C--   Enable some old macro conventions for backward compatibility

C=== IO related options ===
C--   Flag used to indicate whether Fortran formatted write
C     and read are threadsafe. On SGI the routines can be thread
C     safe, on Sun it is not possible - if you are unsure then
C     undef this option.

C--   Flag used to indicate whether Binary write to Local file (i.e.,
C     a different file for each tile) and read are thread-safe.

C--   Flag to turn off the writing of error message to ioUnit zero

C--   Alternative formulation of BYTESWAP, faster than
C     compiler flag -byteswapio on the Altix.

C--   Flag to turn on old default of opening scratch files with the
C     STATUS='SCRATCH' option. This method, while perfectly FORTRAN-standard,
C     caused filename conflicts on some multi-node/multi-processor platforms
C     in the past and has been replace by something (hopefully) more robust.

C--   Flag defined for eeboot_minimal.F, eeset_parms.F and open_copy_data_file.F
C     to write STDOUT, STDERR and scratch files from process 0 only.
C WARNING: to use only when absolutely confident that the setup is working
C     since any message (error/warning/print) from any proc <> 0 will be lost.

C=== MPI, EXCH and GLOBAL_SUM related options ===
C--   Flag turns off MPI_SEND ready_to_receive polling in the
C     gather_* subroutines to speed up integrations.

C--   Control MPI based parallel processing
CXXX We no longer select the use of MPI via this file (CPP_EEOPTIONS.h)
CXXX To use MPI, use an appropriate genmake2 options file or use
CXXX genmake2 -mpi .
CXXX #undef  1

C--   Control use of communication that might overlap computation.
C     Under MPI selects/deselects "non-blocking" sends and receives.
C--   Control use of communication that is atomic to computation.
C     Under MPI selects/deselects "blocking" sends and receives.

C--   Control XY periodicity in processor to grid mappings
C     Note: Model code does not need to know whether a domain is
C           periodic because it has overlap regions for every box.
C           Model assume that these values have been
C           filled in some way.

C--   disconnect tiles (no exchange between tiles, just fill-in edges
C     assuming locally periodic subdomain)

C--   Always cumulate tile local-sum in the same order by applying MPI allreduce
C     to array of tiles ; can get slower with large number of tiles (big set-up)

C--   Alternative way of doing global sum without MPI allreduce call
C     but instead, explicit MPI send & recv calls. Expected to be slower.

C--   Alternative way of doing global sum on a single CPU
C     to eliminate tiling-dependent roundoff errors. Note: This is slow.

C=== Other options (to add/remove pieces of code) ===
C--   Flag to turn on checking for errors from all threads and procs
C     (calling S/R STOP_IF_ERROR) before stopping.

C--   Control use of communication with other component:
C     allow to import and export from/to Coupler interface.

C--   Activate some pieces of code for coupling to GEOS AGCM


CBOP
C     !ROUTINE: CPP_EEMACROS.h
C     !INTERFACE:
C     include "CPP_EEMACROS.h"
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP_EEMACROS.h
C     *==========================================================*
C     | C preprocessor "execution environment" supporting
C     | macros. Use this file to define macros for  simplifying
C     | execution environment in which a model runs - as opposed
C     | to the dynamical problem the model solves.
C     *==========================================================*
CEOP


C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C     Flag used to indicate which flavour of multi-threading
C     compiler directives to use. Only set one of these.
C     USE_SOLARIS_THREADING  - Takes directives for SUN Workshop
C                              compiler.
C     USE_KAP_THREADING      - Takes directives for Kuck and
C                              Associates multi-threading compiler
C                              ( used on Digital platforms ).
C     USE_IRIX_THREADING     - Takes directives for SGI MIPS
C                              Pro Fortran compiler.
C     USE_EXEMPLAR_THREADING - Takes directives for HP SPP series
C                              compiler.
C     USE_C90_THREADING      - Takes directives for CRAY/SGI C90
C                              system F90 compiler.






C--   Define the mapping for the _BARRIER macro
C     On some systems low-level hardware support can be accessed through
C     compiler directives here.

C--   Define the mapping for the BEGIN_CRIT() and  END_CRIT() macros.
C     On some systems we simply execute this section only using the
C     master thread i.e. its not really a critical section. We can
C     do this because we do not use critical sections in any critical
C     sections of our code!

C--   Define the mapping for the BEGIN_MASTER_SECTION() and
C     END_MASTER_SECTION() macros. These are generally implemented by
C     simply choosing a particular thread to be "the master" and have
C     it alone execute the BEGIN_MASTER..., END_MASTER.. sections.

CcnhDebugStarts
C      Alternate form to the above macros that increments (decrements) a counter each
C      time a MASTER section is entered (exited). This counter can then be checked in barrier
C      to try and detect calls to BARRIER within single threaded sections.
C      Using these macros requires two changes to Makefile - these changes are written
C      below.
C      1 - add a filter to the CPP command to kill off commented _MASTER lines
C      2 - add a filter to the CPP output the converts the string N EWLINE to an actual newline.
C      The N EWLINE needs to be changes to have no space when this macro and Makefile changes
C      are used. Its in here with a space to stop it getting parsed by the CPP stage in these
C      comments.
C      #define IF ( a .EQ. 1 ) THEN  IF ( a .EQ. 1 ) THEN  N EWLINE      CALL BARRIER_MS(a)
C      #define ENDIF    CALL BARRIER_MU(a) N EWLINE        ENDIF
C      'CPP = cat $< | $(TOOLSDIR)/set64bitConst.sh |  grep -v '^[cC].*_MASTER' | cpp  -traditional -P'
C      .F.f:
C      $(CPP) $(DEFINES) $(INCLUDES) |  sed 's/N EWLINE/\n/' > $@
CcnhDebugEnds

C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working
C     set size. However, on vector CRAY systems this degrades
C     performance.
C- Note: global_sum/max macros were used to switch to  JAM routines (obsolete);
C  in addition, since only the R4 & R8 S/R are coded, GLOBAL RS & RL macros
C  enable to call the corresponding R4 or R8 S/R.



C- Note: a) exch macros were used to switch to  JAM routines (obsolete)
C        b) exch R4 & R8 macros are not practically used ; if needed,
C           will directly call the corrresponding S/R.

C--   Control use of JAM routines for Artic network (no longer supported)
C     These invoke optimized versions of "exchange" and "sum" that
C     utilize the programmable aspect of Artic cards.
CXXX No longer supported ; started to remove JAM routines.
CXXX #ifdef LETS_MAKE_JAM
CXXX #define CALL GLOBAL_SUM_R8 ( a, b) CALL GLOBAL_SUM_R8_JAM ( a, b)
CXXX #define CALL GLOBAL_SUM_R8 ( a, b ) CALL GLOBAL_SUM_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RS ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RL ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RS ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RL ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #endif

C--   Control use of "double" precision constants.
C     Use d0 where it means REAL*8 but not where it means REAL*16

C--   Substitue for 1.D variables
C     Sun compilers do not use 8-byte precision for literals
C     unless .Dnn is specified. CRAY vector machines use 16-byte
C     precision when they see .Dnn which runs very slowly!

C--   Set the format for writing processor IDs, e.g. in S/R eeset_parms
C     and S/R open_copy_data_file. The default of I9.9 should work for
C     a long time (until we will use 10e10 processors and more)



CBOP
C     !ROUTINE: EXCH1_UV_R4_CUBE

C     !INTERFACE:
      SUBROUTINE EXCH1_UV_R4_CUBE(
     U                 Uarray, Varray,
     I                 withSigns,
     I                 myOLw, myOLe, myOLs, myOLn, myNz,
     I                 exchWidthX, exchWidthY,
     I                 cornerMode, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH1_UV_R4_CUBE
C     | o Forward-mode edge exchanges for R4 vector on CS config.
C     *==========================================================*
C     | Controlling routine for exchange of XY edges of an array
C     | distributed in X and Y. The routine interfaces to
C     | communication routines that can use messages passing
C     | exchanges, put type exchanges or get type exchanges.
C     |  This allows anything from MPI to raw memory channel to
C     | memmap segments to be used as a inter-process and/or
C     | inter-thread communiation and synchronisation
C     | mechanism.
C     | Notes --
C     | 1. Some low-level mechanisms such as raw memory-channel
C     | or SGI/CRAY shmem put do not have direct Fortran bindings
C     | and are invoked through C stub routines.
C     | 2. Although this routine is fairly general but it does
C     | require nSx and nSy are the same for all innvocations.
C     | There are many common data structures ( myByLo,
C     | westCommunicationMode, mpiIdW etc... ) tied in with
C     | (nSx,nSy). To support arbitray nSx and nSy would require
C     | general forms of these.
C     | 3. Exchanges on the cube of vector quantities need to be
C     | paired to allow rotations and sign reversal to be applied
C     | consistently between vector components as they rotate.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     *==========================================================*
C     | The design here supports a three-dimensional model grid
C     | with indices I,J and K. The three-dimensional domain
C     | is comprised of nPx*nSx blocks (or tiles) of size sNx
C     | along the first (left-most index) axis, nPy*nSy blocks
C     | of size sNy along the second axis and one block of size
C     | Nr along the vertical (third) axis.
C     | Blocks/tiles have overlap regions of size OLx and OLy
C     | along the dimensions that are subdivided.
C     *==========================================================*
C     \ev
C
C     Voodoo numbers controlling data layout:
C     sNx :: Number of X points in tile.
C     sNy :: Number of Y points in tile.
C     OLx :: Tile overlap extent in X.
C     OLy :: Tile overlap extent in Y.
C     nSx :: Number of tiles per process in X.
C     nSy :: Number of tiles per process in Y.
C     nPx :: Number of processes to use in X.
C     nPy :: Number of processes to use in Y.
C     Nx  :: Number of points in X for the full domain.
C     Ny  :: Number of points in Y for the full domain.
C     Nr  :: Number of points in vertical direction.
CEOP
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr
      PARAMETER (
     &           sNx =  90,
     &           sNy =  40,
     &           OLx =   3,
     &           OLy =   3,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   1,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =  15)

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Controls the sizing of exch
C                routine buffers.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )

CBOP
C     !ROUTINE: EEPARAMS.h
C     !INTERFACE:
C     include "EEPARAMS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | EEPARAMS.h                                               |
C     *==========================================================*
C     | Parameters for "execution environemnt". These are used   |
C     | by both the particular numerical model and the execution |
C     | environment support routines.                            |
C     *==========================================================*
CEOP

C     ========  EESIZE.h  ========================================

C     MAX_LEN_MBUF  :: Default message buffer max. size
C     MAX_LEN_FNAM  :: Default file name max. size
C     MAX_LEN_PREC  :: Default rec len for reading "parameter" files

      INTEGER MAX_LEN_MBUF
      PARAMETER ( MAX_LEN_MBUF = 512 )
      INTEGER MAX_LEN_FNAM
      PARAMETER ( MAX_LEN_FNAM = 512 )
      INTEGER MAX_LEN_PREC
      PARAMETER ( MAX_LEN_PREC = 200 )

C     MAX_NO_THREADS  :: Maximum number of threads allowed.
CC    MAX_NO_PROCS    :: Maximum number of processes allowed.
CC    MAX_NO_BARRIERS :: Maximum number of distinct thread "barriers"
      INTEGER MAX_NO_THREADS
      PARAMETER ( MAX_NO_THREADS =  4 )
c     INTEGER MAX_NO_PROCS
c     PARAMETER ( MAX_NO_PROCS   =  70000 )
c     INTEGER MAX_NO_BARRIERS
c     PARAMETER ( MAX_NO_BARRIERS = 1 )

C     Particularly weird and obscure voodoo numbers
C     lShare :: This wants to be the length in
C               [148]-byte words of the size of
C               the address "window" that is snooped
C               on an SMP bus. By separating elements in
C               the global sum buffer we can avoid generating
C               extraneous invalidate traffic between
C               processors. The length of this window is usually
C               a cache line i.e. small O(64 bytes).
C               The buffer arrays are usually short arrays
C               and are declared REAL ARRA(lShare[148],LBUFF).
C               Setting lShare[148] to 1 is like making these arrays
C               one dimensional.
      INTEGER cacheLineSize
      INTEGER lShare1
      INTEGER lShare4
      INTEGER lShare8
      PARAMETER ( cacheLineSize = 256 )
      PARAMETER ( lShare1 =  cacheLineSize )
      PARAMETER ( lShare4 =  cacheLineSize/4 )
      PARAMETER ( lShare8 =  cacheLineSize/8 )

CC    MAX_VGS  :: Maximum buffer size for Global Vector Sum
c     INTEGER MAX_VGS
c     PARAMETER ( MAX_VGS = 8192 )

C     ========  EESIZE.h  ========================================

C     Symbolic values
C     precXXXX :: precision used for I/O
      INTEGER precFloat32
      PARAMETER ( precFloat32 = 32 )
      INTEGER precFloat64
      PARAMETER ( precFloat64 = 64 )

C     Real-type constant for some frequently used simple number (0,1,2,1/2):
      Real*8     zeroRS, oneRS, twoRS, halfRS
      PARAMETER ( zeroRS = 0.0D0 , oneRS  = 1.0D0 )
      PARAMETER ( twoRS  = 2.0D0 , halfRS = 0.5D0 )
      Real*8     zeroRL, oneRL, twoRL, halfRL
      PARAMETER ( zeroRL = 0.0D0 , oneRL  = 1.0D0 )
      PARAMETER ( twoRL  = 2.0D0 , halfRL = 0.5D0 )

C     UNSET_xxx :: Used to indicate variables that have not been given a value
      Real*8  UNSET_FLOAT8
      PARAMETER ( UNSET_FLOAT8 = 1.234567D5 )
      Real*4  UNSET_FLOAT4
      PARAMETER ( UNSET_FLOAT4 = 1.234567E5 )
      Real*8     UNSET_RL
      PARAMETER ( UNSET_RL     = 1.234567D5 )
      Real*8     UNSET_RS
      PARAMETER ( UNSET_RS     = 1.234567D5 )
      INTEGER UNSET_I
      PARAMETER ( UNSET_I      = 123456789  )

C     debLevX  :: used to decide when to print debug messages
      INTEGER debLevZero
      INTEGER debLevA, debLevB,  debLevC, debLevD, debLevE
      PARAMETER ( debLevZero=0 )
      PARAMETER ( debLevA=1 )
      PARAMETER ( debLevB=2 )
      PARAMETER ( debLevC=3 )
      PARAMETER ( debLevD=4 )
      PARAMETER ( debLevE=5 )

C     SQUEEZE_RIGHT      :: Flag indicating right blank space removal
C                           from text field.
C     SQUEEZE_LEFT       :: Flag indicating left blank space removal
C                           from text field.
C     SQUEEZE_BOTH       :: Flag indicating left and right blank
C                           space removal from text field.
C     PRINT_MAP_XY       :: Flag indicating to plot map as XY slices
C     PRINT_MAP_XZ       :: Flag indicating to plot map as XZ slices
C     PRINT_MAP_YZ       :: Flag indicating to plot map as YZ slices
C     commentCharacter   :: Variable used in column 1 of parameter
C                           files to indicate comments.
C     INDEX_I            :: Variable used to select an index label
C     INDEX_J               for formatted input parameters.
C     INDEX_K
C     INDEX_NONE
      CHARACTER*(*) SQUEEZE_RIGHT
      PARAMETER ( SQUEEZE_RIGHT = 'R' )
      CHARACTER*(*) SQUEEZE_LEFT
      PARAMETER ( SQUEEZE_LEFT = 'L' )
      CHARACTER*(*) SQUEEZE_BOTH
      PARAMETER ( SQUEEZE_BOTH = 'B' )
      CHARACTER*(*) PRINT_MAP_XY
      PARAMETER ( PRINT_MAP_XY = 'XY' )
      CHARACTER*(*) PRINT_MAP_XZ
      PARAMETER ( PRINT_MAP_XZ = 'XZ' )
      CHARACTER*(*) PRINT_MAP_YZ
      PARAMETER ( PRINT_MAP_YZ = 'YZ' )
      CHARACTER*(*) commentCharacter
      PARAMETER ( commentCharacter = '#' )
      INTEGER INDEX_I
      INTEGER INDEX_J
      INTEGER INDEX_K
      INTEGER INDEX_NONE
      PARAMETER ( INDEX_I    = 1,
     &            INDEX_J    = 2,
     &            INDEX_K    = 3,
     &            INDEX_NONE = 4 )

C     EXCH_IGNORE_CORNERS :: Flag to select ignoring or
C     EXCH_UPDATE_CORNERS    updating of corners during an edge exchange.
      INTEGER EXCH_IGNORE_CORNERS
      INTEGER EXCH_UPDATE_CORNERS
      PARAMETER ( EXCH_IGNORE_CORNERS = 0,
     &            EXCH_UPDATE_CORNERS = 1 )

C     FORWARD_SIMULATION
C     REVERSE_SIMULATION
C     TANGENT_SIMULATION
      INTEGER FORWARD_SIMULATION
      INTEGER REVERSE_SIMULATION
      INTEGER TANGENT_SIMULATION
      PARAMETER ( FORWARD_SIMULATION = 0,
     &            REVERSE_SIMULATION = 1,
     &            TANGENT_SIMULATION = 2 )

C--   COMMON /EEPARAMS_L/ Execution environment public logical variables.
C     eeBootError    :: Flags indicating error during multi-processing
C     eeEndError     :: initialisation and termination.
C     fatalError     :: Flag used to indicate that the model is ended with an error
C     debugMode      :: controls printing of debug msg (sequence of S/R calls).
C     useSingleCpuIO :: When useSingleCpuIO is set, MDS_WRITE_FIELD outputs from
C                       master MPI process only. -- NOTE: read from main parameter
C                       file "data" and not set until call to INI_PARMS.
C     useSingleCpuInput :: When useSingleCpuInput is set, EXF_INTERP_READ
C                       reads forcing files from master MPI process only.
C                       -- NOTE: read from main parameter file "data"
C                          and defaults to useSingleCpuInput = useSingleCpuIO
C     printMapIncludesZeros  :: Flag that controls whether character constant
C                               map code ignores exact zero values.
C     useCubedSphereExchange :: use Cubed-Sphere topology domain.
C     useCoupler     :: use Coupler for a multi-components set-up.
C     useNEST_PARENT :: use Parent Nesting interface (pkg/nest_parent)
C     useNEST_CHILD  :: use Child  Nesting interface (pkg/nest_child)
C     useNest2W_parent :: use Parent 2-W Nesting interface (pkg/nest2w_parent)
C     useNest2W_child  :: use Child  2-W Nesting interface (pkg/nest2w_child)
C     useOASIS       :: use OASIS-coupler for a multi-components set-up.
      COMMON /EEPARAMS_L/
c    &  eeBootError, fatalError, eeEndError,
     &  eeBootError, eeEndError, fatalError, debugMode,
     &  useSingleCpuIO, useSingleCpuInput, printMapIncludesZeros,
     &  useCubedSphereExchange, useCoupler,
     &  useNEST_PARENT, useNEST_CHILD,
     &  useNest2W_parent, useNest2W_child, useOASIS,
     &  useSETRLSTK, useSIGREG
      LOGICAL eeBootError
      LOGICAL eeEndError
      LOGICAL fatalError
      LOGICAL debugMode
      LOGICAL useSingleCpuIO
      LOGICAL useSingleCpuInput
      LOGICAL printMapIncludesZeros
      LOGICAL useCubedSphereExchange
      LOGICAL useCoupler
      LOGICAL useNEST_PARENT
      LOGICAL useNEST_CHILD
      LOGICAL useNest2W_parent
      LOGICAL useNest2W_child
      LOGICAL useOASIS
      LOGICAL useSETRLSTK
      LOGICAL useSIGREG

C--   COMMON /EPARAMS_I/ Execution environment public integer variables.
C     errorMessageUnit    :: Fortran IO unit for error messages
C     standardMessageUnit :: Fortran IO unit for informational messages
C     maxLengthPrt1D :: maximum length for printing (to Std-Msg-Unit) 1-D array
C     scrUnit1      :: Scratch file 1 unit number
C     scrUnit2      :: Scratch file 2 unit number
C     eeDataUnit    :: Unit # for reading "execution environment" parameter file
C     modelDataUnit :: Unit number for reading "model" parameter file.
C     numberOfProcs :: Number of processes computing in parallel
C     pidIO         :: Id of process to use for I/O.
C     myBxLo, myBxHi :: Extents of domain in blocks in X and Y
C     myByLo, myByHi :: that each threads is responsble for.
C     myProcId      :: My own "process" id.
C     myPx          :: My X coord on the proc. grid.
C     myPy          :: My Y coord on the proc. grid.
C     myXGlobalLo   :: My bottom-left (south-west) x-index global domain.
C                      The x-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     myYGlobalLo   :: My bottom-left (south-west) y-index in global domain.
C                      The y-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     nThreads      :: No. of threads
C     nTx, nTy      :: No. of threads in X and in Y
C                      This assumes a simple cartesian gridding of the threads
C                      which is not required elsewhere but that makes it easier
C     ioErrorCount  :: IO Error Counter. Set to zero initially and increased
C                      by one every time an IO error occurs.
      COMMON /EEPARAMS_I/
     &  errorMessageUnit, standardMessageUnit, maxLengthPrt1D,
     &  scrUnit1, scrUnit2, eeDataUnit, modelDataUnit,
     &  numberOfProcs, pidIO, myProcId,
     &  myPx, myPy, myXGlobalLo, myYGlobalLo, nThreads,
     &  myBxLo, myBxHi, myByLo, myByHi,
     &  nTx, nTy, ioErrorCount
      INTEGER errorMessageUnit
      INTEGER standardMessageUnit
      INTEGER maxLengthPrt1D
      INTEGER scrUnit1
      INTEGER scrUnit2
      INTEGER eeDataUnit
      INTEGER modelDataUnit
      INTEGER ioErrorCount(MAX_NO_THREADS)
      INTEGER myBxLo(MAX_NO_THREADS)
      INTEGER myBxHi(MAX_NO_THREADS)
      INTEGER myByLo(MAX_NO_THREADS)
      INTEGER myByHi(MAX_NO_THREADS)
      INTEGER myProcId
      INTEGER myPx
      INTEGER myPy
      INTEGER myXGlobalLo
      INTEGER myYGlobalLo
      INTEGER nThreads
      INTEGER nTx
      INTEGER nTy
      INTEGER numberOfProcs
      INTEGER pidIO

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     Uarray      :: (u-type) Array with edges to exchange.
C     Varray      :: (v-type) Array with edges to exchange.
C     withSigns   :: sign of Uarray,Varray depends on orientation
C     myOLw,myOLe :: West  and East  overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C                    Note --
C                    1. In theory one could have a send width and
C                    a receive width for each face of each tile. The only
C                    restriction would be that the send width of one
C                    face should equal the receive width of the sent to
C                    tile face. Dont know if this would be useful. I
C                    have left it out for now as it requires additional
C                    bookeeping.
C     cornerMode  :: Flag indicating whether corner updates are needed.
C     myThid      :: my Thread Id number

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      Real*4     Uarray( 1-myOLw:sNx+myOLe,
     &                1-myOLs:sNy+myOLn,
     &                myNz, nSx, nSy )
      Real*4     Varray( 1-myOLw:sNx+myOLe,
     &                1-myOLs:sNy+myOLn,
     &                myNz, nSx, nSy )
      LOGICAL withSigns
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     theSimulationMode :: Holds working copy of simulation mode
C     theCornerMode     :: Holds working copy of corner mode
C     I,J,K             :: Loop and index counters
C     bl,bt,bn,bs,be,bw :: tile indices
C     negOne, Utmp,Vtmp :: Temps used in swapping and rotating vectors
c     INTEGER theSimulationMode
c     INTEGER theCornerMode
      INTEGER I,J,K, repeat
      INTEGER bl,bt,bn,bs,be,bw
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      Real*4 negOne, Utmp, Vtmp

C     == Statement function ==
C     tilemod :: Permutes indices to return neighboring tile index
C                on six face cube.
      INTEGER tilemod
      tilemod(I)=1+mod(I-1+6,6)
CEOP

c     theSimulationMode = FORWARD_SIMULATION
c     theCornerMode     = cornerMode

c     IF ( simulationMode.EQ.REVERSE_SIMULATION ) THEN
c       WRITE(msgBuf,'(A)')'EXCH1_UV_R4_CUBE: AD mode not implemented'
c       CALL PRINT_ERROR( msgBuf, myThid )
c       STOP 'ABNORMAL END: EXCH1_UV_R4_CUBE: no AD code'
c     ENDIF
      IF ( sNx.NE.sNy .OR.
     &     nSx.NE.6 .OR. nSy.NE.1 .OR.
     &     nPx.NE.1 .OR. nPy.NE.1 ) THEN
        WRITE(msgBuf,'(2A)') 'EXCH1_UV_R4_CUBE: Wrong Tiling'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'EXCH1_UV_R4_CUBE: ',
     &   'works only with sNx=sNy & nSx=6 & nSy=nPx=nPy=1'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: EXCH1_UV_R4_CUBE: Wrong Tiling'
      ENDIF

      negOne = 1.
      IF (withSigns) negOne = -1.

C     For now tile<->tile exchanges are sequentialised through
C     thread 1. This is a temporary feature for preliminary testing until
C     general tile decomposistion is in place (CNH April 11, 2001)
      CALL BAR2( myThid )
      IF ( myThid .EQ. 1 ) THEN

       DO repeat=1,2

       DO bl = 1, 5, 2

        bt = bl
        bn=tilemod(bt+2)
        bs=tilemod(bt-1)
        be=tilemod(bt+1)
        bw=tilemod(bt-2)

        DO K = 1,myNz

C        Tile Odd:Odd+2 [get] [North<-West]
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Uarray(J,sNy+I,K,bt,1) = negOne*Varray(I,sNy+2-J,K,bn,1)
          ENDDO
         ENDDO
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Varray(J,sNy+I,K,bt,1) = Uarray(I,sNy+1-J,K,bn,1)
          ENDDO
         ENDDO
C        Tile Odd:Odd-1 [get] [South<-North]
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Uarray(J,1-I,K,bt,1) = Uarray(J,sNy+1-I,K,bs,1)
          ENDDO
         ENDDO
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Varray(J,1-I,K,bt,1) = Varray(J,sNy+1-I,K,bs,1)
          ENDDO
         ENDDO
C        Tile Odd:Odd+1 [get] [East<-West]
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Uarray(sNx+I,J,K,bt,1) = Uarray(I,J,K,be,1)
          ENDDO
         ENDDO
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Varray(sNx+I,J,K,bt,1) = Varray(I,J,K,be,1)
          ENDDO
         ENDDO
C        Tile Odd:Odd-2 [get] [West<-North]
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Uarray(1-I,J,K,bt,1) = Varray(sNx+1-J,sNy+1-I,K,bw,1)
          ENDDO
         ENDDO
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Varray(1-I,J,K,bt,1) = negOne*Uarray(sNx+2-J,sNy+1-I,K,bw,1)
          ENDDO
         ENDDO

        ENDDO

        bt = bl+1
        bn=tilemod(bt+1)
        bs=tilemod(bt-2)
        be=tilemod(bt+2)
        bw=tilemod(bt-1)

        DO K = 1,myNz

C        Tile Even:Even+1 [get] [North<-South]
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Uarray(J,sNy+I,K,bt,1) = Uarray(J,I,K,bn,1)
          ENDDO
         ENDDO
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Varray(J,sNy+I,K,bt,1) = Varray(J,I,K,bn,1)
          ENDDO
         ENDDO
C        Tile Even:Even-2 [get] [South<-East]
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Uarray(J,1-I,K,bt,1) = negOne*Varray(sNx+1-I,sNy+2-J,K,bs,1)
          ENDDO
         ENDDO
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Varray(J,1-I,K,bt,1) = Uarray(sNx+1-I,sNy+1-J,K,bs,1)
          ENDDO
         ENDDO
C        Tile Even:Even+2 [get] [East<-South]
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Uarray(sNx+I,J,K,bt,1) = Varray(sNx+1-J,I,K,be,1)
          ENDDO
         ENDDO
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Varray(sNx+I,J,K,bt,1) = negOne*Uarray(sNx+2-J,I,K,be,1)
          ENDDO
         ENDDO
C        Tile Even:Even-1 [get] [West<-East]
         DO J = 1,sNy
          DO I = 1,exchWidthX
           Uarray(1-I,J,K,bt,1) = Uarray(sNx+1-I,J,K,bw,1)
          ENDDO
         ENDDO
         DO J = 1,sNy+1
          DO I = 1,exchWidthX
           Varray(1-I,J,K,bt,1) = Varray(sNx+1-I,J,K,bw,1)
          ENDDO
         ENDDO

        ENDDO

       ENDDO

C-     Add one valid uVel,vVel value next to the corner, that allows
C       to compute vorticity on a wider stencil (e.g., vort3(0,1) & (1,0))
       DO bt = 1,6
        DO K = 1,myNz
C      SW corner:
          Uarray(0,0,K,bt,1)=Varray(1,0,K,bt,1)
          Varray(0,0,K,bt,1)=Uarray(0,1,K,bt,1)
C      NW corner:
          Uarray(0,sNy+1,K,bt,1)= negOne*Varray(1,sNy+2,K,bt,1)
          Varray(0,sNy+2,K,bt,1)= negOne*Uarray(0,sNy,K,bt,1)
C      SE corner:
          Uarray(sNx+2,0,K,bt,1)= negOne*Varray(sNx,0,K,bt,1)
          Varray(sNx+1,0,K,bt,1)= negOne*Uarray(sNx+2,1,K,bt,1)
C      NE corner:
          Uarray(sNx+2,sNy+1,K,bt,1)=Varray(sNx,sNy+2,K,bt,1)
          Varray(sNx+1,sNy+2,K,bt,1)=Uarray(sNx+2,sNy,K,bt,1)
        ENDDO
       ENDDO

C      Fix degeneracy at corners
       IF (.FALSE.) THEN
c      IF (withSigns) THEN
        DO bt = 1, 6
         DO K = 1,myNz
C         Top left
          Utmp=0.5*(Uarray(1,sNy,K,bt,1)+Uarray(0,sNy,K,bt,1))
          Vtmp=0.5*(Varray(0,sNy+1,K,bt,1)+Varray(0,sNy,K,bt,1))
          Varray(0,sNx+1,K,bt,1)=(Vtmp-Utmp)*0.70710678
          Utmp=0.5*(Uarray(1,sNy+1,K,bt,1)+Uarray(2,sNy+1,K,bt,1))
          Vtmp=0.5*(Varray(1,sNy+1,K,bt,1)+Varray(1,sNy+2,K,bt,1))
          Uarray(1,sNy+1,K,bt,1)=(Utmp-Vtmp)*0.70710678
C         Bottom right
          Utmp=0.5*(Uarray(sNx+1,1,K,bt,1)+Uarray(sNx+2,1,K,bt,1))
          Vtmp=0.5*(Varray(sNx+1,1,K,bt,1)+Varray(sNx+1,2,K,bt,1))
          Varray(sNx+1,1,K,bt,1)=(Vtmp-Utmp)*0.70710678
          Utmp=0.5*(Uarray(sNx+1,0,K,bt,1)+Uarray(sNx,0,K,bt,1))
          Vtmp=0.5*(Varray(sNx,1,K,bt,1)+Varray(sNx,0,K,bt,1))
          Uarray(sNx+1,0,K,bt,1)=(Utmp-Vtmp)*0.70710678
C         Bottom left
          Utmp=0.5*(Uarray(1,1,K,bt,1)+Uarray(0,1,K,bt,1))
          Vtmp=0.5*(Varray(0,1,K,bt,1)+Varray(0,2,K,bt,1))
          Varray(0,1,K,bt,1)=(Vtmp+Utmp)*0.70710678
          Utmp=0.5*(Uarray(1,0,K,bt,1)+Uarray(2,0,K,bt,1))
          Vtmp=0.5*(Varray(1,1,K,bt,1)+Varray(1,0,K,bt,1))
          Uarray(1,0,K,bt,1)=(Utmp+Vtmp)*0.70710678
C         Top right
          Utmp=0.5*(Uarray(sNx+1,sNy,K,bt,1)+Uarray(sNx+2,sNy,K,bt,1))
          Vtmp=0.5*(Varray(sNx+1,sNy+1,K,bt,1)+Varray(sNx+1,sNy,K,bt,1))
          Varray(sNx+1,sNy+1,K,bt,1)=(Vtmp+Utmp)*0.70710678
          Utmp=0.5*(Uarray(sNx+1,sNy+1,K,bt,1)+Uarray(sNx,sNy+1,K,bt,1))
          Vtmp=0.5*(Varray(sNx,sNy+1,K,bt,1)+Varray(sNx,sNy+2,K,bt,1))
          Uarray(sNx+1,sNy+1,K,bt,1)=(Utmp+Vtmp)*0.70710678
         ENDDO
        ENDDO
       ENDIF

       ENDDO

      ENDIF
      CALL BAR2(myThid)

      RETURN
      END
