












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


C CPP options file for EXCH2 package


C ... W2_USE_E2_SAFEMODE description ...

C Debug mode option:

C Use only exch2_R1_cube (and avoid calling exch2_R2_cube)

C Fill null regions (face-corner halo regions) with e2FillValue_RX (=0)
C notes: for testing (allow to check that results are not affected)

C Process Global Cumulated-Sum using a Tile x Tile (x 2) Matrix
C notes: should be faster (vectorise) but storage of this matrix might
C        become an issue on large set-up (with many tiles)


CBOP 0
C !ROUTINE: EXCH2_RECV_R41

C !INTERFACE:
      SUBROUTINE EXCH2_RECV_R41(
     I       thisTile, nN,
     I       e2BufrRecSize,
     I       iBufr,
     O       e2Bufr1_R4,
     I       commSetting, myThid )

C !DESCRIPTION:
C     Scalar field (1 component) Exchange:
C     Receive into buffer exchanged data from the source Process.
C     buffer data will be used to fill in the tile-edge overlap region.

C !USES:
      IMPLICIT NONE

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
CBOP
C     !ROUTINE: EESUPPORT.h
C     !INTERFACE:
C     include "EESUPPORT.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | EESUPPORT.h                                              |
C     *==========================================================*
C     | Support data structures for the MITgcm UV ``execution    |
C     | environment'' code. This data should be private to the   |
C     | execution environment routines. Data which needs to be   |
C     | accessed directly by a numerical model goes in           |
C     | EEPARAMS.h.                                              |
C     *==========================================================*
CEOP

C     ERROR_HEADER        - String which prefixes error messages
      CHARACTER*(*) ERROR_HEADER
      PARAMETER ( ERROR_HEADER = ' *** ERROR ***' )
C     PROCESS_HEADER      - String which prefixes processor number
      CHARACTER*(*) PROCESS_HEADER
      PARAMETER ( PROCESS_HEADER = 'PID.TID' )

C     MAX_NUM_COMM_MODES - Maximum number of communication modes
C     COMM_NONE       - No edge communication
C     COMM_MSG        - Use messages to communicate edges
C     COMM_PUT        - Use put to communicate edges
C     COMM_GET        - Use get to communicate edges
C     Note - commName holds an identifying name for each communication
C            mode. The COMM_ parameters are used to index commName
C            so the COMM_ parameters need to be in the range
C            1 : MAX_NUM_COMM_MODES.
      INTEGER MAX_NUM_COMM_MODES
      PARAMETER ( MAX_NUM_COMM_MODES = 4 )
      INTEGER COMM_NONE
      PARAMETER ( COMM_NONE   =   1 )
      INTEGER COMM_MSG
      PARAMETER ( COMM_MSG    =   2 )
      INTEGER COMM_PUT
      PARAMETER ( COMM_PUT    =   3 )
      INTEGER COMM_GET
      PARAMETER ( COMM_GET    =   4 )
      COMMON /EESUPP_COMMNAME/ commName
      CHARACTER*10 commName(MAX_NUM_COMM_MODES)

C     Tile identifiers
C     Tiles have a number that is unique over the global domain.
C     A tile that is not there has its number set to NULL_TILE
      INTEGER NULL_TILE
      PARAMETER ( NULL_TILE = -1 )


C--   COMMON /EESUPP_C/ Execution environment support character variables
C     myProcessStr - String identifying my process number
      COMMON /EESUPP_C/ myProcessStr
      CHARACTER*128 myProcessStr

C--   COMMON /EESUPP_L/ Execution environment support logical variables
C     initMPError - Flag indicating error during multi-processing
C                   initialisation.
C     finMPError  - Flag indicating error during multi-processing
C                   termination.
C     ThError     - Thread detected an error.
C     usingMPI    - Flag controlling use of MPI routines. This flag
C                   allows either MPI or threads to be used in a
C                   shared memory environment which can be a useful
C                   debugging/performance analysis tool.
C     usingSyncMessages - Flag that causes blocking communication to be used
C                         if possible. When false non-blocking EXCH routines
C                         will be used if possible.
C     notUsingXPeriodicity - Flag indicating no X/Y boundary wrap around
C     notUsingYPeriodicity   This affects the communication routines but
C                            is generally ignored in the numerical model
C                            code.
C     threadIsRunning, threadIsComplete - Flags used to check for correct behaviour
C                                         of multi-threaded code.
C                                         threadIsRunning is used to check that the
C                                         threads we need are running. This catches the
C                                         situation where a program eedata file has nTthreads
C                                         greater than the setenv PARALLEL or NCPUS variable.
C                                         threadIsComplete is used to flag that a thread has
C                                         reached the end of the model. This is used as a check to
C                                         trap problems that might occur if one thread "escapes"
C                                         the main.F master loop. This should not happen
C                                         if the multi-threading compilation tools works right.
C                                         But (see for example KAP) this is not always the case!
      COMMON /EESUPP_L/ thError, threadIsRunning, threadIsComplete,
     & allMyEdgesAreSharedMemory, usingMPI, usingSyncMessages,
     & notUsingXPeriodicity, notUsingYPeriodicity
      LOGICAL thError(MAX_NO_THREADS)
      LOGICAL threadIsRunning(MAX_NO_THREADS)
      LOGICAL threadIsComplete(MAX_NO_THREADS)
      LOGICAL allMyEdgesAreSharedMemory(MAX_NO_THREADS)
      LOGICAL usingMPI
      LOGICAL usingSyncMessages
      LOGICAL notUsingXPeriodicity
      LOGICAL notUsingYPeriodicity

C--   COMMON /EESUPP_I/ Parallel support integer globals
C     pidW   -  Process  ID of neighbor to West
C     pidE   -           ditto             East
C     pidN   -           ditto             North
C     pidS   -           ditto             South
C              Note: pid[XY] is not necessairily the UNIX
C                    process id - it is just an identifying
C                    number.
C     myPid  - My own process id
C     nProcs - Number of processes
C     westCommunicationMode  - Mode of communication for each tile face
C     eastCommunicationMode
C     northCommunicationMode
C     southCommunicationMode
C     bi0   - Low cartesian tile index for this process
C     bj0     Note - In a tile distribution with holes bi0 and bj0
C                    are not useful. Neighboring tile indices must
C                    be derived some other way.
C     tileNo       - Tile identification number for my tile and
C     tileNo[WENS]   my N,S,E,W neighbor tiles.
C     tilePid[WENS] - Process identification number for
C                     my N,S,E,W neighbor tiles.
C     nTx, nTy    - No. threads in X and Y. This assumes a simple
C                   cartesian gridding of the threads which is not
C                   required elsewhere but that makes it easier.
      COMMON /EESUPP_I/
     & myPid, nProcs, pidW, pidE, pidN, pidS,
     & tileCommModeW,  tileCommModeE,
     & tileCommModeN,  tileCommModeS,
     & tileNo, tileNoW, tileNoE, tileNoS, tileNoN,
     &  tilePidW, tilePidE, tilePidS, tilePidN,
     &  tileBiW, tileBiE, tileBiS, tileBiN,
     & tileBjW, tileBjE, tileBjS, tileBjN,
     & tileTagSendW, tileTagSendE, tileTagSendS, tileTagSendN,
     & tileTagRecvW, tileTagRecvE, tileTagRecvS, tileTagRecvN
      INTEGER myPid
      INTEGER nProcs
      INTEGER pidW
      INTEGER pidE
      INTEGER pidN
      INTEGER pidS
      INTEGER tileCommModeW ( nSx, nSy )
      INTEGER tileCommModeE ( nSx, nSy )
      INTEGER tileCommModeN ( nSx, nSy )
      INTEGER tileCommModeS ( nSx, nSy )
      INTEGER tileNo( nSx, nSy )
      INTEGER tileNoW( nSx, nSy )
      INTEGER tileNoE( nSx, nSy )
      INTEGER tileNoN( nSx, nSy )
      INTEGER tileNoS( nSx, nSy )
      INTEGER tilePidW( nSx, nSy )
      INTEGER tilePidE( nSx, nSy )
      INTEGER tilePidN( nSx, nSy )
      INTEGER tilePidS( nSx, nSy )
      INTEGER tileBiW( nSx, nSy )
      INTEGER tileBiE( nSx, nSy )
      INTEGER tileBiN( nSx, nSy )
      INTEGER tileBiS( nSx, nSy )
      INTEGER tileBjW( nSx, nSy )
      INTEGER tileBjE( nSx, nSy )
      INTEGER tileBjN( nSx, nSy )
      INTEGER tileBjS( nSx, nSy )
      INTEGER tileTagSendW( nSx, nSy )
      INTEGER tileTagSendE( nSx, nSy )
      INTEGER tileTagSendN( nSx, nSy )
      INTEGER tileTagSendS( nSx, nSy )
      INTEGER tileTagRecvW( nSx, nSy )
      INTEGER tileTagRecvE( nSx, nSy )
      INTEGER tileTagRecvN( nSx, nSy )
      INTEGER tileTagRecvS( nSx, nSy )

C--   Include MPI standard Fortran header file
! -*- fortran -*-
!
! Copyright (c) 2004-2006 The Trustees of Indiana University and Indiana
!                         University Research and Technology
!                         Corporation.  All rights reserved.
! Copyright (c) 2004-2005 The University of Tennessee and The University
!                         of Tennessee Research Foundation.  All rights
!                         reserved.
! Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
!                         University of Stuttgart.  All rights reserved.
! Copyright (c) 2004-2005 The Regents of the University of California.
!                         All rights reserved.
! Copyright (c) 2006-2014 Cisco Systems, Inc.  All rights reserved.
! $COPYRIGHT$
!
! Additional copyrights may follow
!
! $HEADER$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do ***not*** copy this file to the directory where your Fortran
! fortran application is compiled unless it is absolutely necessary!  Most
! modern Fortran compilers now support the -I command line flag, which
! tells the compiler where to find .h files (specifically, this one).  For
! example:
!
!      shell$ mpifort foo.f -o foo -I$OMPI_HOME/include
!
! will probably do the trick (assuming that you have set OMPI_HOME
! properly).
!
! That being said, OMPI's "mpifort" wrapper compiler should
! automatically include the -I option for you.  The following command
! should be equivalent to the command listed above:
!
!      shell$ mpifort foo.f -o foo
!
! You should not copy this file to your local directory because it is
! possible that this file will be changed between versions of Open MPI.
! Indeed, this mpif.h is incompatible with the mpif.f of other
! implementations of MPI.  Using this mpif.h with other implementations
! of MPI, or with other versions of Open MPI will result in undefined
! behavior (to include incorrect results, segmentation faults,
! unexplainable "hanging" in your application, etc.).  Always use the
! -I command line option instead (or let mpifort do it for you).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      include 'mpif-config.h'
      include 'mpif-constants.h'
      include 'mpif-handles.h'
      include "mpif-io-constants.h"
      include "mpif-io-handles.h"
      include 'mpif-externals.h'
      include 'mpif-sentinels.h'
      include 'mpif-sizeof.h'

C--   COMMON /EESUPP_MPI_I/ MPI parallel support integer globals
C     mpiPidW   - MPI process id for west neighbor.
C     mpiPidE   - MPI process id for east neighbor.
C     mpiPidN   - MPI process id for north neighbor.
C     mpiPidS   - MPI process id for south neighbor.
C     mpiPidNW  - MPI process id for northwest neighbor.
C     mpiPidNE  - MPI process id for northeast neighbor.
C     mpiPidSW  - MPI process id for southwest neighbor.
C     mpiPidSE  - MPI process id for southeast neighbor.
C     mpiPidIO  - MPI process to use for IO.
C     mpiNprocs - No. of MPI processes.
C     mpiMyId   - MPI process id of me.
C     mpiComm   - MPI communicator to use.
C     mpiPx     - My MPI proc. grid X coord
C     mpiPy     - My MPI proc. grid Y coord
C     mpiXGlobalLo - My bottom-left (south-west) x-coordinate in
C                    global domain.
C     mpiYGlobalLo - My bottom-left (south-west) y-coordinate in
C                    global domain.
C     mpiTypeXFaceBlock_xy_r4  - Primitives for communicating edge
C     mpiTypeXFaceBlock_xy_r8    of a block.
C     mpiTypeYFaceBlock_xy_r4    XFace is used in east-west transfer
C     mpiTypeYFaceBlock_xy_r8    YFace is used in nrth-south transfer
C     mpiTypeXFaceBlock_xyz_r4   xy is used in two-dimensional arrays
C     mpiTypeXFaceBlock_xyz_r8   xyz is used with three-dimensional arrays
C     mpiTypeYFaceBlock_xyz_r4   r4 is used for real*4 data
C     mpiTypeYFaceBlock_xyz_r8   r8 is used for real*8 data
C     mpiTypeXFaceThread_xy_r4  - Composites of the above primitives
C     mpiTypeXFaceThread_xy_r8    for communicating edges of all blocks
C     mpiTypeYFaceThread_xy_r4    owned by a thread.
C     mpiTypeYFaceThread_xy_r8
C     mpiTypeXFaceThread_xyz_r4
C     mpiTypeXFaceThread_xyz_r8
C     mpiTypeYFaceThread_xyz_r4
C     mpiTypeYFaceBlock_xyz_r8
C     mpiTagE       - Tags are needed to mark requests when MPI is running
C     mpiTagW         between multithreaded processes or when the same process.
C     mpiTagS         is a neighbor in more than one direction. The tags ensure that
C     mpiTagN         a thread will get the message it is looking for.
C     mpiTagSW        The scheme adopted is to tag messages according to
C     mpiTagSE        the direction they are travelling. Thus a message
C     mpiTagNW        travelling east is tagged mpiTagE. However, in a
C     mpiTagNE        multi-threaded environemnt several messages could
C                     be travelling east from the same process at the
C                     same time. The tag is therefore modified to
C                     be mpiTag[EWS...]*nThreads+myThid. This requires that
C                     each thread also know the thread ids of its "neighbor"
C                     threads.
      COMMON /EESUPP_MPI_I/
     & mpiPidW,  mpiPidE,  mpiPidS,  mpiPidN,
     & mpiPidSE, mpiPidSW, mpiPidNE, mpiPidNW,
     & mpiPidIo, mpiMyId, mpiNProcs, mpiComm,
     & mpiPx, mpiPy, mpiXGlobalLo, mpiYGlobalLo,
     & mpiTypeXFaceBlock_xy_r4, mpiTypeXFaceBlock_xy_r8,
     & mpiTypeYFaceBlock_xy_r4, mpiTypeYFaceBlock_xy_r8,
     & mpiTypeXFaceBlock_xyz_r4, mpiTypeXFaceBlock_xyz_r8,
     & mpiTypeYFaceBlock_xyz_r4, mpiTypeYFaceBlock_xyz_r8,
     & mpiTypeXFaceThread_xy_r4, mpiTypeXFaceThread_xy_r8,
     & mpiTypeYFaceThread_xy_r4, mpiTypeYFaceThread_xy_r8,
     & mpiTypeXFaceThread_xyz_r4, mpiTypeXFaceThread_xyz_r8,
     & mpiTypeYFaceThread_xyz_r4, mpiTypeYFaceThread_xyz_r8,
     & mpiTagE, mpiTagW, mpiTagN, mpiTagS,
     & mpiTagSE, mpiTagSW, mpiTagNW, mpiTagNE

      INTEGER mpiPidW
      INTEGER mpiPidE
      INTEGER mpiPidS
      INTEGER mpiPidN
      INTEGER mpiPidSW
      INTEGER mpiPidSE
      INTEGER mpiPidNW
      INTEGER mpiPidNE
      INTEGER mpiPidIO
      INTEGER mpiMyId
      INTEGER mpiNProcs
      INTEGER mpiComm
      INTEGER mpiPx
      INTEGER mpiPy
      INTEGER mpiXGlobalLo
      INTEGER mpiYGlobalLo
      INTEGER mpiTypeXFaceBlock_xy_r4
      INTEGER mpiTypeXFaceBlock_xy_r8
      INTEGER mpiTypeYFaceBlock_xy_r4
      INTEGER mpiTypeYFaceBlock_xy_r8
      INTEGER mpiTypeXFaceBlock_xyz_r4
      INTEGER mpiTypeXFaceBlock_xyz_r8
      INTEGER mpiTypeYFaceBlock_xyz_r4
      INTEGER mpiTypeYFaceBlock_xyz_r8
      INTEGER mpiTypeXFaceThread_xy_r4(MAX_NO_THREADS)
      INTEGER mpiTypeXFaceThread_xy_r8(MAX_NO_THREADS)
      INTEGER mpiTypeYFaceThread_xy_r4(MAX_NO_THREADS)
      INTEGER mpiTypeYFaceThread_xy_r8(MAX_NO_THREADS)
      INTEGER mpiTypeXFaceThread_xyz_r4(MAX_NO_THREADS)
      INTEGER mpiTypeXFaceThread_xyz_r8(MAX_NO_THREADS)
      INTEGER mpiTypeYFaceThread_xyz_r4(MAX_NO_THREADS)
      INTEGER mpiTypeYFaceThread_xyz_r8(MAX_NO_THREADS)
      INTEGER mpiTagNW
      INTEGER mpiTagNE
      INTEGER mpiTagSW
      INTEGER mpiTagSE
      INTEGER mpiTagW
      INTEGER mpiTagE
      INTEGER mpiTagN
      INTEGER mpiTagS

C--   COMMON /MPI_FULLMAP_I/ holds integer arrays of the full list of MPI process
C     mpi_myXGlobalLo :: List of all processors bottom-left X-index in global domain
C     mpi_myYGlobalLo :: List of all processors bottom-left Y-index in global domain
C                        Note: needed for mpi gather/scatter routines & singleCpuIO.
      COMMON /MPI_FULLMAP_I/
     &        mpi_myXGlobalLo, mpi_myYGlobalLo
      INTEGER mpi_myXGlobalLo(nPx*nPy)
      INTEGER mpi_myYGlobalLo(nPx*nPy)

C MPI communicator describing this model realization
      COMMON /MPI_COMMS/
     &        MPI_COMM_MODEL
      INTEGER MPI_COMM_MODEL

CBOP
C    !ROUTINE: W2_EXCH2_SIZE.h
C    !INTERFACE:
C    include W2_EXCH2_SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | W2_EXCH2_SIZE.h
C     | Declare size of Wrapper2-Exch2 arrays
C     *==========================================================*
C     | Expected to be modified for unconventional configuration
C     | (e.g., many blank-tiles) or specific topology.
C     *==========================================================*
CEOP

C---   Size of Tiling topology structures
C  W2_maxNbFacets   :: Maximum number of Facets (also and formerly called
C                   :: "domains" or "sub-domains") of this topology.
C  W2_maxNeighbours :: Maximum number of neighbours any tile has.
C  W2_maxNbTiles    :: Maximum number of tiles (active+blank) in this topology
C  W2_ioBufferSize  :: Maximum size of Single-CPU IO buffer.
       INTEGER W2_maxNbFacets
       INTEGER W2_maxNeighbours
       INTEGER W2_maxNbTiles
       INTEGER W2_ioBufferSize
       INTEGER W2_maxXStackNx
       INTEGER W2_maxXStackNy
       INTEGER W2_maxYStackNx
       INTEGER W2_maxYStackNy

C---   Default values :
C      (suitable for 6-face Cube-Sphere topology, compact global I/O format)
C      W2_maxNbTiles = Nb of active tiles (=nSx*nSy*nPx*nPy) + Max_Nb_BlankTiles
C      default assume a large Max_Nb_BlankTiles equal to Nb of active tiles
C      resulting in doubling the tile number.
       PARAMETER ( W2_maxNbFacets = 10 )
       PARAMETER ( W2_maxNeighbours = 8 )
       PARAMETER ( W2_maxNbTiles = nSx*nSy*nPx*nPy * 2 )
       PARAMETER ( W2_ioBufferSize = W2_maxNbTiles*sNx*sNy )
       PARAMETER ( W2_maxXStackNx = W2_maxNbTiles*sNx )
       PARAMETER ( W2_maxXStackNy = W2_maxNbTiles*sNy )
       PARAMETER ( W2_maxYStackNx = W2_maxNbTiles*sNx )
       PARAMETER ( W2_maxYStackNy = W2_maxNbTiles*sNy )

C- Note: Overestimating W2_maxNbFacets and, to less extent, W2_maxNeighbours
C        have no or very little effects on memory footprint.
C        overestimated W2_maxNbTiles does not have large effect, except
C        through ioBufferSize (if related to, as here).
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: W2_EXCH2_TOPOLOGY.h
C     !INTERFACE:
C     #include W2_EXCH2_TOPOLOGY.h

C     !DESCRIPTION:
C     *==========================================================*
C     | W2_EXCH2_TOPOLOGY.h
C     | o Header defining tile exchange and mapping for W2_EXCH2
C     *==========================================================*
C     | 1rst part holds the full topology structure (same for all
C     |  process) and is independent of tile-processor repartition
C     |  (needs W2_EXCH2_SIZE.h to be included before)
C     | 2nd part (put in this header for convenience) holds
C     |   Tile Ids and is function of tile-process repartition
C     |  (needs SIZE.h to be included before)
C     *==========================================================*
CEOP

C---   Parameters for enumerating directions
       INTEGER W2_NORTH, W2_SOUTH, W2_EAST, W2_WEST
       PARAMETER ( W2_NORTH = 1 )
       PARAMETER ( W2_SOUTH = 2 )
       PARAMETER ( W2_EAST  = 3 )
       PARAMETER ( W2_WEST  = 4 )

C---   Topology data structures
C      exch2_global_Nx   :: Global-file domain length.
C      exch2_global_Ny   :: Global-file domain height.
C      exch2_xStack_Nx   :: Length of domain used for north/south OBCS.
C      exch2_xStack_Ny   :: Height of domain used for north/south OBCS.
C      exch2_yStack_Nx   :: Length of domain used for east/west OBCS.
C      exch2_yStack_Ny   :: Height of domain used for east/west OBCS.
C---   Tiling and Exch data structures
C      exch2_nTiles      :: Number of tiles in this topology
C      exch2_myFace      :: Face number for each tile (used for I/O).
C      exch2_mydNx       :: Face size in X for each tile (for I/O).
C      exch2_mydNy       :: Face size in Y for each tile (for I/O).
C      exch2_tNx         :: Size in X for each tile.
C      exch2_tNy         :: Size in Y for each tile.
C      exch2_tBasex      :: Tile offset in X within its sub-domain (cube face)
C      exch2_tBasey      :: Tile offset in Y within its sub-domain (cube face)
C      exch2_txGlobalo   :: Tile base X index within global index space.
C      exch2_tyGlobalo   :: Tile base Y index within global index space.
C      exch2_txXStackLo  :: Tile base X index within N/S OBCS index space.
C      exch2_tyXStackLo  :: Tile base Y index within N/S OBCS index space.
C      exch2_txYStackLo  :: Tile base X index within E/W OBCS index space.
C      exch2_tyYStackLo  :: Tile base Y index within E/W OBCS index space.
C      exch2_isWedge     :: 1 if West  is at domain edge, 0 if not.
C      exch2_isNedge     :: 1 if North is at domain edge, 0 if not.
C      exch2_isEedge     :: 1 if East  is at domain edge, 0 if not.
C      exch2_isSedge     :: 1 if South is at domain edge, 0 if not.
C      exch2_nNeighbours :: Tile neighbour entries count.
C      exch2_neighbourId :: Tile number for each neighbour entry.
C      exch2_opposingSend:: Neighbour entry in target tile send
C                        :: which has this tile and neighbour as its target.
C      exch2_pij(:,n,t)  :: Matrix which applies to target-tile indices to get
C                        :: source-tile "t" indices, for neighbour entry "n".
C      exch2_oi(n,t)     :: Source-tile "t" X index offset in target
C                        :: to source connection (neighbour entry "n").
C      exch2_oj(n,t)     :: Source-tile "t" Y index offset in target
C                        :: to source connection (neighbour entry "n").
       INTEGER exch2_global_Nx
       INTEGER exch2_global_Ny
       INTEGER exch2_xStack_Nx
       INTEGER exch2_xStack_Ny
       INTEGER exch2_yStack_Nx
       INTEGER exch2_yStack_Ny
       INTEGER exch2_nTiles
       INTEGER exch2_myFace( W2_maxNbTiles )
       INTEGER exch2_mydNx( W2_maxNbTiles )
       INTEGER exch2_mydNy( W2_maxNbTiles )
       INTEGER exch2_tNx( W2_maxNbTiles )
       INTEGER exch2_tNy( W2_maxNbTiles )
       INTEGER exch2_tBasex( W2_maxNbTiles )
       INTEGER exch2_tBasey( W2_maxNbTiles )
       INTEGER exch2_txGlobalo(W2_maxNbTiles )
       INTEGER exch2_tyGlobalo(W2_maxNbTiles )
       INTEGER exch2_txXStackLo(W2_maxNbTiles )
       INTEGER exch2_tyXStackLo(W2_maxNbTiles )
       INTEGER exch2_txYStackLo(W2_maxNbTiles )
       INTEGER exch2_tyYStackLo(W2_maxNbTiles )
       INTEGER exch2_isWedge( W2_maxNbTiles )
       INTEGER exch2_isNedge( W2_maxNbTiles )
       INTEGER exch2_isEedge( W2_maxNbTiles )
       INTEGER exch2_isSedge( W2_maxNbTiles )
       INTEGER exch2_nNeighbours( W2_maxNbTiles )
       INTEGER exch2_neighbourId(  W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_opposingSend( W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_neighbourDir( W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_pij(4,W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_oi (  W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_oj (  W2_maxNeighbours, W2_maxNbTiles )

       COMMON /W2_EXCH2_TOPO_I/
     &        exch2_global_Nx, exch2_global_Ny,
     &        exch2_xStack_Nx, exch2_xStack_Ny,
     &        exch2_yStack_Nx, exch2_yStack_Ny,
     &        exch2_nTiles,
     &        exch2_myFace, exch2_mydNx, exch2_mydNy,
     &        exch2_tNx, exch2_tNy,
     &        exch2_tBasex, exch2_tBasey,
     &        exch2_txGlobalo,exch2_tyGlobalo,
     &        exch2_txXStackLo,exch2_tyXStackLo,
     &        exch2_txYStackLo,exch2_tyYStackLo,
     &        exch2_isWedge, exch2_isNedge,
     &        exch2_isEedge, exch2_isSedge,
     &        exch2_nNeighbours, exch2_neighbourId,
     &        exch2_opposingSend, exch2_neighbourDir,
     &        exch2_pij,
     &        exch2_oi, exch2_oj

C---   Exchange execution loop data structures
C      exch2_iLo,iHi(n,t) :: X-index range of this tile "t" halo-region
C                         :: to be updated with neighbour entry "n".
C      exch2_jLo,jHi(n,t) :: Y-index range of this tile "t" halo-region
C                         :: to be updated with neighbour entry "n".
       INTEGER exch2_iLo( W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_iHi( W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_jLo( W2_maxNeighbours, W2_maxNbTiles )
       INTEGER exch2_jHi( W2_maxNeighbours, W2_maxNbTiles )
       COMMON /W2_EXCH2_HALO_SPEC/
     &        exch2_iLo, exch2_iHi,
     &        exch2_jLo, exch2_jHi

C---   Cumulated Sum operator
C      W2_tMC1, W2_tMC2 :: tile that holds Missing Corners (=f1.NW,f2.SE)
C      W2_cumSum_facet(1,f1,f2) :: cum-sum at facet f2 origin function of
C                                  facet f1 X-increment
C      W2_cumSum_facet(2,f1,f2) :: cum-sum at tile f2 origin function of
C                                  facet f1 Y-increment
C      W2_cumSum_tiles(1,t1,t2) :: cum-sum at tile t2 origin function of
C                                  tile t1 X-increment
C      W2_cumSum_tiles(2,t1,t2) :: cum-sum at tile t2 origin function of
C                                  tile t1 Y-increment
       INTEGER W2_tMC1, W2_tMC2
       INTEGER W2_cumSum_facet( 2, W2_maxNbFacets,W2_maxNbFacets)
       COMMON /W2_CUMSUM_TOPO_I/
     &        W2_tMC1, W2_tMC2,
     &        W2_cumSum_facet

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   COMMON /W2_MAP_TILE2PROC/ mapping of tiles to process:
C     get W2 tile Id from process Id + subgrid indices (bi,bj) or the reverse
C     (tile ids are no longer a simple function of process and subgrid indices).
C
C     W2_tileProc(tN) :: Rank of process owning tile tN (filled at run time).
C     W2_tileIndex(tN):: local subgrid index of tile tN
C     W2_tileRank(tN) :: rank of tile tN in full-tile list (without blank)
C     W2_myTileList   :: list of tiles owned by this process
C     W2_procTileList :: same as W2_myTileList, but contains
C                        information for all processes
      INTEGER W2_tileProc ( W2_maxNbTiles )
      INTEGER W2_tileIndex( W2_maxNbTiles )
c     INTEGER W2_tileRank ( W2_maxNbTiles )
      INTEGER W2_myTileList ( nSx,nSy )
      INTEGER W2_procTileList(nSx,nSy,nPx*nPy )
      COMMON /W2_MAP_TILE2PROC/
     &        W2_tileProc,
     &        W2_tileIndex,
c    &        W2_tileRank,
     &        W2_myTileList, W2_procTileList

C--   COMMON /W2_EXCH2_COMMFLAG/ EXCH2 character Flag for type of communication
      CHARACTER W2_myCommFlag( W2_maxNeighbours, nSx, nSy )
      COMMON /W2_EXCH2_COMMFLAG/ W2_myCommFlag

C--   COMMON /EXCH2_FILLVAL_RX/ real type filling value used by EXCH2
C     e2FillValue_RX :: filling value for null regions (facet-corner
C                    :: halo regions)
      Real*8 e2FillValue_RL
      Real*8 e2FillValue_RS
      Real*4 e2FillValue_R4
      Real*8 e2FillValue_R8
      COMMON /EXCH2_FILLVAL_RL/ e2FillValue_RL
      COMMON /EXCH2_FILLVAL_RS/ e2FillValue_RS
      COMMON /EXCH2_FILLVAL_R4/ e2FillValue_R4
      COMMON /EXCH2_FILLVAL_R8/ e2FillValue_R8

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     thisTile      :: receiving tile Id. number
C     nN            :: Neighbour entry that we are processing
C     e2BufrRecSize :: Number of elements in each entry of e2Bufr1_R4
C     iBufr         :: number of buffer elements to transfer
C     e2Bufr1_R4    :: Data transport buffer array. This array is used in one of
C                   :: two ways. For PUT communication the entry in the buffer
C                   :: associated with the source for this receive (determined
C                   :: from the opposing_send index) is read.
C                   :: For MSG communication the entry in the buffer associated
C                   :: with this neighbor of this tile is used as a receive
C                   :: location for loading a linear stream of bytes.
C     commSetting   :: Mode of communication used to exchange with this neighbor
C     myThid        :: my Thread Id. number

      INTEGER thisTile, nN
      INTEGER e2BufrRecSize
      INTEGER iBufr
      Real*4     e2Bufr1_R4( e2BufrRecSize )
      CHARACTER commSetting
      INTEGER myThid
CEOP

C !LOCAL VARIABLES:
C     == Local variables ==
C     soT     :: Source tile Id. number
C     oNb     :: Opposing send record number
      INTEGER soT
      INTEGER oNb

C     MPI setup
      INTEGER theTag, theType
      INTEGER sProc, tProc
      INTEGER mpiStatus(MPI_STATUS_SIZE), mpiRc

      soT = exch2_neighbourId(nN, thisTile )
      oNb = exch2_opposingSend(nN, thisTile )

C     Handle receive end data transport according to communication mechanism
C     between source and target tile
      IF ( commSetting .EQ. 'M' ) THEN
C      Setup MPI stuff here
       theTag =  (soT-1)*W2_maxNeighbours + oNb
       tProc = W2_tileProc(thisTile)-1
       sProc = W2_tileProc(soT)-1
       theType = MPI_REAL
       CALL MPI_Recv( e2Bufr1_R4, iBufr, theType, sProc,
     &                theTag, MPI_COMM_MODEL, mpiStatus, mpiRc )
      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
