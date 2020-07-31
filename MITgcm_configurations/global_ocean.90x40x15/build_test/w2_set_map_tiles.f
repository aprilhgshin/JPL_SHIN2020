












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


C--  File w2_set_map_tiles.F:
C--   Contents
C--   o W2_SET_MAP_TILES :: Set tiles and IO mapping
C--   o FIND_GCD_N       :: Returns the Greatest Common Divisor

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: W2_SET_MAP_TILES

C !INTERFACE:
      SUBROUTINE W2_SET_MAP_TILES( myThid )

C     !DESCRIPTION:
C     Set-up tiles mapping and IO global mapping

C     !USES:
      IMPLICIT NONE

C      Tile topology settings data structures
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
C     !ROUTINE: W2_EXCH2_PARAMS.h
C     !INTERFACE:
C     #include W2_EXCH2_PARAMS.h

C     !DESCRIPTION:
C     *==========================================================*
C     | W2_EXCH2_PARAMS.h
C     | o Header file defining WRAPPER2-EXCH2 topology parameters
C     *==========================================================*
CEOP

C--   COMMON /W2_EXCH2_PARM_I/ Integer valued parameters used by W2-EXCH2
C     preDefTopol :: pre-defined Topology selector:
C                 :: = 0 : topology defined from processing "data.exch2";
C                 :: = 1 : simple, single facet topology;
C                 :: = 2 : customized topology (w2_set_myown_facets)
C                 :: = 3 : 6-face Cube (3 face-dims: nRed, nGreen, nBlue).
C     nFacets     :: Number of facets (or domains) in this topology
C     facet_dims  :: facet pair of dimensions (n1x,n1y, n2x,n2y, ...)
C     nTiles      :: Number of tiles in this topology <- moved to W2_EXCH2_TOPOLOGY.h
C     nBlankTiles :: Number of "Blank-Tiles"
C     blankList   :: List of "Blank-Tiles" (non active)
C--
C     W2_mapIO    :: select option for global-IO mapping:
C             =-1 :: old format, put facets 1 after the other in the X dir.;
C                 :: this is not necessary "compact".
C             = 1 :: compact format, 1 facet after the other (mostly in Y dir.)
C                 :: but needs to fold some facet (domain) if too large
C             = 0 :: compact format (= 1 long line), one facet after the other.
C--
C     W2_oUnit    :: output fortran Unit for printing Std messages
C     W2_printMsg :: select option for information messages printing
C             < 0 :: open and print to "w2_tile_topology" log file
C             = 0 :: print the minimum, to StdOut
C             = 1 :: no duplicated print from different processes (only Master)
C             = 2 :: all processes do print (used to check).
      INTEGER preDefTopol
      INTEGER nFacets
      INTEGER facet_dims(2*W2_maxNbFacets)
      INTEGER nBlankTiles
      INTEGER blankList(W2_maxNbTiles)
      INTEGER W2_mapIO
      INTEGER W2_oUnit, W2_printMsg
      COMMON /W2_EXCH2_PARM_I/
     &        preDefTopol,
     &        nFacets, facet_dims,
     &        nBlankTiles, blankList,
     &        W2_mapIO,
     &        W2_oUnit, W2_printMsg

C--   COMMON /W2_EXCH2_PARM_L/ Logical valued parameters used by W2-EXCH2
C     W2_useE2ioLayOut :: =T: use Exch2 global-IO Layout; =F: use model default
      LOGICAL W2_useE2ioLayOut
      COMMON /W2_EXCH2_PARM_L/
     &        W2_useE2ioLayOut

C--   COMMON /W2_EXCH2_BUILD_I/ Integer variables used to build topology
C     facet_owns  :: Range of tiles this facet "owns"
C     facet_pij   \  ::
C     facet_oi     } :: indices correspondence matrix (facet_pij) & offsets:
C     facet_oj    /  ::
C-with:  suffix "so" for indices of source facet j ;
C        suffix "tg" for indices of target facet jj= INT(facet_link(i,j))
C      pij(:,i,j) : matrix which gives so indices when applied to tg indices
C        iso = pij(1)*itg + pij(2)*jtg + oi
C        jso = pij(3)*itg + pij(4)*jtg + oj
C-----
      INTEGER facet_owns(2,W2_maxNbFacets)
      INTEGER facet_pij(4,4,W2_maxNbFacets)
      INTEGER facet_oi(4,W2_maxNbFacets)
      INTEGER facet_oj(4,W2_maxNbFacets)
      COMMON /W2_EXCH2_BUILD_I/
     &        facet_owns,
     &        facet_pij, facet_oi, facet_oj

C--   COMMON /W2_EXCH2_PARM_R/ Real*4 valued parameters used by W2-EXCH2
C--   topology defined from processing "data.exch2" (preDefTopol=0):
C     facet_link  :: Face-Edge connectivity map:
C       facet_link(i,j)=XX.1 : face(j)-edge(i) (i=1,2,3,4 <==> N,S,E,W)
C       is connected to Northern edge of face "XX" ; similarly,
C       = XX.2 : to Southern.E, XX.3 = Eastern.E, XX.4 = Western.E of face "XX".
      Real*4  facet_link( 4, W2_maxNbFacets )
      COMMON /W2_EXCH2_PARM_R/ facet_link

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

C     !INPUT PARAMETERS:
C     myThid  :: my Thread Id number
C               (Note: not relevant since threading has not yet started)
      INTEGER myThid

C     !FUNCTIONS:
      INTEGER  FIND_GCD_N
      EXTERNAL FIND_GCD_N

C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf     :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER tNx, tNy, fNx, fNy, nbPts, fBaseX, fBaseY
      INTEGER nbTx, nbTy
      INTEGER j, ii, k, tId, tx, ty
      INTEGER divide, nnx(W2_maxNbFacets)
      INTEGER errCnt, tCnt
      LOGICAL tileIsActive, prtFlag
CEOP

C     Set-up tiles mapping and IO global mapping
      WRITE(msgBuf,'(2A)') 'W2_SET_MAP_TILES:',
     &       ' tile mapping within facet and global Map:'
      CALL PRINT_MESSAGE( msgBuf, W2_oUnit, SQUEEZE_RIGHT, myThid )
      prtFlag = ABS(W2_printMsg).GE.2
     &       .OR. ( W2_printMsg .NE.0 .AND. myProcId.EQ.0 )

      tNx = sNx
      tNy = sNy
C--   Check that tile dims divide facet dims
      errCnt = 0
      tCnt = 0
      nbPts = 0
      DO j=1,nFacets
        fNx = facet_dims(2*j-1)
        fNy = facet_dims( 2*j )
        nbTx = fNx/tNx
        nbTy = fNy/tNy
        IF ( nbTx*tNx .NE. fNx ) THEN
          WRITE(msgBuf,'(A,I3,2(A,I7))') 'Facet',j,
     &      ' : X-size=', fNx, ' not multiple of sNx=', tNx
          CALL PRINT_ERROR( msgBuf, myThid )
          errCnt = errCnt + 1
        ENDIF
        IF ( nbTy*tNy .NE. fNy ) THEN
          WRITE(msgBuf,'(A,I3,2(A,I7))') 'Facet',j,
     &      ' : Y-size=', fNy, ' not multiple of sNy=', tNy
          CALL PRINT_ERROR( msgBuf, myThid )
          errCnt = errCnt + 1
        ENDIF
        facet_owns(1,j) = tCnt+1
        tCnt = tCnt + nbTx*nbTy
        facet_owns(2,j) = tCnt
        nbPts = nbPts + fNx*fNy
      ENDDO
      IF ( errCnt.GT.0 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &   ' W2_SET_MAP_TILES: found', errCnt, ' Fatal errors'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R W2_SET_MAP_TILES'
      ENDIF
C--   Check that domain size and (SIZE.h + blankList) match:
      IF ( tCnt.NE.exch2_nTiles ) THEN
        WRITE(msgBuf,'(A,I6,A)')
     &   'W2_SET_MAP_TILES: Domain Total # of tiles =', tCnt, ' does'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A,I6)')
     &   'W2_SET_MAP_TILES:  not match (SIZE.h+blankList)=',exch2_nTiles
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R W2_SET_MAP_TILES'
      ENDIF

      IF ( W2_mapIO.EQ.1 ) THEN
C--   Compact IO map (mostly in Y dir): search for Greatest Common Divisor
C     of all x-size (faster to apply GCD to Nb of Tiles in X):
        k = 0
        nnx(1) = 0
        DO j=1,nFacets
C     skip empty facet
          IF ( facet_dims(2*j-1).GT.0 ) THEN
            k = k + 1
            nnx(k) = facet_dims(2*j-1)/tNx
          ENDIF
        ENDDO
        divide = FIND_GCD_N( nnx, k )
        W2_mapIO = divide*tNx
        WRITE(msgBuf,'(A,2(I5,A))') ' W2_mapIO =', W2_mapIO,
     &                              ' (=', divide, '*sNx)'
        CALL PRINT_MESSAGE( msgBuf, W2_oUnit, SQUEEZE_RIGHT, myThid )
      ENDIF

C--   Global Map size:
C     facets stacked in x direction
      exch2_xStack_Nx = 0
      exch2_xStack_Ny = 0
      DO j=1,nFacets
        exch2_xStack_Nx =   exch2_xStack_Nx + facet_dims(2*j-1)
        exch2_xStack_Ny = MAX( exch2_xStack_Ny, facet_dims(2*j) )
      ENDDO
C     facets stacked in y direction
      exch2_yStack_Nx = 0
      exch2_yStack_Ny = 0
      DO j=1,nFacets
        exch2_yStack_Nx = MAX( exch2_yStack_Nx, facet_dims(2*j-1) )
        exch2_yStack_Ny =   exch2_yStack_Ny + facet_dims(2*j)
      ENDDO
      IF ( W2_mapIO.EQ.-1 ) THEN
        exch2_global_Nx = exch2_xStack_Nx
        exch2_global_Ny = exch2_xStack_Ny
      ELSEIF ( W2_mapIO.EQ.0 ) THEN
        exch2_global_Nx = nbPts
        exch2_global_Ny = 1
      ELSE
        exch2_global_Nx = W2_mapIO
        exch2_global_Ny = nbPts/W2_mapIO
      ENDIF
      WRITE(msgBuf,'(A,2(A,I8))') ' Global Map (IO):',
     &  ' X-size=', exch2_global_Nx, ' , Y-size=', exch2_global_Ny
      CALL PRINT_MESSAGE( msgBuf, W2_oUnit, SQUEEZE_RIGHT, myThid )

C--   Set tiles mapping within facet (sub-domain) and within Global Map
      WRITE(msgBuf,'(2A)') 'W2_SET_MAP_TILES:',
     &       ' tile offset within facet and global Map:'
      CALL PRINT_MESSAGE( msgBuf, W2_oUnit, SQUEEZE_RIGHT, myThid )
      tId = 0
      nbPts = 0
      fBaseX = 0
      fBaseY = 0
      DO j=1,nFacets
        fNx = facet_dims(2*j-1)
        fNy = facet_dims( 2*j )
        nbTx = fNx/tNx
        nbTy = fNy/tNy
        WRITE(W2_oUnit,'(A,I3,2(A,I6),A,I5,2(A,I4),A)')
     &    '- facet', j, ' : X-size=', fNx, ' , Y-size=', fNy,
     &    ' ;', nbTx*nbTy, ' tiles (Tx,Ty=', nbTx,',',nbTy,')'
c       CALL PRINT_MESSAGE( msgBuf, W2_oUnit, SQUEEZE_RIGHT, myThid )
        DO ty=1,nbTy
         DO tx=1,nbTx
          tId = tId + 1
C--   Tags blank tile by removing facet # (exch2_myFace) but keeps its location
          tileIsActive = .TRUE.
          DO k=1,nBlankTiles
           IF ( blankList(k).EQ.tId ) tileIsActive = .FALSE.
          ENDDO
          IF ( tileIsActive ) exch2_myFace(tId) = j
          exch2_mydNx ( tId ) = fNx
          exch2_mydNy ( tId ) = fNy
          exch2_tNx   ( tId ) = tNx
          exch2_tNy   ( tId ) = tNy
          exch2_tBasex( tId ) = (tx-1)*tNx
          exch2_tBasey( tId ) = (ty-1)*tNy
C--   Global IO Mappings
C       these are for OBCS (vertical slices)
          exch2_txXStackLo( tId ) = 1 + exch2_tBasex(tId) + fBaseX
          exch2_tyXStackLo( tId ) = 1 + exch2_tBasey(tId)
          exch2_txYStackLo( tId ) = 1 + exch2_tBasex(tId)
          exch2_tyYStackLo( tId ) = 1 + exch2_tBasey(tId) + fBaseY
C       and these for global files (3d files/horizontal 2d files)
          IF ( W2_mapIO.EQ.-1 ) THEN
C-        Old format
            exch2_txGlobalo( tId ) = 1 + exch2_tBasex(tId) + fBaseX
            exch2_tyGlobalo( tId ) = 1 + exch2_tBasey(tId)
          ELSEIF ( W2_mapIO.EQ.0 ) THEN
C-        Compact format = 1 long line
            ii = nbPts + exch2_tBasex(tId) + exch2_tBasey(tId)*fNx
            exch2_txGlobalo( tId ) = 1 + ii
            exch2_tyGlobalo( tId ) = 1
          ELSE
C         Compact format: piled in the Y direction
            ii = nbPts + exch2_tBasex(tId) + exch2_tBasey(tId)*fNx
            exch2_txGlobalo( tId ) = 1 + MOD(ii,W2_mapIO)
            exch2_tyGlobalo( tId ) = 1 + ii/W2_mapIO
          ENDIF
          IF ( prtFlag )
     &    WRITE(W2_oUnit,'(A,I5,3(A,I3),2A,2I5,2A,2I8)') '  tile',tId,
     &    ' on facet', exch2_myFace(tId),' (',tx,',',ty,'):',
     &         ' offset=', exch2_tBasex(tId), exch2_tBasey(tId),' ;',
     &    ' on Glob.Map=', exch2_txGlobalo(tId),exch2_tyGlobalo(tId)
         ENDDO
        ENDDO
        fBaseX = fBaseX + fNx
        fBaseY = fBaseY + fNy
        nbPts = nbPts + fNx*fNy
      ENDDO

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: FIND_GCD_N

C     !INTERFACE:
      INTEGER FUNCTION FIND_GCD_N( fldList, nFld )

C     !DESCRIPTION:
C     *==========================================================*
C     | FUNCTION FIND_GCD_N
C     | o Find the Greatest Common Divisor of N integers
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     !INPUT PARAMETERS:
C     fldList :: list of integers to search for GCD
C     nFLd    :: length of the input integer list.
      INTEGER nFLd
      INTEGER fldList(nFld)

C     !LOCAL VARIABLES:
      INTEGER mnFld, divide
      INTEGER j, ii
      LOGICAL flag
      LOGICAL localDBg
CEOP
      PARAMETER ( localDBg = .FALSE. )
c     PARAMETER ( localDBg = .TRUE. )

      mnFld = fldList(1)
      DO j=1,nFld
        mnFld = MIN( mnFld, fldList(j) )
      ENDDO
      IF (localDBg) WRITE(0,'(A,I8)') 'FIND_GCD_N: mnFld=',mnFld

      IF (mnFld.GT.1 ) THEN
        divide = 1
        ii = 2
        DO WHILE ( ii.LE.mnFld )
          IF (localDBg) WRITE(0,'(A,I8)') '     GCD  : try',ii
          flag = .TRUE.
          DO j=1,nFld
            flag = flag.AND.(MOD(fldList(j),ii).EQ.0 )
          ENDDO
          IF ( flag ) THEN
            divide = divide*ii
            DO j=1,nFld
              fldList(j) = fldList(j)/ii
            ENDDO
            IF (localDBg) WRITE(0,'(A,I8)')
     &                              'FIND_GCD_N: com.fact=',ii
            mnFld = mnFld/ii
          ELSE
            ii = ii+2
            IF (ii.EQ.4) ii=3
          ENDIF
        ENDDO
C-    Put back the original Nb:
        IF (localDBg) WRITE(0,'(10I8)') (fldList(j),j=1,nFld)
        DO j=1,nFld
          fldList(j) = fldList(j)*divide
        ENDDO
      ELSE
        divide = MAX( 0, mnFld )
      ENDIF

      FIND_GCD_N = divide

      RETURN
      END
