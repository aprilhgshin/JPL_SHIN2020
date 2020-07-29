C
C     *===============================================================*
C     | USER_INPUT.h Contains user input for open boundary extraction.
C     *===============================================================*
C
C     nOB_mask  :: Number of open boundary masks.
C     nOB_fld  :: Number of fields.
C     ob_fnames :: Filenames of each open boundary mask
C     ob_hydrogSalt :: logical value of whether salinity field should be outputted
C     ob_hydrogTheta :: logical value of whether theta field should be outputted
C     ob_bathymetry :: logical value of whether bathymetry field should be outputted

      INTEGER, PARAMETER :: nOB_mask =  1
      INTEGER, PARAMETER :: nOB_fld = 2

      LOGICAL ob_hydrogSalt
      LOGICAL ob_hydrogTheta
      LOGICAL ob_bathymetry
      PARAMETER (
     &           ob_bathymetry = .TRUE.,
     &           ob_hydrogSalt = .FALSE.,
     &           ob_hydrogTheta = .TRUE.)
