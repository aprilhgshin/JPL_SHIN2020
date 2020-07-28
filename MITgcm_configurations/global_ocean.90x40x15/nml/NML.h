C ======================================================================
C NML Parameters for namelist demo
C ======================================================================


      INTEGER dimension_x, dimension_y, dimension_z
      CHARACTER*10 nml_fields(4)
      REAL*8 nml_map(10,10)


      COMMON / NML_PARMS_I /
     &  dimension_x, dimension_y, dimension_z

      COMMON / NML_PARMS_C /
     &  nml_fields

      COMMON / NML_PARMS_R /
     &  nml_map

CEOF
