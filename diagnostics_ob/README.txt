
---------------------------------------------------------------------------------------------
The following files in the diagnostics_ob package are used for the Open Boundary Extraction:
---------------------------------------------------------------------------------------------


diagnostics_ob_readparms.F
- Read data.diagnostics_ob input parameters

diagnostics_ob_init_fixed.F
- Initialize variables fixed during the model's time loop

DIAGNOSTICS_OB.h
- Declare all variables in COMMON blocks

USER_INPUT.h
- Additional input parameters

ob_extract_output.F
- Contains the OB_EXTRACT_OUTPUT and some supporting subroutines that executes entire Open Boundary Extraction Program:
      - MASTER_PROC_TASKS
      - APPEND_OB_VALUE2D
      - APPEND_OB_VALUE3D
      - OB_PASS_RL_to_R8
      - PLOT_TEST_GLOBAL_OB
      - PRINT_INT_ARR, PRINT_FLOAT_ARR
      - PLOT_GLO_FIELD_XYRL


diag_ob_prepare_subFieldOnMasks.F
- Contains all supporting subroutines that extract field values on open boundary points, accumulates the extracted values, and compute time averaged fields.
      - SET_TO_INI_STATE
      - SET_SUBFIELDS
      - CUMULATE_FLD_TIME
      - TIMEAVE_OB_FLD

diagnostics_ob_write_bin.F
- Writes final output to binary file with subroutines:
      - WRITE_GLOBAL_BIN
      - CREATE_FILENAME
      - OB_R8_to_R4

print_ob_output.py
- python test script to compare outputs with expected outputs



Other files are included for model convention.
