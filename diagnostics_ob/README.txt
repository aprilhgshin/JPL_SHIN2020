
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
- Contains the OB_EXTRACT_OUTPUT subroutine that executes entire Open Boundary Extraction Program

diag_ob_prepare_subFieldOnMasks.F
- Contains all subroutines that extract field values on open boundary points, accumulates the extracted values, and compute time averaged fields.

diagnostics_ob_write_bin.F
- Writes final output to binary file

print_ob_output.py
- python test script to compare outputs with expected outputs



Other files are included for model convention.
