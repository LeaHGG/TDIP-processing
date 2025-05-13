# TDIP-processing
TDIP processing by stretched exponential

Starting files = mat files for logging TDIP data and text files (.tx2 extension) for surface and cross-borehole TDIP data. All data are found in the “input data” folder.

1) SE FITS 

SEfit_from_mat_logging.m calls .mat files in Data repository\Input data\Nesjavellir\
SEfit_from_Tx2_surface.m calls .tx2 files in Data repository\Input data\Krafla\
SEfit_from_Tx2_XB.m calls .tx2 files in Data repository\Input data\Hvedemarken\
The three scripts are relatively similar. They load data and carry out the six steps of the processing procedure described in the paper. They plot the results and save the results as a .mat file which can be used for the next step (Debye decomposition)

2) DEBYE DECOMPOSITION

GEMGAS_RTD_v5.m calls .mat files (Output data from scripts above)
Krafla_RTD_v4.m calls .mat files (Output data from scripts above)
Hvede_RTD_v4.m calls .mat files (Output data from scripts above)



