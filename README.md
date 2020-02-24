# OSIRIS_Model
 
The OSIRIS marine ecosystem model is designed to be relatively fast and simple, and to allow examination of the influence of multiple interacting and simulataneous stressors on marine ecosystems. Running many model variants (using the facility for Latin-hypercube sampling) allows uncertainty in the results to be quantified.
The model was written in Matlab, version R2019b, by Richard Bailey, and runs on Mac OS, Windows, Linux.

## Quickstart:
Once all folders and files are downloaded, in Matlab open and run “OSIRIS_v1_0_0.m”.
Within “OSIRIS_v1_0_0.m” there are 4 parameters to choose (`Prefix, n_LHS, plot_network, plot_timeseries`) and notes on these are provided in comments within the code. The model downloaded from this repository is set up to run for 30 years, running 2 model variants (`n_LHS=2`). The first LHS run always uses the central parameter value, with no uncertainty, and all other model runs use parameters sampled within the uncertainties. The time series plotted (when `n_LHS>1` are for the last model run). To immediately see results with no uncertainty (in forcing or parameters), set `n_LHS=1`. 

## Model inputs and outputs:
Parameters are specified in the excel sheet “OSIRIS_Model_configuration_v1.0.0.xlsm”. Model input files can be generated using this sheet (or can be generated manually by the user). The visual basic script (“saveCSVfile()” in module 3) which produces the input files contains details of which data (within the excel sheet tabs) are reproduced in which files. Notes on using the sheet are provided in the sheet itself. Once created, these files are stored in the “config_files” directory. (A note for Mac OS users: the first time configuration files are creeated using “OSIRIS_Model_configuration_v1.0.0.xlsm”, file permission is requested independently for each of the 31 files, which slows down the process considerably. This is only needed the first time, however, or if the model is moved to a new directory.)

Model outputs are found in the “Output_files” directory, with one file per model run, prefixed by the name set as `Prefix`. Each files contains a time series (columns, in order) for each node state. The time resolution is specified by the user in “OSIRIS_Model_configuration_v1.0.0.xlsm”. 

## Reference
The full description of the model can be found in:
Bailey, R.M., van der Grient, J.M.A. (2020)  OSIRIS: a model for integrating the effects of multiple stressors on marine ecosystems. Journal of Theoretical Biology (accepted).
