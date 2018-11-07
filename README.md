This collection of files is an implementation in R of a simulation and bias correction of wind power generation in Brazil
as well as validation and investigation of correlations between El Nino and La Nina events as well as water inflows into Brazilian 
hydropower plants and wind power

Details about the simulation and data can be found in this Master Thesis: https://doi.org/10.5281/zenodo.1471221
Results are found here: https://doi.org/10.5281/zenodo.1470357

The starting file, where the simulation is performed is "RScript_16.R". For it to work, also the files "MERRA_data.R" and "functions_16.R" 
(the first is for download, conversion and loading of MERRA-2 Data and the second contains all the functions used in "RScript_16.R")

After simulation, the results are analysed. First they are prepared for analysis, which happens in the files "prepare_for_analysis_states.R" 
and "prepare_for_analysis_subsystems.R" (spatially aggregated for the states and subsystems + Brazil), the statistical analysis happens in 
the files "analysis_statistics_states.R" and "analysis_statistics_subsystems.R". For the analysis also some functions need to be loaded
which are found in the file "functions_analysis.R".

Some exemplary plots of the results are shown in "plotting.R".

The correlations with El Nino and La Nina events are examined with the code in the files "analysis_nino_7.R" and "analysis_nino7states.R" 
(again for subsystems + Brazil and the states). Also Linear Models are tested with and without including El Nino and La Nina events, which
happens in the file "LinearModels_BNESSt.R". For this also the file "functions_LM.R" is needed.

Finally, also correlations between water inflows and simulated wind power generation are examined in the file "analysis_water3.R".
