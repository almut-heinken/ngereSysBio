
% Please run this script in order to reproduce all simulations.

currentDir = pwd;
initCobraToolbox
solverOK = changeCobraSolver('ibm_cplex','all');
cd(currentDir)

addpath([pwd filesep 'Scripts'])

createModels
computeFluxes
exportFluxes
doStatisticalAnalysis
calculateCorrelations