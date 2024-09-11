
% Please run this script in order to reproduce all simulations.

currentDir = pwd;
initCobraToolbox
solverOK = changeCobraSolver('ibm_cplex','all');
cd(currentDir)

addpath([pwd filesep 'Scripts'])

% Optional: create sample-specific models
%% Note: due to differeces in gap-filling, newly created models
createModelsRetina

% Start from here to reproduce the data from the article
computeFluxesRetina
exportFluxesRetina
doStatisticalAnalysisRetina
calculateCorrelationsRetina