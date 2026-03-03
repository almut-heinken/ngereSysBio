
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
solverOK=changeCobraSolver('ibm_cplex','MILP');

addpath('Scripts')

% run simulations
runModelBuilding_brain
runSimulations_brain
exportFluxFeatures_brain
runStatistics_brain
plotMinMaxFluxes_brain