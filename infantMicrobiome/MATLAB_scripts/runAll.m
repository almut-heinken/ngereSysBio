
rootDir = pwd;

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

cd(rootDir)

addpath(genpath([rootDir filesep 'Scripts']))
addpath([rootDir filesep 'inputFiles'])

% perform gapfilling with HMO module and testing
runRefinement
plotUsedHMOs

% normalize abundances
createAbundanceFiles

% build microbiome models and perform simulations
runMgPipe_Infants
runMgPipe_Mothers

% analyze the generated data
statisticalAnalysisInfants
statisticalAnalysisInfantsMothers
getSignificantFeatures
getVDvsCSDReactions
exportResults
plotInfantsMothers
exportContributionSummaryTables
