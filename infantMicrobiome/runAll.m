
% Execute this script to run all simulations and data analyses.

% Create a variable containing the current path (required)
rootDir = pwd;

% Start the simulations
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
cd(rootDir)

addpath(genpath([rootDir filesep 'MATLAB_scripts']))
addpath([rootDir filesep 'inputFiles'])

%% OPTIONAL %%
% Refine 289 additional draft reconstructions
runReconstruction
% NOTE: Due to some randomness in automated gap-filling in DEMETER, the 
% resulting refined reconstructions may differ slightly from the previously  
% generated reconstructions published with the paper and used in the 
% following analyses. 

%% REPRODUCTION OF SIMULATIONS AND ANALYSES
% To extract the version of the 289 refined reconstructions used in the
% paper:
unzip('refinedReconstructionsPublished')

% To retrieve AGORA2 from the Virtual Metabolic Human website:
mkdir('AGORA2')
agora2_info = readInputTableForPipeline('AGORA2_infoFile.xlsx');
for i=2:size(agora2_info,1)
    websave(['AGORA2' filesep agora2_info{i,1} '.mat'],['https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/individual_reconstructions/' agora2_info{i,1} '.mat'])
end

% perform gapfilling with HMO module and testing
runRefinement
plotUsedHMOs

% normalize abundances retrieved from Busi et al, ISME J Comm 2021
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
