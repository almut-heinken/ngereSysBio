
% Run the following scripts in order to reproduce the simulations of IECM
% phenotypes and patient fibroblasts.

addpath('simulations')

% add one-carbon metabolism/ vitamin B12 metabolism pathways to Recon3D/
% fibroblast reconstruction
prepareReconstructions

% simulate gene KO for IECMs and vitamin B12 deficiency 
simulateGeneKOs_IEMs

mapGeneLociToRecon3D

% Build sample-specific models
% NOTE: This step is time-consuming. It can be skipped by instead using the
% online deposited versions of the models by moving onto the next step.
% Also, the built models may differ slightly from the deposited versions.
createSampleSpecificModels_MMUT

% extract deposited sample-specific models
unzip('OutputModels')

getModelStats_MMUT

interrogateConstrainedModels_MMUT

exportFluxes_MMUT

testMetabolites_MMUT

statisticalAnalysisFluxes_MMUT

