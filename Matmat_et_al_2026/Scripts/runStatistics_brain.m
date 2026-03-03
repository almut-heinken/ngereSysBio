
% perform statistical analysis

mkdir('Statistics_brain_model')
samples = readInputTableForPipeline('Data/SampleAnnotations.xlsx');

% minimal and maximal fluxes
load(['Results_brain_model' filesep 'MinFluxesByModel.mat'])
Statistics = performStatisticalAnalysis(fluxes,samples);
writetable(cell2table(Statistics),['Statistics_brain_model' filesep 'Statistics_min_Fluxes.csv'],'WriteVariableNames',false)

load(['Results_brain_model' filesep 'MaxFluxesByModel.mat'])
Statistics = performStatisticalAnalysis(fluxes,samples);
writetable(cell2table(Statistics),['Statistics_brain_model' filesep 'Statistics_max_Fluxes.csv'],'WriteVariableNames',false)
