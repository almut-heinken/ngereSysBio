
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

addpath('Scripts')
addpath([pwd filesep 'input'])

mkdir([pwd filesep 'MicrobiomeModels'])
% run mgPipe workflow

% number of cores dedicated to parallelization
numWorkers = 4;

% define a folder where results will be saved
resPath = [pwd filesep 'MicrobiomeModels'];

% use pan-species models
panPath = [pwd filesep 'panSpeciesModels'];

% path to and name of the file with abundance information.
abunFilePath = [pwd filesep 'normalizedCoverage.csv'];

%% custom "Spanish" diets
% path to the simulated diet
dietFilePath = [pwd filesep 'input' filesep 'InputDiets_DIFAMEM.txt'];

% run model building and simulations
[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(panPath, abunFilePath, true, 'resPath', resPath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);

% exclude metabolites that are not produced or are all the same
mkdir([pwd filesep 'MicrobiomeResults'])
delArray=[];
cnt=1;
for i=2:size(netSecretionFluxes,1)
    if abs(sum(cell2mat(netSecretionFluxes(i,2:end))))<0.000001 || all(round(cell2mat(netSecretionFluxes(i,2:end)),3) == round(netSecretionFluxes{i,2},3))
        delArray(cnt)=i;
        cnt=cnt+1;
    end
end
netSecretionFluxes(delArray,:)=[];
writetable(cell2table(netSecretionFluxes),[pwd filesep 'MicrobiomeResults' filesep 'CD_netSecretionFluxes.csv'],'writeVariableNames',false);

% export matched metabolite names for plots
database = loadVMHDatabase;
metNames = netSecretionFluxes(:,1);
metNames(1,2:3) = {'Metabolite name','Subsystem'};
for i=2:size(metNames,1)
    met = strrep(metNames{i,1},'EX_','');
    met = strrep(met,'[fe]','');
    metNames{i,2} = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};
    metNames{i,3} = database.metabolites{find(strcmp(database.metabolites(:,1),met)),13};
end
writetable(cell2table(metNames),[pwd filesep 'MicrobiomeResults' filesep 'CD_metNames.csv'],'writeVariableNames',false);

