
% Run microbiome modeling pipeline consisting of creating an input file
% with relative strain abundances, creation of microbiome models,
% interrogation of microbiome models through simulations, and plotting of
% computed metabolite uptake and secretion by the communities.

%% Gut microbiome

mkdir([rootDir filesep 'Modeling_MaternalGutMicrobiomes'])

%% run mgPipe workflow
% define a folder where results will be saved
resPath = [rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'MicrobiomeModels'];

modPath = [rootDir filesep 'panSpeciesModels'];

% path to and name of the file with abundance information.
abunFilePath = [rootDir filesep 'inputFiles' filesep 'normalizedCoverage_MaternalGutMicrobiome.csv'];

% path to the simulated diet
dietFilePath = 'AverageEuropeanDiet.txt';

% number of cores dedicated to parallelization
numWorkers = 4;

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(modPath, abunFilePath, true, 'resPath', resPath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);
writetable(cell2table(netSecretionFluxes),[resPath filesep 'AE_diet.csv'],'writeVariableNames',false)

%% compute microbe secretion for metabolites of interest
% vitamins, SCFAs
mets = {'ac','ppa','but','lac_L','fol','adocbl','ribflv','thm','nac','btn','pydx','pnto_R'};
constrModPath = [resPath filesep 'Diet'];
contrPath = [rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions'];

% make changes to the model to prevent false prediction of riboflavin
dInfo = dir(constrModPath);
fileList={dInfo.name};
fileList=fileList';
fileList(~contains(fileList(:,1),'.mat'),:)=[];
for i=1:length(fileList)
    load([constrModPath filesep fileList{i}])
    rxns=model.rxns(find(contains(model.rxns,'IEX_rbflvrd[u]tr')));
    model=changeRxnBounds(model,rxns,0,'u');
    save([constrModPath filesep fileList{i}],'model')
end

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath, 'metList', mets, 'numWorkers', numWorkers,'resultsFolder',contrPath);
save([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'minFluxes.mat'],'minFluxes')
cell2csv([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)

% minFluxes = secretion
minFluxes=readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Microbe_Secretion.csv']);
for i=2:size(minFluxes,1)
    for j=2:size(minFluxes,2)
        if minFluxes{i,j} < 0
            minFluxes{i,j} = abs(minFluxes{i,j});
        elseif minFluxes{i,j} > 0
            minFluxes{i,j} = 0;
        end
    end
end
cell2csv([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)
