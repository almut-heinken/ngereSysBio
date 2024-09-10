

% Run microbiome modeling pipeline consisting of creating an input file
% with relative strain abundances, creation of microbiome models,
% interrogation of microbiome models through simulations, and plotting of
% computed metabolite uptake and secretion by the communities.

mkdir([rootDir filesep 'Modeling_COSMIC'])
% run mgPipe workflow

% number of cores dedicated to parallelization
numWorkers = 4;

% define a folder where results will be saved
resPath = [rootDir filesep 'Modeling_COSMIC' filesep 'MicrobiomeModels'];

modPath = [rootDir filesep 'AGORA2_gapfilled'];

cd([rootDir filesep 'inputFiles'])
% create pan-species models
panPath = [rootDir filesep 'panSpeciesModels'];
createPanModels(modPath,panPath,'Species',numWorkers,[rootDir filesep 'inputFiles' filesep 'expanded_AGORA2_infoFile.xlsx'])
cd(rootDir)

% test pan-species models
[notGrowing,biomassFluxes] = plotBiomassTestResults(panPath, 'Species', 'numWorkers',numWorkers);
[tooHighATP,atpFluxes] = plotATPTestResults(panPath, 'Species', 'numWorkers',numWorkers);

% path to and name of the file with abundance information.
abunFilePath = [rootDir filesep 'inputFiles' filesep 'normalizedCoverage_COSMIC.csv'];

% path to the file with characteristics of the study participants
infoFilePath = [rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv'];

% path to the simulated diet
dietFilePath = [rootDir filesep 'inputFiles' filesep 'inputDiets_COSMIC.txt'];

% lower the required minimum biomass production
lowerBMBound = 0.2;

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(panPath, abunFilePath, true, 'resPath', resPath, 'infoFilePath', infoFilePath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers, 'lowerBMBound', lowerBMBound);
writetable(cell2table(netSecretionFluxes),[resPath filesep 'Infant_diet.csv'],'writeVariableNames',false)

%% determine taxon-metabolite correlations
mkdir([rootDir filesep 'Modeling_COSMIC' filesep 'Correlations'])
taxInfo = [rootDir filesep 'inputFiles' filesep 'AGORA2_infoFile.xlsx'];
fluxPath = [rootDir filesep 'Modeling_COSMIC' filesep 'MicrobiomeModels' filesep 'Infant_diet.csv'];
corrMethod = 'Spearman';
% first with net secretion
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, taxInfo, corrMethod);

% export the results
taxa=fieldnames(FluxCorrelations);
for i=1:length(taxa)
    cnt=1;
    delArray=[];
    for j=2:size(FluxCorrelations.(taxa{i}),2)
        if ~any(abs(cell2mat(FluxCorrelations.(taxa{i})(2:end,j))) > 0.8)
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    FluxCorrelations.(taxa{i})(:,delArray)=[];

    cnt=1;
    delArray=[];
    for j=2:size(FluxCorrelations.(taxa{i}),1)
        if ~any(abs(cell2mat(FluxCorrelations.(taxa{i})(j,2:end))) > 0.8)
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    FluxCorrelations.(taxa{i})(delArray,:)=[];
    if i>2
        [C,I]=setdiff(TaxonomyInfo.(taxa{i})(:,1),FluxCorrelations.(taxa{i})(2:end,1),'stable');
        TaxonomyInfo.(taxa{i})(I(2:end),:)=[];
    end

    FluxCorrelations.(taxa{i})(find(strcmp(FluxCorrelations.(taxa{i})(:,1),'')),:)=[];
    writetable(cell2table(FluxCorrelations.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'FluxCorrelations_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    writetable(cell2table(PValues.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'PValues_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    if i>1
        writetable(cell2table(TaxonomyInfo.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'TaxonomyInfo_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    end
end

%% for comparison with adults, compute fluxes on the standard Western diet

dietFilePath = 'AverageEuropeanDiet.txt';

delete([resPath filesep 'simRes.mat'],[resPath filesep 'intRes.mat'])

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(panPath, abunFilePath, true, 'resPath', resPath, 'infoFilePath', infoFilePath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);
writetable(cell2table(netSecretionFluxes),[resPath filesep 'AE_diet.csv'],'writeVariableNames',false)

%% compute microbe secretion for metabolites of interest
% vitamins, SCFAs
mets = {'ac','ppa','but','lac_L','fol','adocbl','ribflv','thm','nac','btn','pydx','pnto_R'};
constrModPath = [resPath filesep 'Diet'];
contrPath = [rootDir filesep 'Modeling_COSMIC' filesep 'Contributions'];

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
save([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'minFluxes.mat'],'minFluxes')
cell2csv([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)

% minFluxes = secretion
minFluxes=readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv']);
for i=2:size(minFluxes,1)
    for j=2:size(minFluxes,2)
        if minFluxes{i,j} < 0
            minFluxes{i,j} = abs(minFluxes{i,j});
        elseif minFluxes{i,j} > 0
            minFluxes{i,j} = 0;
        end
    end
end
cell2csv([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)
