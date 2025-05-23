
% interrogate the RNA-seq constrained models:
% first get get the model statistics and distribution of reactions,
% metabolites, and genes
% then perform flux variability analysis to get the allowed flux through all reactions
% in the sample-specific models

clear all

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

numWorkers = 4;

% start parallel pool
if ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

numWorkers=4;

mkdir('Results')

modelFolder = 'OutputModels';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList(find(strncmp(modelList(:,1),'._',2)),:)=[];

% sample flux spaces
% flux={};
minFluxes = {};
maxFluxes = {};

info = readInputTableForPipeline(['input' filesep 'Metadata.csv']);
dCol = find(strcmp(info(1,:),'Disease_state'));
mutCol = find(strcmp(info(1,:),'mut_category'));
geneCol = find(strcmp(info(1,:),'GeneVariant'));

for i=1:length(modelList)
    i
    model = readCbModel([modelFolder filesep modelList{i}]);
    
    % silence reactions-should not have flux
    drugRxns = {'ACMPthc','ACMPdt','MERACMPtep','MERACMPthc','CRVSM1hr','CVM1GLUChc','ACMPGLUTtep','ACMPGLUTthc','ACMPGLUtep','ACMPGLUthc','CRVSM1hr','CVM1GLUChc','ATVLACGLCURhc','ATVLAChr'};
    model = changeRxnBounds(model,drugRxns,0,'b');
    
    % implement medium
    model=defineMediumDMEM(model);
    
    model = changeRxnBounds(model,'DM_dna[n]',0.001,'l');
    model = changeRxnBounds(model,'DM_dna5mtc[n]',0.001,'l');
    
    % enforce biomass production
    model = changeRxnBounds(model,'biomass_maintenance',0.01,'l');
    
    % check growth
    model = changeObjective(model,'biomass_maintenance');
    
    % implement gene/protein defects
    findSamp = find(strcmp(info(:,1),strrep(modelList{i},'.mat','')));
    if strcmp(info{findSamp,dCol},'MMA')
        if strcmp(info{findSamp,mutCol},'MUT0')
            if ~isempty(find(contains(model.genes,'4594.1')))
                [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '4594.1',0);
            end
        elseif strcmp(info{findSamp,mutCol},'MUT-')
            if ~isempty(find(contains(model.genes,'4594.1')))
                % assume a bit of remaining flux-20%
                [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '4594.1',0.2);
            end
        elseif strcmp(info{findSamp,mutCol},'NA')
            % account for other variants if possible
            if strcmp(info{findSamp,geneCol},'ACSF3')
                if ~isempty(find(contains(model.genes,'197322.1')))
                    [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '197322.1',0);
                end
            elseif strcmp(info{findSamp,geneCol},'MMAA')
                if ~isempty(find(contains(model.genes,'166785.1')))
                    [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '166785.1',0);
                end
            elseif strcmp(info{findSamp,geneCol},'MMAB')
                if ~isempty(find(contains(model.genes,'326625.1')))
                    [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '326625.1',0);
                end
            elseif strcmp(info{findSamp,geneCol},'SUCLA2')
                if ~isempty(find(contains(model.genes,'8803.1')))
                    [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '8803.1',0);
                end
            elseif strcmp(info{findSamp,geneCol},'TCN2')
                if ~isempty(find(contains(model.genes,'197322.1')))
                    [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '197322.1',0);
                end
            end
        end
    end
    
    % compute minimal and maximal fluxes
    [minFlux,maxFlux] = fastFVA(model,0,'max','ibm_cplex',model.rxns, 'S');
    minFluxes{i} = minFlux;
    maxFluxes{i} = maxFlux;
end

save(['Results' filesep 'MinFluxes'],'minFluxes')
save(['Results' filesep 'MaxFluxes'],'maxFluxes')
