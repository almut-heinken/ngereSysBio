
% Predict each model's potential to synthesize all internal metabolites.

clear all

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox
end
solver = CBT_LP_SOLVER;

numWorkers=12;
if numWorkers>0 && ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

info = readInputTableForPipeline(['input' filesep 'Metadata.csv']);
dCol = find(strcmp(info(1,:),'Disease_state'));
mutCol = find(strcmp(info(1,:),'mut_category'));
geneCol = find(strcmp(info(1,:),'GeneVariant'));

modelFolder = 'OutputModels';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList(find(strncmp(modelList(:,1),'._',2)),:)=[];

% get all metabolites
mets = {};
for i=1:length(modelList)
    load([modelFolder filesep modelList{i}]);
    mets = union(mets,model.mets);
end

% go through objectives
Table = [{''} 
    mets];
% loop through models
for i=1:length(modelList)
    Table{1,i+1} = strrep(modelList{i},'.mat','');
    model = readCbModel([modelFolder filesep modelList{i}]);
    Table(2:end,i+1) = {0};
    
    % implement medium
    model=defineMediumDMEM(model);
    
    model = changeRxnBounds(model,'DM_dna[n]',0.001,'l');
    model = changeRxnBounds(model,'DM_dna5mtc[n]',0.001,'l');
    
    % enforce biomass production
    model = changeRxnBounds(model,'biomass_maintenance',0.01,'l');
    
    % implement gene/protein defects
    findSamp = find(strcmp(info(:,1),strrep(modelList{i},'.mat','')));
    if strcmp(info{findSamp,dCol},'MMA')
        if strcmp(info{findSamp,mutCol},'MUT0')
            if ~isempty(find(contains(model.genes,'4594.1')))
                [model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, '4594.1',0);
            end
        elseif strcmp(info{findSamp,mutCol},'MUT-')
            if ~isempty(find(contains(model.genes,'4594.1')))
                % assume a bit of remaining flux
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
    
    metsTemp = {};
    parfor m=1:length(model.mets)
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        
        modelMet = model;
        modelMet = addSinkReactions(modelMet,model.mets{m});
        modelMet = changeObjective(modelMet,['sink_' model.mets{m}]);
        solution = optimizeCbModel(modelMet, 'max');
        metsTemp{m} = solution.f;
    end
    
    for m=1:length(model.mets)
        findMet = find(strcmp(Table(:,1),model.mets{m}));
        Table{findMet,i+1} =  metsTemp{m};
    end
end
% remove fluxes that are zero or all the same
delArray=[];
cnt=1;
for j=2:size(Table,1)
    if abs(sum(cell2mat(Table(j,2:end))))<0.00001 || all(cell2mat(Table(j,2:end)) == Table{j,2})
        delArray(cnt)=j;
        cnt=cnt+1;
    end
end
Table(delArray,:) = [];
writetable(cell2table(Table),[pwd filesep 'Results' filesep 'MetaboliteFluxes.csv'],'writeVariableNames',false)
