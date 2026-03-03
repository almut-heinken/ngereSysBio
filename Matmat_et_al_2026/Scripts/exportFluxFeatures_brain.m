
% first get the total list of reactions with descriptions and subsystems
modelFolder = 'Brain_specific_models';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

rxns = {};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    rxns = union(rxns,model.rxns);
end

Table = {'Reaction_ID','Reaction_Description','Subsystem'};
Table(2:length(rxns)+1,1) = rxns;
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    for j=2:size(Table,1)
        if find(strcmp(Table{j,1},model.rxns))
            rxnID = find(strcmp(model.rxns,Table{j,1}));
            Table{j,2} = model.rxnNames{rxnID,1};
            Table{j,3} = model.subSystems{rxnID,1};
        end
    end
end

models={};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    models{i}=model;
end

samples = readInputTableForPipeline('Data/SampleAnnotations.xlsx');

% get minimal and maximal fluxes
load(['Results_brain_model' filesep 'MinFluxes.mat'])
fluxes = {};
fluxes{1} = 'Reaction';
fluxes(2:length(rxns)+1,1) = rxns;
for k=2:size(samples,1)
    fluxes{1,k} = samples{k,1};
end
for j=2:size(fluxes,1)
    for k=2:size(samples,1)
        modelInd = find(strcmp(modelList,[samples{k,1} '.mat']));
        model=models{modelInd};
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,k} = minFluxes{modelInd}(findRxn,1);
        else
            fluxes{j,k} = 0;
        end
    end
end
writetable(cell2table(fluxes),['Results_brain_model' filesep 'MinFluxes.csv'],'writeVariableNames',false)
save(['Results_brain_model' filesep 'MinFluxesByModel.mat'],'fluxes')

load(['Results_brain_model' filesep 'MaxFluxes.mat'])
fluxes = {};
fluxes{1} = 'Reaction';
fluxes(2:length(rxns)+1,1) = rxns;
for k=2:size(samples,1)
    fluxes{1,k} = samples{k,1};
end
for j=2:size(fluxes,1)
    for k=2:size(samples,1)
        modelInd = find(strcmp(modelList,[samples{k,1} '.mat']));
        model=models{modelInd};
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,k} = maxFluxes{modelInd}(findRxn,1);
        else
            fluxes{j,k} = 0;
        end
    end
end
writetable(cell2table(fluxes),['Results_brain_model' filesep 'MaxFluxes.csv'],'writeVariableNames',false)
save(['Results_brain_model' filesep 'MaxFluxesByModel.mat'],'fluxes')
