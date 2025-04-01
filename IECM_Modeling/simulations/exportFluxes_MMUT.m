
clear all

%% get minimal and maximal fluxes-sample-specific models
% combined all into one table
rxns = {};
% first get the total list of reactions with descriptions and subsystems
modelFolder =[pwd filesep  'OutputModels'];
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

for i=1:length(modelList)
    load([modelFolder filesep modelList{i}]);
    rxns = union(rxns,model.rxns);
end

% get minimal and maximal fluxes
fluxes = {'Reaction_ID','Reaction_Description','Subsystem'};
fluxes(2:length(rxns)+1,1) = rxns;
cnt=4;
load([pwd filesep 'Results' filesep 'MinFluxes.mat'])

for i=1:length(modelList)
    load([modelFolder filesep modelList{i}]);
    fluxes{1,cnt} = strrep(modelList{i},'.mat','');
    for j=2:size(fluxes,1)
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,2} = model.rxnNames{findRxn,1};
            fluxes{j,3} = model.subSystems{findRxn,1};
            fluxes{j,cnt} = minFluxes{i}(findRxn,1);
        else
            fluxes{j,cnt} = 0;
        end
    end
    cnt=cnt+1;
end
writetable(cell2table(fluxes),[pwd filesep 'Results' filesep 'MinFluxes.csv'],'writeVariableNames',false)

fluxes = {'Reaction_ID','Reaction_Description','Subsystem'};
fluxes(2:length(rxns)+1,1) = rxns;
cnt=4;
load([pwd filesep 'Results' filesep 'MaxFluxes.mat'])
% first get the total list of reactions with descriptions and subsystems
for i=1:length(modelList)
    load([modelFolder filesep modelList{i}]);
    fluxes{1,cnt} = strrep(modelList{i},'.mat','');
    for j=2:size(fluxes,1)
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,2} = model.rxnNames{findRxn,1};
            fluxes{j,3} = model.subSystems{findRxn,1};
            fluxes{j,cnt} = maxFluxes{i}(findRxn,1);
        else
            fluxes{j,cnt} = 0;
        end
    end
    cnt=cnt+1;
end
writetable(cell2table(fluxes),[pwd filesep 'Results' filesep 'MaxFluxes.csv'],'writeVariableNames',false)

% combine minimal and maximal fluxes into one table
minFluxes = readInputTableForPipeline([pwd filesep 'Results' filesep 'MinFluxes.csv']);
for i=2:size(minFluxes,1)
    minFluxes{i,1} = [minFluxes{i,1} '_min'];
end
maxFluxes = readInputTableForPipeline([pwd filesep 'Results' filesep 'MaxFluxes.csv']);
for i=2:size(maxFluxes,1)
    maxFluxes{i,1} = [maxFluxes{i,1} '_max'];
end
fluxes = minFluxes;
fluxes(size(minFluxes,1)+1:(size(minFluxes,1)+size(maxFluxes,1))-1,:) = maxFluxes(2:end,:);

Subsystems = fluxes(:,1:3);
fluxes(:,2:3)=[];
writetable(cell2table(fluxes),['Results' filesep 'AllFluxes.csv'],'WriteVariableNames',false)
writetable(cell2table(Subsystems),['Results' filesep 'AllFluxSubsystems.csv'],'WriteVariableNames',false)
