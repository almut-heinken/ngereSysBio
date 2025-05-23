
%% get minimal and maximal fluxes-sample-specific models
% combined all into one table
oriModel = readCbModel('iMM1865.xml');

% first get the total list of reactions with descriptions and subsystems

modelFolder = [pwd filesep 'Models'];
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList(find(strncmp(modelList(:,1),'._',2)),:)=[];

rxns = {};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    rxns = union(rxns,model.rxns);
end

% get minimal and maximal fluxes
fluxes = {'Reaction_ID','Reaction_Description','Subsystem'};
fluxes(2:length(rxns)+1,1) = rxns;

cnt=4;
load([pwd filesep 'Metabolic_Flux_Results' filesep 'MinFluxes.mat'])

for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    fluxes{1,cnt} = strrep(modelList{i},'.mat','');
    for j=2:size(fluxes,1)
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,2} = model.rxnNames{findRxn,1};
            fluxes{j,cnt} = minFluxes{i}(findRxn,1);
        else
            fluxes{j,cnt} = 0;
        end
    end
    cnt=cnt+1;
end

% add subsystems
for j=2:size(fluxes,1)
    rxnInd = find(strcmp(oriModel.rxns(:,1),fluxes{j,1}));
    fluxes{j,3}=oriModel.subSystems{rxnInd,1};
end

% remove fluxes that are all zero or all the same
delArray=[];
cnt=1;
for j=2:size(fluxes,1)
    if sum(abs(cell2mat(fluxes(j,4:end))))<0.000001 || all(cell2mat(fluxes(j,4:end)) == fluxes{j,4})
        delArray(cnt)=j;
        cnt=cnt+1;
    end
end
fluxes(delArray,:)=[];
writetable(cell2table(fluxes),[pwd filesep 'Metabolic_Flux_Results' filesep 'MinFluxes.csv'],'writeVariableNames',false)

% max fluxes
fluxes = {'Reaction_ID','Reaction_Description','Subsystem'};
fluxes(2:length(rxns)+1,1) = rxns;
cnt=4;

load([pwd filesep 'Metabolic_Flux_Results' filesep 'MaxFluxes.mat'])
% first get the total list of reactions with descriptions and subsystems
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    fluxes{1,cnt} = strrep(modelList{i},'.mat','');
    for j=2:size(fluxes,1)
        findRxn = find(strcmp(model.rxns,fluxes{j,1}));
        if ~isempty(findRxn)
            fluxes{j,2} = model.rxnNames{findRxn,1};
            fluxes{j,cnt} = maxFluxes{i}(findRxn,1);
        else
            fluxes{j,cnt} = 0;
        end
    end
    cnt=cnt+1;
end
% add subsystems
for j=2:size(fluxes,1)
    rxnInd = find(strcmp(oriModel.rxns(:,1),fluxes{j,1}));
    fluxes{j,3}=oriModel.subSystems{rxnInd,1};
end

% remove fluxes that are all zero or all the same
delArray=[];
cnt=1;
for j=2:size(fluxes,1)
    if sum(abs(cell2mat(fluxes(j,4:end))))<0.000001 || mean(cell2mat(fluxes(j,4:end))) == fluxes{j,4}
        delArray(cnt)=j;
        cnt=cnt+1;
    end
end
fluxes(delArray,:)=[];
writetable(cell2table(fluxes),[pwd filesep 'Metabolic_Flux_Results' filesep 'MaxFluxes.csv'],'writeVariableNames',false)

% combine minimal and maximal fluxes

minFluxes = readInputTableForPipeline([pwd filesep 'Metabolic_Flux_Results' filesep 'MinFluxes.csv']);
for j=2:size(minFluxes,1)
    minFluxes{j,1} = [minFluxes{j,1} '_min'];
end
maxFluxes = readInputTableForPipeline([pwd filesep 'Metabolic_Flux_Results' filesep 'MaxFluxes.csv']);
for j=2:size(maxFluxes,1)
    maxFluxes{j,1} = [maxFluxes{j,1} '_max'];
end
fluxes = minFluxes;
fluxes(size(minFluxes,1)+1:(size(minFluxes,1)+size(maxFluxes,1))-1,:) = maxFluxes(2:end,:);

writetable(cell2table(fluxes),[pwd filesep 'Metabolic_Flux_Results' filesep 'AllFluxes.csv'],'writeVariableNames',false)
