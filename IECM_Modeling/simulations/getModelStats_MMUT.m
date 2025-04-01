
modelFolder = 'OutputModels';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList(find(strncmp(modelList(:,1),'._',2)),:)=[];

% get reaction and metabolite presence
rxnList={};
metList={};

% get reconstruction statistics:
% reactions, metabolites, and genes
stats={};
stats{1,1}='Model_ID';
stats{1,2}='Reactions';
stats{1,3}='Metabolites';
stats{1,4}='Genes';

for i=1:length(modelList)
    stats{i+1,1} = strrep(modelList{i},'.mat','');
    model = readCbModel([modelFolder filesep modelList{i}]);

    % get reactions and metabolites
    rxnList = union(rxnList,model.rxns);
    metList = union(metList,model.mets);

    % Number of reactions, metabolites, and genes
    stats{i+1,2}=length(model.rxns);
    stats{i+1,3}=length(model.mets);
    stats{i+1,4}=length(model.genes);
end
save(['Results' filesep 'Model_statistics'],'stats')

% statistical analysis
samples = readInputTableForPipeline(['input' filesep 'Metadata.csv']);

% Control vs. MMA
col = find(strcmp(samples(1,:),'Disease_state'));
[Statistics,significantFeatures] = performStatisticalAnalysis(stats',samples,'stratification',samples{1,col});

% Control vs. MUT0
% remove all other cases
statsRed = stats';
otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
otherSamples = samples(find(strcmp(samples(:,col),'MUT-')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
col = find(strcmp(samples(1,:),'mut_category'));
[Statistics,significantFeatures] = performStatisticalAnalysis(statsRed,samples,'stratification',samples{1,col});

% Control vs. MUT-
% remove all other cases
statsRed = stats';
otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
otherSamples = samples(find(strcmp(samples(:,col),'MUT0')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
col = find(strcmp(samples(1,:),'mut_category'));
[Statistics,significantFeatures] = performStatisticalAnalysis(statsRed,samples,'stratification',samples{1,col});

% Control vs. other MMA
col = find(strcmp(samples(1,:),'Group'));
statsRed=stats';
% remove all other cases
otherSamples = samples(find(strcmp(samples(:,col),'MUT-')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
otherSamples = samples(find(strcmp(samples(:,col),'MUT0')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(statsRed,samples,'stratification',samples{1,col});

% MUT0 vs. MUT-
col = find(strcmp(samples(1,:),'mut_category'));
statsRed=stats';
% remove all other cases
otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
[C,I]=intersect(statsRed(1,:),otherSamples);
statsRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(statsRed,samples,'stratification',samples{1,col});


Table2a = {'','Reactions','Non-unique metabolites','Genes'};

Table2a{2,1} = 'All';
for i=2:size(stats,2)
    Table2a{2,i} = [num2str(mean(cell2mat(stats(2:end,i)))) ' +/-' num2str(std(cell2mat(stats(2:end,i))))];
end

% count separately for groups
col = find(strcmp(samples(1,:),'Group'));

statsRed = stats;
otherSamples = samples(find(~strcmp(samples(:,col),'Control')),1);
[C,I]=intersect(statsRed(:,1),otherSamples,'stable');
statsRed(I(2:end),:)=[];
Table2a{3,1} = 'Controls';
for i=2:size(statsRed,2)
    Table2a{3,i} = [num2str(mean(cell2mat(statsRed(2:end,i)))) ' +/-' num2str(std(cell2mat(statsRed(2:end,i))))];
end

statsRed = stats;
otherSamples = samples(find(~strcmp(samples(:,col),'MUT0')),1);
[C,I]=intersect(statsRed(:,1),otherSamples,'stable');
statsRed(I(2:end),:)=[];
Table2a{4,1} = 'MUT0';
for i=2:size(statsRed,2)
    Table2a{4,i} = [num2str(mean(cell2mat(statsRed(2:end,i)))) ' +/-' num2str(std(cell2mat(statsRed(2:end,i))))];
end

statsRed = stats;
otherSamples = samples(find(~strcmp(samples(:,col),'MUT-')),1);
[C,I]=intersect(statsRed(:,1),otherSamples,'stable');
statsRed(I(2:end),:)=[];
Table2a{5,1} = 'MUT-';
for i=2:size(statsRed,2)
    Table2a{5,i} = [num2str(mean(cell2mat(statsRed(2:end,i)))) ' +/-' num2str(std(cell2mat(statsRed(2:end,i))))];
end

statsRed = stats;
otherSamples = samples(find(~strcmp(samples(:,col),'MMA_Other')),1);
[C,I]=intersect(statsRed(:,1),otherSamples,'stable');
statsRed(I(2:end),:)=[];
Table2a{6,1} = 'MMA (other)';
for i=2:size(statsRed,2)
    Table2a{6,i} = [num2str(mean(cell2mat(statsRed(2:end,i)))) ' +/-' num2str(std(cell2mat(statsRed(2:end,i))))];
end
