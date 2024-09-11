
% statistical analysis for the groups
oriModel = readCbModel('iMM1865.xml');

database=loadVMHDatabase;

samples = readInputTableForPipeline([pwd filesep 'Data' filesep 'Samples.csv']);
samples(:,2) = strrep(samples(:,2),'hyperoxia','HO');
samples(:,2) = strrep(samples(:,2),'normoxia','NO');

mkdir([pwd filesep 'Metabolic_Flux_Statistics'])
 
%% perform statistical analysis for lipidomics

data = readInputTableForPipeline([pwd filesep 'Data' filesep 'Plasma_lipidomics.csv']);
% summarize duplicate entries
[uniqueA,i,j] = unique(data(:,1));
n  = accumarray(j(:),1);
Dupes=uniqueA(find(n>1));
delArray=[];
cnt=1;
for i=1:length(Dupes)
    indexToDupes = find(strcmp(data(:,1),Dupes{i}));
    for j=2:length(indexToDupes)
        for k=2:size(data,2)
            data{indexToDupes(1),k}=data{indexToDupes(1),k}+data{indexToDupes(j),k};
            delArray(cnt,1)=indexToDupes(j);
            cnt=cnt+1;
        end
    end
end
data(delArray,:)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(data,samples);

writetable(cell2table(Statistics),[pwd filesep 'Statistics_lipidomics.csv'],'WriteVariableNames',false)
writetable(cell2table(significantFeatures),[pwd filesep 'significantFeatures_lipidomics.csv'],'WriteVariableNames',false)

%% compare fluxes between groups

mkdir([pwd filesep 'Figure_S3'])

Fluxes = readInputTableForPipeline([pwd filesep 'Metabolic_Flux_Results' filesep 'AllFluxes.csv']);
subsystems = Fluxes(:,1:3);
for j=2:size(subsystems,1)
    if strncmp(subsystems{j,3},'Transport',9)
        subsystems{j,3}='Transport';
    end
end
Fluxes(:,2:3)=[];

subs = unique(subsystems(2:end,3));
subs(find(strcmp(subs,'')))=[];

table = {'Reactions by subsystem','Higher in BALBcByJ HO vs. NO','Higher in BALBcByJ NO vs. HO','Higher in C57BL6J HO vs. NO','Higher in C57BL6J NO vs. HO','Higher in BALBcByJ NO vs. C57BL6J NO','Higher in C57BL6J NO vs. BALBcByJ NO','Higher in BALBcByJ HO vs. C57BL6J HO','Higher in C57BL6J HO vs. BALBcByJ HO'};
table(2:length(subs)+1,1) = subs;
table(2:end,2:end) = {0};

% BALBcByJ HO vs. NO
fluxesRed = Fluxes;
remSamples = samples(find(contains(samples(:,2),'C57BL6J')),1);
[C,I] = intersect(fluxesRed(1,:),remSamples);
fluxesRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(fluxesRed,samples);
% add reaction descriptions and subsystems
Statistics(1,12:14) = {'VMH_ID','Description','Subsystem'};
for i=2:size(Statistics,1)
    Statistics{i,12} = Statistics{i,1};
    findRxn = find(strcmp(subsystems(:,1),Statistics{i,1}));
    Statistics{i,13} = subsystems{findRxn,2};
    Statistics{i,14} = subsystems{findRxn,3};
end
Statistics(:,12) = strrep(Statistics(:,12),'_min','');
Statistics(:,12) = strrep(Statistics(:,12),'_max','');
writetable(cell2table(Statistics),[pwd filesep 'Metabolic_Flux_Statistics' filesep 'BALBcByJ_HO_vs_NO.csv'],'WriteVariableNames',false)

Statistics(1,:)=[];
% remove results that are not initially significant
Statistics(find(cell2mat(Statistics(:,3))>=0.05),:)=[];

for j=2:size(Statistics,1)
    findSub = subsystems{find(strcmp(subsystems(:,1),Statistics{j,1})),3};
    if ~isempty(findSub)
        if abs(Statistics{j,8})>abs(Statistics{j,10})
            table{find(strcmp(table(:,1),findSub)),2} = table{find(strcmp(table(:,1),findSub)),2}+1;
        elseif abs(Statistics{j,10})>abs(Statistics{j,8})
            table{find(strcmp(table(:,1),findSub)),3} = table{find(strcmp(table(:,1),findSub)),3}+1;
        end
    end
end

% C57BL6J HO vs. NO
fluxesRed = Fluxes;
remSamples = samples(find(contains(samples(:,2),'BALBcByJ')),1);
[C,I] = intersect(fluxesRed(1,:),remSamples);
fluxesRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(fluxesRed,samples);
% add reaction descriptions and subsystems
Statistics(1,12:14) = {'VMH_ID','Description','Subsystem'};
for i=2:size(Statistics,1)
    Statistics{i,12} = Statistics{i,1};
    findRxn = find(strcmp(subsystems(:,1),Statistics{i,1}));
    Statistics{i,13} = subsystems{findRxn,2};
    Statistics{i,14} = subsystems{findRxn,3};
end
Statistics(:,12) = strrep(Statistics(:,12),'_min','');
Statistics(:,12) = strrep(Statistics(:,12),'_max','');
writetable(cell2table(Statistics),[pwd filesep 'Metabolic_Flux_Statistics' filesep 'C57BL6J_HO_vs_NO.csv'],'WriteVariableNames',false)

Statistics(1,:)=[];
% remove results that are not initially significant
Statistics(find(cell2mat(Statistics(:,3))>=0.05),:)=[];

for j=2:size(Statistics,1)
    findSub = subsystems{find(strcmp(subsystems(:,1),Statistics{j,1})),3};
    if ~isempty(findSub)
        if abs(Statistics{j,8})>abs(Statistics{j,10})
            table{find(strcmp(table(:,1),findSub)),4} = table{find(strcmp(table(:,1),findSub)),4}+1;
        elseif abs(Statistics{j,10})>abs(Statistics{j,8})
            table{find(strcmp(table(:,1),findSub)),5} = table{find(strcmp(table(:,1),findSub)),5}+1;
        end
    end
end

% BALBcByJ NO vs. C57BL6J NO
fluxesRed = Fluxes;
remSamples = samples(find(contains(samples(:,2),' HO')),1);
[C,I] = intersect(fluxesRed(1,:),remSamples);
fluxesRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(fluxesRed,samples);
% add reaction descriptions and subsystems
Statistics(1,12:14) = {'VMH_ID','Description','Subsystem'};
for i=2:size(Statistics,1)
    Statistics{i,12} = Statistics{i,1};
    findRxn = find(strcmp(subsystems(:,1),Statistics{i,1}));
    Statistics{i,13} = subsystems{findRxn,2};
    Statistics{i,14} = subsystems{findRxn,3};
end
Statistics(:,12) = strrep(Statistics(:,12),'_min','');
Statistics(:,12) = strrep(Statistics(:,12),'_max','');
writetable(cell2table(Statistics),[pwd filesep 'Metabolic_Flux_Statistics' filesep 'BALBcByJ_NO_vs_C57BL6J_NO.csv'],'WriteVariableNames',false)

Statistics(1,:)=[];
% remove results that are not initially significant
Statistics(find(cell2mat(Statistics(:,3))>=0.05),:)=[];

for j=2:size(Statistics,1)
    findSub = subsystems{find(strcmp(subsystems(:,1),Statistics{j,1})),3};
    if ~isempty(findSub)
        if abs(Statistics{j,8})>abs(Statistics{j,10})
            table{find(strcmp(table(:,1),findSub)),6} = table{find(strcmp(table(:,1),findSub)),6}+1;
        elseif abs(Statistics{j,10})>abs(Statistics{j,8})
            table{find(strcmp(table(:,1),findSub)),7} = table{find(strcmp(table(:,1),findSub)),7}+1;
        end
    end
end

% BALBcByJ HO vs. C57BL6J HO
fluxesRed = Fluxes;
remSamples = samples(find(contains(samples(:,2),' NO')),1);
[C,I] = intersect(fluxesRed(1,:),remSamples);
fluxesRed(:,I)=[];
[Statistics,significantFeatures] = performStatisticalAnalysis(fluxesRed,samples);
% add reaction descriptions and subsystems
Statistics(1,12:14) = {'VMH_ID','Description','Subsystem'};
for i=2:size(Statistics,1)
    Statistics{i,12} = Statistics{i,1};
    findRxn = find(strcmp(subsystems(:,1),Statistics{i,1}));
    Statistics{i,13} = subsystems{findRxn,2};
    Statistics{i,14} = subsystems{findRxn,3};
end
Statistics(:,12) = strrep(Statistics(:,12),'_min','');
Statistics(:,12) = strrep(Statistics(:,12),'_max','');
writetable(cell2table(Statistics),[pwd filesep 'Metabolic_Flux_Statistics' filesep 'BALBcByJ_HO_vs_C57BL6J_HO.csv'],'WriteVariableNames',false)

Statistics(1,:)=[];
% remove results that are not initially significant
Statistics(find(cell2mat(Statistics(:,3))>=0.05),:)=[];

for j=2:size(Statistics,1)
    findSub = subsystems{find(strcmp(subsystems(:,1),Statistics{j,1})),3};
    if ~isempty(findSub)
        if abs(Statistics{j,8})>abs(Statistics{j,10})
            table{find(strcmp(table(:,1),findSub)),8} = table{find(strcmp(table(:,1),findSub)),8}+1;
        elseif abs(Statistics{j,10})>abs(Statistics{j,8})
            table{find(strcmp(table(:,1),findSub)),9} = table{find(strcmp(table(:,1),findSub)),9}+1;
        end
    end
end

% remove zero columns
delArray=[];
cnt=1;
for j=2:size(table,1)
    if abs(sum(cell2mat(table(j,2:end))))==0
        delArray(cnt)=j;
        cnt=cnt+1;
    end
end
table(delArray,:)=[];

% plot rank ordered for each case
for j=2:size(table,2)
    
    pies = table(2:end,1);
    data = cell2mat(table(2:end,j));
    pies(find(data==0),:)=[];
    data(find(data==0),:)=[];
    [C,I]=sort(data,'descend');
    data=data(I);
    pies=pies(I);
    % reduce number of categories
    if length(data)>12
        summed=sum(data(12:end));
        data(12:end)=[];
        data(12,1)=summed;
        pies(12:end)=[];
        pies{12,1}='Others';
    end
    data=data';
    pies=pies';
    
    % define colors
    cols = hsv(length(data));
    
    f=figure('Renderer', 'painters', 'Position', [10 100 900 600]);
    donut(data,pies,cols,'pie')
    grid off
    axis off
    legend('Location','eastoutside')
    title(table{1,j})
    set(gca, 'FontSize', 14)
    hold on
    f.Children(1).Position=[0.5296    0.3196    0.4156    0.3958];
    f.Children(2).Position=[0.0744    0.1100    0.4438    0.8150];
    filename = strrep(table{1,j},' ','_');
    filename = strrep(filename,'.','');
    print([pwd filesep 'Figure_S3' filesep filename],'-dpng','-r300')
end
close all
