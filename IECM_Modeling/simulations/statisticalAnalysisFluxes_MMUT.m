
clear all

%% Reproduce statistical analyses in the study
% (i) Computed metabolic fluxes, (ii) predcited internally produced
% metabolites
%%
samples = readInputTableForPipeline(['input' filesep 'Metadata.csv']);
model = readCbModel('fibroblast_CBL.mat');

datasets = {'Fluxes','Metabolites'};
for d=1:length(datasets)
    mkdir(['Statistics_' datasets{d}])
    
    if strcmp(datasets{d},'Fluxes')
        data = readInputTableForPipeline(['Results' filesep 'AllFluxes.csv']);
        Subsystems = readInputTableForPipeline(['Results' filesep 'AllFluxSubsystems.csv']);
        Subsystems(find(strncmp(Subsystems(:,3),'Transport',9)),3) = {'Transport'};
        Subsystems(find(strncmp(Subsystems(:,3),'Exchange',8)),3) = {'Exchange/demand reaction'};
        subsystems = unique(Subsystems(2:end,3));
    elseif strcmp(datasets{d},'Metabolites')
        data = readInputTableForPipeline(['Results' filesep 'MetaboliteFluxes.csv']);
    end
    
    % remove fluxes that are zero or all the same
    delArray=[];
    for i=2:size(data,1)
        if abs(sum(cell2mat(data(i,2:end))))<0.00001 || length(find(abs(cell2mat(data(i,2:end)) - mean(cell2mat(data(i,2:end))))<0.00001))==length(data(i,2:end))
            delArray(length(delArray)+1)=i;
        end
    end
    data(delArray,:)=[];

    %% Statistical analyses comparing fluxes in mut0/mut-/other MMA/controls
    
    % create overview tables
    table_Sig = vertcat({''},subsystems);
    colCnt=2;
    statCol = 4;
    sigRxns = {};
    
    % Control vs. mut0
    col = find(strcmp(samples(1,:),'mut_category'));
    dataAnalysis=data;
    % remove all other cases
    otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    otherSamples = samples(find(strcmp(samples(:,col),'mut-')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    % remove zero fluxes
    delArray=[];
    cnt=1;
    for j=2:size(dataAnalysis,1)
        if abs(sum(cell2mat(dataAnalysis(j,2:end))))<0.00001
            delArray(cnt)=j;
            cnt=cnt+1;
        end
    end
    dataAnalysis(delArray,:)=[];
    [Statistics,significantFeatures] = performStatisticalAnalysis(dataAnalysis,samples,'stratification',samples{1,col});
    sigRxns = union(sigRxns,significantFeatures(2:end,1));
    
    % add reaction/ metabolites names
    if strcmp(datasets{d},'Fluxes')
        for i=2:size(Statistics,1)
            findRxn = find(strcmp(Subsystems(:,1),Statistics{i,1}));
            Statistics{i,2} = Subsystems{findRxn,2};
        end
    elseif strcmp(datasets{d},'Metabolites')
        for i=2:size(Statistics,1)
            met = model.metNames{find(strcmp(model.mets,dataAnalysis{i,1}))};
            Statistics{i,2} = met;
        end
    end
    writetable(cell2table(Statistics),['Statistics_' datasets{d} filesep 'Statistics_Control_MUT0.csv'],'WriteVariableNames',false)
    
    if strcmp(datasets{d},'Fluxes')
        % get enriched reactions-significant
        diffReactions={};
        diffReactions{1} = subsystems;
        diffReactions{1}(:,2) = {0};
        diffReactions{2} = subsystems;
        diffReactions{2}(:,2) = {0};
        diffReactions{3} = subsystems;
        diffReactions{3}(:,2) = {0};
        diffReactions{4} = subsystems;
        diffReactions{4}(:,2) = {0};
        
        % get enriched reactions-significant and initially significant
        for i=2:size(Statistics,1)
            if Statistics{i,statCol} < 0.05
                sub = Subsystems{find(strcmp(Subsystems(:,1),Statistics{i,1})),3};
                findSub = find(strcmp(subsystems(:,1),sub));
                if abs(Statistics{i,8})>abs(Statistics{i,10})
                    diffReactions{1}{findSub,2} = diffReactions{1}{findSub,2}+1;
                elseif abs(Statistics{i,10})>abs(Statistics{i,8})
                    diffReactions{2}{findSub,2} = diffReactions{2}{findSub,2}+1;
                end
            end
        end
        
        % fill in tables
        table_Sig{1,colCnt} = 'Higher in controls vs. mut0';
        table_Sig(2:end,colCnt) = diffReactions{1}(:,2);
        table_Sig{1,colCnt+1} = 'Higher in mut0 vs. controls';
        table_Sig(2:end,colCnt+1) = diffReactions{2}(:,2);
        colCnt = colCnt+2;
    end
    
    % Control vs. mut-
    col = find(strcmp(samples(1,:),'mut_category'));
    dataAnalysis=data;
    % remove all other cases
    otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    otherSamples = samples(find(strcmp(samples(:,col),'mut0')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    % remove zero fluxes
    delArray=[];
    cnt=1;
    for j=2:size(dataAnalysis,1)
        if abs(sum(cell2mat(dataAnalysis(j,2:end))))<0.00001
            delArray(cnt)=j;
            cnt=cnt+1;
        end
    end
    dataAnalysis(delArray,:)=[];
    [Statistics,significantFeatures] = performStatisticalAnalysis(dataAnalysis,samples,'stratification',samples{1,col});
    sigRxns = union(sigRxns,significantFeatures(2:end,1));
    
    % add reaction/ metabolites names
    if strcmp(datasets{d},'Fluxes')
        for i=2:size(Statistics,1)
            findRxn = find(strcmp(Subsystems(:,1),Statistics{i,1}));
            Statistics{i,2} = Subsystems{findRxn,2};
        end
    elseif strcmp(datasets{d},'Metabolites')
        for i=2:size(Statistics,1)
            met = model.metNames{find(strcmp(model.mets,dataAnalysis{i,1}))};
            Statistics{i,2} = met;
        end
    end
    writetable(cell2table(Statistics),['Statistics_' datasets{d} filesep 'Statistics_Control_MUT-.csv'],'WriteVariableNames',false)
    
    if strcmp(datasets{d},'Fluxes')
        % get enriched reactions-significant
        diffReactions={};
        diffReactions{1} = subsystems;
        diffReactions{1}(:,2) = {0};
        diffReactions{2} = subsystems;
        diffReactions{2}(:,2) = {0};
        diffReactions{3} = subsystems;
        diffReactions{3}(:,2) = {0};
        diffReactions{4} = subsystems;
        diffReactions{4}(:,2) = {0};
        
        % get enriched reactions-significant and initially significant
        for i=2:size(Statistics,1)
            if Statistics{i,statCol} < 0.05
                sub = Subsystems{find(strcmp(Subsystems(:,1),Statistics{i,1})),3};
                findSub = find(strcmp(subsystems(:,1),sub));
                 if abs(Statistics{i,8})>abs(Statistics{i,10})
                    diffReactions{1}{findSub,2} = diffReactions{1}{findSub,2}+1;
                elseif abs(Statistics{i,10})>abs(Statistics{i,8})
                    diffReactions{2}{findSub,2} = diffReactions{2}{findSub,2}+1;
                end
            end
        end
        
        % fill in tables
        table_Sig{1,colCnt} = 'Higher in controls vs. mut-';
        table_Sig(2:end,colCnt) = diffReactions{1}(:,2);
        table_Sig{1,colCnt+1} = 'Higher in mut- vs. controls';
        table_Sig(2:end,colCnt+1) = diffReactions{2}(:,2);
        colCnt = colCnt+2;
    end
    
    % Controls vs. other MMA
    col = find(strcmp(samples(1,:),'Group'));
    dataAnalysis=data;
    % remove all other cases
    otherSamples = samples(find(strcmp(samples(:,col),'mut-')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    otherSamples = samples(find(strcmp(samples(:,col),'mut0')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    % remove zero fluxes
    delArray=[];
    cnt=1;
    for j=2:size(dataAnalysis,1)
        if abs(sum(cell2mat(dataAnalysis(j,2:end))))<0.00001
            delArray(cnt)=j;
            cnt=cnt+1;
        end
    end
    dataAnalysis(delArray,:)=[];
    [Statistics,significantFeatures] = performStatisticalAnalysis(dataAnalysis,samples,'stratification',samples{1,col});
    sigRxns = union(sigRxns,significantFeatures(2:end,1));
    
    % add reaction/ metabolites names
    if strcmp(datasets{d},'Fluxes')
        for i=2:size(Statistics,1)
            findRxn = find(strcmp(Subsystems(:,1),Statistics{i,1}));
            Statistics{i,2} = Subsystems{findRxn,2};
        end
    elseif strcmp(datasets{d},'Metabolites')
        for i=2:size(Statistics,1)
            met = model.metNames{find(strcmp(model.mets,dataAnalysis{i,1}))};
            Statistics{i,2} = met;
        end
    end
    writetable(cell2table(Statistics),['Statistics_' datasets{d} filesep 'Statistics_Control_MMA_other.csv'],'WriteVariableNames',false)
    
    if strcmp(datasets{d},'Fluxes')
        % get enriched reactions-significant
        diffReactions={};
        diffReactions{1} = subsystems;
        diffReactions{1}(:,2) = {0};
        diffReactions{2} = subsystems;
        diffReactions{2}(:,2) = {0};
        diffReactions{3} = subsystems;
        diffReactions{3}(:,2) = {0};
        diffReactions{4} = subsystems;
        diffReactions{4}(:,2) = {0};
        
        % get enriched reactions-significant and initially significant
        for i=2:size(Statistics,1)
            if Statistics{i,statCol} < 0.05
                sub = Subsystems{find(strcmp(Subsystems(:,1),Statistics{i,1})),3};
                findSub = find(strcmp(subsystems(:,1),sub));
                 if abs(Statistics{i,8})>abs(Statistics{i,10})
                    diffReactions{1}{findSub,2} = diffReactions{1}{findSub,2}+1;
                elseif abs(Statistics{i,10})>abs(Statistics{i,8})
                    diffReactions{2}{findSub,2} = diffReactions{2}{findSub,2}+1;
                end
            end
        end
        
        % fill in tables
        table_Sig{1,colCnt} = 'Higher in controls vs. other MMA';
        table_Sig(2:end,colCnt) = diffReactions{1}(:,2);
        table_Sig{1,colCnt+1} = 'Higher in other MMA vs. controls';
        table_Sig(2:end,colCnt+1) = diffReactions{2}(:,2);
    end
    colCnt = colCnt+2;
    
    % mut0 vs. mut-
    col = find(strcmp(samples(1,:),'mut_category'));
    dataAnalysis=data;
    % remove all other cases
    otherSamples = samples(find(strcmp(samples(:,col),'NA')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    otherSamples = samples(find(strcmp(samples(:,col),'Control')),1);
    [C,I]=intersect(dataAnalysis(1,:),otherSamples);
    dataAnalysis(:,I)=[];
    % remove zero fluxes
    delArray=[];
    cnt=1;
    for j=2:size(dataAnalysis,1)
        if abs(sum(cell2mat(dataAnalysis(j,2:end))))<0.00001
            delArray(cnt)=j;
            cnt=cnt+1;
        end
    end
    dataAnalysis(delArray,:)=[];
    [Statistics,significantFeatures] = performStatisticalAnalysis(dataAnalysis,samples,'stratification',samples{1,col});
    sigRxns = union(sigRxns,significantFeatures(2:end,1));
    
    % add reaction/ metabolites names
    if strcmp(datasets{d},'Fluxes')
        for i=2:size(Statistics,1)
            findRxn = find(strcmp(Subsystems(:,1),Statistics{i,1}));
            Statistics{i,2} = Subsystems{findRxn,2};
        end
    elseif strcmp(datasets{d},'Metabolites')
        for i=2:size(Statistics,1)
            met = model.metNames{find(strcmp(model.mets,dataAnalysis{i,1}))};
            Statistics{i,2} = met;
        end
    end
    writetable(cell2table(Statistics),['Statistics_' datasets{d} filesep 'Statistics_' samples{1,col} '.csv'],'WriteVariableNames',false)
    
    if strcmp(datasets{d},'Fluxes')
        % get enriched reactions-significant
        diffReactions={};
        diffReactions{1} = subsystems;
        diffReactions{1}(:,2) = {0};
        diffReactions{2} = subsystems;
        diffReactions{2}(:,2) = {0};
        diffReactions{3} = subsystems;
        diffReactions{3}(:,2) = {0};
        diffReactions{4} = subsystems;
        diffReactions{4}(:,2) = {0};
        
        % get enriched reactions-significant and initially significant
        for i=2:size(Statistics,1)
            if Statistics{i,statCol} < 0.05
                sub = Subsystems{find(strcmp(Subsystems(:,1),Statistics{i,1})),3};
                findSub = find(strcmp(subsystems(:,1),sub));
                 if abs(Statistics{i,8})>abs(Statistics{i,10})
                    diffReactions{1}{findSub,2} = diffReactions{1}{findSub,2}+1;
                elseif abs(Statistics{i,10})>abs(Statistics{i,8})
                    diffReactions{2}{findSub,2} = diffReactions{2}{findSub,2}+1;
                end
            end
        end
        
        % fill in tables
        table_Sig{1,colCnt} = 'Higher in mut- vs. mut0';
        table_Sig(2:end,colCnt) = diffReactions{1}(:,2);
        table_Sig{1,colCnt+1} = 'Higher in mut0 vs. mut-';
        table_Sig(2:end,colCnt+1) = diffReactions{2}(:,2);
        
        % remove empty rows and columns
        delArraySig = [];
        for i=2:size(table_Sig,1)
            if sum(sum(cell2mat(table_Sig(i,2:end))))==0
                delArraySig(length(delArraySig)+1) = i;
            end
        end
        table_Sig(delArraySig,:) = [];
        delArraySig = [];
        for i=2:size(table_Sig,2)
            if sum(sum(cell2mat(table_Sig(2:end,i))))==0
                delArraySig(length(delArraySig)+1) = i;
            end
        end
        table_Sig(:,delArraySig) = [];
        writetable(cell2table(table_Sig),'Table_2.csv','WriteVariableNames',false)
    end
    
    %% Create boplots of fluxes
    % all reactions that were signficant in at least one comparison.
    if strcmp(datasets{d},'Fluxes')
        dataAnalysis=data;
        [C,I] = setdiff(dataAnalysis(:,1),sigRxns,'stable');
        dataAnalysis(I(2:end),:)=[];
        
        mkdir([pwd filesep 'Boxplots_' datasets{d}])
        col = find(strcmp(samples(1,:),'Group'));
        cats = {'Control','MMA_Other','mut-','mut0'};
        for i=2:size(dataAnalysis,1)
            table = {};
            cnt=1;
            for j=2:size(dataAnalysis,2)
                findSamp = find(strcmp(samples(:,1),dataAnalysis{1,j}));
                table{cnt,1}=samples{findSamp,col};
                table{cnt,2}=dataAnalysis{i,j};
                cnt=cnt+1;
            end
            
            table=cell2table(table,'VariableNames',{'Group','Value'});
            table.Group=categorical(table.Group,cats);
            
            f=figure;
            boxchart(table.Group,table.Value)
            hold on
            swarmchart(table.Group,table.Value)
            set(gca,'TickLabelInterpreter','none')
            ylabel('mmol * g dry weight-1 * hr-1')
            % get plot title
            rxn = strrep(dataAnalysis{i,1},'_min','');
            rxn = strrep(rxn,'_max','');
            title(rxn)
            set(gca,'FontSize',18)
            ax=gca;
            ax.Title.Interpreter = 'none';
            print([pwd filesep 'Boxplots_' datasets{d} filesep dataAnalysis{i,1}],'-dpng','-r300')
            close all
        end
    end
    
    %% Find fluxes that differered between absence and presence of symptoms
    
    if strcmp(datasets{d},'Metabolites')
        mkdir([pwd filesep 'Symptoms_' datasets{d}])
    end
    
    % compare symptoms separately for MMUT/MMA_Other
    group = {'MMUT','MMA_Other'};
    col = find(strcmp(samples(1,:),'MMA_type'));
    
    for g=1:length(group)
        mkdir(['Statistics_' datasets{d} filesep group{g}])
        
        % create overview
        dataAnalysis=data;
        table = dataAnalysis(:,1);
        cnt=2;
        
        subs={};
        diffReactions={};
        statCol = 4;
        for j=8:size(samples,2)
            table{1,cnt} = ['Lower in ' samples{1,j}];
            table{1,cnt+1} = ['Higher in ' samples{1,j}];
            table(2:end,cnt:cnt+1) = {0};
            
            dataAnalysis=data;
            excludeSamples = samples(find(~strcmp(samples(:,col),group{g})),1);
            [C,I]=intersect(dataAnalysis(1,:),excludeSamples(2:end));
            dataAnalysis(:,I)=[];
            % remove fluxes that are zero or all the same
            delArray=[];
            for i=2:size(dataAnalysis,1)
                if abs(sum(cell2mat(dataAnalysis(i,2:end))))<0.00001 || length(find(abs(cell2mat(data(i,2:end)) - mean(cell2mat(data(i,2:end))))<0.00001))==length(data(i,2:end))
                    delArray(length(delArray)+1)=i;
                end
            end
            dataAnalysis(delArray,:)=[];

            % only analyze symptoms with a certain number of cases
            aSamples = length(find(strcmp(samples(:,j),'Absent/ NR')));
            pSamples = length(find(strcmp(samples(:,j),'Present')));
            if aSamples>8 && pSamples>8
                [Statistics,significantFeatures] = performStatisticalAnalysis(dataAnalysis,samples,'stratification',samples{1,j});
                writetable(cell2table(Statistics),['Statistics_' datasets{d} filesep group{g} filesep 'Statistics_' samples{1,j} '.csv'],'WriteVariableNames',false)
                
                if strcmp(datasets{d},'Metabolites')
                    if size(significantFeatures,1)>1
                        % create boxplots for symptoms
                        for k=2:size(significantFeatures,1)
                            plotTable = {};
                            pcnt=1;
                            for l=2:size(significantFeatures,2)
                                findSamp = find(strcmp(samples(:,1),significantFeatures{1,l}));
                                plotTable{pcnt,1}=samples{findSamp,j};
                                plotTable{pcnt,2}=significantFeatures{k,l};
                                pcnt=pcnt+1;
                            end
                            plotTable(:,1) = strrep(plotTable(:,1),'Present',samples{1,j});
                            cats = unique(plotTable(:,1));
                            
                            plotTable=cell2table(plotTable,'VariableNames',{'Group','Value'});
                            plotTable.Group=categorical(plotTable.Group,cats);
                            
                            figure
                            boxchart(plotTable.Group,plotTable.Value)
                            hold on
                            swarmchart(plotTable.Group,plotTable.Value)
                            set(gca,'TickLabelInterpreter','none')
                            ylabel('mmol * g dry weight-1 * hr-1')
                            % get plot title
                            met = model.metNames{find(strcmp(model.mets,significantFeatures{k,1}))};
                            set(gca,'FontSize',14)
                            if strcmp(group{g},'MMUT')
                                h=title({'MMUT patients'},{met});
                            elseif strcmp(group{g},'MMA_Other')
                                h=title({'Other MMA patients'},{met});
                            end
                            set(h,'interpreter','none')
                            print([pwd filesep 'Symptoms_' datasets{d} filesep group{g} '_' samples{1,j} '_' significantFeatures{k,1}],'-dpng','-r300')
                            close all
                        end
                    end
                end
                
                % get enriched reactions
                for i=2:size(Statistics,1)
                    if Statistics{i,statCol} < 0.05
                        table{i,1} = Statistics{i,1};
                        if abs(Statistics{i,8})>abs(Statistics{i,10})
                            table{i,cnt} = 1;
                        elseif abs(Statistics{i,10})>abs(Statistics{i,8})
                            table{i,cnt+1} = 1;
                        end
                    end
                end
            end
            cnt = cnt+2;
        end
        
        % export up/downregulated reactions
        % first remove empty rows/columns
        delArray = [];
        for j=2:size(table,1)
            if sum(cell2mat(table(j,2:end)))==0
                delArray(length(delArray)+1,1) = j;
            end
        end
        table(delArray,:) = [];
        delArray = [];
        for j=2:size(table,2)
            if sum(cell2mat(table(2:end,j)))==0
                delArray(length(delArray)+1,1) = j;
            end
        end
        table(:,delArray) = [];
        % only keep one entry per reaction
        table(:,1) = strrep(table(:,1),'_min','');
        table(:,1) = strrep(table(:,1),'_max','');
        Subsystems(:,1) = strrep(Subsystems(:,1),'_min','');
        Subsystems(:,1) = strrep(Subsystems(:,1),'_max','');
        [uniqueA,i,j] = unique(table(:,1));
        n  = accumarray(j(:),1);
        Dupes=uniqueA(find(n>1));
        delArray=[];
        for i=1:length(Dupes)
            indexToDupes = find(strcmp(table(:,1),Dupes{i}));
                delArray(i,1)=indexToDupes(2);
        end
        table(delArray,:)=[];
        writetable(cell2table(table),['Statistics_' datasets{d} filesep group{g} '_Symptoms.csv'],'WriteVariableNames',false)
        
        if strcmp(datasets{d},'Fluxes')
            exportSubs = {'Reaction ID','Subsystem'};
            for i=2:size(table,1)
                exportSubs{i,1} = table{i,1};
                exportSubs{i,2} = Subsystems{find(strcmp(Subsystems(:,1),table{i,1})),3};
            end
            
            if length(unique(exportSubs(2:end,2)))>15
                [C,subs] = groupcounts(exportSubs(2:end,2));
                [A,I] = sort(C,'descend');
                subs = subs(I);
                for j = 16:length(subs)
                    exportSubs(2:end,2) = strrep(exportSubs(2:end,2),subs{j},'Others');
                end
            end
            writetable(cell2table(exportSubs),['Statistics_' datasets{d} filesep group{g} '_Symptoms_Subsystems.csv'],'WriteVariableNames',false)
        end
    end
end
