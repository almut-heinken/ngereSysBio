
% create Table S7a-b: reaction abundance and presence enriched and 
% depleted in VD vs. CSD

metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);
timePoints = {'5 days','1 month','6 months','1 year'};

cd('inputFiles')
database=loadVMHDatabase;
cd(rootDir)

datasets = {
    'Reaction_abundance','Table_S7a'
    'Reaction_presence','Table_S7b'
    };

for i=1:length(datasets)
    statTable = {'Subsystem'};
    subs = unique(database.reactions(:,11));
    statTable(2:length(subs)+1,1) = subs;
    cnt = 2;

    for j=1:length(timePoints)
        stats=readInputTableForPipeline([rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{i,1} '_' strrep(timePoints{j},' ','_') '_birthMode.csv']);
        stats(find(strcmp(stats(:,1),'biomassPan')),:)=[];
        % get the reactions that are significantly different after
        % correction for FDR, and only before correction for FDR
        % 1st: higher in VD
        % 2nd: higher in CSD
        rxnsListAfterFDR=cell(1,2);
        rxnsListBeforeFDR=cell(1,2);
        for k=2:size(stats,1)
            if stats{k,7}<0.05
                % significant after FDR
                if stats{k,2}>stats{k,4}
                    rxnsListAfterFDR{1}{end+1,1}=stats{k,1};
                elseif stats{k,4}>stats{k,2}
                    rxnsListAfterFDR{2}{end+1,1}=stats{k,1};
                end
            end
            if stats{k,6}<0.05
                % significant before FDR
                if stats{k,2}>stats{k,4}
                    rxnsListBeforeFDR{1}{end+1,1}=stats{k,1};
                elseif stats{k,4}>stats{k,2}
                    rxnsListBeforeFDR{2}{end+1,1}=stats{k,1};
                end
            end
        end
        % get subsystems
        for k=1:size(rxnsListAfterFDR{1},1)
            rxnsListAfterFDR{1}{k,2}=database.reactions{find(strcmp(database.reactions(:,1),rxnsListAfterFDR{1}{k,1})),11};
        end
        for k=1:size(rxnsListAfterFDR{2},1)
            rxnsListAfterFDR{2}{k,2}=database.reactions{find(strcmp(database.reactions(:,1),rxnsListAfterFDR{2}{k,1})),11};
        end
        for k=1:size(rxnsListBeforeFDR{1},1)
            rxnsListBeforeFDR{1}{k,2}=database.reactions{find(strcmp(database.reactions(:,1),rxnsListBeforeFDR{1}{k,1})),11};
        end
        for k=1:size(rxnsListBeforeFDR{2},1)
            rxnsListBeforeFDR{2}{k,2}=database.reactions{find(strcmp(database.reactions(:,1),rxnsListBeforeFDR{2}{k,1})),11};
        end
        % fill in the table
         statTable{1,cnt} = [timePoints{j} ', higher in VD, before correction for FDR'];
        statTable(2:end,cnt) = {'0'};
        for k=1:size(rxnsListBeforeFDR{1},1)
            findSub = find(strcmp(statTable(:,1),rxnsListBeforeFDR{1}{k,2}));
            statTable{findSub,cnt} = num2str(str2double(statTable{findSub,cnt}) + 1);
        end
        cnt = cnt+1;
        statTable{1,cnt} = [timePoints{j} ', higher in CSD, before correction for FDR'];
        statTable(2:end,cnt) = {'0'};
        for k=1:size(rxnsListBeforeFDR{2},1)
            findSub = find(strcmp(statTable(:,1),rxnsListBeforeFDR{2}{k,2}));
            statTable{findSub,cnt} = num2str(str2double(statTable{findSub,cnt}) + 1);
        end
        cnt = cnt+1;
        statTable{1,cnt} = [timePoints{j} ', higher in VD, after correction for FDR'];
        statTable(2:end,cnt) = {'0'};
        for k=1:size(rxnsListAfterFDR{1},1)
            findSub = find(strcmp(statTable(:,1),rxnsListAfterFDR{1}{k,2}));
            statTable{findSub,cnt} = num2str(str2double(statTable{findSub,cnt}) + 1);
        end
        cnt = cnt+1;
        statTable{1,cnt} = [timePoints{j} ', higher in CSD, after correction for FDR'];
        statTable(2:end,cnt) = {'0'};
        for k=1:size(rxnsListAfterFDR{2},1)
            findSub = find(strcmp(statTable(:,1),rxnsListAfterFDR{2}{k,2}));
            statTable{findSub,cnt} = num2str(str2double(statTable{findSub,cnt}) + 1);
        end
        cnt = cnt+1;
    end
    % remove zeros
    delArray = [];
    for j=2:size(statTable,1)
        if sum(str2double(statTable(j,2:end)))==0
            delArray(length(delArray)+1,1) = j;
        end
    end
    statTable(delArray,:) = [];
    writetable(cell2table(statTable),[datasets{i,2} '.csv'],'writeVariableNames',false)
end

