
% export predicted features and fluxes for infant gut microbiomes

cd('inputFiles')
database=loadVMHDatabase;
cd(rootDir)

%% define datasets
datasets = {
    ['MicrobiomeModels' filesep 'ReactionAbundance.csv'],'ReactionAbundance.csv','Reaction_abundance_'
    ['MicrobiomeModels' filesep 'ReactionPresence.csv'],'ReactionPresence.csv','Reaction_presence_'
    ['MicrobiomeModels' filesep 'Infant_diet.csv'],'NetSecretion.csv','Net_secretion_'
    };

%% export statistically significant features by time point
for i=1:size(datasets,1)
    data = readInputTableForPipeline([pwd filesep 'Modeling_COSMIC' filesep datasets{i,1}]);

    % remove small contributions
    delArray = [];
    for j=1:size(data,1)
        delArray(j,1) = sum(cell2mat(data(j,2:end)));
    end
    data(find(delArray<0.00001),:) = [];

    % limit to statistically significant
    stats = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{i,3} 'byTimePoint.csv']);
    stats(1,:) = [];
    sigFeats = stats(find(cell2mat(stats(:,11))<0.05),1);

    if i==3
        % translate metabolites
        for j=2:size(data,1)
            data{j,1} = strrep(data{j,1},'EX_','');
            data{j,1} = strrep(data{j,1},'[fe]','');
            data{j,1} = database.metabolites{find(strcmp(database.metabolites(:,1),data{j,1})),2};
        end
    end
    
    if i==3
        % first export all data
        writetable(cell2table(data),['Exported_results' filesep 'ForSupplement_Infants_' datasets{i,2}],'writeVariableNames',false)
    end
    
    % remove nonsignificant features
    [~,I] = setdiff(data(:,1),sigFeats,'stable');
    data(I(2:end),:) = [];

    if i==3
        % keep 50 most abundant
        delArray = [];
        for j=1:size(data,1)
            delArray(j,1) = sum(cell2mat(data(j,2:end)));
        end
        [I,J] = sort(delArray,'descend');
        data(J(51:end),:) = [];
    end

    % add annotation for reactions
    if i==1
        annotation = {'Reaction','Subsystem_detailed','Subsystem'};
        for j=2:size(data,1)
            annotation{j,1} = data{j,1};
            if ~strncmp(data{j,1},'bio',3)
                findRxn = find(strcmp(database.reactions(:,1),data{j,1}));
                annotation{j,2} = database.reactions{findRxn,11};
                annotation{j,3} = database.reactions{findRxn,12};
            end
        end
        writetable(cell2table(annotation),[rootDir filesep 'Exported_results' filesep 'ByTimePoint_ReactionAbundanceSubsystems.csv'],'writeVariableNames',false)
    end
    if i==2
        annotation = {'Reaction','Subsystem_detailed','Subsystem'};
        for j=2:size(data,1)
            annotation{j,1} = data{j,1};
            if ~strncmp(data{j,1},'bio',3)
                findRxn = find(strcmp(database.reactions(:,1),data{j,1}));
                annotation{j,2} = database.reactions{findRxn,11};
                annotation{j,3} = database.reactions{findRxn,12};
            end
        end
        writetable(cell2table(annotation),[rootDir filesep 'Exported_results' filesep 'ByTimePoint_ReactionPresenceSubsystems.csv'],'writeVariableNames',false)
    end
    writetable(cell2table(data),[rootDir filesep 'Exported_results' filesep 'ByTimePoint_' datasets{i,2}],'writeVariableNames',false)
end

%% export statistically significant features by group
for i=1:size(datasets,1)
    data = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep datasets{i,1}]);

    % remove small contributions
    delArray = [];
    for j=1:size(data,1)
        delArray(j,1) = sum(cell2mat(data(j,2:end)));
    end
    data(find(delArray<0.00001),:) = [];

    % limit to statistically significant in at least one comparison
    timePoints = {'5_days','1_month','6_months','1_year'};
    sigFeats = {};
    for j=1:length(timePoints)
        stats = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{i,3} timePoints{j} '_birthMode.csv']);
        stats(1,:) = [];
        sig = find(cell2mat(stats(:,7))<0.05);
        % sig = find(cell2mat(stats(:,6))<0.05);
        sigFeats = union(sigFeats,stats(sig,1));
    end

    if i==4
        % translate metabolites
        for j=2:size(data,1)
            data{j,1} = strrep(data{j,1},'EX_','');
            data{j,1} = strrep(data{j,1},'[fe]','');
            data{j,1} = database.metabolites{find(strcmp(database.metabolites(:,1),data{j,1})),2};
        end
    end

    % remove nonsignificant features
    [~,I] = setdiff(data(:,1),sigFeats,'stable');
    data(I(2:end),:) = [];

    if i==3 || i==4
        % keep 50 most abundant
        delArray = [];
        for j=1:size(data,1)
            delArray(j,1) = sum(cell2mat(data(j,2:end)));
        end
        [I,J] = sort(delArray,'descend');
        data(J(51:end),:) = [];
    end

    if i==4
        % scale by row
        for j=2:size(data,1)
            maxval = max(cell2mat(data(j,2:end)));
            for k=2:size(data,2)
                data{j,k} = data{j,k}/maxval;
            end
        end
    end

    % add annotation for reactions
    if i==1
        annotation = {'Reaction','Subsystem_detailed','Subsystem'};
        for j=2:size(data,1)
            annotation{j,1} = data{j,1};
            if ~strncmp(data{j,1},'bio',3)
                findRxn = find(strcmp(database.reactions(:,1),data{j,1}));
                annotation{j,2} = database.reactions{findRxn,11};
                annotation{j,3} = database.reactions{findRxn,12};
            end
        end
        writetable(cell2table(annotation),[rootDir filesep 'Exported_results' filesep 'ByGroup_ReactionAbundanceSubsystems.csv'],'writeVariableNames',false)
    end
    if i==2
        annotation = {'Reaction','Subsystem_detailed','Subsystem'};
        for j=2:size(data,1)
            annotation{j,1} = data{j,1};
            if ~strncmp(data{j,1},'bio',3)
                findRxn = find(strcmp(database.reactions(:,1),data{j,1}));
                annotation{j,2} = database.reactions{findRxn,11};
                annotation{j,3} = database.reactions{findRxn,12};
            end
        end
        writetable(cell2table(annotation),[rootDir filesep 'Exported_results' filesep 'ByGroup_ReactionPresenceSubsystems.csv'],'writeVariableNames',false)
    end
    writetable(cell2table(data),[rootDir filesep 'Exported_results' filesep 'ByGroup_' datasets{i,2}],'writeVariableNames',false)
end

%% export infants with adults combined
datasets = {
    ['MicrobiomeModels' filesep 'ReactionAbundance.csv'],'ReactionAbundance.csv','Reaction_abundance_'
    ['MicrobiomeModels' filesep 'ReactionPresence.csv'],'ReactionPresence.csv','Reaction_presence_'
    ['MicrobiomeModels' filesep 'SubsystemAbundance.csv'],'SubsystemAbundance.csv','Subsystem_abundance_'
    ['MicrobiomeModels' filesep 'AE_diet.csv'],'NetSecretion.csv','Net_secretion_'
    ['Contributions' filesep 'Microbe_Secretion.csv'],'Microbe_Secretion.csv','Microbe_secretion_'
    };

for i=1:size(datasets,1)
    dataInfants = readInputTableForPipeline([pwd filesep 'Modeling_COSMIC' filesep datasets{i,1}]);
    dataGutMaternal = readInputTableForPipeline([pwd filesep 'Modeling_MaternalGutMicrobiomes' filesep datasets{i,1}]);

    combFeats=union(dataInfants(2:end,1),dataGutMaternal(2:end,1));
    data = [{''}
        combFeats];

    for j=1:size(combFeats,1)
        dataTmp = [];
        findRxn = find(strcmp(dataInfants(:,1),combFeats{j,1}));
        if ~isempty(findRxn)
            dataTmp = cell2mat(dataInfants(findRxn,2:end));
        else
            dataTmp = zeros(1,size(dataInfants,2)-1);
        end

         % add the maternal samples
        findInM = find(strcmp(dataGutMaternal(:,1),combFeats{j,1}));
        if ~isempty(findInM)
            dataM = cell2mat(dataGutMaternal(findInM,2:end));
        else
            dataM = zeros(1,size(dataGutMaternal,2)-1);
        end
        dataTmp = [dataTmp,dataM];

        data{j+1,1} = combFeats{j};
        data(j+1,2:length(dataTmp)+1) = num2cell(dataTmp);
    end
    modNames = dataInfants(1,2:end);
    modNames = [modNames,dataGutMaternal(1,2:end)];
    data(1,2:length(dataTmp)+1) = modNames;

    % remove very small contributions
    delArray = [];
    for j=1:size(data,1)
        delArray(j,1) = sum(cell2mat(data(j,2:end)));
    end
    data(find(delArray<0.00001),:) = [];
    
    if i==4
        % translate metabolites
        for j=2:size(data,1)
            data{j,1} = strrep(data{j,1},'EX_','');
            data{j,1} = strrep(data{j,1},'[fe]','');
            data{j,1} = database.metabolites{find(strcmp(database.metabolites(:,1),data{j,1})),2};
        end
    end
    
     % export all data for supplemental tables
    writetable(cell2table(data),['Exported_results' filesep 'ForSupplement_All_' datasets{i,2}],'writeVariableNames',false)
end

%% summarize statistical results in one table (S10,S11,S13)

files = {
    'Net_secretion_5_days_birthMode.csv'
    'Net_secretion_1_month_birthMode.csv'
    'Net_secretion_6_months_birthMode.csv'
    'Net_secretion_1_year_birthMode.csv'
    'Net_secretion_5_days_antibiotics.csv'
    'Net_secretion_1_month_antibiotics.csv'
    'Net_secretion_6_months_antibiotics.csv'
    'Net_secretion_1_year_antibiotics.csv'
    'Net_secretion_byTimePoint.csv'
    };

stats = {'','VD vs. CSD five days','VD vs. CSD one month','VD vs. CSD six months','VD vs. CSD one year','No antibiotics vs. antibiotics five days','No antibiotics vs. antibiotics one month','No antibiotics vs. antibiotics six months','No antibiotics vs. antibiotics one year','By time point','VD vs. CSD five days','VD vs. CSD one month','VD vs. CSD six months','VD vs. CSD one year','No antibiotics vs. antibiotics five days','No antibiotics vs. antibiotics one month','No antibiotics vs. antibiotics six months','No antibiotics vs. antibiotics one year','By time point'};
for i=1:length(files)
    statistics = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_COSMIC' filesep files{i}]);
    if i==1
        stats(2:length(statistics),1) = statistics(2:end,1);
    end
    if i<9
        stats(2:end,i+1) = statistics(2:end,6);
        stats(2:end,i+10) = statistics(2:end,7);
    else
        stats(2:end,i+1) = statistics(2:end,10);
        stats(2:end,i+10) = statistics(2:end,11);
    end
end
delArray = [];
cnt=1;
for i=2:size(stats,1)
    if any(isnan(cell2mat(stats(i,2:end))))
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
stats(delArray,:) = [];
writetable(cell2table(stats),['Exported_results' filesep 'Table_S10.csv'],'writeVariableNames',false)

files = {
    'Subsystem_abundance_5_days.csv'
    'Subsystem_abundance_1_month.csv'
    'Subsystem_abundance_6_months.csv'
    'Subsystem_abundance_1_year.csv'
    };

stats = {'','Maternal gut vs. infant gut five days','Maternal gut vs. infant gut one month','Maternal gut vs. infant gut six months','Maternal gut vs. infant gut one year'};
for i=1:length(files)
    statistics = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_with_maternalGut_COSMIC' filesep files{i}]);
    if i==1
        stats(2:length(statistics),1) = statistics(2:end,1);
    end
    stats(2:end,i+1) = statistics(2:end,7);
end
delArray = [];
cnt=1;
for i=2:size(stats,1)
    if any(isnan(cell2mat(stats(i,2:end))))
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
stats(delArray,:) = [];
writetable(cell2table(stats),['Exported_results' filesep 'Table_S11.csv'],'writeVariableNames',false)

files = {
    'Net_secretion_AD_5_days.csv'
    'Net_secretion_AD_1_month.csv'
    'Net_secretion_AD_6_months.csv'
    'Net_secretion_AD_1_year.csv'
    };

stats = {'','Maternal gut vs. infant gut five days','Maternal gut vs. infant gut one month','Maternal gut vs. infant gut six months','Maternal gut vs. infant gut one year'};
for i=1:length(files)
    statistics = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_with_maternalGut_COSMIC' filesep files{i}]);
    if i==1
        stats(2:length(statistics),1) = statistics(2:end,1);
    end
    stats(2:end,i+1) = statistics(2:end,7);
end
delArray = [];
cnt=1;
for i=2:size(stats,1)
    if all(isnan(cell2mat(stats(i,2:end))))
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
stats(delArray,:) = [];
writetable(cell2table(stats),['Exported_results' filesep 'Table_S13.csv'],'writeVariableNames',false)
