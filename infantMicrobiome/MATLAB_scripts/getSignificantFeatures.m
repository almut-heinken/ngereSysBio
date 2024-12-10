
% get an overview of all significant features between VD and CSD for Table 1

datasets = {'Net_secretion','Reaction_abundance','Reaction_presence','Subsystem_abundance','Microbe_secretion'};

% create a table with statistics overview
table={};
analyses = {'byTimePoint','5_days_birthMode','1_month_birthMode','6_months_birthMode','1_year_birthMode','5_days_antibiotics','1_month_antibiotics','6_months_antibiotics','1_year_antibiotics'};

for i=1:length(datasets)
    table{1,i+1}=datasets{i};
    for j=1:length(analyses)
        table{j+1,1}=analyses{j};
        stats=readInputTableForPipeline(['Statistical_analysis_COSMIC' filesep datasets{i} '_' analyses{j} '.csv']);
        % find number of significant features after FDR correction
        sigFeats=length(find(cell2mat(stats(2:end,end))<0.05));
        % find number of initially significant features
        inSigFeats=length(find(cell2mat(stats(2:end,end-1))<0.05));
        table{j+1,i+1}=[num2str(sigFeats) '(' num2str(inSigFeats-sigFeats) ')' '/' num2str(size(stats,1)-1)];
        cnt=cnt+1;
    end
end
writetable(cell2table(table),['Table_1' '.csv'],'writeVariableNames',false)
