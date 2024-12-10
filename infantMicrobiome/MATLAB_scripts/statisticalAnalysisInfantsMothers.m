
% perform statistical analyses for gut microbiomes from infants and maternal
% microbiomes

metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);
timePoints = unique(metadata(2:end,4));

cd('inputFiles')
database = loadVMHDatabase;
cd ..

%% perform for infant and maternal gut samples

datasets = {
    ['MicrobiomeModels' filesep 'AE_diet.csv'],'Net_secretion_AD'
    ['MicrobiomeModels' filesep 'ReactionAbundance.csv'],'Reaction_abundance'
    ['MicrobiomeModels' filesep 'ReactionPresence.csv'],'Reaction_presence'
    ['MicrobiomeModels' filesep 'SubsystemAbundance.csv'],'Subsystem_abundance'
    ['Contributions' filesep 'Microbe_Secretion.csv'],'Microbe_secretion'
    };

mkdir([rootDir filesep 'Statistical_analysis_with_maternalGut_COSMIC'])

% export significant features

for d=1:size(datasets)
    for i=1:length(timePoints)
        data = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep datasets{d,1}]);
        if d==1
            % remove net secretion fluxes that are zero
            delArray=[];
            cnt=1;
            for j=2:size(data,1)
                if abs(sum(cell2mat(data(j,2:end))))<0.000001 || all(cell2mat(data(j,2:end)) == data{j,2})
                    delArray(cnt)=j;
                    cnt=cnt+1;
                end
            end
            data(delArray,:) = [];
        end
        dataM = readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep datasets{d,1}]);
        if d==1
            % remove net secretion fluxes that are zero
            delArray=[];
            cnt=1;
            for j=2:size(dataM,1)
                if abs(sum(cell2mat(dataM(j,2:end))))<0.000001 || all(cell2mat(dataM(j,2:end)) == dataM{j,2})
                    delArray(cnt)=j;
                    cnt=cnt+1;
                end
            end
            dataM(delArray,:) = [];
        end
        feats=union(data(2:end,1),dataM(2:end,1));
        
        findSamp = find(strcmp(metadata(:,4),timePoints{i}));
        samples = metadata(findSamp,1);
        [~,I] = intersect(data(1,:),samples);
        
        % get table
        statistics = {'Feature','Mean Infants','SD Infants','Mean Adults','SD Adults','p-value','after FDR','Decision'};
        
        for j=1:size(feats,1)
            findRxn = find(strcmp(data(:,1),feats{j,1}));
            if ~isempty(findRxn)
                infantsTest = cell2mat(data(findRxn,I));
            else
                infantsTest = zeros(1,length(I));
            end
            % add the adult samples
            findInM = find(strcmp(dataM(:,1),feats{j,1}));
            if ~isempty(findInM)
                adultsTest = cell2mat(dataM(findInM,2:end));
            else
                adultsTest = zeros(1,size(dataM,2)-1);
            end
            [p,h,stats] = ranksum(infantsTest,adultsTest);

            if d==1
                feats{j,1} = strrep(feats{j,1},'EX_','');
                feats{j,1} = strrep(feats{j,1},'[fe]','');
                feats{j,1} = database.metabolites{find(strcmp(database.metabolites(:,1),feats{j,1})),2};
            end
            
            statistics{j+1,1} = feats{j,1};
            statistics{j+1,2} = mean(infantsTest);
            statistics{j+1,3} = std(infantsTest);
            statistics{j+1,4} = mean(adultsTest);
            statistics{j+1,5} = std(adultsTest);
            statistics{j+1,6} = p;
        end
        pAverages=cell2mat(statistics(2:end,6));
        fdr = mafdr(pAverages,'BHFDR', true);
        statistics(2:end,7)=num2cell(fdr);
        writetable(cell2table(statistics),[rootDir filesep 'Statistical_analysis_with_maternalGut_COSMIC' filesep datasets{d,2} '_' strrep(timePoints{i},' ','_') '.csv'],'WriteVariableNames',false)
    end
end
