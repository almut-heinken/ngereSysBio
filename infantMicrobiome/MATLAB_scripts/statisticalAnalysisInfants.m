
% perform statistical analyses for gut microbiomes from infants delivered
% vaginally and by Cesarian section separate by time point, as well as
% antibiotics vs. no antibiotics

rootDir = pwd;

metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);
timePoints = {'5 days','1 month','6 months','1 year'};

datasets = {
    ['MicrobiomeModels' filesep 'Infant_diet.csv'],'Net_secretion'
    ['MicrobiomeModels' filesep 'ReactionAbundance.csv'],'Reaction_abundance'
    ['MicrobiomeModels' filesep 'ReactionPresence.csv'],'Reaction_presence'
    ['MicrobiomeModels' filesep 'SubsystemAbundance.csv'],'Subsystem_abundance'
    ['Contributions' filesep 'Microbe_Secretion.csv'],'Microbe_secretion'
    };

cd('inputFiles')
database = loadVMHDatabase;
cd ..

mkdir([rootDir filesep 'Statistical_analysis_COSMIC'])
mkdir([rootDir filesep 'Metabolite_Plots_COSMIC'])

for d=1:size(datasets,1)
    data = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep datasets{d,1}]);
    data(1,:) = strrep(data(1,:),'microbiota_model_diet_','');
    for i=1:length(timePoints)
        % analyze VD and CSD samples by time point
        statistics = {'Feature','Mean VD','SD VD','Mean CSD','SD CSD','p-value','after FDR'};

        for j=2:size(data,1)
            findSamp = find(strcmp(metadata(:,4),timePoints{i}));
            samples = metadata(findSamp,1);
            strat = metadata(findSamp,2);
            findVD = samples(find(strcmp(strat,'VD')),1);
            findCSD = samples(find(strcmp(strat,'CSD')),1);
            [~,I] = intersect(data(1,:),findVD);
            dataVD = cell2mat(data(j,I));
            [~,I] = intersect(data(1,:),findCSD);
            dataCSD = cell2mat(data(j,I));
            [p,h,stats] = ranksum(dataVD,dataCSD);
            if d==1
                met = strrep(data{j,1},'EX_','');
                met = strrep(met,'[fe]','');
                met = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};

                statistics{j,1} = met;
            else
                statistics{j,1} = data{j,1};
            end
            statistics{j,2} = mean(dataVD);
            statistics{j,3} = std(dataVD);
            statistics{j,4} = mean(dataCSD);
            statistics{j,5} = std(dataCSD);
            statistics{j,6} = p;
        end
        pAverages=cell2mat(statistics(2:end,6));
        fdr = mafdr(pAverages,'BHFDR', true);
        statistics(2:end,7)=num2cell(fdr);
        writetable(cell2table(statistics),[rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{d,2} '_' strrep(timePoints{i},' ','_') '_birthMode.csv'],'WriteVariableNames',false)
        
        % analyze antiboitcs vs. no antibiotics by time point
        statistics = {'Feature','Mean None','SD None','Mean Antibiotics','SD Antibiotics','p-value','after FDR'};

        for j=2:size(data,1)
            findSamp = find(strcmp(metadata(:,4),timePoints{i}));
            samples = metadata(findSamp,1);
            strat = metadata(findSamp,5);
            findNA = samples(find(strcmp(strat,'None')),1);
            findA = samples(find(strcmp(strat,'Antibiotics')),1);
            [~,I] = intersect(data(1,:),findNA);
            dataNA = cell2mat(data(j,I));
            [~,I] = intersect(data(1,:),findA);
            dataA = cell2mat(data(j,I));
            [p,h,stats] = ranksum(dataNA,dataA);
            if d==1
                met = strrep(data{j,1},'EX_','');
                met = strrep(met,'[fe]','');
                met = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};

                statistics{j,1} = met;
            else
                statistics{j,1} = data{j,1};
            end
            statistics{j,2} = mean(dataNA);
            statistics{j,3} = std(dataNA);
            statistics{j,4} = mean(dataA);
            statistics{j,5} = std(dataA);
            statistics{j,6} = p;
        end
        pAverages=cell2mat(statistics(2:end,6));
        fdr = mafdr(pAverages,'BHFDR', true);
        statistics(2:end,7)=num2cell(fdr);
        writetable(cell2table(statistics),[rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{d,2} '_' strrep(timePoints{i},' ','_') '_antibiotics.csv'],'WriteVariableNames',false)
    end

    % analyze differences between time points
    data = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep datasets{d,1}]);

    % create the table
    statistics = {'Feature','Mean 5 days','SD 5 days','Mean 1 month','SD 1 month','Mean 6 months','SD 6 months','Mean 1 year','SD 1 year','p-value','after FDR'};

    for i=2:size(data,1)
        find5d = metadata(find(strcmp(metadata(:,4),'5 days')),1);
        find1m = metadata(find(strcmp(metadata(:,4),'1 month')),1);
        find6m = metadata(find(strcmp(metadata(:,4),'6 months')),1);
        find1y = metadata(find(strcmp(metadata(:,4),'1 year')),1);

        group={};

        [~,I] = intersect(data(1,:),find5d);
        data5d = cell2mat(data(i,I));
        for j=1:length(I)
            group{length(group)+1,1}='5 days';
        end
        [~,I] = intersect(data(1,:),find1m);
        data1m = cell2mat(data(i,I));
        for j=1:length(I)
            group{length(group)+1,1}='1 month';
        end
        [~,I] = intersect(data(1,:),find6m);
        data6m = cell2mat(data(i,I));
        for j=1:length(I)
            group{length(group)+1,1}='6 months';
        end
        [~,I] = intersect(data(1,:),find1y);
        data1y = cell2mat(data(i,I));

        for j=1:length(I)
            group{length(group)+1,1}='1 year';
        end

        dataAll=[data5d,data1m,data6m,data1y];
        [p,ANOVATAB] = kruskalwallis(dataAll,group,'off');

        if d==1
            met = strrep(data{i,1},'EX_','');
            met = strrep(met,'[fe]','');
            met = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};
            statistics{i,1} = met;
        else
            statistics{i,1} = data{i,1};
        end
        statistics{i,2} = mean(data5d);
        statistics{i,3} = std(data5d);
        statistics{i,4} = mean(data1m);
        statistics{i,5} = std(data1m);
        statistics{i,6} = mean(data6m);
        statistics{i,7} = std(data6m);
        statistics{i,8} = mean(data1y);
        statistics{i,9} = std(data1y);
        statistics{i,10} = p;
    end
    pAverages=cell2mat(statistics(2:end,10));
    fdr = mafdr(pAverages,'BHFDR', true);
    statistics(2:end,11)=num2cell(fdr);
    writetable(cell2table(statistics),[rootDir filesep 'Statistical_analysis_COSMIC' filesep datasets{d,2} '_byTimePoint.csv'],'WriteVariableNames',false)
    
    % plot the results
    if d==1
        for l=2:size(data,1)
            if sum(cell2mat(data(l,2:end)))>0.0001 % && statistics{l,6} < 0.05
                table={};
                cnt=1;
                for j=1:length(timePoints)
                    findSamp = find(strcmp(metadata(:,4),timePoints{j}));
                    samples = metadata(findSamp,1);
                    strat = metadata(findSamp,2);
                    findVD = samples(find(strcmp(strat,'VD')),1);
                    findCSD = samples(find(strcmp(strat,'CSD')),1);
                    [~,I] = intersect(data(1,:),findVD);
                    dataVD = cell2mat(data(l,I));
                    [~,I] = intersect(data(1,:),findCSD);
                    dataCSD = cell2mat(data(l,I));
                    
                    for k=1:length(dataVD)
                        table{cnt,1}='VD';
                        table{cnt,2}=timePoints{j};
                        table{cnt,3}=dataVD(k);
                        cnt=cnt+1;
                    end
                    for k=1:length(dataCSD)
                        table{cnt,1}='CSD';
                        table{cnt,2}=timePoints{j};
                        table{cnt,3}=dataCSD(k);
                        cnt=cnt+1;
                    end
                end
                
                table=cell2table(table,'VariableNames',{'TimePoint','Group','Value'});
                table.Group=categorical(table.Group,timePoints);
                
                figure
                boxchart(table.Group,table.Value,'GroupByColor',table.TimePoint)
                legend('Location','best')
                if d==1 || d>5
                    ylabel('mmol/person/day')
                elseif d==3
                    ylabel('Relative presence')
                else
                    ylabel('Relative abundance')
                end
                if d==1
                    met = strrep(data{l,1},'EX_','');
                    met = strrep(met,'[fe]','');
                    met = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};
                    h=title(met);
                else
                    h=title(data{l,1});
                end
                set(h,'interpreter','none')
                xticklabels({'5 days','1 month','6 months','1 year'})
                set(gca,'FontSize',18)
                filename = strrep(data{l,1},' ','_');
                filename = strrep(filename,'(','_');
                filename = strrep(filename,')','_');
                filename = strrep(filename,',','_');
                filename = strrep(filename,'/','_');
                filename = strrep(filename,'[','_');
                filename = strrep(filename,']','');
                print([rootDir filesep 'Metabolite_Plots_COSMIC' filesep datasets{d,2} '_' filename],'-dpng','-r300')
                close all
            end
        end
    end
end
