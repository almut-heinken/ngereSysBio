
metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);
timePoints = {'5 days','1 month','6 months','1 year','Maternal gut'};

cd([rootDir filesep 'inputFiles'])
database = loadVMHDatabase;
cd ..

groups = {'5 days','1 month','6 months','1 year','dummy','Maternal gut'};

% plot metabolites resolved by time point and group
mkdir([rootDir filesep 'Metabolite_Plots_With_Mothers'])

fluxes = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep 'MicrobiomeModels' filesep 'AE_diet.csv']);
fluxesM = readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'MicrobiomeModels' filesep 'AE_diet.csv']);
mets = union(fluxes(2:end,1),fluxesM(2:end,1));

for i=2:size(mets,1)
    met = strrep(mets{i,1},'EX_','');
    met = strrep(met,'[fe]','');
    met = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};
    table={};
    cnt=1;
    for j=1:length(timePoints)
        findSamp = find(strcmp(metadata(:,4),timePoints{j}));
        samples = metadata(findSamp,1);
        strat = metadata(findSamp,2);
        %         findVD = samples(find(strcmp(strat,'VD')),1);
        %         findCSD = samples(find(strcmp(strat,'CSD')),1);
        findMet = find(strcmp(fluxes(:,1),mets{i,1}));
        if ~isempty(findMet)
            %             [~,I] = intersect(fluxes(1,:),findVD);
            %             dataVD = cell2mat(fluxes(findMet,I));
            %             [~,I] = intersect(fluxes(1,:),findCSD);
            %             dataCSD = cell2mat(fluxes(findMet,I));
            [~,I] = intersect(fluxes(1,:),samples);
            data = cell2mat(fluxes(findMet,I));
        else
            %             [~,I] = intersect(fluxes(1,:),findVD);
            %             dataVD = zeros(1,length(I));
            %             [~,I] = intersect(fluxes(1,:),findCSD);
            %             dataCSD = zeros(1,length(I));
            data = zeros(1,length(findSamp));
        end

        for k=1:length(data)
            table{cnt,1}=timePoints{j};
            table{cnt,2}=data(k);
            cnt=cnt+1;
        end
%         for k=1:length(dataVD)
%             table{cnt,1}='VD';
%             table{cnt,2}=timePoints{j};
%             table{cnt,3}=dataVD(k);
%             cnt=cnt+1;
%         end
%         for k=1:length(dataCSD)
%             table{cnt,1}='CSD';
%             table{cnt,2}=timePoints{j};
%             table{cnt,3}=dataCSD(k);
%             cnt=cnt+1;
%         end
    end
    % add the adult samples + dummy values for plot
    findMet = find(strcmp(fluxesM(:,1),mets{i,1}));
    if ~isempty(findMet)
        dataM = cell2mat(fluxesM(findMet,2:end));
    else
        dataM = zeros(1,size(fluxesM,2)-1);
    end
    for j=1:length(dataM)
        table{cnt,1}='Maternal gut';
        table{cnt,2}=dataM(j);
        cnt=cnt+1;
    end
 
    if sum(data)>0 || sum(dataM)>0
        
        %         table=cell2table(table,'VariableNames',{'Group','TimePoint','Value'});
        table=cell2table(table,'VariableNames',{'Group','Value'});
        table.Group=categorical(table.Group,groups);
        %         table.TimePoint=categorical(table.TimePoint,groups);

        figure
        boxchart(table.Value,'GroupByColor',table.Group)
%         boxchart(table.TimePoint,table.Value,'GroupByColor',table.Group)
        legend('5 days','1 month','6 months','1 year','Maternal gut','Location','best')
%         legend('Adults','CSD','VD','Location','best')
%         xticklabels({'5 days','1 month','6 months','1 year','Maternal gut'})
        ylabel('mmol/person/day')
        xticklabels(met)
        set(gca,'TickLabelInterpreter','none');
        set(gca,'FontSize',16)
        print([rootDir filesep 'Metabolite_Plots_With_Mothers' filesep 'Fluxes_' mets{i,1}],'-dpng','-r300')
        close all
    end
end

% export infant and maternal gut microbiome subsystem abundances for Figure
% S2

metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);

% export combined datasets
datasets = {
    ['MicrobiomeModels' filesep 'SubsystemAbundance.csv'],'SubsystemAbundance.csv','Subsystem_abundance_'
    };

mkdir([rootDir filesep 'R_plots'])

for i=1:size(datasets,1)
    dataInfants = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep datasets{i,1}]);
    dataMothers = readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep datasets{i,1}]);

    combFeats=union(dataInfants(2:end,1),dataMothers(2:end,1));
    data = [{''}
        combFeats];

    for j=1:size(combFeats,1)
        findRxn = find(strcmp(dataInfants(:,1),combFeats{j,1}));
        if ~isempty(findRxn)
            dataTmp = cell2mat(dataInfants(findRxn,2:end));
        else
            dataTmp = zeros(1,size(dataInfants,2)-1);
        end
        % add the adult samples
        findInM = find(strcmp(dataMothers(:,1),combFeats{j,1}));
        if ~isempty(findInM)
            dataM = cell2mat(dataMothers(findInM,2:end));
        else
            dataM = zeros(1,size(dataMothers,2)-1);
        end
        dataTmp = [dataTmp,dataM];
        data{j+1,1} = combFeats{j};
        data(j+1,2:length(dataTmp)+1) = num2cell(dataTmp);
    end
    modNames = dataInfants(1,2:end);
    modNames = [modNames,dataMothers(1,2:end)];
    data(1,2:length(dataTmp)+1) = modNames;
    
    if i==2
        % translate metabolites
        for j=2:size(data,1)
            data{j,1} = strrep(data{j,1},'EX_','');
            data{j,1} = strrep(data{j,1},'[fe]','');
            data{j,1} = database.metabolites{find(strcmp(database.metabolites(:,1),data{j,1})),2};
        end
    end

    timePoints = {'5_days','1_month','6_months','1_year'};
    sigFeats = {};
    for j=1:length(timePoints)
        stats = readInputTableForPipeline([rootDir filesep 'Statistical_analysis_with_maternalGut_COSMIC' filesep datasets{i,3} timePoints{j} '.csv']);
        stats(1,:) = [];
        sig = find(cell2mat(stats(:,7))<0.001);
        sigFeats = union(sigFeats,stats(sig,1));
    end
     
    % remove nonsignificant features
    [~,I] = setdiff(data(:,1),sigFeats,'stable');
    data(I(2:end),:) = [];
    
    % keep 50 highest
%     delArray = [];
%     for j=1:size(data,1)
%         delArray(j,1) = sum(cell2mat(data(j,2:end)));
%     end
%     [I,J] = sort(delArray,'descend');
%     J(find(J==1))=[];
%     data(J(51:end),:) = [];

    writetable(cell2table(data),[rootDir filesep 'Exported_results' filesep 'InfantsMothers_' datasets{i,2}],'writeVariableNames',false)
end

