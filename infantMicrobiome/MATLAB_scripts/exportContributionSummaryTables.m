
%% get microbes that contribute metabolites of interest in each microbiome

metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);

mets = {'ac','ppa','but','lac_L','fol','adocbl','ribflv','thm','nac','btn','pydx','pnto_R'};

% reference genomes
data =  readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv']);
data(:,1) = regexprep(data(:,1),'pan','','once');
data(1,:) = strrep(data(1,:),'microbiota_model_diet_','');
dataM =  readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Microbe_Secretion.csv']);
dataM(1,:) = strrep(dataM(1,:),'microbiota_model_diet_','');

timePoints = unique(metadata(2:end,3));

%% extract contributions for infants-summary tables

for i=1:length(mets)
    I = find(endsWith(data(:,1),['_' mets{i}]));
    dataMet = data(1,:);
    dataMet(2:length(I)+1,:) = data(I,:);
    for j=2:size(dataMet,1)
        findStr = strfind(dataMet{j,1},['_' mets{i}]);
        dataMet{j,1} = dataMet{j,1}(1:findStr(end)-1);
    end
    writetable(cell2table(dataMet),[rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Contributions_' mets{i} '.csv'],'writeVariableNames',false)
end

%% export contributions for adults-summary tables

for i=1:length(mets)
    I = find(endsWith(dataM(:,1),['_' mets{i}]));
    dataMet = dataM(1,:);
    dataMet(2:length(I)+1,:) = dataM(I,:);
    for j=2:size(dataMet,1)
        findStr = strfind(dataMet{j,1},['_' mets{i}]);
        dataMet{j,1} = dataMet{j,1}(1:findStr(end)-1);
    end
    writetable(cell2table(dataMet),[rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Contributions_' mets{i} '.csv'],'writeVariableNames',false)
end

%% export summary table of contributions by metabolite and taxon
mkdir([rootDir filesep 'ContributionsByMetabolite'])
metadata = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv']);

% read contributions to metabolites
microbeInfo = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'expanded_AGORA2_infoFile.xlsx']);
microbeInfo(:,1:4)=[];
microbeInfo(:,1)=strrep(microbeInfo(:,1),'[','');
microbeInfo(:,1)=strrep(microbeInfo(:,1),' ','_');
microbeInfo(:,1)=strrep(microbeInfo(:,1),']','');
microbeInfo(:,1)=strrep(microbeInfo(:,1),'(','_');
microbeInfo(:,1)=strrep(microbeInfo(:,1),')','');
microbeInfo(:,1)=strrep(microbeInfo(:,1),'/','_');
microbeInfo(:,1)=strrep(microbeInfo(:,1),'-','_');
microbeInfo(:,1)=strrep(microbeInfo(:,1),'.','');

% then get a summary on different taxon levels for each metabolite
taxLevels={'Genus','Family','Order','Class','Phylum'};
for i=1:length(mets)
    % infants
    dataMet = readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Contributions_' mets{i}]);
    % adults
    dataMetM = readInputTableForPipeline([rootDir filesep 'Modeling_MaternalGutMicrobiomes' filesep 'Contributions' filesep 'Contributions_' mets{i}]);

    for j=1:length(taxLevels)
        % prepare the table
        summary = {'Taxon','5 days','1 month','6 months','1 year','Maternal gut'};
        % first get all taxa in infants/adults
        taxCol=find(strcmp(microbeInfo(1,:),taxLevels{j}));
        taxa=unique(microbeInfo(2:end,taxCol));
        taxa(find(contains(taxa,'unclassified')),:)=[];
        taxa(find(strcmp(taxa,'')),:)=[];
        summary(2:length(taxa)+1,1)=taxa;
        summary(2:end,2:end)={0};
        % first summarize for infants
        % summarize the average contribution of species/strains belonging to each taxon
        for k=2:size(summary,1)
            % first infants
            rep=microbeInfo(find(strcmp(microbeInfo(:,taxCol),summary{k,1})),1);
            [C,I]=intersect(dataMet(:,1),rep);
            for l=1:length(I)
                for m=2:5
                    % get the samples for each time point
                    samp=metadata(find(strcmp(metadata(:,4),summary{1,m})),1);
                    [~,findSamp]=intersect(dataMet(1,:),samp);
                    summary(k,m)=num2cell(cell2mat(summary(k,m)) + mean(cell2mat(dataMet(I(l),findSamp))));
                end
            end
        end
        % then adults
        % summarize the average contribution of species/strains belonging to each taxon
        for k=2:size(summary,1)
            rep=microbeInfo(find(strcmp(microbeInfo(:,taxCol),summary{k,1})),1);
            [C,I]=intersect(dataMetM(:,1),rep);
            for l=1:length(I)
                summary(k,6)=num2cell(cell2mat(summary(k,6)) + mean(cell2mat(dataMetM(I(l),2:end))));
            end
        end
        % remove taxa with no/small contributions
        cnt=1;
        delArray=[];
        for k=2:size(summary,1)
            if sum(cell2mat(summary(k,2:end)))<sum(sum(cell2mat(summary(2:end,2:end))))/100 || sum(cell2mat(summary(k,2:end)))<0.000001
                delArray(cnt)=k;
                cnt=cnt+1;
            end
        end
        summary(delArray,:)=[];
        writetable(cell2table(summary),[rootDir filesep 'ContributionsByMetabolite' filesep mets{i} '_' taxLevels{j} '.csv'],'writeVariableNames',false)

        if j<4
            % create microbe annotation
            taxonomy={'Taxon','Phylum'};
            taxonomy(2:size(summary,1),1)=summary(2:end,1);
            phylCol=find(strcmp(microbeInfo(1,:),'Phylum'));
            for k=2:size(taxonomy,1)
                taxonomy{k,2}=microbeInfo{find(strcmp(microbeInfo(:,taxCol),taxonomy{k,1})),phylCol};
            end
            cell2csv([rootDir filesep 'ContributionsByMetabolite' filesep mets{i} '_' taxLevels{j} '_annotations.csv'],taxonomy)
        end
    end
end
