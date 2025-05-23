
% calculate correlations with plasma lipidomics

samples = readInputTableForPipeline([pwd filesep 'Data' filesep 'Samples.csv']);

% calculate correlations and create plots
mkdir([pwd filesep 'Correlations'])

numWorkers=4;
if ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

%% lipidomics

% only use significant features
lipidomics = readInputTableForPipeline([pwd filesep 'significantFeatures_lipidomics.csv']);

% load all fluxes
data = readInputTableForPipeline([pwd filesep 'Metabolic_Flux_Results' filesep 'AllFluxes.csv']);
data(:,2:3)=[];

% remove fluxes that are all zeros
% remove fluxes with too small flux or that are nearly all the
% same
delArray=[];
cnt=1;
for j=2:size(data,1)
    if abs(sum(cell2mat(data(j,2:end))))<0.001
        delArray(cnt)=j;
        cnt=cnt+1;
    end
    for k=2:size(data,2)
        if length(find(cell2mat(data(j,2:end)) ~= data{j,k}))<=5
            delArray(cnt)=j;
            cnt=cnt+1;
        end
    end
end
delArray = unique(delArray);
data(delArray,:)=[];

% remove samples without data and vice versa
[C,I]=setdiff(data(1,:),lipidomics(1,:),'stable');
data(:,I(2:end))=[];
[C,I]=setdiff(lipidomics(1,:),data(1,:),'stable');
lipidomics(:,I(2:end))=[];

% correlations with lipidomics
correlations = {'Feature'};
pValues = {'Feature'};

sampFlux=data(1,:);
% loop through the datapoints one by one
for j=2:size(data,1)
    correlations{1,j}=data{j,1};
    pValues{1,j}=data{j,1};
    j
    dataTmp=data(j,:);
    % loop through the clinical parameters one by one
    % to parallelize
    corrTmp=zeros(1,size(lipidomics,1));
    pTmp=zeros(1,size(lipidomics,1));
    parfor k=2:size(lipidomics,1)
        % determine correlations
        lipidomicsTmp=lipidomics(k,:);
        corrDat = cell2mat(dataTmp(2:end))';
        for l=2:size(lipidomicsTmp,2)
            samp=find(strcmp(sampFlux(1,:),lipidomics{1,l}));
            corrDat(samp-1,2)=lipidomicsTmp{1,l};
        end
        % exclude samples with missing values
        corrDat(find(isnan(corrDat(:,2))),:)=[];
        % calculate the Spearman correlations
        [RHO, PVAL] = corr(corrDat(:, 1), corrDat(:, 2), 'type', 'Spearman');
        if isnan(RHO)
            RHO=0;
        end
        if isnan(PVAL)
            PVAL=1;
        end
        corrTmp(k)=RHO;
        pTmp(k)=PVAL;
    end
    for k=2:size(lipidomics,1)
        if j==2
            correlations{k,1} = lipidomics{k,1};
            pValues{k,1} = lipidomics{k,1};
        end
        correlations{k,j} = num2str(corrTmp(k));
        pValues{k,j} =num2str(pTmp(k));
    end
end

% perform correction for multiple testing
vals = [];
cnt=1;
for i=2:size(pValues,1)
    for j=2:size(pValues,2)
        vals(cnt,1)=str2double(pValues{i,j});
        cnt=cnt+1;
    end
end
fdr = mafdr(vals);
cnt=1;
for i=2:size(pValues,1)
    for j=2:size(pValues,2)
        pValues{i,j}=fdr(cnt,1);
        cnt=cnt+1;
    end
end

writetable(cell2table(correlations),[pwd filesep 'Correlations' filesep 'Correlations_Lipidomics_Fluxes'],'FileType','text','writeVariableNames',false,'Delimiter','tab')
writetable(cell2table(pValues),[pwd filesep 'Correlations' filesep 'pValues_Lipidomics_Fluxes'],'FileType','text','writeVariableNames',false,'Delimiter','tab')

delArray = [];
cnt=1;
for i=2:size(pValues,1)
    go = 0;
    for j=2:size(pValues,2)
        if pValues{i,j}<0.05
             go = 1;
        end
    end
    if go==0
        delArray(cnt,1) = i;
        cnt=cnt+1;
    end
end
correlations(delArray,:) = [];
pValues(delArray,:) = [];

delArray = [];
cnt=1;
for i=2:size(pValues,2)
    go = 0;
    for j=2:size(pValues,1)
        if pValues{j,i}<0.05
             go = 1;
        end
    end
    if go==0
        delArray(cnt,1) = i;
        cnt=cnt+1;
    end
end
correlations(:,delArray) = [];
pValues(:,delArray) = [];

annotations = correlations(1,:)';
annotations{1,2} = 'Subsystem';

for i=2:size(annotations,1)
    rxn = annotations{i,1};
    rxn = strrep(rxn,'_min','');
    rxn = strrep(rxn,'_max','');
    findRxn = find(strcmp(oriModel.rxns,rxn));
    annotations{i,2} = oriModel.subSystems{findRxn,1}{1};
end

for i=2:size(annotations,1)
    I = find(contains(annotations(:,2),annotations{i,2}));
    if length(I)<=5
        findRxn = find(strcmp(database.reactions(:,11),annotations{i,2}));
        if ~isempty(findRxn)
            annotations{i,2} = database.reactions{findRxn,12};
        else
            annotations{i,2} = 'Others';
        end
    end
    if strncmp(annotations{i,2},'Transport',9)
        annotations{i,2}='Transport';
    end
end
correlations(1,:) = strrep(correlations(1,:),'_min','_reverse');
correlations(1,:) = strrep(correlations(1,:),'_max','_forward');
annotations(:,1) = strrep(annotations(:,1),'_min','_reverse');
annotations(:,1) = strrep(annotations(:,1),'_max','_forward');

writetable(cell2table(correlations),[pwd filesep 'Correlations' filesep 'HighCorrelations_Lipidomics_Fluxes'],'FileType','text','writeVariableNames',false,'Delimiter','tab')
writetable(cell2table(annotations),[pwd filesep 'Correlations' filesep 'HighCorrelations_Lipidomics_Fluxes_Annotations'],'FileType','text','writeVariableNames',false,'Delimiter','tab')
