
% run integration of WT vs. MTR mutant proteomic data into mouse
% reconstruction of brain

%% load reconstruction
modelB=readCbModel('Data/iBrain674_Mm.xml');
modelB=creategrRulesField(modelB);
modelB=changeObjective(modelB,'n_Ex_A_Macro');
% get flux consistent subset
[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, model, modelB] = findFluxConsistentSubset(modelB);

% get some reaction annotations in the model
vmhDB = readInputTableForPipeline('Data/recon-store-reactions-1.tsv');
info = readInputTableForPipeline('Data/d0mo00135j2.xlsx');
info(:,1)=strrep(info(:,1),'-','_');
modelB.subSystems = cell(length(modelB.rxns),1);
for i=1:length(modelB.rxns)
    rxn=find(strcmp(info(:,1),strrep(modelB.rxns{i},'n_','')));
    ECnumbers = strsplit(info{rxn,6},'|');
    ECnumbers = strrep(ECnumbers,'EC-','');
    ECnumbers = strrep(ECnumbers,' ','');
    % find Ec numbers in global reconstruction
    [C,I] = intersect(vmhDB(:,5),ECnumbers);
    if ~isempty(I) && ~isempty(C{1})
%         modelB.rxnNames{i,1}=vmhDB{I(1),2};
        modelB.subSystems{i,1}=vmhDB{I(1),4};
    else
        modelB.subSystems{i,1} = info{rxn,2};
    end
end

%% create mapping from proteins onto gene loci
% download latest mouse assembly from NCBI
websave('GCF_000001635.27_GRCm39_feature_table.txt.gz','https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_feature_table.txt.gz')
gunzip('GCF_000001635.27_GRCm39_feature_table.txt.gz')

geneMapping = readInputTableForPipeline('GCF_000001635.27_GRCm39_feature_table.txt');

%% create input parameters for data integration
%% brain-specific model

mkdir('Brain_specific_models')

proteome = readInputTableForPipeline('Data/Proteomic_data.xlsx');
proteome(find(cellfun(@isempty,proteome(:,3))),:)=[];

% map gene symbols to gene IDs
genesMapped = {'GeneSymbol','GeneID'};
cnt=2;
for i=2:size(proteome,1)
    genesMapped{cnt,1} = proteome{i,3};
    findProt = find(strcmp(geneMapping(:,15),proteome{i,3}));
    if ~isempty(findProt)
        genesMapped{cnt,2}=geneMapping{findProt(1),16};
        cnt=cnt+1;
    end
end

% adapt the input table with the proteome data
proteomeTable = proteome;
[~,I] = setdiff(proteomeTable(:,3),genesMapped(:,1),'stable');
proteomeTable(I(2:end),:) = [];
proteomeTable(:,1:2) = [];

for i=2:size(proteomeTable,1)
    proteomeTable{i,1} = genesMapped{find(strcmp(genesMapped(:,1),proteomeTable{i,1})),2};
end
delArray=[];
cnt=1;
for i=2:size(proteomeTable,1)
    if isempty(find(ismember(str2double(modelB.genes),proteomeTable{i,1})))
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
proteomeTable(delArray,:) = [];

p_values=cell2mat(proteomeTable(2:end,2));
proteomeTable(:,2:3) = [];

for i=2:size(proteomeTable,1)
    proteomeTable{i,1}=num2str(proteomeTable{i,1});
end

for i=2:size(genesMapped,1)
    genesMapped{i,2}=num2str(genesMapped{i,2});
end

% loop through the samples and create a sample-specific model with two
% methods

samples=readInputTableForPipeline('Data/SampleAnnotations.xlsx');

% use INIT and then EFlux
% there may be better methods but let us use this for now
for i=2:size(proteomeTable,2)
    % process the proteomic data
    expressionData = struct;
    expressionData.gene = proteomeTable(2:end,1);
    expressionData.value = cell2mat(proteomeTable(2:end,i));
    %     expressionData.sig = p_values;
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(modelB, expressionData);

    % create options variable for INIT
    options=struct;
    options.solver='INIT';
    options.weights=expressionRxns;
    options.weights(find(isnan(options.weights)))=0;
    % add objective function inclusion
    options.weights(find(strcmp(modelB.rxns(:,1),'n_Ex_A_Macro')),1)=max(options.weights);
    % add Mtr reaction
    %% maybe not best way!!!
    if strcmp(samples{find(strcmp(samples(:,1),proteomeTable{1,i})),2},'WT')
        options.weights(find(strcmp(modelB.rxns(:,1),'n_RA_meth4_c')),1)=100;
        options.weights(find(strcmp(modelB.rxns(:,1),'n_RN_meth4_c')),1)=100;
    elseif strcmp(samples{find(strcmp(samples(:,1),proteomeTable{1,i})),2},'KO')
        options.weights(find(strcmp(modelB.rxns(:,1),'n_RA_meth4_c')),1)=20;
        options.weights(find(strcmp(modelB.rxns(:,1),'n_RN_meth4_c')),1)=20;
    end
    %%

    model = createTissueSpecificModel(modelB, options);

    % now scale reactions with EFlux
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData);
    expression=struct;
    expression.target=model.rxns;
    expression.value=expressionRxns;
    expression.value(find(isnan(expression.value)))=-1;
    expression.value(find(strcmp(expression.target(:,1),'n_Ex_A_Macro')),1)=-1;
    % add Mtr reaction
    if strcmp(samples{find(strcmp(samples(:,1),proteomeTable{1,i})),2},'WT')
        expression.value(find(strcmp(expression.target,'n_RA_meth4_c')),1)=100;
        expression.value(find(strcmp(expression.target,'n_RN_meth4_c')),1)=100;
    elseif strcmp(samples{find(strcmp(samples(:,1),proteomeTable{1,i})),2},'KO')
        expression.value(find(strcmp(expression.target,'n_RA_meth4_c')),1)=20;
        expression.value(find(strcmp(expression.target,'n_RN_meth4_c')),1)=20;
    end
    expression.preprocessed=true;
    % modified version of original e-Flux method
    model = relaxedApplyEFluxConstraints(model, expression);
    writeCbModel(model,'format','mat','fileName',['Brain_specific_models' filesep proteomeTable{1,i}])
end
