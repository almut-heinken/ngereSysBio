
% integrate fecal metabolic data available for 16 patients

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
solverOK=changeCobraSolver('ibm_cplex','QP');

% read fecal metabolomics data
metabolome = readInputTableForPipeline(['input' filesep 'Fecal_Metabolome.csv']);

% read metadata to match the samples
metadata = readInputTableForPipeline(['input' filesep 'SampleMetadataPatients.csv']);

% first map units
for i=2:size(metabolome,1)
    if ~isempty(find(contains(metabolome{i,1},'μg/g')))
        for j=2:size(metabolome,2)
            metabolome{i,j} = metabolome{i,j}/1000;
        end
    end
end

% translate names of metabolites to corresponding model IDs
% mapping = {'CA (μg/g dry feces)','cholate';'3β-DCA (μg/g dry feces)','NA';'DCA (μg/g dry feces)','dchac';'LCA (μg/g dry feces)','HC02191';'7-KLCA (μg/g dry feces)','NA';'GCDCA (μg/g dry feces)','dgchol';'GDCA (μg/g dry feces)','NA';'HCA (μg/g dry feces)','NA';'UDCA (μg/g dry feces)','HC02194';'TLCA (μg/g dry feces)','NA';'5-cholenic (μg/g dry feces)','NA';'iLCA (μg/g dry feces)','NA';'β-MCA (μg/g dry feces)','NA';'3-DHCA (μg/g dry feces)','3dhchol';'MDCA (μg/g dry feces)','NA';'CDCA (μg/g dry feces)','C02528';'Acetic acid (mg/g dry faeces)','ac';'Propionic acid (mg/g dry faeces)','ppa';'Butyric acid (mg/g dry faeces)','but';'Isobutyric acid (mg/g dry faeces)','isobut';'Valeric acid (mg/g dry faeces)','M03134';'Isovaleric acid (mg/g dry faeces)','isoval';'2-methylbutyric acid (mg/g dry faeces)','NA'};
% without some metabolites that cannot always be produced
mapping = {'CA (μg/g dry feces)','cholate';'3β-DCA (μg/g dry feces)','NA';'DCA (μg/g dry feces)','dchac';'LCA (μg/g dry feces)','NA';'7-KLCA (μg/g dry feces)','NA';'GCDCA (μg/g dry feces)','NA';'GDCA (μg/g dry feces)','NA';'HCA (μg/g dry feces)','NA';'UDCA (μg/g dry feces)','HC02194';'TLCA (μg/g dry feces)','NA';'5-cholenic (μg/g dry feces)','NA';'iLCA (μg/g dry feces)','NA';'β-MCA (μg/g dry feces)','NA';'3-DHCA (μg/g dry feces)','3dhchol';'MDCA (μg/g dry feces)','NA';'CDCA (μg/g dry feces)','C02528';'Acetic acid (mg/g dry faeces)','ac';'Propionic acid (mg/g dry faeces)','ppa';'Butyric acid (mg/g dry faeces)','but';'Isobutyric acid (mg/g dry faeces)','isobut';'Valeric acid (mg/g dry faeces)','NA';'Isovaleric acid (mg/g dry faeces)','isoval';'2-methylbutyric acid (mg/g dry faeces)','NA'};
for i=2:size(metabolome,1)
    metabolome{i,1} = mapping{find(strcmp(mapping(:,1),metabolome{i,1})),2};
end

% remove metabolites not in models
metabolome(find(strcmp(metabolome(:,1),'NA')),:) = [];

% load the models, implement constraints, and save new versions
mkdir([pwd filesep 'MetabolomeConstrainedModels'])
mkdir([pwd filesep 'MetabolomeConstrainedModels' filesep 'Models'])
mkdir([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes'])

for i=2:size(metabolome,2)
    modelToLoad = metadata{find(strcmp(metadata(:,7),metabolome{1,i})),1};
    model = readCbModel([pwd filesep 'MicrobiomeModels' filesep 'Diet' filesep 'microbiota_model_diet_' modelToLoad '.mat']);
    % formulate a pseudo-reaction for metabolomics data
    form = '';
    for j=2:size(metabolome,1)-1
        form = [form num2str(metabolome{j,i}) ' ' metabolome{j,1} '[u] + '];
    end
    form = [form num2str(metabolome{end,i}) ' ' metabolome{end,1} '[u] -> metabolome[u]'];
    model = addReaction(model,'metabolomeReaction',form);
    model = addDemandReaction(model,'metabolome[u]');
    model = changeObjective(model,'metabolomeReaction');
    model = changeRxnBounds(model,'communityBiomass',1,'u');

    FBA=optimizeCbModel(model,'max','1e-6')
    save([pwd filesep 'MetabolomeConstrainedModels' filesep 'Models' filesep modelToLoad '.mat'],'model')
    save([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes' filesep 'Fluxes_' modelToLoad '.mat'],'FBA')
end

% extract the fluxes
% first get all reactions
rxns = {};
for i=2:size(metabolome,2)
    modelToLoad = metadata{find(strcmp(metadata(:,7),metabolome{1,i})),1};
    model = readCbModel([pwd filesep 'MetabolomeConstrainedModels' filesep 'Models' filesep modelToLoad '.mat']);
    rxns = union(model.rxns,rxns);
end

fluxTable = vertcat('Reaction',rxns);

for i=2:size(metabolome,2)
    modelToLoad = metadata{find(strcmp(metadata(:,7),metabolome{1,i})),1};
    model = readCbModel([pwd filesep 'MetabolomeConstrainedModels' filesep 'Models' filesep modelToLoad '.mat']);
    flux = load([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes' filesep 'Fluxes_' modelToLoad '.mat']);
    fluxTable(1,i) = {modelToLoad};
    for j=2:size(fluxTable,1)
        findRxn = find(strcmp(model.rxns,fluxTable{j,1}));
        if ~isempty(findRxn)
            fluxTable{j,i} = flux.FBA.x(findRxn);
        else
            fluxTable{j,i} = 0;
        end
    end
end

delArray = [];
for i=2:size(fluxTable,1)
    if abs(sum(cell2mat(fluxTable(i,2:end)))) < 0.000000001
        delArray(length(delArray)+1) = i;
    end
end
fluxTable(delArray,:) = [];

% extract only the contributions
mappedMets = unique(mapping(:,2));
mappedMets(find(strcmp(mappedMets,'NA'))) = [];
for i=1:length(mappedMets)
    mappedMets{i,1} = ['IEX_' mappedMets{i,1} '[u]tr'];
end

delArray = [];
for i=2:size(fluxTable,1)
    if ~any(contains(fluxTable{i,1},mappedMets))
        delArray(length(delArray)+1) = i;
    end
end
fluxTable(delArray,:) = [];
writetable(cell2table(fluxTable),[pwd filesep 'MicrobiomeResults' filesep 'MetabolomeConstrainedExchangeFluxes.csv'],'writeVariableNames',false)

database = loadVMHDatabase;

fluxAnnotations = fluxTable(:,1);
fluxAnnotations(1,2:3) = {'Species','Metabolite'};
for i=2:size(fluxTable,1)
    sp = strsplit(fluxTable{i,1},'_IEX_');
    met = strrep(sp{2},'[u]tr','');
    fluxAnnotations{i,2} = regexprep(sp{1},'pan','','once');
    fluxAnnotations{i,3} = database.metabolites{find(strcmp(database.metabolites(:,1),met)),2};
    fluxAnnotations{i,3}(1) = upper(fluxAnnotations{i,3}(1));
end
writetable(cell2table(fluxAnnotations),[pwd filesep 'MicrobiomeResults' filesep 'MetabolomeConstrainedExchangeFluxes_Annotation.csv'],'writeVariableNames',false)
