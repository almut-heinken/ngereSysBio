function [model,flux] = add1cmPathways(model,species)
% Adds missing one-carbon metabolism reactions and genes to the human or 
% mouse reconstruction.

% remove the more generic versions of reactions/ replace with correct gene
% RE2493C replaced with METS_R
toRemove = {'METS','MMMm','CBLATm','CBL2OR','RE2493C','AQCOBALt','EX_adpcbl[e]','UNK3r','HMR_9550'};
model = removeRxns(model,toRemove);

% define the reactions with associated gene to add
reactions_genes = {
    'CCBLR','25974.1 AND 27249.1','67096 AND 109129','Vitamin B12 metabolism' % CblC and CblD-form complex
    'METS_FORM','25974.1 AND 27249.1 AND 4548.1','67096 AND 109129 AND 238505','Methionine and cysteine metabolism' % CblC and CblD and CblG
    'METS_CAT','25974.1 AND 27249.1 AND 4548.1','67096 AND 109129 AND 238505','Methionine and cysteine metabolism' % CblC and CblD and CblG
    'METS_R_REG','25974.1 AND 27249.1 AND 4548.1 AND 4552.1','67096 AND 109129 AND 238505 AND 210009','Methionine and cysteine metabolism' % CblC and CblD and CblG and CblE
    'COBOXs','','','Vitamin B12 metabolism' % spontaneous
    'CYANt','','','Transport, extracellular' % allowing metabolite to leave system
    'EX_cyan[e]','','','Exchange/demand reaction' % allowing metabolite to leave system
    'CBL2OR','','','Vitamin B12 metabolism' % seems to be spontaneous -> https://www.genome.jp/dbget-bin/www_bget?ec:1.16.1.3
    'MMMm_cbl','4594.1','17850','Valine, leucine, and isoleucine metabolism' % MUT
    'CBL_REL_MUT','166785.1','109136','Methionine and cysteine metabolism' % CblA
    'CBLATm','326625.1','77697','Vitamin B12 metabolism' % CblB
    'CBL2ATm','326625.1','77697','Vitamin B12 metabolism' % CblB
    'EX_C06453[e]','','','Exchange/demand reaction' % allowing metabolite to enter system
    'EX_ccbl[e]','','','Exchange/demand reaction' % allowing metabolite to enter system
    'EX_hoxocbl[e]','','','Exchange/demand reaction' % allowing metabolite to enter system
    'EX_aqcobal[e]','','','Exchange/demand reaction' % allowing metabolite to enter system
    'CBL1tm','25974.1 AND 27249.1','67096 AND 109129','Transport, mitochondrial' % CblC and CblD-assumed
    'CBL2tm','25974.1 AND 27249.1','67096 AND 109129','Transport, mitochondrial' % CblC and CblD-form complex
    'C06453tly','51293.1 AND 6948.1','54219 AND 21452','Transport, lysosomal' % CD320 and TCN2
    'C06453tl','55788.1 AND 5826.1','68421 AND 19300','Transport, lysosomal' % CblF and CblJ
    'CCBLtly','51293.1 AND 6948.1','54219 AND 21452','Transport, lysosomal' % CD320 and TCN2
    'CCBLtl','55788.1 AND 5826.1','68421 AND 19300','Transport, lysosomal' % CblF and CblJ
    'HOXOCBLtly','51293.1 AND 6948.1','54219 AND 21452','Transport, lysosomal' % CD320 and TCN2
    'HOXOCBLtl','55788.1 AND 5826.1','68421 AND 19300','Transport, lysosomal' % CblF and CblJ
    'AQCOBALtly','51293.1 AND 6948.1','54219 AND 21452','Transport, lysosomal' % CD320 and TCN2
    'AQCOBALtl','55788.1 AND 5826.1','68421 AND 19300','Transport, lysosomal' % CblF and CblJ
    'HOXOCBLC_C','','','Vitamin B12 metabolism' % gene unknown
    'PPA2m','','','Miscellaneous' % from Recon3D
    'CBLTDe','','','Transport, extracellular' % from Recon3D
    'HMR_7160','5980.1 OR 23649.1 OR 10721.1 OR 5424.1 OR 5423.1 OR 10714.1 OR 5427.1 OR 11201.1 OR 5422.1 OR 5425.1 OR 5985.1 OR 11044.1 OR 56655.1 OR 51426.1 OR 27434.1 OR 353497.1 OR 10514.1 OR 5428.1','(19714 or 18969 or 77782 or 18971 or 18970 or 67967 or 18974 or 26447 or 18968 or 18972 or 72151 or 210106 or 66979 or 27015 or 54125 or 272158 or 18432 or 18975 or 59001 or 56626 or 80905 or 69745 or 18973 or 50776 or 214627)','Nucleotide metabolism' % to account for DNA methylation
    'HMR_8639','','','Transport, nuclear' % to account for DNA methylation
    'DNAMTn','1787.3 OR 1788.2 OR 1786.1 OR 1787.4 OR 1787.6 OR 1788.3 OR 1789.4 OR 1789.2 OR 1789.1 OR 1787.5 OR 1789.3 OR 1787.1 OR 1788.4 OR 1788.1 OR 1787.2','(13434 or 13435 or 13433 or 13436)','Methionine and cysteine metabolism' % to account for DNA methylation
    'AMETtn','','','Transport, nuclear' % to account for DNA methylation
    'AHCYStn','','','Transport, nuclear' % to account for DNA methylation
    'DM_dna5mtc[n]','','','Exchange/demand reaction' % needed for model functionality
    'DM_con_adocbl[m]','','','Exchange/demand reaction' % needed for model functionality
    'DM_con_C06453[c]','','','Exchange/demand reaction' % needed for model functionality
    'DM_dna[n]','','','Exchange/demand reaction' % needed for model functionality
    'PPCOACm','5096.1 and 5095.1','66904 and 110821','Valine, leucine, and isoleucine metabolism' % to capture propionate metabolism
    };

% define metabolites to add
metabolites = {'adocbl[m]';'adocbl[c]';'adocbl[e]';'ahcys[n]';'amet[n]';'aqcobal[c]';'aqcobal[l]';'C06453[c]';'C06453[e]';'C06453[l]';'cbl1[m]';'cbl1[c]';'ccbl[c]';'ccbl[e]';'ccbl[l]';'con_adocbl[m]';'con_C06453[c]';'dna[c]';'dna[n]';'dna5mtc[n]';'hoxocbl[e]';'hoxocbl[l]';'hoxocbl[c]';'mets[c]';'mets_C06453[c]';'mets_cbl1[c]';'mets_cbl2[c]';'mutApoE[m]';'mutHoloE[m]';'pppi[m]'};

database = loadVMHDatabase;

% add metabolites
model = removeMetabolites(model, metabolites, false);
for i=1:length(metabolites)
    met = strsplit(metabolites{i},'[');
    findMet = find(strcmp(database.metabolites(:,1),met{1}));
    model = addMetabolite(model,metabolites{i},'metName',database.metabolites{findMet,2},'metFormula',database.metabolites{findMet,4},'Charge',str2double(database.metabolites{findMet,5}));
end

% add reactions
for j = 1:length(reactions_genes)
    formula = database.reactions{ismember(database.reactions(:, 1), reactions_genes{j,1}), 3};
    if strcmp(species,'Human')
        model = addReaction(model, reactions_genes{j,1}, 'reactionFormula', formula, 'geneRule', reactions_genes{j,2},'subSystem', reactions_genes{j,4});
    elseif strcmp(species,'Mouse')
        model = addReaction(model, reactions_genes{j,1}, 'reactionFormula', formula, 'geneRule', reactions_genes{j,3},'subSystem', reactions_genes{j,4});
    end
    rxnID = find(strcmp(model.rxns,reactions_genes{j,1}));
    if ~ischar(model.subSystems{rxnID,1})
        model.subSystems{rxnID,1} = model.subSystems{rxnID,1}{1};
    end
    % add reaction description
    model.rxnNames{end,1} = database.reactions{ismember(database.reactions(:, 1), reactions_genes{j,1}), 2};
end

% test if reactions can carry flux
BlockedReaction = findBlockedReaction(model,'L2');
flux = setdiff(reactions_genes(:,1),BlockedReaction);

end
