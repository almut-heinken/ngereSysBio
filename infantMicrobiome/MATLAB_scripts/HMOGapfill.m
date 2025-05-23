function [model, addedRxns] = HMOGapfill(model, microbeID, database, inputDataFolder)
% Gap-fills carbon source utilization pathways in a microbial reconstruction based on
% experimental evidence.
%
% USAGE
%    [model, addedRxns] = HMOGapfill(model, microbeID, database)
%
% INPUT
% model             COBRA model structure
% microbeID:        ID of the reconstructed microbe that serves as the
%                   reconstruction name and to identify it in input tables
% database          rBioNet reaction database containing min. 3 columns:
%                   Column 1: reaction abbreviation, Column 2: reaction
%                   name, Column 3: reaction formula.
% inputDataFolder   Folder with experimental data and database files
%                   to load
%
% OUTPUT
% model             COBRA model structure refined through experimental data
%                   for carbon sources
% addedRxns         List of reactions that were added during refinement
%
% .. Authors:
% Almut Heinken and Rola Shaaban, 02/2023

addedRxns={};

HMOTable = readInputTableForPipeline([inputDataFolder filesep 'HMOTable.txt']);
HMOTable(:,find(strncmp(HMOTable(1,:), 'Ref',3)):end)=[];

mInd = find(ismember(HMOTable(:, 1), microbeID));
if isempty(mInd)
    warning(['Microbe ID not found in HMO data table: ', microbeID])
else
    
    % pathway, rxns to add
    gapfillAdd = {
        '2-Fucosyllactose', {'EX_2fuclac(e)', '2FUCLAC_FUCASEe', 'EX_fuc_L(e)', 'EX_lcts(e)'}
        '3-Fucosyllactose', {'EX_3fuclac(e)', '3FUCLAC_FUCASEe', 'EX_fuc_L(e)', 'EX_lcts(e)'}
        '3-Sialyllactose', {'EX_3slac(e)', '3SLAC_SIASEe', 'EX_acnam(e)', 'EX_lcts(e)'}
        '6-Sialyllactose', {'EX_6slac(e)', '6SLAC_SIASEe', 'EX_acnam(e)', 'EX_lcts(e)'}
        'Lacto-N-difucohexaose', {'EX_lacndfuchx(e)', 'LACNDFUCHX_FUCASEe', 'EX_fuc_L(e)', 'EX_lacnfucpt(e)'}
        'Lacto-N-tetraose', {'EX_lacnttr(e)'}
        'Lacto-N-neotetraose', {'EX_lacnnttr(e)', 'LACNNTTR_DEGe', 'EX_gal(e)', 'EX_lntri_ii(e)'}
        'Lactodifucotetraose', {'EX_lacdfucttr(e)', 'LACDFUCTTR_FUCASEe', 'EX_fuc_L(e)', 'EX_2fuclac(e)'}
        'Lacto-N-biose', {'EX_lctnb(e)','GAL6PI'}
        'Galacto-N-biose', {'EX_glctnb(e)','GAL6PI'}
        'Lacto-N-fucopentaose I', {'EX_lacnfucpt_i(e)', 'LACNFUCPT_I_FUCASEe', 'EX_lacnttr(e)', 'EX_fuc_L(e)'}
        'Lacto-N-fucopentaose V', {'EX_lacnfucpt_v(e)', 'LACNFUCPT_V_FUCASEe', 'EX_lacnttr(e)', 'EX_fuc_L(e)'}
        'Lacto-N-fucopentaose II', {'EX_lacnfucpt_ii(e)', 'LACNFUCPT_II_FUCASEe', 'EX_lacnttr(e)', 'EX_fuc_L(e)'}
        'Lacto-N-fucopentaose III', {'EX_lacnfucpt_iii(e)', 'LACNFUCPT_III_FUCASEe', 'EX_lacnttr(e)', 'EX_fuc_L(e)'}
        'Disialyllacto-N-tetraose', {'EX_dslacnttr(e)', 'DSLACNTTR_SIASEe', 'EX_acnam(e)', 'EX_lacnttr(e)'}
        'Lacto-N-triose I', {'EX_lntri_i(e)', 'LNTRI_I_DEGe', 'EX_acgam(e)', 'EX_lcts(e)'}
        'Lacto-N-triose II', {'EX_lntri_ii(e)', 'LNTRI_II_DEGe', 'EX_acgam(e)', 'EX_lcts(e)'}
        'Lacto-N-hexaose', {'EX_lacnhx(e)', 'LACNHX_DEGe', 'EX_lcts(e)', 'EX_gal(e)', 'EX_lctnb(e)', 'EX_acgam(e)'}
        'Difucosyllacto-N-hexaose a', {'EX_dflnh_a(e)', 'DFLNH_a_FUCASEe', 'EX_fuc_L(e)'}
        'Monofucosyllacto-N-hexaose', {'EX_fuclachx(e)', 'FUCLACHX_FUCASEe', 'EX_fuc_L(e)'}
        'Sialyllacto-N-tetraose a', {'EX_neulacnttr(e)', 'NEULACNTTR_SIASEe', 'EX_acnam(e)', 'EX_lacnttr(e)'}
        'Sialyllacto-N-tetraose b', {'EX_neulacnttr_b(e)', 'NEULACNTTR_B_SIASEe', 'EX_acnam(e)', 'EX_lacnttr(e)'}
        'Sialyllacto-N-tetraose c', {'EX_neulacnttr_c(e)', 'NEULACNTTR_C_SIASEe', 'EX_acnam(e)', 'EX_lacnttr(e)'}
        'Lacto-N-neohexaose', {'EX_lacnnhx(e)', 'LACNNHX_DEGe', 'EX_gal(e)', 'EX_acgam(e)', 'EX_lcts(e)'}
        'Lacto-N-difucopentaose II', {'EX_lacndfucpt_ii(e)', 'LACNDFUCPT_II_FUCASEe', 'EX_lacnttr(e)', 'EX_fuc_L(e)'}
        'Difucosyllactose', {'EX_dfuclac(e)', 'DFUCLAC_FUCASEe', 'EX_fuc_L(e)', 'EX_lcts(e)'}
        'Disialyllactose', {'EX_dslac(e)', 'DSLAC_SIASEe', 'EX_acnam(e)', 'EX_lcts(e)'}
        'Monofucosylmonosialyllacto-N-hexaose', {'EX_fucneulacnhx(e)', 'FUCNEULACNHX_DEGe', 'EX_acnam(e)', ...
        'EX_acgam(e)', 'EX_fuc_L(e)', 'EX_gal(e)', 'EX_lacnttr(e)', 'EX_lacnnttr(e)'}
        % simple sugars and disaccharides
        'Lactose', {'EX_lcts(e)', 'GALK', 'UGLT', 'UDPG4E', 'PGMT', ...
        'HEX1', 'PGI', 'PFK', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'}
        'N-acetyllactosamine', {'EX_naclacn(e)', '6PGALSZ3', 'GAL6PI', ...
        'PFK_2', 'ACGAMK', 'AGDC', 'G6PDA', 'TGBPA'}
        'N-acetylgalactosamine', {'ACGALK3', 'AGDC2', 'GALAM6PDA', 'PFK_2', 'TGBPA'}
        'N-acetylneuraminic acid', {'EX_acnam(e)', 'ACNAMabc', 'ACNML', 'AMANK', ...
        'AMANAPEr', 'AGDC', 'G6PDA'}
        'D-glucose', {'EX_glc_D(e)', 'GLCabc', 'HEX1', 'PGI', 'PFK', 'TPI', ...
        'GAPD', 'PGK', 'PGM', 'ENO', 'PYK'}
        'N-acetylglucosamine', {'EX_acgam(e)', 'ACGAMtr2', 'ACGAMK', 'AGDC', 'G6PDA'}
        'L-fucose', {'EX_fuc_L(e)', 'FUCt2_1', 'FCI', 'FCLK', 'FCLPA', 'LCARS', ...
        'EX_12ppd_S(e)', '12PPDt'}
        };
    
    % find the pathways to add
    pathways=HMOTable(1,2:end);
    if contains(version, '(R202') % for Matlab R2020a and newer
        cSources = pathways(find(cell2mat(HMOTable(mInd, 2:end)) == 1));
    else
        cSources = pathways(find(str2double(HMOTable(mInd, 2:end)) == 1));
    end
    
    gapfillAddConditional = {
        % c-source, condition, add reaction(s)
        'Lacto-N-biose', 'any(ismember(model.rxns, ''AMANK'')) && any(ismember(model.rxns, ''AGDC''))', {'AMANAPEr'}
        'Lacto-N-biose', 'strncmp(''Bifidobacterium_longum'', microbeID, 22)', {'ACGAMK'}
        'Lacto-N-biose', '~strncmp(''Lactobacillus'', microbeID, 13)', {'LNBP'}
        'Galacto-N-biose', '~strncmp(''Lactobacillus'', microbeID, 13)', {'GLCNBPTc'}
        'Lacto-N-biose', 'strncmp(''Lactobacillus'', microbeID, 13)', {'LCTNBpts', 'LCTNBHc'}
        'Galacto-N-biose', 'strncmp(''Lactobacillus'', microbeID, 13)', {'GLCTNBpts', 'GLCNBHc'}
        'N-acetyllactosamine', 'strncmp(''Lactobacillus'', microbeID, 13)', {'NACLACNpts'}
        'N-acetyllactosamine', 'strncmp(''Bifidobacterium'', microbeID, 15)', {'NACLACNabc', 'NACLACNK'}
        'Lacto-N-biose', '~strncmp(''Bifidobacterium_longum_infantis'', microbeID, 31) && ~strncmp(''Bifidobacterium_longum_subsp_infantis'', microbeID, 37) && ~strncmp(''Lactobacillus'', microbeID, 13)', {'EX_lctnb(e)', 'LCTNBabc'}
        'Galacto-N-biose', 'strncmp(''Bifidobacterium'', microbeID, 15) || strncmp(''Roseburia'', microbeID, 9)', {'EX_glctnb(e)', 'GLCTNBabc'}
        'Lactose', 'any(ismember(model.rxns, ''LACpts''))', {'6PGALSZ', 'GAL6PI', 'PFK_2', 'TGBPA'}
        'Lactose', '~any(ismember(model.rxns, ''LACpts''))', {'LACZ','LCTSabc'}
        'Lactose', 'strncmp(''Bifidobacterium_bifidum'', microbeID, 23)', {'LACZe', 'EX_gal(e)', 'EX_glc_D(e)'}
        'Lacto-N-hexaose', 'strncmp(''Bifidobacterium_longum_infantis'', microbeID, 31) || strncmp(''Bifidobacterium_longum_subsp_infantis'', microbeID, 37)', {'LACNHX_DEGc', 'LACNHXabc'}
        'Lacto-N-triose II', 'strncmp(''Bifidobacterium_breve'', microbeID, 21) || strncmp(''Bifidobacterium_longum'', microbeID, 22)', {'LNTRI_II_DEGc'}
        'Lacto-N-tetraose', '~strncmp(''Bifidobacterium_breve'', microbeID, 21) && ~strncmp(''Bifidobacterium_longum_infantis'', microbeID, 31) && ~strncmp(''Bifidobacterium_longum_subsp_infantis'', microbeID, 37)', {'LACNTTR_DEGe', 'EX_lctnb(e)','EX_lcts(e)'}
        'Lacto-N-tetraose', 'strncmp(''Bifidobacterium_breve'', microbeID, 21)', {'LACNTTR_DEG_BB', 'LACNTTRabc'}
        'Lacto-N-tetraose', 'strncmp(''Bifidobacterium_longum_infantis'', microbeID, 31) || strncmp(''Bifidobacterium_longum_subsp_infantis'', microbeID, 37)', {'LACNTTR_DEG_BI', 'EX_lntri_ii(e)'}
        'Lacto-N-neotetraose', 'strncmp(''Bifidobacterium_breve'', microbeID, 21) || strncmp(''Bifidobacterium_longum'', microbeID, 22)', {'LACNNTTRabc', 'LACNNTTR_DEGc'}
        };
    
    % go through carbon sources
    for i = 1:length(cSources)
        % add pathway reactions
        addRxns = gapfillAdd{find(ismember(gapfillAdd(:, 1), cSources{i})), 2};
        for j = 1:length(addRxns)
            if ~any(ismember(model.rxns, addRxns{j}))
                formula = database.reactions{ismember(database.reactions(:, 1), addRxns{j}), 3};
                model = addReaction(model, addRxns{j}, 'reactionFormula', formula, 'geneRule', '');
                addedRxns{length(addedRxns)+1,1} = addRxns{j};
            end
        end
        % add conditional reactions
        if any(ismember(gapfillAddConditional(:, 1), cSources{i}))
            conditions = find(ismember(gapfillAddConditional(:, 1), cSources{i}));
            for k = 1:length(conditions)
                if eval(gapfillAddConditional{conditions(k), 2})
                    addRxns = gapfillAddConditional{conditions(k), 3};
                    addRxns = setdiff(addRxns,model.rxns);
                    for j = 1:length(addRxns)
                        formula = database.reactions{ismember(database.reactions(:, 1), addRxns{j}), 3};
                        model = addReaction(model, addRxns{j}, 'reactionFormula', formula, 'geneRule', 'CarbonSourceGapfill');
                        addedRxns{length(addedRxns)+1,1} = addRxns{j};
                    end
                end
            end
        end
    end
    
    % add reactions with genes
    [model,addAnnRxns,updateGPRCnt]=refineGenomeAnnotation(model,microbeID,database,inputDataFolder);
    addedRxns = union(addedRxns,addAnnRxns);
    model = rebuildModel(model,database);
end

end
