
% Create personalized retina models using the rFASTCORMICS and E-Flux
% algorithms.

%% get reconstruction
websave('41598_2020_63235_MOESM2_ESM.zip','https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-63235-w/MediaObjects/41598_2020_63235_MOESM2_ESM.zip')
unzip('41598_2020_63235_MOESM2_ESM.zip')

%% load reconstruction
model=readCbModel('iMM1865.xml');
model=creategrRulesField(model);
model=changeObjective(model,'BIOMASS_reaction');
% create consistent model
A = fastcc_4_rfastcormics(model, 1e-4,0);
model = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));

% make unlikely reversible reactions irreversible
unlIrrRxns = {'ALAARGCYSr';'ALAARGGLYr';'ALAASNLEUr';'ALAGLYLYSr';'ALAHISALAr';'ALALYSTHRr';'ARGALAALAr';'ARGALAPHEr';'ARGALATHRr';'ARGARGr';'ARGARGLYSr';'ARGARGMETr';'ARGCYSGLYr';'ARGCYSSERr';'ARGGLUGLUr';'ARGGLUPROr';'ARGGLYGLYr';'ARGHISTHRr';'ARGLEUPHEr';'ARGLYSASPr';'ARGPHEARGr';'ARGPROMETr';'ARGPROTHRr';'ARGSERSERr';'ARGTYRVALr';'ARGVALCYSr';'ARGVALTRPr';'ASNASNARGr';'ASNCYSCYSr';'ASNMETPROr';'ASNPHEASPr';'ASNPHECYSr';'ASNTYRGLYr';'ASNTYRPHEr';'ASNTYRTHRr';'ASPALAARGr';'ASPASNGLUr';'ASPGLUr';'ASPGLUPROr';'ASPGLUTRPr';'ASPHISCYSr';'ASPHISPROr';'ASPLYSGLUr';'ASPLYSHISr';'ASPMETASPr';'ASPPROLYSr';'ASPVALASNr';'CYSASNMETr';'CYSASPPHEr';'CYSCYSr';'CYSGLNMETr';'CYSGLUHISr';'CYSGLUTRPr';'CYSLEUTHRr';'CYSSERMETr';'CYSTYRASNr';'GLNASNGLNr';'GLNHISHISr';'GLNHISLYSr';'GLNLYSLYSr';'GLNLYSTRPr';'GLNPROGLUr';'GLNTRPGLUr';'GLNTYRLEUr';'GLUARGLEUr';'GLUASNLEUr';'GLUGLUr';'GLUILELYSr';'GLULEUr';'GLUMETr';'GLUMETHISr';'GLUTHRr';'GLUTHRLYSr';'GLUTRPALAr';'GLYHISASNr';'GLYHISLYSr';'GLYLYSCYSr';'GLYLYSPHEr';'GLYTYRLYSr';'GLYVALHISr';'HISARGCYSr';'HISARGSERr';'HISASPr';'HISCYSCYSr';'HISGLNALAr';'HISGLUr';'HISGLUGLNr';'HISGLYLYSr';'HISHISLYSr';'HISLYSALAr';'HISLYSGLUr';'HISLYSILEr';'HISLYSTHRr';'HISLYSVALr';'HISMETr';'HISMETGLNr';'HISPHEARGr';'HISPROLYSr';'HISTRPHISr';'ILEARGILEr';'ILEASNHISr';'ILEASPr';'ILEGLNGLUr';'ILEGLYARGr';'ILEPROLYSr';'ILESERARGr';'ILETRPTYRr';'LEUALAARGr';'LEUASNASPr';'LEUASPLYSr';'LEULEUTRPr';'LEUPROr';'LEUPROARGr';'LEUSERTRPr';'LEUTRPr';'LEUTRPARGr';'LEUTYRTYRr';'LEUVALr';'LYSARGLEUr';'LYSCYSHISr';'LYSGLNPHEr';'LYSGLUGLUr';'LYSLYSLYSr';'LYSPHEILEr';'LYSTRPARGr';'LYSTYRILEr';'LYSVALPHEr';'LYSVALTRPr';'METARGLEUr';'METASNTYRr';'METGLNTYRr';'METGLYARGr';'METHISLYSr';'METMETILEr';'METPHEARGr';'METTRPPHEr';'PHEASNMETr';'PHEASPr';'PHEGLNPHEr';'PHELEUr';'PHELEUASPr';'PHELEUHISr';'PHELYSALAr';'PHELYSPROr';'PHEPHEr';'PHEPHEASNr';'PHEPHETHRr';'PHEPROARGr';'PHESERTRPr';'PHETHRLYSr';'PHETRPLEUr';'PHETYRr';'PHETYRGLNr';'PHETYRLYSr';'PROARGASPr';'PROARGCYSr';'PROASNCYSr';'PROCYSr';'PROGLNPROr';'PROGLULYSr';'PROHISr';'PROHISTYRr';'PROLEUARGr';'PROLYSPROr';'PROPHEr';'PROPROARGr';'PROPROPROr';'PROTRPLYSr';'PROTRPTHRr';'PROVALGLNr';'SERARGALAr';'SERARGTRPr';'SERCYSARGr';'SERGLYGLUr';'SERLYSHISr';'SERPHELYSr';'SERTRPHISr';'THRARGTYRr';'THRASNTYRr';'THRGLNGLUr';'THRGLNTYRr';'THRHISHISr';'THRILEARGr';'THRMETARGr';'THRPHEARGr';'THRSERARGr';'THRTHRARGr';'THRTYRMETr';'TRPALAPROr';'TRPARGALAr';'TRPASPASPr';'TRPGLNGLNr';'TRPGLUGLYr';'TRPGLULEUr';'TRPGLUPROr';'TRPGLUTYRr';'TRPGLYLEUr';'TRPGLYPHEr';'TRPGLYVALr';'TRPHISMETr';'TRPILELYSr';'TRPILETRPr';'TRPLEUVALr';'TRPLYSr';'TRPMETARGr';'TRPMETVALr';'TRPPHEr';'TRPPROGLYr';'TRPPROLEUr';'TRPPROVALr';'TRPSERTYRr';'TRPTHRGLUr';'TRPTHRILEr';'TRPTHRTYRr';'TRPTYRGLNr';'TRPTYRTYRr';'TRPVALASPr';'TYRALAr';'TYRALAPHEr';'TYRARGGLUr';'TYRARGSERr';'TYRASPARGr';'TYRCYSGLYr';'TYRCYSTHRr';'TYRGLUr';'TYRLEUARGr';'TYRPHETYRr';'TYRTHRr';'TYRTRPPHEr';'TYRTYRr';'TYRVALMETr';'VALARGGLYr';'VALHISASNr';'VALLEUPHEr';'VALLYSTYRr';'VALPHEARGr';'VALPROTRPr';'VALSERARGr';'VALTRPPHEr';'VALTRPVALr';'VALVALr';'TRPGLYASPr'};
model = changeRxnBounds(model,unlIrrRxns,0,'l');

%% create mapping from genes onto gene loci
% download latest mouse assembly from NCBI
websave('GCF_000001635.27_GRCm39_feature_table.txt.gz','https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_feature_table.txt.gz')
gunzip('GCF_000001635.27_GRCm39_feature_table.txt.gz')

geneMapping = readInputTableForPipeline('GCF_000001635.27_GRCm39_feature_table.txt');

%% define media
retina_mets = {'2pg[e]';'34hpp[e]';'35cgmp[e]';'3aib[e]';'3pg[e]';'4abut[e]';'4hpro_LT[e]';'4pyrdx[e]';'5aop[e]';'5htrp[e]';'5oxpro[e]';'ac[e]';'acac[e]';'ach[e]';'ade[e]';'adn[e]';'adp[e]';'adpac[e]';'adrnl[e]';'akg[e]';'ala__L[e]';'allop[e]';'alltn[e]';'amp[e]';'anth[e]';'arg__L[e]';'asn__L[e]';'asp__L[e]';'atp[e]';'avite1[e]';'avite2[e]';'bgly[e]';'bhb[e]';'btn[e]';'but[e]';'C02470[e]';'C02528[e]';'ca2[e]';'chol[e]';'chsterol[e]';'cit[e]';'citr__L[e]';'cl[e]';'cmp[e]';'co[e]';'co2[e]';'creat[e]';'crn[e]';'crtn[e]';'csn[e]';'cu2[e]';'cyst__L[e]';'cytd[e]';'dcmp[e]';'dgchol[e]';'dhap[e]';'dmgly[e]';'dtmp[e]';'duri[e]';'fdp[e]';'fe2[e]';'fe3[e]';'fol[e]';'for[e]';'fru[e]';'fum[e]';'g1p[e]';'gal[e]';'gam[e]';'gchola[e]';'gdp[e]';'glc__D[e]';'glcur[e]';'gln__L[e]';'glu__L[e]';'glutar[e]';'gly[e]';'glyald[e]';'glyb[e]';'glyc__R[e]';'glyc[e]';'glyc3p[e]';'gmp[e]';'gsn[e]';'gthox[e]';'gthrd[e]';'gtp[e]';'gudac[e]';'h[e]';'h2co3[e]';'h2o[e]';'h2o2[e]';'HC00319[e]';'HC00900[e]';'hco3[e]';'hcys__L[e]';'hgentis[e]';'his__L[e]';'hista[e]';'hLkynr[e]';'hom__L[e]';'homoval[e]';'hpdca[e]';'hxan[e]';'i[e]';'ile__L[e]';'imp[e]';'inost[e]';'ins[e]';'k[e]';'kynate[e]';'L2aadp[e]';'lac__L[e]';'lcts[e]';'leu__L[e]';'Lkynr[e]';'lnlc[e]';'lnlnca[e]';'lys__L[e]';'M01966[e]';'mal__L[e]';'melatn[e]';'meoh[e]';'met__L[e]';'methsucc[e]';'mev__R[e]';'mg2[e]';'mhista[e]';'nac[e]';'ncam[e]';'nh4[e]';'no[e]';'no2[e]';'normete__L[e]';'o2[e]';'oaa[e]';'oh1[e]';'orn__D[e]';'orot[e]';'orot5p[e]';'oxa[e]';'oxyp[e]';'pchol_hs[e]';'pep[e]';'phe__L[e]';'phpyr[e]';'pi[e]';'pnto__R[e]';'ppa[e]';'ppi[e]';'pro__L[e]';'prostge2[e]';'prpp[e]';'pydam[e]';'pydx[e]';'pydx5p[e]';'pydxn[e]';'pyr[e]';'q10[e]';'q10h2[e]';'retinol[e]';'ribflv[e]';'sbt__D[e]';'sel[e]';'ser__L[e]';'so3[e]';'so4[e]';'succ[e]';'sucr[e]';'taur[e]';'thf[e]';'thm[e]';'thr__L[e]';'trp__L[e]';'tsul[e]';'ttdca[e]';'tyr__L[e]';'udp[e]';'ura[e]';'urate[e]';'uri[e]';'val__L[e]';'xan[e]';'xtsn[e]';'zn2[e]'};
% retina metabolites based on https://doi.org/10.1074/jbc.M115.698985

retina_medium = {'EX_2pg_e';'EX_34hpp_e';'EX_35cgmp_e';'EX_3aib_e';'EX_3pg_e';'EX_4abut_e';'EX_4hpro_LT_e';'EX_4pyrdx_e';'EX_5aop_e';'EX_5htrp_e';'EX_5oxpro_e';'EX_ac_e';'EX_acac_e';'EX_ach_e';'EX_ade_e';'EX_adn_e';'EX_adp_e';'EX_adpac_e';'EX_adrnl_e';'EX_akg_e';'EX_ala__L_e';'EX_allop_e';'EX_alltn_e';'EX_amp_e';'EX_anth_e';'EX_arg__L_e';'EX_asn__L_e';'EX_asp__L_e';'EX_atp_e';'EX_avite1_e';'EX_avite2_e';'EX_bgly_e';'EX_bhb_e';'EX_btn_e';'EX_but_e';'EX_C02470_e';'EX_C02528_e';'EX_ca2_e';'EX_chol_e';'EX_chsterol_e';'EX_cit_e';'EX_citr__L_e';'EX_cl_e';'EX_cmp_e';'EX_co_e';'EX_co2_e';'EX_creat_e';'EX_crn_e';'EX_crtn_e';'EX_csn_e';'EX_cu2_e';'EX_cyst__L_e';'EX_cytd_e';'EX_dcmp_e';'EX_dgchol_e';'EX_dhap_e';'EX_dmgly_e';'EX_dtmp_e';'EX_duri_e';'EX_fdp_e';'EX_fe2_e';'EX_fe3_e';'EX_fol_e';'EX_for_e';'EX_fru_e';'EX_fum_e';'EX_g1p_e';'EX_gal_e';'EX_gam_e';'EX_gchola_e';'EX_gdp_e';'EX_glc__D_e';'EX_glcur_e';'EX_gln__L_e';'EX_glu__L_e';'EX_glutar_e';'EX_gly_e';'EX_glyald_e';'EX_glyb_e';'EX_glyc__R_e';'EX_glyc_e';'EX_glyc3p_e';'EX_gmp_e';'EX_gsn_e';'EX_gthox_e';'EX_gthrd_e';'EX_gtp_e';'EX_gudac_e';'EX_h_e';'EX_h2co3_e';'EX_h2o_e';'EX_h2o2_e';'EX_HC00319_e';'EX_HC00900_e';'EX_hco3_e';'EX_hcys__L_e';'EX_hgentis_e';'EX_his__L_e';'EX_hista_e';'EX_hLkynr_e';'EX_hom__L_e';'EX_homoval_e';'EX_hpdca_e';'EX_hxan_e';'EX_i_e';'EX_ile__L_e';'EX_imp_e';'EX_inost_e';'EX_ins_e';'EX_k_e';'EX_kynate_e';'EX_L2aadp_e';'EX_lac__L_e';'EX_lcts_e';'EX_leu__L_e';'EX_Lkynr_e';'EX_lnlc_e';'EX_lnlnca_e';'EX_lys__L_e';'EX_M01966_e';'EX_mal__L_e';'EX_melatn_e';'EX_meoh_e';'EX_met__L_e';'EX_methsucc_e';'EX_mev__R_e';'EX_mg2_e';'EX_mhista_e';'EX_nac_e';'EX_ncam_e';'EX_nh4_e';'EX_no_e';'EX_no2_e';'EX_normete__L_e';'EX_o2_e';'EX_oaa_e';'EX_oh1_e';'EX_orn__D_e';'EX_orot_e';'EX_orot5p_e';'EX_oxa_e';'EX_oxyp_e';'EX_pchol_hs_e';'EX_pep_e';'EX_phe__L_e';'EX_phpyr_e';'EX_pi_e';'EX_pnto__R_e';'EX_ppa_e';'EX_ppi_e';'EX_pro__L_e';'EX_prostge2_e';'EX_prpp_e';'EX_pydam_e';'EX_pydx_e';'EX_pydx5p_e';'EX_pydxn_e';'EX_pyr_e';'EX_q10_e';'EX_q10h2_e';'EX_retinol_e';'EX_ribflv_e';'EX_sbt__D_e';'EX_sel_e';'EX_ser__L_e';'EX_so3_e';'EX_so4_e';'EX_succ_e';'EX_sucr_e';'EX_taur_e';'EX_thf_e';'EX_thm_e';'EX_thr__L_e';'EX_trp__L_e';'EX_tsul_e';'EX_ttdca_e';'EX_tyr__L_e';'EX_udp_e';'EX_ura_e';'EX_urate_e';'EX_uri_e';'EX_val__L_e';'EX_xan_e';'EX_xtsn_e';'EX_zn2_e'};

%% loop through the different data files

% read in and normalize the RNA sequencing data
normData = readInputTableForPipeline([pwd filesep 'Data' filesep 'Retina_RNAseq.csv']);

% map gene symbols to gene IDs
genesMapped = {'GeneSymbol','GeneID'};
cnt=2;
for i=2:size(normData,1)
    genesMapped{cnt,1} = normData{i,1};
    findGene = find(strcmp(geneMapping(:,15),normData{i,1}));
    if ~isempty(findGene)
        genesMapped{cnt,2}=geneMapping{findGene(1),16};
        cnt=cnt+1;
    end
end

for i=2:size(genesMapped,1)
    genesMapped{i,2} = num2str(genesMapped{i,2});
end

[C,I] = setdiff(genesMapped(:,2),model.genes,'stable');
genesMapped(I(2:end),:) = [];

% proceed with the normalized and translated datag
data = cell2mat(normData(2:end,2:end));
discretized = discretize_FPKM(data, normData(1,2:end),1);

biomass_rxn = {'BIOMASS_reaction'}; %find the corresponding biomass reaction in your model, here it is for recon 2

epsilon = 1e-4; %avoid small number errors

already_mapped_tag = 0;

consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples

unpenalizedSystems = {'Transport, endoplasmic reticular'
    'Transport, extracellular'
    'Transport, golgi apparatus'
    'Transport, mitochondrial'
    'Transport, peroxisomal'
    'Transport, lysosomal'
    'Transport, nuclear'};

subs = {};
for i=1:length(model.subSystems)
    try
        subs{i,1} = model.subSystems{i,1}{1};
    catch
        subs{i,1} = model.subSystems{i,1};
    end
end
unpenalized = model.rxns(ismember(subs,unpenalizedSystems));

optional_settings.unpenalized = unpenalized;
% make sure one-carbon metabolism is included, some other universal reactions
optional_settings.func = {'BIOMASS_reaction'}; % forced additional reactions into the  model

% make sure diet exchanges are kept
optional_settings.medium = retina_mets;
optional_settings.func = union(optional_settings.func,retina_medium);

growth = {};

%% Create models
mkdir([pwd filesep 'Models'])

samples = normData(1,2:end);
for i = 1:length(samples) %for each sample[modelS, A_keep] = fastcormics_RNAseq(model, discretized(:,i), normData(2:end,1), genesMapped, ...
        biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
    modelS = updateGenes(modelS);
    
    
    % scale fluxes with eFlux
    expressionData = struct;
    expressionData.gene = normData(2:end,1);
    for j=1:length(expressionData.gene)
        findGene = find(strcmp(genesMapped(:,1),expressionData.gene{j}));
        if ~isempty(findGene)
            expressionData.gene{j} = num2str(genesMapped{findGene,2});
        end
    end
    expressionData.value = cell2mat(normData(2:end,i+1));
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(modelS, expressionData);
    expression=struct;
    expression.target=modelS.rxns;
    expression.value=expressionRxns;
    expression.value(find(isnan(expression.value)))=-1;
    % exclude reactions for which constraints should not be changed
    unchRxns = {'BIOMASS_reaction','DM_atp_c_'};
    [C,I] = intersect(expression.target(:,1),unchRxns);
    expression.value(I,1)=-1;
    expression.preprocessed=true;
    % modified version of original e-Flux method
    modelS = relaxedApplyEFluxConstraints(modelS, expression);
    
    modelS = rmfield(modelS,'subSystems');
    modelS.subSystems = cell(length(modelS.rxns),1);
    
    modelS=changeObjective(modelS,'BIOMASS_reaction');
    
    % test growth
    modelTest = changeRxnBounds(modelS,retina_medium,-1,'l');

    % some metabolites need higher uptake rates
    modelTest = changeRxnBounds(modelTest,{'EX_glc__D_e'},-10,'l');
    modelTest = changeRxnBounds(modelTest,{'EX_o2_e';'EX_h2o_e';'EX_co2_e';'EX_hco3_e'},-100,'l');
    % some metabolites need to be included
    modelTest = changeRxnBounds(modelTest,{'EX_idl_hs_e','EX_ldl_hs_e','EX_hdl_hs_e'},-1,'l');
    
    FBA=optimizeCbModel(modelTest,'max');
    growth{i,1} = samples{i};
    growth{i,2} = FBA.f;
    % save model
    writeCbModel(modelS,'format','mat','fileName',[pwd filesep 'Models' filesep samples{i} '.mat'])
end
save([pwd filesep 'Growth_retina_medium.mat'],'growth')
