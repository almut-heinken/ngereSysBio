

% interrogate the RNA-seq constrained models:
% first get get the model statistics and distribution of reactions,
% metabolites, and genes
% then sample the flux space to get the allowed flux through all reactions
% in the sample-specific models

numWorkers=4;
if ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

oriModel=readCbModel('iMM1865.xml');

minFluxes = {};
maxFluxes = {};

retina_medium = {'EX_2pg_e';'EX_34hpp_e';'EX_35cgmp_e';'EX_3aib_e';'EX_3pg_e';'EX_4abut_e';'EX_4hpro_LT_e';'EX_4pyrdx_e';'EX_5aop_e';'EX_5htrp_e';'EX_5oxpro_e';'EX_ac_e';'EX_acac_e';'EX_ach_e';'EX_ade_e';'EX_adn_e';'EX_adp_e';'EX_adpac_e';'EX_adrnl_e';'EX_akg_e';'EX_ala__L_e';'EX_allop_e';'EX_alltn_e';'EX_amp_e';'EX_anth_e';'EX_arg__L_e';'EX_asn__L_e';'EX_asp__L_e';'EX_atp_e';'EX_avite1_e';'EX_avite2_e';'EX_bgly_e';'EX_bhb_e';'EX_btn_e';'EX_but_e';'EX_C02470_e';'EX_C02528_e';'EX_ca2_e';'EX_chol_e';'EX_chsterol_e';'EX_cit_e';'EX_citr__L_e';'EX_cl_e';'EX_cmp_e';'EX_co_e';'EX_co2_e';'EX_creat_e';'EX_crn_e';'EX_crtn_e';'EX_csn_e';'EX_cu2_e';'EX_cyst__L_e';'EX_cytd_e';'EX_dcmp_e';'EX_dgchol_e';'EX_dhap_e';'EX_dmgly_e';'EX_dtmp_e';'EX_duri_e';'EX_fdp_e';'EX_fe2_e';'EX_fe3_e';'EX_fol_e';'EX_for_e';'EX_fru_e';'EX_fum_e';'EX_g1p_e';'EX_gal_e';'EX_gam_e';'EX_gchola_e';'EX_gdp_e';'EX_glc__D_e';'EX_glcur_e';'EX_gln__L_e';'EX_glu__L_e';'EX_glutar_e';'EX_gly_e';'EX_glyald_e';'EX_glyb_e';'EX_glyc__R_e';'EX_glyc_e';'EX_glyc3p_e';'EX_gmp_e';'EX_gsn_e';'EX_gthox_e';'EX_gthrd_e';'EX_gtp_e';'EX_gudac_e';'EX_h_e';'EX_h2co3_e';'EX_h2o_e';'EX_h2o2_e';'EX_HC00319_e';'EX_HC00900_e';'EX_hco3_e';'EX_hcys__L_e';'EX_hgentis_e';'EX_his__L_e';'EX_hista_e';'EX_hLkynr_e';'EX_hom__L_e';'EX_homoval_e';'EX_hpdca_e';'EX_hxan_e';'EX_i_e';'EX_ile__L_e';'EX_imp_e';'EX_inost_e';'EX_ins_e';'EX_k_e';'EX_kynate_e';'EX_L2aadp_e';'EX_lac__L_e';'EX_lcts_e';'EX_leu__L_e';'EX_Lkynr_e';'EX_lnlc_e';'EX_lnlnca_e';'EX_lys__L_e';'EX_M01966_e';'EX_mal__L_e';'EX_melatn_e';'EX_meoh_e';'EX_met__L_e';'EX_methsucc_e';'EX_mev__R_e';'EX_mg2_e';'EX_mhista_e';'EX_nac_e';'EX_ncam_e';'EX_nh4_e';'EX_no_e';'EX_no2_e';'EX_normete__L_e';'EX_o2_e';'EX_oaa_e';'EX_oh1_e';'EX_orn__D_e';'EX_orot_e';'EX_orot5p_e';'EX_oxa_e';'EX_oxyp_e';'EX_pchol_hs_e';'EX_pep_e';'EX_phe__L_e';'EX_phpyr_e';'EX_pi_e';'EX_pnto__R_e';'EX_ppa_e';'EX_ppi_e';'EX_pro__L_e';'EX_prostge2_e';'EX_prpp_e';'EX_pydam_e';'EX_pydx_e';'EX_pydx5p_e';'EX_pydxn_e';'EX_pyr_e';'EX_q10_e';'EX_q10h2_e';'EX_retinol_e';'EX_ribflv_e';'EX_sbt__D_e';'EX_sel_e';'EX_ser__L_e';'EX_so3_e';'EX_so4_e';'EX_succ_e';'EX_sucr_e';'EX_taur_e';'EX_thf_e';'EX_thm_e';'EX_thr__L_e';'EX_trp__L_e';'EX_tsul_e';'EX_ttdca_e';'EX_tyr__L_e';'EX_udp_e';'EX_ura_e';'EX_urate_e';'EX_uri_e';'EX_val__L_e';'EX_xan_e';'EX_xtsn_e';'EX_zn2_e'};

% Starting simulations

mkdir([pwd filesep 'Metabolic_Flux_Results'])

modelFolder = [pwd filesep 'Models'];
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList(find(strncmp(modelList(:,1),'._',2)),:)=[];

for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    
    % add subsystems
    for j=1:length(model.rxns)
        rxnID = find(strcmp(oriModel.rxns,model.rxns{j}));
        model.subSystems{j} = oriModel.subSystems{rxnID};
    end
    
    % constrain unlikely drug metabolism reactions
    drugRxns = {'2HATVACIDhc','2HATVACIDteb','2HATVACIDthc','2HATVLAChc','2HATVLACteb','2HATVLACthc','OXYP1CONJ','OXYPR1tehv','OXYPtepv','OXYPthc','ALLOPOXDhep','ALLOPtepvb','CRVSATPtu','CRVSM1hr','CRVSM23hr','CRVSthc','CRVStu','CVM1GLUChc','CVM23GLUChc','LVACLAChep','LVSTACIDhep','MERACMPtep','MERACMPthc','PTVSTATPtu','PTVSTGLUChc','PTVSTLAChc','PTVSThc','PTVSTtu','PVSATPtu','PVSHtu','PVSOATPtu','RSVATPtu','RSVtu','SMVACIDhep','SMVGLUChep','SMVHYDROhep','SMVLAChep','ATVACIDhr','ATVETHGLUChc','ATVLACGLCURhc','ATVLAChr','ACMPGLUTtep','ACMPGLUTthc','ACMPGLUtep','ACMPGLUthc','6EPSteb','6EPVSthc','3HPVSteb','3HPVSthc','3ISPVSteb','3ISPVSthc','6EPSteb','6EPVSthc','4HATVACIDhc','4HATVACIDteb','4HATVACIDthc','4HATVLAChc','4HATVLACteb','4HATVLACthc','ATVACIDOATPtu','ATVACIDtdu','OXYPtepv','OXYPthc','2HATVACIDtep','2HATVLACtep','35DHPVStep','35DHPVSthc','3HPVSTETtev','3HPVStep','3ISPVStep','4HATVACIDtep','4HATVLACtep','56DHPVStev','56EPPVStev','6EPVStep','7HPVStev','ACMPdt','ALLOPtepvb','ATVACIDMCTtu','ATVACIDOATPtu','ATVLACThc','ATVLACitr','ATVLACtdhc','CRVS1tev','CRVSM23tev','CRVSM24tev','DSPVStev','FVStep','FVStu','OXYPtepv','OXYPthc','PTVSTtep','PVSGLUCtev','PVStep','RSVLACitr','RSVtev','SMVthep','SMVtv','TRIPVStev'};
    model = changeRxnBounds(model,drugRxns,0,'b');
    
    % implement medium
    model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),0,'l');
    model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'SK_',3))),0,'b');
    model = changeRxnBounds(model,retina_medium,-1,'l');

    % some metabolites need higher uptake rates
    model = changeRxnBounds(model,{'EX_glc__D_e'},-10,'l');
    model = changeRxnBounds(model,{'EX_o2_e';'EX_h2o_e';'EX_co2_e';'EX_hco3_e'},-100,'l');
    % some metabolites need to be included
    model = changeRxnBounds(model,{'EX_idl_hs_e','EX_ldl_hs_e','EX_hdl_hs_e'},-1,'l');
    % enforce biomass production
    model = changeRxnBounds(model,'BIOMASS_reaction',0.01,'l');
    
    % compute minimal and maximal fluxes
    [minFlux,maxFlux] = fastFVA(model,99,'max','ibm_cplex',model.rxns, 'S');
    minFluxes{i} = minFlux;
    maxFluxes{i} = maxFlux;
end
save([pwd filesep 'Metabolic_Flux_Results' filesep 'MinFluxes'],'minFluxes')
save([pwd filesep 'Metabolic_Flux_Results' filesep 'MaxFluxes'],'maxFluxes')

% get reconstruction statistics:
% reactions, metabolites, and genes
stats={};
stats{1,1}='Model_ID';
stats{1,2}='Reactions';
stats{1,3}='Metabolites';
stats{1,4}='Genes';

for i=1:length(modelList)
    stats{i+1,1} = strrep(modelList{i},'.mat','');
    model = readCbModel([modelFolder filesep modelList{i}]);
    
    % Number of reactions, metabolites, and genes
    stats{i+1,2}=length(model.rxns);
    stats{i+1,3}=length(model.mets);
    stats{i+1,4}=length(model.genes);
end
save([pwd filesep 'Metabolic_Flux_Results' filesep 'Model_statistics'],'stats')
