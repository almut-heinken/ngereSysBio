
% integrate blood metabolome available for 34 patients

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
solverOK=changeCobraSolver('ibm_cplex','QP');

% load the human reconstruction
modelHuman = readCbModel('input/Recon3D_CBL.mat');
modelHuman = changeRxnBounds(modelHuman,'biomass_reaction',1,'b');

% rename duplicate metabolite 2-hydroxyglutarate
modelHuman = removeRxns(modelHuman,{'HMR_0718','HMR_0719'});
modelHuman = addReaction(modelHuman,'HMR_0718','reactionFormula','akg[c] + fadh2[c] -> fad[c] + 2hydog[c]',...
    'subSystem','Butanoate metabolism','geneRule','79944.1');
modelHuman = addReaction(modelHuman,'HMR_0719','reactionFormula','akg[c] + h[c] + nadh[c] -> 2hydog[c] + nad[c]',...
    'subSystem','Butanoate metabolism','geneRule','26227.1');

% adapt some subsystem names
edit_subs = {'R_group_phosphotase_1','Fatty acid oxidation';'HMR_1403','Fatty acid oxidation';'r0390','Starch and sucrose metabolism';'r0407','Starch and sucrose metabolism';'r0408','Starch and sucrose metabolism';'r0409','Starch and sucrose metabolism';'RE2848C','Fatty acid oxidation';'r0604','Valine, leucine, and isoleucine metabolism';'r1167','Fatty acid oxidation';'r1177','Fatty acid oxidation';'r1418','Fatty acid oxidation';'r1466','Fatty acid oxidation';'CATp','ROS detoxification';'RE2526C','Tyrosine metabolism';'r0683','Fatty acid oxidation';'r0627','Glycerophospholipid metabolism'};
for i=1:length(edit_subs)
    findRxn = find(strcmp(modelHuman.rxns,edit_subs{i,1}));
    modelHuman.subSystems{findRxn,1} = edit_subs{i,2};
end

% prepare the reconstruction
% silence drug eactions-should not have flux
drugRxns = {'12HTACRhr';'12HTACRtu';'12HTACRtep';'1331TAALThr';'1331TACRhr';'1331TACRteb';'1331TACRtev';'13DMThr';'13DMTtu';'13DMTtep';'14HMDZALThr';'14HMDZhr';'14MDZtev';'1513DTALThr';'1513DTACRhr';'1513TACRtu';'1513TACRtep';'1531TACRhr';'1531TACRteb';'1531TACRtev';'1531TALThr';'15DMThr';'15DMTtu';'15DMTtep';'1HIBUPGLUC_Sthv';'1HIBUP_SGLUhep';'1HIBUP_Sthv';'1HMDGLUChr';'1HMDZGLUChc';'1OHMDZhr';'1OHMDZtep';'2HATVLAChc';'2HATVACIDGLUChr';'2HATVACIDGLUCteb';'2HATVACIDOXDhc';'2HATVACIDhc';'2HATVACIDteb';'2HATVACIDtep';'2HATVACIDthc';'2HATVLACGLUChr';'2HATVLACGLUCteb';'2HATVLACOXDhc';'2HATVLACteb';'2HATVLACtep';'2HATVLACthc';'2HIBUPGLUC_Sthv';'2HIBUP_Rthv';'2HIBUP_SGLUhep';'2HIBUP_Sthv';'31DMThr';'31DMTtu';'31DMTtep';'35DHPVShc';'35DHPVStep';'35DHPVSthc';'35DSMVhep';'35DSMVteb';'3HIBUPGLUC_Sthv';'3HIBUP_Rthv';'3HIBUP_SGLUhep';'3HIBUP_Sthv';'3HLVSTAChep';'3HLVSTACtbc';'3HPVSTETCOAhcm';'3HPVSTETCOAhcx';'3HPVSTETteb';'3HPVSTETtev';'3HPVShc';'3HPVSteb';'3HPVStep';'3HPVSthc';'3HSMVACIDhep';'3HSMVACIDteb';'3HSMVhep';'3ISPVShc';'3ISPVSteb';'3ISPVStep';'3ISPVSthc';'3MEACMPhc';'3OHACMPhr';'3OHACMPtev';'4BHGLZABCt';'4BHGLZhr';'4BHGLZtev';'4HATVACIDOXDhc';'4HATVACIDhc';'4HATVACIDteb';'4HATVACIDtep';'4HATVACIDthc';'4HATVLACOXDhc';'4HATVLAChc';'4HATVLACteb';'4HATVLACtep';'4HATVLACthc';'4HMDGLUCtev';'4HMDZGLUChr';'4OHMDZhr';'4OHMDZtev';'56DHPVShc';'56DHPVSteb';'56DHPVStev';'56EPPVShc';'56EPPVSteb';'56EPPVStev';'5OHFVSGLUhc';'5OHFVSGLUtev';'5OHFVShc';'5OHFVSteb';'6AHGLZABCt';'6AHGLZhr';'6AHGLZtev';'6BHGLZABCt';'6BHGLZGLCABCt';'6BHGLZGLChr';'6BHGLZGLCtev';'6BHGLZhr';'6BHGLZtev';'6CSMVACIDhep';'6CSMVACIDteb';'6CSMVhep';'6EPSteb';'6EPVStep';'6EPVShc';'6EPVSthc';'6HLVSTAChep';'6HLVSTthep';'6HMSMVACIDhep';'6HMSMVACIDteb';'6HMSMVhep';'6HSMVACIDhep';'6HSMVACIDteb';'6HSMVhep';'6MELACAChep';'6MELVACtbc';'6MELVSTthep';'6MSMVhep';'6OHFVSGLUhc';'6OHFVSGLUtev';'6OHFVShc';'6OHFVSteb';'7AHGLZABCt';'7AHGLZhr';'7AHGLZtev';'7BHGLZABCt';'7BHGLZGLCABCt';'7BHGLZGLChr';'7BHGLZGLCtev';'7BHGLZhr';'7BHGLZtev';'7HPVShc';'7HPVSteb';'7HPVStev';'ALLOP2tu';'ACMPtu';'ACMPGLUChr';'ACMPGLUTdt';'ACMPGLUTtep';'ACMPGLUTthc';'ACMPGLUtep';'ACMPGLUthc';'ACMPShc';'ACMPdt';'ACMPthc';'ALLOP1tu';'ALLOPOXDhep';'ALLOPtepvb';'AM19CSALThr';'AM19CShr';'AM19CSteb';'AM1A4NCShc';'AM1A4NCSteb';'AM1ACCShr';'AM1ACCStev';'AM1ACShr';'AM1ACSteb';'AM1ACStep';'AM1ALCShr';'AM1ALCSteb';'AM1ALCStep';'AM1C4N9CShc';'AM1C4N9CSteb';'AM1C9CShr';'AM1C9CSteb';'AM1C9CStev';'AM1CCShr';'AM1CCSteb';'AM1CCStev';'AM1CGLChr';'AM1CGLCteb';'AM1CSAhr';'AM1CSAtep';'AM4N9CShc';'AM4N9CShr';'AM4N9CStev';'AM4NC9CSteb';'AM4NCShr';'AM4NCSteb';'AM4NCStep';'AM9CSAhr';'AM9CSAteb';'AM9CSAtep';'ATVACIDMCTtu';'ATVACIDOATPtu';'ATVACIDhc';'ATVACIDhr';'ATVACIDtdu';'ATVACIDtu';'ATVACYLGLUChc';'ATVETHGLUChc';'ATVLACGLCURhc';'ATVLACThc';'ATVLACh2r';'ATVLAChc';'ATVLAChr';'ATVLACtdhc';'ATVLACtu';'Am19CStev';'Am1CSAteb';'CARBIBUP_SGLUthv';'CARIBUP_Rthv';'CARIBUP_SGLUhep';'CARIBUP_Sthv';'CRGLZABCt';'CRGLZhr';'CRGLZtev';'CRVS1M24hc';'CRVS1tev';'CRVS23M24hc';'CRVSATPthc';'CRVSATPtu';'CRVSM1SPhc';'CRVSM1hc';'CRVSM1hr';'CRVSM1teb';'CRVSM22hc';'CRVSM23hc';'CRVSM23hr';'CRVSM23teb';'CRVSM23tev';'CRVSM24teb';'CRVSM24tev';'CRVSM31hc';'CRVStu';'CRVSthc';'CSASULPhc';'CSASULPteb';'CSASULPtev';'CSAtd';'CSAtu';'CVM1GLUChc';'CVM23GLUChc';'CYSACMPAChc';'CYSAMPtev';'DELACCRVSM23hc';'DEOXFVShc';'DEOXFVStev';'DESFVShc';'DESFVSteb';'DHGLZABCt';'DHGLZhc';'DHGLZtev';'DSPVShc';'DSPVSteb';'DSPVStev';'EPOXTAChr';'EPOXTACteb';'EPOXTACtev';'FVSGLUChc';'FVSTETGLUhc';'FVSTETGLUtev';'FVSTETtev';'FVShc';'FVSteb';'FVStep';'FVStu';'GLC3MEACPhr';'GLC3MEACPtev';'GLZABCteb';'GLZtd';'GTACMPhr';'GTACMPtev';'IBUPGLUCtchep';'IBUPGLUCtpvb';'IBUPGT_HEP';'IBUP_RASCL1hep';'IBUP_RCYP2hep';'IBUP_RCYP3hep';'IBUP_RCYPCARhep';'IBUP_Rshep';'IBUP_Rtdhep';'IBUP_Rtdu';'IBUP_SACOT2';'IBUP_SCONJhep';'IBUP_SCYP1hep';'IBUP_SCYP2hep';'IBUP_SCYP3hep';'IBUP_SCYPCARhep';'IBUP_Stbc';'IBUP_Stdhep';'IBUP_Stdu';'ISOLVSTAChep';'ISOLVSTtbc';'LST4EXPTDhc';'LST4EXPhr';'LST4EXPthc';'LSTN1GLUChr';'LSTN1GLUCtev';'LSTNtu';'LSTNM1hr';'LSTNM1tev';'LSTNM2hr';'LSTNM2tev';'LSTNM4hr';'LSTNM4tev';'LSTNM5hr';'LSTNM5tev';'LSTNM7TDhc';'LSTNM7hr';'LSTNM7thc';'LSTNRATt';'LSTNtd';'LVACLAChep';'LVSTACIDhep';'LVSTACIDtu';'LVSTACOXD6Hhep';'LVSTACOXD6MEhep';'LVSTOXD3Hhep';'LVSTOXD6Hhep';'LVSTOXD6METhep';'LVSTPGPtu';'LVSTtu';'MDZGLCtev';'MDZtd';'MDZtu';'MERACMPtep';'MERACMPthc';'MHGLZABCt';'MHGLZhr';'MHGLZtev';'NDERSVhc';'NDERSVteb';'NFDACOXDhc';'NFDACtep';'NFDDMEThr';'NFDLAChc';'NFDLACtep';'NFDNPYtep';'NFDOHtep';'NFDOXDhc';'NFDtd';'OXYP1CONJ';'OXYP2CONJ';'OXYPR1tehv';'OXYPR7tehv';'OXYPthc';'OXYPtepv';'PROFVSCOAhc';'PROFVShc';'PROFVStev';'PTVSTATPtu';'PTVSTGLUChc';'PTVSTLAChc';'PTVSTLACtev';'PTVSTM13hr';'PTVSTM3eb';'PTVSTM3hc';'PTVSThc';'PTVSTtep';'PTVSTtu';'PVSATPtu';'PVSGLUChc';'PVSGLUCteb';'PVSGLUCtev';'PVSHtu';'PVSOATPtu';'PVStep';'RSVATPtu';'RSVGLUChc';'RSVLAChv';'RSVLACteb';'RSVSPONhc';'RSVhc';'RSVtev';'RSVtu';'S3MEACMPhc';'S3MEACMPtev';'SMVACIDATPteb';'SMVACIDhep';'SMVACIDtev';'SMVtu';'SMVGLUCLAChep';'SMVGLUChep';'SMVHYDROhep';'SMVLAChep';'SMVtv';'SMVthep';'STACMPhc';'STACMPtev';'SULPACMPtev';'TACRDtsc';'TACRtu';'TAURIBUP_Sthv';'THRFVShc';'THRFVStev';'THSACMPhr';'TLACFVShc';'TLACFVStev';'TMDM1OATt';'TMDM1hr';'TMDM3OATt';'TMDM3hr';'TMDM5OATt';'TMDM5hr';'TMDOATPtsc';'TMDOATtev';'TMDOATthc';'TMDtd';'TRIPVShc';'TRIPVSteb';'TRIPVStev';'TSACGLUCtev';'TSACMGLUChr';'TSACMSULhc';'TSACMSULtev';'12HTACRitr';'13HTACRitr';'14HMDZitr';'1513TACRitr';'1531TACRitr';'1HIBUP_Sitr';'1HIBUPGLUitr';'1HMDGLUCitr';'2HATVACIDGLUCitr';'2HATVLACGLUCitr';'2HIBUP_Ritr';'2HIBUP_Sitr';'2HIBUPGLUC_Sitr';'35DHPVSitr';'35DSMVitr';'3HIBUP_Ritr';'3HIBUPGLUC_Sitr';'3HLVSTitr';'3HPVSitr';'3HPVSCOAitm';'3HPVSCOAitx';'3HPVSTETCOAitm';'3HPVSTETCOAitx';'3HSMVitr';'3ISPVSitr';'3MEACMPitr';'3OHACMPitr';'4BHGLZitr';'4HATVACIDitr';'4HATVLACitr';'4HMDGLUCitr';'4OHMDZitr';'56DHPVSitr';'56EPPVSitr';'5OHFVSitr';'5OHFVSGLUitr';'6AHGLZitr';'6BHGLZGLCitr';'6CSMVitr';'6HLVSTitr';'6HLVSTACIDitr';'6HMSMVitr';'6HSMVitr';'6MELVACIDitr';'6MELVSTitr';'6OHFVSitr';'6OHFVSGLUitr';'7AHGLZitr';'7BHGLZGLCitr';'7HPVSitr';'ACMPitr';'ACMPGLUitr';'AM19CSitr';'AM1ACCSitr';'AM1ACSitr';'AM1ALCSitr';'AM1C9CSitr';'AM1CGLCitr';'AM4N9CSitr';'AM4NCSitr';'CARIBUP_Sitr';'CARIBUPGLU_Sitr';'CRGLZitr';'CRVSitr';'CRVSM22itr';'CRVSM24itr';'CRVSM31itr';'CSAitr';'DEOXFVSitx';'DESFVSitr';'DSPVSitr';'EPOXTACitr';'FVSitx';'FVSCOAitx';'FVSTETitr';'FVSTETGLUitr';'GLC3MEACPitr';'GLZitr';'GTACMPitr';'IBUP_Ritr';'IBUP_Sitr';'IBUPGLUCitr';'LSTN1GLUCitr';'LSTNitr';'LSTNM1itr';'LSTNM2itr';'LSTNM4itr';'LSTNM5itr';'LSTNM7itr';'LVSTitr';'LVSTACIDitr';'MDZitr';'MDZGLCitr';'NDERSVitr';'NFDNPYitr';'NFDOHitr';'PROFVSCOAitx';'PTVSTLACitr';'PTVSTM13itr';'PTVSTM13te';'PTVSTM3itr';'PVSitr';'PVSGLUCitr';'RSVLACitr';'TACRitr';'THSACMPitr';'TLACFVSitr';'TMACMPitr';'TMDitr';'TMDM1itr';'TMDM3itr';'TMDM5itr';'TSACMGLUCitr';'3HPVSCOAhc';'3HPVSTEThc';'ACMPGLUTTRsc';'FVSCOAhc';'MDZGLChr';'TMACMPhr';'1OHMDZitr';'CYSACMPitr';'NFDACitr';'ACMPGLUTitr';'NAPQIhr';'H2O2itr';'UDPRIBc';'GLYitr';'PAPSitr';'PAPitr';'13DMTitr';'15DMTitr';'ATVACIDitr';'ATVLACitr';'31DMTitr';'SMVitr';'6BHGLZitr';'7BHGLZitr';'AM1CCSitr';'LST4EXPitr';'MHGLZitr';'RSVitr';'TRIPVSitr';'SMVACIDitr';'PTVSTitr';'3HIBUP_Sitr';'FVSitr';'AM1CSAitr';'AM9CSAitr';'2HATVACIDitr'};
modelHuman = changeRxnBounds(modelHuman,drugRxns,0,'b');

% close sink reactions
modelHuman = changeRxnBounds(modelHuman,modelHuman.rxns(find(strncmp(modelHuman.rxns,'sink_',5))),0,'b');
modelHuman=changeRxnBounds(modelHuman,modelHuman.rxns(find(strncmp(modelHuman.rxns,'DM_',3))),0,'l');

% implement a generic diet
modelHuman=changeRxnBounds(modelHuman,modelHuman.rxns(find(strncmp(modelHuman.rxns,'EX_',3))),0,'l');
diet = {'EX_fru[e]',10;'EX_glc_D[e]',10;'EX_lcts[e]',10;'EX_malt[e]',10;'EX_sucr[e]',10;'EX_strch1[e]',10;'EX_ttdca[e]',10;'EX_hdca[e]',10;'EX_ocdca[e]',10;'EX_arach[e]',10;'EX_octa[e]',10;'EX_glyc[e]',10;'EX_chsterol[e]',10;'EX_hdcea[e]',10;'EX_ocdcea[e]',10;'EX_arachd[e]',10;'EX_lnlc[e]',10;'EX_lnlnca[e]',10;'EX_ala_L[e]',10;'EX_arg_L[e]',10;'EX_glu_L[e]',10;'EX_gly[e]',10;'EX_ile_L[e]',10;'EX_leu_L[e]',10;'EX_lys_L[e]',10;'EX_met_L[e]',10;'EX_pro_L[e]',10;'EX_asn_L[e]',10;'EX_asp_L[e]',10;'EX_his_L[e]',10;'EX_phe_L[e]',10;'EX_ser_L[e]',10;'EX_thr_L[e]',10;'EX_trp_L[e]',10;'EX_tyr_L[e]',10;'EX_gln_L[e]',10;'EX_val_L[e]',10;'EX_cys_L[e]',10;'EX_ca2[e]',1;'EX_C06453[e]',1;'EX_ccbl[e]',1;'EX_oxocbl[e]',1;'EX_aqcobale]',1;'EX_hxan[e]',1;'EX_thymd[e]',1;'EX_avite1[e]',1;'EX_yvite1[e]',1;'EX_btn[e]',1;'EX_fol[e]',1;'EX_thf[e]',1;'EX_5mthf[e]',1;'EX_ncam[e]',1;'EX_lipoate[e]',1;'EX_pydxn[e]',1;'EX_ribflv[e]',1;'EX_thm[e]',1;'EX_pnto_R[e]',1;'EX_pheme[e]',1;'EX_retn[e]',1;'EX_inost[e]',1;'EX_chol[e]',1;'EX_fe2[e]',1;'EX_fe3[e]',1;'EX_h2s[e]',1;'EX_k[e]',1;'EX_mg2[e]',1;'EX_mn2[e]',1;'EX_na1[e]',1;'EX_zn2[e]',1;'EX_so4[e]',1;'EX_pi[e]',10;'EX_o2[e]',100;'EX_co2[e]',100;'EX_hco3[e]',100;'EX_h2o[e]',1000;'EX_adp[e]',2;'EX_creat[e]',1;'EX_11_cis_retfa[e]',1;'EX_11_M01570[e]',1};

for i=1:size(diet,1)
    modelHuman = changeRxnBounds(modelHuman,diet{i,1},-(diet{i,2}),'l');
end

%% read one-carbon and TCA cycle metabolome in µg/mL
metabolome = readInputTableForPipeline(['input' filesep 'Metabolome_1CM.csv']);

% translate names of metabolites to corresponding model IDs
mapping = {'Lactate','lac_L';'Choline','chol';'Glycine','gly';'Cystathionine','cyst_L';'Serine','ser_L';'Betaine','glyb';'Cysteine','cys_L';'S-adenosylmethionine','amet';'Homocysteine','hcys_L';'Methionine','met_L';'Glutathione','gthrd';'S-adenosylhomocysteine','ahcys';'Methyltetrahydrofolate','5mthf';'Dimethylglycine','dmgly';'Glutamate','glu_L';'Malate','mal_L';'2-Hydroxyglutarate','2hydog';'Succinate','succ';'Propionate','ppa';'Isocitrate','icit';'Citrate','cit';'Fumarate','fum';'Oxaloacetate','oaa';'2-Ketoglutarate','akg';'Pyruvate','pyr'};
for i=2:size(metabolome,1)
    metabolome{i,1} = mapping{find(strcmp(mapping(:,1),metabolome{i,1})),2};
end

% implement constraints and save fluxes
mkdir([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes_1CM'])

for i=2:size(metabolome,2)
    i
    model = modelHuman;
    % formulate a pseudo-reaction for metabolomics data
    form = '';
    for j=2:size(metabolome,1)-1
        form = [form num2str(metabolome{j,i}) ' ' metabolome{j,1} '[e] + '];
    end
    form = [form num2str(metabolome{end,i}) ' ' metabolome{end,1} '[e] -> metabolome[e]'];
    model = addReaction(model,'metabolomeReaction',form);
    model = addDemandReaction(model,'metabolome[e]');
    model = changeObjective(model,'metabolomeReaction');

    FBA=optimizeCbModel(model,'max','1e-6')
    save([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes_1CM' filesep 'Fluxes_' metabolome{1,i} '.mat'],'FBA')
end

% export the fluxes
fluxTable = vertcat('Reaction',modelHuman.rxns,'metabolomeReaction','DM_metabolome[e]');


for i=2:size(metabolome,2)
    fluxTable{1,i} = metabolome{1,i};
    load([pwd filesep 'MetabolomeConstrainedModels' filesep 'Fluxes_1CM' filesep 'Fluxes_' metabolome{1,i} '.mat'])
    fluxTable(2:end,i) = num2cell(FBA.x);
end

writetable(cell2table(fluxTable),[pwd filesep 'MetabolomeResults' filesep 'Metabolome_1CM_ConstrainedFluxes.csv'],'writeVariableNames',false)

%% export annotations
annotations = {'Reaction','ReactionName','Subsystem'};
annotations(2:length(modelHuman.rxns)+1,1) = modelHuman.rxns;
annotations(2:length(modelHuman.rxns)+1,2) = modelHuman.rxnNames;
annotations(2:length(modelHuman.rxns)+1,3) = modelHuman.subSystems;
writetable(cell2table(annotations),[pwd filesep 'MetabolomeResults' filesep 'HumanReactionAnnotations.csv'],'writeVariableNames',false)

