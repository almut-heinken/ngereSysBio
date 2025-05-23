
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% gapfill fibroblast reconstruction and Recon3D with cobalamin metabolism
% pathways

vmhDB = readInputTableForPipeline(['input' filesep 'recon-store-reactions-1.tsv']);

% make unlikely reversible reactions irreversible
unlIrrRxns = {'ALAARGCYSr';'ALAARGGLYr';'ALAASNLEUr';'ALAGLYLYSr';'ALAHISALAr';'ALALYSTHRr';'ARGALAALAr';'ARGALAPHEr';'ARGALATHRr';'ARGARGr';'ARGARGLYSr';'ARGARGMETr';'ARGCYSGLYr';'ARGCYSSERr';'ARGGLUGLUr';'ARGGLUPROr';'ARGGLYGLYr';'ARGHISTHRr';'ARGLEUPHEr';'ARGLYSASPr';'ARGPHEARGr';'ARGPROMETr';'ARGPROTHRr';'ARGSERSERr';'ARGTYRVALr';'ARGVALCYSr';'ARGVALTRPr';'ASNASNARGr';'ASNCYSCYSr';'ASNMETPROr';'ASNPHEASPr';'ASNPHECYSr';'ASNTYRGLYr';'ASNTYRPHEr';'ASNTYRTHRr';'ASPALAARGr';'ASPASNGLUr';'ASPGLUr';'ASPGLUPROr';'ASPGLUTRPr';'ASPHISCYSr';'ASPHISPROr';'ASPLYSGLUr';'ASPLYSHISr';'ASPMETASPr';'ASPPROLYSr';'ASPVALASNr';'CYSASNMETr';'CYSASPPHEr';'CYSCYSr';'CYSGLNMETr';'CYSGLUHISr';'CYSGLUTRPr';'CYSLEUTHRr';'CYSSERMETr';'CYSTYRASNr';'GLNASNGLNr';'GLNHISHISr';'GLNHISLYSr';'GLNLYSLYSr';'GLNLYSTRPr';'GLNPROGLUr';'GLNTRPGLUr';'GLNTYRLEUr';'GLUARGLEUr';'GLUASNLEUr';'GLUGLUr';'GLUILELYSr';'GLULEUr';'GLUMETr';'GLUMETHISr';'GLUTHRr';'GLUTHRLYSr';'GLUTRPALAr';'GLYHISASNr';'GLYHISLYSr';'GLYLYSCYSr';'GLYLYSPHEr';'GLYTYRLYSr';'GLYVALHISr';'HISARGCYSr';'HISARGSERr';'HISASPr';'HISCYSCYSr';'HISGLNALAr';'HISGLUr';'HISGLUGLNr';'HISGLYLYSr';'HISHISLYSr';'HISLYSALAr';'HISLYSGLUr';'HISLYSILEr';'HISLYSTHRr';'HISLYSVALr';'HISMETr';'HISMETGLNr';'HISPHEARGr';'HISPROLYSr';'HISTRPHISr';'ILEARGILEr';'ILEASNHISr';'ILEASPr';'ILEGLNGLUr';'ILEGLYARGr';'ILEPROLYSr';'ILESERARGr';'ILETRPTYRr';'LEUALAARGr';'LEUASNASPr';'LEUASPLYSr';'LEULEUTRPr';'LEUPROr';'LEUPROARGr';'LEUSERTRPr';'LEUTRPr';'LEUTRPARGr';'LEUTYRTYRr';'LEUVALr';'LYSARGLEUr';'LYSCYSHISr';'LYSGLNPHEr';'LYSGLUGLUr';'LYSLYSLYSr';'LYSPHEILEr';'LYSTRPARGr';'LYSTYRILEr';'LYSVALPHEr';'LYSVALTRPr';'METARGLEUr';'METASNTYRr';'METGLNTYRr';'METGLYARGr';'METHISLYSr';'METMETILEr';'METPHEARGr';'METTRPPHEr';'PHEASNMETr';'PHEASPr';'PHEGLNPHEr';'PHELEUr';'PHELEUASPr';'PHELEUHISr';'PHELYSALAr';'PHELYSPROr';'PHEPHEr';'PHEPHEASNr';'PHEPHETHRr';'PHEPROARGr';'PHESERTRPr';'PHETHRLYSr';'PHETRPLEUr';'PHETYRr';'PHETYRGLNr';'PHETYRLYSr';'PROARGASPr';'PROARGCYSr';'PROASNCYSr';'PROCYSr';'PROGLNPROr';'PROGLULYSr';'PROHISr';'PROHISTYRr';'PROLEUARGr';'PROLYSPROr';'PROPHEr';'PROPROARGr';'PROPROPROr';'PROTRPLYSr';'PROTRPTHRr';'PROVALGLNr';'SERARGALAr';'SERARGTRPr';'SERCYSARGr';'SERGLYGLUr';'SERLYSHISr';'SERPHELYSr';'SERTRPHISr';'THRARGTYRr';'THRASNTYRr';'THRGLNGLUr';'THRGLNTYRr';'THRHISHISr';'THRILEARGr';'THRMETARGr';'THRPHEARGr';'THRSERARGr';'THRTHRARGr';'THRTYRMETr';'TRPALAPROr';'TRPARGALAr';'TRPASPASPr';'TRPGLNGLNr';'TRPGLUGLYr';'TRPGLULEUr';'TRPGLUPROr';'TRPGLUTYRr';'TRPGLYLEUr';'TRPGLYPHEr';'TRPGLYVALr';'TRPHISMETr';'TRPILELYSr';'TRPILETRPr';'TRPLEUVALr';'TRPLYSr';'TRPMETARGr';'TRPMETVALr';'TRPPHEr';'TRPPROGLYr';'TRPPROLEUr';'TRPPROVALr';'TRPSERTYRr';'TRPTHRGLUr';'TRPTHRILEr';'TRPTHRTYRr';'TRPTYRGLNr';'TRPTYRTYRr';'TRPVALASPr';'TYRALAr';'TYRALAPHEr';'TYRARGGLUr';'TYRARGSERr';'TYRASPARGr';'TYRCYSGLYr';'TYRCYSTHRr';'TYRGLUr';'TYRLEUARGr';'TYRPHETYRr';'TYRTHRr';'TYRTRPPHEr';'TYRTYRr';'TYRVALMETr';'VALARGGLYr';'VALHISASNr';'VALLEUPHEr';'VALLYSTYRr';'VALPHEARGr';'VALPROTRPr';'VALSERARGr';'VALTRPPHEr';'VALTRPVALr';'VALVALr';'TRPGLYASPr';'r0643';'r0571';'UDPGLCAter'};

%% fibroblast
% Retrieved from https://doi.org/10.1111/febs.15292
model = readCbModel('models/fibroblast.mat');
model = changeRxnBounds(model,'biomass_maintenance',0,'l');

% add missing subsystems and reaction descriptions
for i=1:length(model.rxns)
    if strcmp(model.subSystems{i,1},'')
        findRxn = find(strcmp(vmhDB(:,1),model.rxns{i,1}));
        model.subSystems{i,1} = vmhDB{findRxn,4};
        if ~ischar(model.subSystems{i,1})
            model.subSystems{i,1} = model.subSystems{i,1}{1};
        end
    end
    if strcmp(model.rxns{i,1},model.rxnNames{i,1})
        findRxn = find(strcmp(vmhDB(:,1),model.rxns{i,1}));
        if ~isempty(findRxn)
            model.rxnNames{i,1} = vmhDB{findRxn,2};
        end
    end
end

[model,fbFlux] = add1cmPathways(model,'Human');

model = changeRxnBounds(model,unlIrrRxns,0,'l');

% enforce DNA synthesis and DNA methylation (uncoupled)
model = changeRxnBounds(model,'DM_dna[n]',1,'l');
model = changeRxnBounds(model,'DM_dna5mtc[n]',1,'l');

% writeCbModel removes coupling constraints
model = changeObjective(model,'biomass_maintenance');
save('fibroblast_CBL.mat','model')
writeCbModel(model,'fileName','fibroblast_CBL.mat','format','mat')

%% Recon3D flux-consistent model
% Retrieved from https://www.vmh.life/#downloadview
model = readCbModel('models/Recon3DModel_301.mat');

% add missing subsystems and reaction descriptions
for i=1:length(model.rxns)
    if strcmp(model.subSystems{i,1},'')
        findRxn = find(strcmp(vmhDB(:,1),model.rxns{i,1}));
        model.subSystems{i,1} = vmhDB{findRxn,4};
        if ~ischar(model.subSystems{i,1})
            model.subSystems{i,1} = model.subSystems{i,1}{1};
        end
    end
    if strcmp(model.rxns{i,1},model.rxnNames{i,1})
        findRxn = find(strcmp(vmhDB(:,1),model.rxns{i,1}));
        if ~isempty(findRxn)
            model.rxnNames{i,1} = vmhDB{findRxn,2};
        end
    end
end

[model,Recon3DFlux] = add1cmPathways(model,'Human');

model = changeRxnBounds(model,unlIrrRxns,0,'l');

% close sink reactions
model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'sink_',5))),0,'b');

% enforce DNA synthesis and DNA methylation (uncoupled)
model = changeRxnBounds(model,'DM_dna[n]',1,'l');
model = changeRxnBounds(model,'DM_dna5mtc[n]',1,'l');

model = changeObjective(model,'biomass_maintenance');

writeCbModel(model,'fileName','Recon3DModel_CBL.mat','format','mat')
writeCbModel(model,'fileName','Recon3DModel_CBL.xml','format','sbml')
