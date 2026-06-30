
initCobraToolbox

solverOK=changeCobraSolver('ibm_cplex','LP');

draftFolder = [pwd filesep 'Draft_models'];

infoFilePath = [pwd filesep 'New_Strains.xlsx'];

inputDataFolder = [pwd filesep 'InputData'];

% [infoFilePath,inputDataFolder] = prepareInputData(infoFilePath,'inputDataFolder', inputDataFolder);

% after that, the files with experimental data were edited

numWorkers = 4;

reconVersion = 'DIFAMEM';

refinedFolder = [pwd filesep 'refinedReconstructions'];

translatedDraftsFolder = [pwd filesep 'translatedDraftReconstructions'];

runDemeter(draftFolder, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion, 'refinedFolder', refinedFolder, 'translatedDraftsFolder', translatedDraftsFolder);

testResultsFolder = [pwd filesep 'TestResults'];

[~,curationReport] = runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'translatedDraftsFolder', translatedDraftsFolder, 'numWorkers', numWorkers, 'testResultsFolder', testResultsFolder);

% needs to run 2x
[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,inputDataFolder,infoFilePath,reconVersion,'numWorkers',numWorkers);


% create SBML files
modelList = {'Alistipes_inops_627.mat';'Amedibacterium_intestinale_JCM_30884.mat';'Blautia_caecimuris_DSM_29492.mat';'Blautia_stercoris_3_YM_SP_D4_24_mj.mat';'Bradyrhizobium_liaoningense_CCBAU_05525.mat';'Brevundimonas_nasdae_JCM_11415.mat';'Burkholderia_glumae_ATCC_33617.mat';'Chryseobacterium_lactis_KC_1864.mat';'Citrobacter_europaeus_67A.mat';'Comamonas_kerstersii_J29.mat';'Coprobacter_secundus_subsp_similis_2CBH44.mat';'Enorma_phocaeensis_Marseille_P3242.mat';'Enteroscipio_rubneri_ResAG_96.mat';'Harryflintia_acetispora_DSM_100433.mat';'Hydrobacter_penzbergensis_DSM_25353.mat';'Marseilla_massiliensis_An824.mat';'Megasphaera_hexanoica_MH.mat';'Methylobacillus_flagellatus_KT.mat';'Parafannyhessea_umbonata_DSM_22619.mat';'Parasutterella_secunda_An562.mat';'Pseudarthrobacter_polychromogenes_CGMCC_1_1927.mat';'Pseudomonas_songnenensis_NEAU_ST5_5.mat';'Sanguibacteroides_justesenii_OUH_969102.mat';'Shewanella_baltica_OS678.mat';'Slackia_isoflavoniconvertens_DSM_22006.mat'};
mkdir('SBML_Files');

for i=1:length(modelList)
    model = readCbModel([refinedFolder filesep modelList{i}]);
    writeCbModel(model,'format','sbml','fileName',['SBML_Files' filesep strrep(modelList{i},'.mat','')])
end

