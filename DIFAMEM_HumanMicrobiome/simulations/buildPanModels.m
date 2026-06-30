
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

cd([pwd filesep 'Additional_strains'])

agoraPath = [pwd filesep 'refinedReconstructions'];

% create pan-species models for new strains
panPath = [pwd filesep 'panSpeciesModelsNew'];
createPanModels(agoraPath,panPath,'Species','New_Strains.xlsx',4)
cd ..

% copy to folder with existing AGORA2 pan-species models
copyfile(panPath, [pwd filesep 'PanSpeciesModels'])

% existing strains
agoraPath =  [pwd filesep 'AGORA2_corrected'];
infoFile = [pwd filesep 'input' filesep 'AGORA2_infoFile.xlsx'];
numWorkers = 12;

% create pan-species models
panPath = [pwd filesep 'panSpeciesModels'];
cd([pwd filesep 'Additional_strains'])
createPanModels(agoraPath,panPath,'Species',infoFile,numWorkers)
cd ..

% create pan-genus models
panPath = [pwd filesep 'panGenusModels'];
cd([pwd filesep 'Additional_strains'])
createPanModels(agoraPath,panPath,'Genus',infoFile,numWorkers)
cd ..

% workaround to build Escherichia/Shigella model
infoFile = [pwd filesep 'input' filesep 'AGORA2_infoFile.xlsx'];

panPath = [pwd filesep 'Escherichia_Shigella'];
cd([pwd filesep 'Additional_strains'])
createPanModels(agoraPath,panPath,'Genus',infoFile,numWorkers)
cd ..

