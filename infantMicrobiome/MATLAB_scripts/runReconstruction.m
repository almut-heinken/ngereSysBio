
% Runs DEMETER to refine and test 288 additional draft reconstructions
% generated for this study.

unzip('draftReconstructions.zip')

draftFolder = [rootDir filesep 'draftReconstructions'];

infoFilePath = [rootDir filesep 'inputReconstruction' filesep 'infoFile_Kraken.csv'];

inputDataFolder = [rootDir filesep 'inputReconstruction'];

numWorkers = 8;

reconVersion = 'Additional_Recons';

runDemeter(draftFolder, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion);
refinedFolder = [rootDir filesep 'refinedReconstructions'];
translatedDraftsFolder = [rootDir filesep 'translatedDraftReconstructions'];

[~,curationReport] = runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'translatedDraftsFolder', translatedDraftsFolder, 'numWorkers', numWorkers);
testResultsFolder = [rootDir filesep 'TestResults'];

[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,inputDataFolder,infoFilePath,reconVersion,'numWorkers',numWorkers);
