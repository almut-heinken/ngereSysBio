
% rebuild AGORA2 with corrected reaction

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

mkdir('AGORA2_corrected')

% rebuild with updated reaction database
cd('Additional_strains')
database=loadVMHDatabase;
cd ..

changed = {};

% first VMH version AGORA2
agora2Folder = 'C:\Users\Almut\Universite de Lorraine\Systems Biology NGERE - Documents\AGORA2_VMH\mat_files';

dInfo = dir(agora2Folder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

for i=1:length(modelList)
    i
    model=readCbModel([agora2Folder filesep modelList{i}]);
    if ~isempty(find(strcmp(model.rxns,'DHPHEOGAT')))
        model=rebuildModel(model,database);
        changed{i} = modelList{i};
    end
    writeCbModel(model,'format','mat','fileName',['AGORA2_corrected' filesep modelList{i}])
end

% then the expansion
agora2Folder = 'C:\Users\Almut\Universite de Lorraine\Systems Biology NGERE - Documents\AGORA2_Expansion_InfantMicrobiome\refinedReconstructions';

dInfo = dir(agora2Folder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

for i=1:length(modelList)
    i
    model=readCbModel([agora2Folder filesep modelList{i}]);
    if ~isempty(find(strcmp(model.rxns,'DHPHEOGAT')))
        model=rebuildModel(model,database);
        changed{i} = modelList{i};
    end
    writeCbModel(model,'format','mat','fileName',['AGORA2_corrected' filesep modelList{i}])
end

save('corrected_AGORA2','changed')

% test if everything still works
[notGrowing,Biomass_fluxes] = plotBiomassTestResults('AGORA2_corrected', 'AGORA2','numWorkers', 4);
[tooHighATP,ATP_fluxes] = plotATPTestResults('AGORA2_corrected', 'AGORA2','numWorkers', 4);
