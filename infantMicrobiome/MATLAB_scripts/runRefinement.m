

cd('inputFiles')
database = loadVMHDatabase;
cd ..

% AGORA2
% to get reconstructions: download from https://vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/zipped/
% then place all reconstructions in a folder named "AGORA2"

cd ..
agora2Folder = [rootDir filesep 'AGORA2'];
dInfo = dir(agora2Folder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

cd(rootDir)

% get all reactions and metabolites
totalContent=struct;

agora2GFFolder = [rootDir filesep 'AGORA2_gapfilled'];
mkdir(agora2GFFolder)
Gapfills = struct;

% parfor loop to reduce time
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox
end
solver = CBT_LP_SOLVER;
numWorkers=8;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(numWorkers)
end
environment = getEnvironment();

totalContent.('AGORA2')={
    'Reactions','Metabolites'
    '',''
    };

steps=100;

for i=1:steps:length(modelList)
    if length(modelList)-i>=steps-1
        endPnt=steps-1;
    else
        endPnt=length(modelList)-i;
    end
    
    modelsTmp = {};
    GapfillsTmp = {};
    allRxns={};
    allMets={};
    
    parfor j=i:i+endPnt
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        % load the model
        modelR = readCbModel([agora2Folder filesep modelList{j}]);
        
        % remove drug reactions
        %         remRxns = intersect(modelR.rxns,drugRxns);
        %         modelR = removeRxns(modelR,remRxns);
        
        microbeID = strrep(modelList{j},'.mat','');
        % run the HMO module function
        [modelR, addedRxns] = HMOGapfill(modelR, microbeID, database, [rootDir filesep 'inputFiles']);
        allRxns{j}=modelR.rxns;
        allMets{j}=modelR.mets;

        GapfillsTmp{j} = addedRxns;
        modelsTmp{j}=modelR;
    end
    % save the data
    for j=i:i+endPnt
        model=modelsTmp{j};
        % save model
        writeCbModel(model,'format','mat', 'fileName', [agora2GFFolder filesep modelList{j}])
        microbeID = strrep(modelList{j},'.mat','');
        if ~isempty(GapfillsTmp{j})
            Gapfills.(microbeID) = GapfillsTmp{j};
        end

        totalContent.('AGORA2'){2,1}=union(totalContent.('AGORA2'){2,1},allRxns{j});
        totalContent.('AGORA2'){2,2}=union(totalContent.('AGORA2'){2,2},allMets{j});
    end
end
save('AGORA2_Gapfills','Gapfills');

% test if this caused any futile cycles or abolished growth
[notGrowing,Biomass_fluxes] = plotBiomassTestResults(agora2GFFolder, 'AGORA2', 'numWorkers', numWorkers);
tooHighATP = plotATPTestResults(agora2GFFolder, 'AGORA2', 'numWorkers', numWorkers);

% get the additionally 289 created reconstructions
inputFolder = [pwd filesep 'refinedReconstructionsPublished'];

dInfo = dir(inputFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
for i=1:length(modelList)
    copyfile([inputFolder filesep modelList{i}],[agora2GFFolder filesep modelList{i}])
    model=readCbModel([agora2GFFolder filesep modelList{i}]);
    totalContent.('AGORA2'){2,1}=union(totalContent.('AGORA2'){2,1},model.rxns);
    totalContent.('AGORA2'){2,2}=union(totalContent.('AGORA2'){2,2},model.mets);
end

% test the gap-filled reconstructions
hmoData = readInputTableForPipeline('HMOTable.txt');
hmoReactions = readInputTableForPipeline('HMO_reactions.txt');

Flux = {};
NoFlux = {};
for i=2:size(hmoData,1)
    i
    Flux{i-1,1} = hmoData{i,1};
    NoFlux{i-1,1} = hmoData{i,1};
    model = readCbModel([agora2GFFolder filesep hmoData{i,1} '.mat']);
    rxns = intersect(model.rxns,hmoReactions(:,1));
    if ~isempty(rxns)
        try
            [minFlux,maxFlux] = fastFVA(model,0,'max','ibm_cplex',rxns);
        catch
            [minFlux,maxFlux] = fluxVariability(model,0,'max',rxns);
        end
        fluxRxns = {};
        fluxRxns = union(fluxRxns,rxns(minFlux<-0.00001));
        fluxRxns = union(fluxRxns,rxns(maxFlux>0.00001));
        noFluxRxns = setdiff(rxns,fluxRxns);
        Flux(i-1,2:length(fluxRxns)+1) = fluxRxns;
        NoFlux(i-1,2:length(noFluxRxns)+1) = noFluxRxns;
    end
end
save('Flux','Flux');
save('NoFlux','NoFlux');

% calculate number of added reactions
addedRxns=[];
cnt=1;

load('AGORA2_Gapfills.mat')
fn=fieldnames(Gapfills);
for i=1:length(fn)
    addedRxns(cnt,1)=length(Gapfills.(fn{i}));
    cnt=cnt+1;
end
avAdded=mean(addedRxns);
stdAdded=std(addedRxns);
AGORA2Added=length(fn);
