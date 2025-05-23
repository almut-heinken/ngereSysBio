
% simulate gene knockouts for the basic model for the different types of
% IECMs

clear all

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

iecms = {
    'WT',''
    'No_B12',''
    'CblA','166785.1'
    'CblB','326625.1'
    'CblC','25974.1'
    'CblD','27249.1'
    'CblE','4552.1'
    'CblF','55788.1'
    'CblG','4548.1'
    'CblJ','5826.1'
    'MMUT','4594.1'
    'CD320','51293.1'
    'TCN2','6948.1'
    };

mkdir([pwd filesep 'Results'])

model = readCbModel('Recon3DModel_CBL.mat');
forms = printRxnFormula(model);

% close sink reactions
model=changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'DM_',3))),0,'l');
model=changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'sink_',5))),0,'l');

flux = {};
minFluxes = {};
maxFluxes = {};

numWorkers=4;
if ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

for i=1:size(iecms,1)
    modelC = model;
    
    % simulate vitamin B12 deficiency
    if i==2
        modelC=changeRxnBounds(modelC,{'EX_aqcobal[e]','EX_ccbl[e]','EX_C06453[e]','EX_hoxocbl[e]','EX_adocbl[e]'},0,'l')';
    end
    
    % simulate the KO
    if i>2
        [modelC, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(modelC, iecms{i,2});
    end
    
    % compute minimal and maximal fluxes
    [minFlux,maxFlux] = fastFVA(modelC,99,'max','ibm_cplex',model.rxns, 'S');
    
    minFluxes{i} = minFlux;
    maxFluxes{i} = maxFlux;
    
    save([pwd filesep 'Results' filesep 'MinimalFluxes_KOs.mat'],'minFluxes')
    save([pwd filesep 'Results' filesep 'MaximalFluxes_KOs.mat'],'maxFluxes')
end

%% export fluxes

subsystems = {'Reaction','Subsystem','Formula'};
fluxes = ['Reaction',iecms(:,1)'];
for i=1:size(iecms,1)
    cnt=2;
    % minimal fluxes
    for j=1:length(model.rxns)
        fluxes{cnt,1} = [model.rxns{j,1}];
        subsystems{cnt,1} = [model.rxns{j,1}];
        subsystems{cnt,2} = model.subSystems{j,1};
        subsystems{cnt,3} = forms{j,1};
        fluxes{cnt,i+1} = minFluxes{i}(j,1);
        cnt=cnt+1;
    end
    
    % maximal fluxes
    for j=1:length(model.rxns)
        fluxes{cnt,1} = [model.rxns{j,1}];
        subsystems{cnt,1} = [model.rxns{j,1}];
        subsystems{cnt,2} = model.subSystems{j,1};
        subsystems{cnt,3} = forms{j,1};
        fluxes{cnt,i+1} = maxFluxes{i}(j,1);
        cnt=cnt+1;
    end
end
cell2csv([pwd filesep 'Results' filesep 'All_KO_Fluxes.csv'],fluxes)

% export altered fluxes-at least 5% difference
cnt = 1;
delArray = [];
for i=2:size(fluxes,1)
    wtFlux = fluxes{i,2};
    if ~any(cell2mat(fluxes(i,2:end))/wtFlux<0.95) ||~any(abs(cell2mat(fluxes(i,2:end)))>0.00001)
        delArray(cnt,1) = i;
        cnt = cnt+1;
    end
end
fluxes(delArray,:) = [];
subsystems(delArray,:) = [];

% scale to WT flux
for i=2:size(fluxes,1)
    wtFlux = fluxes{i,2};
    for j=2:size(fluxes,2)
       fluxes{i,j} = fluxes{i,j}/wtFlux;
    end
end
% remove reactions with constant minimal and maximal flux
[uniqueA,i,j] = unique(fluxes(:,1));
n  = accumarray(j(:),1);
Dupes=uniqueA(find(n>1));
indexToDupes = find(strcmp(fluxes(:,1),Dupes));
fluxes(indexToDupes(2),:)=[];
subsystems(indexToDupes(2),:) = [];

cell2csv([pwd filesep 'Results' filesep 'Altered_KO_Fluxes.csv'],fluxes)
writetable(cell2table(subsystems),[pwd filesep 'Results' filesep 'Altered_KO_Subsystems.csv'],'WriteVariableNames',false)
