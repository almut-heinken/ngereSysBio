
% run simulations

mkdir('Results_brain_model')

modelFolder = 'Brain_specific_models';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

% perform flux variability analysis
minFluxes = {};
maxFluxes = {};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);

    % limit uptake of metabolites
    model=changeRxnBounds(model,model.rxns(find(contains(model.rxns,'_Ex_'))),-1,'l');
    model=changeRxnBounds(model,{'n_Ex_A_CO2','n_Ex_N_CO2','n_Ex_A_O2','n_Ex_N_O2'},-10,'l');
    [minFlux,maxFlux] = fastFVA(model,99,'max','ibm_cplex',model.rxns, 'S');
    minFluxes{i} = minFlux;
    maxFluxes{i} = maxFlux;
end
save(['Results_brain_model' filesep 'MinFluxes'],'minFluxes')
save(['Results_brain_model' filesep 'MaxFluxes'],'maxFluxes')

