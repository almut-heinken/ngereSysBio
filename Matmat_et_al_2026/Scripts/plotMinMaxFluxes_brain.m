
% plot and analyze the entire flux space to identify altered reactions
model = readCbModel('Data/iBrain674_Mm.xml');

% get some reaction annotations in the model
vmhDB = readInputTableForPipeline('Data/recon-store-reactions-1.tsv');
info = readInputTableForPipeline('Data/d0mo00135j2.xlsx');
info(:,1)=strrep(info(:,1),'-','_');
model.subSystems = cell(length(model.rxns),1);
for i=1:length(model.rxns)
    rxn=find(strcmp(info(:,1),strrep(model.rxns{i},'n_','')));
    ECnumbers = strsplit(info{rxn,6},'|');
    ECnumbers = strrep(ECnumbers,'EC-','');
    ECnumbers = strrep(ECnumbers,' ','');
    % find Ec numbers in global reconstruction
    [C,I] = intersect(vmhDB(:,5),ECnumbers);
    if ~isempty(I) && ~isempty(C{1})
        model.subSystems{i,1}=vmhDB{I(1),4};
    else
        model.subSystems{i,1} = info{rxn,2};
    end
end

Table = {'Reaction_ID','Reaction_Description','Subsystem'};
for i=1:length(model.rxns)
    Table{i,1} = model.rxns{i,1};
    Table{i,2} = model.rxnNames{i,1};
    Table{i,3} = model.subSystems{i,1};
end

modelFolder = 'Brain_specific_models';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

rxns = {};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    rxns = union(rxns,model.rxns);
end

[C,IA] = setdiff(Table(:,1),rxns);
Table(IA,:) = [];

load(['Results_brain_model' filesep 'MinFluxes'])
load(['Results_brain_model' filesep 'MaxFluxes'])

% plot fluxes by reaction
mkdir('MinMaxFluxes_plots_brain_model')

% define colors and samples
samples = readInputTableForPipeline('Data/SampleAnnotations.xlsx');
sampWT = samples(find(strcmp(samples(:,2),'WT')),1);
sampKO = samples(find(strcmp(samples(:,2),'KO')),1);

models={};
for i=1:length(modelList)
    model = readCbModel([modelFolder filesep modelList{i}]);
    models{i}=model;
end

cols = [0.8 0.4 0.4
    0.4 0.4 0.8
    0.8 0.4 0.4
    0.4 0.4 0.8
    ];

for i=2:size(Table,1)
    data = [];
    groupIdx=[];
    cnt=1;
    for j=1:size(sampWT,1)
        findMod = find(strcmp(modelList,[sampWT{j,1} '.mat']));
        model = models{findMod};
        rxnID = find(strcmp(model.rxns,Table{i,1}));
        if ~isempty(rxnID)
            data(cnt,1) = minFluxes{findMod}(rxnID,1);
            groupIdx(cnt,1) = 1;
            cnt=cnt+1;
            data(cnt,1) = maxFluxes{findMod}(rxnID,1);
            groupIdx(cnt,1) = 3;
            cnt=cnt+1;
        else
            data(cnt,1) = 0;
            groupIdx(cnt,1) = 1;
            cnt=cnt+1;
            data(cnt,1) = 0;
            groupIdx(cnt,1) = 3;
            cnt=cnt+1;
        end

    end
    for j=1:size(sampKO,1)
        findMod = find(strcmp(modelList,[sampKO{j,1} '.mat']));
        model = models{findMod};
        rxnID = find(strcmp(model.rxns,Table{i,1}));
        if ~isempty(rxnID)
            data(cnt,1) = minFluxes{findMod}(rxnID,1);
            groupIdx(cnt,1) = 2;
            cnt=cnt+1;
            data(cnt,1) = maxFluxes{findMod}(rxnID,1);
            groupIdx(cnt,1) = 4;
            cnt=cnt+1;
        else
            data(cnt,1) = 0;
            groupIdx(cnt,1) = 2;
            cnt=cnt+1;
            data(cnt,1) = 0;
            groupIdx(cnt,1) = 4;
            cnt=cnt+1;
        end

    end

    % only plot the reactions with flux
    if any(abs(data)>0.0001)
        f=figure;
        plot_box_scatter(data, groupIdx, 1:4, {'b','r','b','r'}, {'o','o','o','o'} , 1);
        xticklabels({'WT, minimal flux','KO, minimal flux','WT, maximal flux','KO, maximal flux'})
         ylim([min(data)*1.1,max(data)*1.1])
        set(gca,'TickLabelInterpreter','none')
        h=title(Table{i,1});
        set(h, 'Interpreter', 'none')
        hold on
        ylabel('mmol * g dry weight-1 * hr-1')
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),cols(j,:),'FaceAlpha',.5);
        end
        lcols=[0 0 1
            1 0 0];
        for j = 1:size(lcols,1)
            p(j) = patch(NaN, NaN, lcols(j,:));
        end
        legend(p,{'WT','KO'},'Location','northwest')

        f.Renderer='painters';
        print(['MinMaxFluxes_plots_brain_model' filesep Table{i,1}],'-dpng','-r300')
        close all
    end
end

