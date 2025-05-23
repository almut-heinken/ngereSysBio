
% plot HMOs utilized by each AGORA2 microbe

mkdir([rootDir filesep 'Exported_results'])

hmos = {
    '2-Fucosyllactose', {'EX_2fuclac(e)'}
    '3-Fucosyllactose', {'EX_3fuclac(e)'}
    '3-Sialyllactose', {'EX_3slac(e)'}
    '6-Sialyllactose', {'EX_6slac(e)'}
    'Lacto-N-difucohexaose', {'EX_lacndfuchx(e)'}
    'Lacto-N-tetraose', {'EX_lacnttr(e)'}
    'Lacto-N-neotetraose', {'EX_lacnnttr(e)'}
    'Lactodifucotetraose', {'EX_lacdfucttr(e)'}
    'Lacto-N-biose', {'EX_lctnb(e)'}
    'Galacto-N-biose', {'EX_glctnb(e)'}
    'Lacto-N-fucopentaose I', {'EX_lacnfucpt_i(e)'}
    'Lacto-N-fucopentaose V', {'EX_lacnfucpt_v(e)'}
    'Lacto-N-fucopentaose II', {'EX_lacnfucpt_ii(e)'}
    'Lacto-N-fucopentaose III', {'EX_lacnfucpt_iii(e)'}
    'Disialyllacto-N-tetraose', {'EX_dslacnttr(e)'}
    'Lacto-N-triose I', {'EX_lntri_i(e)'}
    'Lacto-N-triose II', {'EX_lntri_ii(e)'}
    'Lacto-N-hexaose', {'EX_lacnhx(e)'}
    'Difucosyllacto-N-hexaose a', {'EX_dflnh_a(e)'}
    'Monofucosyllacto-N-hexaose', {'EX_fuclachx(e)'}
    'Sialyllacto-N-tetraose a', {'EX_neulacnttr(e)'}
    'Sialyllacto-N-tetraose b', {'EX_neulacnttr_b(e)'}
    'Sialyllacto-N-tetraose c', {'EX_neulacnttr_c(e)'}
    'Lacto-N-neohexaose', {'EX_lacnnhx(e)'}
    'Lacto-N-difucopentaose II', {'EX_lacndfucpt_ii(e)'}
    'Difucosyllactose', {'EX_dfuclac(e)'}
    'Disialyllactose', {'EX_dslac(e)'}
    'Monofucosylmonosialyllacto-N-hexaose', {'EX_fucneulacnhx(e)'}
    };

agora2GFFolder = [rootDir filesep 'AGORA2_gapfilled'];
dInfo = dir(agora2GFFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];
modelList=strrep(modelList,'.mat','');

hmoData = readInputTableForPipeline('HMOTable.txt');
[C,I]=setdiff(modelList,hmoData(:,1));
modelList(I,:)=[];

infoFile=readInputTableForPipeline('expanded_AGORA2_infoFile.xlsx');
hmosUsed = {'Strain'};

% loop through AGORA2 species and retrieve HMO utilization
for i=1:length(modelList)
    model=readCbModel([agora2GFFolder filesep modelList{i} '.mat']);
    hmosUsed{i+1,1}=modelList{i};
    for j=1:size(hmos,1)
        hmosUsed{1,j+1}=hmos{j,1};
        if ~isempty(find(strcmp(model.rxns,hmos{j,2})))
            model=changeObjective(model,hmos{j,2});
            FBA=optimizeCbModel(model,'min');
            if FBA.f <-0.000001
            hmosUsed{i+1,j+1}=1;
            else
                hmosUsed{i+1,j+1}=0;
            end
        else
            hmosUsed{i+1,j+1}=0;
        end
    end
end
% remove strains with no utilization
noUse = [];
cnt=1;
for i=2:size(hmosUsed,1)
    if sum(cell2mat(hmosUsed(i,2:end)))==0
        noUse(cnt,1)=i;
        cnt=cnt+1;
    end
end
hmosUsed(noUse,:)=[];
cell2csv([rootDir filesep 'Exported_results' filesep 'HMOsUsed.csv'],hmosUsed)

% export strain taxonomy
hmosAnnotation = infoFile;
[C,I] = setdiff(hmosAnnotation(:,1),hmosUsed(:,1),'stable');
hmosAnnotation(:,5) = strrep(hmosAnnotation(:,5),'[','');
hmosAnnotation(:,5) = strrep(hmosAnnotation(:,5),']','');
hmosAnnotation(I(2:end),:)=[];
hmosAnnotation(:,11:end)=[];
cell2csv([rootDir filesep 'Exported_results' filesep 'HMOsAnnotation.csv'],hmosAnnotation)
