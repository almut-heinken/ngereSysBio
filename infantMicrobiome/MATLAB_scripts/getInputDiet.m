% get an input file with the feeding for each sample

% path to the file with characteristics of the study participants
infoFilePath = [rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv'];

infoFile = readInputTableForPipeline(infoFilePath);
fCol = find(strcmp(infoFile(1,:),'Feeding (milk)'));
inputDiets = {};
for i=2:size(infoFile,1)
    inputDiets{i-1,1} = infoFile{i,1};
    if strcmp(infoFile{i,fCol},'Breast')
        if contains(infoFile{i,1},'V3')
            inputDiets{i-1,2} = 'breastmilk_5days';
        elseif contains(infoFile{i,1},'V4')
            inputDiets{i-1,2} = 'breastmilk_1month';
        elseif contains(infoFile{i,1},'V5')
            inputDiets{i-1,2} = 'breastmilk_6months';
        elseif contains(infoFile{i,1},'V6')
            inputDiets{i-1,2} = 'breastmilk_1year';
        end
    elseif strcmp(infoFile{i,fCol},'Formula')
        if contains(infoFile{i,1},'V3')
            inputDiets{i-1,2} = 'formulaFed_5days';
        elseif contains(infoFile{i,1},'V4')
            inputDiets{i-1,2} = 'formulaFed_1month';
        elseif contains(infoFile{i,1},'V5')
            inputDiets{i-1,2} = 'formulaFed_6months';
        elseif contains(infoFile{i,1},'V6')
            inputDiets{i-1,2} = 'formulaFed_1year';
        end
    elseif strcmp(infoFile{i,fCol},'Combined')
        if contains(infoFile{i,1},'V3')
            inputDiets{i-1,2} = 'breastmilk_5days';
        elseif contains(infoFile{i,1},'V4')
            inputDiets{i-1,2} = 'breastmilk_1month';
        elseif contains(infoFile{i,1},'V5')
            inputDiets{i-1,2} = 'formulaFed_6months';
        elseif contains(infoFile{i,1},'V6')
            inputDiets{i-1,2} = 'formulaFed_1year';
        end
    end
end
writetable(cell2table(inputDiets),'InputDiets','filetype','text','Delimiter','tab','writeVariableNames',false);
