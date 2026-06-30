
% map the reads to taxonomy and create input file for mgPipe

data = readInputTableForPipeline(['input' filesep 'feature-table.txt']);
taxonomy = readInputTableForPipeline(['input' filesep 'taxonomy.tsv']);

% remove controls
metadata = readInputTableForPipeline(['input' filesep 'sample-metadata.txt']);
controls = metadata(find(strcmp(metadata(:,2),'Control')),1);
[C,I] = intersect(data(1,:),controls);
data(:,I) = [];

for i=2:size(data,1)
    findID = find(strcmp(taxonomy(:,1),data{i,1}));
    data{i,1} = taxonomy{findID,2};
end
cell2csv('taxFeatureTable.csv',data)

% normalize the coverage to remove species below the cutoff
[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage('taxFeatureTable.csv',0);

% read mapping to AGORA2
mapping = readInputTableForPipeline(['input' filesep 'mapping2AGORA2.xlsx']);

for i=2:size(normalizedCoverage,1)
    normalizedCoverage{i,1} = mapping{find(strcmp(mapping(:,2),normalizedCoverage{i,1})),3};
end

normalizedCoverage(find(cellfun(@isempty,normalizedCoverage(:,1))),:) = [];

for i=2:size(normalizedCoverage,1)
    normalizedCoverage{i,1} = ['pan' normalizedCoverage{i,1}];
end

writetable(cell2table(normalizedCoverage),[pwd filesep 'normalizedCoverage.csv'],'writeVariableNames',false)

% calculate the total captured coverage
coverage = [];
for i=2:size(normalizedCoverage,2)
    coverage(i-1,1) = sum(cell2mat(normalizedCoverage(2:end,i)));
end

mean(coverage)
