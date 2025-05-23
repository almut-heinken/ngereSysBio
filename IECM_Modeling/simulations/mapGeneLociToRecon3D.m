
clear all

% read in normalized and batch corrected data
data = readInputTableForPipeline(['input' filesep 'RNA_seq_data_RPKM.csv']);
data(:,1:3)=[];
data(:,3)=[];
data(1,:) = strrep(data(1,:),'/scratch/permanent/PIPELINEOUTPUT/projects/MMA_PHRT_R_BAK20190806/','');
for i=3:size(data,2)
colname = strsplit(data{1,i},'_');
data{1,i} = colname{1};
end

% read in gene symbol to Recon3D identifier mapping, downloaded from vmh.life
geneMapping = readInputTableForPipeline(['input' filesep 'recon-store-genes-1.tsv']);

genesMapped = {'GeneSymbol','ENTREZ','GeneNameRecon3D','Ensembl'};
cnt=2;

for i=2:size(data,1)
    ensembl = strsplit(data{i,1} ,'.');
    findGene = find(strcmp(geneMapping(:,12),ensembl{1}));
    if ~isempty(findGene)
        if isempty(find(strcmp(genesMapped(:,1),geneMapping{findGene(1),2})))
            genesMapped{cnt,1} = geneMapping{findGene(1),2};
            genesMapped{cnt,3} = geneMapping{findGene(1),1};
            genesMapped{cnt,4} = ensembl{1};
            entrez = strsplit(genesMapped{cnt,3} ,'.');
            genesMapped{cnt,2} = entrez{1};
            cnt=cnt+1;
        end
    end
end
genesMapped(:,2:3) = strrep(genesMapped(:,2:3),' ','');

% manually add the additionally mapped Cbl metabolism genes
cblGenes={
    'MMACHC','25974','25974.1',''
    'MMAB','326625','326625.1',''
    'MMADHC','27249','27249.1',''
    'CD320','51293','51293.1',''
    'TCN2','6948','6948.1',''
    'LMBR1','55788','55788.1',''
    'ABCD4','5826','5826.1',''
    };

genesMapped = vertcat(genesMapped,cblGenes);

save('Mapping_Genes_Recon3D','genesMapped');

