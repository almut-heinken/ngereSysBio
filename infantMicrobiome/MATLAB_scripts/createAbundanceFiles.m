
% Normalize infant gut microbiome data retrieved form Busi et al, ISME J
% Comm 2021

% Read in the mapping from the abundance data to AGORA2 taxa (created
% semi-automatically beforehand)
taxonMapping = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Kraken_taxa_mapped_to_AGORA2.csv']);

[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage([rootDir filesep 'inputFiles' filesep 'mappedInfantGutMicrobiome.csv'],0.001);
% exclude unmapped taxa from the normalized coverage
[C,I] = setdiff(normalizedCoverage(:,1),taxonMapping(:,1),'stable');
normalizedCoverage(I(2:end),:) = [];
for i=2:size(normalizedCoverage,1)
    normalizedCoverage{i,1} = taxonMapping{find(strcmp(taxonMapping(:,1),normalizedCoverage{i,1})),4};
    normalizedCoverage{i,1} = ['pan' normalizedCoverage{i,1}];
end
% summarize duplicate entries
[uniqueA,i,j] = unique(normalizedCoverage(:,1));
n  = accumarray(j(:),1);
Dupes=uniqueA(find(n>1));
delArray=[];
cnt=1;
for i=1:length(Dupes)
    indexToDupes = find(strcmp(normalizedCoverage(:,1),Dupes{i}));
    for j=2:length(indexToDupes)
        for k=2:size(normalizedCoverage,2)
                normalizedCoverage{indexToDupes(1),k}=normalizedCoverage{indexToDupes(1),k}+normalizedCoverage{indexToDupes(j),k};
        end
        delArray(cnt,1)=indexToDupes(j);
        cnt=cnt+1;
    end
end
normalizedCoverage(delArray,:)=[];
cell2csv('inputFiles/normalizedCoverage_COSMIC.csv',normalizedCoverage)


% normalize mapped mothers' microbiome data

% remove unmapped species after normalization
toRem = {'Bacteroides_dorei_vulgatus';'Bifidobacterium_catenulatum_Bifidobacterium_pseudocatenulatum_complex';'Enterobacter_hormaechei_cloacae';'butyrate_producing_bacterium_SS3_4';'Megasphaera_genomosp_type_1';'butyrate_producing_bacterium_SS3_4';'Streptococcus_sp_M143';'Streptococcus_sp_oral_taxon_056';'Klebsiella_variicola_pneumoniae';'Streptococcus_sp_oral_taxon_071';'Neisseria_gonorrhoeae';'Actinomyces_sp_oral_taxon_448';'Rahnella_sp_Y9602';'Enterobacter_sp_638';'Neisseria_sicca_macacae';'Propionibacterium_acnes_humerusii';'Actinomyces_sp_oral_taxon_170';'Actinomyces_sp_oral_taxon_180';'Burkholderia_ambifaria';'Burkholderia_dolosa';'Burkholderia_multivorans';'Serratia_plymuthica_odorifera';'unknown_SpeciesCluster_of_Oscillibacter';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Oscillibacter';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Lachnospiraceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Lachnospiraceae';'butyrate_producing_bacterium';'unknown_SpeciesCluster_of_Ruminococcaceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Veillonella';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Lachnospiraceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Ruminococcus';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Oscillibacter';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Oscillibacter';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Acidaminococcaceae';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Oscillibacter';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Coprobacillus';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Collinsella';'unknown_SpeciesCluster_of_Eubacterium';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Oscillospiraceae';'unclassified_Fusobacterium';'unknown_SpeciesCluster_of_Sutterella';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Sutterellaceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Peptostreptococcaceae';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Eubacterium';'unknown_SpeciesCluster_of_Coriobacteriaceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Lachnospiraceae';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Prevotella';'unknown_SpeciesCluster_of_Alistipes';'unknown_SpeciesCluster_of_Acidaminococcaceae';'unknown_SpeciesCluster_of_Ruminococcus';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Brachyspira';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Eubacterium';'unknown_SpeciesCluster_of_unknown_kingdom';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Lachnospiraceae';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroidales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteroides';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteroidetes';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Firmicutes';'unknown_SpeciesCluster_of_Roseburia';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridium';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Bacteria';'unknown_SpeciesCluster_of_Clostridia';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales';'unknown_SpeciesCluster_of_Clostridiales'};

[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage('inputFiles/Maternal_gut_microbiome.csv',0.001);
for i=2:size(normalizedCoverage,1)
    for j=2:size(normalizedCoverage,2)
    normalizedCoverage{i,j} = num2str(normalizedCoverage{i,j});
    end
end
[C,I] = intersect(normalizedCoverage,toRem);
normalizedCoverage(I,:) = [];

for i=2:size(normalizedCoverage,1)
    normalizedCoverage{i,1} = ['pan' normalizedCoverage{i,1}];
end

cell2csv('inputFiles/normalizedCoverage_MaternalGutMicrobiome.csv',normalizedCoverage)
