

% Run microbiome modeling pipeline consisting of creating an input file
% with relative strain abundances, creation of microbiome models,
% interrogation of microbiome models through simulations, and plotting of
% computed metabolite uptake and secretion by the communities.

mkdir([rootDir filesep 'Modeling_COSMIC'])
% run mgPipe workflow

% number of cores dedicated to parallelization
numWorkers = 4;

% define a folder where results will be saved
resPath = [rootDir filesep 'Modeling_COSMIC' filesep 'MicrobiomeModels'];
modPath = [rootDir filesep 'AGORA2_gapfilled'];

cd([rootDir filesep 'inputFiles'])
% create pan-species models
% first define list of species to build
SpeciesToBuild = {'Achromobacter insuavis';'Achromobacter sp. 2789STDY5608615';'Achromobacter sp. 2789STDY5608621';'Achromobacter xylosoxidans';'Actinomyces bouchesdurhonensis';'Actinomyces graevenitzii';'Actinomyces ihuae';'Actinomyces sp. HMSC035G02';'Actinomyces sp. HPA0247';'Actinomyces sp. ICM39';'Actinomyces sp. ICM47';'Actinomyces sp. ICM54';'Actinomyces sp. ICM58';'Actinomyces sp. S6-Spd3';'Actinomyces sp. oral taxon 172';'Actinomyces sp. oral taxon 181';'Actinomyces sp. ph3';'Actinomyces urogenitalis';'Actinotignum timonense';'Adlercreutzia equolifaciens';'Agathobaculum butyriciproducens';'Akkermansia muciniphila';'Alistipes communis';'Alistipes finegoldii';'Alistipes onderdonkii';'Alistipes putredinis';'Alistipes shahii';'Alistipes sp. AM16-43';'Alistipes sp. HGB5';'Anaerobutyricum hallii';'Anaerococcus hydrogenalis';'Anaerococcus jeddahensis';'Anaerococcus nagyae';'Anaerococcus obesiensis';'Anaerococcus rubeinfantis';'Anaerococcus sp. HMSC065G05';'Anaerococcus sp. HMSC068A02';'Anaerococcus sp. HMSC075B03';'Anaerococcus vaginalis';'Anaerostipes caccae';'Anaerostipes hadrus';'Anaerostipes sp. AF04-45';'Anaerotignum lactatifermentans';'Anaerotruncus Ruminococcaceae bacterium AF10-16';'Anaerotruncus Ruminococcaceae bacterium TF06-43';'Anaerotruncus colihominis';'Atlantibacter hermannii';'Atopobium sp. BS2';'Atopobium sp. ICM42b';'Bacteroides acidifaciens';'Bacteroides caccae';'Bacteroides caecimuris';'Bacteroides cellulosilyticus';'Bacteroides clarus';'Bacteroides dorei';'Bacteroides eggerthii';'Bacteroides faecis';'Bacteroides finegoldii';'Bacteroides fragilis';'Bacteroides helcogenes';'Bacteroides intestinalis';'Bacteroides ovatus';'Bacteroides rodentium';'Bacteroides sp. 14_A';'Bacteroides sp. 1_1_14';'Bacteroides sp. 1_1_30';'Bacteroides sp. 2_1_33B';'Bacteroides sp. 2_2_4';'Bacteroides sp. 3_1_33FAA';'Bacteroides sp. 3_1_40A';'Bacteroides sp. 4_1_36';'Bacteroides sp. AF14-46';'Bacteroides sp. AF17-1';'Bacteroides sp. AM37-9';'Bacteroides sp. D20';'Bacteroides sp. D22';'Bacteroides stercoris';'Bacteroides thetaiotaomicron';'Bacteroides timonensis';'Bacteroides uniformis';'Bacteroides xylanisolvens';'Barnesiella intestinihominis';'Bifidobacterium adolescentis';'Bifidobacterium angulatum';'Bifidobacterium animalis';'Bifidobacterium bifidum';'Bifidobacterium breve';'Bifidobacterium catenulatum';'Bifidobacterium dentium';'Bifidobacterium longum';'Bifidobacterium pseudocatenulatum';'Bifidobacterium scardovii';'Bifidobacterium sp. MSTE12';'Bifidobacterium sp. N4G05';'Bifidobacterium sp. N5G01';'Bilophila sp. 4_1_30';'Bilophila wadsworthia';'Blautia Lachnospiraceae bacterium 6_1_37FAA';'Blautia Lachnospiraceae bacterium AM10-38';'Blautia Lachnospiraceae bacterium AM21-21';'Blautia Lachnospiraceae bacterium AM25-17';'Blautia Lachnospiraceae bacterium AM26-1LB';'Blautia Lachnospiraceae bacterium AM48-27BH';'Blautia Lachnospiraceae bacterium GAM79';'Blautia Lachnospiraceae bacterium NLAE-zl-G231';'Blautia Lachnospiraceae bacterium OF11-28';'Blautia hansenii';'Blautia massiliensis';'Blautia obeum';'Blautia producta';'Blautia schinkii';'Blautia sp. AF13-16';'Blautia sp. AF14-40';'Blautia sp. AF17-9LB';'Blautia sp. AF19-1';'Blautia sp. AF22-5LB';'Blautia sp. AF26-2';'Blautia sp. AF32-4BH';'Blautia sp. AF34-10';'Blautia sp. AM23-13AC';'Blautia sp. AM29-29';'Blautia sp. AM42-2';'Blautia sp. AM47-4';'Blautia sp. BCRC 81119';'Blautia sp. KGMB01111';'Blautia sp. KLE 1732';'Blautia sp. Marseille-P2398';'Blautia sp. Marseille-P3201T';'Blautia sp. N6H1-15';'Blautia sp. OF03-13';'Blautia sp. OF03-15BH';'Blautia sp. OF09-25XD';'Blautia sp. OM05-6';'Blautia sp. SF-50';'Blautia sp. SG-772';'Blautia sp. aa_0143';'Blautia wexlerae';'Burkholderiales bacterium 1_1_47';'Butyricicoccus Clostridiaceae bacterium AF18-31LB';'Butyricicoccus Clostridiaceae bacterium AF29-16BH';'Butyricicoccus Clostridiaceae bacterium AF42-6';'Butyricicoccus Clostridiaceae bacterium AM27-36LB';'Butyricicoccus Clostridiaceae bacterium OF09-1';'Butyricicoccus Clostridiaceae bacterium OM02-2AC';'Butyricicoccus Clostridiaceae bacterium OM08-6BH';'Butyricicoccus Clostridiaceae bacterium TF01-6';'Butyricicoccus pullicaecorum';'Butyricicoccus sp. AF10-3';'Butyricicoccus sp. AM05-1';'Butyricicoccus sp. AM27-36';'Butyricicoccus sp. AM28-25';'Butyricicoccus sp. AM42-5AC';'Butyricicoccus sp. GAM44';'Campylobacter concisus';'Campylobacter gracilis';'Campylobacter jejuni';'Campylobacter sp. AAUH-44UCsig-a';'Candidatus Microthrix parvicella';'Candidatus Stoquefichus Firmicutes bacterium AF16-15';'Candidatus Stoquefichus Firmicutes bacterium AF22-6AC';'Candidatus Stoquefichus Firmicutes bacterium AF25-13AC';'Candidatus Stoquefichus Firmicutes bacterium AF36-19BH';'Candidatus Stoquefichus Firmicutes bacterium AF36-3BH';'Candidatus Stoquefichus Firmicutes bacterium AM55-24TS';'Candidatus Stoquefichus Firmicutes bacterium OM08-11AC';'Candidatus Stoquefichus Firmicutes bacterium TM09-10';'Capnocytophaga sp. oral taxon 329';'Citrobacter amalonaticus';'Citrobacter farmeri';'Citrobacter freundii';'Citrobacter koseri';'Citrobacter portucalensis';'Citrobacter sp. 30_2';'Citrobacter sp. CtB712';'Citrobacter youngae';'Clostridiales bacterium 1_7_47FAA';'Clostridioides difficile';'Clostridium baratii';'Clostridium botulinum';'Clostridium butyricum';'Clostridium carnis';'Clostridium celatum';'Clostridium chauvoei';'Clostridium cocleatum';'Clostridium disporicum';'Clostridium leptum';'Clostridium neonatale';'Clostridium nigeriense';'Clostridium paraputrificum';'Clostridium perfringens';'Clostridium phoceensis';'Clostridium saudiense';'Clostridium septicum';'Clostridium sp. 1_1_41A1FAA';'Clostridium sp. 7_2_43FAA';'Clostridium sp. 7_3_54FAA';'Clostridium sp. AF18-27';'Clostridium sp. AF27-2AA';'Clostridium sp. AF34-10BH';'Clostridium sp. AF34-13';'Clostridium sp. AF36-18BH';'Clostridium sp. AM22-11AC';'Clostridium sp. AM27-31LB';'Clostridium sp. AM30-24';'Clostridium sp. AM32-2';'Clostridium sp. AM33-3';'Clostridium sp. AM34-11AC';'Clostridium sp. AM34-9AC';'Clostridium sp. AT4';'Clostridium sp. ATCC BAA-442';'Clostridium sp. C105KSO14';'Clostridium sp. CL-2';'Clostridium sp. DSM 8431';'Clostridium sp. HGF2';'Clostridium sp. HMb25';'Clostridium sp. KLE 1755';'Clostridium sp. L2-50';'Clostridium sp. Marseille-P8228';'Clostridium sp. ND2';'Clostridium sp. OF03-18AA';'Clostridium sp. SN20';'Clostridium sp. SS2/1';'Clostridium sp. TF08-15';'Clostridium sp. TM06-18';'Clostridium tertium';'Collinsella aerofaciens';'Collinsella intestinalis';'Coprobacillus cateniformis';'Coprobacillus sp. 3_3_56FAA';'Coprobacillus sp. 8_2_54BFAA';'Coprobacillus sp. AF09-1A';'Coprobacillus sp. D7';'Coprococcus catus';'Coprococcus comes';'Coprococcus eutactus';'Coprococcus sp. AF18-48';'Coprococcus sp. AF27-8';'Coprococcus sp. AM27-12LB';'Coprococcus sp. HPP0048';'Coprococcus sp. HPP0074';'Corynebacterium aurimucosum';'Corynebacterium kroppenstedtii';'Corynebacterium propinquum';'Corynebacterium pseudodiphtheriticum';'Corynebacterium pseudogenitalium';'Corynebacterium tuberculostearicum';'Culturomica massiliensis';'Cutibacterium avidum';'Cutibacterium granulosum';'Dakarella massiliensis';'Dermabacter hominis';'Dermabacter sp. HFH0086';'Desulfovibrio piger';'Desulfovibrio sp. 3_1_syn3';'Dialister invisus';'Dielma fastidiosa';'Dolosigranulum pigrum';'Dorea formicigenerans';'Dorea longicatena';'Dorea sp. AGR2135';'Dorea sp. Marseille-P4003';'Dorea sp. Marseille-P4042';'Drancourtella massiliensis';'Dysgonomonas gadei';'Dysgonomonas mossii';'Eggerthella lenta';'Eggerthella sinensis';'Eggerthella sp. 1_3_56FAA';'Eikenella corrodens';'Eisenbergiella massiliensis';'Eisenbergiella sp. OF01-20';'Eisenbergiella tayi';'Enorma massiliensis';'Enterobacter asburiae';'Enterobacter cloacae';'Enterobacter hormaechei';'Enterobacter lignolyticus';'Enterobacter sp. Bisph2';'Enterobacter sp. NFR05';'Enterocloster aldenensis';'Enterocloster asparagiformis';'Enterocloster bolteae';'Enterocloster citroniae';'Enterocloster clostridioformis';'Enterococcus avium';'Enterococcus canintestini';'Enterococcus casseliflavus';'Enterococcus durans';'Enterococcus faecalis';'Enterococcus faecium';'Enterococcus gallinarum';'Enterococcus gilvus';'Enterococcus sp. 5B3_DIV0040';'Enterococcus sp. 6D12_DIV0197';'Enterococcus sp. HMSC05C03';'Erysipelatoclostridium innocuum';'Erysipelatoclostridium ramosum';'Erysipelotrichaceae bacterium 21_3';'Erysipelotrichaceae bacterium 2_2_44A';'Erysipelotrichaceae bacterium 5_2_54FAA';'Erysipelotrichaceae bacterium 6_1_45';'Escherichia coli';'Eubacterium limosum';'Eubacterium rectale';'Eubacterium sp. AF22-9';'Eubacterium sp. AM46-8';'Eubacterium sp. AM49-13BH';'Eubacterium sp. SB2';'Eubacterium ventriosum';'Eubacterium yurii';'Ezakiella Firmicutes bacterium AF12-30';'Faecalibacterium prausnitzii';'Faecalibacterium sp. AF28-13AC';'Faecalibacterium sp. OF04-11AC';'Faecalimonas umbilicata';'Finegoldia magna';'Flavonifractor plautii';'Fournierella massiliensis';'Fusicatenibacter saccharivorans';'Fusobacterium nucleatum';'Fusobacterium periodonticum';'Gardnerella vaginalis';'Gemella haemolysans';'Gemella morbillorum';'Gemella sanguinis';'Gemmiger formicilis';'Gordonibacter pamelaeae';'Gordonibacter sp. An230';'Gordonibacter sp. An232A';'Gordonibacter urolithinfaciens';'Granulicatella adiacens';'Granulicatella elegans';'Haemophilus aegyptius';'Haemophilus haemolyticus';'Haemophilus influenzae';'Haemophilus parahaemolyticus';'Haemophilus parainfluenzae';'Haemophilus paraphrohaemolyticus';'Haemophilus pittmaniae';'Haemophilus sp. CCUG 60358';'Haemophilus sp. CCUG 66565';'Haemophilus sp. HMSC061E01';'Haemophilus sp. HMSC068C11';'Haemophilus sp. HMSC71H05';'Haemophilus sp. oral taxon 036';'Haemophilus sputorum';'Holdemania filiformis';'Hungatella Clostridiaceae bacterium AF31-3BH';'Hungatella effluvii';'Hungatella hathewayi';'Intestinibacillus sp. Marseille-P4005';'Intestinibacillus sp. Marseille-P6563';'Intestinibacter bartlettii';'Intestinimonas Clostridiales bacterium AM23-16LB';'Intestinimonas Clostridiales bacterium CCNA10';'Intestinimonas Clostridiales bacterium Choco116';'Intestinimonas Clostridiales bacterium KLE1615';'Intestinimonas Clostridiales bacterium Marseille-P5551';'Intestinimonas Clostridiales bacterium VE202-01';'Intestinimonas Clostridiales bacterium VE202-03';'Intestinimonas Clostridiales bacterium VE202-18';'Intestinimonas butyriciproducens';'Isoptericola variabilis';'Klebsiella aerogenes';'Klebsiella grimontii';'Klebsiella michiganensis';'Klebsiella oxytoca';'Klebsiella pneumoniae';'Klebsiella sp. WCHKl090539';'Klebsiella variicola';'Kluyvera intermedia';'Kluyvera sp. Nf5';'Lachnoclostridium edouardi';'Lachnoclostridium lavalense';'Lachnoclostridium phocaeense';'Lachnoclostridium sp. An14';'Lachnoclostridium sp. An169';'Lachnoclostridium sp. An298';'Lachnoclostridium sp. An76';'Lachnospira eligens';'Lachnospira pectinoschiza';'Lachnospiraceae bacterium 1_4_56FAA';'Lachnospiraceae bacterium 5_1_57FAA';'Lachnospiraceae bacterium 5_1_63FAA';'Lachnospiraceae bacterium 6_1_63FAA';'Lachnospiraceae bacterium 7_1_58FAA';'Lachnospiraceae bacterium 9_1_43BFAA';'Lacrimispora indolis';'Lacrimispora saccharolytica';'Lacticaseibacillus casei';'Lacticaseibacillus paracasei';'Lacticaseibacillus rhamnosus';'Lactiplantibacillus plantarum';'Lactobacillus delbrueckii';'Lactobacillus gasseri';'Lactobacillus johnsonii';'Lactobacillus sp. HMSC068B07';'Lactococcus lactis';'Lancefieldella parvula';'Lawsonella clevelandensis';'Lawsonibacter Clostridiales bacterium VE202-07';'Lawsonibacter Clostridiales bacterium VE202-16';'Lawsonibacter asaccharolyticus';'Ligilactobacillus salivarius';'Limosilactobacillus fermentum';'Limosilactobacillus reuteri';'Mediterraneibacter glycyrrhizinilyticus';'Megasphaera micronuciformis';'Merdibacter Firmicutes bacterium AM59-13';'Metakosakonia sp. MRY16-398';'Methanobrevibacter smithii';'Mogibacterium diversum';'Mycobacterium tuberculosis';'Negativicoccus massiliensis';'Neisseria cinerea';'Neisseria flavescens';'Neisseria meningitidis';'Neisseria subflava';'Neobitarella massiliensis';'Niameybacter massiliensis';'Odoribacter splanchnicus';'Oscillibacter Oscillospiraceae bacterium VE202-24';'Oscillibacter sp. 1-3';'Oscillibacter sp. PC13';'Oscillibacter sp. PEA192';'Oxobacter Clostridiales bacterium VE202-26';'Paeniclostridium sordellii';'Parabacteroides distasonis';'Parabacteroides merdae';'Parabacteroides sp. D13';'Parascardovia denticolens';'Parasutterella excrementihominis';'Parvibacter Coriobacteriaceae bacterium CHKCI002';'Pediococcus acidilactici';'Pediococcus pentosaceus';'Peptoniphilus grossensis';'Peptoniphilus harei';'Peptostreptococcus anaerobius';'Phascolarctobacterium faecium';'Phascolarctobacterium succinatutens';'Phocaeicola coprocola';'Phocaeicola plebeius';'Phocaeicola sartorii';'Phocaeicola vulgatus';'Phytobacter sp. SCO41';'Phytobacter ursingii';'Prevotella buccae';'Prevotella copri';'Propionibacterium freudenreichii';'Propionibacterium sp. HGH0353';'Proteus mirabilis';'Pseudoflavonifractor capillosus';'Pseudoflavonifractor sp. An184';'Pseudomonas aeruginosa';'Pseudomonas yamanorum';'Raoultella ornithinolytica';'Raoultella planticola';'Raoultella terrigena';'Robinsoniella peoriensis';'Romboutsia sp. Marseille-P6047';'Romboutsia timonensis';'Roseburia faecis';'Roseburia hominis';'Roseburia intestinalis';'Roseburia inulinivorans';'Rothia mucilaginosa';'Rothia sp. HMSC061D12';'Rothia sp. HMSC061E04';'Rothia sp. HMSC062F03';'Rothia sp. HMSC064F07';'Rothia sp. HMSC065C03';'Rothia sp. HMSC068F09';'Rothia sp. HMSC072B03';'Rothia sp. HMSC072B04';'Rothia sp. HMSC072E10';'Rubneribacter badeniensis';'Ruminococcaceae bacterium D16';'Ruminococcus bicirculans';'Ruminococcus bromii';'Ruminococcus faecis';'Ruminococcus flavefaciens';'Ruminococcus sp. 5_1_39BFAA';'Ruminococcus sp. AF14-5';'Ruminococcus sp. AF17-12';'Ruminococcus sp. AF20-12LB';'Ruminococcus sp. AF21-42';'Ruminococcus sp. AF31-8BH';'Ruminococcus sp. AF41-9';'Ruminococcus sp. AF46-10NS';'Ruminococcus sp. AM16-34';'Ruminococcus sp. AM26-12LB';'Ruminococcus sp. AM41-10BH';'Ruminococcus sp. AM42-11';'Ruminococcus sp. AM46-18';'Ruminococcus sp. DSM 100440';'Ruminococcus sp. OF02-6';'Ruminococcus sp. OF03-6AA';'Ruminococcus sp. OM02-16LB';'Ruminococcus sp. OM04-4AA';'Ruminococcus sp. OM08-9BH';'Ruminococcus sp. SR1/5';'Ruminococcus sp. TF06-23';'Ruminococcus sp. TF11-2AC';'Ruthenibacterium lactatiformans';'Salmonella enterica';'Sanguibacter keddieii';'Scardovia wiggsiae';'Schaalia odontolytica';'Schaalia turicensis';'Sedimentibacter Tissierellia bacterium S5-A11';'Sellimonas intestinalis';'Serratia marcescens';'Shigella sonnei';'Staphylococcus aureus';'Staphylococcus capitis';'Staphylococcus epidermidis';'Staphylococcus haemolyticus';'Staphylococcus hominis';'Staphylococcus lugdunensis';'Staphylococcus pasteuri';'Staphylococcus saprophyticus';'Staphylococcus sp. HMSC070D05';'Streptococcus agalactiae';'Streptococcus anginosus';'Streptococcus australis';'Streptococcus constellatus';'Streptococcus cristatus';'Streptococcus equinus';'Streptococcus gallolyticus';'Streptococcus gordonii';'Streptococcus infantarius';'Streptococcus infantis';'Streptococcus intermedius';'Streptococcus lutetiensis';'Streptococcus macedonicus';'Streptococcus mitis';'Streptococcus oralis';'Streptococcus parasanguinis';'Streptococcus pasteurianus';'Streptococcus peroris';'Streptococcus pneumoniae';'Streptococcus pseudopneumoniae';'Streptococcus salivarius';'Streptococcus sanguinis';'Streptococcus sp. 1171_SSPC';'Streptococcus sp. 263_SSPC';'Streptococcus sp. 400_SSPC';'Streptococcus sp. A12';'Streptococcus sp. ACS2';'Streptococcus sp. C150';'Streptococcus sp. CCH5-D3';'Streptococcus sp. CCH8-H5';'Streptococcus sp. F0442';'Streptococcus sp. FDAARGOS_192';'Streptococcus sp. HMSC057G03';'Streptococcus sp. HMSC061E03';'Streptococcus sp. HMSC064H09';'Streptococcus sp. HMSC065C01';'Streptococcus sp. HMSC065E03';'Streptococcus sp. HMSC072C09';'Streptococcus sp. HMSC072G04';'Streptococcus sp. HMSC073D05';'Streptococcus sp. HMSC074F05';'Streptococcus sp. HMSC076C09';'Streptococcus sp. HMSC078D09';'Streptococcus sp. HMSC078H03';'Streptococcus sp. HMSC078H12';'Streptococcus sp. HMSC10E12';'Streptococcus sp. HMSC34B10';'Streptococcus sp. HPH0090';'Streptococcus sp. I-P16';'Streptococcus sp. M334';'Streptococcus sp. SK140';'Streptococcus sp. SK643';'Streptococcus sp. SR4';'Streptococcus sp. UMB0029';'Streptococcus sp. UMB1385';'Streptococcus sp. bf_0095';'Streptococcus suis';'Streptococcus thermophilus';'Streptococcus timonensis';'Streptococcus vestibularis';'Streptococcus viridans';'Subdoligranulum sp. 4_3_54A2FAA';'Subdoligranulum sp. AF14-43';'Subdoligranulum sp. AM16-9';'Subdoligranulum sp. APC924/74';'Subdoligranulum sp. OF01-18';'Subdoligranulum variabile';'Sutterella sp. KLE1602';'Sutterella wadsworthensis';'Terrisporobacter glycolicus';'Trueperella pyogenes';'Turicibacter sanguinis';'Tyzzerella nexilis';'Tyzzerella sp. Marseille-P3062';'Unknown Burkholderiales bacterium';'Unknown Candidatus Saccharibacteria bacterium';'Vallitalea Clostridia bacterium UC51-1D1';'Vallitalea Clostridia bacterium UC51-1D10';'Vallitalea Clostridia bacterium UC51-2G4';'Vallitalea Clostridia bacterium UC51-2H11';'Vallitalea Clostridia bacterium UC51-2H6';'Varibaculum cambriense';'Varibaculum sp. Marseille-P2802';'Veillonella atypica';'Veillonella caviae';'Veillonella dispar';'Veillonella infantium';'Veillonella parvula';'Veillonella ratti';'Veillonella rogosae';'Veillonella seminalis';'Veillonella sp. 3_1_44';'Veillonella sp. 6_1_27';'Veillonella sp. ACP1';'Veillonella sp. AF13-2';'Veillonella sp. AF36-20BH';'Veillonella sp. AF42-16';'Veillonella sp. AM51-8BH';'Veillonella sp. HPA0037';'Veillonella sp. ICM51a';'Veillonella sp. S13053-19';'Veillonella sp. S13054-11';'Veillonella sp. T11011-6';'Veillonella sp. T14073-2';'Veillonella sp. T34266-5';'Veillonella sp. oral taxon 158';'Veillonella sp. oral taxon 780';'Veillonella tobetsuensis';'Xylanimonas cellulosilytica';'[Clostridium] methoxybenzovorans';'[Clostridium] saccharogumia';'[Clostridium] scindens';'[Clostridium] spiroforme';'[Clostridium] sporosphaeroides';'[Clostridium] symbiosum';'[Enterobacter] lignolyticus';'[Eubacterium] rectale';'[Eubacterium] siraeum';'[Eubacterium] yurii';'[Ruminococcus] gnavus';'[Ruminococcus] lactaris';'[Ruminococcus] torques'};
panPath = [rootDir filesep 'panSpeciesModels'];
createPanModels(modPath,panPath,'Species',[rootDir filesep 'inputFiles' filesep 'expanded_AGORA2_infoFile.xlsx'],numWorkers,SpeciesToBuild)
cd(rootDir)

% test pan-species models
[notGrowing,biomassFluxes] = plotBiomassTestResults(panPath, 'Species', 'numWorkers',numWorkers);
[tooHighATP,atpFluxes] = plotATPTestResults(panPath, 'Species', 'numWorkers',numWorkers);

% path to and name of the file with abundance information.
abunFilePath = [rootDir filesep 'inputFiles' filesep 'normalizedCoverage_COSMIC.csv'];

% path to the file with characteristics of the study participants
infoFilePath = [rootDir filesep 'inputFiles' filesep 'Sample_metadata.csv'];

% path to the simulated diet
dietFilePath = [rootDir filesep 'inputFiles' filesep 'inputDiets_COSMIC.txt'];

% lower the required minimum biomass production
lowerBMBound = 0.2;

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(panPath, abunFilePath, true, 'resPath', resPath, 'infoFilePath', infoFilePath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers, 'lowerBMBound', lowerBMBound);
writetable(cell2table(netSecretionFluxes),[resPath filesep 'Infant_diet.csv'],'writeVariableNames',false)

%% determine taxon-metabolite correlations
mkdir([rootDir filesep 'Modeling_COSMIC' filesep 'Correlations'])
taxInfo = [rootDir filesep 'inputFiles' filesep 'expanded_AGORA2_infoFile.xlsx'];
fluxPath = [rootDir filesep 'Modeling_COSMIC' filesep 'MicrobiomeModels' filesep 'Infant_diet.csv'];
corrMethod = 'Spearman';
% first with net secretion
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, taxInfo, corrMethod);

% export the results
taxa=fieldnames(FluxCorrelations);
for i=1:length(taxa)
    % first correct for multiple testing
    cnt=1;
    pvals = [];
    for j=2:size(PValues.(taxa{i}),1)
        for k=2:size(PValues.(taxa{i}),2)
            pvals(cnt,1)=PValues.(taxa{i}){j,k};
            cnt=cnt+1;
        end
    end
    pvals = mafdr(pvals,'BHFDR', true);
    cnt=1;
    for j=2:size(PValues.(taxa{i}),1)
        for k=2:size(PValues.(taxa{i}),2)
            PValues.(taxa{i}){j,k} = pvals(cnt,1);
            cnt=cnt+1;
        end
    end
    % then remove lesser correlations
    cnt=1;
    delArray=[];
    for j=2:size(FluxCorrelations.(taxa{i}),2)
        if ~any(abs(cell2mat(FluxCorrelations.(taxa{i})(2:end,j))) > 0.8)
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    FluxCorrelations.(taxa{i})(:,delArray)=[];

    cnt=1;
    delArray=[];
    for j=2:size(FluxCorrelations.(taxa{i}),1)
        if ~any(abs(cell2mat(FluxCorrelations.(taxa{i})(j,2:end))) > 0.8)
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    FluxCorrelations.(taxa{i})(delArray,:)=[];
    if i>2
        [C,I]=setdiff(TaxonomyInfo.(taxa{i})(:,1),FluxCorrelations.(taxa{i})(2:end,1),'stable');
        TaxonomyInfo.(taxa{i})(I(2:end),:)=[];
    end

    FluxCorrelations.(taxa{i})(find(strcmp(FluxCorrelations.(taxa{i})(:,1),'')),:)=[];
    writetable(cell2table(FluxCorrelations.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'FluxCorrelations_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    writetable(cell2table(PValues.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'PValues_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    if i>1
        writetable(cell2table(TaxonomyInfo.(taxa{i})),[rootDir filesep 'Modeling_COSMIC' filesep 'Correlations' filesep 'TaxonomyInfo_Secretion_' taxa{i} '.csv'],'writeVariableNames',false)
    end
end

%% for comparison with adults, compute fluxes on the standard Western diet

dietFilePath = 'AverageEuropeanDiet.txt';

delete([resPath filesep 'simRes.mat'],[resPath filesep 'intRes.mat'])

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] =  initMgPipe(panPath, abunFilePath, true, 'resPath', resPath, 'infoFilePath', infoFilePath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);
writetable(cell2table(netSecretionFluxes),[resPath filesep 'AE_diet.csv'],'writeVariableNames',false)

%% compute microbe secretion for metabolites of interest
% vitamins, SCFAs
mets = {'ac','ppa','but','lac_L','fol','adocbl','ribflv','thm','nac','btn','pydx','pnto_R'};
constrModPath = [resPath filesep 'Diet'];
contrPath = [rootDir filesep 'Modeling_COSMIC' filesep 'Contributions'];

% make changes to the model to prevent false prediction of riboflavin
dInfo = dir(constrModPath);
fileList={dInfo.name};
fileList=fileList';
fileList(~contains(fileList(:,1),'.mat'),:)=[];
for i=1:length(fileList)
    load([constrModPath filesep fileList{i}])
    rxns=model.rxns(find(contains(model.rxns,'IEX_rbflvrd[u]tr')));
    model=changeRxnBounds(model,rxns,0,'u');
    save([constrModPath filesep fileList{i}],'model')
end

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath, 'metList', mets, 'numWorkers', numWorkers,'resultsFolder',contrPath);
save([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'minFluxes.mat'],'minFluxes')
cell2csv([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)

% minFluxes = secretion
minFluxes=readInputTableForPipeline([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv']);
for i=2:size(minFluxes,1)
    for j=2:size(minFluxes,2)
        if minFluxes{i,j} < 0
            minFluxes{i,j} = abs(minFluxes{i,j});
        elseif minFluxes{i,j} > 0
            minFluxes{i,j} = 0;
        end
    end
end
cell2csv([rootDir filesep 'Modeling_COSMIC' filesep 'Contributions' filesep 'Microbe_Secretion.csv'],minFluxes)
