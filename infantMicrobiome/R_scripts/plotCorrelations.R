
library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")
  
#### plot species correlations #####
fpath = paste("Modeling_COSMIC/Correlations/FluxCorrelations_Secretion_Species.csv",sep="")
data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
fpath = paste("Modeling_COSMIC/Correlations/TaxonomyInfo_Secretion_Species.csv",sep="")
taxonomy=read.csv(fpath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum","Class")]
rowsToKeep <- intersect(rownames(taxonomy), rownames(data))
taxonomy = taxonomy[rowsToKeep,, drop = FALSE]
  
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
palette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(palette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)
  
mycol <- c(brewer.pal(6,"Blues")[6:1],"white",brewer.pal(9,"Reds")[1:9])
  
column_ha = HeatmapAnnotation(df = NULL, Phylum=as.character(taxonomy$Phylum),
                         Class=as.character(taxonomy$Class),
                         col=list(Phylum =PhylumAnn.cols,Class =ClassAnn.cols),
                         annotation_name_gp= gpar(fontsize = 14),
                         annotation_legend_param = list(Phylum = list(title_gp = gpar(fontsize = 14),
                                                                      labels_gp = gpar(fontsize = 12)),
                                                        Class = list(title_gp = gpar(fontsize = 14),
                                                                     labels_gp = gpar(fontsize = 12)))
)
  
png("R_plots/Figure_S1.png", width = 18, height =16, units = 'in', res = 300)
  
t <- data.matrix(t(data))
Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,name = "Correlation",column_names_side = c("top"),
          row_names_max_width = unit(14, "cm"),column_names_max_height = unit(14, "cm"),#column_names_rot = 45,
          show_column_names = TRUE,column_names_gp = gpar(fontsize = 12),top_annotation=column_ha, 
          heatmap_legend_param = list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12))
)
  
dev.off()
  
setwd("..")
  