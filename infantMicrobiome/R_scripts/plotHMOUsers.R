
library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

### HMO utilization by strain
data = read.csv("Exported_results/HMOsUsed.csv",header=T,row.names=1, check.names = FALSE)
taxonomy=read.csv("Exported_results/HMOsAnnotation.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum","Species")]

getPalette = colorRampPalette(brewer.pal(8, "Accent"))
palette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(palette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
SpeciesPalette = getPalette(length(levels(taxonomy$Species)))
SpeciesAnn.cols <- c(SpeciesPalette)
names(SpeciesAnn.cols) <- levels(taxonomy$Species)

mycol <- c("white",brewer.pal(5,"Blues")[1:5])

column_ha = HeatmapAnnotation(df = NULL, Phylum=as.character(taxonomy$Phylum),
                              Species=as.character(taxonomy$Species),
                              col=list(Phylum =PhylumAnn.cols,Species =SpeciesAnn.cols),
                              annotation_name_side = "right",
                              annotation_name_gp= gpar(fontsize = 14),
                              annotation_legend_param = list(Phylum = list(title_gp = gpar(fontsize = 14),
                                                                           labels_gp = gpar(fontsize = 14)),
                                                             Species = list(title_gp = gpar(fontsize = 14),
                                                                          labels_gp = gpar(fontsize = 14)))
)

png("R_plots/Figure_1b.png", width = 16, height =8, units = 'in', res = 300)

t <- data.matrix(t(data))
Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = FALSE, top_annotation=column_ha, 
        row_names_max_width = unit(12, "cm"),column_names_max_height = unit(12, "cm"),border = "grey",
        show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14))
)
dev.off()
