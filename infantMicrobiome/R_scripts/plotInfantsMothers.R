
library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")
setwd("..")

mycol <- c("white",brewer.pal(9,"Reds")[1:9])

### metadata infants
metadata = read.csv("inputFiles/Sample_metadata_with_adults.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)[c("Group","Time point")]

#### subsystem abundance ####
data=read.csv("Exported_results/InfantsMothers_SubsystemAbundance.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
  
# scale to data to rows
data=t(scale(t(data)))
  
dataReduced <- data
colsToKeep <- intersect(colnames(data),rownames(metadata))
dataReduced = dataReduced[,colsToKeep, drop = FALSE]
  
rowsToKeep <- intersect(rownames(metadata),colnames(data))
metadata = metadata[rowsToKeep,, drop = FALSE]
  
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(metadata$"Time point")))
ann.cols <- c(palette)
names(ann.cols) <- levels(metadata$"Time point")
column_ha = HeatmapAnnotation(df = NULL, TimePoint=as.character(metadata$"Time point"),
                                col=list(TimePoint =ann.cols))

dataReduced <- as.matrix(dataReduced)

ht = Heatmap(dataReduced,col=mycol,show_column_names = FALSE,top_annotation=column_ha,name = "Relative abundance",cluster_columns = FALSE,
               show_row_names = TRUE, width = unit(10, "cm"), row_names_max_width = unit(8, "cm"),row_names_gp = gpar(fontsize = 12))

png("R_plots/Figure_S2.png", width = 14, height = 6, units = 'in', res = 300)

draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")
  
dev.off()
