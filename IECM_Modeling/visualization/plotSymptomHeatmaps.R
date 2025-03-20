
library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

mycol <- c("white",brewer.pal(9,"Blues")[1:9])

### Fluxes

#### All symptoms MMUT fluxes ####
data=read.csv("Statistics_Fluxes/MMUT_Symptoms.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

#### read reaction subsystems
subsystems=read.csv("Statistics_Fluxes/MMUT_Symptoms_Subsystems.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)
rowsToKeep <- intersect(rownames(subsystems), rownames(data))
subsystems = subsystems[rowsToKeep,, drop = FALSE]

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(subsystems$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(subsystems$Subsystem)

column_ha = HeatmapAnnotation(df = NULL, col=list(Subsystem=SubsystemAnn.cols),
                       Subsystem=as.character(subsystems$Subsystem))

png("Figure_4a.png", width = 14, height = 5, units = 'in', res = 300)
t <- t(data.matrix(data))
Heatmap(t,col=mycol,show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE,
        name = "Significant reactions", show_heatmap_legend = FALSE, top_annotation=column_ha,
        show_row_names =TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 12),
        row_names_max_width = unit(10, "cm"))
dev.off()

#### All symptoms other MMA ####
data=read.csv("Statistics_Fluxes/MMA_Other_Symptoms.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

#### read reaction subsystems
subsystems=read.csv("Statistics_Fluxes/MMA_Other_Symptoms_Subsystems.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)
rowsToKeep <- intersect(rownames(subsystems), rownames(data))
subsystems = subsystems[rowsToKeep,, drop = FALSE]

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(subsystems$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(subsystems$Subsystem)

column_ha = HeatmapAnnotation(df = NULL, col=list(Subsystem=SubsystemAnn.cols),
                              Subsystem=as.character(subsystems$Subsystem))

png("Figure_4b.png", width = 14, height = 5, units = 'in', res = 300)
t <- t(data.matrix(data))
Heatmap(t,col=mycol,show_column_names = TRUE, cluster_columns = TRUE, cluster_rows = TRUE,border = TRUE,
        name = "Significant reactions", show_heatmap_legend = FALSE, top_annotation=column_ha,
        show_row_names =TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10)
)
dev.off()

### Metabolites

data=read.csv("Statistics_Metabolites/MMUT_Symptoms.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

png("Figure_S2.png", width = 10, height = 4, units = 'in', res = 300)
t <- t(data.matrix(data))
Heatmap(t,col=mycol,show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE,
        name = "Significant reactions", show_heatmap_legend = FALSE,
        show_row_names =TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 12),
        row_names_max_width = unit(10, "cm"))
dev.off()

data=read.csv("Statistics_Metabolites/MMA_Other_Symptoms.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

png("Figure_S3.png", width = 10, height = 4, units = 'in', res = 300)
t <- t(data.matrix(data))
Heatmap(t,col=mycol,show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE,
        name = "Significant reactions", show_heatmap_legend = FALSE,
        show_row_names =TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 12),
        row_names_max_width = unit(10, "cm"))
dev.off()
