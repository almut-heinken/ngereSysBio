
## plot IEM data ####

library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

### Fluxes KOs ####
data=read.csv("Results/Altered_KO_Fluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# scale to data to rows
#data=t(scale(t(data)))

#### read reaction subsystems
annotations=read.csv("Results/Altered_KO_Subsystems.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)[c("Subsystem")]

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(annotations$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(annotations$Subsystem)

column_ha = HeatmapAnnotation(df = NULL, Subsystem=as.character(annotations$Subsystem),col=list(Subsystem=SubsystemAnn.cols),
                       annotation_legend_param = list(title_gp = gpar(fontsize = 12), 
                                                      labels_gp = gpar(fontsize = 12)))

mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(9,"YlOrRd")[2:4],brewer.pal(9,"Reds")[6:9])

png("Figure_2.png", width = 12, height = 8, units = 'in', res = 300)
t <- data.matrix(t(data))
ht <- Heatmap(t,col=mycol,show_column_names = TRUE, cluster_columns = TRUE, cluster_rows = TRUE,
        name = "Flux", top_annotation=column_ha,border = TRUE,
        show_row_names = TRUE, row_names_max_width = unit(10, "cm"),
        column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),legend_gp = gpar(fontsize = 12))
        )
draw(ht,heatmap_legend_side = "right", annotation_legend_side = "bottom")
dev.off()
