
### plot significant correlations between reaction fluxes and lipidomics

library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)

### Lipidomics
#### fluxes ####
data=read.delim("Correlations/HighCorrelations_Lipidomics_Fluxes.txt",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
data <- data.matrix(data)

#### read reaction subsystems
annotations=read.delim("Correlations/HighCorrelations_Lipidomics_Fluxes_Annotations.txt",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)[c("Subsystem")]

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(annotations$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(annotations$Subsystem)

row_ha = rowAnnotation(df = NULL, Subsystem=as.character(annotations$Subsystem),col=list(Subsystem=SubsystemAnn.cols))

mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(9,"Reds")[1:9])

png("Figure_6a.png", width = 16, height = 20, units = 'in', res = 300)
t <- t(data.matrix(data))
Heatmap(t,col=mycol,show_column_names = FALSE, show_row_names = TRUE,
        cluster_columns = TRUE, cluster_rows = TRUE,
        name = "Correlation", right_annotation=row_ha,row_names_gp = gpar(fontsize = 10)
)
dev.off()
