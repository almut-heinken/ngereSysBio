
## Plot the reaction presence that differed significantly between time points

library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")
dir.create("R_plots")

########## infants only ########
metadata = read.csv("inputFiles/Sample_metadata.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)[c("Group","Time point")]

##### reaction presence ####
data=read.csv("Exported_results/ByTimePoint_ReactionPresence.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

#### read reaction subsystems
annotations=read.csv("Exported_results/ByTimePoint_ReactionPresenceSubsystems.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)[c("Subsystem")]
  
png("R_plots/Figure_2a.png", width = 16, height = 8, units = 'in', res = 300)
  
# first time point
metadataReduced <- metadata
rowsToKeep <- which(as.character(metadataReduced[,2])=="5 days")
metadataReduced = metadataReduced[rowsToKeep,, drop = FALSE]
  
dataReduced <- data
colsToKeep <- intersect(colnames(data),rownames(metadataReduced))
dataReduced = dataReduced[,colsToKeep, drop = FALSE]
  
metadataReduced$Group <- as.factor(metadataReduced$Group)
column_ha = HeatmapAnnotation(df = NULL, Group=as.character(metadataReduced$Group),
                                col=list(Group = c("VD" = "blue","CSD" = "red")))
  
mycol <- c("white",brewer.pal(9,"Reds")[1:9])
  
t <- data.matrix(dataReduced)
h1 <-Heatmap(t,col=mycol,show_column_names = FALSE,column_title = "Five days",
               top_annotation=column_ha,name = "Five days",show_row_names = FALSE)
  
# second time point
metadataReduced <- metadata
rowsToKeep <- which(as.character(metadataReduced[,2])=="1 month")
metadataReduced = metadataReduced[rowsToKeep,, drop = FALSE]
  
dataReduced <- data
colsToKeep <- intersect(colnames(data),rownames(metadataReduced))
dataReduced = dataReduced[,colsToKeep, drop = FALSE]
  
metadataReduced$Group <- as.factor(metadataReduced$Group)
column_ha = HeatmapAnnotation(df = NULL, Group=as.character(metadataReduced$Group),
                                col=list(Group = c("VD" = "blue","CSD" = "red")))
  
mycol <- c("white",brewer.pal(9,"Oranges")[1:9])
  
t <- data.matrix(dataReduced)
h2 <-Heatmap(t,col=mycol,show_column_names = FALSE,column_title = "One month",
               top_annotation=column_ha,name = "One month",show_row_names = FALSE)
  
# third time point
metadataReduced <- metadata
rowsToKeep <- which(as.character(metadataReduced[,2])=="6 months")
metadataReduced = metadataReduced[rowsToKeep,, drop = FALSE]
  
dataReduced <- data
colsToKeep <- intersect(colnames(data),rownames(metadataReduced))
dataReduced = dataReduced[,colsToKeep, drop = FALSE]
  
metadataReduced$Group <- as.factor(metadataReduced$Group)
column_ha = HeatmapAnnotation(df = NULL, Group=as.character(metadataReduced$Group),
                                col=list(Group = c("VD" = "blue","CSD" = "red")))
  
mycol <- c("white",brewer.pal(9,"Purples")[1:9])
  
t <- data.matrix(dataReduced)
h3 <-Heatmap(t,col=mycol,show_column_names = FALSE,column_title = "Six months",
               top_annotation=column_ha,name = "Six months",show_row_names = FALSE)
  
# fourth time point
metadataReduced <- metadata
rowsToKeep <- which(as.character(metadataReduced[,2])=="1 year")
metadataReduced = metadataReduced[rowsToKeep,, drop = FALSE]
  
dataReduced <- data
colsToKeep <- intersect(colnames(data),rownames(metadataReduced))
dataReduced = dataReduced[,colsToKeep, drop = FALSE]
  
metadataReduced$Group <- as.factor(metadataReduced$Group)
column_ha = HeatmapAnnotation(df = NULL, Group=as.character(metadataReduced$Group),
                                col=list(Group = c("VD" = "blue","CSD" = "red")))
  
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
palette = getPalette(length(levels(annotations$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(annotations$Subsystem)
  
row_ha = rowAnnotation(df = NULL, Subsystem=as.character(annotations$Subsystem),col=list(Subsystem =SubsystemAnn.cols))
  
mycol <- c("white",brewer.pal(9,"Blues")[1:9])
  
t <- data.matrix(dataReduced)
h4 <-Heatmap(t,col=mycol,show_column_names = FALSE,column_title = "One year",
               top_annotation=column_ha,name = "One year",
               right_annotation=row_ha,show_row_names = FALSE)
  
ht_list <- h1+h2+h3+h4
  
draw(ht_list,ht_gap = unit(0.2,"cm"))
  
dev.off()
