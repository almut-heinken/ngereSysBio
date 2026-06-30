
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

library(RColorBrewer)
library(openxlsx)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(cleanepi)
library(viridis)
library(dplyr)
library(ggthemes)

## Microbiome exchange fluxes predicted through FBA
# load fluxes
fluxes = read.csv("MicrobiomeResults/MetabolomeConstrainedExchangeFluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
annotation = read.csv("MicrobiomeResults/MetabolomeConstrainedExchangeFluxes_Annotation.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

## load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# remove samples with no data
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]

# read statistical analysis
stats_res = read.csv("StatisticalAnalysis/FecalMetabolomeModels_Statistical_Results.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

sigFeats <- vector()
for (i in 1:ncol(stats_res)) {
  for (j in 1:nrow(stats_res)) {
    if (stats_res[j,i] <0.05) {
      sigFeats <- union(sigFeats,rownames(stats_res)[j])
    }
  }
}
rowsToKeep <- intersect(rownames(fluxes),sigFeats)
fluxes = fluxes[sigFeats,, drop = FALSE]

rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]
colsToKeep <- intersect(colnames(fluxes),rownames(metadata))
fluxes = fluxes[,colsToKeep, drop = FALSE]

column_ha = HeatmapAnnotation(df = NULL, Treatment_Group=as.character(metadata$"Treatment group"),
                              Time_Point=as.character(metadata$"Time point"),
                              col=list(Treatment_Group =c("HMP"="red","LMP"="blue","Placebo"="green"),
                                       Time_Point =c("Basal"="yellow","Post"="purple")),
                              annotation_legend_param = list(labels_gp = gpar(fontsize = 10)))

mycol <- c(brewer.pal(9,"Spectral"))

scaledFluxes <- fluxes
scaledFluxes=t(scale(t(scaledFluxes)))
scaledFluxes <- as.matrix(scaledFluxes)
rownames(scaledFluxes) <- str_replace_all(rownames(scaledFluxes),'_IEX','')
rownames(scaledFluxes) <- str_replace_all(rownames(scaledFluxes),'\\[u\\]tr','')
rownames(scaledFluxes) <- str_replace_all(rownames(scaledFluxes),'pan','')
png("MicrobiomeResults/MetabolomeConstrainedExchanges.png", width = 10, height = 6, units = 'in', res = 300)
Heatmap(col=mycol,scaledFluxes,show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE,
        name = "Flux",top_annotation=column_ha,#right_annotation = row_ha,
        show_row_names = TRUE,row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8)
)
dev.off()

# Plot as boxplots
plotnames <- rownames(fluxes)

plist <- list()

for (i in 1:nrow(fluxes)) {
  # get microbe and metabolite name
  an_ind = which(rownames(annotation) == rownames(fluxes)[i], arr.ind=TRUE)
  microbe = str_replace(annotation$Species[an_ind],"_"," ")
  met = str_to_title(annotation$Metabolite[an_ind])
  
  plotdata <- as.data.frame(t(fluxes[i,]))
  for (m in 1:ncol(fluxes)) {
    row_ind = which(rownames(metadata) == colnames(fluxes)[m], arr.ind=TRUE)
    plotdata[m,2] <- metadata$"Treatment group"[row_ind]
    plotdata[m,3] <- metadata$"Time point"[row_ind]
  }
  colnames(plotdata) <- c("Flux","Treatment.group","Time.point")
  
  bp = ggplot(data = plotdata, aes(x=Treatment.group, y=Flux, fill=Time.point))
  plist[[i]] = bp +
    geom_hline(yintercept = 0, color = "grey", linewidth = 0.5) +
    geom_boxplot(width=0.6, position = position_dodge2(preserve = "single")) +
    #scale_y_continuous(limits=quantile(plotdata$Flux, c(0, 0.95))) +
    geom_jitter(color="black", size=0.2, alpha=0.6) +
    theme_classic() + coord_flip() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    ggtitle(paste(met,"exchange by",microbe)) + ylab("Flux (mmol/sample/day)") +
    theme(axis.title.y = element_text(hjust = 0.5,size = 10,
                                      color = "black")) +
    theme(axis.title.x = element_text(size = 10,
                                      color = "black")) +
    xlab("") + theme(plot.title = element_text(hjust = 0.5,size = 10,
                                               color = "black",face = "bold"))
}

for (i in 1:nrow(fluxes)) {
  imagePath =  paste("MetabolomeConstrainedModels/Boxplots","/Plot_Exch_FBA_",plotnames[i],".png",sep = "")
  png(imagePath, width = 5, height = 3, units = 'in', res = 300)
  plot(plist[[i]])
  dev.off()
}


## Human models constrained with one-carbon and TCA metabolome
# load fluxes
fluxes = read.csv("MetabolomeResults/Metabolome_1CM_ConstrainedFluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
annotation = read.csv("MetabolomeResults/HumanReactionAnnotations.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

## load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# rename metadata row names
rownames(metadata) <- metadata$Metabolome_sample

# plot for baseline and post
treatmentGroups <- c("Placebo","LMP","HMP")

# read statistical analysis
stats_res = read.csv("StatisticalAnalysis/Metabolome_1CM_Models_Statistical_Results.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

sigFeats <- vector()
for (i in 1:ncol(stats_res)) {
  for (j in 1:nrow(stats_res)) {
    if (stats_res[j,i] <0.05) {
      sigFeats <- union(sigFeats,rownames(stats_res)[j])
    }
  }
}
rowsToKeep <- intersect(rownames(fluxes),sigFeats)
fluxes = fluxes[sigFeats,, drop = FALSE]

# Plot as boxplots
plotnames <- rownames(fluxes)

dir.create("MetabolomeResults/1CM_Boxplots")

plist <- list()

for (i in 1:nrow(fluxes)) {
  # get reaction name
  an_ind = which(rownames(annotation) == rownames(fluxes)[i], arr.ind=TRUE)
  title = paste(annotation$ReactionName[an_ind],"\n",annotation$Subsystem[an_ind])
  
  plotdata <- as.data.frame(t(fluxes[i,]))
  for (m in 1:ncol(fluxes)) {
    row_ind = which(rownames(metadata) == colnames(fluxes)[m], arr.ind=TRUE)
    plotdata[m,2] <- metadata$"Treatment group"[row_ind]
    plotdata[m,3] <- metadata$"Time point"[row_ind]
  }
  colnames(plotdata) <- c("Flux","Treatment.group","Time.point")
  
  bp = ggplot(data = plotdata, aes(x=Treatment.group, y=Flux, fill=Time.point))
  plist[[i]] = bp +
    geom_hline(yintercept = 0, color = "grey", linewidth = 0.5) +
    geom_boxplot(width=0.6, position = position_dodge2(preserve = "single")) +
    scale_y_continuous(limits=quantile(plotdata$Flux, c(0, 0.95))) +
    geom_jitter(color="black", size=0.2, alpha=0.6) +
    theme_classic() + coord_flip() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    ggtitle(title) + ylab("Flux (mmol/g dry weight/hour)") +
    theme(axis.title.y = element_text(hjust = 0.5,size = 10,
                                      color = "black")) +
    theme(axis.title.x = element_text(size = 10,
                                      color = "black")) +
    xlab("") + theme(plot.title = element_text(hjust = 0.5,size = 10,
                                               color = "black",face = "bold"))
}

for (i in 1:nrow(fluxes)) {
  imagePath =  paste("MetabolomeResults/1CM_Boxplots/",plotnames[i],".png",sep = "")
  png(imagePath, width = 6, height = 3, units = 'in', res = 300)
  plot(plist[[i]])
  dev.off()
}

# plot reactions with higher and lower flux by subsystem
subs <- as.data.frame(matrix(nrow=length(rownames(fluxes)),ncol=1))
for (i in 1:nrow(fluxes)) {
  row_ind = which(rownames(annotation) == rownames(fluxes)[i], arr.ind=TRUE)
  subs[i,] <- as.character(annotation$Subsystem[row_ind])
}
unSubs <- sort(unique(subs$V1))

plist <- list()

for (t in 1:length(treatmentGroups)) {
  enrRxnsTable <- as.data.frame(matrix(0,nrow=length(unSubs),ncol=2))
  rownames(enrRxnsTable) <- unSubs
  colnames(enrRxnsTable) <- c("Lower flux post-intervention","Higher flux post-intervention ")
  # get for each reaction the significance and whether it is higher or lower post-intervention
  for (i in 1:nrow(fluxes)) {
    # find subsystem in results table
    sub_ind = which(rownames(enrRxnsTable) == subs$V1[i], arr.ind=TRUE)
    # find in stats table
    row_ind = which(rownames(stats_res) == rownames(fluxes)[i], arr.ind=TRUE)
    # extract the fluxes for treatment groups if significant
    if (stats_res[row_ind,t] < 0.05) {
      samp <- metadata[metadata$"Group" == paste(treatmentGroups[t],"_Basal",sep=""),]
      colsToKeep <- intersect(colnames(fluxes),rownames(samp))
      fluxes_basal = fluxes[i,colsToKeep, drop = FALSE]
      samp <- metadata[metadata$"Group" == paste(treatmentGroups[t],"_Post",sep=""),]
      colsToKeep <- intersect(colnames(fluxes),rownames(samp))
      fluxes_post = fluxes[i,colsToKeep, drop = FALSE]
      if  (mean(abs(as.numeric(fluxes_basal))) > mean(abs(as.numeric(fluxes_post)))) {
        enrRxnsTable[sub_ind,1] = enrRxnsTable[sub_ind,1]+1
      } else {
        enrRxnsTable[sub_ind,2] = enrRxnsTable[sub_ind,2]+1
      }
    }
  }
  # create bar plots
  # filter out rows with all zeros
  enrRxnsTable <- enrRxnsTable[rowSums(enrRxnsTable[])>0,]
  if (length(rownames(enrRxnsTable))>3) {
    # filter out exchange and transport
    enrRxnsTable <- enrRxnsTable %>% filter(!grepl('Exchange', rownames(enrRxnsTable)))
    enrRxnsTable <- enrRxnsTable %>% filter(!grepl('Transport', rownames(enrRxnsTable)))
    data_plot <- data.frame(matrix(nrow=length(rownames(enrRxnsTable))*2,ncol=0))
    cnt <- 1
    for (i in 1:length(rownames(enrRxnsTable))) {
      data_plot$Subsystem[cnt] <- rownames(enrRxnsTable)[i]
      data_plot$Scenario[cnt] <- "Upregulated"
      data_plot$Flux[cnt] <- enrRxnsTable$"Higher flux post-intervention"[i]
      cnt = cnt+1
    }
    for (i in 1:length(rownames(enrRxnsTable))) {
      data_plot$Subsystem[cnt] <- rownames(enrRxnsTable)[i]
      data_plot$Scenario[cnt] <- "Downregulated"
      data_plot$Flux[cnt] <- enrRxnsTable$"Lower flux post-intervention"[i]
      cnt = cnt+1
    }
    
    plist[[1]] = ggplot(data_plot, aes(y=Flux, x = fct_rev(Subsystem), fill=Subsystem)) + 
      geom_bar(position="dodge",stat = "identity") +
      facet_wrap(~Scenario) +
      theme_classic() + 
      coord_flip() + theme(legend.position="none") +
      ggtitle(paste(treatmentGroups[t],"group")) + 
      xlab("") + ylab("Reaction fluxes") +
      theme(plot.title = element_text(hjust = 0.5,size = 12,color = "black",face = "bold"))
    imagePath =  paste("MetabolomeResults","/Metabolome_1CM_",treatmentGroups[t],".png",sep = "")
    
    if (length(rownames(enrRxnsTable))>20) {
      h <- 6
      w <- 8
    } else {
      h <- 4
      w <- 5
    }
    png(imagePath, width = w, height = h, units = 'in', res = 300)
    plot(plist[[1]])
    dev.off()
  }
}
