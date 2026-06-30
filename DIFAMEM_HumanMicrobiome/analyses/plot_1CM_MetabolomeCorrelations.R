
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

library(RColorBrewer)
library(openxlsx)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(cleanepi)

## calculate and plot correlations between significant reaction fluxes and the metabolome
fluxes = read.csv("MicrobiomeResults/CD_netSecretionFluxes.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

## load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# load metabolite names
metnames <- read.csv("MicrobiomeResults/CD_metNames.csv", header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# excluding control
samp <- metadata[metadata$"Time point" != "Control",]
colsToKeep <- intersect(colnames(fluxes),rownames(samp))
fluxes = fluxes[,colsToKeep, drop = FALSE]
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]

# perform statistical analysis for baseline and post
sigFeats <- character()
timePoints <- c("Basal","Post")

stats_res <- data.frame(matrix(ncol = length(timePoints), nrow = 0))
colnames(stats_res) <- timePoints

for (i in 1:length(timePoints)) {
  # extract the data by time point
  dataRed <- fluxes
  samp <- metadata[metadata$"Time point" == timePoints[i] | metadata$"Time point" == "Control",]
  rowsToKeep <- intersect(rownames(metadata),rownames(samp))
  metadataRed = metadata[rowsToKeep,, drop = FALSE]
  colsToKeep <- intersect(colnames(fluxes),rownames(samp))
  dataRed = dataRed[,colsToKeep, drop = FALSE]
  
  for (j in 1:nrow(dataRed)) {
    # run the test
    sdata <- data.frame(t(dataRed[j,]))
    sdata$"Treatment group" <- paste(metadataRed$"Treatment group")
    colnames(sdata) <- c("Flux","TreatmentGroup")
    one.way <- aov(Flux ~ TreatmentGroup, data = sdata)
    # fill in the results
    if  (!is.na(summary(one.way)[[1]][["Pr(>F)"]][1])) {
      stats_res[j,i] <- summary(one.way)[[1]][["Pr(>F)"]][1]
    } else {
      stats_res[j,i] <- 1
    }
  }
}

### correct for multiple testing
for (i in 1:length(timePoints)) {
  stats_res[,i] <- p.adjust(stats_res[,i], method = "fdr")
}

### read in metabolome
metabolome = read.csv("input/Metabolome_1CM.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

for (i in 1:ncol(metabolome)) {
  row_ind = which(metadata$Metabolome_sample == colnames(metabolome)[i], arr.ind=TRUE)
  if (length(row_ind)==1) {
  colnames(metabolome)[i] <- rownames(metadata)[row_ind]
  }
}

# remove metabolome samples with no corresponding flux sample and vice versa
colsToKeep <- intersect(colnames(metabolome), colnames(fluxes))
metabolome = metabolome[,colsToKeep, drop = FALSE]
metabolome <- droplevels(metabolome)
colsToKeep <- intersect(colnames(fluxes), colnames(metabolome))
fluxes = fluxes[,colsToKeep, drop = FALSE]
fluxes <- droplevels(fluxes)

## remove fluxes that are constant after removing some samples
fluxes = t(fluxes)
fluxes <- remove_constants(fluxes, cutoff = 1)
fluxes = t(fluxes)

rowsToKeep <- intersect(rownames(metnames), rownames(fluxes))
metnames = metnames[rowsToKeep,, drop = FALSE]
metnames <- droplevels(metnames)

## create the correlations matrix
correlations <- data.frame(matrix(ncol = length(rownames(metabolome)),nrow = length(rownames(fluxes))))
rownames(correlations) <- rownames(fluxes)
colnames(correlations) <- rownames(metabolome)
pvalues <- correlations

for (i in 1:nrow(fluxes)) {
  data <- data.frame(matrix(nrow = length(colnames(metabolome)), ncol = 2))
  data[,1] <- fluxes[i,]
  for (j in 1:nrow(metabolome)) {
    for (k in 1:ncol(fluxes)) {
    # get corresponding metabolome sample
    col_ind = which(colnames(metabolome) == colnames(fluxes)[k], arr.ind=TRUE)
    data[k,2] = metabolome[j,col_ind]
    }
    corr <- cor.test(as.numeric(data[,1]), data[,2], method = 'spearman')
    correlations[i,j] <- corr$estimate
    pvalues[i,j] <- corr$p.value
  }
}

# adjust for multiple testing
pvaluesFDR <- pvalues
for (i in 1:ncol(pvalues)) {
  pvaluesFDR[,i] <- p.adjust(pvalues[,i], method = "fdr")
}
write.csv(correlations, "MicrobiomeResults/Flux_Metabolome_Correlations.csv", row.names=TRUE)
write.csv(pvalues, "MicrobiomeResults/Flux_Metabolome_pValues.csv", row.names=TRUE)
write.csv(pvaluesFDR, "MicrobiomeResults/Flux_Metabolome_pValues_afterFDR.csv", row.names=TRUE)

# create scatter plots of significant correlations
dir.create("MicrobiomeMetabolomeCorrelations")
plist <- list()
cnt <- 1
for (i in 1:nrow(fluxes)) {
  for (j in 1:nrow(metabolome)) {
    if (pvaluesFDR[i,j] < 0.05) {
      data <- data.frame(matrix(nrow = length(colnames(metabolome)), ncol = 2))
      data[,1] <- fluxes[i,]
      rownames(data) <- colnames(fluxes)
      for (k in 1:ncol(fluxes)) {
        # get corresponding metabolome sample
        col_ind = which(colnames(metabolome) == colnames(fluxes)[k], arr.ind=TRUE)
        data[k,2] = metabolome[j,col_ind]
      }
      colnames(data) <- c("Flux","Metabolite")
      data$Treatment_group <- metadata$"Treatment group"
      data$Time_point <- metadata$"Time point"
      
      # make plot
      plist[[cnt]] = ggplot(data = data, aes(x=Flux, y=Metabolite)) +
        geom_point(aes(color=Treatment_group,shape=Time_point)) +
        theme_classic() +
        ggtitle("Microbial metabolite flux-metabolome correlation",paste("Spearman correlation", round(correlations[i,j],2),",","p-value", round(pvaluesFDR[i,j],3))) +
        theme(axis.title.y = element_text(hjust = 0.5,size = 12,
                                          color = "black")) +
        theme(axis.title.x = element_text(size = 12,
                                          color = "black")) +
        xlab(paste(metnames$"Metabolite name"[i],"(mmol/person/day)")) + ylab(paste(rownames(metabolome)[j],"(µmol/L)")) + 
        theme(plot.title = element_text(hjust = 0.5,size = 14,color = "black",face = "bold"))
      cnt = cnt + 1
    }
  }
}

cnt <- 1
for (i in 1:nrow(fluxes)) {
  for (j in 1:nrow(metabolome)) {
    if (pvaluesFDR[i,j] < 0.05) {
      metname <- str_split_1(rownames(metabolome)[j], " ")
      imagePath =  paste("MicrobiomeMetabolomeCorrelations/",rownames(fluxes)[i],"_",metname[1],".png",sep = "")
      png(imagePath, width = 6, height = 4, units = 'in', res = 300)
      plot(plist[[cnt]])
      dev.off()
      cnt = cnt + 1
    }
  }
}

# drop all rows that were not significant
ns <- rep(NA,times=nrow(pvaluesFDR))
for (i in 1:nrow(pvaluesFDR)) {
  bools <- pvaluesFDR[i,] >= 0.05
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[ns == FALSE,]
pvaluesFDR = pvaluesFDR[ns == FALSE,]

# drop all columns that were not significant
ns <- rep(NA,times=ncol(pvaluesFDR))
for (i in 1:ncol(pvaluesFDR)) {
  bools <- pvaluesFDR[,i] >= 0.05
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[,ns == FALSE]
pvaluesFDR = pvaluesFDR[,ns == FALSE]

## plot correlations
rowsToKeep <- intersect(rownames(metnames), rownames(correlations))
metnames = metnames[rowsToKeep,, drop = FALSE]
metnames <- droplevels(metnames)
rownames(correlations) <- metnames$"Metabolite name"

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
palette = getPalette(length(levels(metnames$Subsystem)))
SubsystemAnn.cols <- c(palette)
names(SubsystemAnn.cols) <- levels(metnames$Subsystem)

row_ha = rowAnnotation(df = NULL, Subsystem=as.character(metnames$Subsystem),col=list(Subsystem=SubsystemAnn.cols),
                              annotation_legend_param = list(labels_gp = gpar(fontsize = 10)))

mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(9,"Reds")[1:9])

png("MicrobiomeResults/Flux_Metabolome_Correlations.png", width = 5, height = 9.5, units = 'in', res = 300)
t <- data.matrix(correlations)
ht <- Heatmap(t,col=mycol,show_column_names = TRUE, cluster_columns = TRUE, cluster_rows = TRUE,
        name = "Correlations",row_names_max_width = unit(10, "cm"), 
        show_row_names = TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(at = seq(-0.6, 0.6))
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
