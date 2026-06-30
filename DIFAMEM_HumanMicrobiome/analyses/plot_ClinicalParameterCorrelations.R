
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

library(RColorBrewer)
library(openxlsx)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(cleanepi)

## calculate and plot correlations between significant reaction fluxes and the clinicalData
fluxes = read.csv("MicrobiomeResults/CD_netSecretionFluxes.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

## load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# load metabolite names
metnames <- read.csv("MicrobiomeResults/CD_metNames.csv", header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

### read in clinical parameters
clinicalData = read.csv("input/ClinicalParameters.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# define features to analyze
feats <- c("Bilirubin","Gamma-glutamyl transferase","Aspartate aminotransferase","Glucose","Total Cholesterol","Alkaline Phosphatase","Triglycerides","A2-Macroglobulin","Apolipoprotein A-1","Haptoglobin","SPT Wheal area (mm2)","sIgE Pru p 3 sIgE (kU/L)","DBPCFC [Pru p 3] (μg/mL)")
colsToKeep <- intersect(colnames(clinicalData), feats)
clinicalData = clinicalData[,colsToKeep, drop = FALSE]
clinicalData <- data.frame(t(clinicalData))

# remove clinicalData samples with no corresponding flux sample and vice versa
colsToKeep <- intersect(colnames(clinicalData), colnames(fluxes))
clinicalData = clinicalData[,colsToKeep, drop = FALSE]
clinicalData <- droplevels(clinicalData)
colsToKeep <- intersect(colnames(fluxes), colnames(clinicalData))
fluxes = fluxes[,colsToKeep, drop = FALSE]
fluxes <- droplevels(fluxes)

## remove fluxes that are constant after removing some samples
fluxes = t(fluxes)
fluxes <- remove_constants(fluxes, cutoff = 1)
fluxes = t(fluxes)

# remove metadata entries not in data
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]

rowsToKeep <- intersect(rownames(metnames), rownames(fluxes))
metnames = metnames[rowsToKeep,, drop = FALSE]
metnames <- droplevels(metnames)

## create the correlations matrix
correlations <- data.frame(matrix(ncol = length(rownames(clinicalData)),nrow = length(rownames(fluxes))))
rownames(correlations) <- rownames(fluxes)
colnames(correlations) <- rownames(clinicalData)
pvalues <- correlations

for (i in 1:nrow(fluxes)) {
  data <- data.frame(matrix(nrow = length(colnames(clinicalData)), ncol = 2))
  data[,1] <- fluxes[i,]
  for (j in 1:nrow(clinicalData)) {
    for (k in 1:ncol(fluxes)) {
      # get corresponding clinicalData sample
      col_ind = which(colnames(clinicalData) == colnames(fluxes)[k], arr.ind=TRUE)
      data[k,2] = clinicalData[j,col_ind]
    }
    data <- na.omit(data)
    corr <- cor.test(as.numeric(data[,1]), data[,2], method = 'spearman')
    if (!is.na(corr$p.value)) {
      correlations[i,j] <- corr$estimate
      pvalues[i,j] <- corr$p.value }
    else {
      correlations[i,j] <- 0
      pvalues[i,j] <- 1
    }
  }
}

# adjust for multiple testing
pvaluesFDR <- pvalues
for (i in 1:ncol(pvalues)) {
  pvaluesFDR[,i] <- p.adjust(pvalues[,i], method = "fdr")
}


write.csv(correlations, "MicrobiomeResults/Flux_ClinicalData_Correlations.csv", row.names=TRUE)
write.csv(pvalues, "MicrobiomeResults/Flux_ClinicalData_pValues.csv", row.names=TRUE)
write.csv(pvaluesFDR, "MicrobiomeResults/Flux_ClinicalData_pValues_afterFDR.csv", row.names=TRUE)

# create scatter plots of significant correlations
dir.create("MicrobiomeResults/MicrobiomeClinicalDataCorrelations")
plist <- list()
cnt <- 1
for (i in 1:nrow(fluxes)) {
  for (j in 1:nrow(clinicalData)) {
    if (pvalues[i,j] < 0.05 & abs(correlations[i,j]) > 0.35) {
      data <- data.frame(matrix(nrow = length(colnames(clinicalData)), ncol = 2))
      data[,1] <- fluxes[i,]
      rownames(data) <- colnames(fluxes)
      for (k in 1:ncol(fluxes)) {
        # get corresponding clinical data sample
        col_ind = which(colnames(clinicalData) == colnames(fluxes)[k], arr.ind=TRUE)
        data[k,2] = clinicalData[j,col_ind]
      }
      colnames(data) <- c("Flux","Metabolite")
      data$Treatment_group <- metadata$"Treatment group"
      data$Time_point <- metadata$"Time point"
      
      # make plot
      plist[[cnt]] = ggplot(data = data, aes(x=Flux, y=Metabolite)) +
        geom_point(aes(color=Treatment_group,shape=Time_point)) +
        theme_classic() +
        ggtitle(paste("Spearman correlation", round(correlations[i,j],2),",","p-value after adjusting for false discovery rate", round(pvaluesFDR[i,j],3))) +
        theme(axis.title.y = element_text(hjust = 0.5,size = 10,
                                          color = "black")) +
        theme(axis.title.x = element_text(size = 10,
                                          color = "black")) +
        xlab(paste(metnames$"Metabolite name"[i],"(mmol/person/day)")) + ylab(paste(rownames(clinicalData)[j])) + 
        theme(plot.title = element_text(hjust = 0.2,size = 10,color = "black",face = "bold"))
      cnt = cnt + 1
    }
  }
}

cnt <- 1
for (i in 1:nrow(fluxes)) {
  for (j in 1:nrow(clinicalData)) {
    if (pvalues[i,j] < 0.05 & abs(correlations[i,j]) > 0.35) {
      metname <- str_split_1(rownames(clinicalData)[j], " ")
      imagePath =  paste("MicrobiomeResults/MicrobiomeClinicalDataCorrelations/",rownames(fluxes)[i],"_",metname[1],".png",sep = "")
      png(imagePath, width = 6, height = 4, units = 'in', res = 300)
      plot(plist[[cnt]])
      dev.off()
      cnt = cnt + 1
    }
  }
}

# drop all rows that were not significant
ns <- rep(NA,times=nrow(pvalues))
for (i in 1:nrow(pvalues)) {
  bools <- pvalues[i,] >= 0.05
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[ns == FALSE,]
pvalues = pvalues[ns == FALSE,]

ns <- rep(NA,times=nrow(pvalues))
for (i in 1:nrow(pvalues)) {
  bools <- abs(correlations[i,]) < 0.35
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[ns == FALSE,]
pvalues = pvalues[ns == FALSE,]

# drop all columns that were not significant
ns <- rep(NA,times=ncol(pvalues))
for (i in 1:ncol(pvalues)) {
  bools <- pvalues[,i] >= 0.05
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[,ns == FALSE]
pvalues = pvalues[,ns == FALSE]

ns <- rep(NA,times=ncol(pvalues))
for (i in 1:ncol(pvalues)) {
  bools <- abs(correlations[,i]) < 0.35
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
correlations= correlations[,ns == FALSE]
pvalues = pvalues[,ns == FALSE]

## plot correlations
rowsToKeep <- intersect(rownames(metnames), rownames(correlations))
metnames = metnames[rowsToKeep,, drop = FALSE]
metnames <- droplevels(metnames)
rownames(correlations) <- metnames$"Metabolite name"

mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(6,"Reds")[1:6])

png("MicrobiomeResults/MicrobialMetabolite_ClinicalData_Correlations.png", width = 6, height = 10, units = 'in', res = 300)
t <- data.matrix(correlations)
ht <- Heatmap(t,col=mycol,show_column_names = TRUE, cluster_columns = TRUE, cluster_rows = TRUE,
        name = "Correlations",column_names_max_height =  unit(10, "cm"),
        show_row_names = TRUE,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10)
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
