
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

library(tidyr)
library(tidyverse)

dir.create("StatisticalAnalysis")

# load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# load metabolite names
metnames <- read.csv("MicrobiomeResults/CD_metNames.csv", header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
  
# perform statistical analysis for baseline and post
sigFeats <- character()
treatmentGroups <- c("Placebo","LMP","HMP")

### net secretion fluxes
fluxes=read.csv("MicrobiomeResults/CD_netSecretionFluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# remove samples with no data
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]
  
stats_res <- data.frame(matrix(ncol = 6, nrow = 0))

for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i] <- treatmentGroups[i]
  # extract the fluxes by treatment group
  fluxesRed <- fluxes
  samp <- metadata[metadata$"Treatment group" == treatmentGroups[i],]
  colsToKeep <- intersect(colnames(fluxes),rownames(samp))
  fluxesRed = fluxesRed[,colsToKeep, drop = FALSE]
  
  for (j in 1:nrow(fluxesRed)) {
    # run the test
    sdata <- data.frame(t(fluxesRed[j,]))
    for (k in 1:nrow(sdata)) {
      # get corresponding time point
      row_ind = which(rownames(metadata) == rownames(sdata)[k], arr.ind=TRUE)
      sdata[k,2] = metadata$"Time point"[row_ind]
    }
    colnames(sdata) <- c("Flux","TimePoint")
    wt <- wilcox.test(Flux ~ TimePoint, data = sdata)
    # fill in the results
    if  (!is.na(wt$p.value)) {
      stats_res[j,i] <- wt$p.value
    } else {
      stats_res[j,i] <- 1
    }
  }
}

### correct for multiple testing
for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i+3] <- paste(colnames(stats_res)[i],"after FDR")
  stats_res[,i+3] <- p.adjust(stats_res[,i], method = "fdr")
}

rownames(stats_res) <- metnames$"Metabolite name"
  
# export table
write.csv(stats_res, "StatisticalAnalysis/Net_Secretion_Statistical_Results.csv", row.names=TRUE)

#######
### analyze and plot model statistics
data=read.csv("MicrobiomeModels/ModelStatistics.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)
data <- t(data)

stats_res <- data.frame(matrix(ncol = length(treatmentGroups), nrow = length(rownames(data))))
rownames(stats_res) <- rownames(data)

# load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# remove samples with no data
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]

for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i] <- treatmentGroups[i]
  # extract the data by treatment group
  dataRed <- data
  samp <- metadata[metadata$"Treatment group" == treatmentGroups[i],]
  colsToKeep <- intersect(colnames(data),rownames(samp))
  dataRed = dataRed[,colsToKeep, drop = FALSE]
  
  for (j in 1:nrow(dataRed)) {
    # run the test
    sdata <- data.frame(dataRed[j,])
    for (k in 1:nrow(sdata)) {
      # get corresponding time point
      row_ind = which(rownames(metadata) == rownames(sdata)[k], arr.ind=TRUE)
      sdata[k,2] = metadata$"Time point"[row_ind]
    }
    colnames(sdata) <- c("Flux","TimePoint")
    wt <- wilcox.test(Flux ~ TimePoint, data = sdata)
    # fill in the results
    if  (!is.na(wt$p.value)) {
      stats_res[j,i] <- wt$p.value
    } else {
      stats_res[j,i] <- 1
    }
  }
}

### correct for multiple testing
for (i in length(treatmentGroups)) {
  stats_res[,i] <- p.adjust(stats_res[,i], method = "fdr")
}

# export table
write.csv(stats_res, "StatisticalAnalysis/Model_Features_Statistical_Results.csv", row.names=TRUE)

plist <- list()
for (i in 1:nrow(data)) {
  plotdata <- data.frame(
    Time.point=metadata[c("Time point")],
    Treatment.group=metadata[c("Treatment group")],
    Flux=data[i,])
  
  plist[[i]] = ggplot(data = plotdata, aes(x=fct_rev(Treatment.group), y=Flux, fill=Time.point)) +
    geom_boxplot(width=0.6, position = position_dodge2(preserve = "single")) +
    #facet_wrap(~Time.point, scale="free") +
    scale_y_continuous(limits = quantile(plotdata$Flux, c(0, 0.95))) +
    geom_jitter(color="black", size=0.2, alpha=0.6) +
    theme_classic() + coord_flip() +
    theme(plot.title = element_text(size=12),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    #ylim(0,NA) +
    ggtitle(rownames(data)[i]) +
    theme(axis.title.x=element_blank()) +
    xlab("") + theme(plot.title = element_text(hjust = 0.2,size = 12,
                                               color = "black",face = "bold"))
}

for (i in 1:nrow(data)) {
  imagePath =  paste("MetaboliteBoxplots","/Plot_",rownames(data)[i],".png",sep = "")
  png(imagePath, width = 6, height = 3, units = 'in', res = 300)
  plot(plist[[i]])
  dev.off()
}
