
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

library(tidyr)
library(tidyverse)
library(viridis)

## load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

### read in data
data = read.csv("input/Metabolome_1CM.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# rename metadata row names
rownames(metadata) <- metadata$Metabolome_sample

treatmentGroups <- c("Placebo","LMP","HMP")

stats_res <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i] <- treatmentGroups[i]
  # extract the data by treatment group
  dataRed <- data
  samp <- metadata[metadata$"Treatment group" == treatmentGroups[i],]
  colsToKeep <- intersect(colnames(data),rownames(samp))
  dataRed = dataRed[,colsToKeep, drop = FALSE]
  
  for (j in 1:nrow(dataRed)) {
    # run the test
    sdata <- data.frame(t(dataRed[j,]))
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
for (i in 1:3) {
  stats_res[,i] <- p.adjust(stats_res[,i], method = "fdr")
}
rownames(stats_res) <- rownames(data)

# export table
write.csv(stats_res, "StatisticalAnalysis/MetabolomicData_Statistical_Results.csv", row.names=TRUE)

## plot data

dir.create("Metabolome_Boxplots")

plist <- list()

for (i in 1:nrow(data)) {
  plotdata <- data.frame(
    Flux=t(data[i,]))
  for (k in 1:nrow(plotdata)) {
    # get corresponding time point
    row_ind = which(rownames(metadata) == rownames(plotdata)[k], arr.ind=TRUE)
    plotdata[k,2] = metadata$"Time point"[row_ind]
    plotdata[k,3] = metadata$"Treatment group"[row_ind]
  }
  colnames(plotdata) <- c("Flux","Time.point","Treatment.group")
  
  plist[[i]] = ggplot(data = plotdata, aes(x=fct_rev(Treatment.group), y=Flux, fill=Time.point)) +
    geom_boxplot(width=0.6, position = position_dodge2(preserve = "single")) +
    #facet_wrap(~Time.point, scale="free") +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.2, alpha=0.6) +
    theme_classic() + coord_flip() +
    theme(plot.title = element_text(size=12),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    #ylim(0,NA) +
    ggtitle(rownames(data)[i]) + ylab("Concentration (µmol/L)") +
    theme(axis.title.y = element_text(hjust = 0.5,size = 12,
                                      color = "black")) +
    theme(axis.title.x = element_text(size = 12,
                                      color = "black")) +
    xlab("") + theme(plot.title = element_text(hjust = 0.2,size = 12,
                                               color = "black",face = "bold"))
}

for (i in 1:nrow(data)) {
  imagePath =  paste("Metabolome_Boxplots/",rownames(data)[i],".png",sep = "")
  png(imagePath, width = 6, height = 3, units = 'in', res = 300)
  plot(plist[[i]])
  dev.off()
}

