
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

# Libraries
library(ggplot2)
library(ggalign)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggthemes)

# Load predicted microbial metabolite production dataset
fluxes <- read.csv("MicrobiomeResults/CD_netSecretionFluxes.csv", header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# load metadata
metadata = read.csv("input/SampleMetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# remove samples with no data
rowsToKeep <- intersect(rownames(metadata),colnames(fluxes))
metadata = metadata[rowsToKeep,, drop = FALSE]

dir.create("MetaboliteBoxplots")

# load metabolite names
metnames <- read.csv(paste("MicrobiomeResults/","CD_metNames.csv", sep = ""), header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# plot metabolites
fluxes <- t(fluxes)

# Plot
plotnames <- colnames(fluxes)

plist <- list()

for (i in 1:ncol(fluxes)) {
  plotdata <- data.frame(
    Time.point=metadata[c("Time point")],
    Treatment.group=metadata[c("Treatment group")],
    Flux=fluxes[,i])
  
  plist[[i]] = ggplot(data = plotdata, aes(x=Treatment.group, y=Flux, fill=Time.point)) +
    geom_boxplot(width=0.6, position = position_dodge2(preserve = "single")) +
    scale_y_continuous(limits = quantile(plotdata$Flux, c(0, 0.95))) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.2, alpha=0.6) +
    theme_classic() + coord_flip() +
    theme(plot.title = element_text(size=12),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    ggtitle(metnames[i,'Metabolite name']) + ylab("Flux (mmol/sample/day)") +
    theme(axis.title.y = element_text(hjust = 0.5,size = 12,
                                      color = "black")) +
    theme(axis.title.x = element_text(size = 12,
                                      color = "black")) +
    xlab("") + theme(plot.title = element_text(hjust = 0.5,size = 12,
                                               color = "black",face = "bold"))
}

for (i in 1:ncol(fluxes)) {
  imagePath =  paste("MetaboliteBoxplots","/Plot_",plotnames[i],".png",sep = "")
  png(imagePath, width = 5, height = 3, units = 'in', res = 300)
  plot(plist[[i]])
  dev.off()
}
