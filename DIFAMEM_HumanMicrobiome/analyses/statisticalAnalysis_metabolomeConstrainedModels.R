
currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

# load metadata
metadata = read.csv("input/SamplemetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# perform statistical analysis for baseline and post
treatmentGroups <- c("Placebo","LMP","HMP")

## fluxes for microbiome models constrained by fecal metabolome
fluxes=read.csv("MicrobiomeResults/MetabolomeConstrainedExchangeFluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# drop all rows that were below the flux cutoff
ns <- rep(NA,times=nrow(fluxes))
for (i in 1:nrow(fluxes)) {
  bools <- abs(fluxes[i,]) < 0.000001
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
fluxes = fluxes[ns == FALSE,]

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
rownames(stats_res) <- rownames(fluxes)
  
### correct for multiple testing
for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i+3] <- paste(colnames(stats_res)[i],"after FDR")
  stats_res[,i+3] <- p.adjust(stats_res[,i], method = "fdr")
}

# export table
exportfile <- paste("StatisticalAnalysis/FecalMetabolomeModels_Statistical_Results.csv", sep = "")
write.csv(stats_res, exportfile, row.names=TRUE)


#### Human model fluxes

# human reaction annotations
annotation = read.csv("MetabolomeResults/HumanReactionAnnotations.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

## fluxes for human models constrained by targeted metabolome
fluxes=read.csv("MetabolomeResults/Metabolome_1CM_ConstrainedFluxes.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)

# load metadata
metadata = read.csv("input/SamplemetadataPatients.csv",header=T,row.names=1, check.names = FALSE, stringsAsFactors=T)

# rename metadata row names
rownames(metadata) <- metadata$Metabolome_sample

# drop all rows that were below the flux cutoff
ns <- rep(NA,times=nrow(fluxes))
for (i in 1:nrow(fluxes)) {
  bools <- abs(fluxes[i,]) < 0.000001
  if (any(bools == FALSE)) {
    ns[i] = FALSE
  } else {
    ns[i] = TRUE
  }
}
fluxes = fluxes[ns == FALSE,]

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
rownames(stats_res) <- rownames(fluxes)

### correct for multiple testing
for (i in 1:length(treatmentGroups)) {
  colnames(stats_res)[i+3] <- paste(colnames(stats_res)[i],"after FDR")
  stats_res[,i+3] <- p.adjust(stats_res[,i], method = "fdr")
}

# export table
exportfile <- paste("StatisticalAnalysis/Metabolome_1CM_Models_Statistical_Results.csv", sep = "")
write.csv(stats_res, exportfile, row.names=TRUE)
