
## plot species-metabolite contributions summarized on genus level

library(RColorBrewer)
library(ComplexHeatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")

####### have metabolite contributions in one plot ##########
##### different taxon levels ######

taxa=c("Genus");

### short-chain fatty acids ###

iw=c(9);
ih=c(16);

for (i in 1:length(taxa)) {
  
  ####
  fpath = paste("ContributionsByMetabolite/","ac_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  if (i>2){
    apath = paste("ContributionsByMetabolite/","ac_",taxa[i],"_annotations.csv",sep="")
    annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  }
  
  mycol <- c("white",brewer.pal(9,"Reds")[3:9])
  t <- data.matrix(data)
  
  if (i>2){
    row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
    ht1 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(10, "cm"),column_names_max_height = unit(10, "cm"),name = "Acetate flux",
                  row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,show_column_names = TRUE,column_names_side="top",
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Acetate")
    
  } else {
    
    ht1 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(10, "cm"),column_names_max_height = unit(10, "cm"),name = "Acetate flux",
                  row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",show_column_names = TRUE,column_names_side="top",
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Acetate")
  }
  
  ###
  fpath = paste("ContributionsByMetabolite/","ppa_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","ppa_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]

  mycol <- c("white",brewer.pal(9,"Oranges")[3:9])
  t <- data.matrix(data)
  

  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht2 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(10, "cm"),column_names_max_height = unit(10, "cm"),name = "Propionate flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Propionate")
 
  ### 
  fpath = paste("ContributionsByMetabolite/","but_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","but_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Purples")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht3 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(10, "cm"),column_names_max_height = unit(10, "cm"),name = "Butyrate flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Butyrate")

  ### 
  fpath = paste("ContributionsByMetabolite/","lac_L_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","lac_L_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Greens")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht4 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(10, "cm"),column_names_max_height = unit(10, "cm"),name = "L-lactate flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "L-lactate")

  ht_list = ht1 %v% ht2 %v% ht3 %v% ht4
  ipath = paste("R_plots/","Figure_4b.png",sep="")
  png(ipath, width = iw[i], height = ih[i], units = 'in', res = 300)
  draw(ht_list,ht_gap = unit(0.5,"cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

### vitamins ###

iw=c(8);
ih=c(16);

for (i in 1:length(taxa)) {
  
  ####
  fpath = paste("ContributionsByMetabolite/","adocbl_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)

  apath = paste("ContributionsByMetabolite/","adocbl_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Reds")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht1 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "AdoCbl flux",
                  row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,show_column_names = TRUE,column_names_side="top",
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "AdoCbl")

  ###
  fpath = paste("ContributionsByMetabolite/","btn_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","btn_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]

  mycol <- c("white",brewer.pal(9,"Oranges")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht2 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Biotin flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Biotin")
  
  ### 
  fpath = paste("ContributionsByMetabolite/","fol_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","fol_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
 
  mycol <- c("white",brewer.pal(9,"YlOrBr")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht3 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Folate flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Folate")

  ###
  fpath = paste("ContributionsByMetabolite/","nac_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","nac_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Purples")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht4 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Nicotinic acid flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Nicotinic acid")
  
  ####
  fpath = paste("ContributionsByMetabolite/","pnto_R_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","pnto_R_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Blues")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht5 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Pantothenate flux",
                  row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,show_column_names = TRUE,column_names_side="top",
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Pantothenate")

  ###
  fpath = paste("ContributionsByMetabolite/","pydx_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","pydx_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]

  mycol <- c("white",brewer.pal(9,"GnBu")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht6 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Pyridoxal flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Pyridoxal")

  ### 
  fpath = paste("ContributionsByMetabolite/","ribflv_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","ribflv_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
  
  mycol <- c("white",brewer.pal(9,"Greens")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht7 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Riboflavin flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Riboflavin")

  ###
  fpath = paste("ContributionsByMetabolite/","thm_",taxa[i],".csv",sep="")
  data = read.csv(fpath,header=T,row.names=1, check.names = FALSE)
  
  apath = paste("ContributionsByMetabolite/","thm_",taxa[i],"_annotations.csv",sep="")
  annotations=read.csv(apath,header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Phylum")]
 
  mycol <- c("white",brewer.pal(9,"Greys")[3:9])
  t <- data.matrix(data)
  
  row_ha = rowAnnotation(df = NULL, Phylum=as.character(annotations$Phylum),annotation_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),show_annotation_name = FALSE,
                           col=list(Phylum =c("Actinobacteria"="#E41A1C","Bacteroidetes"="#3A86A5","Firmicutes"="#658E67","Fusobacteria"="#CB6651","Proteobacteria"="#FFD422","Synergistetes"="#B47229","Verrucomicrobia"="#F781BF")))
    
  ht8 <-Heatmap(t,col=mycol,show_row_names = TRUE,show_heatmap_legend = TRUE,
                  row_names_max_width = unit(12, "cm"),column_names_max_height = unit(10, "cm"),name = "Thiamin flux",
                  show_column_names = FALSE,row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),cluster_columns=FALSE,cluster_rows=FALSE,
                  border = "grey",right_annotation=row_ha,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),row_title = "Thiamin")

  
  ht_list = ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8
  ipath = paste("R_plots/","Figure_4c.png",sep="")
  png(ipath, width = iw[i], height = ih[i], units = 'in', res = 300)
  draw(ht_list,ht_gap = unit(0.5,"cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
}
