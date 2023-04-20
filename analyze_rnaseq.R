
library(stringr)
library(Biostrings)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(gplots)
library(reshape2)
library(ggrepel)
library(sqldf)
library(stringr)

dat_gpcr <- read.csv("../gpcr_final_list converted3.csv")
dat_gpcr$highlight[is.na(dat_gpcr$highlight)] <- FALSE


dat_gpcr_red <- dat_gpcr[dat_gpcr$includefig1,]

dat_gpcr_red$genetype_red <- dat_gpcr_red$genetype
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"sr")] <- "Chemo receptors"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"str")] <- "Chemo receptors"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"DUF11")] <- "DUF11*"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"DUF6")] <- "DUF6*"
dat_gpcr_red$genetype_red[dat_gpcr_red$genetype_red %in% c("Light receptors")] <- "Neuropeptide receptors"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"Frizzled")] <- "Frizzled/Taste"
unique(dat_gpcr_red$genetype_red)


n <- 36
qual_col_pals_ct = brewer.pal.info[brewer.pal.info$category == 'qual',]
large_col_vector_ct = unlist(mapply(brewer.pal, qual_col_pals_ct$maxcolors, rownames(qual_col_pals_ct)))



##############################################################################
####################### Expression analysis - plotting ####################### 
##############################################################################


############## Compare expression in the L4
l4tpm <- qs::qread("L4.all.TPM.raw.qs")

#### Find agreeing set of genes
l4tpm_datpcr <- dat_gpcr_red[dat_gpcr_red$highlight=="TRUE",]
l4tpm <- l4tpm[rownames(l4tpm) %in% l4tpm_datpcr$wbid,]
#mean(dat_gpcr_red$wbid %in% rownames(l4tpm)) #8% lost
rownames(l4tpm_datpcr) <- l4tpm_datpcr$wbid
l4tpm_datpcr <- l4tpm_datpcr[rownames(l4tpm),]


#normalize each gene; use pseudo log scale
for(i in 1:nrow(l4tpm)){
  l4tpm[i,] <- 1e6*l4tpm[i,]/sum(l4tpm[i,])
}
l4tpm <- log10(1+l4tpm)


if(FALSE){
  ## Heatmap of ALL genes
  heatmap.2(l4tpm, scale = "none", col = bluered(100), 
            trace = "none", density.info = "none")
}


##### order cells by amount, then primarily by annotated cell type or gene type
l4tpm_celltype <- read.csv("../ce_celltypes_annot.csv")
rownames(l4tpm_celltype) <- l4tpm_celltype$Cell
l4tpm_celltype$SimpleDescription <- factor(
  l4tpm_celltype$SimpleDescription, 
  levels=c("other", "interneuron","motor neuron","sensory neuron"))
l4tpm_celltype <- l4tpm_celltype[order(l4tpm_celltype[colnames(l4tpm),]$SimpleDescription, colSums(l4tpm)),]
l4tpm <- l4tpm[,rownames(l4tpm_celltype)]

l4tpm_celltype$color <- c("black","red","green","blue")[l4tpm_celltype$SimpleDescription]


##### Order genes by type
#genetype_l4 <- dat_gpcr[,c("wbid","genetype")]
theorder <- order(l4tpm_datpcr$genetype)  
l4tpm_datpcr <- l4tpm_datpcr[theorder,]
l4tpm <- l4tpm[theorder,]


if(TRUE){
  l4tpm_genetype <- data.frame(row.names=unique(l4tpm_datpcr$genetype_red))
  l4tpm_genetype$color <- large_col_vector_ct[1:nrow(l4tpm_genetype)]
} else {
  l4tpm_genetype <- data.frame(row.names=unique(l4tpm_datpcr$genetype))
  l4tpm_genetype$color <- large_col_vector_ct[1:nrow(l4tpm_genetype)]
}



#order cells by expression level, then clustering
png("/home/mahogny/umea/project/changchun/gpcr/outfigs/rnaseq_heatmap.png", width = 1000, height = 1000)
heatmap.2(l4tpm, scale = "none", 
          col = hcl.colors(50),
          dendrogram="none",
          ColSideColors=l4tpm_celltype$color, 
          RowSideColors=l4tpm_genetype[l4tpm_datpcr$genetype_red,], #l4tpm_genetype$color, 
          #RowSideColors=l4tpm_genetype[l4tpm_datpcr$genetype,], #l4tpm_genetype$color, 
          Rowv=FALSE,
          Colv=FALSE,
          labRow=FALSE, labCol=FALSE,
          xlab="Cells",ylab="Genes",
          trace = "none", density.info = "none")
dev.off()

######## For legend
l4tpm_genetype
unique(l4tpm_celltype$SimpleDescription)
table(l4tpm_genetype[l4tpm_datpcr$genetype_red,])







##############################################################################
################ Expression analysis - most similar genes ####################
##############################################################################



######## Find closest gene for each
l4tpm <- qs::qread("L4.all.TPM.raw.qs")

#### Find agreeing set of genes
l4tpm_datpcr <- dat_gpcr_red#[dat_gpcr_red$highlight=="TRUE",]
l4tpm <- l4tpm[rownames(l4tpm) %in% l4tpm_datpcr$wbid,]
#mean(dat_gpcr_red$wbid %in% rownames(l4tpm)) #8% lost
rownames(l4tpm_datpcr) <- l4tpm_datpcr$wbid
l4tpm_datpcr <- l4tpm_datpcr[rownames(l4tpm),]


#normalize each gene; use pseudo log scale
for(i in 1:nrow(l4tpm)){
  l4tpm[i,] <- 1e6*l4tpm[i,]/sum(l4tpm[i,])
}
l4tpm <- log10(1+l4tpm)



rnaseq_cor <- cor(t(l4tpm))
rnaseq_mostsimilar <- NULL
for(i in 1:nrow(rnaseq_cor)){
  list_mostsimilar <- l4tpm_datpcr[order(rnaseq_cor[i,], decreasing = TRUE)[1:4],]$genename
  list_mostsimilar <- do.call(paste, c(as.list(list_mostsimilar), sep = ","))
  
  rnaseq_mostsimilar <- rbind(
    rnaseq_mostsimilar, 
    data.frame(from=l4tpm_datpcr$genename[i], to=list_mostsimilar))
}
rownames(rnaseq_mostsimilar) <- rnaseq_mostsimilar$from
rnaseq_mostsimilar$wbid <- l4tpm_datpcr$wbid
rnaseq_mostsimilar$genetype <- l4tpm_datpcr$genetype

#rnaseq_mostsimilar <- rnaseq_mostsimilar[order(rnaseq_mostsimilar$from),]
#rnaseq_mostsimilar["DUF-621",]

rnaseq_mostsimilar[rnaseq_mostsimilar$genetype=="DUF621",]

write.csv(rnaseq_mostsimilar, "newout//list_mostsimilar_rnaseq.csv")



