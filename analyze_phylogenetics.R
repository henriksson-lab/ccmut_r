################################################################################
## This script is for comparing GENES based on their sequence, and their      ##
## expression pattern                                                         ##
################################################################################


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
library(ape)
library(igraph)



##################### much moved to finaltreeanalysis/



##############################################################################
###################### output for clustalo ###################################
##############################################################################


if(FALSE){
  fileConn <- file("finalclustalanalysis/gpcr_list.fasta")
  writeLines(as.character(t(cbind(sprintf(">%s",dat_gpcr_red$wbid),dat_gpcr_red$aaseq))), fileConn)
  close(fileConn)
  # clustalo -i gpcr_list.fasta -o gpcr_list.aligned.fasta --threads=4
  # gpcr_list.aligned.fasta
}

##############################################################################
################## Compute distances from alignment ##########################
##############################################################################

myTree <- ape::read.tree(file="finalclustalanalysis/veryfasttree.dnd")


dat_gpcr <- read.csv("finalclustalanalysis/gpcr_final_list converted3.csv")
dat_gpcr$highlight[is.na(dat_gpcr$highlight)] <- FALSE

rownames(dat_gpcr) <- dat_gpcr$wbid
dat_gpcr_red <- dat_gpcr[myTree$tip.label,]

dat_gpcr_red$genetype_red <- dat_gpcr_red$genetype
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"sr")] <- "Chemo receptors"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"str")] <- "Chemo receptors"
dat_gpcr_red$genetype_red[dat_gpcr_red$genetype_red %in% c("Light receptors")] <- "Neuropeptide receptors"
dat_gpcr_red$genetype_red[str_starts(dat_gpcr_red$genetype_red,"Frizzled")] <- "Frizzled/Taste"
unique(dat_gpcr_red$genetype_red)


n <- 36
qual_col_pals_ct = brewer.pal.info[brewer.pal.info$category == 'qual',]
large_col_vector_ct = unlist(mapply(brewer.pal, qual_col_pals_ct$maxcolors, rownames(qual_col_pals_ct)))


possible_srh <- c(
  "F07C3.16",
  "F18E3.10",
  "F59B1.4",
  "F59B1.6",
  "T03E6.8",
  "T03E6.9",
  "W06G6.15",
  "Y70C5A.2",
  "Y70C5B.2",
  "F07B10.4",
  "F07B10.7",
  "F49C5.9",
  "Y62H9A.1")
dat_gpcr_red$genetype[dat_gpcr_red$genename %in% possible_srh] <- "NEW_SRH"

sort(unique(dat_gpcr_red$genetype))

############################
############################
############################
############################


#use_genetype <- dat_gpcr_red$genetype
use_genetype <- dat_gpcr_red$genetype_red



onegraph <- igraph::graph_from_data_frame(data.frame(myTree$edge), directed=FALSE)  
onegraph <- igraph::set_vertex_attr(onegraph, "genetype", 1:length(use_genetype), use_genetype)
alldist <- igraph::distances(onegraph)

# Need to get the right order
alldist <- alldist[as.character(1:nrow(alldist)),as.character(1:nrow(alldist))]

# Labels of vertices are whatever the labelled closest vertex is
all_vertex_labels <- NULL
for(i in 1:nrow(alldist)){
  closest_genetypes <- V(onegraph)$genetype[order(alldist[i,], decreasing = FALSE)]
  closest_genetypes <- closest_genetypes[!is.na(closest_genetypes)]
  all_vertex_labels <- c(all_vertex_labels, closest_genetypes[1])
}
length(all_vertex_labels)

# Comparison
data.frame(
  detected=all_vertex_labels[1:length(myTree$tip.label)],
  actual=use_genetype)

# For edges, pick the first edge, simply. This is not ideal but good enough
all_edge_labels <- all_vertex_labels[as.integer(myTree$edge[,1])]


############################
############################
############################
############################


#For each branch, set labels of closest neighbours. This will prune many...
new_vertex_labels <- rep("", max(myTree$edge))
for(onetype in unique(use_genetype)){
  allthenodes <- which(all_vertex_labels==onetype)
  thenode <- allthenodes[round(length(allthenodes)/2)]
  print(c(thenode,onetype))
  new_vertex_labels[thenode] <- onetype
}
myTree$node.label <- new_vertex_labels[-(1:myTree$Nnode)]


## Set up color table
if(FALSE){
  l4tpm_genetype <- data.frame(row.names=unique(use_genetype))
  l4tpm_genetype$color <- large_col_vector_ct[1:nrow(l4tpm_genetype)]
  #write.csv(l4tpm_genetype,"gpcr_list.treecolor.csv")
} else {
  l4tpm_genetype <- read.csv("finalclustalanalysis/gpcr_list.treecolor.csv", row.names = "X")
}
l4tpm_genetype


myTree$tip.color <- l4tpm_genetype[use_genetype,]
use_genetype

#tips are coded 1...n  in labels; not edges!
numLeaf <- nrow(dat_gpcr_red)

## Which are the final edges?
myTree$edge.isfinal <- rep(FALSE, nrow(myTree$edge))
myTree$edge.isfinal[myTree$edge[,2] <= numLeaf] <- TRUE
myTree$edge.isfinal.highlight <- myTree$edge[,2] %in% which(dat_gpcr_red$highlight)

## Edge lengths
myTree$edge.length <- 0.01 + 0.3*myTree$edge.isfinal.highlight
#myTree$edge.length <- 0.1*myTree$edge.length*(!myTree$edge.isfinal.highlight) + 0.1*myTree$edge.isfinal.highlight

## Edge colors
myTree$edge.color <- l4tpm_genetype[all_edge_labels,]

# Tip labels
myTree_label <- myTree
myTree$tip.label<- use_genetype
myTree$tip.label[duplicated(myTree$tip.label)] <- ""


pdf("newout/phylo1.pdf")
plot(myTree, 
           "u",
           edge.width = 0.3,
           use.edge.length = TRUE,
           tip.color = myTree$tip.color, 
           edge.color = myTree$edge.color, 
           show.tip.label = TRUE,
           cex = 0.45)
dev.off()
  




myTree_temp <- myTree_label
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="Frizzled/Taste"])] <- ""
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="Adhesion-type GPCRs"])] <- ""
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="Metabotropic neurotransmitter receptors" ])] <- ""
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="GPR180" ])] <- ""
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="C24A8.6" ])] <- ""
#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="TPRA1" ])] <- ""

#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype_red=="DUF1182" ])] <- ""
myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype=="NEW_SRH" ])] <- ""


#myTree_temp$tip.label[!(myTree_temp$tip.label %in% dat_gpcr_red$wbid[dat_gpcr_red$genetype=="str" ])] <- ""


dat_gpcr_red$wbid[dat_gpcr_red$genetype=="Frizzled/Taste"]
unique(dat_gpcr_red$genetype_red)#=="Frizzled/Taste"

pdf("newout/phylo1_label.pdf")
plot(myTree_temp, 
     "u",
     edge.width = 0.3,
     use.edge.length = TRUE,
     tip.color = myTree$tip.color, 
     edge.color = myTree$edge.color, 
     show.tip.label = TRUE,
     cex = 0.45)
dev.off()





myTree_temp <- myTree_label
myTree_temp$tip.label <- dat_gpcr_red$genename
myTree_temp$tip.label[!dat_gpcr_red$highlight] <- ""


pdf("newout/phylo1_highlight.pdf")
plot(myTree_temp, 
     "u",
     edge.width = 0.3,
     use.edge.length = TRUE,
     tip.color = myTree$tip.color, 
     edge.color = myTree$edge.color, 
     show.tip.label = TRUE,
     cex = 0.45)
dev.off()
