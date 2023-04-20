 
library(stringr)
library(ComplexHeatmap)
library(umap)
library(patchwork)
library(gridExtra)
library(ggcorrplot)
library(gplots)
library(reshape2)
library(reshape2)

################################################################################
####################### load data ##############################################
################################################################################

##### Load chemotaxis data
chemotaxis_dat <- read.csv("chemotaxis/chemotaxis_data.csv",as.is = TRUE)
colnames(chemotaxis_dat)
chemotaxis_meta <- read.csv("chemotaxis/chemotaxis_meta.csv")
chemotaxis_meta$condname <- sprintf("cond%s",1:nrow(chemotaxis_meta))
rownames(chemotaxis_meta) <- chemotaxis_meta$condname

##### Load phenotype data: avoidance, resistance, sensitive
list_phenotype <- read.csv("chemotaxis/list_phenotype.csv")

##### Load pathogen phenotypes
pathogen_avoidance <- read.csv("chemotaxis/pathogen_response_v3.csv")


map_strain_type <- read.csv("chemotaxis/mapping_strain_category.csv")  #For clustering
list_neuropep <- map_strain_type$ccid[map_strain_type$genetype=="neuropeptide muants"]
list_neuropep_recep <- map_strain_type$ccid[map_strain_type$genetype=="Neuropeptide receptor"]



################################################################################
####################### chemotaxis violins #####################################
################################################################################


allviolin <- NULL
for(i in 1:7){
  allviolin <- rbind(
    allviolin,
    data.frame(
      ccid=chemotaxis_dat$ccid, 
      condname=chemotaxis_meta$condition[i],
      compound=chemotaxis_meta$compound[i],
      val=chemotaxis_dat[,chemotaxis_meta$condname[i]])
  )  
}

p <- ggplot(allviolin, aes(x=condname, y=val, fill=compound, color=compound)) + 
  coord_flip()+
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.3), size=0.5, color="black") 
p
ggsave("newout/chemotaxis_violin_dil.pdf", height = 3)


###
allviolin <- NULL
for(i in 8:14){
  allviolin <- rbind(
    allviolin,
    data.frame(
      ccid=chemotaxis_dat$ccid, 
      condname=chemotaxis_meta$condition[i],
      compound=chemotaxis_meta$compound[i],
      val=chemotaxis_dat[,chemotaxis_meta$condname[i]])
  )  
}

p <- ggplot(allviolin, aes(x=condname, y=val, fill=compound, color=compound)) + 
  coord_flip()+
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.3), size=0.5, color="black") 
p
ggsave("newout/chemotaxis_violin_undil.pdf", height = 3)




################################################################################
####################### pathogen response - general ############################
################################################################################


### Simple correlation analysis
if(FALSE){
  corr_pathogen <- cor(pathogen_avoidance[,-1], method = "spearman")
  ggcorrplot(corr_pathogen)
}

## See how similar avoidance and improved measure is
if(FALSE){
  ggplot(pathogen_avoidance, aes(occupancy_change,avoidance_index_change))+ geom_point()
  thelm <- lm(occupancy_change ~ avoidance_index_change, data=pathogen_avoidance)
  pathogen_avoidance$combined_avoidance <- pathogen_avoidance$avoidance_index_change - pathogen_avoidance$occupancy_change/15
  
  
  ggplot(pathogen_avoidance, aes(avoidance_index_change,survival_change_2))+ geom_point() |
    ggplot(pathogen_avoidance, aes(combined_avoidance,survival_change_2))+ geom_point()
}



#cor.test(pathogen_avoidance$survival_change_2, pathogen_avoidance$avoidance_index_change)
# 0.2826092 %, pval < 2.74e-7

#Occupancy change reflects how many animals stay on the pathogen lawn compared to WT. 
#The higher value means more animals stay on the pathogen.

#Avoidance changes reflects how many animals leave the pathogen lawn compared to WT.
#The lower value means more animals leave the pathogen.

#We just need to use one of the values between occupancy and avoidance.
#But I am not sure which one highlights the difference better.

################################################################################
################# pathogen occupancy vs survival, scatter ######################
################################################################################

############# Test: correlation occupancy_change vs survival change
cor.test(pathogen_avoidance$occupancy_change, pathogen_avoidance$survival_change_2)
thelm <- lm(data=pathogen_avoidance, survival_change_2 ~ occupancy_change)
anova(thelm)  #p < 2.74e-7


#Supplementary Figure 4A, correlation plot, the sensitive and resistant strains indicated, but I saw the avoidance strains on the plot.


list_strain_3e_sensitive <- c(
  "CHS 1021",
  "CHS 1024",
  "CHS 1057",
  "CHS 1062",
  "CHS 1120",
  "CHS 1157",
  "CHS 1178",
  "CHS 1179",
  "CHS 1256",
  "CHS 10009",
  "CHS 10088",
  "CHS 10091",
  "CHS 1666"
)

list_strain_3g_resistant <- c(
  "CHS 1025",
  "CHS 1040",
  "CHS 1101",
  "CHS 1103",
  "CHS 1237",
  "CHS 1270",
  "CHS 10149",
  "CHS 10119"
)

pathogen_avoidance$lab_res_sens <- as.character(pathogen_avoidance$ccid)
pathogen_avoidance$col_res_sens<-"black"
pathogen_avoidance$col_res_sens[pathogen_avoidance$ccid %in% list_strain_3e_sensitive] <- "red"
pathogen_avoidance$col_res_sens[pathogen_avoidance$ccid %in% list_strain_3g_resistant] <- "blue"
pathogen_avoidance$lab_res_sens[pathogen_avoidance$col_res_sens=="black"] <- ""

#pathogen_avoidance$lab[abs(pathogen_avoidance$survival_change_2) < 0.3 & pathogen_avoidance$occupancy_change<5 & !(pathogen_avoidance$ccid %in% list_phenotype$ccid)] <- ""

ggplot(pathogen_avoidance, aes(occupancy_change,survival_change_2,label=lab_res_sens)) + 
  geom_point(size=1, color=pathogen_avoidance$col_res_sens) + 
  xlim(-1,20) +
  geom_smooth(method=lm, level=0.95, color="black") +
  xlab("Occupancy change") +
  ylab("Mean survival change") +
  geom_text( size=2, color=pathogen_avoidance$col_res_sens, nudge_y=0.01)
ggsave("newout/supfig 4a pathogen_occupancy.pdf")  


cor.test(pathogen_avoidance$survival_change_2, pathogen_avoidance$occupancy_change)
# 0.2826092 %, pval < 2.74e-7



############# Test: correlation occupancy_change vs survival change
cor.test(pathogen_avoidance$avoidance_index_change, pathogen_avoidance$survival_change_2)
#28% p < 4.52e-7
#thelm <- lm(data=pathogen_avoidance, survival_change_2 ~ avoidance_index_change)
#anova(thelm)  #p < 2.74e-7

ggplot(pathogen_avoidance, aes(avoidance_index_change,survival_change_2,label=lab_res_sens)) + 
  geom_point(size=1, color=pathogen_avoidance$col_res_sens) + 
  geom_smooth(method=lm, level=0.95, color="black") +
  xlab("Avoidance change") +
  ylab("Mean survival change") +
  geom_text( size=2, color=pathogen_avoidance$col_res_sens, nudge_y=0.01)
#ggsave("newout/supfig 4a pathogen_avoidance.pdf")  
ggsave("newout/supfig 4a pathogen_avoidance.pdf")  


################################################################################
################# pathogen avoidance, violin ###################################  Figure 3a
################################################################################

thefit <- density(pathogen_avoidance$avoidance_index_change)
pathogen_avoidance$jitter_advanced <- NA
for(i in 1:nrow(pathogen_avoidance)){
  pathogen_avoidance$jitter_advanced[i] <- thefit$y[which.min(abs(thefit$x - pathogen_avoidance$avoidance_index_change[i]))]/60
  pathogen_avoidance$jitter_advanced[i] <- rnorm(1, mean=0, sd=pathogen_avoidance$jitter_advanced[i])
}

pathogen_avoidance$label_here <- as.character(pathogen_avoidance$ccid)
list_strains_here <- c(
  "CHS 1020",
  "CHS 1021",
  "CHS 1024",
  "CHS 1057",
  "CHS 1103",
  "CHS 1127",
  "CHS 1157",
  "CHS 1223",
  "CHS 1256",
  "CHS 10009",
  "CHS 1666"
  )
pathogen_avoidance$label_here[!(pathogen_avoidance$ccid %in% list_strains_here)] <- ""
pathogen_avoidance$col_here <- "red"
pathogen_avoidance$col_here[pathogen_avoidance$label_here==""] <- "black"


pathogen_avoidance$const_0 <- 0 # "a"
pathogen_avoidance$const_0 <- factor(pathogen_avoidance$const_0)
ggplot(pathogen_avoidance, aes(x=const_0,y=avoidance_index_change)) + 
  geom_violin(fill="gray", color="gray")+
  scale_fill_manual(values=c("red", "black", "blue")) + 
  scale_color_manual(values=c("gray"))+
  geom_point(mapping=aes(x=1+jitter_advanced, avoidance_index_change), size=0.5, color=pathogen_avoidance$col_here)+
  ylab("Avoidance index change")+xlab("") + 
  #geom_text_repel(aes(x=const_0,y=avoidance_index_change,label=label_here), size=2)#, nudge_x=-.05)
  geom_text(aes(x=const_0,y=avoidance_index_change,label=label_here), size=2, nudge_x=-.15)+
 theme_classic()
ggsave("newout/fig3a.pdf", width = 3, height = 3)


################################################################################
################# pathogen avoidance, violin ###################################  Figure 3a - with occupancy
################################################################################

thefit <- density(pathogen_avoidance$occupancy_change)
pathogen_avoidance$jitter_advanced <- NA
for(i in 1:nrow(pathogen_avoidance)){
  pathogen_avoidance$jitter_advanced[i] <- thefit$y[which.min(abs(thefit$x - pathogen_avoidance$occupancy_change[i]))]/60
  pathogen_avoidance$jitter_advanced[i] <- rnorm(1, mean=0, sd=pathogen_avoidance$jitter_advanced[i])
}

pathogen_avoidance$label_here <- as.character(pathogen_avoidance$ccid)
list_strains_here <- c(
  "CHS 1020",
  "CHS 1021",
  "CHS 1024",
  "CHS 1057",
  "CHS 1103",
  "CHS 1127",
  "CHS 1157",
  "CHS 1223",
  "CHS 1256",
  "CHS 10009",
  "CHS 1666"
)
pathogen_avoidance$label_here[!(pathogen_avoidance$ccid %in% list_strains_here)] <- ""
pathogen_avoidance$col_here <- "red"
pathogen_avoidance$col_here[pathogen_avoidance$label_here==""] <- "black"


pathogen_avoidance$const_0 <- 0 # "a"
pathogen_avoidance$const_0 <- factor(pathogen_avoidance$const_0)
ggplot(pathogen_avoidance, aes(x=const_0,y=occupancy_change)) + 
  geom_violin(fill="gray", color="gray")+
  scale_fill_manual(values=c("red", "black", "blue")) + 
  scale_color_manual(values=c("gray"))+
  geom_point(mapping=aes(x=1+jitter_advanced, occupancy_change), size=0.5, color=pathogen_avoidance$col_here)+
  ylab("Avoidance index change")+xlab("") + 
  #geom_text_repel(aes(x=const_0,y=avoidance_index_change,label=label_here), size=2)#, nudge_x=-.05)
  geom_text(aes(x=const_0,y=occupancy_change,label=label_here), size=2, nudge_x=-.15)+
  theme_classic()
ggsave("out/fig3a occupancy.pdf", width = 3, height = 3)



################################################################################
################# pathogen survival, violin ####################################  Figure 3d
################################################################################

thefit <- density(pathogen_avoidance$survival_change_2)
pathogen_avoidance$jitter_advanced <- NA
for(i in 1:nrow(pathogen_avoidance)){
  pathogen_avoidance$jitter_advanced[i] <- thefit$y[which.min(abs(thefit$x - pathogen_avoidance$survival_change_2[i]))]/30
  pathogen_avoidance$jitter_advanced[i] <- rnorm(1, mean=0, sd=pathogen_avoidance$jitter_advanced[i])
}

pathogen_avoidance$label_here <- as.character(pathogen_avoidance$ccid)
list_strains_here <- c(
  "CHS 1057",
  "CHS 1157",
  "CHS 1062",
  "CHS 1021",
  "CHS 10091",
  "CHS 1178",
  "CHS 1024",
  "CHS 1120",
  "CHS 10009",
  "CHS 1666",
  "CHS 1256",
  "CHS 1179",
  "CHS 10088",
  "CHS 1040",
  "CHS 1101",
  "CHS 1237",
  "CHS 1025",
  "CHS 1270",
  "CHS 10149",
  "CHS 10119",
  "CHS 1103"
)
pathogen_avoidance$label_here[!(pathogen_avoidance$ccid %in% list_strains_here)] <- ""
pathogen_avoidance$col_here <- "red"
pathogen_avoidance$col_here[pathogen_avoidance$label_here==""] <- "black"

pathogen_avoidance$const_0 <- 0 # "a"
pathogen_avoidance$const_0 <- factor(pathogen_avoidance$const_0)
ggplot(pathogen_avoidance, aes(x=const_0,y=survival_change_2)) + 
  geom_violin(fill="gray", color="gray")+
  scale_fill_manual(values=c("red", "black", "blue")) + 
  scale_color_manual(values=c("gray"))+
  geom_point(mapping=aes(x=1+jitter_advanced, survival_change_2), size=0.5, color=pathogen_avoidance$col_here)+
  ylab("Mean survival change")+xlab("") + 
  geom_text(aes(x=const_0,y=survival_change_2,label=label_here), size=2, nudge_x=-.15)+
  theme_classic()
ggsave("newout/fig3d.pdf", width = 3, height = 3)







################################################################################
########### correlation: pathogen avoidance vs chemotaxis ######################
################################################################################

if(FALSE){
  dat_path_chemo <- merge(chemotaxis_dat, pathogen_avoidance)
  colnames(dat_path_chemo) <- c("ccid","genotype",chemotaxis_meta$condition,"occupancy_change","avoidance_index_change","survival_change_2")
  dat_path_chemo <- dat_path_chemo[colnames(dat_path_chemo)!="undiluted 2-nonanone"]
  
  corr_path_chemo <- cor(dat_path_chemo[,-(1:2)], method = "spearman")
  for(i in 1:nrow(corr_path_chemo)){
    corr_path_chemo[i,i] <- 0
  }
  ggcorrplot(corr_path_chemo)
}




################################################################################
#################### Violin: chemotaxis,  ######################################
################################################################################

chemotaxis_long <- chemotaxis_dat[,-(1:2)]
rownames(chemotaxis_long) <- chemotaxis_dat$ccid
chemotaxis_long <- melt(as.matrix(chemotaxis_long))
colnames(chemotaxis_long) <- c("ccid","condname","value")
chemotaxis_meta$condname
chemotaxis_long <- merge(chemotaxis_meta, chemotaxis_long)

#for each diluted condition, label top 5
tocolor <- NULL
for(curcond in chemotaxis_meta$condname[chemotaxis_meta$diluted]){
  chemotaxis_long <- chemotaxis_long[order(chemotaxis_long$value),]
  tocolor <- rbind(
    tocolor,
    head(chemotaxis_long[chemotaxis_long$diluted & chemotaxis_long$condname==curcond,],n=5))
}
for(curcond in chemotaxis_meta$condname[!chemotaxis_meta$diluted]){
  chemotaxis_long <- chemotaxis_long[order(chemotaxis_long$value, decreasing = TRUE),]
  tocolor <- rbind(
    tocolor,
    head(chemotaxis_long[!chemotaxis_long$diluted & chemotaxis_long$condname==curcond,],n=5))
}
tocolor <- tocolor[,c("ccid","condname")]
tocolor$label <- as.character(tocolor$ccid)
chemotaxis_long <- merge(chemotaxis_long, tocolor, all.x=TRUE)
chemotaxis_long$label[is.na(chemotaxis_long$label)]<-""


chemotaxis_long$col_here <- "black"
chemotaxis_long$col_here[chemotaxis_long$label!=""] <- "red"
chemotaxis_long$shortname <- factor(chemotaxis_long$shortname, levels=c("DA","PZ","TMT","IAA","BZ","BU","PD",       "NON","OCT"))

ggplot(chemotaxis_long[chemotaxis_long$diluted,], aes(x=shortname,y=value, fill=shortname)) + 
  geom_violin(scale = "width") + 
  geom_jitter(width=0.2, color=chemotaxis_long$col_here[chemotaxis_long$diluted]) + 
  geom_text(aes(x=shortname,y=value, fill=shortname, label=label), nudge_x = 0.4) +
  xlab("Diluted")+ylab("")
ggsave("out/fig5b.pdf")

ggplot(chemotaxis_long[!chemotaxis_long$diluted,], aes(x=shortname,y=value, fill=shortname)) + 
  geom_violin(scale = "width") + 
  geom_jitter(width=0.2, color=chemotaxis_long$col_here[!chemotaxis_long$diluted]) + 
  geom_text(aes(x=shortname,y=value, fill=shortname, label=label), nudge_x = 0.4) +
  xlab("Undiluted")+ylab("")
ggsave("newout/fig5d.pdf")



################################################################################
#################### UMAP: chemotaxis data #####################################
################################################################################


## Merge and rescale
all_pheno <- chemotaxis_dat[,-2]
all_pheno[,-1] <- scale(all_pheno[,-1])
all_pheno <- all_pheno[,colnames(all_pheno)!="cond13"]

#UMAP or heatmap to summarize the pathogen and chemotaxis data as in the paper 
#(https://www.science.org/doi/10.1126/sciadv.ade1249).


all_pheno <- merge(all_pheno, map_strain_type)

###################### chemotax
### umap based on sequence
custom.config <- umap.defaults
custom.config$random_state <- 124
#custom.config$n_neighbors <- 20
#custom.config$min_dist  <- 0.3  #from 0.1 to 1
umap.dist <- umap(all_pheno[,str_starts(colnames(all_pheno),"cond")], config=custom.config)
all_pheno$chemotax_umap_x <- umap.dist$layout[,1]
all_pheno$chemotax_umap_y <- umap.dist$layout[,2]

### Vanilla
all_pheno$genetype_neuro <- all_pheno$genetype   #<- all_pheno$ccid %in% list_neuropep
all_pheno$genetype_neuro[!(all_pheno$ccid %in% c(list_neuropep, list_neuropep_recep))] <- "other type"
# ggplot(all_pheno, aes(umap_x, umap_y, label=ccid, color=genetype_neuro))  + geom_point() + 
#   scale_color_manual(values=c("#FF0000", "#00FF00", "#000000"))
# ggsave("umap_chemotaxANDpathogen.pdf")
ggplot(all_pheno, aes(chemotax_umap_x, chemotax_umap_y, label=ccid, color=genetype_neuro))  + geom_point()+ 
  scale_color_manual(values=c("#FF0000", "#00FF00", "#000000"))
ggsave("newout/umap_chemotaxONLY.pdf")

#How are neuropeptide muants different from the rest?

alltests <- NULL
for(i in 2:17){
  onet <- t.test(
    all_pheno[all_pheno$ccid %in% list_neuropep,i],
    all_pheno[!(all_pheno$ccid %in% list_neuropep),i]
  )
  alltests <- rbind(alltests, data.frame(condname=colnames(all_pheno)[i], pval=onet$p.value))
}
alltests

alltests <- merge(alltests,chemotaxis_meta,all=TRUE)
alltests <- alltests[order(alltests$pval),]
alltests

ggplot(all_pheno, aes(umap_x, umap_y, label=ccid, color=genetype_neuro))  + geom_point() + 
  scale_color_manual(values=c("#FF0000", "#00FF00", "#000000"))




################################################################################
#################### Heatmap: chemotaxis data ##################################
################################################################################

# many options https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/

## Merge and rescale
all_pheno <- chemotaxis_dat[,-2]  #merge(chemotaxis_dat[,-2], pathogen_avoidance[,-2])
all_pheno[,-1] <- scale(all_pheno[,-1], scale = FALSE, center = FALSE)   #likely not needed
all_pheno <- merge(all_pheno,map_strain_type)
rownames(all_pheno) <- all_pheno$ccid
all_pheno <- all_pheno[,-1]

colnames(all_pheno)[colnames(all_pheno) %in% rownames(chemotaxis_meta)] <- chemotaxis_meta[colnames(all_pheno)[colnames(all_pheno) %in% rownames(chemotaxis_meta)],]$condition

all_pheno$genetype <- factor(all_pheno$genetype, levels = c("Chemoreceptor","Adhesion-type GPCR","Frizzled_Taste2","GPR180","Light receptors","Metabotropic neurotransmitter receptors","neuropeptide muants","Neuropeptide receptor"))
all_pheno <- all_pheno[order(all_pheno$genetype),]
rownames(all_pheno)<- sprintf("%s  %s",rownames(all_pheno), all_pheno$genetype)
all_pheno_genetype <- all_pheno$genetype
all_pheno <- all_pheno[,!(colnames(all_pheno) %in% "genetype")]






png("newout/heatmap_fig4c.png", width = 1600, height = 800)
#pdf("newout/heatmap_fig4c.pdf", width = 100, height = 15)

ha <- HeatmapAnnotation(bar = all_pheno_genetype,
                        show_legend = FALSE,
                        col = list(bar = c(
                          "Adhesion-type GPCR"="black",                      
                          "Chemoreceptor"="blue",                           
                          "Frizzled_Taste2"="black",                        
                          "GPR180"="black",                                  
                          "Light receptors"="black",                         
                          "Metabotropic neurotransmitter receptors"="black",
                          "neuropeptide muants"="red",                     
                          "Neuropeptide receptor"="green"              
                        )))
ht <- Heatmap(t(as.matrix(all_pheno)), 
              name = "", 
              column_title = "", row_title = "",
              cluster_columns = FALSE,
              top_annotation = ha,
              #show_heatmap_legend =TRUE,
              show_heatmap_legend =FALSE,
              row_names_gp = gpar(fontsize = 30) 
)
draw(ht, padding = unit(c(2, 20, 20, 200), "mm"))
dev.off()





################################################################################
####################### chemotaxis similarity matrix ###########################
################################################################################

c_chemotaxis_dat <- chemotaxis_dat
rownames(c_chemotaxis_dat) <- c_chemotaxis_dat$ccid
c_chemotaxis_dat <- c_chemotaxis_dat[,-(1:2)]
c_chemotaxis_dat <- scale(c_chemotaxis_dat)

similarity_receptor_ligand <- merge(
  data.frame(gene_receptor=map_strain_type$ccid[map_strain_type$genetype=="Neuropeptide receptor"]),
  data.frame(gene_ligand=map_strain_type$ccid[map_strain_type$genetype=="neuropeptide muants"])
)

similarity_receptor_ligand$cor <- NA
similarity_receptor_ligand$pval <- NA
for(i in 1:nrow(similarity_receptor_ligand)){
  v1 <- c_chemotaxis_dat[similarity_receptor_ligand$gene_receptor[i],]
  v2 <- c_chemotaxis_dat[similarity_receptor_ligand$gene_ligand[i],]
  testval <- cor.test(as.double(v1),as.double(v2))
  similarity_receptor_ligand$cor[i] <- testval$estimate  
  similarity_receptor_ligand$pval[i] <- testval$p.value  
}
similarity_receptor_ligand <- similarity_receptor_ligand[order(similarity_receptor_ligand$pval),]
similarity_receptor_ligand$quartile <- (1:nrow(similarity_receptor_ligand))/nrow(similarity_receptor_ligand)
similarity_receptor_ligand

write.csv(similarity_receptor_ligand, "newout/corr_ligand_receptor.csv")

if(FALSE){
  similarity_receptor_ligand[
    "CHS 10063"==similarity_receptor_ligand$gene_ligand & 
      "CHS 1025"==similarity_receptor_ligand$gene_receptor,]
  
  similarity_receptor_ligand[
    "CHS 10091"==similarity_receptor_ligand$gene_ligand & 
      "CHS 1178"==similarity_receptor_ligand$gene_receptor,]
}

unmelted <- acast(similarity_receptor_ligand, gene_receptor~gene_ligand, value.var = "cor")
pdf("newout/corr_ligand_receptor.pdf")
ht <- Heatmap(
  unmelted,#as.matrix(unm), 
  name = "Correlation",
  column_title = "Receptor", 
  row_title = "Ligand",
  column_names_gp = gpar(fontsize = 5), 
  row_names_gp = gpar(fontsize = 5) 
)
draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()







################################################################################
####################### chemotaxis AND survival/occ similarity matrix ##########
################################################################################


rownames(pathogen_avoidance) <- pathogen_avoidance$ccid
c_chemotaxis_dat <- cbind(
  chemotaxis_dat,
  pathogen_avoidance[chemotaxis_dat$ccid,c("occupancy_change","survival_change_2")]
)
rownames(c_chemotaxis_dat) <- c_chemotaxis_dat$ccid
c_chemotaxis_dat <- c_chemotaxis_dat[,-(1:2)]
c_chemotaxis_dat <- c_chemotaxis_dat[!is.na(c_chemotaxis_dat$occupancy_change),]
c_chemotaxis_dat <- scale(c_chemotaxis_dat)


similarity_receptor_ligand <- merge(
  data.frame(gene_receptor=map_strain_type$ccid[map_strain_type$genetype=="Neuropeptide receptor"]),
  data.frame(gene_ligand=map_strain_type$ccid[map_strain_type$genetype=="neuropeptide muants"])
)

similarity_receptor_ligand$cor <- NA
similarity_receptor_ligand$pval <- NA
for(i in 1:nrow(similarity_receptor_ligand)){
  v1 <- c_chemotaxis_dat[similarity_receptor_ligand$gene_receptor[i],]
  v2 <- c_chemotaxis_dat[similarity_receptor_ligand$gene_ligand[i],]
  testval <- cor.test(as.double(v1),as.double(v2))
  similarity_receptor_ligand$cor[i] <- testval$estimate  ###good!
  similarity_receptor_ligand$pval[i] <- testval$p.value  ###good!
}
similarity_receptor_ligand <- similarity_receptor_ligand[order(similarity_receptor_ligand$pval),]
similarity_receptor_ligand$quartile <- (1:nrow(similarity_receptor_ligand))/nrow(similarity_receptor_ligand)
similarity_receptor_ligand

write.csv(similarity_receptor_ligand, "newout/corr_ligand_receptor_withsurvivalocc.csv")

if(FALSE){
  similarity_receptor_ligand[
    "CHS 10063"==similarity_receptor_ligand$gene_ligand & 
      "CHS 1025"==similarity_receptor_ligand$gene_receptor,]
  
  similarity_receptor_ligand[
    "CHS 10091"==similarity_receptor_ligand$gene_ligand & 
      "CHS 1178"==similarity_receptor_ligand$gene_receptor,]
}

unmelted <- acast(similarity_receptor_ligand, gene_receptor~gene_ligand, value.var = "cor")
pdf("newout/corr_ligand_receptor_withsurvivalocc.pdf")
ht <- Heatmap(
  unmelted,#as.matrix(unm), 
  name = "Correlation",
  column_title = "Receptor", 
  row_title = "Ligand",
  column_names_gp = gpar(fontsize = 5), 
  row_names_gp = gpar(fontsize = 5) 
)
draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()
