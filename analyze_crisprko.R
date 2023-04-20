################################################################################
## This script is for checking the alignments and figuring out if it worked out
################################################################################

library(stringr)
library(ggplot2)
library(plotly)
library(reshape2)


################################################################################
#### read mapping wbid - genesymbol
################################################################################

map_wbid_genesym_new <- read.csv("finalstrainannotation/genecoord_20230331.csv.gz",sep="\t")
colnames(map_wbid_genesym_new) <- c("wbid","transcript","genestart","geneend","chr","genesym")
map_wbid_genesym <- unique(map_wbid_genesym_new[,c("wbid","genesym")])


################################################################################
## Statistics on insertions/deletions
################################################################################

## Read stuff
tf <- read.table("finalstrainannotation/BOWTIEmut/CHS_1001_R1.fastq.gz.bowtie.bam.mut")  #CHS_1001_R1.fastq.gz.out.bam.mut
nr <- nrow(tf)
af <- list.files("finalstrainannotation/BOWTIEmut/",pattern = "*mut")
nummut <- length(af)
mat <- matrix(nrow = nr, ncol=nummut)
rownames(mat) <- tf[,2]
colnames(mat) <- str_replace_all(substr(af, 1,8),"_"," ")
#all(colnames(mat)==map_ngi_ccid$ccid) 
#colnames(mat) <- map_ngi_ccid$ccid
for(i in 1:length(af)){
  tf <- read.table(sprintf("finalstrainannotation/BOWTIEmut/%s",af[i]))
  mat[,i] <- tf[,3]
}

######## split the stats table
#reads del   ins   fine 
mut_cnt <- mat[tf$V1=="reads",]
mut_del <- mat[tf$V1=="del",]
mut_ins <- mat[tf$V1=="ins",]
mut_fine <- mat[tf$V1=="fine",]
mut_totc <- mat[tf$V1=="totc",]   #what is diff from cnt??

numgene <- nrow(mut_cnt)  #1695  

#del: number of CigarOperator.D
#ins: number of CigarOperator.I


mut_long <- cbind(
  melt(mut_cnt),
  melt(mut_fine)[,3],
  melt(mut_del)[,3],
  melt(mut_ins)[,3])
colnames(mut_long) <- c("wbid","ccid","cnt","fine","del","ins")

##Calculate rations
mut_long$pfine <- mut_long$fine/mut_long$cnt
mut_long$pdel  <- mut_long$del/mut_long$cnt
mut_long$pins  <- mut_long$ins/mut_long$cnt
mut_long$ratio <- log10(mut_long$pdel/mut_long$pfine)
mut_long$ratio[is.na(mut_long$ratio)] <- 0


mut_long$log_cnt <- log10(mut_long$cnt+1)



################################################################################
## Which genes were intended to be edited?
################################################################################

alledit <- unique(read.csv("finalstrainannotation/list_geneedit_strain_gene_NEW.csv")[,-1])
dim(alledit) #1747 entries
alledit$edited<-TRUE

dim(alledit) #1747


##################
mut_long <- merge(mut_long,alledit,all=TRUE)
mut_long$edited[is.na(mut_long$edited)] <- FALSE

#A silly label. remove?
mut_long$pcr <- "pcr-wt"
mut_long$pcr[mut_long$edited] <- "pcr-ko"  

#dim(mut_long)
#sum(mut_long$edited)  #1739  genes




################################################################################
## Add secondary MSA analysis
################################################################################

#dim(mut_long)
msa_diff <- read.csv("finalstrainannotation/msa_diff.csv", sep="\t")
mut_long <- merge(msa_diff,mut_long, all=TRUE)


### First attempt at annotation
mut_long$annotation <- "Inconclusive"
mut_long$annotation[!is.na(mut_long$numRE) & mut_long$numRE>0] <- "Edited"

table(mut_long$annotation)
table(mut_long$annotation[mut_long$edited])


################################################################################
## Which genes to include in the counts?
################################################################################

list_include_genes <- read.csv("finalstrainannotation/genes_for_fig1.csv")
#genomeseqeditfig1$wbid


################################################################################
## Add in manual annotation
################################################################################


manual_annot <- read.csv("finalstrainannotation/manual_gpcr_annot.csv")
dim(mut_long)
mut_long <- merge(mut_long, manual_annot, all=TRUE)

mut_long$annotation[!is.na(mut_long$manual_annotation)] <- mut_long$manual_annotation[!is.na(mut_long$manual_annotation)]

#table(mut_long$manual_annotation[!is.na(manual_annot$manual_annotation)])
#mut_long$manual_annotation[!is.na(manual_annot$manual_annotation)]
#mut_long$manual_annotat

table(mut_long$manual_annotation[!is.na(manual_annot$manual_annotation)])

#table(manual_annot$annotation)

table(mut_long$annotation)  #8 orig
table(mut_long$annotation[mut_long$edited])

################################################################################
## How much due to coverage?
################################################################################

#median(mut_long$cnt, na.rm = TRUE)
mut_long$annotation[mut_long$annotation=="Inconclusive" & (is.na(mut_long$cnt) | mut_long$cnt<10)] <- "Low coverage"


table(mut_long$annotation)
table(mut_long$annotation[mut_long$edited])


################################################################################
## Figure out annotation
################################################################################


################## volcano included in fig1
mut_long$color_by_onechs <- mut_long$pcr
mut_long$color_by_onechs[mut_long$ccid=="CHS 1050" & mut_long$edited] <- "y_CHS 1050"
mut_long$color_by_onechs[mut_long$color_by_onechs=="pcr-ko"] <- "x_ko"
mut_long <- mut_long[order(mut_long$color_by_onechs),]
ggplot(mut_long, aes(ratio, log_cnt, color=color_by_onechs)) + geom_point() + 
  ylab("Log(Counts)") + 
  xlab("Log10(p_del / p_fine)") +
  scale_color_manual(values=c("#999999", "#E69F00", "black"))
ggsave("newout/pdel_volcano_onehighlight.png")
ggsave("newout/pdel_volcano_onehighlight.pdf")




################################################################################
## Pie charts
################################################################################


sub_mut_long <- mut_long[mut_long$wbid %in% list_include_genes$wbid & mut_long$edited,]
dim(sub_mut_long)
have_num_genes <- length(unique(sub_mut_long$wbid))
have_num_genes

expected_num_genes_to_have <- length(unique(list_include_genes$wbid))  #possibly -1 because of lat-1?
expected_num_genes_to_have

length(alledit$wbid %in% list_include_genes)  #1747 targets

#########################
######################### pie chart for each target
#########################

stat_targeted <- as.data.frame(table(sub_mut_long$annotation))
colnames(stat_targeted) <- c("group","cnt")
stat_targeted$group_cnt <- sprintf("%s (%s)", stat_targeted$group, stat_targeted$cnt)
ggplot(stat_targeted, aes(x="", y=cnt, fill=group_cnt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() 
ggsave("newout/pie_pertarget.pdf", width = 4, height = 4)

sum(stat_targeted$cnt)   #1701



#########################
######################### pie chart for each gene
#########################

pie_geneedit <- data.frame(wbid=unique(sub_mut_long$wbid))
pie_geneedit$status <- "Inconclusive"
for(i in 1:nrow(pie_geneedit)){
  statforone <- unique(sub_mut_long$annotation[sub_mut_long$wbid==pie_geneedit$wbid[i]])
  if(length(statforone)==1){
    pie_geneedit$status[i] <- statforone
  } else {
    
    if("Edited" %in% statforone){
      pie_geneedit$status[i] <- "Edited"
    }
    print(statforone)
  }
}



stat_targeted <- as.data.frame(table(pie_geneedit$status))
colnames(stat_targeted) <- c("group","cnt")
stat_targeted$group_cnt <- sprintf("%s (%s)", stat_targeted$group, stat_targeted$cnt)
stat_targeted

sum(stat_targeted$cnt)  #1628

#stat_targeted <- stat_targeted[stat_targeted$group!="Unedited",]

ggplot(stat_targeted, aes(x="", y=cnt, fill=group_cnt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() 
ggsave("newout/pie_pergene.pdf", width = 4, height = 4)

sum(stat_targeted$cnt)









################################################################################
## Compare each gene to controls... t-test?
################################################################################

# 
# mut_long$color_by_onechs <- mut_long$pcr
# mut_long$color_by_onechs[mut_long$ccid=="CHS 1050"] <- "CHS 1050"
# alledit[alledit$ccid=="CHS 1050",]
# 
# #genomeseqeditfig1 <- read.csv("../Genome sequencing analysis for fig1.wbid.csv")
# #all(genomeseqeditfig1$wbid %in% alledit$wbid)
# genomeseqeditfig1$wbid[!(genomeseqeditfig1$wbid %in% alledit$wbid)]   #lat-1 is not in there. but it was never edited anyway (essential)
# sub_mut_long <- mut_long[mut_long$wbid %in% genomeseqeditfig1$wbid & mut_long$edited,]
###

minratio <- min(mut_long$ratio[!is.infinite(mut_long$ratio)])

sub_mut_long$pval <- NA
for(i in 1:nrow(sub_mut_long)){
  print(i)
  the_wbid <- sub_mut_long$wbid[i]

  neg_ratios <- mut_long$ratio[mut_long$wbid==the_wbid & !mut_long$edited]
  neg_ratios <- neg_ratios[!is.infinite(neg_ratios)]
  neg_ratios[is.infinite(neg_ratios)] <- minratio
  sub_mut_long$pval[i] <- 1-pnorm(sub_mut_long$ratio[i], mean=mean(neg_ratios), sd=sd(neg_ratios))
}


sub_mut_long$log_pval <- -log10(sub_mut_long$pval)

hist(-log10(sub_mut_long$pval))

#### p-values 
ggplot(sub_mut_long, aes(x=log_pval, y=log_cnt)) + geom_point() + xlab("-Log10 p-value") + ylab("Log10 Counts")
ggsave("newout/pval_of_targets.pdf")


p<-ggplot(sub_mut_long, aes(x=log_pval)) + 
  geom_histogram(color="black", fill="gray") 
p






################################################################################
## random crap
################################################################################

#these are the RE sites to look for
#taaaagcttataaataa      taaaagctttaa      taactcgagtaa taagaattcataaataa      taagaattctaa      taaggatcctaa      taatctagataa 


sub_mut_long <- mut_long[mut_long$wbid %in% list_include_genes$wbid & mut_long$edited,]
dim(sub_mut_long)
have_num_genes <- length(unique(sub_mut_long$wbid))
have_num_genes

expected_num_genes_to_have <- length(unique(list_include_genes$wbid))  #possibly -1 because of lat-1?
expected_num_genes_to_have



#sub_mut_long <- mut_long[mut_long$wbid %in% genomeseqeditfig1$wbid & mut_long$edited,]

#sum(sub_mut_long$pdel[sub_mut_long$edited] + sub_mut_long$pins[sub_mut_long$edited] >0.1, na.rm = TRUE)
#length(unique(sub_mut_long$wbid))  #1620 

#missing_genes <- should_have_genes[!(should_have_genes %in% sub_mut_long$wbid)]

#map_wbid_genesym[map_wbid_genesym$wbid %in% missing_genes,]
#Missing genes
#717   WBGene00000478   cfz-2  ... never mentioned in slurm out!   missing in BED and fa
# 4116  WBGene00003238   mig-1  not mentioned!
# 10496 WBGene00007432 srxa-18   #is in the bed-file... should be in CHS 1226. we have CHS 1226. but bam file is missing??
# 13802 WBGene00009707 srxa-16
# 18272 WBGene00012846 srxa-14
# 18273 WBGene00012847 srxa-15
# 32977 WBGene00044115 srxa-12
# 32978 WBGene00044116 srxa-13

#alledit[alledit$wbid %in% missing_genes,]
#ChS 1091, CHS 1226 missing. especially the latter. need to upload to AE??




