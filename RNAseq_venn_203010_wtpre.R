library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(made4)
library(VennDiagram)

# this script calculates log2fc with comparing each genotype to its pre condition and plot heatmaps

# PRE is WT 

# first part is like RNAseq_rep_PCA_190115, putting all the counts in one long dataset

setwd("C:/Siavash/Codes/RNAseq/data")
counts_all <- read.table("count.tsv")
counts_all <- as.data.frame(counts_all)

marx_list2 <- c("gene","3-7-E","2-3-W","4-1-E","10-1-E","3-1-E","4-3-W", "8-1-E","1-1-W","3-3-E","7-3-W","4-2-W","7-3-E","12-7-E","11-1-E",
                "11-3-W","3-7-W","11-2-E","7-2-W","10-9-E","12-9-E","10-3-E","7-1-E","1-5-W","3-9-E","11-1-W","4-1-W","11-3-E","10-7-W",
                "8-3-W","3-3-W","13-1-E","10-1-W","1-1-E","2-9-E","3-9-W","2-3-E","13-3-W","3-1-W","10-3-W","1-5-E","7-1-W","8-2-E","8-3-E",
                "4-2-E","2-9-W","2-7-E","13-1-W","10-7-E","12-9-W","13-3-E","8-1-W","7-2-E","10-9-W","12-7-W","4-3-E","11-2-W","2-1-W",
                "2-1-E", "8-2-W", "2-7-W")

colnames(counts_all) <- marx_list2
counts_all <- counts_all[-1,] # remove name of treaments 

###################new 

key_id_full <- read.csv("key_id_full2.csv")
key_id_full$Sample <- gsub('.','-',key_id_full$Sample,fixed = T)
counts_all_m <- melt(counts_all,id.vars  ="gene",variable.name = "Sample",value.name ="counts")
#counts_all_m <- gather(counts_all, gene, counts) # alternative 
counts_all_m2 <- merge(counts_all_m,key_id_full,by="Sample")

# now smaller dataset 
mext_counts <- read_csv(file="C:/Siavash/Codes/RNAseq/data/m_extorquens_counts_JB.csv") %>%
  dplyr::rename(feat_id = gene) %>%
  gather(sample, counts,2:22) 

#makes matrix for DESeq2 and changes names to Marx lab nomenclature
sample <- data.frame(sample = c("13E","13W","142W","47E","47W","48E","48W","49E","49W","77E","77W","78E","78W","79E","79W","87E","87W","88E","88W",  "89E","89W"))
marx_list <- c("e-pre","wt-pre","wt-f-720","e-40","wt-40","e-k-40","wt-k-40","e-f-40","wt-f-40","e-180","wt-180","e-k-180","wt-k-180","e-f-180","wt-f-180","e-360","wt-360","e-k-360","wt-k-360","e-f-360","wt-f-360")
name_change <- data.frame(sample)
name_change$marx <- marx_list
#write_csv(name_change, path="key_id.csv")
key_id <- read.csv("key_id2.csv")
colnames(mext_counts)[1:2] <- c("gene","Sample")
small_data <- merge(mext_counts,key_id,by="Sample")
all_data <- rbind(small_data,counts_all_m2) # this dataset has everything in long format 
all_data <- all_data[,-1]
colnames(all_data)[3:6] <- c('genotype','treatment','time','replicate')

all_data$genotype <- as.character(all_data$genotype)
all_data$treatment <- as.character(all_data$treatment)
all_data$time <- as.character(all_data$time)

# a function to cast wt_pre and e_pre based on each replicate's counts 
casting <- function(data){
  return(dcast(data,gene+genotype+treatment+time~replicate,value.var = "counts"))
}

# wt_pre_short and e_pre_short are reference treatments 
wt_pre <- subset(all_data,genotype=='WT' & time=='pre' & treatment=='pre')
wt_pre_short <- as.matrix(casting(wt_pre))

e_pre <- subset(all_data,genotype=='efgA' & time=='pre' & treatment=='pre')
e_pre_short <- casting(e_pre)

# this is a martix (short data frame) that keeps all normalized counts  
significance <- data.frame(gene=unique(all_data$gene))

# all_data2 keeps all the log2FC in long format for plotting change of log2fc in sepcific genes over the time 
all_data2 <- subset(all_data,replicate==1)
all_data2 <- all_data2[,-6]
all_data2$counts <- 0
all_data2$padj <- 0
all_data2 <- all_data2[,c(1,3,4,5,2,6)]

for (g in unique(all_data$genotype)){
  for (tr in unique(all_data$treatment)[-1]){
    for (ti in c('5','20','40','180','360')){
      if(dim(subset(all_data,genotype==g & treatment==tr & time==ti))[1]>0){
        sample <- subset(all_data,genotype==g & treatment==tr & time==ti)
        sample_short <- casting(sample)
        
        PRE <- wt_pre_short[,5:7] # THIS LINE BEEN CHANGED, WE HAVE ONLY ONE PRE  
        colnames(PRE) <- paste(unique(wt_pre$genotype),unique(wt_pre$time),sort(as.numeric(unique(wt_pre$replicate))),sep = '.')
        
        df <- sample_short[,5:dim(sample_short)[2]]
        colnames(df) <- paste(unique(sample$genotype),unique(sample$treatment),unique(sample$time),sort(as.numeric(unique(sample$replicate))),sep = '.')
        df <- cbind(PRE,df)
        df <-apply(df, 2,as.numeric)
        condition <- c(rep(1,3),rep(2,length(unique(sample$replicate))))
        colData <- data.frame(condition)
        rownames(colData) <- colnames(df)
        colData$condition <- as.factor(colData$condition)
        
        #normalized counts, adds gene descriptions, saves
        mext_dds <- DESeqDataSetFromMatrix(countData = df, colData = colData, design = ~ condition)
        dds <- DESeq(mext_dds)
        res <- results(dds)
        res <- as.data.frame(res)
        all_data2$counts[all_data2$genotype==g & all_data2$treatment==tr & all_data2$time==ti] <- res$log2FoldChange
        all_data2$padj[all_data2$genotype==g & all_data2$treatment==tr & all_data2$time==ti] <- res$padj
        #res$log2FoldChange[is.na(res$padj)] <- 0 # if we want to plot only significant changes
        #res$log2FoldChange[res$padj>=0.001] <- 0 # if we want to plot only significant changes 
        #res$log2FoldChange[res$log2FoldChange>0 & res$log2FoldChange<1] <- 0 # cut of value based on log2FC, keep>2X 
        #res$log2FoldChange[res$log2FoldChange<0 & res$log2FoldChange>-1] <- 0 # cut of value based on log2FC, keep <-2X
        logf <- data.frame(res$padj)
        colnames(logf) <- paste(unique(sample$genotype),unique(sample$treatment),unique(sample$time),sep = '.')
        significance <- cbind(significance,logf)
      }
      else {}
    }
  }
}


# this function extract name of genes that are significantly up/down regualted in each treatment from the signficance dataset
cud <- function(ge,tr,ti,d){
  id <- paste(ge,'.',tr,'.',ti,sep = '')
  if (d =='u') return(as.character(significance$gene[significance[id]>0])) else return(as.character(significance$gene[significance[id]<0]))
}

# this loop makes list of genes names with this format: [genotype][treatment][time]_[up or down]
for (ti in c('5','20','40','180','360')){
  for (ge in c('WT','efgA')){
    for (di in c('u','d')){
      for (tr in c('none','form','kan')){
        if (tr=='kan' & (ti=='5' || ti=='20')){} # we don't have time 5 and 20 with kanamycin 
        else{
          id <- ''
          if (tr=='form') id <- 'f' else if (tr=='kan') id <- 'k'
          if (ge=='WT'){
            if (di=='u') {
              var <- paste('w',id,ti,'_up',sep='')
              assign(var,cud(ge,tr,ti,di))
            }
            else {
              var <- paste('w',id,ti,'_down',sep='')
              assign(var,cud(ge,tr,ti,di))
            }
          } else {
            if (di=='u') {
              var <- paste('e',id,ti,'_up',sep='')
              assign(var,cud(ge,tr,ti,di))
            }
            else {
              var <- paste('e',id,ti,'_down',sep='')
              assign(var,cud(ge,tr,ti,di))
            }
          }
        }
      }
    }
  }
}

# genes that are commonly up-regulated in WT from 5 min to 180 min (before stationary) in no-treatment 
wpre_up <- unique(c(w5_up,w20_up,w40_up,w180_up))

# genes that are commonly down-regulated in WT from 5 min to 180 min (before stationary) in no-treatment 
wpre_down <- unique(c(w5_down,w20_down,w40_down,w180_down))

# genes that are commonly up-regulated in defgA from 5 min to 180 min (before stationary) in no-treatment 
epre_up <- unique(c(e5_up,e20_up,e40_up,e180_up))

# genes that are commonly down-regulated in defgA from 5 min to 180 min (before stationary) in no-treatment 
epre_down <- unique(c(e5_down,e20_down,e40_down,e180_down))


# Venn diagrams 
# fig 3.5 left in dissertation 
overlap_pre <- list(wpre_up,wpre_down,epre_up,epre_down)
names(overlap_pre) <- c("w.pre.u","w.pre.d","e.pre.u","e.pre.d")
venn.diagram(overlap_pre,height=1500,width=2500,filename = "pre_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.5 right in dissertation 
overlap_pre360 <- list(w360_up,w360_down,e360_up,e360_down)
names(overlap_pre360) <- c("w.pre.360.u","w.pre.360.d","e.pre360.u.","e.pre.360.d")
venn.diagram(overlap_pre360,height=1500,width=2500,filename = "pre_360_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.7 left in dissertation  
overlap_k_wt <- list(wk360_up,wk360_down,wk180_up,wk180_down)
names(overlap_k_wt) <- c("w.k.360.u","w.k.360.d","w.k.180.u","w.k.180.d")
venn.diagram(overlap_k_wt,height=1500,width=2500,filename = "w_k_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.7 right in dissertation 
overlap_k_e <- list(ek360_up,ek360_down,ek180_up,ek180_down)
names(overlap_k_e) <- c("e.k.360.u","e.k.360.d","e.k.180.u","e.k.180.d")
venn.diagram(overlap_k_e,height=1500,width=2500,filename = "e_k_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.8 in dissertation 
overlap_k_360 <- list(wk360_up,wk360_down,ek360_up,ek360_down)
names(overlap_k_360) <- c("w.k.360.u","w.k.360.d","e.k.360.u","e.k.360.d")
venn.diagram(overlap_k_360,height=1500,width=2500,filename = "k_venn360.png",fill=c("red","blue","chocolate3","darkblue"))

# fig sup1 in paper 
overlap_k_180 <- list(wk180_up,wk180_down,ek180_up,ek180_down)
names(overlap_k_180) <- c("w.k.180.u","w.k.180.d","e.k.180.u","e.k.180.d")
venn.diagram(overlap_k_180,height=1500,width=2500,filename = "k_venn180.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.10 left in dissertation 
overlap_wf520 <- list(wf5_up,wf5_down,wf20_up,wf20_down)
names(overlap_wf520) <- c("w.f.5.u","w.f.5.d","w.f.20.u","w.f.20.d")
venn.diagram(overlap_wf520,height=1500,width=2500,filename = "wf_520_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.10 right in dissertation 
overlap_ef520 <- list(ef5_up,ef5_down,ef20_up,ef20_down)
names(overlap_ef520) <- c("e.f.5.u","e.f.5.d","e.f.20.u","e.f.20.d")
venn.diagram(overlap_ef520,height=1500,width=2500,filename = "ef_520_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.11 in dissertation 
overlap_f <- list(wf5_up,wf5_down,ef5_up,ef5_down)
names(overlap_f) <- c("w.f.5.u","w.f.5.d","e.f.5.u","e.f.5.d")
venn.diagram(overlap_f,height=1500,width=2500,filename = "f_venn2.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 3.12 in dissertation 
overlap_f_k360 <- list(wk360_up,wk360_down,wf5_up,wf5_down)
names(overlap_f_k360) <- c("w.k.360.u","w.k.360.d","w.f.5.u","w.f.5.d")
venn.diagram(overlap_f_k360,height=1500,width=2500,filename = "f_k360_venn.png",fill=c("red","blue","chocolate3","darkblue"))

# fig 9 for the paper 
overlap_f_k180 <- list(wk180_up,wk180_down,wf5_up,wf5_down)
names(overlap_f_k180) <- c("w.k.180.u","w.k.180.d","w.f.5.u","w.f.5.d")
venn.diagram(overlap_f_k180,height=1500,width=2500,filename = "f_k180_venn.png",fill=c("red","blue","chocolate3","darkblue"))


# looking at the genes 
genes <- read.csv('mext_desc.csv')
# gene_list is the set we want to look at, overlap_pre is an example 
gene_list <- overlap_pre 
gene_list2 <- subset(genes,feat_id %in% unlist(gene_list))
# looking at the intersections
subsets <- venn(gene_list, show.plot=FALSE)

## plotting log2fc over time

# in all_data2 I want 0 as time for pre (in fact it's -45)
all_data2$time[all_data2$time=='pre'] <- 0

## list of candidate genes 
# list of beneficial alleles 
efgA_order <- paste("Mext_",c(4158,4271,"0606",3949,1636,"0925"),sep = "")
names(efgA_order) <- c("efgA","efgA_hom","efgB","def_pdf","def","marR") 

# list of formaldehyde oxidation genes 
formaldehyde_order <- paste("Mext_",c(seq(4137,4151),1834,1829,1831,seq(1824,1827),seq(4581,4582),seq(4404,4406),
                                      c("0389","0390","0391"),seq(2104,2105),"0414",1798,1797),sep = '')


names(formaldehyde_order) <- c(paste("mxa",c('B','H','E','D','L','K','C','A','S','R','I','G','J','F','W'),sep=''),
                               "fae","mtdB","mch",paste("fhc",c('C','D','A','B'),sep = ''),"fdh1B","fdh1A","fdh2C",
                               "fdh2B","fdh2A","fdh3A","fdh3B","fdh3C","fdh4B","fdh4A","ftfL","fch","mtdA")


candidates <- c(efgA_order,formaldehyde_order)

# log2fT is a subset of all_data2 that has only candidate genes in formaldehyde treatment 
log2fT <- subset(all_data2,gene %in% candidates & time!='720' & treatment %in% c('form','pre'))
log2fT$counts <- as.numeric(as.character(log2fT$counts))
log2fT$time <-  as.numeric(as.character(log2fT$time))

# changing locus names to names of proteins 
for (i in 1:dim(log2fT)[1]){
  for (j in 1:length(candidates)){
    if (log2fT$gene[i]==candidates[j]) log2fT$gene[i] <- names(candidates)[j]
  }
}

## for plotting I break log2fT in two datasets
log2fT_efgA <- subset(log2fT,gene %in% names(efgA_order))
gene <- as.factor(log2fT_efgA$gene)
plot_efgA <- ggplot(data=log2fT_efgA, aes(x=time, y=counts, color=gene,group=factor(gene))) + geom_line(size=1) + geom_point(size=1.4) + facet_grid(~genotype)
plot_efgA <- plot_efgA + xlab('time (min)') + ylab('log2FC')
plot_efgA
write.csv(log2fT_efgA,file='log2ft_efgA.csv')

log2fT_form <- subset(log2fT,gene %in% names(formaldehyde_order))
gene <- as.factor(log2fT_form$gene)
plot_form <- ggplot(data=log2fT_form, aes(x=time, y=counts, color=gene,group=factor(gene))) + geom_line(size=1) + geom_point(size=1.4) + facet_grid(~genotype)
plot_form <- plot_form + xlab('time (min)') + ylab('log2FC')
plot_form
write.csv(log2fT_form,file='log2ft_form.csv')

# log2fT is a subset of all_data2 that has only candidate genes in formaldehyde treatment 
log2fT <- subset(all_data2,gene %in% form_growth_order & time!='720' & treatment %in% c('form','pre'))
log2fT$counts <- as.numeric(as.character(log2fT$counts))
log2fT$time <-  as.numeric(as.character(log2fT$time))

## list of formaldehyde growth genes  

form_growth_order <- paste("Mext_",c("0563","0564",1191,1206,
                                     1355,1356,1357,1358,1359,
                                     1360,1361,1362,1363,2114,
                                     2115,3498,3499,3500,3501,
                                     3502,3503,3504,3505,3506,
                                     3507,3508),sep = "")

colnames(gene_list2)[1] <- "gene"
all_data3 <- merge(all_data2,gene_list2)
colnames(all_data3)[5] <- "log2fc"
all_data3 <- subset(all_data3,select=-c(type))
all_data3 <- all_data3[,c(1,7,2:6)]

# looking at genes in formaldehyde treatement 
log2F_form_growth <- subset(all_data3,gene %in% form_growth_order & time!='720' & treatment %in% c('form','pre'))
write.csv(log2F_form_growth,file='log2f_form_growth.csv')

gene <- as.factor(log2F_form_growth$gene)
plot_form_growth <- ggplot(data=log2F_form_growth, aes(x=time, y=log2fc, color=gene,group=factor(gene))) + geom_line(size=1) + geom_point(size=1.4) + facet_grid(~genotype)
plot_form_growth <- plot_form_growth + xlab('time (min)') + ylab('log2FC')
plot_form_growth

# looking at genes in no-treatment 

log2F_form_growth_untreated <- subset(all_data3,gene %in% form_growth_order & time!='720' & treatment %in% c('none','pre'))
write.csv(log2F_form_growth_untreated,file='log2f_form_growth_untreated.csv')

gene <- as.factor(log2F_form_growth_untreated$gene)
plot_form_growth <- ggplot(data=log2F_form_growth_untreated, aes(x=time, y=log2fc, color=gene,group=factor(gene))) + geom_line(size=1) + geom_point(size=1.4) + facet_grid(~genotype)
plot_form_growth <- plot_form_growth + xlab('time (min)') + ylab('log2FC')
plot_form_growth