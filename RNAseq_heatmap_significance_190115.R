library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(made4)

# this script calculates log2fc with comparing each genotype to its pre condition and plot heatmaps

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

for (g in unique(all_data$genotype)){
  for (tr in unique(all_data$treatment)[-1]){
    for (ti in c('5','20','40','180','360')){
      if(dim(subset(all_data,genotype==g & treatment==tr & time==ti))[1]>0){
        sample <- subset(all_data,genotype==g & treatment==tr & time==ti)
        sample_short <- casting(sample)
        if (unique(sample_short$genotype)=='WT') {
          PRE <- wt_pre_short[,5:7]
          colnames(PRE) <- paste(unique(wt_pre$genotype),unique(wt_pre$time),sort(as.numeric(unique(wt_pre$replicate))),sep = '.')
        }
        if (unique(sample_short$genotype)=='efgA'){
          PRE <- e_pre_short[,5:7]
          colnames(PRE) <- paste(unique(e_pre$genotype),unique(e_pre$time),sort(as.numeric(unique(e_pre$replicate))),sep = '.')
        }
        
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
        res$log2FoldChange[is.na(res$padj)] <- 0 # if we want to plot only significant changes
        res$log2FoldChange[res$padj>=0.001] <- 0 # if we want to plot only significant changes 
        res$log2FoldChange[res$log2FoldChange <= 1 & res$log2FoldChange >= -1] <- 0 # 2X cutoff
        logf <- data.frame(res$log2FoldChange)
        colnames(logf) <- paste(unique(sample$genotype),unique(sample$treatment),unique(sample$time),sep = '.')
        significance <- cbind(significance,logf)
      }
      else {}
    }
  }
}

# adding the reference columns 
significance$WT.pre <- 0
significance$efgA.pre <- 0

# list of different treatments for plotting heatmaps
nostress <-  c("WT.pre","WT.none.5","WT.none.20","WT.none.40","WT.none.180","WT.none.360","efgA.pre","efgA.none.5","efgA.none.20","efgA.none.40","efgA.none.180","efgA.none.360")
formaldehyde <-  c("WT.pre","WT.form.5","WT.form.20","WT.form.40","WT.form.180","WT.form.360","efgA.pre","efgA.form.5","efgA.form.20","efgA.form.40","efgA.form.180","efgA.form.360")
kanamycin <-  c("WT.pre","WT.kan.40","WT.kan.180","WT.kan.360","efgA.pre","efgA.kan.40","efgA.kan.180","efgA.kan.360")
kan_form <- c("WT.pre","WT.kan.40","WT.kan.180","WT.kan.360","WT.form.5","WT.form.20","WT.form.40","WT.form.180","WT.form.360")

hmcol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# a function for plotting heatmaps 
plotting_heatmap <- function(columns){
  df <- data.matrix(subset(significance,select = columns))
  # dim 400-310
  # if dendrogram = "row" it shows dendrogram
  return(heatmap.2(df,col=hmcol,
                   scale="none", trace="none", dendrogram = "none",symbreaks = T, 
                   breaks = seq(-3,3,length.out = 101), 
                   keysize=2, srtCol=45, Colv=FALSE,labRow=F,cexCol = 1))
}

# dim 1000-500
plotting_heatmap(nostress)
plotting_heatmap(formaldehyde)
plotting_heatmap(kanamycin)
plotting_heatmap(kan_form)

# making another dataset from significance which doesn't have 0 values 
significance2 <- significance[,-1]
significance3 <- data.frame() # this dataset keep the values 
gene_list <- c() # to make a column of gene list that don't have 0 values

for (i in seq(1,length(significance$gene))){
  if (any(significance2[i,]!= 0) ){
    significance3 <- rbind(significance3,significance2[i,])
    gene_list = c(gene_list,i)
  }
}

# making sure there is no row with all 0 values 
for (i in seq(1,dim(significance3)[1])){
  if (all(significance3[i,] == 0)){
   print(i)
  }
}
significance3$"gene" <- significance$gene[gene_list]
significance3 <- significance3[,c(29,1:28)]


# making another dataset from significance which all L2FC are bigger < -1 or > 1 
significance2 <- significance[,-1]
significance3 <- data.frame() # this dataset keep the values 
gene_list <- c() # to make a column of gene list that don't have 0 values

for (i in seq(1,length(significance$gene))){
  if (any(significance2[i,]) > 1 | any(significance2[i,] < -1)){
    significance3 <- rbind(significance3,significance2[i,])
    gene_list = c(gene_list,i)
  }
}
significance3$"gene" <- significance$gene[gene_list]
significance3 <- significance3[,c(29,1:28)]

