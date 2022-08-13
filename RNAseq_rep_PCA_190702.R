library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
# this script put all the raw counts in one big matrix and plot PCA 
setwd("C:/Siavash/Codes/RNAseq/data")

# reading big raw data 
counts_all <- read.table("count.tsv")
counts_all <- as.data.frame(counts_all)

# passing sample id's
marx_list2 <- c("gene","3-7-E","2-3-W","4-1-E","10-1-E","3-1-E","4-3-W", "8-1-E","1-1-W","3-3-E","7-3-W","4-2-W","7-3-E","12-7-E","11-1-E",
                "11-3-W","3-7-W","11-2-E","7-2-W","10-9-E","12-9-E","10-3-E","7-1-E","1-5-W","3-9-E","11-1-W","4-1-W","11-3-E","10-7-W",
                "8-3-W","3-3-W","13-1-E","10-1-W","1-1-E","2-9-E","3-9-W","2-3-E","13-3-W","3-1-W","10-3-W","1-5-E","7-1-W","8-2-E","8-3-E",
                "4-2-E","2-9-W","2-7-E","13-1-W","10-7-E","12-9-W","13-3-E","8-1-W","7-2-E","10-9-W","12-7-W","4-3-E","11-2-W","2-1-W",
                "2-1-E", "8-2-W", "2-7-W")

# changing the column names to sample id's
colnames(counts_all) <- marx_list2

# removing the first row 
counts_all <- counts_all[-1,] 

# reading key id's of big data 
key_id_full <- read.csv("key_id_full2.csv")

# changing '.' to '-' in id's to merge with marx_list2
key_id_full$Sample <- gsub('.','-',key_id_full$Sample,fixed = T)

# converting counts_all to long format
counts_all_m <- melt(counts_all,id.vars  ="gene",variable.name = "Sample",value.name ="counts")
#counts_all_m <- gather(counts_all, gene, counts) # alternative 

# merging counts_all with id's 
counts_all_m2 <- merge(counts_all_m,key_id_full,by="Sample")

# now smaller dataset from the first Austin's run
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
colnames(all_data)[6] <- "Time"

# removing sample column 
all_data <- all_data[,-1]

# converting all_data to short format 
all_data_c <- dcast(all_data,gene~Genotype+Treatment+Time+Replicate,value.var = "counts")

# removing WT form 720 (it's only one replicate)
all_data_c <- subset(all_data_c,select=-WT_form_720_1)
###################
rownames(all_data_c) <- all_data_c$gene

# removing name of genes
all_data_c <- all_data_c[,-1]
all_data_c <- data.matrix(all_data_c)

# changing names to smaller names 
colnames(all_data_c) <- gsub('_','.',colnames(all_data_c))
colnames(all_data_c) <- gsub('WT','w',colnames(all_data_c))
colnames(all_data_c) <- gsub('efgA','e',colnames(all_data_c))
colnames(all_data_c) <- gsub('form','f',colnames(all_data_c))
colnames(all_data_c) <- gsub('kan','k',colnames(all_data_c))
colnames(all_data_c) <- gsub('.none','',colnames(all_data_c))
colnames(all_data_c) <- gsub('.pre.pre','.pre',colnames(all_data_c))
write.csv(as.data.frame(all_data_c),file='all_counts_rep.csv')

#makes colData
Condition_temp2 <- gsub('.{2}$', '', colnames(all_data_c))
colData2 <- data.frame(Condition_temp2)
rownames(colData2) <- colnames(all_data_c)
colData2$Condition_temp2 <- as.factor(colData2$Condition_temp2)

############################## PCA
mext_dds2 <- DESeqDataSetFromMatrix(countData = all_data_c, colData = colData2, design = ~ Condition_temp2)
vsd <- vst(mext_dds2, blind=FALSE)
plotPCA(vsd, intgroup="Condition_temp2")
pcaData <- plotPCA(vsd, intgroup="Condition_temp2",returnData=TRUE)
# dim 600-500
gg_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=group,label=group)) + geom_text(size=5) + theme(legend.position="none",text = element_text(size=17),axis.text = element_text(size = 17)) 
gg_pca <- gg_pca + xlab("PC1: 74% variance") + ylab("PC2: 11% variance") 
gg_pca

################################ fancy plotting
# a function to extract the time from id's
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

# extracting genotype 
pcaData$genotype <-  substr(pcaData$group,start = 1,stop = 1)

# extracting treatment 
pcaData$treatment <- substr(pcaData$group,start = 3,stop = 3)

# extracting time 
pcaData$time <- numextract(pcaData$group)
pcaData$time[is.na(pcaData$time)==TRUE] <- "pre"

pcaData$treatment[!(pcaData$treatment %in% c('f','k'))] <- "cnt"

# plotting with grouping based on treatment, genotype
gg_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=interaction(treatment,genotype),label=group)) + geom_text(size=5) #+ theme(legend.position="none",text = element_text(size=17),axis.text = element_text(size = 17)) 
gg_pca <- gg_pca + xlab("PC1: 74% variance") + ylab("PC2: 11% variance") 
#gg_pca

# for playing with shape and color 
pcaData$time[pcaData$time=='pre'] <- -45
pcaData$time <- as.numeric(pcaData$time)

gg_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, shape=interaction(genotype,treatment),color=as.factor(time),label=group)) + geom_point(size=5) #+ theme(legend.position="none",text = element_text(size=17),axis.text = element_text(size = 17)) 
gg_pca <- gg_pca + xlab("PC1: 74% variance") + ylab("PC2: 11% variance") + scale_color_manual(values=brewer.pal(length(unique(pcaData$time))+2, 'YlOrRd'), name='time (hours)')
#gg_pca

############ very fancy plotting 
wcnt <- subset(pcaData,genotype=='w' & treatment == 'cnt') 
wcnt <- wcnt[order(wcnt$time),]

wform <- subset(pcaData,genotype=='w' & treatment == 'f')
wform <- wform[order(wform$time),]

wkan <- subset(pcaData,genotype=='w' & treatment == 'k')
wkan <- wkan[order(wkan$time),]

ecnt <- subset(pcaData,genotype=='e' & treatment == 'cnt') 
ecnt <- ecnt[order(ecnt$time),]

eform <- subset(pcaData,genotype=='e' & treatment == 'f')
eform <- eform[order(eform$time),]

ekan <- subset(pcaData,genotype=='e' & treatment == 'k')
ekan <- ekan[order(ekan$time),]

cnt <- rbind(wcnt,ecnt)
cnt <- cnt[order(cnt$time),]

form <- rbind(wform,eform)
form <- form[order(form$time),]

kan <- rbind(wkan,ekan)
kan <- kan[order(kan$time),]

pcaData3 <- rbind(cnt,form,kan)

# we have 6 time points in control 
ccnt <- colorRampPalette(c('white','black')) # blue
ccnt <- ccnt(7)
ccnt <- ccnt[-1]

# need 3 time points in kanamyicn 
ckan <- colorRampPalette(c('white','#08306b')) #blue
ckan <- ckan(4)
ckan <- ckan[-1]

# need 5 time points for formaldehyde 
cform <- colorRampPalette(c('white','red')) 
cform <- cform(6)
cform <- cform[-1]

# I need a factor that has both treatment and time 
pcaData3$treatment[pcaData3$treatment=='cnt'] <- 'untr'
pcaData3$treatment[pcaData3$treatment=='f'] <- 'form'
pcaData3$treatment[pcaData3$treatment=='k'] <- 'kan'
pcaData3$trTime <- paste(pcaData3$treatment,pcaData3$time,sep = '.')
Treatment <- factor(pcaData3$trTime,levels = unique(pcaData3$trTime))
Genotype <- factor(pcaData3$genotype,levels = unique(pcaData3$genotype))

# dim: 600-450                     
gg_pca <- ggplot(pcaData3, aes(x=PC1, y=PC2, color=Treatment, shape=Genotype)) + geom_point(size=5,alpha=0.85)  
gg_pca <- gg_pca + xlab("PC1: 74% variance") + ylab("PC2: 11% variance") + scale_color_manual(values = c(ccnt,cform,ckan))
gg_pca

