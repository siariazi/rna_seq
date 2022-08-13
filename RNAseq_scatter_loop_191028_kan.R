library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggrepel)
library(reshape2)

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
#df <- counts_all[1:20,]

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
colnames(all_data)[6] <- "Time"
#all_data$ counts <- as.numeric(all_data$counts)
all_data <- all_data[,-1]
all_data_c <- dcast(all_data,gene~Genotype+Treatment+Time+Replicate,value.var = "counts")
###################
# remove name of genes
rownames(all_data_c) <- all_data_c$gene
all_data_c <- all_data_c[,-1] 
all_data_c <- data.matrix(all_data_c)
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

#normalized counts, adds gene descriptions, saves
mext_dds2 <- DESeqDataSetFromMatrix(countData = all_data_c, colData = colData2, design = ~ Condition_temp2)
mext_dds_counts2 <- estimateSizeFactors(mext_dds2)
normalized_counts2 <- counts(mext_dds_counts2, normalized = T)
normalized_all_data_c <- data.frame(normalized_counts2)
write.csv(normalized_all_data_c, file="all_normalized_counts_rep.csv")

############## mean of replicates 
normalized_mean <- data.frame(gene=counts_all$gene)
j = 1
for (i in unique(colData2$Condition_temp2)[-18]){ # -18 takes w.f.720 out of the mean
  normalized_mean$i <- apply(normalized_all_data_c[,colnames(normalized_all_data_c)[gsub('.{2}$', '', colnames(normalized_all_data_c))==i]],1,mean)
  j = j +1
  colnames(normalized_mean)[j] <- i
}
#normalized_mean <- normalized_mean[,-19]

normalized_mean$w.f.720 <- all_data_c[,50]
normalized_mean <- normalized_mean[,c("gene","w.pre","w.5","w.40","w.180","w.360","e.pre","e.5","e.20","e.40",
                                      "e.180","e.360","w.f.5","w.f.40","w.f.180","w.f.360","w.f.720","e.f.5","e.f.20","e.f.40","e.f.180",
                                      "e.f.360","w.k.40","w.k.180","w.k.360","e.k.40","e.k.180","e.k.360")]


write.csv(normalized_mean,file = 'normalized_mean.csv')

mext_gene_trans <- read_csv(file="C:/Siavash/Codes/RNAseq/data/mext_desc.csv")
mext_ddsres <- DESeq(mext_dds2)
mext_results <- results(mext_ddsres) #choose conditions to compare from colData
mext_results_df <- as_data_frame(mext_results)
mext_results_df <- cbind(feat_id = rownames(mext_results_df), mext_results_df) #adds mext IDs
rownames(mext_results_df) <- NULL #removes rownames
mext_results_df[,c(1,2,4,5,6,7)] <- NULL #removes statistics
mext_results_df <- cbind(mext_gene_trans,mext_results_df)
colnames(normalized_mean)[1] <- "feat_id"
normalized_mean$feat_id <- as.character(normalized_mean$feat_id) # feat_id is factor

mext_results_df_original <- mext_results_df # because I'm using mext_results_df in the loop I need to keep the orig
A <- "w.f.720" # depends on the condition choose 
A <- "w.pre" # if you run this, run line 112
A <- "e.pre" # if you run this, run line 113

plotting_data <- data.frame()

B_names <- c("w.pre","w.40","w.180","w.360","e.pre","e.40","e.180","e.360", # depends on the condition 
             "w.f.40","w.f.180","w.f.360","e.f.40","e.f.180","e.f.360")
B_names <- c("w.k.40","w.k.180","w.k.360") 
B_names <- c("e.k.40","e.k.180","e.k.360")

for (B in B_names){
  mext_results_df <- mext_results_df_original
  norm_subset_2 <- select(normalized_mean,1,A,B)
  norm_subset_2$mean_counts <- (norm_subset_2[,A] + norm_subset_2[,B]) / 2 
  mext_results_df <- left_join(mext_results_df,norm_subset_2, by = "feat_id")
  mext_results_df <- filter(mext_results_df, mean_counts > 10) #removes genes whose mean counts for the two strains of interest are below 10.
  
  mext_A_v_B_down <- filter(mext_results_df, log2FoldChange < 0)
  file_name <- gsub('.','_',paste(A,"V",B,sep = '.'),fixed = TRUE)
  #write_csv(mext_A_v_B_down, path=paste(file_name,"_down.csV",sep = '')) #change for each
  mext_A_v_B_up <- filter(mext_results_df, log2FoldChange > 0)
  mext_A_v_B_up <- arrange(mext_A_v_B_up, desc(log2FoldChange))
  #write_csv(mext_A_v_B_up, path=paste(file_name,"_up.csV",sep = '')) #change for each
  
  #plot of log2 normalized counts for the two samples compared above. 
  norm_subset <- select(normalized_mean, 1, A, B) 
  mext_ids_large <- filter(mext_results_df, log2FoldChange > 0.1 | log2FoldChange < -0.1) %>%
    select(feat_id)
  mext_ids_large_list <- mext_ids_large$feat_id
  norm_subset_large <- filter(norm_subset, feat_id %in% mext_ids_large_list)
  norm_subset_large$mean_counts <- (norm_subset_large[,A] + norm_subset_large[,B]) / 2 
  norm_subset_large <- filter(norm_subset_large, mean_counts > 10)
  norm_subset_large$mean_counts <- NULL
  norm_subset_large$large <- 1
  norm_subset_large_10 <- as.factor(norm_subset_large$feat_id)
  norm_subset_small <- filter(norm_subset, ! feat_id %in% norm_subset_large_10)
  norm_subset_small$large <- 0
  norm_subset_final <- bind_rows(norm_subset_large, norm_subset_small)
  cols <- c(A,B) 
  norm_subset_final[cols] <- log2(norm_subset_final[cols]) 
  norm_subset_final$"variation" <- B
  colnames(norm_subset_final)[3] <- "treatment"
  plotting_data <- rbind(plotting_data,norm_subset_final)
}

write.csv(plotting_data,file='plot_norm_subset.csv')
plotting_data <- read.csv('plot_norm_subset.csv')

plotting_data$variation <- factor(plotting_data$variation, levels=B_names) # re-ordering the factors

#large <- as.factor(plotting_data$large)
plot_norm_subset <- ggplot(plotting_data, aes_string(x = A, y="treatment", colour = "large")) + 
  geom_point() + facet_wrap(~variation,nrow = 1) + theme(legend.position="none") + labs(y="") 

print(plot_norm_subset)
for_save <- plot_grid(plot_norm_subset)
#save_plot(filename = paste(file_name,".pdf",sep = ''), for_save, base_width = 10, base_height = 6) #change for each

