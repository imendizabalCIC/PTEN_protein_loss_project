################################################################################
################################### xCELL ######################################
################################################################################

#| Code for inference the score of different cell types using xCell.

################################################################################

########################### LIBRARIES AND DATA #################################
suppressMessages(library(xCell))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(VennDiagram))
suppressMessages(library(RColorBrewer))
suppressMessages(library(survival))  
suppressMessages(library(tidyestimate))
suppressMessages(library(survminer)) 
suppressMessages(library(ggfortify))
suppressMessages(library(colorspace))

theme_set(theme_classic())

#| Loading data expression and sample info
exprMatrix <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt",header=TRUE,row.names=1, as.is=TRUE)
data_trait <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt", sep ="\t", header=T)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
################################################################################


################################################################################
#| PROCESSING DATA EXPRESSION PCa PTEN IHC 
################################################################################

#| NOTE: xCells performs single-sample gene set enrichment analysis on the expression
#| data to evaluate the enrichment of 64 immune and stroma cell types

#| Filtering out NA H-score values
data_trait <- data_trait[which(!is.na(data_trait$H_score)),]
exprMatrix <- exprMatrix[, data_trait$AC.basurto]

#| Filtering low counts
filter_counts <- rowSums(exprMatrix>5) >= 0.7*ncol(exprMatrix)
exprMatrix <- exprMatrix[filter_counts,]

#| Converting rownames to HGNC symbol (REQUIREMENT OF XCELL)
exprMatrix$GeneID <- rownames(exprMatrix)
exprMatrix <- merge(exprMatrix, genome_GRCh39.94, by ="GeneID")
exprMatrix <- aggregate(exprMatrix, by=list(exprMatrix$gene_name), mean)
rownames(exprMatrix) <- exprMatrix$Group.1

gene.length <- exprMatrix$Length
exprMatrix <- exprMatrix[,data_trait$AC.basurto]

#| Compurting TPM
x <- exprMatrix / gene.length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )

#| Running xCell
data <-xCellAnalysis(tpm.mat)

#| Confirming that are in the same order
any(colnames(data) == rownames(data_trait))

#| Adding the info from data to the data_trait
data_trait <- cbind(data_trait, t(data))

#| Renaming some columns 
data_trait$PTEN_status <- data_trait$H_score_cut_0
data_trait$PTEN_mrna <- data_trait$PTEN_Exp_log2

#| save data
write.table(data_trait, "Results/Sample_info_table/Xcell_output_PTEN_IHC_147.txt", sep="\t", row.names=T)

#| PTEN status and mRNA levels
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = PTEN_mrna, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="serif", size =5)+
  theme(text=element_text(size=18,  family="serif"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("PTEN expression \n(Normalized log2 counts)") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/PTEN_mRNA_PTEN_protein_presence_vs_intact.pdf", heigh= 3.9, width =4.5)
################################################################################


################################################################################
#| PTEN presence VS INTACT 
################################################################################

data_trait$PTEN_status[which(data_trait$PTEN_status == "PTEN presence")] <- "presence"
data_trait$PTEN_status[which(data_trait$PTEN_status == "PTEN intact")] <- "presence"

#| StromaScore
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = StromaScore, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="sans", size =5)+
  theme(text=element_text(size=18,  family="sans"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("Stromal Score") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Stroma_Score_Results_PTEN_presence_vs_intact.pdf", heigh= 3.9, width =3.1)

#| ImmuneScore
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = ImmuneScore, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="sans", size =5)+
  theme(text=element_text(size=18,  family="sans"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("Immune Score") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Immune_Score_Results_PTEN_presence_vs_intact.pdf",heigh= 3.9, width =3.1)

#| MicroenvironmentScore
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = MicroenvironmentScore, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="sans", size =5)+
  theme(text=element_text(size=18,  family="sans"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("Microenvironment Score") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Microenvironment_score_Results_PTEN_presence_vs_intact.pdf", heigh= 3.9, width =3.1)

#| Fibroblasts
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = Fibroblasts, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="sans", size =5)+
  theme(text=element_text(size=18,  family="sans"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("Fibroblasts") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Fibroblasts_Results_PTEN_presence_vs_intact.pdf", heigh= 3.9, width =3.1)

#| Epithelial cells
ggplot(data_trait[which(!is.na(data_trait$PTEN_status)),], aes(x= PTEN_status, y = `Epithelial cells`, fill =PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),
             position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),
             pch = 21, size=2.4, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=0.5,hjust= -0.5,family="sans", size =5)+
  theme(text=element_text(size=18,  family="sans"), 
        legend.key.size = unit(2, 'cm'),legend.position = "none", 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab("PTEN protein status")+
  ylab("Epithelial cells") +
  scale_fill_manual(values =c("#EFD5B2", "#504667"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Box_plots/xCell_Epithelial_cells_Results_PTEN_presence_vs_intact.pdf", heigh= 3.9, width =3.1)
################################################################################


################################################################################
#| DIFFERENCES BY CELL TYPE UPON PTEN PROTEIN presence
################################################################################
h_data <- data_trait[which(!is.na(data_trait$PTEN_status)),]
h_data <- h_data[,c(39:(39+67))]
h_data
p_value <-  c()
difference <- c()
for (i in 1:(dim(h_data)[2] - 1) ){
  
  wilcox_test <- wilcox.test(h_data[which(h_data$PTEN_status == "presence"), i], h_data[which(h_data$PTEN_status == "intact"), i])
  p_value <- c(p_value, wilcox_test$p.value )
  difference <- c(difference, mean(h_data[which(h_data$PTEN_status == "presence"), i]) - mean(h_data[which(h_data$PTEN_status == "intact"), i]))
  
}

h_k <- data.frame( datasets = colnames(h_data[,-c(68)]),
                   p_value =p_value,
                   difference =difference)


h_k <- h_k[order(h_k$p_value, decreasing = T),]
h_k$log <- -log(h_k$p_value,2)
h_k$difference[which(h_k$difference >= 0)] <- "Positive"
h_k$difference[which(h_k$difference < 0)] <- "Negative"

ggplot(h_k[56:67,], aes(x = log, y = reorder(datasets,log), fill =difference)) +
  geom_bar(stat = "identity", color ="black") + 
  scale_fill_manual(values = c("midnightblue", "darkgoldenrod3")) +
  scale_fill_manual(values = c("black", "white")) +
  theme(text=element_text(size=16,  family="sans"), axis.text.y = element_text(size = 14, face ="bold", color ="black"),
        axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("") +
  xlab("-log2(p-value)") +
  labs(fill = "PTEN protein\npresence - intact")
ggsave("W:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/xCell_Results_PTEN_presence_vs_intact_wilxon_test_last.pdf", heigh= 5.5, width =7.5)


ggplot(h_k[56:67,], aes(x =log , y =reorder(datasets,log), color =difference)) +
  geom_point(size =5)+
  #scale_color_continuous_sequential(palette = "ag_GrnYl")+
  theme(text=element_text(size=16,  family="sans"), 
        axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("") +
  xlab("-log2(p value)")+
  labs(color ="Mean difference") +
  scale_color_manual(values =c("peachpuff3", "black"))
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Bubble_Plots/xCell_PTEN_presence_vs_intact_bubble_last.pdf", heigh =4, width = 7)
################################################################################


################################################################################

h_data <- data_trait

pi3k <- h_data[,c("Mean.expression.PI3K.AKT.mTOR")]

h_data <- h_data[,c(39:(39+66))]

p_value <-  c()
r <- c()

for (i in 1:(dim(h_data)[2]) ){
  
  cor <- cor.test(pi3k, h_data[,c(i)], method="spearman")
  p_value <- c(p_value,cor$p.value)
  r <- c(r, cor$estimate[[1]])
  
}

data_frame <- data.frame(celltypes = colnames(h_data),
                         p_value = p_value,
                         r =r)
data_frame <- data_frame[order(data_frame$p_value),]

data_frame_positive <- data_frame[which(data_frame$r >0),]
data_frame_negative <- data_frame[which(data_frame$r <0),]

ggplot(data_frame_positive[c(1:10),], aes(x =p_value , y =reorder(celltypes,p_value))) +
  geom_point(aes(size =r))+
  #scale_color_continuous_sequential(palette = "ag_GrnYl")+
  theme(text=element_text(size=16,  family="sans"), 
        axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("") +
  xlab("-log2(p value)")+
  labs(color ="Mean difference") 
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Bubble_Plots/xCell_Correlation_PI3K-AKT-mTOR_celltype_PCa_samples_Basurto_Positive_Correlation.pdf", heigh =4, width = 6)

ggplot(data_frame_negative[c(1:10),], aes(x =p_value , y =reorder(celltypes,p_value))) +
  geom_point(aes(size =r))+
  #scale_color_continuous_sequential(palette = "ag_GrnYl")+
  theme(text=element_text(size=16,  family="sans"), 
        axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("") +
  xlab("-log2(p value)")+
  labs(color ="Mean difference")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Bubble_Plots/xCell_Correlation_PI3K-AKT-mTOR_celltype_PCa_samples_Basurto_Negative_Correlation.pdf", heigh =4, width = 6)

################################################################################