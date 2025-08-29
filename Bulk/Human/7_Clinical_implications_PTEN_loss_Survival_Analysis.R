################################################################################
#| Survival analysis on transcriptional signatures
################################################################################
#| Date: 29/08/2025
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script performs survival analysis on transcriptional signatures and 
#| WGCNA module scores (particularly the green module) in prostate cancer 
#| cohorts. It integrates Basurto (H-score/PTEN protein-based) and TCGA 
#| (genomic) datasets to assess prognostic relevance.  
#|
#| Workflow:
#|   1) Load and normalize RNA-seq data (Basurto and TCGA), compute log2 counts 
#|      and z-scores, and merge with WGCNA module and DEG annotations.  
#|   2) Stratify patients into high/low expression groups (quartiles Q1 vs Q4).  
#|   3) Perform Cox proportional hazards regression and Kaplan-Meier analysis 
#|      for:
#|        * Basurto cohort (overall and PTEN loss vs presence).  
#|        * TCGA cohort (overall and PTEN CNA intact vs loss).  
#|   4) Evaluate combined effects of PTEN status and module/signature scores.  
#|   5) Generate survival plots, risk tables, and hazard ratio estimates.  
#|
#| Outputs:
#|   - Kaplan-Meier survival curves (Basurto and TCGA, Q1 vs Q4, stratified by PTEN status and green module scores).  
#|   - Cox model hazard ratios with p-values.  
#|   - Risk tables and confidence intervals (ggsurvfit).  
#|   - Distribution plots of module/signature scores.
################################################################################

################################################################################
#| LIBRARY AND DATA
################################################################################
suppressMessages(library(dplyr))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(viridis))
suppressMessages(library(readxl))
suppressMessages(library(colorspace))
suppressMessages(library(ggpubr))
suppressMessages(library(DESeq2))
suppressMessages(library(glmnet))
suppressMessages(library(e1071))
suppressMessages(library(DBI))
suppressMessages(library(RMySQL))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))
suppressMessages(library(ggcorrplot))
suppressMessages(library(caret))
suppressMessages(library(randomForest))
suppressMessages(library(gam))
suppressMessages(library(tidyverse))
suppressMessages(library(xgboost))
suppressMessages(library(pheatmap))
suppressMessages(library(reshape2))
suppressMessages(library(survival))  
suppressMessages(library(survminer))
suppressMessages(library(ggfortify))

#| Setting working directory
setwd("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Survival_Analysis/")

#| For plots
theme_set(theme_classic())

#| Data pathways 
data.file_WGCNA <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_intramodularinfo_last_try.txt"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Tables/Limma_voom_analysis_DEGS_results.txt"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Images/11_MachineLearning_Regression_Models_to_identify_candidates/"
sample_info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_ERG_Fusion.txt"
tcga.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_merge.txt"

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Sample info with pAKT info
sample_info <- read.table(sample_info.file, sep = "\t", header =T)
sample_info <- sample_info[which(sample_info$Diagnostico == "PCa"),]
sample_info <- sample_info[order(rownames(sample_info)),]
sample_info <- sample_info[which(! (rownames(sample_info) == "AC92")),]

#| Counts data Basurto ALL PCa
counts_data <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt", header=T) 
counts_data <- counts_data[,rownames(sample_info)]
filter_counts <- rowSums(counts_data>5) >= 0.7*ncol(counts_data)
counts_data <- counts_data[filter_counts,]
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ 1)
dds <- estimateSizeFactors(dds)

counts_data_normalized <- counts(dds, normalized=TRUE)
counts_data_normalized <- log(counts_data_normalized + 1, base =2)
counts_data_normalized <- data.frame(counts_data_normalized)
counts_data_normalized$GeneID <- rownames(counts_data_normalized)
counts_data_normalized <- merge(counts_data_normalized, genome_GRCh39.94, by ="GeneID")
counts_data_normalized <- aggregate(counts_data_normalized, by =list(counts_data_normalized$gene_name), mean)
rownames(counts_data_normalized) <- counts_data_normalized$Group.1
counts_data_normalized <- counts_data_normalized[,rownames(sample_info)]

counts_zscore <- scale(t(as.data.frame(counts_data_normalized)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)
counts_zscore <- as.data.frame(counts_zscore)

#| WGCNA data
wgcna_Modules <- read.table(data.file_WGCNA, sep = "\t")
module_colors <- unique(wgcna_Modules$moduleColors)
wgcna_Modules <- merge(wgcna_Modules, genome_GRCh39.94, by ="GeneID")
rownames(wgcna_Modules) <- wgcna_Modules$GeneID

#| DEGs data (Limma Voom)
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")
DEGs <- DEGs_genes[which((DEGs_genes$log2FC_limma_voom > 0 | DEGs_genes$log2FC_limma_voom < 0) & DEGs_genes$padj_limma_voom < 0.05 & !is.na(DEGs_genes$padj_limma_voom)),]

#| Changing name intact by "Presence"
sample_info$PTEN_status <- sample_info$H_score_cut_0
sample_info$PTEN_status[which(sample_info$PTEN_status == "PTEN intact")] <- "PTEN presence"
sample_info$PTEN_status[which(sample_info$PTEN_status == "PTEN presence")] <- "presence"
sample_info$PTEN_status[which(sample_info$PTEN_status == "PTEN loss")] <- "loss"
sample_info$PTEN_status <- factor(sample_info$PTEN_status, levels = c("presence", "loss"))

#| ADDING MODULE INFORMATION
genes_mod_purple <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "purple")], DEGs$gene_name)
genes_mod_green <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)
samples_media_expression_purple <- colMeans(counts_zscore[genes_mod_purple,], na.rm=TRUE)
samples_media_expression_green <- colMeans(counts_zscore[genes_mod_green,], na.rm=TRUE)

sample_info$mod_purple <- samples_media_expression_purple
sample_info$mod_green <- samples_media_expression_green

sample_info$PTEN_protein_status<- as.character(sample_info$PTEN_status)
sample_info$PTEN_Exp_log2 <- as.numeric(scale(sample_info$PTEN_Exp_log2))
sample_info$PTEN_protein_status <- as.numeric(scale(sample_info$H_score))
sample_info$mod_purple_zscore <- as.numeric(scale(samples_media_expression_purple))
sample_info$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| TCGA data
###| TCGA DATA
#| Loading TCGA data already processed
TCGA_data <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_PTEN_CNA_complete.txt", sep ="\t")
TCGA_xcell <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/xCell/xCell_results_TPM_all.txt", sep ="\t")

dim(TCGA_data)
#| Sample ID
TCGA_data$SampleID <- TCGA_data$PATIENT_ID
TCGA_xcell$SampleID <- rownames(TCGA_xcell)

#| Merging these dataframes
TCGA_data <- merge(TCGA_xcell,TCGA_data, by ="SampleID" )

#| Changing values of DFS states
TCGA_data$DFS_STATUS[which(TCGA_data$DFS_STATUS== "0:DiseaseFree")] <- "0"
TCGA_data$DFS_STATUS[which(TCGA_data$DFS_STATUS== "1:Recurred/Progressed")] <- "1"
TCGA_data$DFS_STATUS[which(TCGA_data$DFS_STATUS== "[Not Available]")] <- "NA"
TCGA_data$DFS_MONTHS <- as.numeric(TCGA_data$DFS_MONTHS)
TCGA_data$DFS_STATUS <- as.numeric(TCGA_data$DFS_STATUS)

#| Gleason score
TCGA_data$GLEASON_SCORE <- as.character(TCGA_data$GLEASON_SCORE)

#| Standardizing stromal score
TCGA_data$StromaScore_Zscore <- (TCGA_data$StromaScore - mean(TCGA_data$StromaScore))/sd(TCGA_data$StromaScore)

#| Loading normalized data
counts_data_tcga <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/NormalizedCounts/TCGA_raw_counts_filtered_normalized_log2.txt", sep ="\t")

#| Computing z-scores
counts_zscore_tcga <- scale(t(as.data.frame(counts_data_tcga)))
counts_zscore_tcga <- as.data.frame(counts_zscore_tcga)
counts_zscore_tcga <- t(counts_zscore_tcga)
counts_zscore_tcga <- as.data.frame(counts_zscore_tcga)

#| Taking only samples ID in data info
counts_zscore_tcga <- counts_zscore_tcga[, TCGA_data$SampleID]

#| Confirming same sample ID in sample info and counts
TCGA_data$SampleID == colnames(counts_zscore_tcga)

#| Obtaining GeneIDs from modules
genes_mod_purple <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "purple")], DEGs$GeneID)
genes_mod_green <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)

#| Computing sample mean expression in tcga data 
samples_media_expression_purple <- colMeans(counts_zscore_tcga[genes_mod_purple,], na.rm=TRUE)
samples_media_expression_green <- colMeans(counts_zscore_tcga[genes_mod_green,], na.rm=TRUE)

TCGA_data$mod_purple<- samples_media_expression_purple
TCGA_data$mod_green <- samples_media_expression_green

TCGA_data$mod_purple_zscore <- as.numeric(scale(samples_media_expression_purple))
TCGA_data$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

TCGA_data$PTEN_status <- NA
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "-1")] <- "0"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "2")] <- "0"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "1")] <- "0"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "NP")] <- "0"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "0")] <- "0"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_cna == "-2")] <- "-2"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_status == "-2")] <- "loss"
TCGA_data$PTEN_status[which(TCGA_data$PTEN_status == "0")] <- "presence"
TCGA_data$PTEN_status <- factor(TCGA_data$PTEN_status, levels = c("presence", "loss"))

#| Quantiles of PTEN protein
TCGA_data$PTEN_protein_quantile <- cut(as.numeric(TCGA_data$PTEN_protein),breaks=quantile(as.numeric(TCGA_data$PTEN_protein),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
TCGA_data$PTEN_protein_quantile <- as.character(TCGA_data$PTEN_protein_quantile)
TCGA_data$PTEN_protein_quantile <- factor(TCGA_data$PTEN_protein_quantile, levels = c("4", "3", "2", "1"))

TCGA_data$PTEN_protein_quantile_join <- as.character(TCGA_data$PTEN_protein_quantile)
TCGA_data$PTEN_protein_quantile_join[which(TCGA_data$PTEN_protein_quantile_join == "2")] <- "1"
TCGA_data$PTEN_protein_quantile_join[which(TCGA_data$PTEN_protein_quantile_join == "3")] <- "2"
TCGA_data$PTEN_protein_quantile_join[which(TCGA_data$PTEN_protein_quantile_join == "4")] <- "4,3"
TCGA_data$PTEN_protein_quantile_join[which(TCGA_data$PTEN_protein_quantile_join == "1")] <- "2,1"
TCGA_data$PTEN_protein_quantile_join <- factor(TCGA_data$PTEN_protein_quantile_join, levels = c("4,3", "2,1"))

TCGA_data$PTEN_protein_quantile_percentiles <- cut(TCGA_data$PTEN_protein, 
                                                   breaks = quantile(as.numeric(TCGA_data$PTEN_protein), probs = seq(0, 1, by = 0.1), na.rm = TRUE), 
                                                   include.lowest = TRUE, 
                                                   labels = 1:10)

TCGA_data$PTEN_protein_quantile_percentiles_num <- as.numeric(TCGA_data$PTEN_protein_quantile_percentiles)

levels(TCGA_data$PTEN_protein_quantile_percentiles) <- rev(levels(TCGA_data$PTEN_protein_quantile_percentiles))
levels(TCGA_data$PTEN_protein_quantile_percentiles) <- c("100%", "90%", "80%", "70%", "60%" , "50%", "40%", "30%" , "20%", "10%")

TCGA_data$PTEN_protein_quantile_percentiles_num <- as.numeric(TCGA_data$PTEN_protein_quantile_percentiles)
TCGA_data$PTEN_Exp_log2 <- as.numeric(scale(as.numeric(counts_data_tcga[which(rownames(counts_data_tcga) == "ENSG00000171862"),])))
TCGA_data$PTEN_protein_status <- scale(TCGA_data$PTEN_protein)
TCGA_data$PTEN_protein_status <- as.numeric(TCGA_data$PTEN_protein_status)
TCGA_data$Gleason_score_pieza <- as.numeric(TCGA_data$GLEASON_SCORE)
TCGA_data$DFS.STATUS <- TCGA_data$DFS_STATUS
################################################################################


################################################################################
#| ADDING SIGNATURES EVALUATED
################################################################################

signature_name <- c("TGFB_ACTIVATION", 
                    "F-TBRS", 
                    "T-TBRS", 
                    "Ma-TBRS", 
                    "End-TBRS",
                    "REACTOME_SASP",
                    "SASP_ALIMONTI",
                    "SASP_p53")

#| Signatures
TGFB_ACTIVATION <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[1],".txt", sep =""), sep =" ", header =T)
TGFB_ACTIVATION <- names(TGFB_ACTIVATION)

F_TBRS <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[2],".txt", sep =""), sep =" ", header =T)
F_TBRS <- names(F_TBRS)

T_TBRS <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[3],".txt", sep =""), sep =" ", header =T)
T_TBRS <- names(T_TBRS)

Ma_TBRS <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[4],".txt", sep =""), sep =" ", header =T)
Ma_TBRS <- names(Ma_TBRS)

End_TBRS <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[5],".txt", sep =""), sep =" ", header =T)
End_TBRS <- names(End_TBRS)

REACTOME_SASP <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[6],".txt", sep =""), sep =" ", header =T)
REACTOME_SASP <- names(REACTOME_SASP)

SASP_ALIMONTI <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[7],".txt", sep =""), sep =" ", header =T)
SASP_ALIMONTI <- names(SASP_ALIMONTI)

SASP_p53 <- read.table(paste("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Data/",signature_name[8],".txt", sep =""), sep =" ", header =T)
SASP_p53 <- names(SASP_p53)


####| BASURTO
sample_info$TGFB_ACTIVATION <- colMeans(counts_zscore[TGFB_ACTIVATION,], na.rm=TRUE)
sample_info$F_TBRS <- colMeans(counts_zscore[F_TBRS,], na.rm=TRUE)
sample_info$T_TBRS <- colMeans(counts_zscore[T_TBRS,], na.rm=TRUE)
sample_info$Ma_TBRS <- colMeans(counts_zscore[Ma_TBRS,], na.rm=TRUE)
sample_info$End_TBRS <- colMeans(counts_zscore[End_TBRS,], na.rm=TRUE)
sample_info$REACTOME_SASP <- colMeans(counts_zscore[REACTOME_SASP,], na.rm=TRUE)
sample_info$SASP_ALIMONTI <- colMeans(counts_zscore[SASP_ALIMONTI,], na.rm=TRUE)
sample_info$SASP_p53 <- colMeans(counts_zscore[SASP_p53,], na.rm=TRUE)

sample_info$TGFB_ACTIVATION_zscore <- as.numeric(scale(colMeans(counts_zscore[TGFB_ACTIVATION,], na.rm=TRUE)))
sample_info$F_TBRS_zscore <- as.numeric(scale(colMeans(counts_zscore[F_TBRS,], na.rm=TRUE)))
sample_info$T_TBRS_zscore <- as.numeric(scale(colMeans(counts_zscore[T_TBRS,], na.rm=TRUE)))
sample_info$Ma_TBRS_zscore <- as.numeric(scale(colMeans(counts_zscore[Ma_TBRS,], na.rm=TRUE)))
sample_info$End_TBRS_zscore <- as.numeric(scale(colMeans(counts_zscore[End_TBRS,], na.rm=TRUE)))
sample_info$REACTOME_SASP_zscore <- as.numeric(scale(colMeans(counts_zscore[REACTOME_SASP,], na.rm=TRUE)))
sample_info$SASP_ALIMONTI_zscore <- as.numeric(scale(colMeans(counts_zscore[SASP_ALIMONTI,], na.rm=TRUE)))
sample_info$SASP_p53_zscore <- as.numeric(scale(colMeans(counts_zscore[SASP_p53,], na.rm=TRUE)))

####| TCGA
TCGA_data$TGFB_ACTIVATION <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% TGFB_ACTIVATION)],], na.rm=TRUE)
TCGA_data$F_TBRS <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% F_TBRS)],], na.rm=TRUE)
TCGA_data$T_TBRS <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% T_TBRS)],], na.rm=TRUE)
TCGA_data$Ma_TBRS <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% Ma_TBRS)],], na.rm=TRUE)
TCGA_data$End_TBRS <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% End_TBRS)],], na.rm=TRUE)
TCGA_data$SASP_ALIMONTI <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% SASP_ALIMONTI)],], na.rm=TRUE)
TCGA_data$End_TBRS <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% End_TBRS)],], na.rm=TRUE)
TCGA_data$SASP_p53 <-  colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% SASP_p53)],], na.rm=TRUE)

TCGA_data$TGFB_ACTIVATION_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% TGFB_ACTIVATION)],], na.rm=TRUE)))
TCGA_data$F_TBRS_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% F_TBRS)],], na.rm=TRUE)))
TCGA_data$T_TBRS_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% T_TBRS)],], na.rm=TRUE)))
TCGA_data$Ma_TBRS_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% Ma_TBRS)],], na.rm=TRUE)))
TCGA_data$End_TBRS_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% End_TBRS)],], na.rm=TRUE)))
TCGA_data$SASP_ALIMONTI_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% SASP_ALIMONTI)],], na.rm=TRUE)))
TCGA_data$End_TBRS_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% End_TBRS)],], na.rm=TRUE)))
TCGA_data$SASP_p53_zscore <-  as.numeric(scale(colMeans(counts_zscore_tcga[genome_GRCh39.94$GeneID[which(genome_GRCh39.94$gene_name %in% SASP_p53)],], na.rm=TRUE)))

#| Adding module green
genes_mod_green <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)
counts<- counts_zscore_tcga[,which(colnames(counts_zscore_tcga) %in% data_DFS$SampleID)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
data_DFS$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| Getting the quartiles
data_DFS$quantile <- cut(as.numeric(data_DFS$mod_green_zscore),breaks=quantile(as.numeric(data_DFS$mod_green_zscore),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_DFS$quantile_joint<-NA
data_DFS$quantile_joint[which( (data_DFS$quantile == "1") | (data_DFS$quantile == "2"))]<-"low"
data_DFS$quantile_joint[which( (data_DFS$quantile == "4") | (data_DFS$quantile == "3"))]<-"high"

#| Separating cases
data_DFS$Cases <- NA
data_DFS$Cases[which(data_DFS$PTEN_cna == "0")] <- "Intact"
data_DFS$Cases[which( (data_DFS$PTEN_cna == "-2") & (data_DFS$quantile_joint == "high") )] <- "Loss/High"
data_DFS$Cases[which( (data_DFS$PTEN_cna == "-2") & (data_DFS$quantile_joint == "low") )] <- "Loss/Low"

################################################################################




################################################################################
#| BASURTO: MODULE GREEN Q1 and Q4
################################################################################

tag <- "Module Green Q1_Q4"
data_DFS <- sample_info[which(!is.na(sample_info$DFS.STATUS)),]

#| Adding module green
genes_mod_green <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)
counts<- counts_zscore[,which(colnames(counts_zscore) %in% data_DFS$AC.basurto)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
data_DFS$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| Getting the quartiles
data_DFS$quantile <- cut(as.numeric(data_DFS$mod_green_zscore),breaks=quantile(as.numeric(data_DFS$mod_green_zscore),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_DFS <- data_DFS[which( (data_DFS$quantile == 1) | (data_DFS$quantile == 4)),]
data_DFS$quantile[which(data_DFS$quantile ==1)] <- "High risk"
data_DFS$quantile[which(data_DFS$quantile ==4)] <- "Low risk"

#| Applying cox
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ quantile,  data=data_DFS)
fit_summary <- summary(fit)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~  quantile , data=data_DFS)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 4), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)

# Add a starting row (time = 0) for each strata if missing
strata_starts <- sfit %>%
  group_by(strata) %>%
  summarise(min_time = min(time), .groups = "drop") %>%
  filter(min_time > 0)

starting_rows <- strata_starts %>%
  rowwise() %>%
  mutate(
    n.risk = max(sfit$n.risk[sfit$strata == strata]),
    n.event = 0,
    n.censor = 0,
    surv = 1,
    std.err = 0,
    upper = 1,
    lower = 1,
    time = 0
  ) %>%
  select(time, n.risk, n.event, n.censor, surv, std.err, upper, lower, strata)

# Combine the original and new rows, and sort
sfit <- bind_rows(sfit, starting_rows) %>%
  arrange(strata, time)

ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point(shape=16, size =2) +
  theme(plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Value") +
  scale_color_manual(values =  c( "blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for ",tag, "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave(paste0("Results/Basurto_Kaplan-Meier_plots/",tag,".pdf"), heigh= 4.5, width = 6)

################################################################################



################################################################################
#| BASURTO: MODULE GREEN PTEN LOSS Q1 and Q4
################################################################################

tag <- "PTEN loss-presence Mod Green High-Low_Q1_Q4"
data_DFS <- sample_info[which(!is.na(sample_info$DFS.STATUS)),]
data_DFS <- data_DFS[which(!is.na(data_DFS$H_score)),]

#| Adding module green
genes_mod_green <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)
counts<- counts_zscore[,which(colnames(counts_zscore) %in% data_DFS$AC.basurto)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
data_DFS$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| Getting the quartiles
data_DFS$quantile <- cut(as.numeric(data_DFS$mod_green_zscore),breaks=quantile(as.numeric(data_DFS$mod_green_zscore),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_DFS$quantile_joint<-NA
data_DFS$quantile_joint[which( (data_DFS$quantile == 1) | (data_DFS$quantile == 2))]<-"low"
data_DFS$quantile_joint[which( (data_DFS$quantile == 4) | (data_DFS$quantile == 3))]<-"high"

#| Separating cases
data_DFS$Cases <- NA
data_DFS$Cases[which(data_DFS$PTEN_status == "presence")] <- "Presence"
data_DFS$Cases[which( (data_DFS$PTEN_status == "loss") & (data_DFS$quantile_joint == "high") )] <- "Loss/High"
data_DFS$Cases[which( (data_DFS$PTEN_status == "loss") & (data_DFS$quantile_joint == "low") )] <- "Loss/Low"

data_DFS <- data_DFS[which((data_DFS$quantile == 1) | (data_DFS$quantile == 4)),]

#| Applying cox
fit <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ Cases,  data=data_DFS)
fit_summary <- summary(fit)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~  Cases , data=data_DFS)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 3), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)

# Add a starting row (time = 0) for each strata if missing
strata_starts <- sfit %>%
  group_by(strata) %>%
  summarise(min_time = min(time), .groups = "drop") %>%
  filter(min_time > 0)

starting_rows <- strata_starts %>%
  rowwise() %>%
  mutate(
    n.risk = max(sfit$n.risk[sfit$strata == strata]),
    n.event = 0,
    n.censor = 0,
    surv = 1,
    std.err = 0,
    upper = 1,
    lower = 1,
    time = 0
  ) %>%
  select(time, n.risk, n.event, n.censor, surv, std.err, upper, lower, strata)

# Combine the original and new rows, and sort
sfit <- bind_rows(sfit, starting_rows) %>%
  arrange(strata, time)

ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point(shape=16, size =2) +
  theme(plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Cases") +
  scale_color_manual(values = c("#FFEA46FF", "blue", "black"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for \n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave(paste0("Results/Basurto_Kaplan-Meier_plots/",tag,".pdf"), heigh= 4.5, width = 6)

library(ggsurvfit)
ggsurvfit( survfit(Surv(DFS.TIME, DFS.STATUS) ~  Cases , data=data_DFS)) +
  theme(plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Cases") +
  scale_color_manual(values = viridis(3, option="C"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for \n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))+
  add_confidence_interval() +
  add_risktable()
ggsave(paste0("Results/Basurto_Kaplan-Meier_plots/",tag,"_risk_table.pdf"), heigh= 4.5, width = 6)
################################################################################


################################################################################
#| TCGA: MODULE GREEN Q1 and Q4
################################################################################

tag <- "Module Green Q1_Q4"
data_DFS <- TCGA_data[which(!is.na(TCGA_data$DFS_STATUS)),]

#| Adding module green
genes_mod_green <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)
counts<- counts_zscore_tcga[,which(colnames(counts_zscore_tcga) %in% data_DFS$SampleID)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
data_DFS$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| Getting the quartiles
data_DFS$quantile <- cut(as.numeric(data_DFS$mod_green_zscore),breaks=quantile(as.numeric(data_DFS$mod_green_zscore),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_DFS <- data_DFS[which( (data_DFS$quantile == 1) | (data_DFS$quantile == 4)),]

data_DFS$quantile[which(data_DFS$quantile ==1)] <- "High risk"
data_DFS$quantile[which(data_DFS$quantile ==4)] <- "Low risk"

fit <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ quantile ,  data=data_DFS)
fit_summary <- summary(fit)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~  quantile  , data=data_DFS)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 4), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)

# Add a starting row (time = 0) for each strata if missing
strata_starts <- sfit %>%
  group_by(strata) %>%
  summarise(min_time = min(time), .groups = "drop") %>%
  filter(min_time > 0)

starting_rows <- strata_starts %>%
  rowwise() %>%
  mutate(
    n.risk = max(sfit$n.risk[sfit$strata == strata]),
    n.event = 0,
    n.censor = 0,
    surv = 1,
    std.err = 0,
    upper = 1,
    lower = 1,
    time = 0
  ) %>%
  select(time, n.risk, n.event, n.censor, surv, std.err, upper, lower, strata)

# Combine the original and new rows, and sort
sfit <- bind_rows(sfit, starting_rows) %>%
  arrange(strata, time)

ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point(shape=16, size =2) +
  theme(plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Value") +
  scale_color_manual(values = c( "blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for ",tag, "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave(paste0("Results/TCGA_Kaplan-Meier_plots/",tag,".pdf"), heigh= 4.5, width = 6)
################################################################################



################################################################################
#| TCGA: MODULE GREEN PTEN LOSS Q1 and Q4
################################################################################

#| Genomic loss PTEN intact and loss with module green Low and High
tag <- "PTEN genomic loss-intact Mod Green High-Low_Q1_Q4"
data_DFS <- TCGA_data[which(!is.na(TCGA_data$DFS_STATUS)),]
data_DFS <- data_DFS[ which( (data_DFS$PTEN_cna == "0") | (data_DFS$PTEN_cna == "-2")  ), ]

#| Adding module green
genes_mod_green <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)
counts<- counts_zscore_tcga[,which(colnames(counts_zscore_tcga) %in% data_DFS$SampleID)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
data_DFS$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

#| Getting the quartiles
data_DFS$quantile <- cut(as.numeric(data_DFS$mod_green_zscore),breaks=quantile(as.numeric(data_DFS$mod_green_zscore),na.rm = TRUE), include.lowest=TRUE, labels=FALSE)
data_DFS$quantile_joint<-NA
data_DFS$quantile_joint[which( (data_DFS$quantile == "1") | (data_DFS$quantile == "2"))]<-"low"
data_DFS$quantile_joint[which( (data_DFS$quantile == "4") | (data_DFS$quantile == "3"))]<-"high"

#| Separating cases
data_DFS$Cases <- NA
data_DFS$Cases[which(data_DFS$PTEN_cna == "0")] <- "Intact"
data_DFS$Cases[which( (data_DFS$PTEN_cna == "-2") & (data_DFS$quantile_joint == "high") )] <- "Loss/High"
data_DFS$Cases[which( (data_DFS$PTEN_cna == "-2") & (data_DFS$quantile_joint == "low") )] <- "Loss/Low"

#| fitting cox
fit <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ Cases ,  data=data_DFS)
fit_summary <- summary(fit)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~  Cases  , data=data_DFS)
p_value <- surv_pvalue(sfit)
p_value <- format(round(p_value[[2]], 4), nsmall = 2)
sfit <- fortify(sfit)
sfit$strata <- as.character(sfit$strata)

# Add a starting row (time = 0) for each strata if missing
strata_starts <- sfit %>%
  group_by(strata) %>%
  summarise(min_time = min(time), .groups = "drop") %>%
  filter(min_time > 0)

starting_rows <- strata_starts %>%
  rowwise() %>%
  mutate(
    n.risk = max(sfit$n.risk[sfit$strata == strata]),
    n.event = 0,
    n.censor = 0,
    surv = 1,
    std.err = 0,
    upper = 1,
    lower = 1,
    time = 0
  ) %>%
  select(time, n.risk, n.event, n.censor, surv, std.err, upper, lower, strata)

# Combine the original and new rows, and sort
sfit <- bind_rows(sfit, starting_rows) %>%
  arrange(strata, time)

ggplot(sfit, aes(x=time, y= surv, color =strata)) +
  geom_line() + 
  geom_point(shape=16, size =2) +
  theme(text=element_text(size=16,  family="sans"), 
        plot.title=element_text(size=14, hjust = 0))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Value") +
  scale_color_manual(values = c("black", "#FFEA46FF","blue"  ))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for  \n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave(paste0("Results/TCGA_Kaplan-Meier_plots/",tag,".pdf"), heigh= 4.5, width = 6)


ggsurvfit( survfit(Surv(DFS_MONTHS, DFS_STATUS) ~  Cases , data=data_DFS)) +
  theme(plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"))+ 
  ylim(c(0,1)) + 
  ylab("Recurrence-free survival (RFS)") + 
  xlab("Time to recurrence (months)") +
  labs(color="Cases") +
  scale_color_manual(values = viridis(3, option="C"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_value,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for \n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))+
  add_confidence_interval() +
  add_risktable()
ggsave(paste0("Results/TCGA_Kaplan-Meier_plots/",tag,"_risk_table.pdf"), heigh= 4.5, width = 6)
################################################################################
