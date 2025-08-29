################################################################################
#| MUTIVARIATE COX IN BASURTO AND TCGA COHORTS
################################################################################
#| Date: 29/08/2025
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script performs multivariate modeling to evaluate the prognostic impact 
#| of WGCNA module signatures (particularly the green module) and clinical 
#| covariates in prostate cancer. The analysis is trained on the Basurto cohort 
#| (with H-score-based PTEN protein assessment) and validated on TCGA data 
#| (using PTEN genomic alterations).  
#|
#| Workflow:
#|   1) Load and preprocess Basurto RNA-seq counts:
#|        * Filter, normalize (DESeq2), log2 transform, and z-score.  
#|        * Merge with WGCNA modules and DEGs to extract module-level expression.  
#|        * Annotate samples with PTEN protein status, Gleason score, and PSA.  
#|   2) Load and preprocess TCGA data:
#|        * Clinical, CNA, protein, and normalized RNA-seq data.  
#|        * Compute z-scores and module expression.  
#|        * Define PTEN status by CNA and protein quantiles.   
#|   3) Fit multivariate Cox regression models:
#|        * Model 1: PSA + Gleason + Green module.  
#|        * Model 2: PSA + Gleason.   
#|   4) Derive risk scores, stratify patients into high/low risk groups.  
#|   5) Plot Kaplan-Meier survival curves and calculate hazard ratios (Basurto and TCGA).  
#|   6) Generate forest plots summarizing hazard ratios with confidence intervals.  
#|
#| Outputs:
#|   - Frequency distributions of module scores (Basurto and TCGA).  
#|   - Kaplan-Meier survival plots for training and validation sets.  
#|   - Forest plots of multivariate Cox regression models.  
#|   - Tables with scaled/annotated sample info (PTEN, modules, signatures).  
################################################################################


################################################################################
#| LIBRARY AND DATA
################################################################################
suppressMessages(library(WGCNA))
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
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(broom))
suppressMessages(library(survival))
suppressMessages(library(timeROC))

#| Setting working directory
setwd("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Multivariant_model/")

#| For plots
theme_set(theme_classic())

#| Data pathways 
data.file_WGCNA <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_intramodularinfo_last_try.txt"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Tables/Limma_voom_analysis_DEGS_results.txt"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Images/11_MachineLearning_Regression_Models_to_identify_candidates/"
sample_info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
sample_info_PTEN_pAKT.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted_pAKT.txt"
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
sample_info$PTEN_protein_status[which(sample_info$PTEN_protein_status == "loss")] <- 1
sample_info$PTEN_protein_status[which(sample_info$PTEN_protein_status == "presence")] <- 0
sample_info$PTEN_protein_status <- as.numeric(sample_info$PTEN_protein_status)
sample_info$PTEN_Exp_log2 <- as.numeric(scale(sample_info$PTEN_Exp_log2))
sample_info$PTEN_protein_status <- as.numeric(scale(sample_info$H_score))
sample_info$mod_purple <- as.numeric(scale(samples_media_expression_purple))
sample_info$mod_green <- as.numeric(scale(samples_media_expression_green))

#| Module green distribution
ggplot(sample_info, aes(x=mod_green))+
  geom_histogram(color = "black", binwidth = 1, fill = "orange3")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12, color = "black"),
        axis.title.x  = element_text(size = 14, color = "black"),
        axis.title.y  = element_text(size = 14, color = "black"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, color = "black"))+
  ylab("Frequency")
ggsave("Basurto_Green_Module_Frequency_distribution_197.pdf", width = 5, height = 4.5)


###| TCGA DATA
#| Loading TCGA data already processed
TCGA_data <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/Sample_info/sample_info_TCGA_PTEN_CNA_complete.txt", sep ="\t")
TCGA_xcell <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/xCell/xCell_results_TPM_all.txt", sep ="\t")

#| Sample ID
TCGA_data$SampleID <- TCGA_data$PATIENT_ID
TCGA_xcell$SampleID <- rownames(TCGA_xcell)

#| Merging these dataframes
TCGA_data <- merge(TCGA_xcell,TCGA_data, by ="SampleID" )
names(TCGA_data)
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

TCGA_data$mod_purple <- as.numeric(scale(samples_media_expression_purple))
TCGA_data$mod_green <- as.numeric(scale(samples_media_expression_green))
#| Module green distribution
ggplot(TCGA_data, aes(x=mod_green))+
  geom_histogram(color = "black", binwidth = 1, fill = "orange3")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12, color = "black"),
        axis.title.x  = element_text(size = 14, color = "black"),
        axis.title.y  = element_text(size = 14, color = "black"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, color = "black"))+
  ylab("Frequency")
ggsave("TCGA_Green_Module_Frequency_distribution_ALL.pdf", width = 5, height = 4.5)

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
#| COX MULTIVARIATE REGRESSION ANALYSIS
################################################################################


# ----------- Step 1: Fit Cox model -----------

train_data <- sample_info[which(!is.na(sample_info$DFS.STATUS)),]

genes_mod_green <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)
counts<- counts_zscore[,which(colnames(counts_zscore) %in% train_data$AC.basurto)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
train_data$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))
train_data$PSA <- as.numeric(scale(log2(train_data$PSA +1) ))

cox_model_1 <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ PSA + 
                     Gleason_score_pieza + 
                     mod_green_zscore, data = train_data)

cox_model_2 <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ PSA + 
                       Gleason_score_pieza, data = train_data)

# ----------- Step 2: Predict Risk Score -----------
train_data_1 <- train_data[which(rownames(train_data) %in% names(cox_model_1$y)),]
train_data_1$risk_score <- predict(cox_model_1, type = "lp")

train_data_2 <- train_data[which(rownames(train_data) %in% names(cox_model_2$y)),]
train_data_2$risk_score <- predict(cox_model_2, type = "lp")

# ----------- Step 3: Stratify Risk Groups -----------
train_data_1$risk_group <- ifelse(train_data_1$risk_score > median(train_data_1$risk_score), "High Risk", "Low Risk")
train_data_2$risk_group <- ifelse(train_data_2$risk_score > median(train_data_2$risk_score), "High Risk", "Low Risk")

# ----------- Step 4: Kaplan-Meier Survival Plot -----------

#| Model 1
cox <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ risk_group, data = train_data_1)
fit_summary <- summary(cox)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ risk_group, data = train_data_1)
p_val <- surv_pvalue(sfit)$pval
p_val_formatted <- format(round(p_val, 4), nsmall = 2)
sfit <- surv_summary(sfit)
sfit$strata <- as.character(sfit$risk_group)
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
  scale_color_manual(values = c("blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_val_formatted,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for Multivariate Cox", "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave("Results/Multivariate Cox Regression/training_BASURTO_Kaplan-Meier_Model_with_Green-Mod.pdf",heigh= 4.5, width = 6)

#| Model 2
cox <- coxph(Surv(DFS.TIME, DFS.STATUS) ~ risk_group, data = train_data_2)
fit_summary <- summary(cox)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS.TIME, DFS.STATUS) ~ risk_group, data = train_data_2)
p_val <- surv_pvalue(sfit)$pval
p_val_formatted <- format(round(p_val, 4), nsmall = 2)
sfit <- surv_summary(sfit)
sfit$strata <- as.character(sfit$risk_group)
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
  scale_color_manual(values = c("blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_val_formatted,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for Multivariate Cox", "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave("Results/Multivariate Cox Regression/training_BASURTO_Kaplan-Meier_Model_without_Green-Mod.pdf",heigh= 4.5, width = 6)

# ----------- Step 5: Forest Plot -----------

# Tidy the model
cox_df <- tidy(cox_model_1, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(variable = term,
         HR = estimate,
         lower = conf.low,
         upper = conf.high)
cox_df$variable[which(cox_df$variable == "Gleason_score_pieza")] <- "Gleason Score"
cox_df$variable[which(cox_df$variable == "mod_green")] <- "Module Green"

# Set factor levels to control order
cox_df$variable <- factor(cox_df$variable, levels = rev(cox_df$variable))

# Plot
ggplot(cox_df, aes(x = HR, y = variable, color = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, size = 0.8) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 3)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size =12, color ="black"),
    axis.title.y = element_blank(),
    plot.title=element_text(size=14, hjust = 0),
    text=element_text(size=14,  family="sans"),
    axis.text.x = element_text(size =12, color ="black")) +
  labs(
    x = "Hazard Ratio (log scale)")+
  scale_color_manual(values = rainbow(3))
ggsave("Results/Multivariate Cox Regression/training_BASURTO_Forest-plot.pdf",heigh= 4.5, width = 5)


# ----------- Step 6: Validate on Test Data  -----------
test_data <- TCGA_data[which(!is.na(TCGA_data$DFS_STATUS)),]

genes_mod_green <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)
counts<- counts_zscore_tcga[,which(colnames(counts_zscore_tcga) %in% test_data$SampleID)]
samples_media_expression_green <- colMeans(counts[genes_mod_green,], na.rm=TRUE)
test_data$mod_green_zscore <- as.numeric(scale(samples_media_expression_green))

test_data$PSA <- as.numeric(test_data$PSA_MOST_RECENT_RESULTS)
test_data$PSA <- scale(log2(test_data$PSA +1))
test_data$PSA <- as.numeric(test_data$PSA )

cox_model_1 <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ PSA + 
                     Gleason_score_pieza + 
                     mod_green_zscore , data = test_data)

cox_model_2 <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ PSA + 
                       Gleason_score_pieza, data = test_data)

test_data_1 <- test_data[which(rownames(test_data) %in% names(cox_model_1$y)),]
test_data_2 <- test_data[which(rownames(test_data) %in% names(cox_model_2$y)),]

test_data_1$risk_score <- predict(cox_model_1, newdata = test_data_1, type = "lp")
test_data_2$risk_score <- predict(cox_model_2, newdata = test_data_2, type = "lp")

test_data_1$risk_group <- ifelse(test_data_1$risk_score > median(test_data_1$risk_score), "High Risk", "Low Risk")
test_data_2$risk_group <- ifelse(test_data_2$risk_score > median(test_data_2$risk_score), "High Risk", "Low Risk")

# Model 1
cox <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ risk_group, data = test_data_1)
fit_summary <- summary(cox)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ risk_group, data = test_data_1)
p_val <- surv_pvalue(sfit)$pval
p_val_formatted <- format(round(p_val, 12), nsmall = 2)
sfit <- surv_summary(sfit)
sfit$strata <- as.character(sfit$risk_group)
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
  scale_color_manual(values = c("blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_val_formatted,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for Multivariate Cox", "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave("Results/Multivariate Cox Regression/test_TCGA_Kaplan-Meier_Model_with_Green-Mod.pdf",heigh= 4.5, width = 6)


# Model 2
cox <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ risk_group, data = test_data_2)
fit_summary <- summary(cox)
p_value_cox <- fit_summary$coefficients[1, 5]
coef_gene_cox <- fit_summary$coefficients[1,1]
exp_coef_gene_cox <- fit_summary$coefficients[1,2]

sfit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ risk_group, data = test_data_2)
p_val <- surv_pvalue(sfit)$pval
p_val_formatted <- format(round(p_val, 14), nsmall = 2)
sfit <- surv_summary(sfit)
sfit$strata <- as.character(sfit$risk_group)
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
  scale_color_manual(values = c("blue","#FFEA46FF"))+
  annotate(geom="text", x=20, y=0.05, label=paste("p < ",p_val_formatted,sep=""), color="black", family ="sans", size =4.5) +
  ggtitle(paste("Hazard Ratio for Multivariate Cox", "\n ",round(exp_coef_gene_cox,2), "; p < ", format(p_value_cox,scientific = T, digits =2), sep =""))
ggsave("Results/Multivariate Cox Regression/test_TCGA_Kaplan-Meier_Model_without_Green-Mod.pdf",heigh= 4.5, width = 6)


# Tidy the model
cox_df <- tidy(cox_model_1, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(variable = term,
         HR = estimate,
         lower = conf.low,
         upper = conf.high)
cox_df$variable[which(cox_df$variable == "Gleason_score_pieza")] <- "Gleason Score"
cox_df$variable[which(cox_df$variable == "mod_green")] <- "Module Green"

# Set factor levels to control order
cox_df$variable <- factor(cox_df$variable, levels = rev(cox_df$variable))

# Plot
ggplot(cox_df, aes(x = HR, y = variable, color = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, size = 0.8) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 3)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size =12, color ="black"),
        axis.title.y = element_blank(),
        plot.title=element_text(size=14, hjust = 0),
        text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black")) +
  labs(
    x = "Hazard Ratio (log scale)")+
  scale_color_manual(values = rainbow(3))
ggsave("Results/Multivariate Cox Regression/test_TCGA_Forest-plot.pdf",heigh= 4.5, width = 5)
################################################################################


