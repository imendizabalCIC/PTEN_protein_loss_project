################################################################################
#|           MuSiC for BulkRNAseq and single-cell data
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#| 
#| Description:
#| This script applies the MuSiC framework to deconvolve bulk RNA-seq datasets 
#| (Basurto human cohort and AC-12 mouse model) using Chen et al.'s prostate tumor 
#| single-cell RNA-seq dataset as reference. The analysis estimates cell type 
#| proportions in bulk samples and compares them across PTEN statuses or 
#| experimental conditions.  
#|
#| Workflow:
#|   1) Load and preprocess bulk RNA-seq data:  
#|        * Basurto human cohort.  
#|        * AC-12 mouse RNA-seq (map mouse -> human orthologs).  
#|   2) Load annotated scRNA-seq reference (Chen localized tumors) and prepare as 
#|      a SingleCellExperiment object (RNA assay).  
#|   3) Run MuSiC to estimate cell type proportions in bulk samples using reference 
#|      cluster annotations (luminal, T cells, macrophages, fibroblasts, endothelial, 
#|      mast cells, cycling).  
#|   4) Stratify bulk samples by PTEN status:  
#|        * Human Basurto cohort -> PTEN loss vs PTEN presence.  
#|        * Mouse AC-12 model -> KO vs WT at 3 months and 6 months.  
#|   5) Statistical comparison of estimated cell type proportions between groups 
#|      (boxplots with jittered points, Wilcoxon test p-values).  
#|
#| Outputs:
#|   - Normalized and filtered bulk RNA-seq matrices (human + mouse).  
#|   - Estimated cell type proportions for each bulk sample.  
#|   - Boxplots of cell type proportions by PTEN status (Basurto cohort).  
#|   - Boxplots of cell type proportions by genotype and timepoint (AC-12 mouse).  
#|   - PDF/PNG figures for publication-ready visualization.  
################################################################################


################################################################################
#| LIBRARIES 
################################################################################
suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("Matrix"))
suppressPackageStartupMessages(require("RCurl"))
suppressPackageStartupMessages(require("scales"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("devtools"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("dplyr"))
suppressPackageStartupMessages(require("ensembldb"))
suppressPackageStartupMessages(require("reprex"))
suppressPackageStartupMessages(require("remotes"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("DESeq2"))
suppressPackageStartupMessages(require("viridis"))
suppressPackageStartupMessages(require("multtest"))
suppressPackageStartupMessages(require("metap"))
suppressPackageStartupMessages(require("gprofiler2"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("UCell"))
suppressPackageStartupMessages(library("MuSiC"))
suppressPackageStartupMessages(library("MuSiC2"))
suppressPackageStartupMessages(library("ggpubr"))

#| Setting theme to classic
theme_set(theme_classic())

#| Set working directory
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/SingleCell/"
setwd(workingDir)

################################################################################


################################################################################
#| DATA BASURTO
################################################################################

#| Bulk Basurto
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/Xcell_output.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt"

#| Counts data
counts_data <- read.table(counts.file, header=T) 

#| Sample info
datTraits <- read.table(info.file, sep ="\t", header =T)

#| Filtering BPH and H-score NA
datTraits <- datTraits[which(!is.na(datTraits$H_score_cut_0)),]
counts_data <- counts_data[,rownames(datTraits)]

#| Filtering low counts
filter_counts <- rowSums(counts_data>5) >= 0.7*ncol(counts_data)
counts_data <- counts_data[filter_counts,]

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

counts_data$GeneID <- rownames(counts_data)
counts_data <- merge(counts_data, genome_GRCh39.94, by ="GeneID")
counts_data <- aggregate(counts_data, list(by = counts_data$gene_name), mean)
rownames(counts_data) <- counts_data$by
counts_data <- counts_data[,datTraits$AC.basurto]
counts_data <- round(counts_data)

#| PTEN loss
counts_data_PTEN_loss <- counts_data[, datTraits$AC.basurto[which(datTraits$PTEN_status == "PTEN loss")]]
counts_data_PTEN_loss <- as.matrix(counts_data_PTEN_loss)

#| PTEN intact
counts_data_PTEN_intact <- counts_data[, datTraits$AC.basurto[which(datTraits$PTEN_status == "PTEN intact")]]
counts_data_PTEN_intact <- as.matrix(counts_data_PTEN_intact)

counts_data <- as.matrix(counts_data)
################################################################################


################################################################################
#| MuSiC WITH DATA CHEN LOCALIZED
################################################################################

#| Data filtered by localized samples
loc.combined.sct <- readRDS("Annotation_Manual_Chen_PT.rds")
DefaultAssay(loc.combined.sct) <- "RNA"
loc.combined.sct$cluster_names <- as.character(Idents(loc.combined.sct))
unique(loc.combined.sct$cluster_names)

#| Transforming to singlecell experiment
sce_localized <- as.SingleCellExperiment(loc.combined.sct, assay = "RNA")

# Estimate cell type proportions
Est.prop_localized <- music_prop(bulk.mtx = counts_data, 
                       sc.sce  = sce_localized, 
                       clusters = 'cluster_names', 
                       samples = 'orig.ident', 
                       select.ct = c("4-Luminal", "6-T cell", "1-Luminal", "13-Luminal", "3-Luminal",
                                     "9-Macrophage", "2-Luminal", "10-Myofibroblasts" , "5-Luminal" , 
                                     "8-Intermediate", "15-Luminal", "0-Luminal", "11-Luminal", "7-Endothelial", 
                                     "14-Luminal", "17-Luminal cycling", "16-Endothelial", "12-Mast cell"), 
                       verbose = F)

#| MuSiC estimated proportions
prop.weighted_localized <- as.data.frame(Est.prop_localized$Est.prop.weighted)
prop.weighted_localized$PTEN_status <- NA
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_PTEN_intact))] <- "presence"  #| Changing to presence instead of intact
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_PTEN_loss))] <- "loss"
prop.weighted_localized_data <- stack(prop.weighted_localized)
prop.weighted_localized_data$PTEN_status <- prop.weighted_localized$PTEN_status
prop.weighted_localized_data$values <- as.numeric(prop.weighted_localized_data$values)
prop.weighted_localized_data$PTEN_status_factor <- factor(prop.weighted_localized_data$PTEN_status, level = c("presence", "loss"))

#| Plotting
ggplot(prop.weighted_localized_data[which(prop.weighted_localized_data$ind != "PTEN_status"),], aes( x = reorder(ind, -values), y = values, fill = PTEN_status_factor)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status_factor),
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),
             pch = 21, size=1.5, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =2, vjust = -0.6)+
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = 0.5, face ="bold"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  xlab("")+
  ylab("Cell Type Proportions") +
  scale_fill_manual(values = c("#EFD5B2","#504667")) +
  labs(fill ="PTEN protein status")
#ggsave("Results/Images/MuSiC/Cell_Type_Proportions_all_PTEN_protein_loss_vs_intact_Basurto_Chen_Localized.pdf", height =5, width = 12.5)
ggsave("Results/Images/MuSiC/Cell_Type_Proportions_all_PTEN_protein_loss_vs_intact_Basurto_Chen_Localized_last.pdf", height =5, width = 12.5)

################################################################################


################################################################################
#| DATA MOUSE
################################################################################

#| Sample info
sample <- colnames(counts_data_human_mouse)
condition <- gsub("_[0-9]*","",sample)
month <- gsub("[A-Z]", "", condition) 
sample_info <- data.frame(sample = sample,
                          condition = condition,
                          month = month)
rownames(sample_info) <- sample_info$sample
sample_info_mouse <- sample_info[order(sample_info$sample),]
sample_info_mouse$condition <- factor(sample_info_mouse$condition, levels =c("WT3", "WT6", "KO3", "KO6"))

#| Counts data
counts_data_mouse <- read.table( "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt", sep ="\t")
counts_data_mouse <- counts_data_mouse[,order(colnames(counts_data_mouse))]

#| Filtering low counts
filter_counts <- rowSums(counts_data_mouse>0) >= 0.9*ncol(counts_data_mouse)
counts_data_mouse <- counts_data_mouse[filter_counts,]

#| From mouse gene name to human data
mouse_human <- read.table("X:/DATA_shared/Human_to_Mouse_Gene_name.txt", sep ="\t", header =T )
samples <- colnames(counts_data_mouse)
names(mouse_human)

counts_data_human_mouse <- counts_data_mouse
counts_data_human_mouse$Mouse.gene.stable.ID <- rownames(counts_data_human_mouse)
counts_data_human_mouse <- merge(counts_data_human_mouse, mouse_human, by ="Mouse.gene.stable.ID")
counts_data_human_mouse <- aggregate(counts_data_human_mouse, by =list(counts_data_human_mouse$Gene.name), mean)
rownames(counts_data_human_mouse) <- counts_data_human_mouse$Group.1
counts_data_human_mouse <- counts_data_human_mouse[, samples]
counts_data_human_mouse <- as.data.frame(counts_data_human_mouse)
counts_data_human_mouse<-round(counts_data_human_mouse)

####################| KO3 VS WT3
counts_data_KO3 <- counts_data_human_mouse[, sample_info$sample[which(sample_info$condition == "KO3")]]
counts_data_WT3 <- counts_data_human_mouse[, sample_info$sample[which(sample_info$condition == "WT3")]]

counts_data_KO3 <- as.matrix(counts_data_KO3)
counts_data_WT3 <- as.matrix(counts_data_WT3)

counts_data_KO3_WT3 <- as.matrix(counts_data_human_mouse[,sample_info$sample[which(sample_info$month == "3")]])


####################| KO6 VS WT6
counts_data_KO6 <- counts_data_human_mouse[, sample_info$sample[which(sample_info$condition == "KO6")]]
counts_data_WT6 <- counts_data_human_mouse[, sample_info$sample[which(sample_info$condition == "WT6")]]

counts_data_KO6 <- as.matrix(counts_data_KO6)
counts_data_WT6 <- as.matrix(counts_data_WT6)

counts_data_KO6_WT6 <- as.matrix(counts_data_human_mouse[,sample_info$sample[which(sample_info$month == "6")]])

################################################################################


################################################################################
#| Music Chen Localized MOUSE 6 MONTHS
################################################################################

#| Data filtered by localized samples
loc.combined.sct <- readRDS("Annotation_Manual_Chen_PT.rds")
DefaultAssay(loc.combined.sct) <- "RNA"
loc.combined.sct$cluster_names <- as.character(Idents(loc.combined.sct))

#| Transforming to singlecell experiment
sce_localized <- as.SingleCellExperiment(loc.combined.sct, assay = "RNA")
intersect(rownames(sce_localized), rownames(counts_data_KO6_WT6))

#| Extracting common genes
common_genes <- intersect(rownames(sce_localized), rownames(counts_data_KO6_WT6))
cat("Number of common genes:", length(common_genes), "\n")

#| Subset matrices to common genes
counts_data_KO6_WT6 <- counts_data_KO6_WT6[common_genes, , drop = FALSE]
sce_localized <- sce_localized[common_genes, , drop = FALSE]

#| Estimate cell type proportions
Est.prop_localized <- music_prop(bulk.mtx = counts_data_KO6_WT6, 
                                 sc.sce  = sce_localized, 
                                 clusters = 'cluster_names', 
                                 samples = 'orig.ident', 
                                 select.ct = c("4-Luminal", "6-T cell", "1-Luminal", "13-Luminal", "3-Luminal",
                                               "9-Macrophage", "2-Luminal", "10-Myofibroblasts" , "5-Luminal" , 
                                               "8-Intermediate", "15-Luminal", "0-Luminal", "11-Luminal", "7-Endothelial", 
                                               "14-Luminal", "17-Luminal cycling", "16-Endothelial", "12-Mast cell"), 
                                 verbose = F)

#| MuSiC estimated proportions
prop.weighted_localized <- as.data.frame(Est.prop_localized$Est.prop.weighted)
prop.weighted_localized$PTEN_status <- NA
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_WT6))] <- "WT6"
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_KO6))] <- "KO6"
prop.weighted_localized_data <- stack(prop.weighted_localized)
prop.weighted_localized_data$PTEN_status <- prop.weighted_localized$PTEN_status
prop.weighted_localized_data$values <- as.numeric(prop.weighted_localized_data$values)
prop.weighted_localized_data$PTEN_status_factor <- factor(prop.weighted_localized_data$PTEN_status, level = c("WT6", "KO6"))

#| Plotting
ggplot(prop.weighted_localized_data[which(prop.weighted_localized_data$ind != "PTEN_status"),], aes( x = reorder(ind, -values), y = values, fill = PTEN_status_factor)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status_factor),
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),
             pch = 21, size=1.5, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =2, vjust = -0.6)+
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = 0.5, face ="bold"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  xlab("")+
  ylab("Cell Type Proportions") +
  scale_fill_manual(values = c(c("#ff9700","#6C2EA5"))) +
  labs(fill ="PTEN protein status")
#ggsave("Results/Images/MuSiC/Cell_Type_Proportions_all_PTEN_protein_loss_vs_intact_Basurto_Chen_Localized.pdf", height =5, width = 12.5)
ggsave("Results/Images/MuSiC/Cell_Type_Proportions_all_Pten_KO6_vs_WT6_Mouse_Chen_Localized_last.pdf", height =5, width = 12.5)

################################################################################


################################################################################
#| Music Chen Localized MOUSE 3 MONTHS
################################################################################

#| Data filtered by localized samples
loc.combined.sct <- readRDS("Annotation_Manual_Chen_PT.rds")
DefaultAssay(loc.combined.sct) <- "RNA"
loc.combined.sct$cluster_names <- as.character(Idents(loc.combined.sct))

#| Transforming to singlecell experiment
sce_localized <- as.SingleCellExperiment(loc.combined.sct, assay = "RNA")

#| Extracting common genes
common_genes <- intersect(rownames(sce_localized), rownames(counts_data_KO3_WT3))
cat("Number of common genes:", length(common_genes), "\n")

#| Subset matrices to common genes
counts_data_KO3_WT3 <- counts_data_KO3_WT3[common_genes, , drop = FALSE]
sce_localized <- sce_localized[common_genes, ,drop = FALSE]

#| Estimate cell type proportions
Est.prop_localized <- music_prop(bulk.mtx = counts_data_KO3_WT3, 
                                 sc.sce  = sce_localized, 
                                 clusters = 'cluster_names', 
                                 samples = 'orig.ident', 
                                 select.ct = c("4-Luminal", "6-T cell", "1-Luminal", "13-Luminal", "3-Luminal",
                                               "9-Macrophage", "2-Luminal", "10-Myofibroblasts" , "5-Luminal" , 
                                               "8-Intermediate", "15-Luminal", "0-Luminal", "11-Luminal", "7-Endothelial", 
                                               "14-Luminal", "17-Luminal cycling", "16-Endothelial", "12-Mast cell"), 
                                 verbose = F)

#| MuSiC estimated proportions
prop.weighted_localized <- as.data.frame(Est.prop_localized$Est.prop.weighted)
prop.weighted_localized$PTEN_status <- NA
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_WT3))] <- "WT3"
prop.weighted_localized$PTEN_status[which(rownames(prop.weighted_localized) %in% colnames(counts_data_KO3))] <- "KO3"
prop.weighted_localized_data <- stack(prop.weighted_localized)
prop.weighted_localized_data$PTEN_status <- prop.weighted_localized$PTEN_status
prop.weighted_localized_data$values <- as.numeric(prop.weighted_localized_data$values)


prop.weighted_localized_data$PTEN_status_factor <- factor(prop.weighted_localized_data$PTEN_status, level = c("WT3", "KO3"))
#| Plotting
ggplot(prop.weighted_localized_data[which(prop.weighted_localized_data$ind != "PTEN_status"),], aes( x = reorder(ind, -values), y = values, fill = PTEN_status_factor)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status_factor),
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),
             pch = 21, size=1.5, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =2, vjust = -0.6)+
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = 0.5, face ="bold"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  xlab("")+
  ylab("Cell Type Proportions") +
  scale_fill_manual(values = c(c("#ff9700","#6C2EA5"))) +
  labs(fill ="PTEN protein status")
ggsave("Results/Images/MuSiC/Cell_Type_Proportions_all_Pten_KO3_vs_WT3_Mouse_Chen_Localized_last.pdf", height =5, width = 12.5)

################################################################################

