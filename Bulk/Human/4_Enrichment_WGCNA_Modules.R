################################################################################
#| ENRICHMENT ANALYSIS OF WGCNA MODULES AND DEGs ASSOCIATED WITH PTEN LOSS                             
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script performs enrichment analyses on WGCNA modules significantly 
#| associated with PTEN protein loss (purple and green). It integrates 
#| differential expression results (DEGs) and module membership to identify 
#| functionally relevant pathways, regulators, and biological processes.
#|
#| Workflow:
#|   1) Load WGCNA module assignments, DEGs, normalized counts, and sample traits.  
#|   2) Calculate the overlap (Jaccard index) between DEGs and modules.  
#|   3) Perform enrichment analysis (gprofiler2) for selected modules (purple, green) 
#|      and their DEG intersections across databases (GO, KEGG, Reactome, TF, HPA).  
#|   4) Visualize enrichment results (bubble plots, dot plots).  
#|   5) Score module signatures at the single-cell level using UCell (on annotated scRNA-seq).  
#|   6) Compare gene expression (log2 normalized counts) between PTEN loss vs presence 
#|      samples within modules.  
#|   7) Generate heatmaps of DEGs and enriched gene sets (ComplexHeatmap).  
#|
#| Outputs:
#|   - Jaccard index plots (module-DEG overlap).  
#|   - Enrichment plots for purple and green modules (all terms, REAC, KEGG, TF, HPA).  
#|   - DEG intersection signatures (tables and enrichment maps).  
#|   - UMAP and violin plots showing module signatures in single-cell data.  
#|   - Boxplots comparing module DEG expression by PTEN status.  
#|   - Heatmaps of module DEGs and ECM-related genes.  
################################################################################


################################################################################
#| LIBRARIES AND DATA
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(vctrs ))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(viridis))
suppressMessages(library(readxl))
suppressMessages(library(colorspace))
suppressMessages(library(ggpubr))
suppressMessages(library(DESeq2))
suppressMessages(library(cowplot))
suppressMessages(library(UCell))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(colorspace))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))

#| For plots
theme_set(theme_classic())

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Pathways directories
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/"
data.file_WGCNA <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_intramodularinfo_last_try.txt"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Tables/Limma_voom_analysis_DEGS_results.txt"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Images/7_Enrichment_Analysis_of_Modules/"
tables.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Tables/"
ligands.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Data/Human-2020-Shao-LR-pairs.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt"
sample.file_immuno <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Low_High_PTEN_protein_samples_basurto_PI3K-AKT-mTOR_H_score_0.xlsx"

#| Setting working directory
setwd(workingDir)

#| WGCNA data
wgcna_Modules <- read.table(data.file_WGCNA, sep = "\t")
module_colors <- unique(wgcna_Modules$moduleColors)

wgcna_Modules <- merge(wgcna_Modules, genome_GRCh39.94, by ="GeneID")
rownames(wgcna_Modules) <- wgcna_Modules$GeneID

FC <- 1
FDR <- 0.05

#| DEGs by comparing PTEN loss vs presence Basurto
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")

#| Limma voom
DEGs <- DEGs_genes[which((DEGs_genes$log2FC_limma_voom > 0 | DEGs_genes$log2FC_limma_voom < 0) & DEGs_genes$padj_limma_voom < 0.05 & !is.na(DEGs_genes$padj_limma_voom)),]

#| Loading WGCNA data preprocessed
load(file ="X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Data/AC-45_RNAseq-FFPE_dataInput_last_try.RData")

#| Opening the processed data
loc.combined.sct <- readRDS("X:/irondon/Project_AC-45_RNAseq-FFPE/SingleCell/Results/Data/loc_combined_sct_annotated.rds")

#| Obtaining raw counts
counts_data <- read.table(counts.file, header=T) 
counts_data <- counts_data[,rownames(datTraits)]

#| Filtering 
filter_counts <- rowSums(counts_data>5) >= 0.7*ncol(counts_data)
counts_data <- counts_data[filter_counts,]

#| DESeq2 normalization (media of ratios)
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = datTraits, design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_data_normalized <- counts(dds, normalized=TRUE)
counts_data_normalized <- log(counts_data_normalized + 1, base =2)
counts_data_normalized <- data.frame(counts_data_normalized)
counts_data_normalized$GeneID <- rownames(counts_data_normalized)
counts_data_normalized <- merge(counts_data_normalized, genome_GRCh39.94, by ="GeneID")
counts_data_normalized <- aggregate(counts_data_normalized, by =list(counts_data_normalized$gene_name), mean)
rownames(counts_data_normalized) <- counts_data_normalized$Group.1
counts_data_normalized <- counts_data_normalized[,rownames(datTraits)]

#| Z-score of normalized log2 counts
counts_zscore <- scale(t(as.data.frame(counts_data_normalized)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)
counts_zscore <- as.data.frame(counts_zscore)

#| Sample info of patients 
sample_info_immuno <- read_xlsx(sample.file_immuno)
################################################################################


################################################################################
#|  FUNCTIONS
################################################################################

#| Jaccard index
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

################################################################################


################################################################################
#|  % OF DEGS CONTAINED IN THE DIFFERENT MODULES
################################################################################

DEGs_up <- DEGs_genes[which((DEGs_genes$log2FC_limma_voom > 0) & DEGs_genes$padj_limma_voom < 0.05 & !is.na(DEGs_genes$padj_limma_voom)),]
DEGs_down <- DEGs_genes[which((DEGs_genes$log2FC_limma_voom < 0) & DEGs_genes$padj_limma_voom < 0.05 & !is.na(DEGs_genes$padj_limma_voom)),]

colors_up <- unique(wgcna_Modules$moduleColors[which(wgcna_Modules$GeneID %in% DEGs_up$GeneID)])
colors_down <- unique(wgcna_Modules$moduleColors[which(wgcna_Modules$GeneID %in% DEGs_down$GeneID)])

jaccard_value_up <- c()
jaccard_value_down <- c()

size_up <- c()
size_down <- c()

for (i in 1:length(colors_up)){

  jaccard_value_up <- c(jaccard_value_up,jaccard(DEGs_up$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_up[i])]))
  size_up <- c(size_up, length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors_up[i])]))
  
}

for (i in 1:length(colors_down)){
  
  jaccard_value_down <- c(jaccard_value_down,jaccard(DEGs_down$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_down[i])]))
  size_down <- c(size_down, length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors_down[i])]))
  
}

data_percentage_DGEs_up <- data.frame(modules = colors_up,
                                   jaccard_value = jaccard_value_up,
                                   size = size_up,
                                   seq = jaccard_value_up)
data_percentage_DGEs_down <- data.frame(modules = colors_down,
                                   jaccard_value = -jaccard_value_down,
                                   size = size_down,
                                   seq = jaccard_value_down)

data_percentage_DGEs <- rbind(data_percentage_DGEs_up, data_percentage_DGEs_down)

#| Plotting the results
ggplot(data_percentage_DGEs, aes(x =jaccard_value,y = reorder(modules,seq),  fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21, color ="black")+
  theme(axis.text.x = element_text(vjust = 0.99, hjust=1),text=element_text(size=12,  family="sans"), legend.key.size = unit(0.1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors[order(colors)]) +
  #scale_color_manual(values = colors[order(colors)])+
  ylab("Modules WGCNA (147 PCa samples)") +
  xlab("Jaccard index")+
  guides(fill ="none")+
  labs(size = "Module size") 
ggsave(paste(results.file, "Jaccard_index_WGCNA_DEGs_limma_voom_up_down.pdf", sep =""), heigh =4.5, width = 6 )

#| Plotting the results
ggplot(data_percentage_DGEs_up, aes(x =jaccard_value,y = reorder(modules,jaccard_value),  fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21, color ="black")+
  theme(axis.text.x = element_text(vjust = 0.99, hjust=1),text=element_text(size=12,  family="sans"), legend.key.size = unit(0.1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors_up[order(colors_up)]) +
  #ylab("Modules WGCNA (147 PCa samples)") +
  ylab("") +
  xlab("Jaccard index")+
  guides(fill ="none")+
  labs(size = "Module size")
ggsave(paste(results.file, "Jaccard_index_WGCNA_DEGs_limma_voom_up.pdf", sep =""), height = 3.1, width = 4.3 )
################################################################################


################################################################################
#| Module: PURPLE
################################################################################
mod <- "purple"
wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing = T),]
genes <- wgcna_Modules$GeneID[wgcna_Modules$moduleColors == mod]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "hsapiens", ordered_query = TRUE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_",mod,".pdf", sep =""), height = 4.5, width = 5.6)


#| Top with REAC
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 5)

#| Top with HPA
dat1_filtered <- results[which(results$source == "HPA"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"HPA_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 3.6)

#| Top with TF
dat1_filtered <- results[which(results$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"TF_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4.4)
##############################################################################


################################################################################
#| Module: GREEN
################################################################################

mod <- "green"
wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing = T),]
genes <- wgcna_Modules$GeneID[wgcna_Modules$moduleColors == mod]
genes_name <- wgcna_Modules$gene_name[wgcna_Modules$moduleColors == mod]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "hsapiens", ordered_query = TRUE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_",mod,".pdf", sep =""), height = 4.5, width = 5)


#| Top with REAC
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 3.5)


#| Top with HPA
dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"KEGG_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 3.9)

#| Top with TF
dat1_filtered <- results[which(results$source == "TF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"TF_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4)

#| UCell
signatures <- list()
signatures$mod_green <- genes_name

DefaultAssay(loc.combined.sct) <- "RNA"

#| SCORING SIGNATURES USING UCELL
loc.combined.sct <- AddModuleScore_UCell(loc.combined.sct, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#| Module purple and DEGs label TRUE
p <- FeaturePlot(loc.combined.sct, reduction = "umap", 
                 features = signature.names[1], 
                 ncol = 1, 
                 order = T,
                 keep.scale = "all",
                 min.cutoff = 'q10',
                 label = TRUE,
                 label.color = "black",
                 cols = c("seashell1","cornflowerblue"),
                 label.size = 3.5,
                 repel = T) 
p <- p + theme(text=element_text(size=12,  family="sans"), 
               legend.key.size = unit(1, 'cm'), 
               plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
p + labs(title = "Module Green")
ggsave(file= paste(results.file,"UMAP_Chen_Module_Green.pdf", sep =""), width = 7, height = 5)

##############################################################################


################################################################################
#| Module Purple intersection with DEGs
################################################################################

#| Size of module purple
mod <- "Purple_DEGs"
wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing = T),]
genes <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "purple")], DEGs$GeneID)
genes_name <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "purple")], DEGs$gene_name)
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]

source <- unique(results$source)
source <- c("TF", "GO:CC", "GO:BP","HPA", "REAC")

n <- 6
value <- results[which(results$source == source[1]),c("term_id", "term_name", "p_value", "intersection", "source","intersection_size")]
value <- value[1:n,]
  
for (i in 2:length(source)){
  value2 <- results[which(results$source == source[i]),c("term_id", "term_name", "p_value", "intersection", "source","intersection_size")]
  value2 <- value2[1:6,]
  
  value <- rbind(value,
                 value2)
}

value <- value[which(!is.na(value$term_id)),]

colnames(value) <- c("GO.ID", "Description", "p.Val", "Genes", "source","intersection_size")
value$FDR <- value$p.Val
value$Phenotype <- "+1"
value <- value[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes", "source","intersection_size")]
write.table(value, "Images/7_Enrichment_Analysis_of_Modules/Modules_signature/EnrichmentMap_DEGs_purple.txt", sep = "\t", quote = F, row.names = F)

value$source_name <- NA
value$source_name[which(value$source == "TF")]<- 1
value$source_name[which(value$source == "GO:CC")]<- 2
value$source_name[which(value$source == "GO:BP")]<- 5
value$source_name[which(value$source == "HPA")]<- 3
value$source_name[which(value$source == "REAC")]<- 4

value$Description <- toupper(value$Description)
value$Description_paste <-paste(value$source_name,  -log2(value$p.Val),value$Description , sep = " ")
value <-value[ order(value$Description_paste, decreasing = T),]
value$source_name
value <-value[ reorder( value$Description_paste , -value$source_name), ]

ggplot(value, aes(x = -log2(p.Val), y = reorder( Description_paste , -source_name))) +
  geom_point(aes(size = intersection_size, color = source)) +
  theme(text=element_text(size=9,  family="sans"), 
        legend.key.size = unit(0.6, 'cm'), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold"),
        axis.text.y = element_text(size=9)) +
  background_grid(major = "xy", minor = "none") +
  scale_y_discrete(labels = value$Description) +
  scale_color_manual(values = rainbow(5))+
  ylab("")

"#FF0000" "#CCFF00" "#00FF66" "#0066FF" "#CC00FF"
ggsave("Images/7_Enrichment_Analysis_of_Modules/Enrichment_together_gprofiler_DEGs_Purple.pdf", heigh=7, width =9)

#| Save DEGs in module Purple
write.table(data.frame(genes =genes, genes_name=genes_name), "Images/7_Enrichment_Analysis_of_Modules/Modules_signature/results_DEGs_Purple_enrichment_gprofiler.txt", sep ="\t")

results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_",mod,".pdf", sep =""), height = 4.5, width = 5.3)

#| Top with REAC
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), axis.text.y = element_text(size = 8, face ="bold", color ="black"),
        legend.key.size = unit(0.3, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.11, width = 4.1)


genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Signaling by NTRKs\n (N = 121)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Signaling by Receptor Tyrosine Kinases\n (N = 452)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Signaling by SCF-KIT\n (N = 37)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Signaling by NTRK1 (TRKA)\n (N = 103)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "PI-3K cascade:FGFR3\n (N = 8)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "PI-3K cascade:FGFR1\n (N = 8)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "PI-3K cascade:FGFR4\n (N = 9)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "PI-3K cascade:FGFR2\n (N = 10)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Tie2 Signaling\n (N = 14)")], split  = ",")))]

#| Top with KEGG
dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"KEGG_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 3.5)

#| UCell
signatures <- list()
signatures$mod_Purple_DEGs <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "purple")], DEGs$gene_name)

DefaultAssay(loc.combined.sct) <- "RNA"

#| SCORING SIGNATURES USING UCELL
loc.combined.sct <- AddModuleScore_UCell(loc.combined.sct, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#| Selecting the colors of label
colors <- data.frame(Celltype = unique(loc.combined.sct$annotation),
                cols = rainbow(length(unique(loc.combined.sct$annotation))),
                num = gsub("\\_.*$", "", unique(loc.combined.sct$annotation)))
colors$num <- as.numeric(colors$num)
colors <- colors[order(colors$num, decreasing= F),]

#| UMAP Module purple and DEGs label TRUE
p <- FeaturePlot(loc.combined.sct, reduction = "umap", 
                 features = signature.names[1], 
                 ncol = 1, 
                 order = T,
                 keep.scale = "all",
                 min.cutoff = 'q10',
                 label = TRUE,
                 label.color = colors$cols,
                 cols = c("seashell1","cornflowerblue"),
                 label.size = 3.5,
                 repel = T) 
p <- p + theme(text=element_text(size=12,  family="sans"), 
               legend.key.size = unit(1, 'cm'), 
               plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
p + labs(title = "Module Purple and DEGs") 
ggsave(file= paste(results.file,"UMAP_Chen_Module_Purple_DEGs.pdf", sep =""), width = 6.6, height = 4.5)

#| Violin plot
v <- VlnPlot(loc.combined.sct, 
        features = signature.names[1], 
        sort=T,
        cols = colors$cols)
v <- v + theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 90))
v + labs(title = "Module Purple and DEGs") 
ggsave(file= paste(results.file,"Vln_Chen_Module_Purple_DEGs.pdf", sep =""), width = 8, height = 5)

##############################################################################


################################################################################
#| Module Green intersection with DEGs
################################################################################
mod <- "Green_DEGs"
length(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")])
length(DEGs$GeneID)
length(intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID))

genes <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "green")], DEGs$GeneID)
genes_name <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = rownames(wgcna_Modules), 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]

source <- unique(results$source)
source <- c("TF", "GO:CC", "GO:BP","GO:MF", "KEGG", "REAC")

n <- 6
value <- results[which(results$source == source[1]),c("term_id", "term_name", "p_value", "intersection", "source","intersection_size")]
value <- value[1:n,]

for (i in 2:length(source)){
  value2 <- results[which(results$source == source[i]),c("term_id", "term_name", "p_value", "intersection", "source","intersection_size")]
  value2 <- value2[1:6,]
  
  value <- rbind(value,
                 value2)
}

value <- value[which(!is.na(value$term_id)),]

colnames(value) <- c("GO.ID", "Description", "p.Val", "Genes", "source","intersection_size")
value$FDR <- value$p.Val
value$Phenotype <- "+1"
value <- value[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes", "source","intersection_size")]
write.table(value, "Images/7_Enrichment_Analysis_of_Modules/Modules_signature/EnrichmentMap_DEGs_green.txt", sep = "\t", quote = F, row.names = F)

value$source_name <- NA
value$source_name[which(value$source == "TF")]<- 1
value$source_name[which(value$source == "GO:CC")]<- 2
value$source_name[which(value$source == "GO:BP")]<- 3
value$source_name[which(value$source == "GO:MF")]<- 4
value$source_name[which(value$source == "KEGG")]<- 5
value$source_name[which(value$source == "REAC")]<- 6

value$Description <- toupper(value$Description)
value$Description_paste <-paste(value$source_name,  -log2(value$p.Val),value$Description , sep = " ")
value <-value[ order(value$Description_paste, decreasing = T),]
value$source_name
value <-value[ reorder( value$Description_paste , -value$source_name), ]

ggplot(value, aes(x = -log2(p.Val), y = reorder( Description_paste , -source_name))) +
  geom_point(aes(size = intersection_size, color = source)) +
  theme(text=element_text(size=9,  family="sans"), 
        legend.key.size = unit(0.6, 'cm'), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold"),
        axis.text.y = element_text(size=9)) +
  background_grid(major = "xy", minor = "none") +
  scale_y_discrete(labels = value$Description) +
  scale_color_manual(values = rainbow(6))+
  ylab("")
ggsave("Images/7_Enrichment_Analysis_of_Modules/Enrichment_together_gprofiler_DEGs_Green.pdf", heigh=8, width =9.5)

results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Save DEGs in module Purple
write.table(data.frame(genes =genes, genes_name=genes_name), "Images/7_Enrichment_Analysis_of_Modules/Modules_signature/DEGs_green.txt", sep ="\t")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_",mod,".pdf", sep =""), height = 4.5, width = 5)

#| Top with REAC
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), 
        legend.key.size = unit(0.3, 'cm'),  axis.text.y = element_text(size = 8, face ="bold", color ="black"),
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.11, width = 5.5)

#| Genes contained in the Extracellular matrix enrichment
genes_ECM <- results$intersection[which(results$source == "REAC" & results$term_name == "Extracellular matrix organization")]
genes_ECM <- strsplit(genes_ECM, split =",")
genes_ECM <- genes_ECM[[1]]


genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Extracellular matrix organization\n (N = 248)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Cell-extracellular matrix interactions\n (N = 18)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Localization of the PINCH-ILK-PARVIN complex to focal adhesions\n (N = 4)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "ECM proteoglycans\n (N = 61)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Muscle contraction\n (N = 142)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Signaling by TGF-beta Receptor Complex\n (N = 88)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "RUNX3 regulates YAP1-mediated transcription\n (N = 8)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "TGF-beta receptor signaling activates SMADs\n (N = 44)")], split  = ",")))]
genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% unlist(strsplit(results$intersection[which(results$`Term name` == "Tie2 Signaling\n (N = 14)")], split  = ",")))]


#| Top with KEGG
dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"KEGG_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4)

#| UCell
signatures <- list()
signatures$mod_Green_DEGs <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)

DefaultAssay(loc.combined.sct) <- "RNA"

#| SCORING SIGNATURES USING UCELL
loc.combined.sct <- AddModuleScore_UCell(loc.combined.sct, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#| Selecting the colors of label
colors <- data.frame(Celltype = unique(loc.combined.sct$annotation),
                     cols = rainbow(length(unique(loc.combined.sct$annotation))),
                     num = gsub("\\_.*$", "", unique(loc.combined.sct$annotation)))
colors$num <- as.numeric(colors$num)
colors <- colors[order(colors$num, decreasing= F),]

#| Module green and DEGs label TRUE
p <- FeaturePlot(loc.combined.sct, reduction = "umap", 
                 features = signature.names[1], 
                 ncol = 1, 
                 order = T,
                 keep.scale = "all",
                 min.cutoff = 'q10',
                 label = TRUE,
                 label.color = colors$cols,
                 cols = c("seashell1","cornflowerblue"),
                 label.size = 3.5,
                 repel = T) 
p <- p + theme(text=element_text(size=12,  family="sans"), 
               legend.key.size = unit(1, 'cm'), 
               plot.title=element_text(size=16, hjust = 0.5, face ="bold"))
p + labs(title = "Module Green and DEGs")
ggsave(file= paste(results.file,"UMAP_Chen_Module_Green_DEGs.pdf", sep =""),  width = 6.6, height = 5)

#| Violin plot
v <- VlnPlot(loc.combined.sct, 
             features = signature.names[1], 
             sort=T,
             cols = colors$cols)
v <- v + theme(legend.position = 'none') +theme(axis.text.x = element_text(angle = 90))
v + labs(title = "Module Green and DEGs") 
ggsave(file= paste(results.file,"Vln_Chen_Module_Green_DEGs.pdf", sep =""), width = 8, height = 5)

##############################################################################


################################################################################
#| EXPRESSION UNDER CONDITIONS PTEN PROTEIN LOSS VS PRESENCE
################################################################################

#| Module purple
module_name <- "purple"
wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing = T),]
genes <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "purple")], DEGs$gene_name)

loss_presence_patients <- datTraits
loss_presence_patients <- cbind(loss_presence_patients,t(counts_data_normalized[genes[1:10],]))
loss_presence_patients <- stack(loss_presence_patients[, (dim(datTraits)[2]+1):(dim(loss_presence_patients)[2])])
loss_presence_patients$PTEN_status <- datTraits$PTEN_status

loss_presence_patients$PTEN_status[which(loss_presence_patients$PTEN_status == 1)] <- "loss"
loss_presence_patients$PTEN_status[which(loss_presence_patients$PTEN_status == 0)] <- "presence"

ggplot(loss_presence_patients, aes( x = ind, y = values, fill = PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="serif", size =3)+
  theme(text=element_text(size=15,  family="serif"), legend.key.size = unit(0.8, 'cm'),plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste("DEGs in module ",module_name, sep =""))+
  ylab("Normalized log2 counts") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))
ggsave(paste(results.file ,"Expression_DEGs_in_module_",module_name,"_PTEN_Protein_loss_vs_presence.pdf", sep = ""), heigh=5, width= 9.6)


#| Module green
module_name <- "green"
wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing =T),]
genes <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "green")], DEGs$gene_name)

loss_presence_patients <- datTraits
loss_presence_patients <- cbind(loss_presence_patients,t(counts_data_normalized[genes[1:10],]))
loss_presence_patients <- stack(loss_presence_patients[, (dim(datTraits)[2]+1):(dim(loss_presence_patients)[2])])
loss_presence_patients$PTEN_status <- datTraits$PTEN_status

loss_presence_patients$PTEN_status[which(loss_presence_patients$PTEN_status == 1)] <- "loss"
loss_presence_patients$PTEN_status[which(loss_presence_patients$PTEN_status == 0)] <- "presence"

ggplot(loss_presence_patients, aes( x = ind, y = values, fill = PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="serif", size =3)+
  theme(text=element_text(size=15,  family="serif"), legend.key.size = unit(0.8, 'cm'),plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  xlab(paste("DEGs in module ",module_name, sep =""))+
  ylab("Normalized log2 counts") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))
ggsave(paste(results.file ,"Expression_DEGs_in_module_",module_name,"_PTEN_Protein_loss_vs_presence.pdf", sep = ""), heigh=5, width= 9.6)

##############################################################################


################################################################################
#| HEATMAPS
################################################################################

####| Purple
name <- "green_DEGs"
module_name <- "green"
genes <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == module_name)], DEGs$gene_name)
subset_genes <- counts_zscore[which(rownames(counts_zscore)%in%genes),]
subset_genes <- data.frame(subset_genes)

datTraits$PTEN_status[which(datTraits$PTEN_status == 1)] <- "loss"
datTraits$PTEN_status[which(datTraits$PTEN_status == 0)] <- "presence"

#| Ordering and selecting only PTEN protein samples
datTraits <- datTraits[order(datTraits$PTEN_status),]
subset_genes <- subset(subset_genes, select = rownames(datTraits))

#| Renaming values
sample_col <- subset(datTraits, select =c( "Age", "stromal", "purity", "PTEN_status"))

#| ComplexHeatmap library. Convert gene expression to matrix form
subset_genes_matrix <- as.matrix(subset_genes)

col_fun_Age <- colorRamp2(c(min(sample_col$Age),  max(sample_col$Age)), c("white", "#243e36"))
col_fun_stromal <- colorRamp2(c(min(sample_col$stromal),  max(sample_col$stromal)), c("white", "#2a9d8f"))
col_fun_purity <- colorRamp2(c(min(sample_col$purity),  max(sample_col$purity)), c("white", "#e9c46a"))

#| To change colors 
ha <- HeatmapAnnotation(df = sample_col, col =list(PTEN_status = c(`loss` =  "#EFD5B2", `presence` = "#504667"),
                                                   Age = col_fun_Age,
                                                   stromal = col_fun_stromal,
                                                   purity = col_fun_purity))

pdf(paste(results.file, "Heatmap_", name,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 5.5),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf(paste(results.file, "Heatmap_", name,"_row_columns.pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = T,
             row_names_gp = gpar(fontsize = 5.5),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf(paste(results.file, "Heatmap_", name,"_row.pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = F,
             row_names_gp = gpar(fontsize = 5.5),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()


##################### ECM genes in module Green 
#| Note: Run first the enrichment in "Module Green intersection with DEGs"
subset_genes <- counts_zscore[which(rownames(counts_zscore) %in% genome_GRCh39.94$gene_name[which(genome_GRCh39.94$GeneID %in% genes_ECM)]),]
subset_genes <- data.frame(subset_genes)

datTraits$PTEN_status[which(datTraits$PTEN_status == 1)] <- "loss"
datTraits$PTEN_status[which(datTraits$PTEN_status == 0)] <- "presence"

#| Ordering and selecting only PTEN protein samples
datTraits <- datTraits[order(datTraits$PTEN_status),]
subset_genes <- subset(subset_genes, select = rownames(datTraits))
names(datTraits)
#| Renaming values
sample_col <- subset(datTraits, select =c( "Age", "stromal", "purity", "PTEN_status"))

#| ComplexHeatmap library. Convert gene expression to matrix form
subset_genes_matrix <- as.matrix(subset_genes)

#| To change colors 
ha <- HeatmapAnnotation(df = sample_col, col =list(PTEN_status = c(`loss` =  "#EFD5B2", `presence` = "#504667")))

col_fun_Age <- colorRamp2(c(min(sample_col$Age),  max(sample_col$Age)), c("white", "#243e36"))
col_fun_stromal <- colorRamp2(c(min(sample_col$stromal),  max(sample_col$stromal)), c("white", "#2a9d8f"))
col_fun_purity <- colorRamp2(c(min(sample_col$purity),  max(sample_col$purity)), c("white", "#e9c46a"))

#| To change colors 
ha <- HeatmapAnnotation(df = sample_col, col =list(PTEN_status = c(`loss` =  "#EFD5B2", `presence` = "#504667"),
                                                   Age = col_fun_Age,
                                                   stromal = col_fun_stromal,
                                                   purity = col_fun_purity))

name <- "ECM_genes_green"
pdf(paste(results.file, "Heatmap_", name,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 13),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf(paste(results.file, "Heatmap_", name,"_row_columns.pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = T,
             row_names_gp = gpar(fontsize = 13),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf(paste(results.file, "Heatmap_", name,"_row.pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = F,
             row_names_gp = gpar(fontsize = 13),
             column_names_gp = gpar(fontsize = 1),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()
################################################################################
