################################################################################
#| LIMMA VOOM TO TEST DIFFERENTIALLY EXPRESSED GENES BETWEEN PTEN PROTEIN LOSS AND PRESENCE                    
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script applies the limma-voom pipeline to identify differentially expressed 
#| genes (DEGs) between PTEN protein loss and presence in the AC-45_RNAseq-FFPE 
#| dataset. It integrates preprocessing, DEG detection, visualization, and functional 
#| enrichment to characterize expression changes.
#|
#| Workflow:
#|   1) Load raw count matrix and sample metadata.
#|   2) Filter and normalize counts using edgeR/voom.
#|   3) Fit linear models (limma) including covariates (Age, DV200).
#|   4) Test contrasts for PTEN loss vs presence.
#|   5) Generate DEG tables with log2FC, p-values, and adjusted p-values.
#|   6) Visualize results with volcano plots, barplots (log2FC), and heatmaps.
#|   7) Perform enrichment analysis (GO, KEGG, Reactome) separately for up- and 
#|      down-regulated genes.
#|
#| Outputs:
#|   - Tables/Limma_voom_analysis_DEGS_results.txt (DEG table)
#|   - Volcano and bar plots of DEGs (PDF)
#|   - Heatmaps of DEG expression (PDF, clustered and unclustered versions)
#|   - Enrichment analysis plots for DEGs (GO, KEGG, Reactome)
#|
#| Reference pipeline: 
#|   https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
################################################################################



################################ LIBRARIES #####################################
suppressMessages(library(DESeq2))
suppressWarnings(library(edgeR, quietly = T))
suppressMessages(library(viridis))
suppressMessages(library(VennDiagram))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(apeglm))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(ashr))
suppressMessages(library(calibrate))
suppressMessages(library(vsn))
suppressMessages(library(genefilter))
suppressMessages(library(gplots))
suppressMessages(library(circlize))
suppressMessages(library(ggrepel))
suppressMessages(library(variancePartition))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(limma))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(corrr))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))

#| For plots
theme_set(theme_classic())
################################################################################


################################################################################
#| DATA 
################################################################################
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Limma_Analysis/4_Limma_Voom_PTEN_loss_vs_presence_Cohort/Results/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt"
setwd(dir.proj)

#| Readtable Full counts.
counts_data <- read.table(counts.file, sep = "\t", stringsAsFactors = TRUE)

#| Readtable Sample information
sample_info <- read.table(info.file, sep ="\t")

#| Filtering PTEN loss vs presence
sample_info <- sample_info[which(!is.na(sample_info$H_score_cut_0)),]
counts_data <- counts_data[, sample_info$AC.basurto]

#| Modifying name 
sample_info$PTEN_status <- NA
sample_info$PTEN_status[which(sample_info$H_score_cut_0 == "PTEN loss")] <-"loss"
sample_info$PTEN_status[which(sample_info$H_score_cut_0 == "PTEN presence")] <-"presence"

#| Defining the conditions: PTEN protein loss vs presence
conditions <- factor(t(sample_info["PTEN_status"]))
DV200 <- as.numeric(t(sample_info["DV200_Zscore"]))
Age <- as.numeric(t(sample_info["Edad_Zscore"]))

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)
genome_GRCh39.94 <- subset(genome_GRCh39.94, select =c("GeneID", "gene_name"))

length(sample_info$PTEN_status[which(sample_info$H_score_cut_0 == "PTEN loss")])
length(sample_info$PTEN_status[which(sample_info$H_score_cut_0 == "PTEN presence")])
################################################################################


################################################################################
#| LIMMA VOOM 
################################################################################
#| Create DGEList object
d0 <- DGEList(counts_data)

#| Calculate normalization factors
d0 <- calcNormFactors(d0)

#| Multidimensional scaling (MDS) plot
plotMDS(d0, col = as.numeric(conditions))

#| Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + DV200 + Age + conditions)

#| Voom
y <- voom(d0, mm, plot = T)

#| Filtering a bit more genes
cutoff <- 3.6
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

#| With more filtering
y <- voom(d, mm, plot = T)

#| lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))

#| Comparison between PTEN loss vs presence
contr <- makeContrasts( conditionsloss - conditionspresence, levels = colnames(coef(fit)))

#| Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#| Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) (see https://www.degruyter.com/doi/10.2202/1544-6115.1027)
tmp <- eBayes(tmp)

#| What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
res_dataframe_limma_voom <- top.table
colnames(res_dataframe_limma_voom) <- c("log2FC_limma_voom", "AveExpr", "t", "pvalues_limma_voom", "padj_limma_voom", "B")
res_dataframe_limma_voom$GeneID <- rownames(res_dataframe_limma_voom)
res_dataframe_limma_voom <- merge(res_dataframe_limma_voom, genome_GRCh39.94, by ="GeneID")

write.table(res_dataframe_limma_voom,"Tables/Limma_voom_analysis_DEGS_results.txt", sep = "\t", row.names=T)

res_dataframe_limma_voom <- read.table("Tables/Limma_voom_analysis_DEGS_results.txt", sep = "\t")
################################################################################


################################################################################
#| VOLCANO PLOT 
################################################################################

#| VOLCANO PLOT LIMMA VOOM
res_dataframe_limma_voom$delabel <- NA
res_dataframe_limma_voom$delabel[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1)) | (res_dataframe_limma_voom$log2FC_limma_voom > (log2(1)))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] <- 
  res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1)) | (res_dataframe_limma_voom$log2FC_limma_voom > log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))]

#| Assigning direction to the DEGs
res_dataframe_limma_voom$DEGs_direction <- NA
res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] <- "UP"
res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] <- "DOWN"
res_dataframe_limma_voom$DEGs_direction[which(is.na(res_dataframe_limma_voom$DEGs_direction))] <- "No DEGs"

length(res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))])
length(res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))])

#| Plotting and saving
ggplot(res_dataframe_limma_voom, aes(x = log2FC_limma_voom, y = -log10(padj_limma_voom), color =DEGs_direction, label=delabel )) + 
  geom_point(size =3.5,alpha = 0.4) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=14, hjust = 0.5, face ="bold"),
        legend.position = "top") +
  scale_color_manual(values = c("aquamarine4","darkgrey", "lightgoldenrod3") ,guide = "legend") +
  labs(color ="Direction") +
  xlab("Log2 Fold Change") +
  ylab("-log10(padj value)") #+geom_text_repel(family ="sans", size =0.9,force = 13,min.segment.length = 0.1)
ggsave("Images/Volcano_Plot_LimmaVoom.pdf", height = 4.5, width = 5.5)

?geom_text_repel

length(res_dataframe_limma_voom$gene_name[which( (res_dataframe_limma_voom$padj_limma_voom < 0.05) & (res_dataframe_limma_voom$log2FC_limma_voom < 0))])
length(res_dataframe_limma_voom$gene_name[which( (res_dataframe_limma_voom$padj_limma_voom < 0.05) & (res_dataframe_limma_voom$log2FC_limma_voom > 0))])

#| WATER FALL PLOT LIMMA VOOM
DEGs <- res_dataframe_limma_voom[which((res_dataframe_limma_voom$padj_limma_voom < 0.05) & (!is.na(res_dataframe_limma_voom$padj_limma_voom)) & ((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) | (res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))))),]
DEGs$Direction <- DEGs$log2FC_limma_voom
DEGs$Direction[which(DEGs$Direction < 0)] <- "Down"
DEGs$Direction[which(DEGs$Direction >= 0 & DEGs$Direction != "Down")] <- "Up"

data3 <- DEGs
data3 <- data3[which(!(duplicated(data3$gene_name))),]
data3$gene_name <- factor(data3$gene_name,levels = data3$gene_name[order(data3$log2FC_limma_voom, decreasing = FALSE)])

ggplot(data3, aes(x = gene_name, y = log2FC_limma_voom, fill =Direction)) +
  theme(text=element_text(size=16,  family="sans"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #coord_cartesian(ylim = c(-4.1,4.1)) +
  geom_bar( stat = "identity") +
  xlab("Genes") +
  ylab("Log2 Fold Change") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values =c("aquamarine4", "lightgoldenrod3"))
ggsave("Images/Log2FC_LimmaVoom_PTEN_protein_loss_vs_presence.pdf", height = 4.5, width = 6.6)

#| Top up and down
data3 <- data3[order(data3$log2FC_limma_voom),]
data3$Direction[1:15]
data3$Direction

top_up_down <- rbind(data3[1:17,], data3[(dim(data3)[1] -20):dim(data3)[1],])
ggplot(top_up_down, aes(x = log2FC_limma_voom, y = gene_name, fill =Direction)) +
  theme(text=element_text(size=13,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #coord_cartesian(ylim = c(-5,5)) +
  geom_bar( stat = "identity") +
  xlab("Log2 Fold Change") +
  ylab("Gene name") +
  scale_fill_manual(values =c("aquamarine4", "lightgoldenrod3"))
ggsave("Images/Log2FC_TOP_LimmaVoom_PTEN_protein_loss_vs_presence.pdf", height = 7, width = 8)
################################################################################


############################### HEAT MAP #######################################
counts_zscore <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Normalized_Log2_Counts/FullCounts_Basurto_Filtered_CPM_phenotype_Normalized_Log2.txt")
counts_zscore <- scale(t(as.data.frame(counts_zscore)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)

#| Sorting values
res_dataframe_limma_voom <- res_dataframe_limma_voom[order(res_dataframe_limma_voom$log2FC_limma_voom, decreasing = T),]
up_DEGs <- res_dataframe_limma_voom$GeneID[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] 
up_DEGs_gene_name <- res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] 
res_dataframe_limma_voom <- res_dataframe_limma_voom[order(res_dataframe_limma_voom$log2FC_limma_voom, decreasing = T),]
down_DEGs <- res_dataframe_limma_voom$GeneID[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] 
down_DEGs_gene_name <- res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))] 

#| Expression of Up and Down regulated
subset_genes <- counts_zscore[which(rownames(counts_zscore) %in% c(up_DEGs, down_DEGs)),]
subset_genes <- data.frame(subset_genes)
subset_genes$GeneID <- rownames(subset_genes)
subset_genes <- merge(subset_genes, genome_GRCh39.94, by ="GeneID")

#| Aggregating values to discard duplicated gene names
subset_genes <- aggregate(subset_genes, by=list(subset_genes$gene_name), mean)
rownames(subset_genes) <- subset_genes$Group.1

#| Ordering and selecting only PTEN protein samples
sample_info <- sample_info[order(sample_info$H_score_cut_0),]
subset_genes <- subset(subset_genes, select = sample_info$AC.basurto)

#| Renaming values
sample_col <- subset(sample_info, select =c("Edad", "purity", "stromal", "immune", "PTEN_status" ))
colnames(sample_col) <- c("Age", "Purity", "Stromal", "Immune", "PTEN_status")

overlapped_matrix <- as.matrix(subset_genes)

ha <- HeatmapAnnotation(df = sample_col, col =list(PTEN_status = c(`loss` =  "#EFD5B2", `presence` = "#504667")))

pdf("Images/Heatmap_LimmaVoom_DEGs_PTEN_protein_presence_vs_loss.pdf")
pushViewport(viewport(gp = gpar(fontfamily = "serif")))
ht = Heatmap(overlapped_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = 1.5),
             column_names_gp = gpar(fontsize = 2),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf("Images/Heatmap_LimmaVoom_DEGs_PTEN_protein_presence_vs_loss_columnsclustered.pdf")
pushViewport(viewport(gp = gpar(fontfamily = "serif")))
ht = Heatmap(overlapped_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = T,
             row_names_gp = gpar(fontsize = 1.5),
             column_names_gp = gpar(fontsize = 2),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf("Images/Heatmap_LimmaVoom_DEGs_PTEN_protein_presence_vs_loss_rowsclustered.pdf")
pushViewport(viewport(gp = gpar(fontfamily = "serif")))
ht = Heatmap(overlapped_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = F,
             row_names_gp = gpar(fontsize = 1.5),
             column_names_gp = gpar(fontsize = 2),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

pdf("Images/Heatmap_LimmaVoom_DEGs_PTEN_protein_presence_vs_loss_columns_rows_clustered.pdf")
pushViewport(viewport(gp = gpar(fontfamily = "serif")))
ht = Heatmap(overlapped_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = T,
             row_names_gp = gpar(fontsize = 1.5),
             column_names_gp = gpar(fontsize = 2),
             row_names_side="left")
draw(ht, newpage = FALSE)
popViewport()
dev.off()

#| Extracting the 30 most up-regulated and down regulated
res_dataframe_limma_voom <- res_dataframe_limma_voom[order(res_dataframe_limma_voom$log2FC_limma_voom, decreasing = T),]

#| UP-regulated
up_DEGs_dataframe <- res_dataframe_limma_voom[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(1)) & (res_dataframe_limma_voom$padj_limma_voom < 0.05)),] 
up_DEGs_dataframe <- up_DEGs_dataframe[order(up_DEGs_dataframe$log2FC_limma_voom, decreasing =T),]

up_DEGs <- up_DEGs_dataframe$GeneID[1:30]
up_DEGs_gene_name <- up_DEGs_dataframe$gene_name[1:30]


#| DOWN-regulated
down_DEGs_dataframe <- res_dataframe_limma_voom[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05)),] 
down_DEGs_dataframe <- down_DEGs_dataframe[order(down_DEGs_dataframe$log2FC_limma_voom, decreasing =F),]

down_DEGs <- down_DEGs_dataframe$GeneID[1:30]
down_DEGs_gene_name <- down_DEGs_dataframe$gene_name[1:30]


#| Expression of Up and Down regulated
subset_genes <- counts_zscore[which(rownames(counts_zscore) %in% c(up_DEGs, down_DEGs)),]
subset_genes <- data.frame(subset_genes)
subset_genes$GeneID <- rownames(subset_genes)
subset_genes <- merge(subset_genes, genome_GRCh39.94, by ="GeneID")

#| Aggregating values to discard duplicated gene names
subset_genes <- aggregate(subset_genes, by=list(subset_genes$gene_name), mean)
rownames(subset_genes) <- subset_genes$Group.1

#| Ordering and selecting only PTEN protein samples
sample_info <- sample_info[order(sample_info$H_score_cut_0),]
subset_genes <- subset(subset_genes, select = sample_info$AC.basurto)

#| Renaming values
sample_col <- subset(sample_info, select =c("Edad", "purity", "stromal", "immune", "PTEN_status" ))
sample_col <- sample_col[order(sample_col$PTEN_status, decreasing = T),]

colnames(sample_col) <- c("Age", "Purity", "Stromal", "Immune", "PTEN_status")

subset_genes <- subset_genes[, rownames(sample_col)]

overlapped_matrix <- as.matrix(subset_genes)

ha <- HeatmapAnnotation(df = sample_col, col =list(PTEN_status = c(`loss` =  "#EFD5B2", `presence` = "#504667")))

pdf("Images/Heatmap_LimmaVoom_DEGs_PTEN_protein_presence_vs_loss_rows_clustered_Top_30_DEGs.pdf")
pushViewport(viewport(gp = gpar(fontfamily = "serif")))
ht = Heatmap(overlapped_matrix,
             name="Expression",
             col =colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, 
             cluster_columns = F,
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 2),
             row_names_side="left",
             row_dend_reorder=F)
draw(ht, newpage = FALSE)
popViewport()
dev.off()
################################################################################


################################################################################
#| ENRICHMENT ANALYSIS 
################################################################################

#| Enrichment DEGs UP
genes_DEGs_UP <- res_dataframe_limma_voom$GeneID[which(((res_dataframe_limma_voom$log2FC_limma_voom > (log2(1)))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))]
total_gost_DEGs_UP <- gost(list("Module DEGs Limma Voom" = genes_DEGs_UP), 
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = res_dataframe_limma_voom$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs_UP, interactive = FALSE, capped =FALSE)
results_DEGs_UP <- total_gost_DEGs_UP$result[order(total_gost_DEGs_UP$result$p_value),]
results_DEGs_UP$`Term name` <- paste(results_DEGs_UP$term_name, "\n (N = ",results_DEGs_UP$term_size, ")",sep ="")


#| Top 15 most significant
ggplot(results_DEGs_UP[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")

dat1_filtered <- results_DEGs_UP[which(results_DEGs_UP$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Images/LimmaVoom/UP_enrichment_REAC_LimmaVoom.pdf", height=3.5, width = 4)

dat1_filtered <- results_DEGs_UP[which(results_DEGs_UP$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Images/LimmaVoom/UP_enrichment_KEGG_LimmaVoom.pdf", height=3.5, width = 3.2)

dat1_filtered <- results_DEGs_UP[which(results_DEGs_UP$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Images/LimmaVoom/UP_enrichment_KEGG_LimmaVoom.pdf", height=3.5, width = 3.2)



#| Enrichment DEGs DOWN
genes_DEGs_DOWN <- res_dataframe_limma_voom$GeneID[which(((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(1)))) & (res_dataframe_limma_voom$padj_limma_voom < 0.05))]
total_gost_DEGs_DOWN <- gost(list("Module DEGs Limma Voom" = genes_DEGs_DOWN), 
                           organism = "hsapiens", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = res_dataframe_limma_voom$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost_DEGs_DOWN, interactive = FALSE, capped =FALSE)
results_DEGs_DOWN <- total_gost_DEGs_DOWN$result[order(total_gost_DEGs_DOWN$result$p_value),]
results_DEGs_DOWN$`Term name` <- paste(results_DEGs_DOWN$term_name, "\n (N = ",results_DEGs_DOWN$term_size, ")",sep ="")


#| Top 15 most significant
ggplot(results_DEGs_DOWN[1:13,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="serif"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Images/LimmaVoom/DOWN_enrichment_LimmaVoom.pdf", height=4, width = 5)

dat1_filtered <- results_DEGs_DOWN[which(results_DEGs_DOWN$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = source, y = `Term name`)) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=11,  family="serif"), legend.key.size = unit(0.8, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis() +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Magenta REAC") +
  labs(size= "Intersection size")+
  xlab("Data")
ggsave("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Images/LimmaVoom/UP_enrichment_REAC_LimmaVoom.pdf", height=5.5, width = 5)

################################################################################

