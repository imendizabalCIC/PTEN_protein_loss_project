#| Last change: 12/12/2024
#| Ivana Rondon-Lorefice

################################################################################
#|  SCRIPT TO PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS WITH LIMMA VOOM
################################################################################

#| Pipeline: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

################################################################################
#| LIBRARIES
################################################################################
suppressMessages(library(DESeq2))
suppressWarnings(library(edgeR, quietly = T))
suppressMessages(library(VennDiagram))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ashr))
suppressMessages(library(calibrate))
suppressMessages(library(genefilter))
suppressMessages(library(circlize))
suppressMessages(library(ggrepel))
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(variancePartition))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(corrr))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(VennDiagram))

#| For plots
theme_set(theme_classic()) 
################################################################################


################################################################################
#| DATA
################################################################################
dir.proj <- "X:/irondon/AC-12_RNAseq/04_DEGs/"
counts.file <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt"
setwd(dir.proj)

#| Counts data
counts_data <- read.table(counts.file, sep ="\t")
counts_data <- counts_data[,order(colnames(counts_data))]

#| Filtering low counts (Not recommended to filter for Limma Voom)
#filter_counts <- rowSums(counts_data>0) >= 0.5*ncol(counts_data)
#counts_data <- counts_data[filter_counts,]

#| Sample info
sample <- colnames(counts_data)
condition <- gsub("_[0-9]*","",sample)
month <- gsub("[A-Z]", "", condition) 
sample_info <- data.frame(sample = sample,
                          condition = condition,
                          month = month)
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info[order(sample_info$sample),]
write.table(sample_info, "Data/sample_info_AC-12_RNAseq.txt", sep ="\t")

#| Checking the same order
any(colnames(counts_data) == rownames(sample_info))

#| Checking that there are the same number of samples
any(colnames(counts_data) %in% rownames(sample_info))

#| Normalization
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_data_normalized <- counts(dds, normalized=TRUE)
counts_data_normalized <- log(counts_data_normalized + 1, base =2)

#| Load the genome for merging GeneID
genome_mouse<- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Data/geneInfo.txt", sep ="\t", header =F)
colnames(genome_mouse) <- c("GeneID", "gene_name", "gene_info")

#| Defining the condition
conditions <- factor(t(sample_info["condition"]))
################################################################################



############################### LIMMA VOOM #####################################

#| Create DGEList object
d0 <- DGEList(counts_data)

#| Calculate normalization factors
d0 <- calcNormFactors(d0)

#| Multidimensional scaling (MDS) plot
plotMDS(d0, col = as.numeric(conditions))

#| Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + conditions)

#| Voom
y <- voom(d0, mm, plot = T)

#| Filtering a bit more genes
cutoff <- 2.5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 
dim(d0)

#| With more filtering
y <- voom(d, mm, plot = T)

#| lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(y, mm)
head(coef(fit))

#| Comparison between Western vs Standard Diet
contr <- makeContrasts( conditionsKO6 - conditionsWT6, levels = colnames(coef(fit)))

#| Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#| Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) (see https://www.degruyter.com/doi/10.2202/1544-6115.1027)
tmp <- eBayes(tmp)

#| What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
res_dataframe_limma_voom <- top.table
colnames(res_dataframe_limma_voom) <- c("log2FC_limma_voom", "AveExpr", "t", "pvalues_limma_voom", "padj_limma_voom", "B")
res_dataframe_limma_voom$GeneID <- rownames(res_dataframe_limma_voom)
res_dataframe_limma_voom <- merge(res_dataframe_limma_voom, genome_mouse, by ="GeneID")

#| Saving the table
write.table(res_dataframe_limma_voom,"Results/Tables/Limma_voom_analysis_DEGS_results.txt", sep = "\t", row.names=T)

#| Parameters 
FC <- 2
FDR <- 0.05

dir_tag <- "KO6_vs_WT6"
tag_contrast <- dir_tag

#| VOLCANO PLOT LIMMA VOOM
res_dataframe_limma_voom$delabel <- NA
res_dataframe_limma_voom$delabel[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC)) | (res_dataframe_limma_voom$log2FC_limma_voom > (log2(FC)))) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] <- 
  res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC)) | (res_dataframe_limma_voom$log2FC_limma_voom > log2(FC))) & (res_dataframe_limma_voom$padj_limma_voom < FDR))]

#| Assigning direction to the DEGs
res_dataframe_limma_voom$DEGs_direction <- NA
res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(FC)) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] <- "UP"
res_dataframe_limma_voom$DEGs_direction[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC))) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] <- "DOWN"
res_dataframe_limma_voom$DEGs_direction[which(is.na(res_dataframe_limma_voom$DEGs_direction))] <- "No DEGs"

#| Plotting and saving
ggplot(res_dataframe_limma_voom, aes(x = log2FC_limma_voom, y = -log10(padj_limma_voom), color =DEGs_direction, label=delabel )) + 
  geom_point(size =4,alpha = 0.5) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_color_manual(values = c("darkslategray","darkgrey", "goldenrod3"),guide = "legend") +
  labs(color ="Direction") +
  xlab("Log2 Fold Change") +
  ylab("-log10(padj value)") +
  geom_text_repel(family ="serif") 
ggsave(paste("Results/Images/",dir_tag,"/VolcanoPlot_",tag_contrast,"_LimmaVoom_FC_",FC,"_FDR_",FDR,".pdf", sep =""), height = 5, width = 6)

################################################################################


################################################################################
#| WATER FALL PLOT LIMMA VOOM
################################################################################

DEGs <- res_dataframe_limma_voom[which((res_dataframe_limma_voom$padj_limma_voom < FDR) & (!is.na(res_dataframe_limma_voom$padj_limma_voom)) & ((res_dataframe_limma_voom$log2FC_limma_voom > log2(FC)) | (res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC))))),]
DEGs$Direction <- DEGs$log2FC_limma_voom
DEGs$Direction[which(DEGs$Direction < 0)] <- "Down"
DEGs$Direction[which(DEGs$Direction >= 0 & DEGs$Direction != "Down")] <- "Up"

data3 <- DEGs
data3 <- data3[which(!(duplicated(data3$gene_name))),]
data3$gene_name <- factor(data3$gene_name,levels = data3$gene_name[order(data3$log2FC_limma_voom, decreasing = FALSE)])

ggplot(data3, aes(x = gene_name, y = log2FC_limma_voom, fill =Direction)) +
  theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #coord_cartesian(ylim = c(-4.1,4.1)) +
  geom_bar( stat = "identity") +
  xlab("Gene name") +
  ylab("Log2 Fold Change") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values =c("darkslategray","goldenrod3"))
ggsave(paste("Results/Images/",dir_tag,"/WaterFallPlot_",tag_contrast, "_LimmaVoom_FC_",FC,"_FDR_",FDR,".pdf", sep =""), height = 5, width = 7)


#| Top up and down
data3 <- data3[order(data3$log2FC_limma_voom),]
data3$Direction[1:15]
data3$Direction

top_up_down <- rbind(data3[1:15,], data3[(dim(data3)[1]-15):dim(data3)[1],])
ggplot(top_up_down, aes(x = log2FC_limma_voom, y = gene_name, fill =Direction)) +
  theme(text=element_text(size=13,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  #coord_cartesian(ylim = c(-5,5)) +
  geom_bar( stat = "identity") +
  xlab("Log2 Fold Change") +
  ylab("Gene name") +
  scale_fill_manual(values =c("darkslategray","goldenrod3"))
ggsave(paste("Results/Images/",dir_tag,"/Plot_TOP_",tag_contrast, "_LimmaVoom_FC_",FC,"_FDR_",FDR,".pdf", sep =""), height = 5, width = 7)


################################################################################
#| HEATMAP LIMMA VOOM
################################################################################

counts_zscore <- scale(t(as.data.frame(counts_data_normalized)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)

#| Sorting values
res_dataframe_limma_voom <- res_dataframe_limma_voom[order(res_dataframe_limma_voom$log2FC_limma_voom, decreasing = T),]
up_DEGs <- res_dataframe_limma_voom$GeneID[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(FC)) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] 
up_DEGs_gene_name <- res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom > log2(FC)) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] 
res_dataframe_limma_voom <- res_dataframe_limma_voom[order(res_dataframe_limma_voom$log2FC_limma_voom, decreasing = T),]
down_DEGs <- res_dataframe_limma_voom$GeneID[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC))) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] 
down_DEGs_gene_name <- res_dataframe_limma_voom$gene_name[which((res_dataframe_limma_voom$log2FC_limma_voom < (-log2(FC))) & (res_dataframe_limma_voom$padj_limma_voom < FDR))] 

#| Expression of Up and Down regulated
subset_genes <- counts_zscore[which(rownames(counts_zscore) %in% c(up_DEGs, down_DEGs)),]
subset_genes <- data.frame(subset_genes)
subset_genes$GeneID <- rownames(subset_genes)
subset_genes <- merge(subset_genes, genome_mouse, by ="GeneID")

#| In case of having repeated gene_names, aggregate its expression (with mean)
subset_genes <- aggregate(subset_genes, by=list(subset_genes$gene_name), mean)
rownames(subset_genes) <- subset_genes$Group.1

#| Subsetting
sample_col <- sample_info[which( (sample_info$condition == "KO6") | (sample_info$condition == "WT6")),]
subset_genes <- subset(subset_genes, select = sample_col$sample)
sample_col <- data.frame(condition = sample_col$condition)

#| ComplexHeatmap library. Convert gene expression to matrix form
subset_genes_matrix <- as.matrix(subset_genes)

#| To change colors 
ha <- HeatmapAnnotation(df = sample_col, col =list(condition = c(`KO6` =  "#440154FF", `WT6` = "#FDE725FF")))

FONT <- 0.2
pdf(paste("Results/Images/",dir_tag,"/HeatMap_",tag_contrast, "_LimmaVoom_FC_",FC,"_FDR_",FDR,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col = colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = F,
             row_names_gp = gpar(fontsize = FONT),
             column_names_gp = gpar(fontsize = 10),
             row_names_side="left")
draw(ht, newpage = F)
popViewport()
dev.off()


pdf(paste("Results/Images/",dir_tag,"/HeatMap_",tag_contrast, "_LimmaVoom_rowclustered_FC_",FC,"_FDR_",FDR,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col = colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = F,
             row_names_gp = gpar(fontsize = FONT),
             column_names_gp = gpar(fontsize = 5),
             row_names_side="left")
draw(ht, newpage = F)
popViewport()
dev.off()

pdf(paste("Results/Images/",dir_tag,"/HeatMap_",tag_contrast, "_LimmaVoom_columnsclustered_FC_",FC,"_FDR_",FDR,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col = colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = F, cluster_columns = T,
             row_names_gp = gpar(fontsize = FONT),
             column_names_gp = gpar(fontsize = 5),
             row_names_side="left")
draw(ht, newpage = F)
popViewport()
dev.off()

pdf(paste("Results/Images/",dir_tag,"/HeatMap_",tag_contrast, "_LimmaVoom_columnsclustered_columns_rows_clustered_FC_",FC,"_FDR_",FDR,".pdf", sep =""))
pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ht = Heatmap(subset_genes_matrix,
             name="Expression",
             col = colorRampPalette(c(brewer.pal(4, "Spectral")))(100),
             top_annotation =ha,
             cluster_rows = T, cluster_columns = T,
             row_names_gp = gpar(fontsize = FONT),
             column_names_gp = gpar(fontsize = 5),
             row_names_side="left")
draw(ht, newpage = F)
popViewport()
dev.off()
