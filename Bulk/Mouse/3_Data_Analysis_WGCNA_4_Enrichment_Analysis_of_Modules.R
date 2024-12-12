################################################################################
#####             ENRICHMENT ANALYSIS OF MODULES DE FROM WGCNA             #####                        
################################################################################


################################################################################


############################  LIBRARIES AND DATA  ##############################
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
suppressMessages(library(VennDiagram))
suppressMessages(library(UCell))
suppressMessages(library(homologene))
package.version("homologene")

#| For plots
theme_set(theme_classic())

#| Genome mouse
genome_mouse<- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Data/geneInfo.txt", sep ="\t", header =F)
colnames(genome_mouse) <- c("GeneID", "gene_name", "gene_info")

#| GeneID from mouse to human
mouse_human <- read.table("X:/DATA_shared/Human_to_Mouse_Gene_name.txt", sep ="\t", header =T )

#| Pathways directories
dir.proj <- "X:/irondon/AC-12_RNAseq/07_WGCNA/"
data.file_WGCNA <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Data/3_geneInfo_counts_geneMM_geneS_moduleColors.txt"
data.file_DEGs <- "X:/irondon/AC-12_RNAseq/04_DEGs/Results/Tables/Limma_voom_analysis_DEGS_results.txt"
results.file <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Images/07_Enrichment_Analysis_of_Modules/"
setwd(dir.proj)

#| WGCNA data
wgcna_Modules <- read.table(data.file_WGCNA, sep = "\t")
wgcna_Modules$GeneID <- rownames(wgcna_Modules)
wgcna_Modules <- merge(wgcna_Modules, genome_mouse, by ="GeneID")
module_colors <- unique(wgcna_Modules$moduleColors)

FC <- 2
FDR <- 0.05

#| DEGs by comparing KO6 VS WT6
DEGs_genes <- read.table(data.file_DEGs, sep = "\t")

#| DESEQ2
#DEGs <- DEGs_genes[which( !(is.na(DEGs_genes$padj_deseq2)) & (DEGs_genes$padj_deseq2 < FDR ) & ( (DEGs_genes$log2FC_deseq2 > log2(FC)) | (DEGs_genes$log2FC_deseq2 < (-log2(FC)))) ),]

#| Limma voom
DEGs <- DEGs_genes[which(( (DEGs_genes$log2FC_limma_voom > log2(FC)) | (DEGs_genes$log2FC_limma_voom <  (-log2(FC)))) & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom))),]

#| UP
length(DEGs_genes$gene_name[which( (DEGs_genes$log2FC_limma_voom > log2(FC)) & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom)))])

#| DOWN
length(DEGs_genes$gene_name[which(( (DEGs_genes$log2FC_limma_voom <  (-log2(FC)))) & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom)))])

#| Opening the processed data of single cell
#loc.combined.sct <- readRDS("X:/irondon/Project_AC-45_RNAseq-FFPE/SingleCell/Results/Data/loc_combined_sct_annotated.rds")
################################################################################


############################ PLOT: SIZE OF MODULES #############################
#| Size
size <- c()
for (i in 1:length(module_colors)){
  size <- c(size,length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors ==module_colors[i])]) )
}

length(rownames(wgcna_Modules)[which(wgcna_Modules$moduleColors == "grey")])

#| New dataframe
data_size <- data.frame(modules = module_colors,
                        size = size)
colors <- sort(unique(wgcna_Modules$moduleColors))
ggplot(data_size, aes(x= size, y = reorder(modules,size), fill =modules)) +
  geom_bar(stat="identity", color ="black") +
  theme(axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        text=element_text(size=14,  family="sans"), 
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  xlab("Module size") +
  ylab("Module by color name")
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/5_Module_size.pdf", height = 4.5, width = 4.9)
################################################################################


####################### VOLCANO WITH MODULE COLORS #############################
#| Enrichment over all the modules
modules_color <- sort(module_colors)
dat_Modules_DEGs <- merge(wgcna_Modules, DEGs_genes, by ="GeneID") 
names(dat_Modules_DEGs)

dat_Modules_DEGs$delabel <- NA
dat_Modules_DEGs$delabel[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < FDR))] <- 
  dat_Modules_DEGs$gene_name.x[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < FDR))]

ggplot(dat_Modules_DEGs, aes(x = log2FC_limma_voom, y = -log10(padj_limma_voom), color=moduleColors)) + 
  geom_point(size =3,alpha = 0.5) +
  theme(text=element_text(size=16,  family="sans"), legend.key.size = unit(0.2, 'cm'), 
        legend.position = "none",
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_color_manual(values = modules_color ,guide = "legend") +
  labs(color ="Module color") +
  xlab("Log2 Fold Change")
ggsave(paste(results.file, "VolcanoPlots_with_Module_colors_LimmaVoom.pdf", sep =""), height = 5.5, width = 7)

dat_Modules_DEGs$delabel <- NA
dat_Modules_DEGs$delabel[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < FDR))] <- 
  dat_Modules_DEGs$gene_name.x[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < 0.05))]

dat_Modules_DEGs$moduleColors[which(is.na(dat_Modules_DEGs$delabel))] <- "grey90"
modules_color <- unique(dat_Modules_DEGs$moduleColors[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < 0.05))])
modules_color <- c(modules_color, "grey90")
modules_color <- sort(modules_color)

dat_Modules_DEGs$moduleColors[which(is.na(dat_Modules_DEGs$delabel))] <- "grey90"
modules_color <- unique(dat_Modules_DEGs$moduleColors[which( ( (dat_Modules_DEGs$log2FC_limma_voom < (-log2(FC))) | (dat_Modules_DEGs$log2FC_limma_voom > (log2(FC)))) & (dat_Modules_DEGs$padj_limma_voom < 0.05))])
modules_color <- c(modules_color, "grey90")
modules_color <- sort(modules_color)
ggplot(dat_Modules_DEGs, aes(x = log2FC_limma_voom, y = -log10(padj_limma_voom), color=moduleColors, label=delabel)) + 
  geom_point(size =2.5, alpha =0.4) +
  theme(text=element_text(size=12,  family="sans"), 
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black"),
        legend.position = "none",
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_color_manual(values =  modules_color) +
  labs(color ="Module color") +
  xlab("Log2 Fold Change")+
  ylab("Log(p-adjusted value)")+
  geom_hline(yintercept=(-log10(0.05)), color = "grey50", linetype='dashed', size=0.2)+
  geom_vline(xintercept=(-log2(2)), color = "grey50", linetype='dashed', size=0.2)+
  geom_vline(xintercept=(log2(2)), color = "grey50", linetype='dashed', size=0.2)
ggsave(paste(results.file, "VolcanoPlots_with_Module_colors_LimmaVoom_grey.pdf", sep =""), height = 4.5, width = 5)

################################################################################


################# % OF DEGS CONTAINED IN THE DIFFERENT MODULES #################  
#| Jaccard function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

DEGs_up <- DEGs_genes[which( (DEGs_genes$log2FC_limma_voom > log2(FC))  & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom))),]
DEGs_down <- DEGs_genes[which( (DEGs_genes$log2FC_limma_voom <  (-log2(FC))) & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom))),]

colors_up <- unique(wgcna_Modules$moduleColors[which(wgcna_Modules$GeneID %in% DEGs_up$GeneID)])
colors_down <- unique(wgcna_Modules$moduleColors[which(wgcna_Modules$GeneID %in% DEGs_down$GeneID)])

jaccard_value_up <- c()
jaccard_value_down <- c()

size_up <- c()
size_down <- c()

intersect_up <- c()
intersect_down <- c()

for (i in 1:length(colors_up)){
  jaccard_value_up <- c(jaccard_value_up,jaccard(DEGs_up$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_up[i])]))
  size_up <- c(size_up, length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors_up[i])]))
  intersect_up <- c(intersect_up, length(intersect(DEGs_up$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_up[i])])))
}

for (i in 1:length(colors_down)){
  
  jaccard_value_down <- c(jaccard_value_down,jaccard(DEGs_down$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_down[i])]))
  size_down <- c(size_down, length(wgcna_Modules$moduleColors[which(wgcna_Modules$moduleColors == colors_down[i])]))
  intersect_down <- c(intersect_down, length(intersect(DEGs_down$gene_name, wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == colors_down[i])])))
  
}


data_percentage_DGEs_up <- data.frame(modules = colors_up,
                                      jaccard_value = jaccard_value_up,
                                      size = size_up,
                                      seq = jaccard_value_up,
                                      intersect = intersect_up)
data_percentage_DGEs_down <- data.frame(modules = colors_down,
                                        jaccard_value = -jaccard_value_down,
                                        size = size_down,
                                        seq = jaccard_value_down)

data_percentage_DGEs <- rbind(data_percentage_DGEs_up, data_percentage_DGEs_down)

#| Plotting the results
ggplot(data_percentage_DGEs_up, aes(x =jaccard_value,y = reorder(modules,seq),  fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21, color ="black")+
  theme(axis.text.x = element_text(size = 12, color ="black"),
        axis.text.y = element_text(size = 12, color ="black"),
        text=element_text(size=12,  family="sans"), 
        plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors_up[order(colors_up)]) +
  #scale_color_manual(values = colors[order(colors)])+
  ylab("") +
  xlab("Jaccard index")+
  guides(fill ="none")+
  labs(size = "Module size") 
ggsave(paste(results.file, "Jaccard_index_WGCNA_DEGs_LimmaVoomm_up.pdf", sep =""),  height = 4.5, width = 4.5)

#| Plotting the results
ggplot(data_percentage_DGEs, aes(x =jaccard_value,y = reorder(modules,jaccard_value),  fill = modules)) + 
  geom_point(aes(size = size), alpha = 0.75, shape = 21, color ="black")+
  theme(axis.text.x = element_text(vjust = 0.99, hjust=1),text=element_text(size=12,  family="sans"), legend.key.size = unit(0.1, 'cm'), plot.title=element_text(size=14, hjust = 0.5, face ="bold")) +
  scale_fill_manual(values = colors[order(colors)]) +
  #scale_color_manual(values = colors[order(colors)])+
  ylab("Modules WGCNA (147 PCa samples)") +
  xlab("Jaccard index of DEGs KO6 vs WT6 \n in the WGCNA modules")+
  guides(fill ="none")+
  labs(size = "Module size") #+guides(color ="none")
ggsave(paste(results.file, "Jaccard_index_WGCNA_DEGs_LimmaVoomm.pdf", sep =""), heigh =4.5, width = 5 )
################################################################################


################################### Yellow #####################################
mod <- "yellow"
genes <- wgcna_Modules$GeneID[wgcna_Modules$moduleColors == mod]
genes_name <- wgcna_Modules$gene_name[wgcna_Modules$moduleColors == mod]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "mmusculus", ordered_query = TRUE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = wgcna_Modules$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
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
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 5)

#| Top with HP
dat1_filtered <- results[which(results$source == "HP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"HP_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4.4)


#| Transforming gene names of mouse to human
genes_name_human <- mouse_human$Gene.name[which(mouse_human$Mouse.gene.name %in% genes_name)]

#| UCell
signatures <- list()
signatures$mod_green <- genes_name_human

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
p + labs(title = "Module Yellow")
ggsave(file= paste(results.file,"UMAP_Chen_Module_Yellow.pdf", sep =""), width = 7, height = 5)

##############################################################################


################################# Tan ###################################
mod <- "black"
genes <- wgcna_Modules$GeneID[wgcna_Modules$moduleColors == mod]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "mmusculus", ordered_query = TRUE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = wgcna_Modules$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
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
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 5)

#| Top with KEGG
dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"HP_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4.4)
##############################################################################

################################# lightgreen ###################################
mod <- "black"
genes <- wgcna_Modules$GeneID[wgcna_Modules$moduleColors == mod]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "mmusculus", ordered_query = TRUE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = wgcna_Modules$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
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
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 5)

#| Top with KEGG
dat1_filtered <- results[which(results$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"HP_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 4.4)
##############################################################################


################### Yellow intersection with DEGs  #############################
color <- "yellow"
mod <- "yellow_DEGs"
genes <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == color)], DEGs$GeneID)
genes_name <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == color)], DEGs$gene_name)

total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "mmusculus", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = wgcna_Modules$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

writexl::write_xlsx(results, "Results/Tables/Supplementary_Table12_Module_yellow_DEGs_Mouse.xlsx")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
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
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 4.5, width = 7)

#| Top with KEGG
dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        axis.text.x = element_text(size =12, color = "black"),
        axis.text.y = element_text(size =12, color = "black"),
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"KEGG_Enrichment_module_",mod,".pdf",sep =""), height = 4.5, width = 6.5)

#| Transforming gene names of mouse to their human homologues
library(homologene)
package.version("homologene")

# Get human homologous genes
human_homologues <- homologene(genes_name, inTax = 10090, outTax = 9606)
human_homologues <- human_homologues$`9606`

#| UCell
signatures <- list()
signatures$mod_yellow_degs <- human_homologues

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

#| Module purple and DEGs label TRUE
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
               axis.text.x = element_text(size =12, color = "black"),
               axis.text.y = element_text(size =12, color = "black"),
               plot.title=element_blank())
p 
#ggsave(file= paste(results.file,"UMAP_Chen_Module_Yellow_DEGs.pdf", sep =""),  width = 6.4, height = 4.5)
ggsave(file= paste(results.file,"UMAP_Chen_Module_Yellow_DEGs_last.pdf", sep =""),  width = 6.4, height = 4.5)


#| Violin plot
v <- VlnPlot(loc.combined.sct, 
             features = signature.names[1], 
             sort=T,
             cols = colors$cols)
v <- v + theme(legend.position = 'none',
               axis.text.x = element_text(angle = 60, size =12),
               axis.text.y = element_text(size =12, color ="black"),
               plot.title = element_blank())
v + labs(title = "Module Yellow and DEGs") 
#ggsave(file= paste(results.file,"Vln_Chen_Module_Yellow_DEGs.pdf", sep =""), width = 8, height = 5)
ggsave(file= paste(results.file,"Vln_Chen_Module_Yellow_DEGs_last.pdf", sep =""), width = 8, height = 5)
##############################################################################


################### Turquoise intersection with DEGs  #############################
#| Size of module purple
length(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "turquoise")])
mod <- "Turquoise_DEGs"
genes <- intersect(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "turquoise")], DEGs$GeneID)
genes_name <- intersect(wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "turquoise")], DEGs$gene_name)

length(genes_name)
length(wgcna_Modules$GeneID[which(wgcna_Modules$moduleColors == "turquoise")])
length(DEGs$GeneID)

total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "mmusculus", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE,
                   user_threshold = 0.05, correction_method = "fdr",
                   domain_scope = "custom", custom_bg = wgcna_Modules$GeneID, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)
ggsave(paste(results.file,"ggplot_results_WGCNA_Module_",mod,".pdf", sep =""), height = 5, width = 7)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_",mod,".pdf", sep =""), height = 4.5, width = 5.3)

#| Top with REAC
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[2:11,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=14,  family="sans"),
        axis.text.x = element_text(size =12, color ="black"),
        axis.text.y = element_text(size =12, color ="black")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"REAC_Enrichment_module_",mod,".pdf",sep =""), height = 4.5, width = 8.2)

#| Top with KEGG
dat1_filtered <- results[which(results$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")
ggsave(paste(results.file,"KEGG_Enrichment_module_",mod,".pdf",sep =""), height = 3.5, width = 3.11)

#| Transforming gene names of mouse to their human homologues: Get human homologous genes
human_homologues <- homologene(genes_name, inTax = 10090, outTax = 9606)
human_homologues <- human_homologues$`9606`

#| UCell
signatures <- list()
signatures$mod_tuquoise_degs <- human_homologues

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

#| Module purple and DEGs label TRUE
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
               axis.text.x = element_text(size =12, color = "black"),
               axis.text.y = element_text(size =12, color = "black"),
               plot.title=element_blank())
p 
ggsave(file= paste(results.file,"UMAP_Chen_Module_Turquoise_DEGs.pdf", sep =""),  width = 6.4, height = 4.5)


#| Violin plot
v <- VlnPlot(loc.combined.sct, 
             features = signature.names[1], 
             sort=T,
             cols = colors$cols)
v <- v + theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 90))
v + labs(title = "Module Turquoise and DEGs") 
ggsave(file= paste(results.file,"Vln_Chen_Module_Turquoise_DEGs.pdf", sep =""), width = 8, height = 5)

##############################################################################



##################### % GREEN MOD IN ALL MODULES MICE ##########################

#| WGCNA mouse
wgcna_mouse <- wgcna_Modules
module_colors_mouse <- unique(wgcna_mouse$moduleColors)

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Loading data from Human samples PTEN loss vs intact
wgcna_human <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_intramodularinfo_last_try.txt", sep = "\t")
wgcna_human <- merge(wgcna_human, genome_GRCh39.94, by ="GeneID")
rownames(wgcna_human) <- wgcna_human$GeneID
module_colors_human <- unique(wgcna_human$moduleColors)

#| Extracting genes in module green
mod_green_human <- wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]
mod_green_human <- mouse_human$Mouse.gene.name[which(mouse_human$Gene.name %in% mod_green_human)]

#| Intersecting module green from HUMAN with all the modules found in MICE
intersect_D_list <- list()
intersect_D <- c()

for (i in 1:length(module_colors_mouse)){
  
  #| Saving genes intersection
  intersect_D_list[[i]] <- intersect(mod_green_human, wgcna_mouse$gene_name[which(wgcna_mouse$moduleColors == module_colors_mouse[i])])
  
  #| Saving length of genes intersection
  intersect_D <- c(intersect_D,jaccard(mod_green_human, wgcna_mouse$gene_name[which(wgcna_mouse$moduleColors == module_colors_mouse[i])]))
  
}

dataframe_intersection <- data.frame(moduleColors = module_colors_mouse,
                                     intersect_D = intersect_D)
ggplot(dataframe_intersection, aes(y= reorder(moduleColors,intersect_D), x = intersect_D)) +
  geom_bar(stat="identity", fill ="green", color ="black") +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("Jaccard index: genes intersected with module green human")+
  ylab("WGCNA Modules mouse KO6 vs WT6")
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_Jaccard_intersection_all_modules_MOUSE_module_green_HUMAN.pdf",  height = 6, width = 5.5)

################################################################################



################################################################################
#| CORRELATING HUMAN AND MOUSE INFORMATION COMPUTED WITH WGCNA
################################################################################

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Loading data from Human samples PTEN loss vs intact
wgcna_human <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Tables/3_geneInfo_counts_geneMM_geneS_moduleColors_intramodularinfo_last_try.txt", sep = "\t")
wgcna_human <- merge(wgcna_human, genome_GRCh39.94, by ="GeneID")
rownames(wgcna_human) <- wgcna_human$GeneID
module_colors_human <- unique(wgcna_human$moduleColors)

#| Modules green and purple from human PTEN loss vs intact
mod_green_human <- wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]
mod_purple_human <- wgcna_human$gene_name[which(wgcna_human$moduleColors == "purple")]

#| Module yellow from mouse analysis KO6 vs WT6. Enriched in stromal remodeling
mod_yellow_mouse <- wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "yellow")]
mod_yellow_mouse <- mouse_human$Gene.name[which(mouse_human$Mouse.gene.name %in% mod_yellow_mouse)]

#| Intersecting module yellow from mouse with all the modules found in human
intersect_D_list <- list()
intersect_D <- c()

for (i in 1:length(module_colors_human)){
  
  #| Saving genes intersection
  intersect_D_list[[i]] <- intersect(mod_yellow_mouse, wgcna_human$gene_name[which(wgcna_human$moduleColors == module_colors_human[i])])
  
  #| Saving length of genes intersection
  intersect_D <- c(intersect_D,jaccard(mod_yellow_mouse, wgcna_human$gene_name[which(wgcna_human$moduleColors == module_colors_human[i])]))
  
}

dataframe_intersection <- data.frame(moduleColors = module_colors_human,
                                     intersect_D = intersect_D)

ggplot(dataframe_intersection, aes(y= reorder(moduleColors,intersect_D), x = intersect_D)) +
  geom_bar(stat="identity", fill ="yellow", color ="black") +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("Jaccard index: genes intersected with module yellow mouse")+
  ylab("WGCNA Modules human PTEN protein loss vs intact")
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_Jaccard_intersection_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 6, width = 5.5)

#| Intersection
venn.diagram(
  x = list(mod_yellow_mouse,   
           wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]),
  category.names = c("mod_yellow_mouse" ,  
                     "mod green" ),
  filename = 'Results/Images/07_Enrichment_Analysis_of_Modules/Venn_diagramm_Mod_yellow_mouse_mod_green_human.png',
  output=T,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("yellow", "green"),
  cex = 1,
  fontface = "bold",
  cat.cex = 1,
  cat.pos= c(0,10),
  sub.fontfamily = "sans")

#| MODULE YELLOW FROM MOUSE AND MODULE GREEN FROM INTERSECTION: Which processes are enriched?
genes <- intersect_D_list[[4]]
total_gost <- gost(list("WGCNA Enrichment" = genes), 
                   organism = "hsapiens", 
                   ordered_query = TRUE, 
                   multi_query = FALSE, 
                   significant = TRUE, 
                   exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, 
                   evcodes = TRUE,
                   user_threshold = 0.05, 
                   correction_method = "fdr",
                   domain_scope = "custom", 
                   custom_bg = wgcna_human$gene_name, 
                   numeric_ns = "", 
                   sources = NULL, 
                   as_short_link = FALSE)
gostplot(total_gost, interactive = FALSE, capped =FALSE)

results <- total_gost$result[order(total_gost$result$p_value),]
results$`Term name` <- paste(results$term_name, "\n (N = ",results$term_size, ")",sep ="")

#| Top 10 most significant
ggplot(results[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.6, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+
  ylab("Term name")
ggsave(paste(results.file,"Enrichment_module_yellow_mouse_intersected_module_green_human.pdf", sep =""), height = 3.5, width = 5)

#| REACTOME
dat1_filtered <- results[which(results$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=9,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data") +
  ylab("Term name")

#| HOW ARE THE GENES IN MODULE YELLOW EXPRESS IN HUMAN SC DATA?
#| UMAP: SCORING SIGNATURES USING UCELL
signatures <- list()
signatures$mod_yellow_mouse <- mod_yellow_mouse
loc.combined.sct_filtered <- AddModuleScore_UCell(loc.combined.sct_filtered, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

#| DEGs UP
p <- FeaturePlot(loc.combined.sct_filtered, reduction = "umap", 
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
p 
ggsave(file= "Results/Images/07_Enrichment_Analysis_of_Modules/Chen_UMAP_LOC_DEGs_UP_PI3K-AKT-mTOR_Module_yellow_mouse_KO6_vs_WT6.pdf", width = 6.5, height = 5)




wgcna_Modules <- wgcna_Modules[order(wgcna_Modules$kWithin, decreasing = T),]

mod_yellow_mouse <- wgcna_Modules$gene_name[which(wgcna_Modules$moduleColors == "yellow")]
mod_yellow_mouse <- mouse_human$Gene.name[which(mouse_human$Mouse.gene.name %in% mod_yellow_mouse)]

intersect(ligands$ligand_gene_symbol,intersect(mod_yellow_mouse, wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]))
intersect(ligands$receptor_gene_symbol,intersect(mod_yellow_mouse, wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]))


g <- mod_yellow_mouse[which(mod_yellow_mouse %in% wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")])]

g <- mouse_human$Mouse.gene.name[which(mouse_human$Gene.name %in% g)]

wgcna_Modules$gene_name[which(wgcna_Modules$gene_name %in% g)]



wgcna_human <- wgcna_human[order(wgcna_human$kWithin, decreasing = T),]
g <- mod_yellow_mouse[which(mod_yellow_mouse %in% wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")])]

wgcna_human$gene_name[which(wgcna_human$gene_name %in% g)]

intersect(wgcna_human$gene_name[which(wgcna_human$gene_name %in% g)],DEGs$gene_name)

mod_green_human <- wgcna_human$gene_name[which(wgcna_human$moduleColors == "green")]

mod_yellow_mouse[which(mod_yellow_mouse=="TGFBR3")]
mod_green_human[which(mod_green_human=="TGFB2")]


#| Sample info
sample <- colnames(counts_data)
condition <- gsub("_[0-9]*","",sample)
month <- gsub("[A-Z]", "", condition) 
sample_info <- data.frame(sample = sample,
                          condition = condition,
                          month = month)
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info[order(sample_info$sample),]



#| Signature of the differences

counts.file <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt"
counts_data <- read.table(counts.file, sep ="\t")
counts_data <- counts_data[,order(colnames(counts_data))]
#| Normalization
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_data_normalized <- counts(dds, normalized=TRUE)
counts_data_normalized <- log(counts_data_normalized + 1, base =2)
counts_data_normalized <- as.data.frame(counts_data_normalized)

#| Load the genome for merging GeneID
genome_mouse<- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Data/geneInfo.txt", sep ="\t", header =F)
colnames(genome_mouse) <- c("GeneID", "gene_name", "gene_info")



counts_data_normalized$GeneID <- rownames(counts_data_normalized)
counts_data_normalized <- merge(counts_data_normalized, genome_mouse, by ="GeneID")
counts_data_normalized <- aggregate(counts_data_normalized, by =list(counts_data_normalized$gene_name), mean)
rownames(counts_data_normalized) <- counts_data_normalized$Group.1
counts_data_normalized <- counts_data_normalized[,sample_info$sample]


counts_data_normalized_6 <- counts_data_normalized[,sample_info$month ==6]
counts_zscore <- scale(t(as.data.frame(counts_data_normalized_6)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)

sample_info_6 <- sample_info[which(sample_info$month == "6"),]

#| Mean expression
genes <- intersect_D_list[[14]]
genes <- mouse_human$Mouse.gene.name[which(mouse_human$Gene.name %in% genes)]
genes <- genes[which(genes != "" )]
samples_media_expression <- colMeans(counts_zscore[genes,])
sample_info_6$genes <- samples_media_expression

cor <- wilcox.test(sample_info_6$genes[which(sample_info_6$condition == "WT6")], sample_info_6$genes[which(sample_info_6$condition == "KO6")])
cor$p.value

ggplot(sample_info_6, aes( x = condition , y = genes, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = condition),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="serif", size =4)+
  theme(text=element_text(size=12,  family="serif"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Mean expression of z scores \n from normalized log2 counts DEGs and Purple") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))+
  xlab("PTEN protein status")







#| And human? 

#| Counts data
counts_data <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt", header=T) 
counts_data <- counts_data[,rownames(datTraits)]
filter_counts <- rowSums(counts_data>5) >= 0.7*ncol(counts_data)
counts_data <- counts_data[filter_counts,]
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

#| Computing Z scores
counts_zscore <- scale(t(as.data.frame(counts_data_normalized)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)
counts_zscore <- as.data.frame(counts_zscore)

#| Green
genes <- intersect_D_list[[4]]
samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
datTraits$genes <- samples_media_expression

datTraits$PTEN_status[which(datTraits$PTEN_status == 1)] <- "Loss"
datTraits$PTEN_status[which(datTraits$PTEN_status == 0)] <- "Intact"
ggplot(datTraits, aes( x = PTEN_status, y = genes, fill = PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =4, hjust=-0.5)+
  theme(text=element_text(size=12,  family="sans"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Mean expression of z scores \n from normalized log2 counts") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))+
  xlab("PTEN protein status")
ggsave(file= "Results/Images/07_Enrichment_Analysis_of_Modules/Mean_expression_module_yellow_mouse_intersected_module_green_human.pdf", heigh=4, width= 4.5)

cor <- cor.test(datTraits$genes, datTraits$StromaScore)
p_value <- cor$p.value
r <- cor$estimate[[1]]
ggplot(datTraits, aes(y= genes, x =StromaScore)) +
  geom_point(color ="blue", size =3.5, alpha = 0.6) +
  geom_smooth(method = "lm", color ="black")+
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=10)) +
  xlab("Stromal Score")+
  ylab("Mean expression genes Module Green human \n intersected with Yellow Mice")+
  labs(size ="-log10(p_value)") +
  ggtitle(paste("r = ",round(r),"; p-value =",format(p_value, scientific=T,digits = 2), sep =""))
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/StromaScore_Mod_green_Human_Mod_Yellow_mouse.pdf", width = 4.5, height = 4.5)


#| yellow
genes <- intersect_D_list[[16]]
samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
datTraits$genes <- samples_media_expression

datTraits$PTEN_status[which(datTraits$PTEN_status == 1)] <- "Loss"
datTraits$PTEN_status[which(datTraits$PTEN_status == 0)] <- "Intact"
ggplot(datTraits, aes( x = PTEN_status, y = genes, fill = PTEN_status)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_status),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =4,hjust=-0.5)+
  theme(text=element_text(size=12,  family="sans"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Mean expression of z scores \n from normalized log2 counts") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))+
  xlab("PTEN protein status")
ggsave(file= "Results/Images/07_Enrichment_Analysis_of_Modules/Mean_expression_module_yellow_mouse_intersected_module_yellow_human.pdf", heigh=4, width= 4.5)


p_value <- c()
length_genes <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- intersect_D_list[[i]]
  if (length(genes) ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
    datTraits$genes <- samples_media_expression
    cor <- wilcox.test(datTraits$genes[which(datTraits$PTEN_status == "Loss")], datTraits$genes[which(datTraits$PTEN_status == "Intact")])
    
    p_value <- c(p_value, cor$p.value)
  }
  
}

length(mod_green_human)
module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]

data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                        p_value = p_value,
                        correlation=correlation)

ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = -log10(p_value))) +
  geom_bar(stat="identity", fill ="black", color ="black") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("-Log(p-value) PTEN protein loss vs intact")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_p_value_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 5, width = 5.5)


#| Correlation with the stromal
p_value <- c()
correlation <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- intersect_D_list[[i]]
  if (length(genes) ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
    datTraits$genes <- samples_media_expression
    cor <- cor.test(datTraits$genes , datTraits$stromal, method= "pearson")
    correlation <- c(correlation, cor$estimate[[1]])
    p_value <- c(p_value, cor$p.value)
  }
  
}

module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]
data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                                  p_value = p_value,
                                  correlation=correlation )

data_frame_p_values$value <- NA
data_frame_p_values$value[which(data_frame_p_values$correlation >= 0)] <- "Positive"
data_frame_p_values$value[which(data_frame_p_values$correlation < 0)] <- "Negative"
ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = correlation)) +
  geom_point(aes(size =-log10(p_value), color =value),stat="identity") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("Pearson correlation")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")+
  labs(size ="-log10(p_value)") +
  labs(color ="Direction")+
  scale_color_manual(values= c("red", "blue"))
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_STROMASCORE_CORRELATION_P_value_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 4.5, width = 5.5)



#| Correlation with the immune
p_value <- c()
correlation <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- intersect_D_list[[i]]
  if (length(genes) ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
    datTraits$genes <- samples_media_expression
    cor <- cor.test(datTraits$genes , datTraits$immune, method= "pearson")
    correlation <- c(correlation, cor$estimate[[1]])
    p_value <- c(p_value, cor$p.value)
  }
  
}

module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]
data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                                  p_value = p_value,
                                  correlation=correlation )

data_frame_p_values$value <- NA
data_frame_p_values$value[which(data_frame_p_values$correlation >= 0)] <- "Positive"
data_frame_p_values$value[which(data_frame_p_values$correlation < 0)] <- "Negative"
ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = correlation)) +
  geom_point(aes(size =-log10(p_value), color =value),stat="identity") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("Pearson correlation")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")+
  labs(size ="-log10(p_value)") +
  labs(color ="Direction")+
  scale_color_manual(values= c("red", "blue"))
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_IMMUNESCORE_CORRELATION_P_value_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 4.5, width = 5.5)


#| Correlation with the purity
p_value <- c()
correlation <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- intersect_D_list[[i]]
  if (length(genes) ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(counts_zscore[genes,], na.rm=TRUE)
    datTraits$genes <- samples_media_expression
    cor <- cor.test(datTraits$genes , datTraits$purity, method= "pearson")
    correlation <- c(correlation, cor$estimate[[1]])
    p_value <- c(p_value, cor$p.value)
  }
  
}

module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]
data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                                  p_value = p_value,
                                  correlation=correlation )

data_frame_p_values$value <- NA
data_frame_p_values$value[which(data_frame_p_values$correlation >= 0)] <- "Positive"
data_frame_p_values$value[which(data_frame_p_values$correlation < 0)] <- "Negative"
ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = correlation)) +
  geom_point(aes(size =-log10(p_value), color =value),stat="identity") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("Pearson correlation")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")+
  labs(size ="-log10(p_value)") +
  labs(color ="Direction")+
  scale_color_manual(values= c("red", "blue"))
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/WGCNA_PURITYSCORE_CORRELATION_P_value_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 4.5, width = 5.5)






#| Correlation with the stromal
p_value <- c()
correlation <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- counts_zscore[which(rownames(counts_zscore) %in%  intersect_D_list[[i]]),]
  
  if (dim(genes)[1] ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(genes, na.rm=TRUE)
    xcell_TCGA_2$genes <- samples_media_expression
    
    cor <- cor.test(xcell_TCGA_2$genes , xcell_TCGA_2$StromaScore, method= "pearson")
    correlation <- c(correlation, cor$estimate[[1]])
    p_value <- c(p_value, cor$p.value)
  }
  
}

length(xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="0")])
length(xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="-2")])


module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]

data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                                  p_value = p_value,
                                  correlation = correlation)
data_frame_p_values$value <- NA
data_frame_p_values$value[which(data_frame_p_values$correlation >= 0)] <- "Positive"
data_frame_p_values$value[which(data_frame_p_values$correlation < 0)] <- "Negative"
ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = -log10(p_value))) +
  geom_point(aes(size =correlation, color =value),stat="identity") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("-Log(p-value) PTEN protein loss vs intact")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")+
  labs(size ="r") +
  labs(color ="Direction")+
  scale_color_manual(values= c("red", "blue"))









#| Genes intersected with purple
intersect_D_list[[16]]

#| Genes intersected with green
intersect_D_list[[4]]


datTraits$CLIC4 <- as.numeric(counts_zscore["CLIC4",])
datTraits$PTEN <- as.numeric(counts_zscore["PTEN",])
ggplot(datTraits, aes(x = PTEN, y = CLIC4)) +
  geom_point()


datTraits$PLXDC2 <- as.numeric(counts_zscore["PLXDC2",])
datTraits$PTEN <- as.numeric(counts_zscore["PTEN",])
ggplot(datTraits, aes(x = CLIC4, y = StromaScore )) +
  geom_point()


#| TCGA?


xcell_TCGA <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/xCell/Table_results_TCGA.txt", sep ="\t")
data_phenotype_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostatephenotype")
data_expression_tcga <- dbGetQuery(conn = con, statement = "SELECT * from tcgaprostateexpression")
data_expression_tcga$Gene <- gsub("\\|", "", data_expression_tcga$Gene)

#| Finding low counts
low_counts <- data_expression_tcga$Gene[which(rowSums(data_expression_tcga[,-c(1)] > 0) > 0.3*dim(data_expression_tcga[,-c(1)])[2] )]

#| Filtering low counts counts
data_expression_tcga <- data_expression_tcga[ which(data_expression_tcga$Gene %in% low_counts),]
rownames(data_expression_tcga) <- data_expression_tcga$Gene

pten_cna <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_TCGA/TCGA PTEN info/cna.txt", sep = "\t", header =T)
colnames(pten_cna) <- c("STUDY_ID","PATIENT_ID", "PTEN_cna")
pten_cna$PATIENT_ID <- sub("-01","",pten_cna$PATIENT_ID )
pten_cna <- pten_cna[which(pten_cna$PATIENT_ID %in% rownames(xcell_TCGA)),]
pten_cna <- pten_cna[order(pten_cna$PATIENT_ID),]
xcell_TCGA$PTEN_cna <- pten_cna$PTEN_cna
xcell_TCGA <- xcell_TCGA[order(rownames(xcell_TCGA)),]

rownames(xcell_TCGA) == pten_cna$PATIENT_ID
data_expression_tcga <- data_expression_tcga[, rownames(xcell_TCGA)]

#| Computing z-scores
counts_zscore <- scale(t(as.data.frame(data_expression_tcga)))
counts_zscore <- as.data.frame(counts_zscore)
counts_zscore <- t(counts_zscore)
counts_zscore <- as.data.frame(counts_zscore)


xcell_TCGA_2 <- xcell_TCGA[which((xcell_TCGA$PTEN_cna=="0") | (xcell_TCGA$PTEN_cna=="-2") ),]
#| Correlation with the stromal
p_value <- c()
correlation <- c()
color <- c()
for (i in 1:length(module_colors_human)){
  
  genes <- counts_zscore[which(rownames(counts_zscore) %in%  intersect_D_list[[i]]),]

  if (dim(genes)[1] ==0 ){
    color <- c(color, module_colors_human[i])
  }else{
    samples_media_expression <- colMeans(genes, na.rm=TRUE)
    xcell_TCGA_2$genes <- samples_media_expression
    
    cor <- cor.test(xcell_TCGA_2$genes , xcell_TCGA_2$StromaScore, method= "pearson")
    correlation <- c(correlation, cor$estimate[[1]])
    p_value <- c(p_value, cor$p.value)
  }
  
}

length(xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="0")])
length(xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="-2")])


module_colors_human_2 <- module_colors_human[which(!(module_colors_human %in% color))]

data_frame_p_values <- data.frame(module_colors_human = module_colors_human_2,
                                  p_value = p_value,
                                  correlation = correlation)
data_frame_p_values$value <- NA
data_frame_p_values$value[which(data_frame_p_values$correlation >= 0)] <- "Positive"
data_frame_p_values$value[which(data_frame_p_values$correlation < 0)] <- "Negative"
ggplot(data_frame_p_values, aes(y= reorder(module_colors_human,-log10(p_value)), x = -log10(p_value))) +
  geom_point(aes(size =correlation, color =value),stat="identity") +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.5, 'cm'), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  xlab("-Log(p-value) PTEN protein loss vs intact")+
  ylab("WGCNA Modules HUMAN PTEN protein loss vs intact")+
  labs(size ="r") +
  labs(color ="Direction")+
  scale_color_manual(values= c("red", "blue"))
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/TCGA_STROMASCORE_CORRELATION_P_value_all_modules_HUMAN_module_yellow_MOUSE.pdf",  height = 4.5, width = 5.5)


xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="-2")] <- "Loss"
xcell_TCGA_2$PTEN_cna[which(xcell_TCGA_2$PTEN_cna =="0")] <- "Intact"
ggplot(xcell_TCGA_2, aes( x = PTEN_cna, y = StromaScore, fill = PTEN_cna)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_cna),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =4,hjust=-0.5)+
  theme(text=element_text(size=12,  family="sans"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Stroma Score") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF"))+
  xlab("PTEN genomic status")
ggsave("Results/Images/07_Enrichment_Analysis_of_Modules/TCGA_PTEN_genomic_loss_intact_StromaScore_xcell.pdf",  height = 5, width = 5.5)


#| Green
genes <- counts_zscore[which(rownames(counts_zscore) %in%  intersect_D_list[[4]]),]
samples_media_expression <- colMeans(genes, na.rm=TRUE)
xcell_TCGA$genes <- samples_media_expression


ggplot(data_phenotype_tcga[which( (data_phenotype_tcga$PTEN_cna == "0") | (data_phenotype_tcga$PTEN_cna== "-2") |(data_phenotype_tcga$PTEN_cna== "-1") ),], aes( x = PTEN_cna, y = genes, fill = PTEN_cna)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_cna),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =4, hjust=-0.5)+
  theme(text=element_text(size=12,  family="sans"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Mean expression of z scores \n from normalized log2 counts") +
  scale_fill_manual(values = rainbow(3))+
  xlab("PTEN protein status")
ggsave(file= "Results/Images/07_Enrichment_Analysis_of_Modules/Mean_expression_module_yellow_mouse_intersected_module_green_human_TCGA.pdf", heigh=4, width= 4.5)

#| purple
genes <- counts_zscore[which(rownames(counts_zscore) %in%  intersect_D_list[[16]]),]
samples_media_expression <- colMeans(genes, na.rm=TRUE)
data_phenotype_tcga$genes <- samples_media_expression

ggplot(data_phenotype_tcga[which( (data_phenotype_tcga$PTEN_cna == "0") | (data_phenotype_tcga$PTEN_cna== "-2")),], aes( x = PTEN_cna, y = genes, fill = PTEN_cna)) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
  geom_point(aes(fill = PTEN_cna),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),pch = 21, size=2, alpha=0.6)+
  stat_compare_means(label ="p.format",family="sans", size =4, hjust=-0.5)+
  theme(text=element_text(size=12,  family="sans"), legend.position = "none",
        legend.key.size = unit(0.8, 'cm'),
        plot.title=element_text(size=16, hjust = -1, face ="bold")) +
  ylab("Mean expression of z scores \n from normalized log2 counts") +
  scale_fill_manual(values = c("#504667", "#EFD5B2"))+
  xlab("PTEN protein status")
ggsave(file= "Results/Images/07_Enrichment_Analysis_of_Modules/Mean_expression_module_yellow_mouse_intersected_module_Purple_human_TCGA.pdf", heigh=4, width= 4.5)


#| Plot here the most significant?
