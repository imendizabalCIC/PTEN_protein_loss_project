#| Last change: 12/12/2024
#| Ivana Rondon-Lorefice


################################################################################
#| UCell                 
################################################################################

#| Single cell data from Chen: https://www.nature.com/articles/s41556-020-00613-6
#| This contains 13 samples (8 localized)
#| After filtering localized samples, I have integrated the data using these tutorial: 
#|      https://satijalab.org/seurat/articles/sctransform_vignette.html 
#|      https://satijalab.org/seurat/articles/sctransform_v2_vignette
#|      https://satijalab.org/seurat/articles/integration_introduction.html
#| In the integration, I have applied SCTransform independently. Then followed 
#| the functions pre-defined for SC-data integration. However, take into account
#| the comments from satija about the RNA Assay vs SCTAssay for visualization: 
#|      https://github.com/satijalab/seurat/issues/4082

#| When creating the annotations, I have to use the RNA assay instead of the integrated
#| assay, because in this one, there are a number of feature selected for optimal
#| visualization

#| UCell: Robust and scalable single-cell gene signature scoring

#| To install: https://github.com/carmonalab/UCell

#| Use of UCell and downloading a single-cell data from the GEO repository. 
#| For Chen: GSE141445 

#| https://carmonalab.github.io/UCell_demo/UCell_matrix_vignette.html
#| https://carmonalab.github.io/UCell_demo/UCell_Seurat_vignette.html

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
suppressPackageStartupMessages(require("UCell"))
suppressPackageStartupMessages(require("viridis"))

#| Typically, the set.seed() would be placed at the beginning of your script. In this way the selected random number would be applied to any function that uses pseudorandom numbers in its algorithm.
set.seed(123456)

#| To check that all packages have been loaded
sessionInfo()

#| Set working directory
workingDir <- "X:/irondon/Project_AC-45_RNAseq-FFPE/SingleCell/"
setwd(workingDir)
################################################################################


################################################################################
#| DATA HUMAN
################################################################################

#| Opening the processed data
loc.combined.sct <- readRDS("Annotation_Manual_Chen_PT.rds")

#| DEGs in Module purple
mod_purple <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Images/7_Enrichment_Analysis_of_Modules/Modules_signature/DEGs_purple.txt", sep ="\t")
mod_purple <- mod_purple$genes_name

#| DEGs in Module green
mod_green <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_intact/Results/Images/7_Enrichment_Analysis_of_Modules/Modules_signature/DEGs_green.txt", sep ="\t")
mod_green <- mod_green$genes_name

#| PI3K-AKT-mTOR signature
pi3k_akt_mtor <- read.table("X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt", sep =",", header= T)
pi3k_akt_mtor <- names(pi3k_akt_mtor)

#| Saving signatures to test in a list
signatures <- list()
signatures$mod_purple <- mod_purple
signatures$mod_green <- mod_green
signatures$pi3k_akt_mtor <- pi3k_akt_mtor

################################################################################


################################################################################
#| HUMAN PROJECTION OF MODULES IN SC-DATA USING UCELL
################################################################################

#| Default assay to evaluate signatures with UCell
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
p + labs(title = "DEGs in module Purple")
ggsave(file= "Results/Images/UMAPPlots/Chen_FeaturePlot_LOC_integrated_modules_purple_DEGs_Label_Human.pdf", width = 7, height = 5)

#| Module green and DEGs label TRUE
p <- FeaturePlot(loc.combined.sct, reduction = "umap", 
                 features = signature.names[2], 
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
p + labs(title = "DEGs in module Green")
ggsave(file= "Results/Images/UMAPPlots/Chen_FeaturePlot_LOC_integrated_modules_green_DEGs_Label_Human.pdf", width = 7, height = 5)


#| PI3K-AKT-mTOR signature per cluster
VlnPlot(loc.combined.sct, 
        features = signature.names[3], 
        sort=T,
        cols = rainbow(21))
ggsave(file= "Results/Images/VlnPlots/Chen_VlnPlot_LOC_integrated_PI3K-AKT-mTOR_Human.pdf", width = 10, height = 5)

#| PI3K-AKT-mTOR expression per samples
VlnPlot(loc.combined.sct, 
        features = signature.names[3],
        group.by = "orig.ident",
        sort = T,
        cols = rainbow(8))
ggsave(file= "Results/Images/VlnPlots/Chen_VlnPlot_LOC_integrated_PI3K-AKT-mTOR_orig.ident_Human.pdf", width = 10, height = 5)
################################################################################



################################################################################
#| DATA MOUSE
################################################################################

#| Genome mouse
genome_mouse <- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Data/geneInfo.txt", sep ="\t", header =F)
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

#| Limma voom
DEGs <- DEGs_genes[which(( (DEGs_genes$log2FC_limma_voom > log2(FC)) | (DEGs_genes$log2FC_limma_voom <  (-log2(FC)))) & (DEGs_genes$padj_limma_voom < FDR) & (!is.na(DEGs_genes$padj_limma_voom))),]

#| Genes in module Yellow intersected with DEGs
mod_yellow <-intersect(wgcna_Modules$gene_name,DEGs$gene_name)

#| Transforming gene names of mouse to human
genes_name_human_mod_yellow <- mouse_human$Gene.name[which(mouse_human$Mouse.gene.name %in% mod_yellow)]

#| Saving signatures to test in a list
signatures <- list()
signatures$mod_yellow <- genes_name_human_mod_yellow

################################################################################


################################################################################
#| MOUSE PROJECTION OF MODULES IN SC-DATA USING UCELL
################################################################################

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
ggsave(file= "Results/Images/UMAPPlots/Chen_FeaturePlot_LOC_integrated_modules_yellow_DEGs_Label_Mouse.pdf", width = 7, height = 5)

################################################################################