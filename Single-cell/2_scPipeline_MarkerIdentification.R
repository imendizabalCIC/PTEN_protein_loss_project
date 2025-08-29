
####################################################################
#| SINGLE CELL MARKER ANNOTATION
####################################################################
#| Date: 12/12/2024
#| Author: Saioa Garcia Longarte
#|
#| Description:
#| This script performs manual cell type annotation of Chen's scRNA-seq prostate 
#| tumor dataset after integration. It identifies cluster-specific marker genes, 
#| compares them with canonical lineage markers, and assigns biological labels to 
#| clusters for downstream analysis.  
#|
#| Workflow:
#|   1) Load integrated Seurat object and project metadata from the log file.  
#|   2) Run FindAllMarkers to identify cluster-enriched genes.  
#|   3) Save top markers per cluster (top 15 and top 100).  
#|   4) Visualize canonical cell type markers via:  
#|        * DotPlots (RNA expression across clusters)  
#|        * UMAP FeaturePlots (marker distribution in 2D space)  
#|        * Violin plots (marker expression by cluster)  
#|   5) Assign cluster identities based on marker expression (e.g. luminal, T cells, 
#|      endothelial, macrophages, fibroblasts, mast cells, cycling cells).  
#|   6) Save annotated Seurat object and generate annotated cluster visualizations 
#|      (UMAP, PCA, DotPlots).  
#|
#| Outputs:
#|   - Text files with top 15 and top 100 marker genes per cluster.  
#|   - PDF/PNG plots of canonical markers (DotPlot, FeaturePlot, ViolinPlot).  
#|   - Annotated UMAP and PCA cluster plots.  
#|   - Annotated Seurat object saved as `Annotation_Manual_<project>.rds`.  
#|
#| Process PT data from Chen

## General pipeline for single-cell data analysis
## Following:
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Guided%20Clustering%20Tutorial%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Introduction%20to%20scRNA-seq%20integration%20(May2022).R 
## -  Single-cell RNA-seq data analysis workshop: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html 

## Definitions of the functions from: https://cloud.r-project.org/web/packages/Seurat/Seurat.pdf (Version 4.1.1 || Date 2022-05-01)
####################################################################


####################################################################
## Libraries
####################################################################

suppressPackageStartupMessages(require("SingleR"))
suppressPackageStartupMessages(require("scCATCH"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("future"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("dplyr"))
suppressPackageStartupMessages(require("readxl"))

###################################################################
## Variable prep
###################################################################

setwd("C:/Users/sgarcia/CIC bioGUNE/Arkaitz group - Documentos/Individual folders/Saioa Garcia/Requests/Ivana/Request(2024_11_18) - Chen PT")
dir.create("4_MarkerIdetification_Figures")

# The variables are stored in the log file from the 2_scPipeline_Integration
log_file <- read.delim("W:/sgarcia/sc_Chen_all/2_scPipeline_Integration_13052022_LOG.log", header = TRUE)

project_name <- log_file$Project_name
samples_v <- log_file$Samples
organism <- log_file$Organism
approach <- log_file$Approach
integration_reduction <- log_file$Integration_reduction
pc_number_integration <- log_file$PC_number_integration
resolution <- log_file$Resolution_find_clusters

if (approach == "Norm_Feature_Scale"){
  DefaultAssay(data_integrated) <- "RNA"
} else if (approach == "SCT") {
  DefaultAssay(data_integrated) <- "SCT"
  data_integrated <- PrepSCTFindMarkers(data_integrated, assay = "SCT")
  saveRDS(data_integrated, paste("1_SeuratObjects/Data_integrated_PrepSCT_", project_name,  ".rds", sep = ""))
}

column_name <- colnames(data_integrated@meta.data)[grep("integrated_snn", colnames(data_integrated@meta.data))]

for (i in column_name){
  data_integrated[[i]] <- NULL
}

###################################################################
## 1. Mannual annotation
###################################################################
data_markers <- FindAllMarkers(data_integrated, slot = "data", assay = approach, min.pct = 0.1, logfc.threshold = 0.25)
dim(data_markers)
table(data_markers$cluster)
top15_markers <- as.data.frame(data_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC))
write.table(top15_markers, file = paste("4_MarkerIdetification_Figures/Annotation_Manual_top15genesFindAllMarkers_", project_name,  ".txt", sep = ""), row.names = FALSE)
top100_markers <- as.data.frame(data_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
write.table(top100_markers, file = paste("4_MarkerIdetification_Figures/Annotation_Manual_top100genesFindAllMarkers_", project_name,  ".txt", sep = ""), row.names = FALSE)

to_plot <- c("CD3D", "CD3E", "CD3G", "PTPRC", # T cells - PTPRC, CD3D, CD3E, and CD3G
             "CD79A", "CD79B", "IGKC", "MS4A1", # B cells - CD79A, IGKC, IGLC2, and IGLC3
             "CD14", "CD68", "CSF1R", "FCGR3A", "LYZ", # macrophages - CD68, CD163, FCGR3A, LYZ,and CSF1R 
             "KIT", "MS4A2", "TPSAB1", "TPSB2",# mast cells - ENPP3, KIT, SLC18A2, and MS4A2
             "COL1A1", "COL1A2", "COL3A1", "DCN", # fibros - aSMA, FN1, and FAP, "LUM", 
             "RGS5", "ACTA2", 
             "TAGLN", "TPM2", # MYOFIBROBLASTS (contractlity)
             "PDGFRB",  # MURAL CELLS
             "CSPG4", # PERICYTES
             "CDH5", "ENG", "PECAM1", "VWF", # endo - PECAM1, ENG, CDH5, and VWF
             "AR", "EPCAM", "KRT5", "KRT8", "KRT14", # epi - EPCAM, AR, KRT5, KRT14, KRT8, and KRT18
             "MUC15", "CYP4B1", "SCGB3A1", "BPIFB1", # CLUB CELLS
             "BIRC5", "CENPF" # cycling
)

DefaultAssay(data_integrated) <- "RNA"

pdf(paste("4_MarkerIdetification_Figures/DotPlot_RNA_", project_name, ".pdf",sep=''), width=15, height=10)
DotPlot(data_integrated, features = to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

# In the features line write the markers for the Cell Type
for (feature in to_plot){
  FeaturePlot(data_integrated,
              reduction = "umap",
              features = feature,
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE,
              keep.scale = "all")
  ggsave(file= paste("4_MarkerIdetification_Figures/Annotation_Manual_FeaturePlot_", feature , "_", project_name, ".png", sep = ""), width = 8, height = 6)

  VlnPlot(data_integrated, features = feature, pt.size = 0) # pt.size = 0 without points
  ggsave(file= paste("4_MarkerIdetification_Figures/Annotation_Manual_VlnPlot_", feature , "_", project_name, ".png", sep = ""), width = 12, height = 6)
}

print("WARNIING: Make the changes in the next variable for the cluster names.")

new.cluster.ids <- c("0" = "0-Luminal",
                     "1" = "1-Luminal", 
                     "2" = "2-Luminal",
                     "3" = "3-Luminal",
                     "4" = "4-Luminal",
                     "5" = "5-Luminal",
                     "6" = "6-T cell",
                     "7" = "7-Endothelial",
                     "8" = "8-Intermediate",
                     "9" = "9-Macrophage",
                     "10" = "10-Myofibroblasts",
                     "11" = "11-Luminal",
                     "12" = "12-Mast cell",
                     "13" = "13-Luminal",
                     "14" = "14-Luminal",
                     "15" = "15-Luminal",
                     "16" = "16-Endothelial",
                     "17" = "17-Luminal cycling")

names(new.cluster.ids) <- levels(data_integrated)
data_annotated <- RenameIdents(data_integrated, new.cluster.ids)
data_annotated$annotation <- data_annotated@active.ident

DimPlot(data_annotated, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.8) 
ggsave(file= paste("4_MarkerIdetification_Figures/Manual annotation_", project_name, ".png", sep = ""), width = 15, height = 10)

DimPlot(data_annotated, 
        reduction = "pca", 
        label = TRUE, 
        pt.size = 0.8) 
ggsave(file= paste("4_MarkerIdetification_Figures/DimPlot_pca_", project_name, ".png", sep = ""), width = 15, height = 15)

pdf(paste("4_MarkerIdetification_Figures/DotPlot_Annotated_", project_name, ".pdf",sep=''), width=15, height=10)
DotPlot(data_annotated, features = to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

saveRDS(data_annotated, paste("1_SeuratObjects/Annotation_Manual_", project_name, ".rds", sep = ""))
