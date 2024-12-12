#| Last change: 12/12/2024
#| Saioa Garcia-Longarte

setwd("C:/Users/sgarcia/CIC bioGUNE/Arkaitz group - Documentos/Individual folders/Saioa Garcia/Requests/Ivana/Request(2024_11_18) - Chen PT")

dir.create(paste("3_Integration_Figures", sep =""))
dir.create(paste("3_Integration_Figures/Integration", sep =""))
dir.create(paste("3_Integration_Figures/No_integration", sep =""))
dir.create(paste("3_Integration_Figures/SCT", sep =""))

# Process PT data from Chen


####################################################################
## Libraries
####################################################################

suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("dplyr"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("purrr"))
suppressPackageStartupMessages(require("future"))
suppressPackageStartupMessages(require("SingleR"))
suppressPackageStartupMessages(require("reshape2"))
source("utils/custom_seurat_functions.R")


###################################################################
## Variable prep
###################################################################
# //////////////////////////////////////////////////////////////////////////////
# Integration type choose CCA or PCA for the integration (options: "rPCA" || "CCA")
print("Remember to choose between rPCA or CCA for the integration.")
integration_reduction <- "CCA"
# //////////////////////////////////////////////////////////////////////////////

project_name <- "Chen_PT"
samples_v <- c("SC159", "SC174", "SC162", "SC176", "SC154", "SC155", "SC156", "SC153")

organism <- "human"

# Approach (options: "Norm_Feature_Scale" || "SCT")
approach <- "SCT"

print(paste("According to your log file, you are currently usign: ", approach, sep = ""))
print("If your approach is Norm_Feature_Scale, please choose rPCA or CCA for integration. 
  CCA it's more accurate when cell types are conserved but could lead to overcorrection when a large proportion of cells are non overlaping across datasets and it is more computationally expensive. 
  (https://satijalab.org/seurat/articles/integration_rpca.html")

samples_v <- unlist(strsplit(samples_v,","))
samples_v <- gsub(" ", "", samples_v)

for (sample in samples_v){
  # We are going to save each of the objects separately 
  RDS <- readRDS(paste("W:/sgarcia/sc_Chen_all/1_SeuratObjects/QC_", sample, "_Chen_all.rds", sep = ""))
  column_name <- colnames(RDS@meta.data)[grep("pANN", colnames(RDS@meta.data))]
  RDS[[column_name]] <- NULL
  assign(sample, RDS)
}

if (length(samples_v) > 1){
  data_combined <- merge(x = get(samples_v[1]), y = sapply(samples_v[-1], get), add.cell.ids = samples_v, project = project_name)
}

###################################################################
## 1. Pre-processing workflow
###################################################################
# There are 2 approaches: 
# 1. Normalization + FindVariableFeatures + ScaleData workflow
# 2. SCT workflow

# Iterating over samples in a dataset
# Since we have several samples in our dataset (from different conditions), we want to keep them as separate objects and transform them as that is what is required for integration.
data_combined <- SetIdent(data_combined, value = "orig.ident")

DefaultAssay(data_combined) <- "RNA"
split_data <- SplitObject(data_combined, split.by = "orig.ident")

# both approaches are going to be performed
##################################################
## Norm_Feature_Scale
##################################################
# https://satijalab.org/seurat/articles/integration_introduction.html

# normalize and identify variable features for each dataset independently
split_data <- lapply(X = split_data, FUN = function(x) {
  ##################################################
  ## Normalizing the data
  ##################################################
  # By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in data_combined[["RNA"]]@data.
  # normalization.method: Method for normalization. There are three: LogNormalize, CLR, RC. Default: normalization.method = "LogNormalize"
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ##################################################
  ## Identification of highly variable features (feature selection)
  ##################################################
  # We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
  
  # Only works for RNA assay
  gene_info_distribution <- summary(Matrix::colSums(x@assays$RNA@counts[,]>0))
  hvg_number <- round(gene_info_distribution[4]+100)
  print(paste("You will calculate FindVariableFeatures with: ", hvg_number, " genes", sep = ""))
  
  # FindVariableFeatures: identifies features that are outliers on a ’mean variability plot’.
  # Selection.method: How to choose top variable features. There are three: vst, mean.var.plot, dispersion. Default: vst
  # nfeatures = 2000 by default
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg_number)
})

##################################################
## SCTransform - https://satijalab.org/seurat/articles/sctransform_vignette.html
##################################################
for (l in 1:length(split_data)){
  DefaultAssay(split_data[[l]]) <- "SCT"
}

# Now we will use a ‘for loop’ to run the SCTransform() on each sample, and regress out mitochondrial expression by specifying in the vars.to.regress argument of the SCTransform() function.
# Before we run this for loop, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:

# SCTransform: This function calls sctransform::vst. The sctransform package is available at https://github.com/ChristophH/sctransform. Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData work-flow. Results are saved in a new assay (named SCT by default) with counts being (corrected) counts, data being log1p(counts), scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay.
# Now, we run the following loop to perform the sctransform on all samples. This may take some time (~10 minutes):
split_data <- lapply(X = split_data, FUN = function(x) {
  x@meta.data$orig.ident[1]
  x <- SCTransform(x, vars.to.regress = c("percent_mt"), vst.flavor = "v2")
  a <- RunPCA(x, npcs = 30, assay = "SCT") %>% 
    RunUMAP(reduction = "pca", dims = 1:30, assay = "SCT") %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.7)
  DimPlot(a, label = T, repel = T) + ggtitle("Unsupervised clustering")
  ggsave(paste("3_Integration_Figures/SCT/Unsupervised clustering_", x@meta.data$orig.ident[1], "_", project_name, ".pdf", sep = ""))
  return(x)
})

##################################################
## UMAP plot before integration
##################################################
# //////////////////////////////////////////////////////////////////////////////
# Write as many as samples in y
non_integrated_data <- merge(x = split_data[[1]], y = c(split_data[[2]], 
                                                        split_data[[3]], 
                                                        split_data[[4]], 
                                                        split_data[[5]],
                                                        split_data[[6]], 
                                                        split_data[[7]], 
                                                        split_data[[8]]),
                             add.cell.ids = samples_v, project = project_name)
# //////////////////////////////////////////////////////////////////////////////

DefaultAssay(non_integrated_data) <- "RNA"

# Only works for RNA assay
gene_info_distribution <- summary(Matrix::colSums(non_integrated_data@assays$RNA@counts[,]>0))
hvg_number <- round(gene_info_distribution[4]+100)
print(paste("You will calculate FindVariableFeatures with: ", hvg_number, " genes", sep = ""))
non_integrated_data <- FindVariableFeatures(non_integrated_data, selection.method = "vst", nfeatures = hvg_number) 
non_integrated_data <- ScaleData(non_integrated_data)                                                              
non_integrated_data <- RunPCA(object = non_integrated_data, assay = "RNA")

# Plot PCA
PCAPlot(non_integrated_data, group.by = "orig.ident", split.by = "orig.ident")  
ggsave(paste("3_Integration_Figures/No_integration/PCAPlot_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

# From https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- non_integrated_data[["pca"]]@stdev / sum(non_integrated_data[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
max_dim_no_integration <- min(co1, co2)

# Plot the elbow plot
ElbowPlot(object = non_integrated_data, ndims = 60)
ggsave(paste("3_Integration_Figures/No_integration/ElbowPlot_", project_name, ".png", sep = ""), width = 10, height = 8, bg="white")

# Run UMAP
non_integrated_data <- RunUMAP(non_integrated_data, dims = 1:max_dim_no_integration, reduction = "pca", min.dist = 0.6, assay = approach)

# VizDimLoadings: Visualize top genes associated with reduction components
VizDimLoadings(non_integrated_data, dims = 1:2, reduction = "pca")
ggsave(paste("3_Integration_Figures/No_integration/VizDimLoadings_top_genes_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Plot UMAP
# DimPlot: Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a cell and it’s positioned based on the cell embeddings determined by the reduction technique. 
DimPlot(non_integrated_data, group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/UMAP_preIntegration_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Plot UMAP split by sample (Side-by-side comparison of clusters)
# Sometimes it’s easier to see whether all of the cells align well if we split the plotting between conditions, which we can do by adding the split.by argument to the DimPlot() function:
DimPlot(non_integrated_data, split.by = "orig.ident", group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/No_integration/No_Integration_UMAP_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

saveRDS(non_integrated_data, file = paste("1_SeuratObjects/", "Data_NON_integrated_", project_name, ".rds", sep = ""))
saveRDS(split_data, file = paste("1_SeuratObjects/", "Data_split_data_", project_name, ".rds", sep = ""))

# //////////////////////////////////////////////////////////////////////////////
# In case of integration, please execute script: 
"0_Scripts/2.3_scPipeline_Integration_2022-09-20.R"
# //////////////////////////////////////////////////////////////////////////////

# Write a Log output
Date_time <- Sys.time()
log_data <- as.data.frame(Date_time, row.names = NULL)
log_data$Project_name <- project_name
log_data$Samples <- do.call(paste, c(as.list(samples_v), sep = ", "))
log_data$Organism <- organism
log_data$Approach <- approach
log_data$Integration_reduction <- integration_reduction
# log_data$PC_number_doublet <- pc_number_doublet
log_data$PC_number_integration <- max_dim
# log_data$PC_number_no_integration <- max_dim_no_integration
log_data$Resolution_find_clusters <- resolution_find_clusters
# log_data$Resolution_find_clusters_NON_integrated <- no_integration_resolution_find_clusters
# log_data$BigData_path <- BigData_path

write.table(log_data, file = "2_scPipeline_Integration_13052022_LOG.log", sep = "\t", row.names = FALSE)
