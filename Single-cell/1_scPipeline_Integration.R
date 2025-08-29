###################################################################
#| SINGLE CELL DATA INTEGRATION
###################################################################
#| Date: 12/12/2024
#| Author: Saioa Garcia Longarte
#|
#| Description:
#| This script integrates scRNA-seq datasets using Seurat, harmonizing batch effects 
#| while preserving biological variation. It supports both standard (LogNormalize + scaling) 
#| and SCTransform-based workflows.
#|
#| Workflow:
#|   1) Select top 3000 variable genes across datasets for integration.
#|   2) Build integration anchors using CCA or RPCA (for Norm_Feature_Scale) 
#|      or SCT-based anchors (for SCT approach).
#|   3) Integrate datasets into a unified Seurat object.
#|   4) Combine doublet/singlet annotations (DF) and visualize distributions
#|      across samples, phenotypes, and clusters.
#|   5) Perform linear dimensionality reduction (PCA, UMAP) on integrated data.
#|   6) Explore elbow plot to guide number of PCs for clustering.
#|   7) Cluster cells at multiple resolutions and assign cluster identities.
#|   8) Visualize UMAP embeddings split by sample, phenotype, DF, and metadata.
#|   9) Assess potential technical variation (UMIs, features, ribo/mito %).
#|  10) Explore PCs driving clusters, cell cycle scores, and DF proportions.
#|
#| Outputs:
#|   - Barplots of DF distributions (by sample, phenotype, clusters).
#|   - PCA and UMAP plots (colored by sample, pheno, DF, clusters).
#|   - Elbow plot to select PCs.
#|   - Cluster distribution tables and visualizations.
#|   - Saved integrated Seurat object for downstream analyses.
###################################################################

###################################################################
## 2. Integration
###################################################################
# First, we need to specify that we want to use all of the 3000 most variable genes identified by SCTransform for the integration. By default, this function only selects the top 2000 genes.
# SelectIntegrationFeatures: Choose the features to use when integrating multiple datasets. This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. It returns the top scoring features by this ranking. 
# https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
integ_features <- SelectIntegrationFeatures(object.list = split_data, nfeatures = 3000) 

if (approach == "Norm_Feature_Scale"){
  for (l in 1:length(split_data)){
    DefaultAssay(split_data[[l]]) <- "RNA"
  }
  
  if(integration_reduction == "CCA"){
    integ_anchors <- FindIntegrationAnchors(object.list = split_data, anchor.features = integ_features, reduction = "cca", dim = 1:50)
    data_integrated <- IntegrateData(anchorset = integ_anchors, dim = 1:50)
  }
  
  if(integration_reduction == "rPCA"){
    split_data <- lapply(X = split_data, FUN = function(x) {
      ##################################################
      ## Scaling the data
      ##################################################
      # We apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
      # The ScaleData() function: 
      #   1. Shifts the expression of each gene, so that the mean expression across cells is 0
      #   2. Scales the expression of each gene, so that the variance across cells is 1 (gives equal weight in downstream analyses, so that highly-expressed 
      #       genes do not dominate)
      #   3. The results of this are stored in s_01_data[["RNA"]]@scale.data
      x <- ScaleData(x, features = integ_features, verbose = FALSE)
      x <- RunPCA(x, features = integ_features, verbose = FALSE)
    })
    
    # Perform CCA, find the best anchors and filter incorrect anchors. (Note: the progress bar in your console will stay at 0%, but know that it is actually running.)
    # FindIntegrationAnchors: Find a set of anchors between a list of Seurat objects. 
    # normalization.method: LogNormalize or SCT
    integ_anchors <- FindIntegrationAnchors(object.list = split_data, normalization.method = "LogNormalize", anchor.features = integ_features, reduction = "rpca")
    
    # Integrate across conditions.
    # IntegrateData: Perform dataset integration using a pre-computed AnchorSet.
    data_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "LogNormalize")
  }
  
  data_integrated <- ScaleData(data_integrated, verbose = FALSE)
  
} else if (approach == "SCT"){
  for (l in 1:length(split_data)){
    DefaultAssay(split_data[[l]]) <- "SCT"
  }
  
  # Now, we need to prepare the SCTransform object for integration.
  # PrepSCTIntegration: Prepare an object list normalized with sctransform for integration.
  # Prepare the SCT list object for integration, To perform integration using the pearson residuals calculated above, we use the PrepSCTIntegration() function after selecting a list of informative features using SelectIntegrationFeatures()
  split_data <- PrepSCTIntegration(object.list = split_data, anchor.features = integ_features)
  split_data <- lapply(X = split_data, FUN = RunPCA, features = integ_features)
  integ_anchors <- FindIntegrationAnchors(object.list = split_data, normalization.method = "SCT", anchor.features = integ_features, reduction = tolower(integration_reduction))
  data_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
}

# We combine all the Doublet/Singlet info in one column
colnames <- colnames(data_integrated@meta.data[grep("DF.", colnames(data_integrated@meta.data))])
data_integrated$DF <- data.frame(col = apply(data_integrated@meta.data[colnames], 1, function(x) toString(na.omit(x))))

for (i in colnames){
  data_integrated[[i]] <- NULL
}

DF_byClsusters <- as.data.frame.matrix(table(data_integrated$seurat_clusters, data_integrated$DF)) 
write.table(DF_byClsusters, paste("3_Integration_Figures/Integration/DF_byClsusters.txt", sep =""), sep = "\t")

data_integrated@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=DF)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")
ggsave(paste("3_Integration_Figures/Integration/DF_Barplot_bySample_", project_name, ".png", sep = ""), width = 10 , height = 8)

data_integrated@meta.data %>% 
  ggplot(aes(x=pheno, fill=DF)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")
ggsave(paste("3_Integration_Figures/Integration/DF_Barplot_byPheno_", project_name, ".png", sep = ""), width = 10 , height = 8)

data_integrated@meta.data %>% 
  ggplot(aes(x=descr, fill=DF)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")
ggsave(paste("3_Integration_Figures/Integration/DF_Barplot_byDescr_", project_name, ".png", sep = ""), width = 10 , height = 8)

## among individuals
metadata <- data_integrated@meta.data
melted <- melt(metadata[which(colnames(metadata) %in% c("orig.ident","DF"))])
tab <- table(melted)
tab <- tab/apply(tab,1,sum)
melted <- melt(tab)

melted %>% 
  ggplot(aes(x=orig.ident, y=value, fill=DF)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")+
  ylab("Fraction")
ggsave(paste("3_Integration_Figures/Integration/DF_proportion_bySample_", project_name, ".png", sep = ""), width = 10 , height = 8)

## by pheno
metadata <- data_integrated@meta.data
melted <- melt(metadata[which(colnames(metadata) %in% c("pheno","DF"))])
tab <- table(melted)
tab <- (tab)/apply(tab,1,sum)
melted <- melt(tab)

melted %>% 
  ggplot(aes(x=pheno, y=value, fill=DF)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")+
  ylab("Fraction")
ggsave(paste("3_Integration_Figures/Integration/DF_proportion_byPheno_", project_name, ".png", sep = ""), width = 10 , height = 8)

####################################################################
## 3. Linear dimensional reduction
####################################################################
DefaultAssay(data_integrated) <- "integrated"

# Run PCA on the scaled data. By default, only the previously determined variable features are used as input, it is possible to chose a certain subset.
# RunPCA: Run a PCA dimensionality reduction. 
data_integrated <- RunPCA(object = data_integrated)

# Plot PCA
PCAPlot(data_integrated, group.by = "orig.ident", split.by = "orig.ident")  
ggsave(paste("3_Integration_Figures/Integration/PCAPlot_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

# From https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- data_integrated[["pca"]]@stdev / sum(data_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
max_dim <- min(co1, co2)

# Run UMAP
data_integrated <- RunUMAP(data_integrated, dims = 1:max_dim, reduction = "pca", min.dist = 0.6)

# VizDimLoadings: Visualize top genes associated with reduction components
VizDimLoadings(data_integrated, dims = 1:2, reduction = "pca")
ggsave(paste("3_Integration_Figures/Integration/VizDimLoadings_top_genes_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Plot UMAP
# DimPlot: Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a cell and it’s positioned based on the cell embeddings determined by the reduction technique. 
DimPlot(data_integrated, group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/Integration/Integration_UMAP_", project_name, ".png", sep = ""), width = 10 , height = 8)

DimPlot(data_integrated, group.by = "pheno")
ggsave(paste("3_Integration_Figures/Integration/Integration_UMAP_byPheno_", project_name, ".png", sep = ""), width = 10 , height = 8)

DimPlot(data_integrated, group.by = "DF")
ggsave(paste("3_Integration_Figures/Integration/Integration_UMAP_byDF_", project_name, ".png", sep = ""), width = 10 , height = 8)

DF_bySeurat_cluster <- as.data.frame.matrix(table(data_integrated$seurat_clusters, data_integrated$DF))
write.table(DF_bySeurat_cluster, paste("3_Integration_Figures/Integration/DF_bySeuratCluster.txt", sep =""), sep = "\t")

table(data_integrated$seurat_clusters, data_integrated$DF, data_integrated$pheno)

# Plot UMAP split by sample (Side-by-side comparison of clusters)
# Sometimes it’s easier to see whether all of the cells align well if we split the plotting between conditions, which we can do by adding the split.by argument to the DimPlot() function:
DimPlot(data_integrated, split.by = "orig.ident", group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/Integration/Integration_UMAP_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

####################################################################
## 4. Clustering
####################################################################
# To overcome the extensive technical noise in the expression of any single gene for scRNA-seq data, Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes, with each PC essentially representing a “metagene”? that combines information across a correlated gene set. Determining how many PCs to include in the clustering step is therefore important to ensure that we are capturing the majority of the variation, or cell types, present in our dataset.

# The elbow plot is a helpful way to determine how many PCs to use for clustering so that we are capturing majority of the variation in the data. The elbow plot visualizes the standard deviation of each PC, and we are looking for where the standard deviations begins to plateau. Essentially, where the elbow appears is usually the threshold for identifying the majority of the variation. However, this method can be quite subjective.

# Let’s draw an elbow plot using the top 60 PCs:
# Plot the elbow plot
ElbowPlot(object = data_integrated, ndims = 60)
ggsave(paste("3_Integration_Figures/Integration/ElbowPlot_", project_name, ".png", sep = ""), width = 10 , height = 8, bg="white")

##################################################
## Cluster the cells 
##################################################
# Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. A nice in-depth description of clustering methods is provided in the SVI Bioinformatics and Cellular Genomics Lab course.

# We will use the FindClusters() function to perform the graph-based clustering. The resolution is an important argument that sets the “granularity” of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.
# Determine the K-nearest neighbor graph
data_integrated <- FindNeighbors(object = data_integrated, dims = 1:max_dim)

# Determine the clusters for various resolutions                                
data_integrated <- FindClusters(object = data_integrated, resolution = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4))
# If we look at the metadata of our Seurat object(seurat_integrated@meta.data), there is a separate column for each of the different resolutions calculated.

# Explore resolutions only works in Rstudio local not in Nextera or Indar
#data_integrated@meta.data %>% 
#  View()
snn_dim <- rep(NA,11)
snn_dim[1] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.4))[1]
snn_dim[2] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.5))[1]
snn_dim[3] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.6))[1]
snn_dim[4] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.7))[1]
snn_dim[5] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.8))[1]
snn_dim[6] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.9))[1]
snn_dim[7] <- dim(table(data_integrated@meta.data$integrated_snn_res.1))[1]
snn_dim[8] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.1))[1]
snn_dim[9] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.2))[1]
snn_dim[10] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.3))[1]
snn_dim[11] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.4))[1]
names(snn_dim) <- seq(0.4,1.4,0.1)
print(snn_dim)

# To choose a resolution to start with, we often pick something in the middle of the range like 0.6 or 0.8. We will start with a resolution of 0.8 by assigning the identity of the clusters using the Idents() function.
# Assign identity of clusters
# //////////////////////////////////////////////////////////////////////////////
resolution_find_clusters <- "integrated_snn_res.0.8"
# //////////////////////////////////////////////////////////////////////////////
Idents(object = data_integrated) <- resolution_find_clusters

# choose how many clusters you want and assign it to seurat_clusters
print("WARNING: Choose the resolution that you want.")
# //////////////////////////////////////////////////////////////////////////////
data_integrated$seurat_clusters <- data_integrated$integrated_snn_res.0.8
# //////////////////////////////////////////////////////////////////////////////

# Plot the UMAP
DimPlot(data_integrated, reduction = "umap", label = TRUE, label.size = 6)
ggsave(paste("3_Integration_Figures/Integration/Clusters_UMAP_", project_name, ".png", sep = ""), width = 10 , height = 8)

##################################################
## Segregation of clusters by sample
##################################################
# Explore the distribution of cells per cluster in each sample:
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(data_integrated, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

print(n_cells)
# View table
#View(n_cells)

#Let’s plot the distribution among clusters using our custom function:
pdf(paste("3_Integration_Figures/Integration/Distribution_Cell_", project_name, ".pdf", sep = ""), width = 100 , height = 100)
plot_integrated_clusters(data_integrated) 
dev.off()

# We can visualize the cells per cluster for each sample using the UMAP: UMAP of cells in each cluster by sample
DimPlot(data_integrated, label = TRUE, split.by = "orig.ident")  + NoLegend()
ggsave(paste("3_Integration_Figures/Integration/Clusters_UMAP_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)
# Generally, we expect to see the majority of the cell type clusters to be present in all conditions; however, depending on the experiment we might expect to see some condition-specific cell types present. These clusters look pretty similar between conditions, which is good since we expected similar cell types to be present in both control and stimulated conditions.

##################################################
## Segregation of clusters by various sources of uninteresting variation 
##################################################
# Explore additional metrics, such as the number of UMIs and genes per cell and mitochondrial gene expression by UMAP. 
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "percent_rb")

FeaturePlot(data_integrated, reduction = "umap", features = metrics, pt.size = 0.4, order = TRUE, min.cutoff = 'q10', label = TRUE)
ggsave(paste("3_Integration_Figures/Integration/Metadata_FeaturePlot_", project_name, ".png", sep = ""), width = 12 , height = 8)

##################################################
## Exploration of the PCs driving the different clusters 
##################################################
# Explore how well the clusters separate by the different PCs. To visualize this information, we need to extract the UMAP coordinate information for the cells along with their corresponding scores for each of the PCs to view by UMAP.
# First, we identify the information we would like to extract from the Seurat object, then, we can use the FetchData() function to extract it.

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:max_dim), "ident", "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(data_integrated, vars = columns)

# In the UMAP plots, the cells are colored by their PC score for each respective principal component.
# Let’s take a quick look at the top PCs:
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(data_integrated, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
# NOTE: choose how many PCs you want to analyse 
map(paste0("PC_", 1:max_dim), function(pc){
  ggplot(pc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), alpha = 0.7) +
    scale_color_gradient(guide = FALSE, low = "grey90",  high = "blue")  +
    geom_text(data=umap_label, aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
ggsave(paste("3_Integration_Figures/Integration/PC_cluster_relation_", project_name, ".png", sep = ""), width = 12 , height = 10)

# check how cell phase behaves in each cluster  
VlnPlot(data_integrated, features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters", ncol = 2, pt.size = 0.1)
ggsave(paste("3_Integration_Figures/Integration/VlnPlot_CC_bySeuratCluster_SCT_", project_name, ".png", sep = ""), width = 12, height = 8)

data_integrated@meta.data %>% 
  ggplot(aes(x=seurat_clusters, fill=DF)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")
ggsave(paste("3_Integration_Figures/Integration/DF_Barplot_byClusters_", project_name, ".png", sep = ""), width = 10 , height = 8)

metadata <- data_integrated@meta.data
melted <- melt(metadata[which(colnames(metadata) %in% c("seurat_clusters","DF"))])
tab <- table(melted)
tab <- tab/apply(tab,1,sum)
melted <- melt(tab)

melted %>% 
  ggplot(aes(x=seurat_clusters, y=value, fill=DF)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("DF")+
  ylab("Fraction")
ggsave(paste("3_Integration_Figures/Integration/DF_proportion_byClusters_", project_name, ".png", sep = ""), width = 10 , height = 8)

data_integrated@meta.data %>% 
  ggplot(aes(x=seurat_clusters, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Cells per cluster")
ggsave(paste("3_Integration_Figures/Integration/NCells_Barplot_bySeuratCluster_", project_name, ".png", sep = ""), width = 10 , height = 8)

data_integrated@meta.data %>% 
  ggplot(aes(x=seurat_clusters, fill=pheno)) + 
  geom_bar() +
  theme_classic() +
  geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Cells per pheno")
ggsave(paste("3_Integration_Figures/Integration/NCells_Barplot_bySeuratCluster_byPheno_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Save integrated seurat object
saveRDS(data_integrated, file = paste("1_SeuratObjects/", "Data_integrated_", project_name, ".rds", sep = ""))
