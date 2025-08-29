################################################################################
#| DATA INPUT AND CLEAN FOR WGCNA ANALYSIS              
################################################################################
#| Date 12/12/2024
#| Author: Ivana Rondon-Lorefice
#|
#| Description:
#| This script performs preprocessing of RNA-seq expression and trait data to 
#| prepare them for Weighted Gene Co-expression Network Analysis (WGCNA). The 
#| workflow includes transformation, filtering, clustering, and saving cleaned 
#| data in a suitable format for downstream network construction.
#|
#| Main steps:
#|   1) Import expression counts and sample trait metadata.
#|   2) Transform categorical variables into binary/numeric representations.
#|   3) Filter low-count genes and apply variance-stabilizing transformation (VST).
#|   4) Detect and remove samples/genes with excessive missing values or outliers.
#|   5) Perform hierarchical clustering of samples to visualize structure and detect outliers.
#|   6) Plot sample dendrograms with trait heatmaps for QC.
#|   7) Save processed expression and trait data in an RData object for WGCNA input.
#|
#| Output:
#|   - Sample clustering plots (PDF).
#|   - Sample dendrogram with trait heatmap (PDF).
#|   - Cleaned expression matrix and trait data (RData file).
################################################################################


################################################################################
#| LIBRARIES AND DATA
################################################################################

suppressMessages(library(WGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(skimr))
suppressMessages(library(DataExplorer))
suppressMessages(library(ggraph))
suppressMessages(library(igraph))
suppressMessages(library(tidyverse))
suppressMessages(library(dendextend))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(DESeq2))

options(stringsAsFactors = FALSE)

#| File directories
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/Xcell_output.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt"
result.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Images/0_DataInput_and_Clean/"
data.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Data/"
tag <- "AC-45_RNAseq-FFPE"

#| Setting working directory
setwd(dir.proj)

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

#| vst transformation
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = datTraits, 
                              design = ~ 1 ) #| No specify the model. Just do vst
dds <- estimateSizeFactors(dds)
vsd.noblind <- vst(dds, blind=FALSE)
vsd.blind <- vst(dds, blind=TRUE)
any(assay(vsd.blind) == assay(vsd.noblind)) #| When the design is not given, it seems that these values are the same
################################################################################


################################################################################
#| DEALING WITH CATEGORICAL VARIABLES   
################################################################################

#|  NOTE:
#| In the case of the variable group (or those which have more than two levels) 
#| it is better to separate then into groups of categorical variable using the 
#| WGCNA function called binarizeCategoricalVariable()
#| binarize_out = binarizeCategoricalVariable(datTraits$Group, includePairwise = TRUE)
#| datTraits <- cbind(datTraits, binarize_out)


#| Assigning 0 to PTEN loss and 1 to PTEN presence CUT 0
datTraits$H_score_cut_0[datTraits$H_score_cut_0 == "PTEN loss"] <- 1
datTraits$H_score_cut_0[datTraits$H_score_cut_0 == "PTEN presence"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN presence CUT 10
datTraits$H_score_cut_10[datTraits$H_score_cut_10 == "PTEN loss"] <- 1
datTraits$H_score_cut_10[datTraits$H_score_cut_10 == "PTEN presence"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN presence CUT 20
datTraits$H_score_cut_20[datTraits$H_score_cut_20 == "PTEN loss"] <- 1
datTraits$H_score_cut_20[datTraits$H_score_cut_20 == "PTEN presence"] <- 0

#| Assigning 0 to PTEN loss and 1 to PTEN presence CUT 30
datTraits$H_score_cut_30[datTraits$H_score_cut_30 == "PTEN loss"] <- 1
datTraits$H_score_cut_30[datTraits$H_score_cut_30 == "PTEN presence"] <- 0

#| Assigning 0 to No info of pAKT and 1 to positive for pAKT 
datTraits$pAKT_positive[which(datTraits$pAKT_positive == "No pAKT")] <- 0
datTraits$pAKT_positive[which(datTraits$pAKT_positive == "pos")] <- 1

#| Converting to numeric values
datTraits$H_score_cut_0 <-  as.numeric(datTraits$H_score_cut_0)
datTraits$H_score_cut_10 <- as.numeric(datTraits$H_score_cut_10)
datTraits$H_score_cut_20 <- as.numeric(datTraits$H_score_cut_20)
datTraits$H_score_cut_30 <- as.numeric(datTraits$H_score_cut_30)
datTraits$pAKT_positive <- as.numeric(datTraits$pAKT_positive)

#| Selecting the columns of interest
datTraits <- subset( datTraits, select = c( Edad, 
                                            PTEN_Exp_log2, 
                                            DV200.value, 
                                            DFS.TIME,
                                            DFS.STATUS,
                                            purity,
                                            stromal,
                                            immune,
                                            Mean.expression.FOXO,
                                            Mean.expression.FOXO.pathway,
                                            Mean.expression.PTEN_loss,
                                            Mean.expression.PI3K.AKT.mTOR,
                                            H_score,
                                            H_score_cut_0,
                                            H_score_cut_10,
                                            H_score_cut_20,
                                            H_score_cut_30,
                                            pAKT_positive,
                                            Macrophages.M1,
                                            Macrophages,
                                            StromaScore,
                                            Smooth.muscle,
                                            NKT,
                                            Epithelial.cells,
                                            Fibroblasts,
                                            MicroenvironmentScore))

colnames(datTraits)[1] <- "Age"
colnames(datTraits)[14] <- "PTEN_status"

################################################################################


################################################################################
#|    DATA CLEAN   
################################################################################

#| We need to transpose the expression data for further analysis.
datExpr <- as.data.frame(t(assay(vsd.noblind)))

#| Checking data for excessive missing values and identification of oulier microarray 
#| samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

#| IF the last statement is TRUE, all genes have passed the cuts. If not, we remove 
#| the offending genes and samples from the data 
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
################################################################################


################################################################################
#| CLUSTERING DENDOGRAMS 
################################################################################

#| Next, we cluster the samples (which is not the same as the clustering genes 
#| that will come later) to see if the are any obvious outliers
sampleTree <- hclust(dist(datExpr)/max(dist(datExpr)), method ="average")

#| Customizing
hcd <- as.dendrogram(sampleTree)
hcd <-highlight_branches_col(hcd,rev(colorRampPalette(c(brewer.pal(4, "Spectral")))(10)))
hcd <- set(hcd,"labels_cex",0.26)

#| Plot the sample tree: Open a graphic output window of size 12 by 9 inches the 
#| user should change the dimensions if the window is too large or too small.
pdf(file = paste(result.file,"0_sampleClustering.pdf", sep =""))
plot(hcd, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 1.4,
     cex=0.1,
     family ="serif", ylim= c(0,1.1))
dev.off()

#| Plot without restriction
pdf(file = paste(result.file,"0_sampleClustering_2.pdf", sep =""))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 1.4,
     cex=0.28,
     family ="serif", ylim= c(0,1.1))
dev.off()
################################################################################


################################################################################
#| CHECKING OUTLIERS
################################################################################

outlier <- F
if (outlier){
  
  #| IF TRUE: We can remove it by hand or by choosing a height cut that will remove 
  #| the offending sample. We can draw a line and determine the cluster under the 
  
  h <- 90
  color <- "red"
  abline(h=h, col=color)
  clust <- cutreeStatic(sampleTree, cutHeight = 90, minSize = 60)
  table(clust)
  
  #| IF clust 1 contains the samples we want to keep we do:
  keepSamples <- (clust == 1) 
  datExpr <- datExpr[keepSamples,]
}
################################################################################


################################################################################
#| PLOTTING DENDOGRAMS AND SAVING RESULTS
################################################################################

#| Before continuing with network construction and module detection, we visualize 
#| how the clinical traits relate to the sample dendrogram. But first we need to 
#| convert traits to a color representation: white means low, read high, grey 
#| means missing entry
collectGarbage()
  
#| Associating numbers to colors
traitColors <- numbers2colors(datTraits, signed = FALSE, colors =colorRampPalette(c(brewer.pal(4, "Spectral")))(100))
traitColors <- numbers2colors(datTraits, signed = FALSE, colors =viridis(100, option ="B"))

#| Plot the sample dendrogram and the colors underneath
pdf(paste(result.file,"0_sample_dendrogram_and_trait_heatmap.pdf", sep =""))
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", 
                    cex.colorLabels = 0.5, 
                    cex.dendroLabels = 0.2, 
                    family ="serif")
dev.off()

#| Saving data
save(datExpr, datTraits, file = paste(data.file, tag, "_dataInput_last_try.RData", sep =""))
################################################################################




