##############################################################################################################
########################  NETWORK CONSTRUCTION OF THE AC-45_RNAseq-FFPE DATA  ################################
##############################################################################################################

#  1) weights = Prior information of weights in the same format (and dimensions) as 'datExp' Default NULL.

#  2) checkMissingData = Logical Parameter. Default TRUE.

#  3) blocks = If given, must be a numeric vector with one entry per column (gene) of exprData giving the number 
#  of the block to which the corresponding gene belongs. DEFAULT NULL.

#  4) maxBlockSize = This is the maximum block size for module detection. If there are too many genes and they 
#  exceed the maxBlockSize the function will divide the modules into blocks that does not exceed the maxBlockSize
#  defined. DEFAULT 5000.

#  5) blockSizePenaltyPower = This is related by how strong must be penalized for exceeding the maximum size. Set 
#  this to a large number is very important because we don't want to exceed the maxblock size, otherwise the
#  precision of obtaining true positive could decrease. DEFAULT 5.

#  6) nPreclusteringCenters = Number of centers for pre-clustering. Larger numbers typically results in better but 
#  slower pre-clustering. DEFAULT as.integer(min(ncol(datExpr)/20, 100*ncol(datExpr)/maxBlockSize)),

#  7) randomSeed = For reproducibility. DEFAULT 54321.

#  8) loadTOM = It may be useful to load previously saved TOM matrices if these have been calculated previously, 
#  since TOM calculation is often the most computationally expensive part of network construction and module 
#  identification. Default FALSE. But how this is affected when the value of the parameters are changed?

#  9) corType = Options "pearson" and "bicor". It is recommended to use "bicor" because it is less sensitive to 
#  outliers. The different with pearson is that bicor it is median-based rather that mean-based.
#  DEFAULT "pearson". 

#  10) maxPOutliers = Only used when 'corType = "bicor"'. Specifies the maximum percentile of data that can be considered
#  outlier on either side of the median. Using 'maxPOutlier = 1'  will disable all weight function broadening
#  using 'maxPOutlier = 0' will give results similar to Pearson Correlation. We strongly recommend using the 
#  argument maxPOutliers = 0.05 or 0.10 whenever the biweight midcorrelation is used. This argument essentially 
#  forces bicor to never regard more than the specified proportion of samples as outliers. DEFAULT 1.

#  11) quickCor = Real number between 0 and 1 that controls the handling of missing data. DEFAULT 0.

#  12) pearsonFallback = This parameter is to specify if one should return to pearson correlation (is bicor is used)
#  when the median absolute deviation is zero. Possible values: "none", "individual" and "all". 
#  DEFAULT "individual".

#  13) cosineCorrelation = Logical. The cosine calculation differs from the standard one in that it does not substract
#  the mean. Cosine similarity does not depend on the magnitudes of the vectors, but only on their angle. It is particularly
#  used in positive space. One advantage of cosine similarity is its low complexity, especially for sparse vectors: 
#  only the non-zero coordinates need to be considered. The fisical interpretation is: The cosine similarity computes 
#  the similarity between two samples, whereas the Pearson correlation coefficient computes the correlation between 
#  two jointly distributed random variables. DEFAULT FALSE.

#  14) power = Soft-thresholding power for network construction. DEFAULT 6.

#  15) networkType = Values: "unsigned", "signed", "signed hybrid". DEFAULT "unsigned".

#  16) replaceMissingAdjacencies = Logical. Should missing values in the calculation of adjacency be replaced by 0?. 
#  DEFAULT FALSE

#  17) TOMType = one of "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2" and "signed Nowick 2". 
#  If "none", adjacency will be used for clustering. DEFAULT "signed". The Topological Overlap Measure (TOM) was 
#  introduced  to make networks less sensitive to spurious connections or to connections missing due to random noise.
#  The central idea of TOM is to count the direct connection strengths as well as connection strengths "mediated" by 
#  shared neighbors. The standard, or "unsigned" TOM assumes that neighbor-mediated connections can always be considered 
#  as "reinforcing" the direct connection. This may not always be the case, and the signed TOM is an attempt to take this 
#  into account. Since the anti-reinforcing connection strengths (practically) cannot occur in signed networks, in signed 
#  networks the signed and unsigned TOM are (practically) identical. (see information in PDF)

#  TOMDenom = A character string specifying the TOM variant to be used. DEFAULT "min".  

#  suppressTOMForZeroAdjacencies = Logical. Should TOM be set to zero for zero adjacencies?. DEFAULT FALSE.

#  suppressNegativeTOM = should the result be set to zero when negative? Negative TOM values can occur when TOMType is 
#  "signed Nowick". DEFAULT FALSE.

#  getTOMs = DEFAULT NULL.

#  saveTOMs = Logical. should the consensus topological overlap matrices for each block be saved and returned? DEFAULT FALSE.

#  saveTOMFileBase = character string containing the file name base for files containing the consensus topological overlaps. 
#  DEFAULT "blockwiseTOM".

#  deepSplit = integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to
#  module splitting, with 0 least and 4 most sensitive. DEFAULT 2. In other words, the higher the deepSplit, the more modules 
#  I'll obtain. And they will be smaller.

#  detectCutHeight = Dendrogram cut height for module detection. Selecting deeper cut heights would result in loss of some 
#  other branch information DEFAULT 0.995.

#  minModuleSize = DEFAULT min(20, ncol(datExpr)/2 )

#  maxCoreScatter = DEFAULT NULL, minGap = DEFAULT NULL

#  maxAbsCoreScatter = DEFAULT NULL, minAbsGap = DEFAULT NULL

#  minSplitHeight = DEFAULT NULL, minAbsSplitHeight = DEFAULT NULL

#  useBranchEigennodeDissim = DEFAULT FALSE,

#  minBranchEigennodeDissim = DEFAULT mergeCutHeight,

#  stabilityLabels = DEFAULT NULL.

#  stabilityCriterion = DEFAULT c("Individual fraction", "Common fraction").

#  minStabilityDissim = DEFAULT NULL.

#  pamStage = TRUE, pamRespectsDendro = DEFAULT TRUE.

#  reassignThreshold = p-value ratio threshold for reassigning genes between modules. DEFAULT 1e-6.

#  minCoreKME = a number between 0 and 1. If a detected module does not have at least minModuleKMESize genes with eigengene
#  connectivity at least minCoreKME, the module is disbanded DEFAULT 0.5. 

#  minCoreKMESize = DEFAULT minModuleSize/3.

#  minKMEtoStay = genes whose eigengene connectivity to their module eigengene is lower than minKMEtoStay are removed from 
#  the module. DEFAULT 0.3.

#  mergeCutHeight = Dendrogram cut height for module merging. DEFAULT 0.15.  

#  impute = Logical. Should imputation be used for module eigengene calculation?. DEFAULT TRUE

#  trapErrors = Logical. should errors in calculations be trapped?. DEFAULT FALSE.

#  numericLabels = Logical. Should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)? DEFAULT FALSE.

#  nThreads = This is a non-negative integer specifying the number of parallel threads to be used by certain parts of correlation 
#  calculations. If zero, the number of online processors will be used if it can be determined dynamically, otherwise correlation
#  calculations will use 2 threads. DEFAULT 0.

#  useInternalMatrixAlgebra = Logical: should WGCNA's own, slow, matrix multiplication be used instead of R-wide BLAS? Only useful
#  for debugging. DEFAULT FALSE.

#  useCorOptionsThroughout = Logical: should correlation options passed to network analysis also be used in calculation of kME?
#  Set to FALSE to reproduce results obtained with WGCNA 1.62 and older. DEFAULT TRUE,

#  verbose = integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose.
#  DEFAULT 0, indent = 0,

###############################################################################################################
################################# POSSIBLE VARIABLES TO CHANGE ################################################
###############################################################################################################

#         *corType (pearson, bicor)  *deepSplit (2,4)  *nThreads (12, 14)   *mergeCutHeight (0.2, 0.25)    
#         *minKMEtoStay     *networkType (signed)    *minCoreKME  *maxPOutliers (0.05,0.1)  * maxBlocksize

#         *Change the cosineCorrelation parameter?     
#         *replaceMissingAdjacencies = FALSE is not TRUE better?

# Q: TOMtype singed, unsinged or singed hybridic? 
# R: if you want to have genes that follow the same pattern then use signed hybrid, otherwise use unsigned, or if
# you want to be in the middle use the "signed" option

# Q: How to choose the mergeCutHeight parameter?
# R: In practice, on datasets of 50-100 samples using 0.2 or 0.15 works fairly well. For fewer samples a larger 
# valuse (0.25 to 0.3) may be warranted. If you want larger modules, increase the value; if you want smaller modules 
# at the risk of having redundant modules (modules with very similar functional annotation and very similar fuzzy
# module membership), you can decrease the value to say 0.10, maybe even lower if you have lots of samples (hundreds
# or more).

###############################################################################################################


start_time <- Sys.time()
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#=============================================================================================================
#
# blockwiseModules: Function to calculate modules and eigengenes from all genes.
#
#=============================================================================================================
# 
data <- load(file="../Results/Data/AC-12_RNAseq_dataInput.RData")
net <- blockwiseModules(
                       # Input data
                       datExpr, 
                       weights = NULL,
                       
                       # Data checking options
                       checkMissingData = TRUE,

                       # Options for splitting data into blocks
                       blocks = NULL,
                       maxBlockSize=20000,    
                       #blockSizePenaltyPower = 7, #The higher the better when there is not computational power
                       #randomSeed = 54321, #Default

                       # Load TOM from previously saved file?
                       loadTOM = FALSE,

                       # Network construction arguments: correlation options
                       corType = "bicor", 
                       maxPOutliers = 0.1,
                       pearsonFallback = "individual",
                       cosineCorrelation =FALSE,

                       # Adjacency function options
                       power = 18,
                       networkType ="signed", 
                       replaceMissingAdjacencies = FALSE,

                       # Topological overlap options
                       TOMType = "signed",  # When we are considering "signed" networks, this parameters does not matter a lot
                       TOMDenom="min",
                       suppressTOMForZeroAdjacencies = FALSE,
                       suppressNegativeTOM = FALSE,

                       # Saving or returning TOM
                       getTOMs = NULL,
                       saveTOMs = TRUE, 
                       saveTOMFileBase = "../Results/Data/2_blockwiseModules_TOM_AC-12_RNAseq",

                       # Basic tree cut options
                       deepSplit = 2, # Could be worthy increasing it when there is a huge amount of genes in our data set?
                       #detectCutHeight = 0.99,
                       minModuleSize =  min(25, ncol(datExpr)/2), 

                       # Module merging options
                       mergeCutHeight = 0.25, 
                       impute = TRUE, 
                       trapErrors = FALSE,

                       reassignThreshold = 1e-6, 
                       
                       # Output options
                       numericLabels = TRUE, 

                        # Advanced tree cut options
                       pamRespectsDendro = FALSE, 

                       # Options controlling behaviour
                       nThreads = 5, 
                       verbvose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(net, MEs, moduleLabels, moduleColors, geneTree, file = "../Results/Data/2_blockwiseModules_AC_12_RNAseq.RData")

end_time <- Sys.time()
ex_time <- end_time - start_time
write.csv(ex_time, "../Results/Data/ex_time_bicor.txt")
