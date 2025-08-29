
##############################################################################################################
#|  NETWORK CONSTRUCTION FOR PROJECT PTEN PROTEIN LOSS (WGCNA)
##############################################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| Build a signed gene co-expression network with WGCNA’s `blockwiseModules` for the
#| AC-45_RNAseq-FFPE cohort (PTEN protein loss vs presence). Uses biweight mid-correlation
#| (bicor) with robust outlier handling and saves module assignments, eigengenes, dendrograms,
#| and TOM files for downstream analyses.
#|
#| Key settings:
#|   - corType = "bicor", maxPOutliers = 0.1, pearsonFallback = "individual"
#|   - networkType = "signed", TOMType = "signed", power = 14
#|   - deepSplit = 1 (conservative module splitting), mergeCutHeight = 0.1 (aggressive merge)
#|   - maxBlockSize = 25,000, minModuleSize = min(25, ncol(datExpr)/2)
#|   - numericLabels = TRUE, nThreads = 5
#|
#| Workflow:
#|   1) Load preprocessed expression matrix (`datExpr`) from the data input/cleaning step.
#|   2) Run `blockwiseModules` with the settings above to construct the signed network.
#|   3) Save: network object (`net`), module eigengenes (MEs), colors/labels, and gene tree.
#|   4) Persist TOMs to disk for reproducibility and later module refinement.
#|   5) Log total execution time.
#|
#| Outputs:
#|   - ../Results/Data/2_blockwiseModules_TOM_AC_45_RNAseq_FFPE_last_try.*   (TOM files)
#|   - ../Results/Data/2_blockwiseModules_AC_45_RNAseq_FFPE_last_try.RData   (net, MEs, labels, colors, tree)
#|   - ../Results/Data/ex_time_bicor_last_try.txt                            (runtime log)
###############################################################################################################


###############################################################################################################
# RUNNING blockwiseModules FUNCTION
###############################################################################################################

#| LOADING WGCNA LIBRARY AND SETTING START SYSTEM TIME
start_time <- Sys.time()
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

data <- load(file="../Results/Data/AC-45_RNAseq-FFPE_dataInput_last_try.RData")
net <- blockwiseModules(
                       # Input data
                       datExpr, 
                       weights = NULL,
                       
                       # Data checking options
                       checkMissingData = TRUE,

                       # Options for splitting data into blocks
                       blocks = NULL,
                       maxBlockSize=25000,    

                       # Load TOM from previously saved file?
                       loadTOM = FALSE,

                       # Network construction arguments: correlation options
                       corType = "bicor", 
                       maxPOutliers = 0.1,
                       pearsonFallback = "individual",
                       cosineCorrelation =FALSE,

                       # Adjacency function options
                       power = 14,
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
                       saveTOMFileBase = "../Results/Data/2_blockwiseModules_TOM_AC_45_RNAseq_FFPE_last_try",

                       # Basic tree cut options
                       deepSplit = 1, 
                       minModuleSize =  min(25, ncol(datExpr)/2), 

                       # Module merging options
                       mergeCutHeight = 0.1, 
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
save(net, MEs, moduleLabels, moduleColors, geneTree, file = "../Results/Data/2_blockwiseModules_AC_45_RNAseq_FFPE_last_try.RData")

end_time <- Sys.time()
ex_time <- end_time - start_time
write.csv(ex_time, "../Results/Data/ex_time_bicor_last_try.txt")
###############################################################################################################