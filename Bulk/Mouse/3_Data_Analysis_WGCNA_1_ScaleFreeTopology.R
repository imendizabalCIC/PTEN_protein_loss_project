################################################################################
#| SCALE-FREE TOPOLOGY ANALYSIS FOR WGCNA MOUSE DATA 
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| Determine the soft-thresholding power (beta) for WGCNA on the AC-12_RNAseq mouse 
#| cohort by assessing scale-free topology fit under a signed network, using both 
#| Pearson and biweight midcorrelation (bicor).
#|
#| Overview:
#|   * Test beta {1-10, 12, 14, ..., 30} with pickSoftThreshold().
#|   * Compute and visualize: scale-free fit (R^2) and mean connectivity vs beta.
#|   * Save plots and the pickSoftThreshold results for downstream network building.
#|
#| Workflow:
#|   1) Load preprocessed expression and traits (from previous data-cleaning step).
#|   2) Run pickSoftThreshold (networkType = "signed") with:
#|        - Pearson correlation
#|        - Bicor (biweight midcorrelation)
#|   3) Plot for each method:
#|        - beta vs scale-free fit (signed R^2; higher is better)
#|        - beta vs mean connectivity (lower connectivity with adequate fit is ideal)
#|   4) Save:
#|        - PDFs of both Pearson and bicor selection curves
#|        - RData with sft_pearson and sft_bicor objects
#|
#| Outputs
#|   - PDFs:
#|       * 1_Pearson_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity.pdf
#|       * 1_Bicor_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity.pdf
#|   - RData:
#|       * <tag>_sft_pearson_bicor.RData  (fit indices for both methods)
#|
#| Background (concise)
#|   WGCNA approximates biological networks with scale-free topology. Adjacency
#|   a_ij = s_ij^beta (soft threshold). Because log(a_ij) = beta log(s_ij), pickSoftThreshold
#|   searches beta that yields high scale-free fit (signed R^2) while retaining a
#|   reasonable mean connectivity. In weighted networks, gene connectivity is the
#|   sum of connection strengths to other nodes (not just neighbor counts).
#|
#| Notes & Recommendations
#|   * Signed vs unsigned networks:
#|       Signed networks preserve correlation sign and are generally preferred
#|       (see Peter Langfelder's guidance). Signed-hybrid and signed are nearly
#|       equivalent if you approximately double beta for signed vs signed-hybrid.
#|   * Correlation choice:
#|       Bicor is robust to outliers; evaluate both Pearson and bicor and select
#|       the beta that achieves adequate signed R^2 (often >= 0.8) with acceptable
#|       mean connectivity for your data.
################################################################################


###############################  LIBRARIES  ####################################
suppressMessages(library(WGCNA))

options(stringsAsFactors = FALSE)
################################################################################



#################################  DATA  #######################################
dir.proj <- "X:/irondon/AC-12_RNAseq/07_WGCNA/"
info.file <- "X:/irondon/AC-12_RNAseq/04_DEGs/Data/sample_info_AC-12_RNAseq.txt"
counts.file <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt"
results.file <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Images/01_ScaleFreeTopology/"
data.file <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Data/"
tag <- "AC-12_RNAseq"

setwd(dir.proj)

#| Loading the data
load(file =  paste(data.file, tag, "_dataInput.RData", sep =""))
################################################################################



#####################  CHOOSING THE SOFT-THRESHOLDING POWER:  ##################

# Choose a set of soft-thresholding powers to test
powers <- c(c(1:10), seq(from=10, to=30, by=2))


#| Call the network topology analysis function
sft_pearson = pickSoftThreshold(datExpr, powerVector = powers, networkType ="signed", verbose = 5)
sft_bicor = pickSoftThreshold(datExpr, powerVector = powers, corFnc = "bicor", networkType ="signed", verbose = 5)

#| Scale-free topology fit index as a function of the soft-thresholding power
#                          Power 
power_plot_pearson = sft_pearson$fitIndices[,1]
power_plot_bicor = sft_bicor$fitIndices[,1]

#                           slope              SFT.R.sq(related to the )
slope_SFT.R.sq_pearson = -sign(sft_pearson$fitIndices[,3])*sft_pearson$fitIndices[,2]
slope_SFT.R.sq_bicor = -sign(sft_bicor$fitIndices[,3])*sft_bicor$fitIndices[,2]

mean.k._pearson = sft_pearson$fitIndices[,5]
mean.k._bicor = sft_bicor$fitIndices[,5]
################################################################################


##########################  PLOTTING PEARSON: ##################################
#R^2
sizeGrWindow(9,5)
pdf(file = paste(results.file,"1_Pearson_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity_last_try.pdf", sep = ""), width = 9, height = 5)
par(mfrow=c(1,2))
cex1= 0.9
plot(power_plot_pearson, slope_SFT.R.sq_pearson, family ="sans",
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main= paste("Scale independence Pearson"))
text(power_plot_bicor, slope_SFT.R.sq_pearson, family="serif",
     labels=powers,cex=cex1,col="red");
abline(h=max(slope_SFT.R.sq_pearson),col="blue")
#MEAN CONNECTIVITY
plot(power_plot_pearson, mean.k._pearson, family ="sans",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Pearson"))
text(power_plot_pearson, mean.k._pearson, family ="sans", labels=powers, cex=cex1,col="red")
dev.off()
################################################################################


###########################  PLOTTING BICOR: ###################################
#R^2
sizeGrWindow(9,5)
pdf(file = paste(results.file, "1_Bicor_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity_last_try.pdf", sep =""), width = 8, height = 4.5)
par(mfrow=c(1,2))
cex1= 0.9
plot(power_plot_bicor, slope_SFT.R.sq_bicor, family ="sans",
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main= paste("Scale independence Bicor"))
text(power_plot_bicor, slope_SFT.R.sq_bicor, family="sans",
     labels=powers,cex=cex1,col="red");
abline(h=max(slope_SFT.R.sq_bicor),col="blue")
#MEAN CONNECTIVITY
plot(power_plot_bicor, mean.k._bicor, family ="sans",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Bicor"))
text(power_plot_bicor, mean.k._bicor, family ="sans", labels=powers, cex=cex1,col="red")
dev.off()
################################################################################


#| Saving the analysis of bicor and pearson soft thresholds
save(sft_pearson, sft_bicor, file = paste(data.file, tag, "_sft_pearson_bicor.RData", sep =""))




