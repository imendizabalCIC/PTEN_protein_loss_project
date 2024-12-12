################################################################################
###### SCALE-FREE TOPOLOGY FOR JOINT ANALYSIS OF DATA PTEN PROTEIN LOSS ########
################################################################################

#| WGCNA assumes that the topology of the network being constructed follows the
#| scale free topology. For instance, it is important to find the correct soft
#| threshold who best fit to the data in consideration. The assumption made in 
#| WGCNA is that an element in the adjacent matrix behaves by the following rule:
#| a_{ij} = (s_{ij})^beta, where beta is the soft-thresholding power and is the
#| parameter searching in this section. In this way, we can construct a weighted
#| network.

#| Because log(a_{ij}) = beta*log(s_{ij})), one can infer the value of beta, where
#| the idea is to try to find the value that fits the best to the scale-free to-
#| pology. Note that the GENE CONNECTIVIY it is different for weighted and un-
#| weighted networks:

#|  - For unweighted networks = number of direct neighbors

#|  - For weighted networks = sum of connection strengths to other nodes

#| Scale Free Topology refers to the frequency distribution of the connectivity k
#| p(k)=proportion of nodes that have connectivity k. So the linear model is fit
#| by Log transforming p(k) and k and look at scatter plots. This adjustment
#| is perform by the pickSoftThreshold function.

#| NOTE: To edit a given function use trace(scaleFreeFitIndex, edit = TRUE)

#| Additionally, here is tested the pearson and the bicor approach. See the comments
#| bellow

#| Call the network topology analysis function
#| To decide the kind of network please read https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/ 
#| In this blog, Peter explains why he thinks taking signed networks in better than
#| unsigned ones. In other words, using the absolute value of the correlation may 
#| obfuscate biologically relevant information, since no distinction is made between
#| gene repression and activation. In contrast, in signed networks the similarity 
#| between genes reflects the sign of the correlation of their expression profiles.

#|  Signed hybridic? https://peterlangfelder.com/2018/11/25/__trashed/ 
#| "Thus, in the end the two signed network variants result in nearly identical 
#| networks as long as (1) the signed network uses twice the soft thresholding 
#| power of the signed hybrid network, and (2) the power is suitable for analysis 
#| of high-throughput, more-variables-than-samples, data: at least 4 for signed 
#| hybrid, and at least 8 for signed networks. In other words, use either but 
#| remember to double the power for the signed, compared to the signed hybrid"

#|  Another decision we have to make is which type of correlation we must use, i.e. 
#| Pearson or Bicor? It is recommended to do the analysis on bicor but it is better
#| to perform the analysis of the softhresold on both correlation functions.

################################################################################

################################################################################
#| LIBRARIES AND DATA
################################################################################

suppressMessages(library(WGCNA))

options(stringsAsFactors = FALSE)

dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Images/1_ScaleFreeTopology"
data.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Data/"

tag <- "AC-45_RNAseq-FFPE"

setwd(dir.proj)

#| Loading the data
load(file =  paste(data.file, tag, "_dataInput_last_try.RData", sep =""))

################################################################################


################################################################################
#| CHOOSING THE SOFT-THRESHOLDING POWER:
################################################################################

# Choose a set of soft-thresholding powers to test
powers <- c(c(1:10), seq(from=10, to=30, by=2))

#| Call the network topology analysis function
sft_pearson = pickSoftThreshold(datExpr, powerVector = powers, networkType ="signed", verbose = 5)
sft_bicor = pickSoftThreshold(datExpr, powerVector = powers, corFnc = "bicor", networkType ="signed", verbose = 5)

#| Scale-free topology fit index as a function of the soft-thresholding power 
power_plot_pearson = sft_pearson$fitIndices[,1]
power_plot_bicor = sft_bicor$fitIndices[,1]

#| slope, SFT.R.sq(related to the )
slope_SFT.R.sq_pearson = -sign(sft_pearson$fitIndices[,3])*sft_pearson$fitIndices[,2]
slope_SFT.R.sq_bicor = -sign(sft_bicor$fitIndices[,3])*sft_bicor$fitIndices[,2]

#| mean k
mean.k._pearson = sft_pearson$fitIndices[,5]
mean.k._bicor = sft_bicor$fitIndices[,5]

################################################################################


################################################################################
#|  PLOTTING PEARSON:
################################################################################

#| R^2
sizeGrWindow(9,5)
pdf(file = paste(results.file,"1_Pearson_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity_last_try.pdf", sep = ""), width = 9, height = 5)
par(mfrow=c(1,2))
cex1= 0.9
plot(power_plot_pearson, slope_SFT.R.sq_pearson, family ="serif",
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main= paste("Scale independence Pearson"))
text(power_plot_bicor, slope_SFT.R.sq_pearson, family="serif",
     labels=powers,cex=cex1,col="red");
abline(h=max(slope_SFT.R.sq_pearson),col="blue")

#| MEAN CONNECTIVITY
plot(power_plot_pearson, mean.k._pearson, family ="serif",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Pearson"))
text(power_plot_pearson, mean.k._pearson, family ="serif", labels=powers, cex=cex1,col="red")
dev.off()

################################################################################


################################################################################
#|  PLOTTING BICOR:
################################################################################

#| R^2
sizeGrWindow(9,5)
pdf(file = paste(results.file, "1_Bicor_SoftThreshold_vs_ScaleFreeTopology_and_MeanConnectivity_last_try.pdf", sep =""), width = 9, height = 5)
par(mfrow=c(1,2))
cex1= 0.9
plot(power_plot_bicor, slope_SFT.R.sq_bicor, family ="serif",
     xlab="Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main= paste("Scale independence Bicor"))
text(power_plot_bicor, slope_SFT.R.sq_bicor, family="serif",
     labels=powers,cex=cex1,col="red");
abline(h=max(slope_SFT.R.sq_bicor),col="blue")

#| MEAN CONNECTIVITY
plot(power_plot_bicor, mean.k._bicor, family ="serif",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Bicor"))
text(power_plot_bicor, mean.k._bicor, family ="serif", labels=powers, cex=cex1,col="red")
dev.off()

################################################################################


#| Saving the analysis of bicor and pearson soft thresholds
save(sft_pearson, sft_bicor, file = paste(data.file, tag, "_sft_pearson_bicor.RData", sep =""))




