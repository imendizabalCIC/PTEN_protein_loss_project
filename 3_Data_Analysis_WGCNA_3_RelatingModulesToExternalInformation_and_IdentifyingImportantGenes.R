################################################################################
#####          RELATING MODULES TO EXTERNAL INFORMATION WITH DATA          #####
################################################################################

#| In this script, we correlate network parameters with traits of interest. Specifically, 
#| we correlate module eigengenes with clinical traits, quantify module membership, 
#| and assess gene significance. We also plot the association of each module with the loss 
#| of PTEN protein to identify those that are highly correlated with this trait.

################################################################################


################################################################################
#|  LIBRARIES  AND DATA
################################################################################
suppressMessages(library(WGCNA))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(corrplot))
options(stringsAsFactors = FALSE)

dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/"
results.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Images/3_Module-trait relationships/"
data.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Data/"
table.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/5_STEP_WGCNA/5_WGCNA_PCa_PTEN_protein_loss_vs_presence/Results/Tables/"
data.file_DEGs <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Wilcoxon_LimmaVoom_EdgeR_DESeq2_DEGs/Tables/Limma_voom_analysis_DEGS_results.txt"
tag <- "AC-45_RNAseq-FFPE"

#| Setting working Dir
setwd(dir.proj)

#| Loading Data
load(file = paste(data.file, "2_blockwiseModules_AC_45_RNAseq_FFPE_last_try.RData", sep =""))
load(file = paste(data.file, tag, "_dataInput_last_try.RData", sep =""))

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| For plots
theme_set(theme_classic())
################################################################################


################################################################################
#|  GENE DENDOGRAMS AND MODULE COLORS
################################################################################

#| Relabel blockwise modules
bwLabels = matchLabels(net$colors, moduleLabels)

#| Number of modules identified
nModules <- length(table(bwLabels))

#| Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

#| Plotting dendogram and module colors
sizeGrWindow(12,9)
pdf(paste(results.file,"3_SingleBlock_geneDendogram_and_ModuleColors_new.pdf", sep =""))
plotDendroAndColors(geneTree,
                    cbind(bwModuleColors),
                    c("Modules"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#| Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

length(moduleColors)
#| Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
################################################################################


################################################################################
#|  ASSOCIATION OF MODULES WITH TRAITS
################################################################################

Association_function <- function(MEs, Trait, nSamples,nModules, bwModuleColors,label,results.file){
  
  #| Correlation of Modules eigengene of every module with the trait
  moduleTraitCor <- cor(MEs, Trait, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  #| Changing rownames 
  rownames(moduleTraitCor) <- gsub("ME", "", rownames(moduleTraitCor))
  rownames(moduleTraitPvalue) <- gsub("ME", "", rownames(moduleTraitPvalue))
  
  if(any(rownames(moduleTraitCor) == rownames(moduleTraitPvalue))){
    
    size <- c()
    for (i in 1:nModules){
      size <- c(size,length(bwModuleColors[which(bwModuleColors == rownames(moduleTraitCor)[i])])) 
    }
    
    #| Creating data frame with values
    log_pvalue_association <- data.frame(pvalue = moduleTraitPvalue ,
                                         Association = moduleTraitCor,
                                         module_size = size,
                                         Modules = rownames(moduleTraitCor))
    
    #| Changing module size to numeric
    log_pvalue_association$module_size <- as.numeric(log_pvalue_association$module_size)
    
    modules <- rownames(moduleTraitCor)
    modules <- as.character(modules)
    modules <- sort(modules)
    
    #| Plot and save
    ggplot(log_pvalue_association, aes(x = Association, y = -log10(pvalue))) + 
      geom_point(aes(size =module_size, fill=Modules), colour="black",pch=21,alpha = 0.8) +
      theme(text=element_text(size=12,  family="sans"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
      scale_fill_manual(values = modules, guide = "legend") +
      labs(size ="Module size") +
      xlab(paste("Association with ",label,sep="")) +
      guides(fill = "none")
    ggsave(paste(results.file,"3_Association_with_", label,"_logPvalue.pdf", sep =""), height = 3.1, width = 4.2)
    
  }else{
    
    print("Check rownames(moduleTraitCor) == rownames(moduleTraitPvalue)")
    
  }
  
}

#| Association with PTEN protein loss
Association_function(MEs, datTraits$PTEN_status, nSamples,nModules,bwModuleColors, "PTEN protein loss", results.file)

#| Association with H-score
Association_function(MEs, datTraits$H_score, nSamples,nModules,bwModuleColors, "H-score", results.file)

#| Association with stromal score (obtained from tidyestimate)
Association_function(MEs, datTraits$stromal, nSamples,nModules,bwModuleColors, "Stromal_score_tidyestimate", results.file)

#| Association with immune score (obtained from tidyestimate)
Association_function(MEs, datTraits$immune, nSamples,nModules,bwModuleColors, "Immune_score_tidyestimate", results.file)

#| Association with purity score (obtained from tidyestimate)
Association_function(MEs, datTraits$purity, nSamples,nModules,bwModuleColors, "Purity_score_tidyestimate", results.file)

################################################################################


################################################################################
#| QUANTIFYING MODULE-TRAIT ASSOCIATIONS
################################################################################

#| Module trait correlation of every variable
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#| Plotting Module-Trait relationships
sizeGrWindow(8,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 15, 2, 0.5), family ="sans");

#| Variables of interest
datTraits_2 <- subset(datTraits, select=c(Age, DV200.value, purity, stromal,immune, H_score,`PTEN_status`,PTEN_Exp_log2, Mean.expression.PI3K.AKT.mTOR,Macrophages ))
colnames(datTraits_2) <- c("Age", "DV200", "Purity", "Stromal", "Immune", "H-score","PTEN protein","PTEN mRNA", "PI3K-AKT-mTOR", "Macrophages" )
moduleTraitCor_2 <- cor(MEs, datTraits_2, use = "p")
moduleTraitPvalue_2 <- corPvalueStudent(moduleTraitCor_2, nSamples)
textMatrix_2 = paste(signif(moduleTraitCor_2, 2), "\n(",
                     signif(moduleTraitPvalue_2, 1), ")", sep = "")
sizeGrWindow(12,11)
pdf(paste(results.file, "3_Module-trait relationships.pdf", sep =""))
labeledHeatmap(Matrix = moduleTraitCor_2,
               xLabels = names(datTraits_2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               textMatrix = textMatrix_2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               cex.lab.y= 0.45,
               cex.lab.x= 0.4,
               cex.lab= 0.55,
               colors  = colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
               main = paste("Module-trait relationships"))
dev.off()

#| Eigengene heatmap
pdf(paste(results.file,"3_Eigengene_heatmap.pdf", sep = ""))
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90,
                      heatmapColors =colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                      letterSubPlots=T ) 
dev.off()


colnames(MEs) <- gsub("ME", "", colnames(MEs))
sampleDists <- dist(t(MEs), method = "euclidean")/max(dist(t(MEs), method = "euclidean"))


pdf(paste(results.file,"3_Eigengene_heatmap_colors_name.pdf", sep = ""), height = 3.1, width = 4.2)
pheatmap(sampleDistMatrix,
         color = colorRampPalette(c("midnightblue","papayawhip" ,"violetred4"), alpha = TRUE)(100),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows = F)
dev.off()


sampleDists <- cor(MEs)
sampleDistMatrix <- as.matrix(sampleDists)
pdf(paste(results.file,"3_Eigengene_heatmap_colors_name_corrplot.pdf", sep = ""), height = 4.5, width = 5)
corrplot(sampleDistMatrix, 
         type ="upper", 
         tl.col = "black",
         col = colorRampPalette(c("midnightblue","papayawhip" ,"violetred4"), alpha = TRUE)(100))
dev.off()
################################################################################


################################################################################
#| GENE SIGNIFICANCE AND MODULE MEMBERSHIP
################################################################################

#|       Module Membership

#| Name of the color modules
modNames <- substring(names(MEs), 3)

#| Gene Module Membership and MM p-value
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use ="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#|       Gene significance

#| Name of the traits
traitNames <- names(datTraits)

#| Gene Trait significance
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

#| Saving the information
geneInfo <- data.frame(t(datExpr),geneModuleMembership, 
                       geneTraitSignificance, moduleColors, 
                       row.names=colnames(datExpr))

#| Summarized table 
write.table(geneInfo, paste(table.file, "3_geneInfo_counts_geneMM_geneS_moduleColors_last_try.txt", sep =""), row.names = T, sep = "\t")
################################################################################


################################################################################
#| GENE SIGNIFICANCE AND MODULE MEMBERSHIP FOR EVERY MODULE AND EVERY TRAIT
################################################################################

#| Creating a function for Gene Significance and Module membership association
GS_MM <- function(module, geneInfo, GS.trait, label, result.file){
  
  info <- geneInfo
  colnames(info)[which(colnames(info) == paste("MM",module,sep=""))] <- "MM"
  colnames(info)[which(colnames(info) == GS.trait)] <- "GS"
  
  correlation <- cor.test(info[which(info$moduleColors ==module), "GS"], info[which(info$moduleColors ==module), "MM"])
  r <- correlation$estimate[[1]]
  p <- correlation$p.value[[1]]
  
  ggplot(info[which(info$moduleColors ==module), ], aes(x= GS, y = MM)) +
    geom_point( size=3.4, color=module, alpha =0.5) +
    geom_smooth(method="lm", color ="black")+
    theme(text=element_text(size=16,  family="serif"), legend.key.size = unit(0.6, 'cm'),  plot.title=element_text(size=16, hjust = 0.5, face ="bold"))+
    xlab(paste("Gene significance for ", label, sep =""))+
    ylab(paste("Module membership in ", module, sep ="")) +
    ggtitle(paste("Module Membership vs. Gene Significance \n", "cor=", round(r,2)," p=", format(p, scientific = TRUE,digits = 2)))
  ggsave(paste(result.file,"/3_GS_MM_",label,"_",module,".pdf", sep =""), height = 5.5, width = 6)
  
}

#| Plots for every module and every trait
for (i in 1:length(names(geneTraitSignificance))){
  
  GS_trait <- names(geneTraitSignificance)[i]
  GS_trait_label <- gsub("GS.","", GS_trait)
  GS_trait_label <- gsub("\\.","-", GS_trait_label)
  
  for (j in 1:length(modNames)){
    
    dir <- paste(results.file, "GS_MM_plots/", modNames[j], sep = "")
    
    dir.create(dir)
    
    GS_MM(modNames[j], geneInfo,GS_trait, GS_trait_label, dir)
  }
}

################################################################################
