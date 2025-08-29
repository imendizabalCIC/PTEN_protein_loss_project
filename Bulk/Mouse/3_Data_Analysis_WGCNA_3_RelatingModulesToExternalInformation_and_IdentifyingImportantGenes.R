################################################################################
#| RELATING MODULES TO EXTERNAL INFORMATION WITH DATA          
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script integrates WGCNA module eigengenes with clinical and molecular traits 
#| from the AC-12_RNAseq mouse cohort (Pten KO6, WT6, KO3, WT3). It quantifies 
#| correlations between modules and external traits, identifies biologically 
#| relevant modules, and explores the relationship between module membership and 
#| gene significance.
#|
#| Workflow:
#|   1) Load WGCNA network results (modules, eigengenes, dendrograms) and trait data.  
#|   2) Recalculate eigengenes and plot dendrograms with module colors.  
#|   3) Correlate module eigengenes with clinical variables (PTEN protein, H-score, 
#|      stromal, immune, purity, etc.).  
#|   4) Visualize results with scatter plots and module-trait heatmaps.  
#|   5) Quantify gene-level associations:  
#|        - Module membership (MM) with respect to eigengenes.  
#|        - Gene significance (GS) with respect to traits.  
#|   6) Plot GS-MM relationships per module and trait to highlight key drivers.  
#|   7) Save summary tables of module membership and gene-trait significance.  
#|
#| Outputs:
#|   - Module-trait association plots (PDF).  
#|   - Heatmaps of module-trait correlations and eigengene networks.  
#|   - GS-MM scatterplots for each module-trait pair.  
#|   - Summarized table: geneInfo with counts, MM, GS, and module colors.  
################################################################################


################################################################################
#|  LIBRARIES AND DATA
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

dir.proj <- "X:/irondon/AC-12_RNAseq/07_WGCNA/"
results.file <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Images/03_Module-trait relationships/"
data.file <- "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Data/"

setwd(dir.proj)

load(file = paste(data.file, "AC-12_RNAseq_dataInput.RData", sep =""))
load(file = paste(data.file, "2_blockwiseModules_AC_12_RNAseq.RData", sep =""))

data_loaded = TRUE
names(datTraits_final)

#| For plots
theme_set(theme_classic())

datTraits <- datTraits_final
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
pdf(paste(results.file,"3_SingleBlock_geneDendogram_and_ModuleColors.pdf", sep =""),height = 4.5, width = 6)
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
  
  # Apply FDR correction
  moduleTraitPvalue_FDR <- apply(moduleTraitPvalue, 2, function(p) p.adjust(p, method = "fdr"))
  
  #| Changing rownames 
  rownames(moduleTraitCor) <- gsub("ME", "", rownames(moduleTraitCor))
  rownames(moduleTraitPvalue_FDR) <- gsub("ME", "", rownames(moduleTraitPvalue_FDR))
  
  if(any(rownames(moduleTraitCor) == rownames(moduleTraitPvalue_FDR))){
    
    size <- c()
    for (i in 1:nModules){
      size <- c(size,length(bwModuleColors[which(bwModuleColors == rownames(moduleTraitCor)[i])])) 
    }
    
    #| Creating data frame with values
    log_pvalue_association <- data.frame(fdr = moduleTraitPvalue_FDR ,
                                         Association = moduleTraitCor,
                                         module_size = size,
                                         Modules = rownames(moduleTraitCor))
    
    #| Changing module size to numeric
    log_pvalue_association$module_size <- as.numeric(log_pvalue_association$module_size)
    
    modules <- rownames(moduleTraitCor)
    modules <- as.character(modules)
    modules <- sort(modules)
    
    #| Plot and save
    ggplot(log_pvalue_association, aes(x = Association, y = -log10(fdr))) + 
      geom_point(aes(size =module_size, fill=Modules), colour="black",pch=21,alpha = 0.8) +
      theme(text=element_text(size=12,  family="sans"), 
            axis.text.x = element_text(size=12, color ="black"),
            axis.text.y = element_text(size=12, color ="black")) +
      scale_fill_manual(values = modules, guide = "legend") +
      labs(size ="Module size") +
      xlab(paste("Association with ",label,sep="")) +
      guides(fill = "none")+
      geom_hline(yintercept=-log10(0.05), linetype='dotted', col = 'grey', linewidth= 0.3)+
      geom_vline(xintercept=0.2, linetype='dotted', col = 'black', linewidth= 0.3)+
      geom_vline(xintercept=-0.2, linetype='dotted', col = 'black', linewidth= 0.3)#+
      #xlim(-0.3,0.3)
    ggsave(paste(results.file,"3_Association_with_", label,"_log10FDR.pdf", sep =""),height = 4.5, width = 4.7)
    
  }else{
    
    print("Check rownames(moduleTraitCor) == rownames(moduleTraitPvalue)")
    
  }
  
}

#| Association with KO6
Association_function(MEs, datTraits$KO6.vs.WT6, nSamples,nModules,bwModuleColors, "KO6", results.file)

#| Association with PTEN protein loss
Association_function(MEs, datTraits$KO3.vs.WT3, nSamples,nModules,bwModuleColors, "KO3", results.file)

################################################################################


################################################################################
#| QUANTIFYING MODULE-TRAIT ASSOCIATIONS
################################################################################

#| Module trait correlation of every variable
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue_FDR <- apply(moduleTraitPvalue, 2, function(p) p.adjust(p, method = "fdr"))
moduleTraitPvalue_FDR <- as.data.frame(moduleTraitPvalue_FDR)

Supplementary_Table4_1 <- as.data.frame(moduleTraitCor[,c("KO6.vs.WT6","KO3.vs.WT3")])
Supplementary_Table4_2 <- as.data.frame(moduleTraitPvalue_FDR[,c("KO6.vs.WT6","KO3.vs.WT3")])

Supplementary_Table4 <- cbind(Supplementary_Table4_1, Supplementary_Table4_2)
Supplementary_Table4$Module_name <- rownames(Supplementary_Table4)
writexl::write_xlsx(Supplementary_Table4, "X:/irondon/AC-12_RNAseq/07_WGCNA/Results/Tables/Supplementary_Table8.xlsx")


#| Plotting Module-Trait relationships
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 15, 2, 0.5), family ="sans");
sizeGrWindow(12,11)
pdf(paste(results.file, "3_Module-trait relationships.pdf", sep =""), height = 4.5, width = 3.9)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.35,
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

#| Heatmap with correlations
colnames(MEs) <- gsub("ME", "", colnames(MEs))
sampleDists <- dist(t(MEs), method = "euclidean")/max(dist(t(MEs), method = "euclidean"))
sampleDistMatrix <- as.matrix(sampleDists)
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
write.table(geneInfo, paste(data.file, "3_geneInfo_counts_geneMM_geneS_moduleColors.txt", sep =""), row.names = T, sep = "\t")
################################################################################


################################################################################
#| GENE SIGNIFICANCE AND MODULE MEMBERSHIP
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


