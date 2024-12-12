#| Last change: 12/12/2024
#| Ivana Rondon-Lorefice

################################################################################
#########  ENRICHMENT ANALYSIS FOR THE COMPARISONS FROM DESEQ2 RESULTS  ########
################################################################################



################################  LIBRARIES  ###################################
suppressMessages(library(dplyr))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))

theme_set(theme_classic())
################################################################################


################################  DATA  ###################################
#| Data directories 
dir.proj <- "X:/irondon/AC-12_RNAseq/05_ENRICHMENT/gprofiler/"
counts.file <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt"
setwd(dir.proj)

#| Counts data
counts_data <- read.table(counts.file, sep ="\t")
counts_data <- counts_data[,order(colnames(counts_data))]

#| Sample info
sample <- colnames(counts_data)
condition <- gsub("_[0-9]*","",sample)
month <- gsub("[A-Z]", "", condition) 
sample_info <- data.frame(sample = sample,
                          condition = condition,
                          month = month)

#| Load the genome data
genome_mouse<- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Data/geneInfo.txt", sep ="\t", header =F)
colnames(genome_mouse) <- c("GeneID", "gene_name", "gene_info")

#| Results from the differential expression analysis
condition_KO3_vs_WT3 <- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Results/Tables/res_KO3_WT3_dataframe.txt", sep ="\t", header = T)
condition_KO6_vs_WT6 <- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Results/Tables/res_KO6_WT6_dataframe.txt", sep ="\t", header = T)
condition_KO6_vs_KO3 <- read.table("X:/irondon/AC-12_RNAseq/04_DEGs/Results/Tables/res_KO6_KO3_dataframe.txt", sep ="\t", header = T)

FC <- 2
FDR <- 0.05

#| DEGs 
condition_KO3_vs_WT3_DEGs <- condition_KO3_vs_WT3$GeneID[which( (condition_KO3_vs_WT3$padj_deseq2 < FDR) & (!is.na(condition_KO3_vs_WT3$padj_deseq2)) & ( (condition_KO3_vs_WT3$log2FC_deseq2 < (-log2(FC))) | (condition_KO3_vs_WT3$log2FC_deseq2 > log2(FC)) ) )]
condition_KO6_vs_WT6_DEGs <- condition_KO6_vs_WT6$GeneID[which( (condition_KO6_vs_WT6$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_WT6$padj_deseq2)) & ( (condition_KO6_vs_WT6$log2FC_deseq2 < (-log2(FC))) | (condition_KO6_vs_WT6$log2FC_deseq2 > log2(FC)) ) )]
condition_KO6_vs_KO3_DEGs <- condition_KO6_vs_KO3$GeneID[which( (condition_KO6_vs_KO3$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_KO3$padj_deseq2)) & ( (condition_KO6_vs_KO3$log2FC_deseq2 < (-log2(1))) | (condition_KO6_vs_KO3$log2FC_deseq2 > log2(1)) ) )]

#| DEGs UP
condition_KO3_vs_WT3_DEGs_UP <- condition_KO3_vs_WT3$GeneID[which( (condition_KO3_vs_WT3$padj_deseq2 < FDR) & (!is.na(condition_KO3_vs_WT3$padj_deseq2)) & ( (condition_KO3_vs_WT3$log2FC_deseq2 > log2(FC)) ) )]
condition_KO6_vs_WT6_DEGs_UP <- condition_KO6_vs_WT6$GeneID[which( (condition_KO6_vs_WT6$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_WT6$padj_deseq2)) & ( (condition_KO6_vs_WT6$log2FC_deseq2 > log2(FC)) ) )]
condition_KO6_vs_KO3_DEGs_UP <- condition_KO6_vs_KO3$GeneID[which( (condition_KO6_vs_KO3$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_KO3$padj_deseq2)) & ( (condition_KO6_vs_KO3$log2FC_deseq2 > log2(1)) ) )]

#| DEGs DOWN
condition_KO3_vs_WT3_DEGs_DOWN <- condition_KO3_vs_WT3$GeneID[which( (condition_KO3_vs_WT3$padj_deseq2 < FDR) & (!is.na(condition_KO3_vs_WT3$padj_deseq2)) & ( (condition_KO3_vs_WT3$log2FC_deseq2 < (-log2(FC))) ) )]
condition_KO6_vs_WT6_DEGs_DOWN <- condition_KO6_vs_WT6$GeneID[which( (condition_KO6_vs_WT6$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_WT6$padj_deseq2)) & ( (condition_KO6_vs_WT6$log2FC_deseq2 < (-log2(FC))) ) )]
condition_KO6_vs_KO3_DEGs_DOWN <- condition_KO6_vs_KO3$GeneID[which( (condition_KO6_vs_KO3$padj_deseq2 < FDR) & (!is.na(condition_KO6_vs_KO3$padj_deseq2)) & ( (condition_KO6_vs_KO3$log2FC_deseq2 < (-log2(1))) ) )]
################################################################################


################################################################################
#| ENRICHEMENT WITH GPROFILER
################################################################################

################################# KO3 vs WT3 #####################################
gost_KO3_vs_WT3 <- gost(list("Condition KO3 vs WT3" = condition_KO3_vs_WT3_DEGs), 
                  organism = "mmusculus", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = TRUE,
                  user_threshold = 0.05, correction_method = "fdr",
                  domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO3_vs_WT3, interactive = FALSE, capped =FALSE)
ggsave("Results/KO3_vs_WT3/Images/ALL_DEGs/ggplot_results_KO3_vs_WT3.pdf", height = 4, width = 7)

results_KO3_vs_WT3 <- gost_KO3_vs_WT3$result[order(gost_KO3_vs_WT3$result$p_value),]
results_KO3_vs_WT3$`Term name` <- paste(results_KO3_vs_WT3$term_name, "\n (N = ",results_KO3_vs_WT3$term_size, ")",sep ="")

write.table(results_KO3_vs_WT3[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO3_vs_WT3/Tables/ggplot_results_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO3_vs_WT3[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/TOP_15_ggplot_results_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 5, width = 5)


#| Top with GO:CC
dat1_filtered <- results_KO3_vs_WT3[which(results_KO3_vs_WT3$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/GO-CC_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with GO:BP
dat1_filtered <- results_KO3_vs_WT3[which(results_KO3_vs_WT3$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/GO-BP_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with GO:MF
dat1_filtered <- results_KO3_vs_WT3[which(results_KO3_vs_WT3$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/GO-MF_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with REAC
dat1_filtered <- results_KO3_vs_WT3[which(results_KO3_vs_WT3$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/REACTOME_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 6)

#| Top with KEGG
dat1_filtered <- results_KO3_vs_WT3[which(results_KO3_vs_WT3$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/ALL_DEGs/KEGG_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 5)


########| UP
gost_KO3_vs_WT3_UP <- gost(list("Condition KO3 vs KT3 UP" = condition_KO3_vs_WT3_DEGs_UP), 
                        organism = "mmusculus", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO3_vs_WT3_UP, interactive = FALSE, capped =FALSE)
ggsave("Results/KO3_vs_WT3/Images/UP_DEGs/ggplot_results_KO3_vs_WT3.pdf", height = 4, width = 7)

results_KO3_vs_WT3_UP <- gost_KO3_vs_WT3_UP$result[order(gost_KO3_vs_WT3_UP$result$p_value),]
results_KO3_vs_WT3_UP$`Term name` <- paste(results_KO3_vs_WT3_UP$term_name, "\n (N = ",results_KO3_vs_WT3_UP$term_size, ")",sep ="")

write.table(results_KO3_vs_WT3_UP[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO3_vs_WT3/Tables/ggplot_results_UP_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO3_vs_WT3_UP[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/TOP_15_ggplot_results_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 5, width = 5.5)


#| Top with GO:CC
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/GO-CC_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with GO:BP
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/GO-BP_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 5.5)

#| Top with GO:MF
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/GO-MF_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with REAC
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/REACTOME_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 6)

#| Top with KEGG
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/UP_DEGs/KEGG_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 5)


########| DOWN

gost_KO3_vs_WT3_DOWN <- gost(list("Condition KO3 vs KT3 DOWN" = condition_KO3_vs_WT3_DEGs_DOWN), 
                           organism = "mmusculus", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO3_vs_WT3_DOWN, interactive = FALSE, capped =FALSE)
ggsave("Results/KO3_vs_WT3/Images/DOWN_DEGs/ggplot_results_KO3_vs_WT3.pdf", height = 4, width = 7)

results_KO3_vs_WT3_DOWN <- gost_KO3_vs_WT3_DOWN$result[order(gost_KO3_vs_WT3_DOWN$result$p_value),]
results_KO3_vs_WT3_DOWN$`Term name` <- paste(results_KO3_vs_WT3_DOWN$term_name, "\n (N = ",results_KO3_vs_WT3_DOWN$term_size, ")",sep ="")

write.table(results_KO3_vs_WT3_DOWN[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO3_vs_WT3/Tables/ggplot_results_DOWN_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO3_vs_WT3_DOWN[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO3_vs_WT3/Images/DOWN_DEGs/TOP_10_ggplot_results_KO3_vs_WT3_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 8)

################################################################################


################################# KO6 vs WT6 #####################################
gost_KO6_vs_WT6 <- gost(list("Condition KO6 vs WT6" = condition_KO6_vs_WT6_DEGs), 
                        organism = "mmusculus", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE,
                        user_threshold = 0.05, correction_method = "fdr",
                        domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO6_vs_WT6, interactive = FALSE, capped =FALSE)
ggsave("Results/KO6_vs_WT6/Images/ALL_DEGs/ggplot_results_KO6_vs_WT6.pdf", height = 4, width = 7)

results_KO6_vs_WT6 <- gost_KO6_vs_WT6$result[order(gost_KO6_vs_WT6$result$p_value),]
results_KO6_vs_WT6$`Term name` <- paste(results_KO6_vs_WT6$term_name, "\n (N = ",results_KO6_vs_WT6$term_size, ")",sep ="")

write.table(results_KO6_vs_WT6[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO6_vs_WT6/Tables/ggplot_results_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO6_vs_WT6[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/TOP_15_ggplot_results_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 5, width = 5)


#| Top with GO:CC
dat1_filtered <- results_KO6_vs_WT6[which(results_KO6_vs_WT6$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/GO-CC_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with GO:BP
dat1_filtered <- results_KO6_vs_WT6[which(results_KO6_vs_WT6$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/GO-BP_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with GO:MF
dat1_filtered <- results_KO6_vs_WT6[which(results_KO6_vs_WT6$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/GO-MF_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with REAC
dat1_filtered <- results_KO6_vs_WT6[which(results_KO6_vs_WT6$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/REACTOME_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 6)

#| Top with KEGG
dat1_filtered <- results_KO6_vs_WT6[which(results_KO6_vs_WT6$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/ALL_DEGs/KEGG_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 5)


########| UP
gost_KO6_vs_WT6_UP <- gost(list("Condition KO6 vs KT6 UP" = condition_KO6_vs_WT6_DEGs_UP), 
                           organism = "mmusculus", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO6_vs_WT6_UP, interactive = FALSE, capped =FALSE)
ggsave("Results/KO6_vs_WT6/Images/UP_DEGs/ggplot_results_KO6_vs_WT6.pdf", height = 4, width = 7)

results_KO6_vs_WT6_UP <- gost_KO6_vs_WT6_UP$result[order(gost_KO6_vs_WT6_UP$result$p_value),]
results_KO6_vs_WT6_UP$`Term name` <- paste(results_KO6_vs_WT6_UP$term_name, "\n (N = ",results_KO6_vs_WT6_UP$term_size, ")",sep ="")

write.table(results_KO6_vs_WT6_UP[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO6_vs_WT6/Tables/ggplot_results_UP_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO6_vs_WT6_UP[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/TOP_15_ggplot_results_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 5, width = 5.5)


#| Top with GO:CC
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "GO:CC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/GO-CC_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 3.8)

#| Top with GO:BP
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "GO:BP"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/GO-BP_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 3.8)

#| Top with GO:MF
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "GO:MF"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/GO-MF_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 4)

#| Top with REAC
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/REACTOME_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 8.1)

#| Top with KEGG
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "KEGG"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/UP_DEGs/KEGG_ggplot_results_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 3.5, width = 5)


########| DOWN

gost_KO6_vs_WT6_DOWN <- gost(list("Condition KO6 vs WT6 DOWN" = condition_KO6_vs_WT6_DEGs_DOWN), 
                             organism = "mmusculus", ordered_query = FALSE, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = TRUE,
                             user_threshold = 0.05, correction_method = "fdr",
                             domain_scope = "custom", custom_bg = condition_KO3_vs_WT3$GeneID, 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO6_vs_WT6_DOWN, interactive = FALSE, capped =FALSE)
ggsave("Results/KO6_vs_WT6/Images/DOWN_DEGs/ggplot_results_KO6_vs_WT6.pdf", height = 4, width = 7)

results_KO6_vs_WT6_DOWN <- gost_KO6_vs_WT6_DOWN$result[order(gost_KO6_vs_WT6_DOWN$result$p_value),]
results_KO6_vs_WT6_DOWN$`Term name` <- paste(results_KO6_vs_WT6_DOWN$term_name, "\n (N = ",results_KO6_vs_WT6_DOWN$term_size, ")",sep ="")

write.table(results_KO6_vs_WT6_DOWN[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection", "Term name")], paste("Results/KO6_vs_WT6/Tables/ggplot_results_DOWN_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".txt", sep =""), sep = "\t", row.names = T)

#| Top 15 most significant
ggplot(results_KO6_vs_WT6_DOWN[1:2,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
ggsave(paste("Results/KO6_vs_WT6/Images/DOWN_DEGs/TOP_2_ggplot_results_KO6_vs_WT6_FC-",FC,"_FDR-",FDR,".pdf", sep = ""), height = 2.5, width = 6)

################################################################################


################################################################################



################################################################################
#| ENRICHMENT OF DEGS THAT DOES NOT OVERLAPPED
################################################################################


res_KO3_WT3_dataframe <- condition_KO3_vs_WT3
res_KO6_WT6_dataframe <- condition_KO6_vs_WT6

#| KO3 vs WT3
DEGs_KO3_vs_WT3 <- res_KO3_WT3_dataframe$gene_name[which((res_KO3_WT3_dataframe$log2FC_deseq2 < (-log2(FC)) | (res_KO3_WT3_dataframe$log2FC_deseq2 > (log2(FC)))) & (res_KO3_WT3_dataframe$padj_deseq2 < FDR))] 
DEGs_UP_KO3_vs_WT3 <- res_KO3_WT3_dataframe$gene_name[which((res_KO3_WT3_dataframe$log2FC_deseq2 > log2(FC)) & (res_KO3_WT3_dataframe$padj_deseq2 < FDR))] 
DEGs_DOWN_KO3_vs_WT3 <- res_KO3_WT3_dataframe$gene_name[which((res_KO3_WT3_dataframe$log2FC_deseq2 < (-log2(FC))) & (res_KO3_WT3_dataframe$padj_deseq2 < FDR))] 

#| KO6 vs WT6
DEGs_KO6_vs_WT6 <- res_KO6_WT6_dataframe$gene_name[which((res_KO6_WT6_dataframe$log2FC_deseq2 < (-log2(FC)) | (res_KO6_WT6_dataframe$log2FC_deseq2 > (log2(FC)))) & (res_KO6_WT6_dataframe$padj_deseq2 < FDR))] 
DEGs_UP_KO6_vs_WT6 <- res_KO6_WT6_dataframe$gene_name[which((res_KO6_WT6_dataframe$log2FC_deseq2 > log2(FC)) & (res_KO6_WT6_dataframe$padj_deseq2 < FDR))] 
DEGs_DOWN_KO6_vs_WT6 <- res_KO6_WT6_dataframe$gene_name[which((res_KO6_WT6_dataframe$log2FC_deseq2 < (-log2(FC))) & (res_KO6_WT6_dataframe$padj_deseq2 < FDR))] 


#| UP DEGS KO3 vs WT3

gost_KO3_vs_WT3_UP <- gost(list("Condition KO3 vs KT3 UP no in KO6 vs WT6" = DEGs_UP_KO3_vs_WT3[which( !(DEGs_UP_KO3_vs_WT3 %in% DEGs_UP_KO6_vs_WT6 ))]), 
                           organism = "mmusculus", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = res_KO3_WT3_dataframe$gene_name, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO3_vs_WT3_UP, interactive = FALSE, capped =FALSE)

results_KO3_vs_WT3_UP <- gost_KO3_vs_WT3_UP$result[order(gost_KO3_vs_WT3_UP$result$p_value),]
results_KO3_vs_WT3_UP$`Term name` <- paste(results_KO3_vs_WT3_UP$term_name, "\n (N = ",results_KO3_vs_WT3_UP$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_KO3_vs_WT3_UP[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")

#| Top with REAC
dat1_filtered <- results_KO3_vs_WT3_UP[which(results_KO3_vs_WT3_UP$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")

#| UP DEGS KO3 vs WT3
gost_KO3_vs_WT3_DOWN <- gost(list("Condition KO3 vs KT3 DOWN" = DEGs_DOWN_KO3_vs_WT3[which( !(DEGs_DOWN_KO3_vs_WT3 %in% DEGs_DOWN_KO6_vs_WT6 ))]), 
                             organism = "mmusculus", ordered_query = FALSE, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = TRUE,
                             user_threshold = 0.05, correction_method = "fdr",
                             domain_scope = "custom", custom_bg = res_KO3_WT3_dataframe$gene_name, 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO3_vs_WT3_DOWN, interactive = FALSE, capped =FALSE)

results_KO3_vs_WT3_DOWN <- gost_KO3_vs_WT3_DOWN$result[order(gost_KO3_vs_WT3_DOWN$result$p_value),]
results_KO3_vs_WT3_DOWN$`Term name` <- paste(results_KO3_vs_WT3_DOWN$term_name, "\n (N = ",results_KO3_vs_WT3_DOWN$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_KO3_vs_WT3_DOWN[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")



#| UP DEGs KO6 vs WT
gost_KO6_vs_WT6_UP <- gost(list("Condition KO6 vs WT6" = DEGs_UP_KO6_vs_WT6[which( !(DEGs_UP_KO6_vs_WT6 %in% DEGs_UP_KO3_vs_WT3 ))]), 
                           organism = "mmusculus", ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = TRUE,
                           user_threshold = 0.05, correction_method = "fdr",
                           domain_scope = "custom", custom_bg = res_KO6_WT6_dataframe$gene_name, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gost_KO6_vs_WT6_UP, interactive = FALSE, capped =FALSE)

results_KO6_vs_WT6_UP <- gost_KO6_vs_WT6_UP$result[order(gost_KO6_vs_WT6_UP$result$p_value),]
results_KO6_vs_WT6_UP$`Term name` <- paste(results_KO6_vs_WT6_UP$term_name, "\n (N = ",results_KO6_vs_WT6_UP$term_size, ")",sep ="")

#| Top 15 most significant
ggplot(results_KO6_vs_WT6_UP[1:15,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  #ggtitle("Enrichment Module Green") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")

#| Top with REAC
dat1_filtered <- results_KO6_vs_WT6_UP[which(results_KO6_vs_WT6_UP$source == "REAC"),]
ggplot(dat1_filtered[1:10,], aes(x = reorder(source, -p_value), y = reorder(`Term name`, -p_value))) + 
  geom_point(aes(size = intersection_size, fill = p_value), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.4, 'cm'), plot.title=element_text(size=12, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "Intersection size")+
  xlab("Data")+ 
  ylab("Term name")
################################################################################
