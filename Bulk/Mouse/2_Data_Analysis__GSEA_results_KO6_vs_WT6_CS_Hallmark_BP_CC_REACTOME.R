################################################################################
########### GSEA RESULTS VISUALIZATION PTEN PROTEIN LOSS VS INTACT #############
################################################################################

#| This script sumarizes the results from GSEA over different databases when comparing
#| PTEN protein loss vs intact over 147 PCa patients from Basurto's cohort

################################################################################


###############################  LIBRARIES  ####################################
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(WGCNA))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))

#| Setwd
setwd("C:/Users/irondon/Desktop/AC-12_RNAseq/05_ENRICHMENT/GSEA/05_KO6_vs_WT6/Results/")

#| For plots
theme_set(theme_classic()) 
################################################################################



#################################  DATA  #######################################
GSEA_HALLMARK <- read.table("HALLMARK.Gsea.1705485804450/gsea_report_for_KO6_1705485804450.tsv", sep = "\t", header = T)
GSEA_BP <- read.table("GOBP.Gsea.1705485950544/gsea_report_for_KO6_1705485950544.tsv", sep = "\t", header = T)
GSEA_CC <- read.table("GOCC.Gsea.1705485964551/gsea_report_for_KO6_1705485964551.tsv", sep = "\t", header = T)
GSEA_MF <- read.table("GOMF.Gsea.1705485996780/gsea_report_for_KO6_1705485996780.tsv", sep = "\t", header = T)
GSEA_REACTOME <- read.table("REACTOME.Gsea.1705485931231/gsea_report_for_KO6_1705485931231.tsv", sep = "\t",  header = T)

GSEA_HALLMARK$NAME <- gsub("HALLMARK_","",GSEA_HALLMARK$NAME)
GSEA_BP$NAME <- gsub("GOBP_","",GSEA_BP$NAME)
GSEA_CC$NAME <- gsub("GOCC_","",GSEA_CC$NAME)
GSEA_MF$NAME <- gsub("GOMF_","",GSEA_MF$NAME)
GSEA_REACTOME$NAME <- gsub("REACTOME_","",GSEA_REACTOME$NAME)

GSEA_HALLMARK$NAME <- gsub("_"," ",GSEA_HALLMARK$NAME)
GSEA_BP$NAME <- gsub("_"," ",GSEA_BP$NAME)
GSEA_CC$NAME <- gsub("_"," ",GSEA_CC$NAME)
GSEA_MF$NAME <- gsub("_"," ",GSEA_MF$NAME)
GSEA_REACTOME$NAME <- gsub("_"," ",GSEA_REACTOME$NAME)

GSEA_HALLMARK$source <- "Hallmark"
GSEA_BP$source <- "GO:BP"
GSEA_CC$source <- "GO:CC"
GSEA_MF$source <- "GO:MF"
GSEA_REACTOME$source <- "REACTOME"

#| Save tables
writexl::write_xlsx(GSEA_HALLMARK, "Hallmark.xlsx")
writexl::write_xlsx(GSEA_BP, "GO-BP.xlsx")
writexl::write_xlsx(GSEA_CC, "GO-CC.xlsx")
writexl::write_xlsx(GSEA_MF, "GO-MF.xlsx")
writexl::write_xlsx(GSEA_REACTOME, "REACTOME.xlsx")
################################################################################



########################## HEAT MAP BY DATASETTS ###############################
vector <- c("HALLMARK", "BP", "CC", "MF", "REACTOME")

n <- 3
data_plot <- data.frame(value = c(rep("PTEN loss", times =n), rep("PTEN loss", times =n), rep("PTEN loss", times =n), rep("PTEN loss", times =n), rep("PTEN loss", times =n)),
                             Name = c(GSEA_HALLMARK$NAME[1:n],GSEA_BP$NAME[1:n], GSEA_CC$NAME[1:n], GSEA_MF$NAME[1:n],REACTOME = GSEA_REACTOME$NAME[1:n]),
                             NES = c(GSEA_HALLMARK$NES[1:n],GSEA_BP$NES[1:n], GSEA_CC$NES[1:n],GSEA_MF$NES[1:n],GSEA_REACTOME$NES[1:n]),
                             dataset = c(rep("HALLMARK", times =n), rep("BP", times =n), rep("CC", times =n), rep("MF", times =n), rep("REACTOME", times =n)))

ggplot(data_plot,aes(reorder(dataset, NES), reorder(Name, NES), color= NES)) + 
  geom_point(size=8) +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(1, 'cm'), 
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"),
        axis.text.y = element_text(size = 8, face ="bold", color ="black"))+
  scale_color_viridis(option="C" ) +
  ylab("Term name") + 
  xlab("Datasets") +
  grids(linetype = "dashed")
ggsave("HALLMARK_BP_CC_MF_REACTOME_dotplot_KO6_WT6.pdf", height = 4.5, width = 6.8)

ggplot(GSEA_HALLMARK[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), 
        legend.key.size = unit(0.2, 'cm'), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C" ) +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("HALLMARK_dotplot_Non_phenotype.pdf", height = 3.11, width = 3.4)

ggplot(GSEA_REACTOME[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C" ) +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("REACTOME_dotplot_Non_phenotype.pdf", height = 3.11, width = 4.7)

ggplot(GSEA_BP[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C" )+
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;BP_dotplot_Non_phenotype.pdf", height = 3.11, width = 6.3)

ggplot(GSEA_CC[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C" ) +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;CC_dotplot_Non_phenotype.pdf", height = 3.11, width = 3.6)

ggplot(GSEA_MF[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="sans"), legend.key.size = unit(0.2, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_viridis(option="C" ) +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;MF_dotplot_Non_phenotype.pdf", height = 3.11, width = 4.8)
################################################################################


#################################  DATA  #######################################
GSEA_HALLMARK <- read.table("HALLMARK.Gsea.1691054162039/gsea_report_for_Responders_1691054162039.tsv", sep = "\t", header = T)
GSEA_BP <- read.table("GO.BP.Gsea.1691054356164/gsea_report_for_Responders_1691054356164.tsv", sep = "\t", header = T)
GSEA_CC <- read.table("GO.CC.Gsea.1691054380670/gsea_report_for_Responders_1691054380670.tsv", sep = "\t", header = T)
GSEA_MF <- read.table("GO.MF.Gsea.1691054432662/gsea_report_for_Responders_1691054432662.tsv", sep = "\t", header = T)
GSEA_REACTOME <- read.table("REACTOME.Gsea.1691054333300/gsea_report_for_Responders_1691054333300.tsv", sep = "\t",  header = T)

GSEA_HALLMARK$NAME <- gsub("HALLMARK_","",GSEA_HALLMARK$NAME)
GSEA_BP$NAME <- gsub("GOBP_","",GSEA_BP$NAME)
GSEA_CC$NAME <- gsub("GOCC_","",GSEA_CC$NAME)
GSEA_MF$NAME <- gsub("GOMF_","",GSEA_MF$NAME)
GSEA_REACTOME$NAME <- gsub("REACTOME_","",GSEA_REACTOME$NAME)

GSEA_HALLMARK$NAME <- gsub("_"," ",GSEA_HALLMARK$NAME)
GSEA_BP$NAME <- gsub("_"," ",GSEA_BP$NAME)
GSEA_CC$NAME <- gsub("_"," ",GSEA_CC$NAME)
GSEA_MF$NAME <- gsub("_"," ",GSEA_MF$NAME)
GSEA_REACTOME$NAME <- gsub("_"," ",GSEA_REACTOME$NAME)

GSEA_HALLMARK$source <- "Hallmark"
GSEA_BP$source <- "GO:BP"
GSEA_CC$source <- "GO:CC"
GSEA_MF$source <- "GO:MF"
GSEA_REACTOME$source <- "REACTOME"
################################################################################



########################## HEAT MAP BY DATASETTS ###############################
vector <- c("HALLMARK", "BP", "CC", "MF", "REACTOME")

data_plot <- data.frame(value = c(rep("PTEN loss", times =5), rep("PTEN loss", times =5), rep("PTEN loss", times =5), rep("PTEN loss", times =5), rep("PTEN loss", times =5)),
                        Name = c(GSEA_HALLMARK$NAME[1:5],GSEA_BP$NAME[1:5], GSEA_CC$NAME[1:5], GSEA_MF$NAME[1:5],REACTOME = GSEA_REACTOME$NAME[1:5]),
                        NES = c(GSEA_HALLMARK$NES[1:5],GSEA_BP$NES[1:5], GSEA_CC$NES[1:5],GSEA_MF$NES[1:5],GSEA_REACTOME$NES[1:5]),
                        dataset = c(rep("HALLMARK", times =5), rep("BP", times =5), rep("CC", times =5), rep("MF", times =5), rep("REACTOME", times =5)))

ggplot(data_plot,aes(reorder(dataset, NES), reorder(Name, NES), color= NES)) + 
  geom_point(size=8) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(1, 'cm'), plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  scale_color_viridis(option="C" ) +
  ylab("Term name") + 
  xlab("Datasets") +
  grids(linetype = "dashed")
ggsave("HALLMARK_BP_CC_MF_REACTOME_dotplot_Responders_phenotype.pdf", height = 7, width = 8.6)


ggplot(GSEA_HALLMARK[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("HALLMARK_dotplot_Responders_phenotype.pdf", height = 4, width = 3.8)

GSEA_HALLMARK$HALLMARK <- "HALLMARK"
ggplot(GSEA_HALLMARK[1:10,], aes(x= HALLMARK, y= reorder(NAME, NES)) ) +
  geom_point( aes( size =NES), fill = rainbow(5)[1], shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'), 
        axis.text.y = element_text(size = 8, color ="black"),
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  ylab("Datasets") +
  xlab("")
ggsave("HALLMARK_bubbleplot.pdf", heigh=3.5, width = 3.8)



ggplot(GSEA_REACTOME[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("REACTOME_dotplot_Responders_phenotype.pdf", height = 4, width = 6.8)

GSEA_REACTOME$REAC <- "REACTOME"
ggplot(GSEA_REACTOME[1:10,], aes(x= REAC, y= reorder(NAME, NES)) ) +
  geom_point( aes( size =NES), fill = rainbow(5)[4], shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'), 
        axis.text.y = element_text(size = 8, color ="black"),
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  ylab("Datasets") +
  xlab("")
ggsave("REACTOME_bubbleplot.pdf", heigh=3.5, width = 5)

ggplot(GSEA_BP[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;BP_dotplot_Responders_phenotype.pdf", height = 4, width = 6)


GSEA_BP$BP <- "GO:BP"
ggplot(GSEA_BP[1:10,], aes(x= BP, y= reorder(NAME, NES)) ) +
  geom_point( aes( size =NES), fill = rainbow(5)[2], shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'), 
        axis.text.y = element_text(size = 8, color ="black"),
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  ylab("Datasets") +
  xlab("")
ggsave("GO-BP_bubbleplot.pdf", heigh=3.5, width = 6.5)



ggplot(GSEA_CC[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;CC_dotplot_Responders_phenotype.pdf", height = 4, width = 4.6)

GSEA_CC$CC <- "GO:CC"
ggplot(GSEA_CC[1:10,], aes(x= CC, y= reorder(NAME, NES)) ) +
  geom_point( aes( size =NES), fill = rainbow(5)[3], shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'), 
        axis.text.y = element_text(size = 8, color ="black"),
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  ylab("Datasets") +
  xlab("")
ggsave("GO-CC_bubbleplot.pdf", heigh=3.5, width = 3.9)


ggplot(GSEA_MF[1:10,], aes(x = reorder(source, NES), y = reorder(NAME, NES))) + 
  geom_point(aes(size = NES, fill = NOM.p.val), alpha = 0.75, shape = 21) +
  theme(text=element_text(size=10,  family="serif"), legend.key.size = unit(0.5, 'cm'), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill= "p value") +
  labs(color ="Module size") +
  labs(size= "NES")+
  xlab("Data") + 
  ylab("Term name")
ggsave("GO;MF_dotplot_Responders_phenotype.pdf", height = 4, width = 5.5)

GSEA_MF$MF <- "GO:MF"
ggplot(GSEA_MF[1:10,], aes(x= MF, y= reorder(NAME, NES)) ) +
  geom_point( aes( size =NES), fill = rainbow(5)[5], shape = 21) +
  theme(text=element_text(size=12,  family="sans"), 
        legend.key.size = unit(0.8, 'cm'), 
        axis.text.y = element_text(size = 8, color ="black"),
        plot.title=element_text(size=13, hjust = 0.5, face ="bold"))+
  ylab("Datasets") +
  xlab("")
ggsave("GO-MF_bubbleplot.pdf", heigh=3.5, width = 5.12)

################################################################################