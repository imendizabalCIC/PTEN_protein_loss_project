################################################################################
#| QUANTIFICATION AND BINARIZATION OF PTEN
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script integrates PTEN immunohistochemistry (IHC) data with RNA-seq 
#| sample information for the AC-45_RNAseq-FFPE cohort. PTEN staining intensities 
#| (H-score) were obtained from a pathologist's evaluation of tumor tissue in 197 
#| patients and merged into the sample_info table.
#|
#| Workflow:
#|   1) Process raw IHC assessment data (TMA matrix + PTEN sheet) and compute 
#|      H-scores using:  
#|        H-score = (0 X %negative) + (1 X %weak) + (2 X %moderate) + (3 X %strong).  
#|   2) Merge H-scores and % tumor values with clinical sample_info metadata.  
#|   3) Explore binary classification thresholds (cut-offs: 0, 10, 20, 30, 200) 
#|      to define PTEN protein status as "presence" vs "intact".  
#|   4) Correlate thresholds with PI3K-AKT-mTOR pathway expression signature 
#|      using Wilcoxon tests.  
#|   5) Generate visualizations:  
#|        - Scatterplots of threshold p-values vs signature expression.  
#|        - Boxplots of PI3K-AKT-mTOR expression across PTEN cut-offs.  
#|        - Barplots of sample proportions across thresholds.  
#|   6) Save updated sample_info with H-scores and PTEN classification.  
#|
#| Outputs:
#|   - Updated sample_info_extracted.txt including PTEN H-scores and thresholds.  
#|   - Boxplots, scatterplots, and barplots summarizing PTEN threshold analysis.  
################################################################################


################################################################################
#| LIBRARIES AND DATA
################################################################################
suppressMessages(library(apeglm))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(rstatix))
suppressMessages(library(viridis))
suppressMessages(library(vsn))
suppressMessages(library(utils))
suppressMessages(library(tidyestimate))
suppressMessages(library(gprofiler2))
suppressMessages(library(ggplot2))
suppressMessages(library(corrr))
suppressMessages(library(readxl))
suppressMessages(library(writexl))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(colorspace))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))

#| Pathways directories
dir.proj <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/"
info.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Sample_info_table/sample_info_extracted.txt"
counts.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Results/Counts/FullCounts_Basurto.txt"
pi3k_akt_mtor.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/PI3K_AKT_MTOR_SIGNALING_PATHWAY_SIGNATURE_GSEA.txt"

#| Setting working directory
setwd(dir.proj)

#| For plots
theme_set(theme_classic()) 

#| Load the genome for merging GeneID. Used the BigData/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt
genome_GRCh39.94 <- read.table("X:/DATA_shared/Genomes/Homo_sapiens.GRCh38.94_geneLength.txt", sep ="\t", header = T)

#| Sample information
sample_info <- read.table(info.file, sep ="\t", header=T)
sample_info$`AC basurto` <- sample_info$AC.basurto

#| Counts data (normalized)
counts_data <- read.table(counts.file, sep ="\t", header=T)
################################################################################


################################################################################
#| PTEN ASSESSMENT
################################################################################
Assessment.file <- "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/4_STEP_DATA_ANALYSIS/4_Data_Basurto/Data/TMA assessment PTEN and pAKT.xlsx"
TMA_matrix_ID <-"Matriz TMA"
PTEN_Assessment.sheet <-"PTEN"

#| Reading Matrix TMA
data_TMA_matrix_ID <- read_excel(Assessment.file, sheet = TMA_matrix_ID)
data_TMA_matrix_ID <- data_TMA_matrix_ID[, 1:8]

#| Reading Assessment data
data_PTEN_Assessment <- read_excel(Assessment.file, sheet =PTEN_Assessment.sheet)

#| Renaming first column
colnames(data_TMA_matrix_ID)[1] <- "TMA ID"

#| Finding TMA IDs
TMA_ID <- data_TMA_matrix_ID$`TMA ID`[grepl("TMA", data_TMA_matrix_ID$`TMA ID`)]
AC_ID <- c()

#| Over all the data
for (i in 1:length(TMA_ID)){
  
  if (i < length(TMA_ID)){
    df <- data_TMA_matrix_ID[(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i])+1):(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i+1])-3),]
  } else {
    df <- data_TMA_matrix_ID[(which(data_TMA_matrix_ID$`TMA ID` == TMA_ID[i])+1):dim(data_TMA_matrix_ID)[1],]
  }
  
  for (j in 0:(dim(df)[2]-1)){
    value <- df[1:dim(df)[1],dim(df)[2]-j]
    AC_ID <- c(AC_ID, value[[1]])
  }
  
}

#| Adding a column with AC ID information to the assessment data
data_PTEN_Assessment$`AC TMA` <- AC_ID

#| Changing AC ID to AC Basurto's ID
data_PTEN_Assessment$`AC basurto` <- gsub("AC-","AC",data_PTEN_Assessment$`AC TMA`)

#| H-score = [(0 x % negative cells) + (1 x %weakly positive cells) + (2 x %moderately positive cells) + (3 x %strongly positive cells)] 
data_PTEN_Assessment$H_score <- (data_PTEN_Assessment$`% negativo`)*0 + (data_PTEN_Assessment$`% débil`)*1 + (data_PTEN_Assessment$`% moderado`)*2 + (data_PTEN_Assessment$`% intenso`)*3
data_PTEN_Assessment$`% tumor numeric` <- as.numeric(data_PTEN_Assessment$`% tumor`)

#| Assigning NA to those value where there was tumor but negative for PTEN in the stroma
data_PTEN_Assessment$H_score[which(data_PTEN_Assessment$estroma == "neg")] <- NA

#| Discarding H-score with NAs
data_PTEN_Assessment2 <- data_PTEN_Assessment[which(!is.na(data_PTEN_Assessment$H_score)),]
H_score_dataframe2 <- aggregate( H_score ~ `AC basurto`, data_PTEN_Assessment2, mean )
H_score_dataframe2$H_score == 0
sample_info$`AC basurto` <- sample_info$AC.basurto

#| Aggregating dataframes
H_score_dataframe <- aggregate( H_score ~ `AC basurto`, data_PTEN_Assessment, mean )
Percentage_tumor <- aggregate( `% tumor numeric` ~ `AC basurto`, data_PTEN_Assessment, mean )

#| Merging dataframes
sample_info <- merge(sample_info, H_score_dataframe, by="AC basurto")
sample_info <- merge(sample_info, Percentage_tumor, by="AC basurto")
sample_info$H_score<- sample_info$H_score.x
rownames(sample_info) <- sample_info$`AC basurto`
################################################################################


################################################################################
#| FINDING H-SCORE THRESHOLD BASED ON SIGNATURES
################################################################################

#| Can we set a threshold of the H-score based on the values of Mean expression FOXO?
sample_info_extracted <- sample_info
threshold <- unique(sample_info_extracted$H_score)
threshold <- threshold[which(!is.na(threshold))]
threshold <- sort(threshold)
p_value_wilcox_PI3K_AKT_mTOR <- c()
difference_wilcox_PI3K_AKT_mTOR <- c()

for (i in 1:(length(threshold)-1)){
  
  #| Creating a new column in the sample_info_extracted
  sample_info_extracted$H_score_threshold <- sample_info_extracted$H_score
  
  #| If values are below the threshold set as "PTEN presence", otherwise set "PTEN intact"
  sample_info_extracted$H_score_threshold[which(sample_info_extracted$H_score_threshold <= threshold[i] & !is.na(sample_info_extracted$H_score_threshold))] <- "presence"
  sample_info_extracted$H_score_threshold[which(sample_info_extracted$H_score_threshold != "presence" & !is.na(sample_info_extracted$H_score_threshold))] <- "intact"
  
  if (i == 1){
    cut_0 <- sample_info_extracted$H_score_threshold
  } else {
    if( i==2){
      cut_10<- sample_info_extracted$H_score_threshold
    } else {
      if(i== 3){
        cut_20 <- sample_info_extracted$H_score_threshold
      }else {
        if(i == 4){
          cut_30 <- sample_info_extracted$H_score_threshold
        }else {
          if(i == 29){
            cut_200 <- sample_info_extracted$H_score_threshold
          }
        }
      }
    }
  }
  
  #| PI3K-AKT-mTOR signature (105 genes)
  test_wilcox_PI3K_AKT_mTOR <- wilcox.test(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "presence")],sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "intact")] )
  p_value_wilcox_PI3K_AKT_mTOR <- c(p_value_wilcox_PI3K_AKT_mTOR, test_wilcox_PI3K_AKT_mTOR$p.value) 
  difference_wilcox_PI3K_AKT_mTOR <- c(difference_wilcox_PI3K_AKT_mTOR, mean(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "presence")]) - mean(sample_info_extracted$Mean.expression.PI3K.AKT.mTOR[which(sample_info_extracted$H_score_threshold == "intact")] )) 
  
}

#| Defining directionality
difference_wilcox_PI3K_AKT_mTOR[which(difference_wilcox_PI3K_AKT_mTOR <0)] <- "Negative"
difference_wilcox_PI3K_AKT_mTOR[which(difference_wilcox_PI3K_AKT_mTOR >=0 & difference_wilcox_PI3K_AKT_mTOR != "Negative")] <- "Positive"


#| DATA FOR PLOTS 
data_frame_p_value_wilcoxon_threshold <- data.frame(threshold =threshold[-length(threshold)],
                                                    p_value_wilcox_PI3K_AKT_mTOR= p_value_wilcox_PI3K_AKT_mTOR,
                                                    difference_wilcox_PI3K_AKT_mTOR=difference_wilcox_PI3K_AKT_mTOR)


#| PI3K-AKT-mTOR signature (105 genes)
ggplot(data_frame_p_value_wilcoxon_threshold, aes(x=threshold, y =p_value_wilcox_PI3K_AKT_mTOR, color = difference_wilcox_PI3K_AKT_mTOR)) +
  geom_point(size =5)+
  theme(text=element_text(size=13,  family="serif"), 
        axis.text=element_text(size=12), 
        plot.title=element_text(size=12, hjust = 0.5, face ="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(size =10)) + 
  ylab("p-values of wilcoxon test \n PI3K-AKT-mTOR pathway signature")+ 
  xlab("threshold of H-score") +
  scale_color_manual(values = c("black", "azure3")) +
  labs(color = "Mean differences PTEN\n presence - intact", size =2)
ggsave("Results/Correlation_plots/FINDING_H_SCORE_THRESHOLD_p_value_PI3K_AKT_mTOR.pdf", height = 4, width = 5.4)


#| Ploting values for different cutoffs

sample_info_extracted$H_score_cut_0 <- cut_0
sample_info_extracted$H_score_cut_10 <- cut_10
sample_info_extracted$H_score_cut_20 <- cut_20
sample_info_extracted$H_score_cut_30 <- cut_30
sample_info_extracted$H_score_cut_200 <- cut_200

#| H_score_cut_0
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_0)),], aes(x = H_score_cut_0, y = Mean.expression.PI3K.AKT.mTOR, fill = H_score_cut_0)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = H_score_cut_0),position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),pch = 21, size=2.6, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=-0.7, hjust=-1,family="sans", size =3, method="wilcox")+
  theme(text=element_text(size=12,  family="sans"), axis.text=element_text(size=12), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold",),
        legend.position = "none") + 
  ylab("Mean expression PI3K-AKT-mTOR\n (Z score of Normalized log2 counts)")+ 
  scale_fill_manual(values =c("#504667", "#EFD5B2")) +
  labs(fill ="PTEN status") +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_0_PI3K-AKT-mTOR_PTEN_intact_presence.pdf", height = 4.2, width = 3.9)

#| H_score_cut_10
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_10)),], aes(x = H_score_cut_10, y = Mean.expression.PI3K.AKT.mTOR, fill = H_score_cut_10)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = H_score_cut_10),position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),pch = 21, size=2.6, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=-0.7, hjust=-1,family="sans", size =3, method="wilcox")+
  theme(text=element_text(size=12,  family="sans"), axis.text=element_text(size=12), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold",),
        legend.position = "none") + 
  ylab("Mean expression PI3K-AKT-mTOR\n (Z score of Normalized log2 counts)")+ 
  scale_fill_manual(values =c("#504667", "#EFD5B2")) +
  labs(fill ="PTEN status") +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_10_PI3K-AKT-mTOR_PTEN_intact_presence.pdf", height = 4.2, width = 3.9)

#| H_score_cut_20
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_20)),], aes(x = H_score_cut_20, y = Mean.expression.PI3K.AKT.mTOR, fill = H_score_cut_20)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = H_score_cut_20),position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),pch = 21, size=2.6, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=-0.7, hjust=-1,family="sans", size =3, method="wilcox")+
  theme(text=element_text(size=12,  family="sans"), axis.text=element_text(size=12), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold",),
        legend.position = "none") + 
  ylab("Mean expression PI3K-AKT-mTOR\n (Z score of Normalized log2 counts)")+ 
  scale_fill_manual(values =c("#504667", "#EFD5B2")) +
  labs(fill ="PTEN status") +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_20_PI3K-AKT-mTOR_PTEN_intact_presence.pdf", height = 4.2, width = 3.9)

#| H_score_cut_30
ggplot(sample_info_extracted[which(!is.na(sample_info_extracted$H_score_cut_30)),], aes(x = H_score_cut_30, y = Mean.expression.PI3K.AKT.mTOR, fill = H_score_cut_30)) +
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
  geom_point(aes(fill = H_score_cut_30),position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7),pch = 21, size=2.6, alpha=0.6)+
  stat_compare_means(label ="p.format", vjust=-0.7, hjust=-1,family="sans", size =3, method="wilcox")+
  theme(text=element_text(size=12,  family="sans"), axis.text=element_text(size=12), 
        plot.title=element_text(size=16, hjust = 0.5, face ="bold",),
        legend.position = "none") + 
  ylab("Mean expression PI3K-AKT-mTOR\n (Z score of Normalized log2 counts)")+ 
  scale_fill_manual(values =c("#504667", "#EFD5B2")) +
  labs(fill ="PTEN status") +
  xlab("PTEN protein status")
ggsave("Results/Box_plots/H_score_cut_30_PI3K-AKT-mTOR_PTEN_intact_presence.pdf", height = 4.2, width = 3.9)

#| What cut it is more suitable to use? H-SCORE AND SAMPLE PROPORTION
fraction_samples_PTEN_presence_intact <- data.frame(group =c("PTEN presence", "PTEN intact", "PTEN presence", "PTEN intact", "PTEN presence", "PTEN intact", "PTEN presence", "PTEN intact"),
                                                Samples = c(length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$H_score_cut_0== "PTEN presence")])/length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0) & sample_info_extracted$H_score_cut_0== "PTEN intact")])/length(sample_info_extracted$H_score_cut_0[which(!is.na(sample_info_extracted$H_score_cut_0))]), 
                                                            length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_10) & sample_info_extracted$H_score_cut_10== "PTEN presence")])/length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_10) & sample_info_extracted$H_score_cut_10== "PTEN intact")])/length(sample_info_extracted$H_score_cut_10[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_20) & sample_info_extracted$H_score_cut_20== "PTEN presence")])/length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_20) & sample_info_extracted$H_score_cut_20== "PTEN intact")])/length(sample_info_extracted$H_score_cut_20[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_30) & sample_info_extracted$H_score_cut_30== "PTEN presence")])/length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_0))]),
                                                            length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_30) & sample_info_extracted$H_score_cut_30== "PTEN intact")])/length(sample_info_extracted$H_score_cut_30[which(!is.na(sample_info_extracted$H_score_cut_0))])),
                                                cut = c("0", "0", "10", "10", "20", "20", "30", "30"))

fraction_samples_PTEN_presence_intact$round <- round(fraction_samples_PTEN_presence_intact$Samples, 2)
ggplot(fraction_samples_PTEN_presence_intact, aes(x =cut, y = Samples, fill =group)) +
  geom_bar(stat = "identity") + scale_fill_manual(values =c("#504667", "#EFD5B2")) +
  theme(text=element_text(size=16,  family="sans"), axis.text=element_text(size=12), plot.title=element_text(size=16, hjust = 0.5, face ="bold")) + 
  ylab("Fraction of samples") +
  xlab("H-score cutoff")+
  labs(fill ="PTEN protein status") +
  geom_text(aes(label=round),position="stack",vjust=1, colour = "white")#+ 
  #ggtitle("Fraction of samples for PTEN presence and PTEN intact changes depending on the \n H-score cut-off (147 PCa samples)")
ggsave("Results/Bar_plots/Fraction_samples_PTEN_presence_PTEN_intact_Different_H-score_cutoff_147_PCa_samples.pdf", heigh= 4.5, width = 5.5)
################################################################################


################################################################################
#| SAVING SAMPLE INFORMATION WITH H-SCORE VALUES
################################################################################
write.table(sample_info_extracted, "Results/Sample_info_table/sample_info_extracted.txt", sep ="\t",row.names=T)
################################################################################

