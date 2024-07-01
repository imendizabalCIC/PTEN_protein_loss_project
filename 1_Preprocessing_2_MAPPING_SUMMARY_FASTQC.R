################################################################################
#                   SUMMARY OF THE FASTQC ANALYSIS
################################################################################

# Once the fastqc software is applied to the fastq files obtained from the RNASeq
#we need to summarize the results obtained on each files (the zip files) regarding
#sequence quality, adapter content, etc..

# This script is written to evaluate the results from the FastQC applied to the
#AC-45_RNAseq-FFPE data

################################################################################


################################################################################
#|  LIBRARIES AND DATA
################################################################################

#| Libraries
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

#| Setting working directory
workingDir = "X:/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/1_STEP_FASTQCs/"
setwd(workingDir)

# List with the names of the files we will process, only keeping the sample name.
Samples <- list.files(path = "X:/DATA_shared/AC-45_RNAseq-FFPE/FASTQs", pattern = "*.fastq")
Samples <- gsub(".fastq.gz", "", Samples)
Samples <- sort(Samples)

# Empty dataframe to load all data
fastqc_summary <- data.frame()

################################################################################


################################################################################
#|  LOOP OVER EVERY SAMPLE AND FORMATING DATAFRAME  
################################################################################

#| With a loop I go through the Samples character list
for(s in 1:length(Samples)){
  
  #| Reading the fastqc summary file, but I have to unzip the directory in the process
  base_table <- read.table(unz(paste("X:/DATA_shared/AC-45_RNAseq-FFPE/FASTQCs/", Samples[s], "_fastqc.zip", sep = ""), paste(Samples[s], "_fastqc/summary.txt", sep = "")), sep = "\t")
  
  #| Keeping the PASS/FAIL/WARN results
  results <- base_table[,1]
  
  #| Adding the name of the sample in the first column of the dataframe (row number is "s")
  fastqc_summary[s,1] <- Samples[s]
  
  #| With a for loop, I print every PASS/FAIL/WARN data of that sample in the next columns
  for(k in 1:nrow(base_table)){
    fastqc_summary[s, k+1] <- results[k]
  }
}

#| Taking the econd column of the base_table that was recollected in the loop. This column has the names of the different studies
col_names <- base_table[,2]

#| Adding the name of the first column to the object
col_names <- c("FastQC Sample", col_names)

#| Applying the colnames to the data frame
colnames(fastqc_summary) <- col_names

################################################################################


################################################################################
#|  PLOTTING  
################################################################################

#| To only plot the studies with the results (without samples), I first have to transpose the data.frame (what will transform it to a matrix, so I retransforme it into a dataframe to edit later) and delete the first column
fastqc_summary_transpose <- as.data.frame(t(fastqc_summary[-1]))

#| Adding a new column with the rownames called "Studies"
fastqc_summary_transpose$Study <- rownames(fastqc_summary_transpose)

#| Now I have to use the melt() function to put the dataframe in long format by the column "Study"
fastqc_summary_molten <- melt(as.data.table(fastqc_summary_transpose), id.vars="Study")

#| Plotting and saving
pdf(file = "FASTQCs_SUMMARY_RESULTS/FASTQsFASTQC_Summary_OVERVIEW.pdf")
ggplot(fastqc_summary_molten, 
       aes(x = Study, 
           fill = value)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("PASS" = "#52bd46", "FAIL" = "#e33a14", "WARN" = "#f5b40f")) +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("FASTQC Summary Overview")
dev.off()

################################################################################


################################################################################
#|  INTERPRETATION  
################################################################################

#| In RNAseq studies, we have to make an special treatment of two parameteres:
#| 
#|      - Per base sequence content: 
#|      Always gives a FAIL for RNA-seq data.This is because the first 10-12 bases 
#|      result from the 'random' hexamer priming that occurs during RNA-seq library 
#|      preparation. This priming is not as random as we might hope giving an 
#|      enrichment in particular bases for these initial nucleotides.
#|      
#|      - Sequence Duplication Levels:
#|      This plot can help identify a low complexity library, which could result 
#|      from too many cycles of PCR amplification or too little starting material. 
#|      For RNA-seq we don't normally do anything to address this in the analysis, 
#|      but if this were a pilot experiment, we might adjust the number of PCR 
#|      cycles, amount of input, or amount of sequencing for future libraries.

#|  OBSERVATIONS FROM THE FASTQC SUMMARY:
#|      - Over-represented sequences: This module will often be triggered when used
#|       to analyse small RNA libraries where sequences are not subjected to random 
#|       fragmentation, and the same sequence may naturally be present in a significant 
#|       proportion of the library.
#|       
#|      - G/C content: give you information about possible contamination. Warnings 
#|      in this module usually indicate a problem with the library. Sharp peaks
#|      on an otherwise smooth distribution are normally the result of a specific 
#|      contaminant (adapter dimers for example), which may well be picked up by 
#|      the overrepresented sequences module. Broader peaks may represent contamination 
#|      with a different species.

################################################################################