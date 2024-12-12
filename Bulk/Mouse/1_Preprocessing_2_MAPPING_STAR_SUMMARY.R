#| Last change: 12/12/2024
#| Ivana Rondon-Lorefice

###############################################################################
############### ANALYSIS OF THE RESULTS OBTAINED WITH STAR ####################
###############################################################################

#| This script creates a dataframe that summarizes the mapping process results,  
#| such as errors in the process, mapped length, percentage of unique mapping, etc...

###############################################################################


################################################################################
#| LIBRARIES AND WORKING DIRECTORY SPECIFICATIONS
################################################################################
library(ggplot2)
workingDir = "X:/irondon/AC-12_RNAseq/02_STAR/"
setwd(workingDir)
################################################################################


################################################################################
#|  DATA INPUT  
################################################################################

#| FASTQs trimmed file directory. In this scripts, trimmed fastqs ends with _out_*.sh
fastqs_dir <- "X:/DATA_shared/AC-12_RNAseq/FASTQs/"

#| Results file 
results.file <- "X:/irondon/AC-12_RNAseq/02_STAR/"

#| Outdir directory
Out_dir <- "X:/irondon/AC-12_RNAseq/02_STAR/"

#| Creating a list with the names of the files we will process, only keeping the sample name.
Samples <- list.files(path = fastqs_dir, pattern = "*_1.fastq")
Samples <- gsub("_1.fastq.gz", "", Samples)

#| Number of job
slurm_number <- Samples

#| Creating a dataframe whose first column is the samples names (from Samples)
mapping_summary <- data.frame(Samples)
################################################################################
  

################################################################################
#|  SLURM.OUT LOOP  
################################################################################
#| This loop runs through the numbers of the slurm IDs. In this files we will see if the mapping successfully finished.
for(i in slurm_number){
  
  #| Reading the slurm.out file
  slurm <- read.table(paste(Out_dir, i, ".out", sep = ""), sep = "\t")
  
  #| Taking the name of the sample from the slurm.out file to know where to save the data in the df
  name <- slurm[9,]
  name <- gsub("Job name: ", "", name)
  
  #| Looking for the string "finished successfully" in the column 27. success = TRUE if I find it
  ending <- slurm[27,]
  success <- grepl("finished successfully", ending)
  
  #| If success = TRUE, I write in the df "Finished successfully". If not, I assume there has been some type of Error.
  if(success){
    mapping_summary$Mapping_Process[mapping_summary$Samples == name] <- "Finished successfully"
  }else{
    mapping_summary$Mapping_Process[mapping_summary$Samples == name] <- "Error"
  }
}
################################################################################


################################################################################
#|  STARLOG.FINAL.OUT LOOP  
################################################################################

#| For this loop I go through the Samples character list
for(k in 1:length(Samples)){
  
  #| I read the STARLog.final.out. I use read.delim since read.table doesn't work in this case
  results <- read.delim(paste(Out_dir, "Outdir_STAR/", Samples[k], "_STARLog.final.out", sep = ""), sep = "\t")
  
  #| Now I take each cell with useful information and create a column in the df
  mapping_summary$Read_Length[mapping_summary$Samples == Samples[k]] <- results[5,2]
  mapping_summary$Read_Number[mapping_summary$Samples == Samples[k]] <- results[7,2]
  mapping_summary$Mapped_Unique[mapping_summary$Samples == Samples[k]] <- results[8,2]
  mapping_summary$Mapped_Length[mapping_summary$Samples == Samples[k]] <- results[9,2]
  mapping_summary$Mapped_MMrate[mapping_summary$Samples == Samples[k]] <- results[16,2]
  mapping_summary$Mapped_Multi[mapping_summary$Samples == Samples[k]] <- results[23,2]
  mapping_summary$Unmapped_MM[mapping_summary$Samples == Samples[k]] <- results[28,2]
  mapping_summary$Unmapped_Short[mapping_summary$Samples == Samples[k]] <- results[30,2]
  mapping_summary$Unmapped_Other[mapping_summary$Samples == Samples[k]] <- results[32,2]
}
################################################################################


################################################################################
#|  PLOTTING  
################################################################################

#| Duplicating the dataframe to avoid modifying the original one
mapping_summary_num <- mapping_summary

#| Saving the name of the columns i will represent. I delete the % and make numeric using lapply
num_col <- c("Mapped_Unique","Mapped_MMrate", "Mapped_Multi", "Unmapped_MM", "Unmapped_Short", "Unmapped_Other")
mapping_summary_num[num_col] <- lapply(mapping_summary_num[num_col], gsub, pattern = "%", replacement = "")
mapping_summary_num[num_col] <- lapply(mapping_summary_num[num_col], as.numeric)

# Keeping only the columns I am interested in
mapping_summary_boxplot <- mapping_summary_num[num_col]

# Representing the data in a boxplot
bp <- ggplot(stack(mapping_summary_boxplot), aes(x = factor(ind, levels = names(mapping_summary_boxplot)), y = values), fill = factor(ind, levels = names(mapping_summary_boxplot))) +geom_boxplot(alpha=0.3) 
bp + theme_classic() + scale_fill_brewer(palette="BuPu")
ggsave(paste(results.file, "2_STAR_Summary_Results.pdf", sep ="" ), height =5, width =7)
################################################################################