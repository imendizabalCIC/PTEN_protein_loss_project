################################################################################
#| SCRIPT TO OBTAIN RAW COUNT MATRIX FROM STAR OUTPUT 
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon
#| 
#| Description:
#| This script generates a raw count matrix from STAR aligner results.
#| Specifically, it extracts read counts from the `ReadsPerGene.out.tab` files 
#| produced by STAR, keeping the correct strandedness column (reverse strand for 
#| SMARTer Stranded Total RNA Pico libraries).
#|
#|
#| Input:
#|   - STAR outputs (*.ReadsPerGene.out.tab) located in:
#|       X:/irondon/AC-12_RNAseq/02_STAR/Outdir_STAR/
#|
#| Output:
#|   - Raw count matrix (ENSEMBL gene IDs as rows, samples as columns), saved at:
#|       X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/AC-12_RNAseq_Raw_Counts.txt
#|
#| Notes on strandedness:
#|   The `ReadsPerGene.out.tab` file has 4 columns:
#|     1) Gene ID
#|     2) Unstranded counts
#|     3) Forward-stranded counts
#|     4) Reverse-stranded counts
#|   Since this library prep is SMARTer Stranded Total RNA Pico (paired-end, 
#|   cDNA synthesis by SMARTScribe RT), we use column 4 (reverse strand).
################################################################################



################################################################################
#| Into the ReadsPerGene.out.tab we can find 4 columns which correspond to different
#| strandedness options:
#|
#| 1) column 1: gene ID
#| 2) column 2: counts for unstranded RNA-seq
#| 3) column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#| 4) column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
#|
#| This information comes from the sequencer. The library used is SMARTer Stranded 
#| Total RNA Pico. The sequencing is paired-end. The 1st strand cDNA synthesis was 
#| performed using SMARTScribe reverse Transcriptase. Hence you would use the 4th 
#| column (Reversely stranded)
strandSpecific <- "stranded"
stranded <- 3 
#|
#|if (strandSpecific == "unstranded") {
#|  stranded <- 1
#|} else if (strandSpecific == "stranded") {
#|  stranded <- 2 # Forward strand
#|} else {
#|  stranded <- 3 # Reversely stranded    
#|} 
#|
workingDir <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/"
setwd(workingDir)
################################################################################

################################################################################
#|  DATA INPUT  
################################################################################

#| Outdir directory from STAR
Out_dir <- "X:/irondon/AC-12_RNAseq/02_STAR/Outdir_STAR/"

#| Selecting Ficheros
listaFicheros <- list.files(path = Out_dir, pattern = "*ReadsPerGene.out.tab")
listaFicheros <- gsub("_STARReadsPerGene.out.tab", "", listaFicheros)

#| Directory to save the counts
counts.dir <- "X:/irondon/AC-12_RNAseq/03_RAW_COUNTS/"

#| Tag name 
tag.name <- "AC-12_RNAseq_Raw_Counts"
################################################################################

################################################################################
#|  LOOP  
################################################################################

#| To be able to run the lines downstream for second time
if (exists("RawCounts")) {
  rm(RawCounts)
} 

#|  LOOP FOR THE ANALYSIS OF ALL THE FILES CONTAINING ReadsPerGene.out.tab  ###
for (file in listaFicheros){
  
  #| Reading the files
  temp <- read.table(paste(Out_dir, file, "_STARReadsPerGene.out.tab", sep=""), sep="\t", header=F)
  
  #| Creating a data frame to store the value of the V1 column and selecting the counts in the column number 1+stranded
  datosTemp <- data.frame(temp$V1, temp[, (1 + stranded)])
  
  #| Deleting the first 4 rows of the data.frame (which contains useless info)
  datosTemp <- datosTemp[-c(1:4), ]
  
  #| Changing the name of the columns to agree with the name of the sample (per patient)
  colnames(datosTemp) <- c("Sample", strsplit(file, "\\.")[[1]][1])
  
  #| Check if the RawCounts file exist or not to create the data frame of the counts
  if (exists("RawCounts")) {
    RawCounts <- merge(RawCounts, datosTemp, by="Sample")
  } else {
    RawCounts <- datosTemp
  }
}

#| Testing if there is any duplicate inside our generated RawCounts
any(duplicated(RawCounts))

#| Assigning the ENSEMBL ID to the row name of the counts matrix
rownames(RawCounts) <- RawCounts$Sample

#| Deleting Sample column
RawCounts <- RawCounts[,-c(1)]

#| Sorting by rownames (Genes)
RawCounts <- RawCounts[order(row.names(RawCounts)),]

#| Sorting by colnames (Samples)
RawCounts <- RawCounts[order(colnames(RawCounts))]

#| Save the count dataframe
write.table(RawCounts, paste(counts.dir, tag.name, ".txt", sep =""), quote=F, row.names=T, sep="\t") 
################################################################################