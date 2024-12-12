#/opt/R/R-4.1.2/bin/R

#| Last change: 12/12/2024
#| Ivana Rondon

###############################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FASTQ FILE AND SUMIT JOBS ########
####### TO PERFORM TRIMMING ON THE FASTQ FILES OF THE SECOND READ      ########
###############################################################################

#| This script creates and execute 198 jobs for trimming corresponding sequences 

#| The raw data analyzed is contained in the folder /vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/
#| This data is the RNAseq from the prostate of 198 patients from Basurto's hospital
#| that was extracted after detecting PCa.

#| The library used is SMART SMARTer(a) Stranded Total RNASeq Kit v2 - Pico Input 
#| Mammalian (Cat. Nos. 634411, 634412, 634413, 63441) in an Illumina NovaSeq 
#| intrument. The output is paired-end sequencing.

#| The preliminar instructions from the sequencer are:
#|   1) When performing paired-end sequencing, the first three nucleotides of 
#|   the second sequencing read (Read 2) are derived from the Pico v2 SMART Adapter.
#|   These three nucleotides must be trimmed prior to mapping
#|   2) Cut the TruSeq adapters
#|   3) Conclusions of the trimming of the Read 1: It is necessary to trim the 3 
#|   nucleotides of Read 1 adjacent of the 3' adapter (because the inserts contains
#|   a size < 150bps

#| In summary: This scripts creates a .sh file and then it is launched to the job 
#| executor in indar to:
#|   1) Cut the TruSeq CD illumina adaptors from the fastq files
#|   2) Cut the first 3 nucleotides of Read 2. To be sure about the quality of our
#| data given that we have inserts with <150bps, it is necessary to cut also the 
#| first 3 nucleotides adjacent to the adapter of Read 1.
#|   1) Perform quality trimming with a cutoff of 10
#|   2) Minimum length size of 20bp. This is chosen based on the distribution of 
#| length for our insert.

#| From the cutadapt documentation (paired-end sequence)
#|     cutadapt [options] -o output 1 -p output2 input1 input2
#| [options]: 
#|     -a: Cut adapters of the 3' of Read 1
#|     -A: Cut adapters of the 3' of Read 2
#|     -m: Defines a threshold of the minimum lenght sequence to keep
#|     -q: Defines a threshold of the quality for each sequence
#|     -u: Cut a specific number of nucleotides of Read 1
#|     -U: Cut a specific number of nucleotides of Read 2

###############################################################################


################################################################################
#|  DATA LIST  
################################################################################

#| List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs", pattern = "*_2.fastq")
samples <- gsub("_2.fastq.gz", "", samples)

################################################################################


################################################################################
#|  VARIABLES  
################################################################################

#| Defining variables
cpu <- "8"
memory <- "1G"
time <- "5:00:00"
partition <- "FAST"

################################################################################


################################################################################
#|  JOBS LOOP  
################################################################################

#| This creates a loop over all the samples and create a .sh file which is then submitted as a job
for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste("TRIMMING_", samples[s], "_2.sh", sep = "") 

  # Forward strand  (Read 1)
  input1 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/", samples[s], "_1.fastq.gz", sep = "")
  input2 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs/", samples[s], "_2.fastq.gz", sep = "")

  # Output directories
  output1_out <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_1_out.fastq.gz", sep = "") 
  output2_out <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_2_out.fastq.gz", sep = "")
  
  output1 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_1.fastq.gz", sep = "") 
  output2 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_2.fastq.gz", sep = "")
  
  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o", samples[s], "_TRIM.out", sep = ""),
    paste("#SBATCH -e", samples[s], "_TRIM.err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    "#SBATCH --mail-user=irondon@cicbiogune.es",
    "#SBATCH --mail-type=FAIL",
    "\n\n",
    "source /share/apps/anaconda3/cic-env", 
    "conda config --add channels defaults",
    "conda config --add channels bioconda",
    "conda config --add channels conda-forge",
    "conda config --set channel_priority strict",
    "conda create -n cutadaptenv cutadapt",
    "conda activate cutadaptenv",
    "\n",
    paste("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o", output1_out,"-p", output2_out, input1, input2, sep = " "),
    paste("cutadapt -q 10 -m 20 -u -3 -U 3 -o", output1, "-p", output2, output1_out, output2_out, sep = " "),
    file = filename, sep = "\n", append = FALSE)
  system(paste("sbatch", filename, sep = " "))
}
################################################################################
