#!/bin/bash

#| Last change: 12/12/2024
#| Ivana Rondon-Lorefice

#################################################################################
#| QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS 
#| READ 1 
#################################################################################

#|  FastQC reads a set of sequence files and produces from each one a quality control report
#| consisting of a number of different modules, each one of which will help to identify a di-
#| fferent potential type of problem in your data. Basically, the command line has to follow:
#| 
#|     fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#|           [-c contaminant file] seqfile1 .. seqfileN
#| where:
#|    -h --help:     Print the help file
#|    -v --version:  Print the version of the program
#|    -o --outdir:   Creates an output folder with the same name where all the outfiles will
#| be stored
#|    --casava:      Files come from raw casava output. Files in the same sample group (diffe-
#| ring only by the group number) will be analysed as a set rather than individually. Files 
#| must have the same names given to them by casava (including being gzipped and ending with 
#| .gz) otherwise they won't be grouped together correctly.
#|    --nano:        Files come from nanopore sequences and are in the fast5 format
#|    -t --threads   Specifies the number of files which can be processed simultaneously.
#| Take into account that each thread will allocate 250MB of memory, so you shouldn't run 
#| more threads than your available memory will cope with, and not more than 6 threads on a 
#| 32 bit machine

##################################################################################


##################################################################################
#|  JOB for FASTQC in Read 1
##################################################################################
#!/bin/sh
#SBATCH --job-name="FASTQC_READ_1"
#SBATCH -o FASTQC_READ1.out
#SBATCH -e FASTQC_READ1.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3GB
#SBATCH --partition=FAST
#SBATCH --time=12:00:00

# Defining the working directory
projectDir="/vols/GPArkaitz_bigdata/irondon/MG-05_TotalRNAseq_processing_handled/1_FASTQCs/"
FASTQs_out="/vols/GPArkaitz_bigdata/DATA_shared/MG-05_TotalRNAseq/FASTQCs_trimmed"
FASTQs_dir="/vols/GPArkaitz_bigdata/DATA_shared/MG-05_TotalRNAseq/FASTQs_trimmed"
cd $projectDir


# Quality control of the sequences from FASTqs files (all the files in the folder)
/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/FastQC/fastqc $FASTQs_dir/*_out_1.fastq.gz -o $FASTQs_out/ 
