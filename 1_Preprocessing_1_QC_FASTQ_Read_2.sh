#!/bin/bash

#################################################################################
### QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS ##
###    READ 2
#################################################################################

# FastQC reads a set of sequence files and produces from each one a quality control report
#consisting of a number of different modules, each one of which will help to identify a di-
#fferent potential type of problem in your data. Basically, the command line has to follow:
#
#    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#          [-c contaminant file] seqfile1 .. seqfileN
#where:
#   -h --help:     Print the help file
#   -v --version:  Print the version of the program
#   -o --outdir:   Creates an output folder with the same name where all the outfiles will
#be stored
#   --casava:      Files come from raw casava output. Files in the same sample group (diffe-
#ring only by the group number) will be analysed as a set rather than individually. Files 
#must have the same names given to them by casava (including being gzipped and ending with 
#.gz) otherwise they won't be grouped together correctly.
#   --nano:        Files come from nanopore sequences and are in the fast5 format
#   -t --threads   Specifies the number of files which can be processed simultaneously.
#Take into account that each thread will allocate 250MB of memory, so you shouldn't run 
#more threads than your available memory will cope with, and not more than 6 threads on a 
#32 bit machine

##################################################################################


#SBATCH --job-name="FASTQC_READ_2"
#SBATCH -o FASTQC_READ2.out
#SBATCH -e FASTQC_READ2.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=3GB
#SBATCH --partition=NORMAL
#SBATCH --time=24:00:00
#SBATCH --exclude=gn[05-09]
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

# Defining the working directory
projectDir="/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq"
numProcessors=8  # max 16 with Indar and 32 with Nextera
FASTQs_out="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQCs"
FASTQs_dir="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed"
cd $projectDir

# Create folder to store the output values. If the folder exist it does not modify the folder
mkdir $FASTQs_out 

# Quality control of the sequences from FASTqs files (all the files in the folder)
fastqc $FASTQs_dir/*_2.fastq.gz -o $FASTQs_out/ -t $numProcessors