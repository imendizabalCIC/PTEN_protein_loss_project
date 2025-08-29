#!/bin/bash

#################################################################################
#| QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS ##
#|    READ 2
#################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This SLURM job runs FastQC on all Read 2 FASTQ files from the AC-45_RNAseq-FFPE 
#| project. FastQC generates quality control reports (HTML + zipped folders) for 
#| each sample, assessing per-base quality, GC content, adapter contamination, and 
#| other metrics. 
#|
#| Workflow:
#|   - Input: Trimmed Read 2 FASTQ files (*.fastq.gz) from $FASTQs_dir
#|   - Tool:  FastQC (run with 8 CPUs)
#|   - Output: Quality reports stored in $FASTQs_out
#|
#| Notes:
#|   - Only Read 2 files are analyzed in this script
#|   - Folder $FASTQs_out is created if it does not already exist.
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

#| Defining the working directory
projectDir="/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq"
numProcessors=8  
FASTQs_out="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQCs"
FASTQs_dir="/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed"
cd $projectDir

#| Create folder to store the output values. If the folder exist it does not modify the folder
mkdir $FASTQs_out 

#| Quality control of the sequences from FASTqs files (all the files in the folder)
fastqc $FASTQs_dir/*_2.fastq.gz -o $FASTQs_out/ -t $numProcessors