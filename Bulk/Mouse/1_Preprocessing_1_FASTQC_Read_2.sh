#!/bin/bash

#################################################################################
#| QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS ##
#|    READ 2
#################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This SLURM job runs FastQC on all Read 2 FASTQ files from the AC-12_RNAseq (mouse)
#| project. FastQC generates quality control reports (HTML + zipped folders) for 
#| each sample, assessing per-base quality, GC content, adapter contamination, and 
#| other metrics. 
#|
#| Workflow:
#|   - Input: Trimmed Read 2 FASTQ files (*.fastq.gz) from $FASTQs_dir
#|   - Tool:  FastQC (run with 2 CPUs)
#|   - Output: Quality reports stored in $FASTQs_out
#|
#| Notes:
#|   - Only Read 2 files are analyzed in this script
#|   - Folder $FASTQs_out is created if it does not already exist.
##################################################################################

#SBATCH --job-name="FASTQC_READ_2"
#SBATCH -o FASTQC_READ2.out
#SBATCH -e FASTQC_READ2.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=3GB
#SBATCH --partition=FAST
#SBATCH --time=12:00:00

# Defining the working directory
projectDir="/vols/GPArkaitz_bigdata/irondon/AC-12_RNAseq/1_FASTQCs/"
FASTQs_out="/vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQCs_trimmed"
FASTQs_dir="/vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQs_trimmed"
cd $projectDir


# Quality control of the sequences from FASTqs files (all the files in the folder)
/vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/FastQC/fastqc $FASTQs_dir/*_out_2.fastq.gz -o $FASTQs_out/ 