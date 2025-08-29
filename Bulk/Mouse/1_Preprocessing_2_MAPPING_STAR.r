#/opt/R/R-4.1.2/bin/R

###############################################################################
####### THIS CODES CREATES A SCRIPT FOR EACH FILE FASTA AND SUMIT JOBS ########
###############################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon-Lorefice
#|
#| Description:
#| This script scans a FASTQ directory, derives sample IDs from paired-end files
#| (<sample>_1.fastq.gz and <sample>_2.fastq.gz), and for each sample:
#|   * Writes a per-sample SLURM shell script that runs STAR alignment
#|   * Submits the job to the Indar cluster via `sbatch`
#|
#| Use case
#| Mapping RNA-seq reads from a prostate cancer cohort of mice to the
#| mouse reference (STAR genome index), producing coordinate-sorted BAMs and
#| per-gene counts.
#|
#| Key STAR options (per job)
#|   --genomeDir <index>                    : STAR genome index
#|   --runThreadN <cpu>                     : threads per sample
#|   --readFilesIn R1 R2                    : gzipped PE FASTQs
#|   --readFilesCommand gunzip -c           : on-the-fly decompression
#|   --outFilterMultimapNmax 1              : keep uniquely mapped reads
#|   --outReadsUnmapped Fastx               : write unmapped mates
#|   --outSAMtype BAM SortedByCoordinate    : sorted BAM output
#|   --twopassMode Basic                    : 2-pass mapping
#|   --limitBAMsortRAM 2000000000           : RAM for BAM sorting (bytes)
#|   --quantMode TranscriptomeSAM GeneCounts: transcriptome BAM + gene counts
#|
#| Cluster settings (per job)
#|   * CPUs: 8
#|   * Memory: 40G
#|   * Time limit: 10:00:00
#|   * Partition: FAST
#|   * Email on FAIL: irondon@cicbiogune.es
#|
#| Inputs
#|   * FASTQs (gz): /vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQs/
#|   * STAR index : /vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Indexes/Mouse_101
#|
#| Outputs
#|   * Per-sample STAR results under:
#|       /vols/GPArkaitz_bigdata/irondon/AC-12_RNAseq/02_STAR/Outdir_STAR/<sample>_STAR*
#|   * Per-sample SLURM scripts: STAR_<sample>.sh
#|   * SLURM logs: <sample>.out / <sample>.err (in working dir)
#|
#| Notes
#|   * Samples are inferred from *_1.fastq.gz; matching *_2.fastq.gz must exist.
#|   * Adjust CPU/memory/time/partition to your queue availability if needed.
################################################################################


################################################################################
#|  DATA INPUT  
################################################################################

#| FASTQs file directory
fastqs_dir <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQs/"

#| FASTQs trimmed file directory. In this scripts, trimmed fastqs ends with _out_*.sh
#fastqs_trimmed_dir <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQs_trimmed"

#| email for reporting failing of jobs
mail <- "irondon@cicbiogune.es"

# List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = "/vols/GPArkaitz_bigdata/DATA_shared/AC-12_RNAseq/FASTQs", pattern = "*_1.fastq")
samples <- gsub("_1.fastq.gz", "", samples)

#| Genome directory
genomedir <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes_Rocky/Indexes/Mouse_101"  

#| Outdir directory
Out_dir <- "/vols/GPArkaitz_bigdata/irondon/AC-12_RNAseq/02_STAR/Outdir_STAR/"

#| Cluster requirements
cpu <- "8"
memory <- "40G"
time <- "10:00:00"
partition <- "FAST"
################################################################################


################################################################################
#|  JOBS LOOP  
################################################################################
for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste("STAR_", samples[s], ".sh", sep = "") 

  # Forward strand
  filedir1 <- paste(fastqs_dir, samples[s], "_1.fastq.gz", sep = "")

  # Reverse strand
  filedir2 <- paste(fastqs_dir, samples[s], "_2.fastq.gz", sep = "")
  
  # Out directory
  outdir <- paste(Out_dir, samples[s], "_STAR", sep = "") 
  
  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o ", samples[s], ".out", sep = ""),
    paste("#SBATCH -e ", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    paste("#SBATCH --mail-user=", mail, sep = ""),
    "#SBATCH --mail-type=FAIL",
    "\n\n",
    "source /opt/ohpc/pub/apps/star/STAR-2.7.10a/cic-env", 
    "\n",
    paste("STAR --genomeDir", genomedir, "--runThreadN", cpu, "--readFilesIn", filedir1, filedir2, "--readFilesCommand gunzip -c --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --limitBAMsortRAM 2000000000 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix", outdir, sep = " "),
    file = filename, sep = "\n", append = FALSE)
  system(paste("sbatch", filename, sep = " "))

}
################################################################################