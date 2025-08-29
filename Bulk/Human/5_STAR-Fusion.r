#/opt/R/R-4.1.2/bin/R

################################################################################
#| STAR-Fusion JOB GENERATOR & SUBMITTER 
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| For each sample ID, this R script auto-generates a per-sample SLURM bash script
#| to run STAR-Fusion on paired-end FASTQ files, then submits the job via `sbatch`.
#| It activates the STAR-Fusion conda env, points to the CTAT genome library, and
#| writes logs and outputs to sample-specific folders.
#|
#| Workflow:
#|   1) Discover sample IDs from FASTQ names (R1/R2 pairs).
#|   2) Compose per-sample paths and output directories.
#|   3) Write a SLURM script with resources and STAR-Fusion command.
#|   4) Submit each job to the scheduler (`system("sbatch <file>")`).
#|
#| Inputs (configure below):
#|   - samples: vector of sample IDs inferred from *_2.fastq.gz.
#|   - FASTQs_directory: folder with trimmed paired FASTQs.
#|   - out_directory: destination for STAR-Fusion outputs.
#|   - genome_lib_dir: CTAT genome library path.
#|   - SLURM resources: memory, time, partition, CPUs, email.
#|
#| Outputs:
#|   - Per-sample SLURM script: 6_STAR-FUSION_<SAMPLE>.sh
#|   - Logs: logs/6_STAR-FUSION_<SAMPLE>.out/.err
#|   - STAR-Fusion output: <out_directory>/<SAMPLE>_starfusion_outdir/
################################################################################


################################################################################
###########################  DATA CONFIGURATION  ###############################
################################################################################

samples <- list.files(path = "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs", pattern = "*_2.fastq")
samples <- gsub("_2.fastq.gz", "", samples)

FATSQs_directory <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/"
out_directory <- "/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/"
genome_lib_dir <- "/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/Scripts/ctat_genome_lib_build_dir"

memory <- "50GB"
time <-"120:00:00"
partition <- "NORMAL"
cpu <- "16"
mail <- "irondon@cicbiogune.es"
################################################################################


################################################################################
###############################  JOBS LOOP  ####################################
################################################################################

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste("6_STAR-FUSION_", samples[s], ".sh", sep = "") 

  # Forward strand
  input_fastq_read1 <- paste(FATSQs_directory, samples[s], "_1.fastq.gz", sep = "")

  # Reverse strand
  input_fastq_read2 <- paste(FATSQs_directory, samples[s], "_2.fastq.gz", sep = "")
  
  # Output bam
  output_dir <- paste(out_directory, samples[s], "_starfusion_outdir", sep = "") 

  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/6_STAR-FUSION_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/6_STAR-FUSION_", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    paste("#SBATCH --mail-user=",mail,sep = ""),
    "#SBATCH --mail-type=FAIL",

    "\n\n",
    "source /opt/ohpc/pub/apps/anaconda3/cic-env", 
    "conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/star_fusion", 
    "export LD_LIBRARY_PATH=/share/apps/star/STAR-2.7.7a/bin/Linux_x86_64/",
    "\n",

    paste("STAR-Fusion --left_fq ", input_fastq_read1, " --right_fq ", input_fastq_read2," --genome_lib_dir ", genome_lib_dir, " --output_dir ", output_dir, " --CPU ",cpu, sep = ""),

    
    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}