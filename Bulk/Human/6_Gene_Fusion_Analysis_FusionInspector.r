#/opt/R/R-4.1.2/bin/R


################################################################################
#|  BATCH FUSIONINSPECTOR JOB GENERATOR & SUBMITTER ###########
################################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| For each sample ID, this R script auto-generates a per-sample SLURM bash script
#| to run FusionInspector on paired-end FASTQ files, then submits the job with
#| `sbatch`. It uses the STAR-Fusion results (abridged predictions) as input and
#| writes logs and outputs to sample-specific folders.
#|
#| Workflow:
#|   1) Read the list of sample IDs (vector `samples`).
#|   2) Compose paths to trimmed FASTQs (_1/_2), STAR-Fusion output dir, and
#|      CTAT genome lib (`genome_lib_dir`).
#|   3) Write a per-sample bash script with SBATCH resources, env activation
#|      (conda), and `FusionInspector` command (with `--vis`).
#|   4) Submit each script to the scheduler via `system("sbatch <file>")`.
#|
#| Inputs (configure in "DATA CONFIGURATION"):
#|   - samples: vector of sample IDs (e.g., "AC137").
#|   - FATSQs_directory: folder with paired FASTQs (trimmed).
#|   - out_directory: where STAR-Fusion and FusionInspector outputs live.
#|   - genome_lib_dir: CTAT genome lib path.
#|   - SLURM resources: memory, time, partition, CPUs, email.
#|
#| Outputs:
#|   - Per-sample SLURM script: 6_FUSION-INSPECTOR_<SAMPLE>.sh
#|   - Logs: logs/6_FUSION-INSPECTOR_<SAMPLE>.out/.err
#|   - FusionInspector output: <out_directory>/<SAMPLE>_fusioninspector_outdir/
################################################################################


################################################################################
#|  DATA CONFIGURATION  
################################################################################

#| Sample ID
samples <- list.files(path = "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs", pattern = "*_2.fastq")
samples <- gsub("_2.fastq.gz", "", samples)

#| Specifying directories
FATSQs_directory <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/"
out_directory <- "/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/"
genome_lib_dir <- "/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/Scripts/ctat_genome_lib_build_dir"

#| SLURM parameters
memory <- "50GB"
time <-"12:00:00"
partition <- "FAST"
cpu <- "12"
mail <- "irondon@cicbiogune.es"
################################################################################


################################################################################
#|  JOBS LOOP  
################################################################################

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste("6_FUSION-INSPECTOR_", samples[s], ".sh", sep = "") 

  # Forward strand
  input_fastq_read1 <- paste(FATSQs_directory, samples[s], "_1.fastq.gz", sep = "")

  # Reverse strand
  input_fastq_read2 <- paste(FATSQs_directory, samples[s], "_2.fastq.gz", sep = "")
  
  # Output bam
  output_dir <- paste(out_directory, samples[s], "_starfusion_outdir", sep = "") 

  # Output Fusion_inspector
  output_dir_fusioninspector <- paste(out_directory, samples[s], "_fusioninspector_outdir", sep = "") 

  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=FUSION-INSPECTOR", samples[s], sep = ""),
    paste("#SBATCH -o logs/6_FUSION-INSPECTOR_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/6_FUSION-INSPECTOR_", samples[s], ".err", sep = ""),
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

    paste("FusionInspector --left_fq ", input_fastq_read1, 
    " --right_fq ", input_fastq_read2,
    " --genome_lib_dir ", genome_lib_dir, 
    " --fusions ", paste0(output_dir, "/star-fusion.fusion_predictions.abridged.tsv"),
    " --output_dir ", output_dir_fusioninspector, 
    " --CPU ",cpu, 
    " --vis", sep = ""),

    
    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}