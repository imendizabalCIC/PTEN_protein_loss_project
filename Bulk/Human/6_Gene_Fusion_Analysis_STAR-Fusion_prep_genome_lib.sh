#!/bin/bash

################################################################################
#|  GENOME CONSTRUCTION FOR FUSION INFERENCE
################################################################################
#| Date: 05/05/2025
#| Author: Ivana Rondon Lorefice
#|
#| Description: This SLURM job script prepares a STAR-Fusion compatible genome 
#| library using the CTAT genome lib builder. It takes the GRCh38 reference genome, 
#| GTF annotations, and fusion annotation library to generate the required 
#| genome indices and databases (Pfam, Dfam) for downstream fusion detection 
#| with STAR-Fusion. The script runs with 16 CPUs and 50 GB of memory on the 
#| NORMAL partition, and logs output/error files for tracking.
################################################################################

#SBATCH --job-name=prep_genome_lib
#SBATCH -o logs/prep_genome_lib_STAR-fusion.out
#SBATCH -e logs/prep_genome_lib_STAR-fusion.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=120:00:00
#SBATCH --partition=NORMAL
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

#| Loading conda environment: initialize Anaconda module
source /opt/ohpc/pub/apps/anaconda3/cic-env

#| Activating STAR-Fusion conda environment
conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/star_fusion 

#| Setting STAR binary path (needed for CTAT genome lib builder)
export LD_LIBRARY_PATH=/share/apps/star/STAR-2.7.7a/bin/Linux_x86_64/

#| Running CTAT genome lib builder:
#|  - --genome_fa: GRCh38 reference FASTA file
#|  - --gtf: filtered GTF annotations (Ensembl GRCh38.94)
#|  - --fusion_annot_lib: curated fusion annotation library (CTAT HumanFusionLib)
#|  - --CPU: number of threads to use (16 in this case)
#|  - --dfam_db: Dfam database for repeat masking (set to human)
#|  - --pfam_db: Pfam database for protein domain annotation (set to current)
#| This step generates the genome resource library required by STAR-Fusion.
/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/ctat-genome-lib-builder/prep_genome_lib.pl --genome_fa /vols/GPArkaitz_bigdata/DATA_shared/Genomes/Sequences/GRCh38/GRCh38_r94.all.fa \
  --gtf /vols/GPArkaitz_bigdata/DATA_shared/Genomes/Annotations/Homo_sapiens.GRCh38.94.filtered.gtf \
  --fusion_annot_lib /vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/6_STEP_STAR_GENE_FUSION/CTAT_HumanFusionLib.v0.1.0.dat.gz \
  --CPU 16 \
  --dfam_db human  \
  --pfam_db current


