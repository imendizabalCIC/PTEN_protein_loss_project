# Exploring PTEN Loss in Prostate Cancer Through Bulk and Single-Cell Transcriptome Analysis

Analysis of the bulk RNA-seq data for Human and Mice are inside the folder "Bulk" where scripts labeled with 
  1. "1_Preprocessing_*"
     This step includes:
      a. Adapter trimming of raw data using Cutadapt.
      b. Quality control with FastQC.
      c. Alignment of FASTQ files using STAR.
      d. Generation of raw counts from STAR output.

  3. "2_PTEN_Quantification.R"
     Computation of the H-score based on immunohistochemistry assessment of PTEN protein aand characterization of PTEN protein loss samples based on activation of the PI3K-AKT-mTOR pathway.

  4. "3_Data_Analysis_*"
      a. Differential expression analysis of PTEN protein loss versus presence using Limma Voom.
      b. Network analysis implementation with WGCNA.

  5. "4_Enrichment_*"
     Enrichment analysis of bulk RNA-seq data using GSEA and g:Profiler on genes significantly associated with PTEN loss.

  6. "5_Celltype_xCell.R" (Human)
     xCell-based deconvolution analysis of bulk RNA-seq data to infer 64 distinct cell types in clinical samples. 
    
