# Transcriptional network analysis of PTEN-protein-deficient prostate tumors reveals robust stromal reprogramming and signs of senescent paracrine communication

Analysis of the bulk RNA-seq data for Human and Mice are inside the folder **"Bulk"**, where scripts labeled with:

1. **"1_Preprocessing_*"**  
   - Adapter trimming of raw data using **Cutadapt**.  
   - Quality control with **FastQC**.  
   - Alignment of FASTQ files using **STAR**.  
   - Generation of raw counts from STAR output.

2. **"2_PTEN_Quantification.R"**  
   - Computation of the H-score based on immunohistochemistry assessment of PTEN protein.  
   - Characterization of PTEN protein loss samples based on activation of the PI3K-AKT-mTOR pathway.

3. **"3_Data_Analysis_*"**  
   - Differential expression analysis of PTEN protein loss versus presence using **Limma Voom** package.  
   - Network analysis implementation with **WGCNA**.

4. **"4_Enrichment_*"**  
   - Enrichment analysis of bulk RNA-seq data using **GSEA** and **g:Profiler** on genes significantly associated with PTEN loss.

5. **"5_Celltype_xCell.R"** (Human)  
   - **xCell**-based deconvolution analysis of bulk RNA-seq data to infer 64 distinct cell types in the clinical samples.

6. **"6_Gene_Fusion_Analysis_"** (Human)  
   - Identification of candidate gene fusions from RNA-Seq data using STAR-Fusion, with downstream quality assessment, refinement, and correction using FusionInspector to reduce false positives.
  
7. **"7_Clinical_implications_PTEN_loss_"** (Human)  
   - Evaluates the prognostic and clinical impact of PTEN loss using green module expression scores in prostate cancer cohorts (Basurto and TCGA). 

The analysis scripts for single-cell RNA-seq data are located in the "Single-cell" folder, organized as follows:

1. **"0-2_scPipeline_*"**
   - Preprocessing, integration, and marker identification of a single-cell data of 8 localized prostate cancer samples from Chen et al. (Nature, 2021) using the **Seurat** package.

2. **"3_sc_UCell_Modules_Human_Mouse.R"**
   - Projection of the differentially expressed genes from green and purple modules (human analysis) and the yellow module (mouse dataset) using the **UCell** package.

3. **"4_sc_MuSiC.R"**
   - Prostate-specific cell type deconvolution analysis of the clinical and murine data with a single-cell data from 8 localized prostate cancer samples from Chen et al. (Nature, 2021) using the **MuSiC** package.

For details on package versions, please refer to the manuscript.

