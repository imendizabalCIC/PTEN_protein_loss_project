# Exploring PTEN Loss in Prostate Cancer Through Bulk and Single-Cell Transcriptome Analysis

Analysis of the bulk RNA-seq data for Human and Mice are inside the folder **"Bulk"**, where scripts labeled with:

1. **"1_Preprocessing_*"**  
   This step includes:  
   - Adapter trimming of raw data using **Cutadapt**.  
   - Quality control with **FastQC**.  
   - Alignment of FASTQ files using **STAR**.  
   - Generation of raw counts from STAR output.

2. **"2_PTEN_Quantification.R"**  
   - Computation of the H-score based on immunohistochemistry assessment of PTEN protein.  
   - Characterization of PTEN protein loss samples based on activation of the PI3K-AKT-mTOR pathway.

3. **"3_Data_Analysis_*"**  
   - Differential expression analysis of PTEN protein loss versus presence using **Limma Voom**.  
   - Network analysis implementation with **WGCNA**.

4. **"4_Enrichment_*"**  
   - Enrichment analysis of bulk RNA-seq data using **GSEA** and **g:Profiler** on genes significantly associated with PTEN loss.

5. **"5_Celltype_xCell.R"** (Human)  
   - **xCell**-based deconvolution analysis of bulk RNA-seq data to infer 64 distinct cell types in clinical samples.
