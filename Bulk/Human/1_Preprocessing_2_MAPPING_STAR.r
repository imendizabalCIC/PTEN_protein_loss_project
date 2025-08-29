#/opt/R/R-4.1.2/bin/R

###############################################################################
#| THIS CODES CREATES A SCRIPT FOR EACH FILE FASTA AND SUMIT JOBS 
###############################################################################
#| Date: 12/12/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This scripts creates a .sh file and then is launched to the job executor in 
#| indar to mapping the RNAseq of 197 patients with prostate cancer that underwent
#| surgery for prostate extraction. Note that the files selected for the mapping 
#| with STAR contains the Forward (*_1 files) and Reverse (*_2 files). 
#|
#|     STAR command line:
#|  STAR --genomeDir genomedir [options] file1 file2 --outFileNamePrefix outdir
#|
#| Parameters used here:
#|    1) --genomeDir: specifies path to the directory where the genome indices 
#|       are stored.
#|    2) --runThreadN: this option defines the number of threads to be used for 
#|       genome generation, it has to be set to the number of available cores on 
#|       the server node. (max 16 with Indar and 36 with Nextera)
#|    3) --readFilesIn: name(s) (with path) of the files containing the sequences 
#|       to be mapped (e.g. RNA-seq FASTQ files). If using Illumina paired-end 
#|       reads, the read1 and read2 files have to be supplied. 
#|    4) --readFilesCommand: UncompressionCommand option. (set gunzip -c)
#|    5) --outFilterMultimapNmax: max number of multiple alignments allowed for 
#|       a read: if exceeded, the read is considered unmapped. Default 20.
#|       COMMENTS FROM THE WEB: --outFilterMultimapNmax 1 limits the output to 
#|       just the uniquely-mapping reads (i.e. reads that confidently map to only 
#|       one locus). If the genes of interest have high sequence similarity, this
#|       option will probably eliminate a large number of reads that map to multiple 
#|       loci. Working with multi-mappers, on the other hand, requires a strategy 
#|       to assign them to specific loci. Note that: if the GeneCounts option is
#|       used, STAR counts only uniquely mapping reads, without regard to the 
#|       --outFilterMultimapNmax value.
#|    6) --outReadsUnmapped: string: output of unmapped and partially mapped 
#|      (i.e. mapped only one mate of a paired end read) reads in separate file(s).
#|      Default None.
#|          * None: no output
#|          * Fastx: output in separate fasta/fastq files, Unmapped.out.mate1/2   
#|    7) --outSAMtype: Default: SAM. Strings: type of SAM/BAM output
#|           1st word:
#|             * BAM: output BAM without sorting
#|             * SAM: output SAM without sorting
#|             * None: no SAM/BAM output
#|           2nd, 3rd:
#|             * Unsorted: standard unsorted
#|             * SortedByCoordinate: sorted by coordinate. This option will allocate 
#|               extra memory for sorting which can be specified by -limitBAMsortRAM.
#|    8) --twopassMode: String 2-pass mapping mode. Default None. Annotated junctions 
#|       will be included in both the 1st and 2nd passes. To run STAR 2-pass mapping 
#|       for each sample separately, use --twopassMode Basic option. STAR will perform 
#|       the 1st pass mapping, then it will automatically extract junctions, insert 
#|       them into the genome index, and, finally, re-map all reads in the 2nd mapping 
#|       pass. 
#|           * None: 1-pass mapping
#|           * Basic: 2-pass mapping, with all 1st pass junction inserted into the 
#|             genome indices on the fly  
#|    9) --limitBAMsortRAM: int>=0: maximum available RAM (bytes) for sorting BAM. 
#|       If =0, it will be set to the genome index size. 0 value can only be used with 
#|      -genomeLoad NoSharedMemory option. Defaul 0.
#|    10) --quantMode: string(s): types of quantification requested.
#|           * -: None
#|           * TranscriptomicSAM: output SAM/BAM alignments to transcriptome into a 
#|             separate file
#|           * GeneCounts: count reads per gene
################################################################################


################################################################################
#|  DATA LIST  
################################################################################

#| List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = "/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs", pattern = "*_1.fastq")
samples <- gsub("_1.fastq.gz", "", samples)

################################################################################


################################################################################
#|  PARAMETERS  
################################################################################

#| Parameters used by SLURM to assign a node to the job
cpu <- "8"
memory <- "40G"
time <- "10:00:00"
partition <- "FAST"

#| Genome directory
genomedir <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes/Indexes/Human_151"   
################################################################################

################################################################################
#|  JOBS LOOP  
################################################################################

for(s in 1:length(samples)){

  #| Name of the file (regarding the sample variable)
  filename <- paste("STAR_", samples[s], ".sh", sep = "") 

  #| Forward strand
  filedir1 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_1.fastq.gz", sep = "")

  #| Reverse strand
  filedir2 <- paste("/vols/GPArkaitz_bigdata/DATA_shared/AC-45_RNAseq-FFPE/FASTQs_trimmed/", samples[s], "_2.fastq.gz", sep = "")
  
  #| Out directory
  outdir <- paste("/vols/GPArkaitz_bigdata/irondon/Project_AC-45_RNAseq-FFPE/RNAseq/2_STEP_STAR/Outdir_STAR/", samples[s], "_STAR", sep = "") 
  
  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o", samples[s], ".out", sep = ""),
    paste("#SBATCH -e", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    "#SBATCH --exclude=gn[05-09]",
    "#SBATCH --mail-user=irondon@cicbiogune.es",
    "#SBATCH --mail-type=FAIL",
    "\n\n",
    "export LD_LIBRARY_PATH=/share/apps/star/STAR-2.7.7a/bin/Linux_x86_64/", 
    "\n",
    paste("STAR --genomeDir", genomedir, "--runThreadN", cpu, "--readFilesIn", filedir1, filedir2, "--readFilesCommand gunzip -c --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --limitBAMsortRAM 2000000000 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix", outdir, sep = " "),
    file = filename, sep = "\n", append = FALSE)
  system(paste("sbatch", filename, sep = " "))

}
################################################################################