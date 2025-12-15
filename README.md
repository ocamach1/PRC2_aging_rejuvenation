# Convergence of aging- and rejuvenation-related epigenetic alterations on PRC2 targets

This repository contains code used in the analysis of the data presented in our paper **Convergence of aging- and rejuvenation-related epigenetic alterations on PRC2 targets** (Oscar Camacho, Michael A. Koldobskiy, Pradeep Reddy, Atharv Oak, Yuxiang Sun, Kenna Sherman, Juan Carlos Izpisua Belmonte, Andrew P. Feinberg).

It includes scripts to:
- **Analyze WGBS data:**
  - Align and process WGBS data from raw fastq files.
  - Run informME pipeline: compute mean methylation level (MML), normalized entropy level (NME) for each WGBS sample, and Jensen-Shannon distance (JSD) between pairs of samples.
  - Rank epigenetically discordant genes when comparing samples (i.e. old vs young) based on JSD.
  - Perform cell type deconvolution with MeDeCom.
  - Perform Student's t-test and Wilcoxon test to compare distributions of MML and NME across groups.
    
- **Analyze ChIP-Seq data:**
  - Align and process ChIP-Seq data from raw fastq files.
  - Perform peak calling with Sicer2 (for H3K27me3).
  - Perform peak calling with EDD (for H3K9me2).
  - Perform differential binding analysis with DiffBind.

- **Analyze RNA-Seq data:**
  - Align and process RNA-Seq data from raw fastq files.
  - Perform differential expression analysis with DESeq2.
    
- **Perform over-representation analysis (ORA) and overlaps**

- **Plots**
