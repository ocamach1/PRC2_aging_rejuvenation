library(tximport)
library(GenomicFeatures)
library(readr)
library('biomaRt')
library(DESeq2)
library("org.Mm.eg.db")

# Upload transcript quantification files (output from Salmon aligner)
file_dirs <- c("O+OSKM-1", "O+OSKM-2", "O+OSKM-3", "O+OSKM-4", "O+OSKM-5","O-1", "O-2", "O-3", "O-4")
files <- file.path("/RNA-seq/analysis/Old+OSKM_vs_Old/quants", file_dirs, "quant.sf")
names(files) <- c("O+OSKM-1", "O+OSKM-2", "O+OSKM-3", "O+OSKM-4", "O+OSKM-5","O-1", "O-2", "O-3", "O-4")

# Import mouse transcriptome file (mm10 refrence genome) and merge transcript-level quantification (output from Salmon aligner) into gene-level quantification
txdb <- makeTxDbFromGFF("gencode.vM23.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
# Run tximport:
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Assign groups to samples
sampleTable <- data.frame(condition = factor(c("O+OSKM","O+OSKM","O+OSKM","O+OSKM","O+OSKM","O","O","O","O")))
rownames(sampleTable) <- colnames(txi$counts)

# Create the DESEQ object
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Setting the desired condition as the reference
dds$condition <- relevel(dds$condition, ref = "O")

# Performing differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)                   
res <- res[order(res$padj), ]

# Add gene names to ENSEMBL IDs
tmp=gsub("\\..*","",row.names(res))
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=tmp,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

write.csv(res, file='DEG_Old+OSKM_vs_Old_results.csv')
