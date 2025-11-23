library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(biomaRt)
library(org.Mm.eg.db)
library(apeglm)

# Upload transcript quantification files (output from Salmon aligner)
file_dirs <- c("Y_1","Y_2","Y_3","Y_4","O_1", "O_2", "O_3", "O_4", "O_5","O_OSKM_1", "O_OSKM_2", "O_OSKM_3", "O_OSKM_4", "O_OSKM_5")
files <- file.path("/RNA-seq/analysis/Old+OSKM_vs_Old/quants", file_dirs, "quant.sf")
names(files) <- c("Y_1","Y_2","Y_3","Y_4","O_1", "O_2", "O_3", "O_4", "O_5","O_OSKM_1", "O_OSKM_2", "O_OSKM_3", "O_OSKM_4", "O_OSKM_5")

# Import mouse transcriptome file (mm10 refrence genome) and merge transcript-level quantification (output from Salmon aligner) into gene-level quantification
txdb <- makeTxDbFromGFF("gencode.vM23.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
# Run tximport:
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Assign groups to samples
sampleTable <- data.frame(condition = factor(c("Y","Y","Y","Y","O","O","O","O","O","O_OSKM","O_OSKM","O_OSKM","O_OSKM","O_OSKM")))
rownames(sampleTable) <- colnames(txi$counts)

# Create the DESEQ object
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

# Generate the DEG results table, specifying the desired contrast
res_raw <- results(dds, contrast=c("condition","O_OSKM","O"),alpha=0.05)

# Applying shrinkage of fold changes
res <- lfcShrink(dds, coef = "condition_O_OSKM_vs_O", type = "apeglm", res=res_raw)
res <- res[order(res$padj), ]
res <- res[!is.na(res$padj),]

# Convert ENSEMBL IDs to gene symbols
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
ens_ids <- gsub("\\..*", "", rownames(res))
bm_res <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
                filters = "ensembl_gene_id",
                values = ens_ids,
                mart = mart)
res$ensembl_gene_id <- ens_ids
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- gsub("\\..*", "", rownames(res_df))
res_annotated <- merge(res_df, bm_res, by = "ensembl_gene_id", all.x = TRUE)
res_annotated <- res_annotated[!duplicated(res_annotated$ensembl_gene_id), ]
res_annotated <- res_annotated[order(res_annotated$padj), ]

write.csv(res_annotated, file='DatasetEV4.csv')
