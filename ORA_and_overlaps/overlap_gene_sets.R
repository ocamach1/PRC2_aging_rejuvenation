# Overlap of gene sets
# Required packages:
# - TxDb.Mmusculus.UCSC.mm10.knownGene
# - hypeR

# Example: Overlapping differentially methylated genes (based on JSD) and genes with differential H3K27me3 during aging.

# Define the Universe for gene set 1 (in this case it is all genes that were part of the JSD analysis):
gene_set1_universe <- JSD_df$Gene
gene_set1_universe <- gene_set1_universe[!is.na(gene_set1_universe)]

# Upload differentially methylated genes (q-value < 0.05 from JSD analysis):
JSD_df <- read_excel("DatasetEV1.xlsx", skip = 2)
JSD_df <- JSD_df[JSD_df$q-value comp. <= 0.05, ]
gene_set1 <- JSD_df$Gene
gene_set1 <- gene_set1[!is.na(gene_set1)]


# Define the Universe for gene set 2 (in this case it is all genes included in the annotations from TxDb.Mmusculus.UCSC.mm10.knownGene, used to map the ChIP peaks):
# Note: since JSD analysis includes all genes in TxDb.Mmusculus.UCSC.mm10.knownGene except the sex chromosomes, genes in sex chromosomes are the only genes
# not shared between the Universes of both gene sets.
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
g <- genes(txdb) 
entrez <- as.character(mcols(g)$gene_id)
sym <- mapIds(org.Mm.eg.db,
              keys     = entrez,
              keytype  = "ENTREZID",
              column   = "SYMBOL",
              multiVals = "first")
mcols(g)$SYMBOL <- unname(sym[entrez])
gene_set2_universe <- g$SYMBOL

# Upload genes with differential H3K27me3 (FDR < 0.01 from differential analyisis with diffbind):
H3K27me3_diffbind <- read.csv("DatasetEV5.xlsx", skip =2)
H3K27me3_diffbind <- H3K27me3_diffbind[H3K27me3_diffbind$FDR <= 0.01, ]
# Get only genes that overlap with H3K27me3 mark (overlap = 0) or are 2000 bp from the gene body
H3K27me3_diffbind <- H3K27me3_diffbind[abs(H3K27me3_diffbind$distance) <= 2000, ]
gene_set2 <- H3K27me3_diffbind$Closest_gene
gene_set2 <- gene_set2[!is.na(gene_set2)]

# Get the intersection from both Universes and restrict each gene set to the common Universe:
U <- intersect(gene_set1_universe, gene_set2_universe)
gene_set1 <- intersect(gene_set1, U)
gene_set2 <- intersect(gene_set2, U)  


# Assess significance of the overlap with Fisher's exact test:
library(hypeR)
overlap <- intersect(gene_set1, gene_set2)
Overlap=length(overlap)
group1=length(gene_set1)
group2=length(gene_set2)
Total=length(U)
phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)

# Calculate odds ratio
A=Overlap
B=group1 - Overlap
C=group2 - Overlap
D=Total - group1 - group2 + Overlap
OR=(A*D)/(B*C)

