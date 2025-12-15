# Over-Representation Analysis (ORA)

library(fgsea)
library(msigdbr)
library(dplyr)
library(tidyr)
library(tibble)
library(orthogene)
library(readxl)


## ORA of Top 500 genes ranked by JSD in aging and OSKM-mediated rejuvenation

# Uploading gene sets from C2 Collection (Curated Gene Sets) in the MSigDB database
msigdb_c2_cgp <- msigdbr(species = "Homo sapiens", category = "C2")
pathways_c2_cgp <- split(msigdb_c2_cgp$gene_symbol, msigdb_c2_cgp$gs_name)
gene_sets <- pathways_c2_cgp

# Reading list of genes ranked by JSD in 'Old vs Young' and 'Old+OSKM vs Old' comparisons
df <- read_excel("DatasetEV1.xlsx", skip = 2)
#df <- read_excel("DatasetEV2.xlsx", skip = 2)

# Removing NA and duplicated genes 
df_clean <- df %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)
dim(df_clean)

# Converting ortholog genes (mouse to human)
gene_df_h <- convert_orthologs(
  gene_df = df_clean,
  gene_input = "Gene",
  input_species = "mouse"
)
gene_df_h <- as.data.frame(gene_df_h) %>%
  rownames_to_column(var = "gene_human")

# Getting the top 500 genes based on JSD
sig_genes <- gene_df_h %>%
  filter(as.numeric(`RANK comp.`) <= 500) %>%
  pull(gene_human)

# Defining the genes in the Universe/background as all genes taken into account in the JSD calculation and ranking
universe_genes <- gene_df_h$gene_human

# Performing ORA against our gene sets of interest
foraRes <- fora(genes=sig_genes, universe=universe_genes, pathways=gene_sets)
foraRes


##############################################################
##############################################################

## ORA of Differentially Expressed Genes in aging and OSKM-mediated rejuvenation

# Uploading gene sets from C2 Collection (Curated Gene Sets) in the MSigDB database
msigdb_c2_cgp <- msigdbr(species = "Homo sapiens", category = "C2")
pathways_c2_cgp <- split(msigdb_c2_cgp$gene_symbol, msigdb_c2_cgp$gs_name)
pathways_c2_cgp

# Subsetting collection to retain only the top 5 gene sets that were enriched with top 500 genes ranked based on JSD
keep_sets <- c(
  "BENPORATH_ES_WITH_H3K27ME3",
  "BENPORATH_EED_TARGETS",
  "BENPORATH_SUZ12_TARGETS",
  "BENPORATH_PRC2_TARGETS",
  "MIKKELSEN_MEF_HCP_WITH_H3K27ME3"
)
subset_gene_sets <- pathways_c2_cgp[keep_sets]


# Reading results of differential expression analysis in 'Old vs Young' and 'Old+OSKM vs Old' comparisons
df <- read_excel("DatasetEV3.xlsx", skip = 2)
#df <- read_excel("DatasetEV4.xlsx", skip = 2)

# Remove NA and duplicated genes before converting orthologs
df_clean <- df %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)
dim(df_clean)

gene_df_h <- convert_orthologs(
  gene_df = df_clean,
  gene_input = "symbol",
  input_species = "mouse"
)

# Get differentially expressed genes. Filter rows with padj < 0.05 and abs(log2FC) > 1.5
sig_genes <- rownames(gene_df_h)[
  !is.na(gene_df_h$padj) &
    !is.na(gene_df_h$log2FoldChange) &
    gene_df_h$padj < 0.05 &
    abs(gene_df_h$log2FoldChange) > 1.5
]

# Defining the genes in the Universe/background as all genes captured by RNA-Seq and taken into account in the differential expression analysis
universe_genes <- rownames(gene_df_h)

foraRes <- fora(genes=sig_genes, universe=universe_genes, pathways=subset_gene_sets)
foraRes
