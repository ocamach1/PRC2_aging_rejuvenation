# Metaplots of JSD and H3K27me3 occupancy over genes

# Example below: Metaplot of JSD and H3K27me3 over differentially methylated genes (based on JSD magnitude) during aging in whole skin.
# JSD corresponds to comparisons of all old samples against one representative young sample in whole skin.
# H3K27me3 occupancy is measured in old and young epidermis.

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)
library(EnrichedHeatmap)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Importing JSD tracks:
JSD_1 <- import("/Users/oscarcamacho/Documents/Feinberg_lab/skin_aging/bwig_files/2025_bwig_files/JSD-5760-con-VS-Y4F-4.bw")
genome(JSD_1) <- "mm10"
colnames(mcols(JSD_1))[1] <- "JSD"
JSD_1$phenoname <- "old-1"

JSD_2 <- import("/Users/oscarcamacho/Documents/Feinberg_lab/skin_aging/bwig_files/2025_bwig_files/JSD-5761-con-VS-Y4F-4.bw")
genome(JSD_2) <- "mm10"
colnames(mcols(JSD_2))[1] <- "JSD"
JSD_2$phenoname <- "old-2"

JSD_3 <- import("/Users/oscarcamacho/Documents/Feinberg_lab/skin_aging/bwig_files/2025_bwig_files/JSD-5762-con-VS-Y4F-4.bw")
genome(JSD_3) <- "mm10"
colnames(mcols(JSD_3))[1] <- "JSD"
JSD_3$phenoname <- "old-3"

JSD_4 <- import("/Users/oscarcamacho/Documents/Feinberg_lab/skin_aging/bwig_files/2025_bwig_files/JSD-5764-con-VS-Y4F-4.bw")
genome(JSD_4) <- "mm10"
colnames(mcols(JSD_4))[1] <- "JSD"
JSD_4$phenoname <- "old-4"

JSD_list <- list(
  old1 = JSD_1,
  old2 = JSD_2,
  old3 = JSD_3,
  old4 = JSD_4)


# Importing H3K27me3 occupancy tracks:
K27_Y1 <- rtracklayer::import("Y1-K27_minus_input_scaled.bw")
genome(K27_Y1) <- "mm10"
colnames(mcols(K27_Y1))[1] <- "K27"
K27_Y1$phenoname <- "K27_Y1"
K27_Y1$group     <- "K27_Y"

K27_Y2 <- rtracklayer::import("Y3-K27_minus_input_scaled.bw")
genome(K27_Y2) <- "mm10"
colnames(mcols(K27_Y2))[1] <- "K27"
K27_Y2$phenoname <- "K27_Y2"
K27_Y2$group     <- "K27_Y"

K27_O1 <- rtracklayer::import("O1-K27_minus_input_scaled.bw")
genome(K27_O1) <- "mm10"
colnames(mcols(K27_O1))[1] <- "K27"
K27_O1$phenoname <- "K27_O1"
K27_O1$group     <- "K27_O"

K27_O2 <- rtracklayer::import("O2-K27_minus_input_scaled.bw")
genome(K27_O2) <- "mm10"
colnames(mcols(K27_O2))[1] <- "K27"
K27_O2$phenoname <- "K27_O2"
K27_O2$group     <- "K27_O"

K27_O3 <- rtracklayer::import("O3-K27_minus_input_scaled.bw")
genome(K27_O3) <- "mm10"
colnames(mcols(K27_O3))[1] <- "K27"
K27_O3$phenoname <- "K27_O3"
K27_O3$group     <- "K27_O"

K27_list <- list(
  K27_Y1 = K27_Y1,
  K27_Y2 = K27_Y2,
  K27_O1 = K27_O1,
  K27_O2 = K27_O2,
  K27_O3 = K27_O3
)

# Upload list of differentially methylated genes by JSD
JSD_df <- read_excel("DatasetEV1.xlsx", skip = 2)
JSD_df <- JSD_df[JSD_df$q-value comp. <= 0.05, ]
gene_list <- JSD_df$Gene
gene_list <- gene_set1[!is.na(gene_set1)]

# Obtaining gene body coordinates for differentially methylated gene set (scaled with ±2 kb flanks)
txdb    <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)
genome(genes_gr) <- "mm10"
autosomes <- paste0("chr", 1:19)
genes_gr  <- keepSeqlevels(genes_gr, autosomes, pruning.mode = "coarse")
seqlevelsStyle(genes_gr) <- seqlevelsStyle(JSD_1)
entrez_ids <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = gene_list,   # your chosen set (genes_neither, etc.)
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
entrez_ids <- unique(na.omit(entrez_ids))
genes_gr   <- genes_gr[names(genes_gr) %in% entrez_ids]

# Build metagene matrices
bin_size <- 50  # 50 bp bins

# JSD matrices (one per sample)
mat_JSD_list <- lapply(JSD_list, function(gr) {
  normalizeToMatrix(
    signal       = gr,
    target       = genes_gr,
    value_column = "JSD",
    extend       = c(2000, 2000),  # -2 kb / +2 kb
    w            = bin_size,
    mean_mode    = "w0",
    smooth       = FALSE,
    empty_value  = NA
  )
})

# H3K27me3 matrices (one per H3K27 sample)
mat_K27_list <- lapply(K27_list, function(gr) {
  gr <- keepSeqlevels(gr, autosomes, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- seqlevelsStyle(genes_gr)
  
  normalizeToMatrix(
    signal       = gr,
    target       = genes_gr,
    value_column = "K27",
    extend       = c(2000, 2000),
    w            = bin_size,
    mean_mode    = "w0",
    smooth       = FALSE,
    empty_value  = NA
  )
})

n_JSD  <- ncol(mat_JSD_list[[1]])
n_K27  <- ncol(mat_K27_list[[1]])
common_bins <- min(n_JSD, n_K27)
mat_JSD_list <- lapply(mat_JSD_list, function(m) m[, seq_len(common_bins)])
mat_K27_list <- lapply(mat_K27_list, function(m) m[, seq_len(common_bins)])

# Define metagene coordinates (upstream, scaled body, downstream)
up_bp    <- 2000
down_bp  <- 2000
up_bins   <- up_bp   / bin_size
down_bins <- down_bp / bin_size
total_bins <- common_bins
body_bins  <- total_bins - up_bins - down_bins
L_up   <- up_bp
L_body <- 4000   # just for labeling the scaled body region
L_down <- down_bp
pos_scaled <- numeric(total_bins)
pos_scaled[1:up_bins] <- seq(-L_up, 0, length.out = up_bins)
pos_scaled[(up_bins + 1):(up_bins + body_bins)] <-
  seq(0, L_body, length.out = body_bins)
pos_scaled[(up_bins + body_bins + 1):total_bins] <-
  seq(L_body, L_body + L_down, length.out = down_bins)

# JSD profiles
JSD_profiles <- lapply(mat_JSD_list, function(m) colMeans(m, na.rm = TRUE))
jsd_df <- purrr::imap_dfr(JSD_profiles, ~ {
  data.frame(
    pos    = pos_scaled,
    value  = .x,
    sample = .y,
    signal = "JSD",
    group  = "Old"
  )
})

jsd_mean_df <- jsd_df %>%
  dplyr::group_by(pos) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(signal = "JSD")

# H3K27me3 profiles
k27_profiles <- purrr::imap(mat_K27_list, ~ colMeans(.x, na.rm = TRUE))
k27_sample_info <- data.frame(
  sample = names(K27_list),
  group  = c("Young", "Young", "Old", "Old", "Old"),
  stringsAsFactors = FALSE
)
k27_df <- purrr::imap_dfr(k27_profiles, ~ {
  data.frame(
    pos    = pos_scaled,
    value  = .x,
    sample = .y,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::left_join(k27_sample_info, by = "sample") %>%
  dplyr::mutate(signal = "H3K27me3")

k27_mean_df <- k27_df %>%
  dplyr::group_by(pos, group) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(signal = "H3K27me3")

# Plot of JSD and H3K27me3 over scaled gene bodies and 2 kb upstream/downstream
k27_colors <- c(
  Young = "#1f78b4",  # blue
  Old   = "#e31a1c"   # red
)

p <- ggplot() +
  ## shaded upstream / body / downstream (shown in both facets)
  annotate("rect",
           xmin = -L_up, xmax = 0,
           ymin = -Inf,  ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect",
           xmin = 0, xmax = L_body,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  annotate("rect",
           xmin = L_body, xmax = L_body + L_down,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(
    data = jsd_df,
    aes(x = pos, y = value, group = sample),
    colour = "grey70",
    size   = 0.4,
    alpha  = 0.8
  ) +
  geom_line(
    data = jsd_mean_df,
    aes(x = pos, y = value),
    colour = "black",
    size   = 0.9
  ) +
  geom_line(
    data = k27_df,
    aes(x = pos, y = value, group = sample, colour = group),
    size   = 0.4,
    alpha  = 0.7
  ) +
  geom_line(
    data = k27_mean_df,
    aes(x = pos, y = value, colour = group),
    size = 1.0
  ) +
  geom_vline(xintercept = 0,      linetype = "dashed") +
  geom_vline(xintercept = L_body, linetype = "dashed") +
  facet_wrap(~ signal, ncol = 1, scales = "free_y") +
  scale_colour_manual(
    values = k27_colors,
    name   = "H3K27me3 group"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  
  # region labels at the top of each panel
  annotate("text", x = -L_up/2,          y = Inf,
           label = "-2 kb upstream",     vjust = 1.3, size = 3.2) +
  annotate("text", x = L_body/2,         y = Inf,
           label = "Gene body (scaled)", vjust = 1.3, size = 3.2) +
  annotate("text", x = L_body + L_down/2, y = Inf,
           label = "+2 kb downstream",   vjust = 1.3, size = 3.2) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.spacing    = unit(0.5, "lines"),
    strip.background = element_rect(fill = "grey95"),
    plot.margin      = margin(5, 10, 5, 10),
    axis.title.x     = element_text(margin = margin(t = 8)),
    axis.title.y     = element_text(margin = margin(r = 8)),
    legend.position  = "top"
  ) +
  labs(
    x = "Metagene coordinate (scaled: −2 kb, gene body expanded, +2 kb)",
    y = "Signal"
  )
