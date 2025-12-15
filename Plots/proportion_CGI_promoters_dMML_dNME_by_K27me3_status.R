# Plots showing the distribution of promoters with increased, unchanged or decreased MML/NME over CpG islands in aging/rejuvenation comparisons, 
# stratifying by H3K27me3 status on promoter (unmarked, unchanged during aging, loss during aging) (Figure EV3D)

# Increase or decrease in MML/NME was defined as showing a difference in condition-level averages exceeding + or - 5%, respectively.

library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(matrixStats)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(effsize)
library(tibble)
library(patchwork)

# Uploaed MML and NME tracks:
MML_1<-import("MML-Y4F-1.bw")
genome(MML_1)<-"mm10"
colnames(mcols(MML_1))[1]<-"MML"
MML_1$phenoname="Y-1"
MML_1$group="Y"

MML_2<-import("MML-Y4F-3.bw")
genome(MML_2)<-"mm10"
colnames(mcols(MML_2))[1]<-"MML"
MML_2$phenoname="Y-2"
MML_2$group="Y"

MML_3<-import("MML-Y4F-4.bw")
genome(MML_3)<-"mm10"
colnames(mcols(MML_3))[1]<-"MML"
MML_3$phenoname="Y-3"
MML_3$group="Y"

MML_4<-import("MML-5760-con.bw")
genome(MML_4)<-"mm10"
colnames(mcols(MML_4))[1]<-"MML"
MML_4$phenoname="O-1"
MML_4$group="O"

MML_5<-import("MML-5761-con.bw")
genome(MML_5)<-"mm10"
colnames(mcols(MML_5))[1]<-"MML"
MML_5$phenoname="O-2"
MML_5$group="O"

MML_6<-import("MML-5762-con.bw")
genome(MML_6)<-"mm10"
colnames(mcols(MML_6))[1]<-"MML"
MML_6$phenoname="O-3"
MML67$group="O"

MML_7<-import("MML-5764-con.bw")
genome(MML_7)<-"mm10"
colnames(mcols(MML_7))[1]<-"MML"
MML_7$phenoname="O-4"
MML_7$group="O"

MML_8<-import("MML-253-dox.bw")
genome(MML_8)<-"mm10"
colnames(mcols(MML_8))[1]<-"MML"
MML_8$phenoname="O+OSKM-1"
MML_8$group="O+OSKM"

MML_9<-import("MML-284-dox.bw")
genome(MML_9)<-"mm10"
colnames(mcols(MML_9))[1]<-"MML"
MML_9$phenoname="O+OSKM-2"
MML_9$group="O+OSKM"

MML_10<-import("MML-285-dox.bw")
genome(MML_10)<-"mm10"
colnames(mcols(MML_10))[1]<-"MML"
MML_10$phenoname="O+OSKM-3"
MML_10$group="O+OSKM"

MML_11<-import("MML-288-dox.bw")
genome(MML_11)<-"mm10"
colnames(mcols(MML_11))[1]<-"MML"
MML_11$phenoname="O+OSKM-4"
MML_11$group="O+OSKM"

MML_12<-import("MML-904-dox.bw")
genome(MML_12)<-"mm10"
colnames(mcols(MML_12))[1]<-"MML"
MML_12$phenoname="O+OSKM-5"
MML_12$group="O+OSKM"


NME_1<-import("NME-Y4F-1.bw")
genome(NME_1)<-"mm10"
colnames(mcols(NME_1))[1]<-"NME"
NME_1$phenoname="Y-1"
NME_1$group="Y"

NME_2<-import("NME-Y4F-3.bw")
genome(NME_3)<-"mm10"
colnames(mcols(NME_3))[1]<-"NME"
NME_3$phenoname="Y-2"
NME_3$group="Y"

NME_3<-import("NME-Y4F-4.bw")
genome(NME_3)<-"mm10"
colnames(mcols(NME_3))[1]<-"NME"
NME_3$phenoname="Y-3"
NME_3$group="Y"

NME_4<-import("NME-5760-con.bw")
genome(NME_4)<-"mm10"
colnames(mcols(NME_4))[1]<-"NME"
NME_4$phenoname="O-1"
NME_4$group="O"

NME_5<-import("NME-5761-con.bw")
genome(NME_5)<-"mm10"
colnames(mcols(NME_5))[1]<-"NME"
NME_5$phenoname="O-2"
NME_5$group="O"

NME_6<-import("NME-5762-con.bw")
genome(NME_6)<-"mm10"
colnames(mcols(NME_6))[1]<-"NME"
NME_6$phenoname="O-3"
NME_6$group="O"

NME_7<-import("NME-5764-con.bw")
genome(NME_7)<-"mm10"
colnames(mcols(NME_7))[1]<-"NME"
NME_7$phenoname="O-4"
NME_7$group="O"

NME_8<-import("NME-253-dox.bw")
genome(NME_8)<-"mm10"
colnames(mcols(NME_8))[1]<-"NME"
NME_8$phenoname="O+OSKM-1"
NME_8$group="O+OSKM"

NME_9<-import("NME-284-dox.bw")
genome(NME_9)<-"mm10"
colnames(mcols(NME_9))[1]<-"NME"
NME_9$phenoname="O+OSKM-2"
NME_9$group="O+OSKM"

NME_10<-import("NME-285-dox.bw")
genome(NME_10)<-"mm10"
colnames(mcols(NME_10))[1]<-"NME"
NME_10$phenoname="O+OSKM-3"
NME_10$group="O+OSKM"

NME_11<-import("NME-288-dox.bw")
genome(NME_11)<-"mm10"
colnames(mcols(NME_11))[1]<-"NME"
NME_11$phenoname="O+OSKM-4"
NME_11$group="O+OSKM"

NME_12<-import("NME-904-dox.bw")
genome(NME_12)<-"mm10"
colnames(mcols(NME_12))[1]<-"NME"
NME_12$phenoname="O+OSKM-5"
NME_12$group="O+OSKM"


MML_list <- list(
  Y_1      = MML_1,
  Y_2      = MML_3,
  Y_3      = MML_4,
  O_1      = MML_5,
  O_2      = MML_6,
  O_3      = MML_7,
  O_4      = MML_8,
  OOSKM_1  = MML_9,
  OOSKM_2  = MML_10,
  OOSKM_3  = MML_11,
  OOSKM_4  = MML_12,
  OOSKM_5  = MML_13
)

NME_list <- list(
  Y_1      = NME_1,
  Y_2      = NME_3,
  Y_3      = NME_4,
  O_1      = NME_5,
  O_2      = NME_6,
  O_3      = NME_7,
  O_4      = NME_8,
  OOSKM_1  = NME_9,
  OOSKM_2  = NME_10,
  OOSKM_3  = NME_11,
  OOSKM_4  = NME_12,
  OOSKM_5  = NME_13
)

sample_info <- tibble(
  sample = names(MML_list),
  group  = c(
    rep("Y", 3),
    rep("O", 4),
    rep("O+OSKM", 5)
  )
)

# Make sure genome and column names are set consistently
MML_list <- lapply(MML_list, function(gr) {
  genome(gr) <- "mm10"
  colnames(mcols(gr))[1] <- "MML"
  gr
})
NME_list <- lapply(NME_list, function(gr) {
  genome(gr) <- "mm10"
  colnames(mcols(gr))[1] <- "NME"
  gr
})

# Define CpG islands within promoters. Promoters are defined as ±2kb around TSS.
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)
genome(genes_gr) <- "mm10"
prom_gr <- promoters(genes_gr, upstream = 2000, downstream = 2000)
autosomes <- paste0("chr", 1:19)
prom_gr   <- keepSeqlevels(prom_gr, autosomes, pruning.mode = "coarse")
seqlevelsStyle(prom_gr) <- seqlevelsStyle(MML_1)
entrez_ids <- names(prom_gr)
sym <- mapIds(
  org.Mm.eg.db,
  keys      = entrez_ids,
  keytype   = "ENTREZID",
  column    = "SYMBOL",
  multiVals = "first"
)
prom_gr$symbol <- unname(sym[entrez_ids])
prom_gr <- prom_gr[!is.na(prom_gr$symbol)]
prom_gr$prom_id <- seq_len(length(prom_gr))

# Restrict promoters to CGI portions
cpg_islands <- import.bed(
  "CpG_islands_mm10.txt", # downloaded from UCSC genome browser
  trackLine = FALSE
)
genome(cpg_islands) <- "mm10"
seqlevelsStyle(cpg_islands) <- seqlevelsStyle(prom_gr)
cpg_islands <- keepSeqlevels(cpg_islands, autosomes, pruning.mode = "coarse")
ov_cgi <- findOverlaps(prom_gr, cpg_islands, ignore.strand = TRUE)
prom_cgi_gr <- pintersect(
  prom_gr[queryHits(ov_cgi)],
  cpg_islands[subjectHits(ov_cgi)],
  drop.nohit.ranges = TRUE
)
prom_cgi_gr$prom_id <- prom_gr$prom_id[queryHits(ov_cgi)]
prom_cgi_gr$symbol  <- prom_gr$symbol[queryHits(ov_cgi)]
prom_meta <- tibble(
  prom_id = prom_cgi_gr$prom_id,
  symbol  = prom_cgi_gr$symbol
) %>%
  distinct()

# Define H3K27me3 status on promoter
ChIP_promoters <- read.csv("DatasetEV5_only_promoters.xlsx, skip=2) # Dataframe mapping H3K27me3 marks to only promoters
prom_peak_df <- ChIP_promoters %>%
  mutate(
    is_loss          = (Fold < 0 & FDR < 0.01),
    is_gain          = (Fold > 0 & FDR < 0.01),
    has_H3K27me3_any = TRUE
  ) %>%
  rename(symbol = gene_symbol)

H3_class_df <- prom_peak_df %>%
  group_by(symbol) %>%
  summarise(
    any_loss = any(is_loss, na.rm = TRUE),
    any_gain = any(is_gain, na.rm = TRUE),
    any_mark = any(has_H3K27me3_any, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    H3_class = case_when(
      any_loss ~ "loss",
      any_mark ~ "marked_no_change",
      TRUE     ~ NA_character_
    )
  ) %>%
  select(symbol, H3_class)


## Calculate MML/NME on CpG islands within promoters
compute_promoter_MML <- function(MML_gr, samp_name) {
  ov <- findOverlaps(prom_cgi_gr, MML_gr, ignore.strand = TRUE)=
  if (length(ov) == 0) {
    return(
      tibble(
        prom_id  = prom_meta$prom_id,
        symbol   = prom_meta$symbol,
        sample   = samp_name,
        n_tiles  = 0L,
        MML_prom = NA_real_
      )
    )
  }
  ov_df <- tibble(
    prom_id = prom_cgi_gr$prom_id[queryHits(ov)],
    MML     = mcols(MML_gr)$MML[subjectHits(ov)]
  )
  summ <- ov_df %>%
    group_by(prom_id) %>%
    summarise(
      n_tiles  = sum(!is.na(MML)),
      MML_prom = ifelse(
        n_tiles >= 5,
        median(MML, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )
  prom_meta %>%
    left_join(summ, by = "prom_id") %>%
    mutate(
      n_tiles  = ifelse(is.na(n_tiles), 0L, n_tiles),
      MML_prom = ifelse(is.na(MML_prom), NA_real_, MML_prom),
      sample   = samp_name
    ) %>%
    dplyr::select(prom_id, symbol, sample, n_tiles, MML_prom)
}

compute_promoter_NME <- function(NME_gr, samp_name) {
  ov <- findOverlaps(prom_cgi_gr, NME_gr, ignore.strand = TRUE)
  if (length(ov) == 0) {
    return(
      tibble(
        prom_id  = prom_meta$prom_id,
        symbol   = prom_meta$symbol,
        sample   = samp_name,
        n_tiles  = 0L,
        NME_prom = NA_real_
      )
    )
  }
  ov_df <- tibble(
    prom_id = prom_cgi_gr$prom_id[queryHits(ov)],
    NME     = mcols(NME_gr)$NME[subjectHits(ov)]
  )
  summ <- ov_df %>%
    group_by(prom_id) %>%
    summarise(
      n_tiles  = sum(!is.na(NME)),
      NME_prom = ifelse(
        n_tiles >= 5,
        median(NME, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )
  
  prom_meta %>%
    left_join(summ, by = "prom_id") %>%
    mutate(
      n_tiles  = ifelse(is.na(n_tiles), 0L, n_tiles),
      NME_prom = ifelse(is.na(NME_prom), NA_real_, NME_prom),
      sample   = samp_name
    ) %>%
    dplyr::select(prom_id, symbol, sample, n_tiles, NME_prom)
}

# Compute condition-level mean MML averaging the sample-level MML across samples
# Sample-level MML is the average MML over CpG islands within promoters in a specific sample
prom_mml <- imap_dfr(MML_list, compute_promoter_MML) %>%
  left_join(sample_info, by = "sample")
prom_summary_MML <- prom_mml %>%
  group_by(prom_id, symbol, group) %>%
  summarise(
    MML_group_mean = mean(MML_prom, na.rm = TRUE),
    n_samples      = sum(!is.na(MML_prom)),
    .groups        = "drop"
  )
prom_wide_MML <- prom_summary_MML %>%
  dplyr::select(prom_id, symbol, group, MML_group_mean) %>%
  pivot_wider(
    names_from  = group,
    values_from = MML_group_mean
  ) %>%
  filter(!is.na(Y),
         !is.na(O),
         !is.na(`O+OSKM`)) %>%
  mutate(
    delta_O_vs_Y     = O - Y,
    delta_O_vs_OOSKM = O - `O+OSKM`
  ) %>%
  left_join(H3_class_df, by = "symbol") %>%
  mutate(
    H3_class = ifelse(is.na(H3_class), "unmarked", H3_class),
    H3_class = factor(H3_class,
                      levels = c("unmarked", "marked_no_change", "loss"))
  ) %>%
  filter(!is.na(H3_class))

# Compute condition-level mean NME averaging the sample-level NME across samples
# Sample-level NME is the average NME over CpG islands within promoters in a specific sample
prom_nme <- imap_dfr(NME_list, compute_promoter_NME) %>%
  left_join(sample_info, by = "sample")

prom_summary_NME <- prom_nme %>%
  group_by(prom_id, symbol, group) %>%
  summarise(
    NME_group_mean = mean(NME_prom, na.rm = TRUE),
    n_samples      = sum(!is.na(NME_prom)),
    .groups        = "drop"
  )

prom_wide_NME <- prom_summary_NME %>%
  dplyr::select(prom_id, symbol, group, NME_group_mean) %>%
  pivot_wider(
    names_from  = group,
    values_from = NME_group_mean
  ) %>%
  filter(!is.na(Y),
         !is.na(O),
         !is.na(`O+OSKM`)) %>%
  mutate(
    delta_O_vs_Y     = O - Y,
    delta_O_vs_OOSKM = O - `O+OSKM`
  ) %>%
  left_join(H3_class_df, by = "symbol") %>%
  mutate(
    H3_class = ifelse(is.na(H3_class), "unmarked", H3_class),
    H3_class = factor(H3_class,
                      levels = c("unmarked", "marked_no_change", "loss"))
  ) %>%
  filter(!is.na(H3_class))


Categorize no change, up down or no_change based on % difference in condition-level mean MML/NME
threshold_up   <- 0.05
threshold_down <- -0.05  
prom_wide_MML <- prom_wide_MML %>%
  mutate(
    # Old vs Young
    cat_OY_MML = case_when(
      delta_O_vs_Y > threshold_up    ~ "up",
      delta_O_vs_Y < threshold_down  ~ "down",
      TRUE                           ~ "no_change"
    ),
    # Old vs Old+OSKM
    cat_OO_MML = case_when(
      delta_O_vs_OOSKM > threshold_up   ~ "up",
      delta_O_vs_OOSKM < threshold_down ~ "down",
      TRUE                              ~ "no_change"
    ),
    # For over-representation (hyper = "up")
    hyper_OY_MML      = cat_OY_MML == "up",
    reversed_OSKM_MML = cat_OO_MML == "up"
  )

prom_wide_NME <- prom_wide_NME %>%
  mutate(
    cat_OY_NME = case_when(
      delta_O_vs_Y > threshold_up    ~ "up",
      delta_O_vs_Y < threshold_down  ~ "down",
      TRUE                           ~ "no_change"
    ),
    cat_OO_NME = case_when(
      delta_O_vs_OOSKM > threshold_up   ~ "up",
      delta_O_vs_OOSKM < threshold_down ~ "down",
      TRUE                              ~ "no_change"
    ),
    hyper_OY_NME      = cat_OY_NME == "up",
    reversed_OSKM_NME = cat_OO_NME == "up"
  )

# 4-panel figure showing no change, up or down in difference in condition-level mean MML/NME in Aging and Rejuvenation comparisons
dir_levels <- c("down", "no_change", "up")
dir_colors <- c(
  down      = "#4C72B0",   # blue
  no_change = "grey85",
  up        = "#DD8452"    # orange/red
)
# MML: Old vs Young:
summary_MML_OY <- prom_wide_MML %>%
  mutate(cat_OY_MML = factor(cat_OY_MML, levels = dir_levels)) %>%
  group_by(H3_class, cat_OY_MML) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(H3_class) %>%
  mutate(
    n_total = sum(n),
    prop    = 100 * n / n_total
  ) %>%
  ungroup()

p_MML_OY <- ggplot(summary_MML_OY,
                   aes(x = H3_class, y = prop,
                       fill = cat_OY_MML)) +
  geom_col() +
  scale_fill_manual(values = dir_colors,
                    name   = "Δ direction\n(Old − Young)") +
  geom_text(
    aes(label = sprintf("%.1f%%\n%d/%d", prop, n, n_total)),
    position = position_stack(vjust = 0.5),
    size     = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 11) +
  labs(
    title = "MML, Old vs Young",
    x     = "H3K27me3 status",
    y     = "Promoters (%)"
  )

#MML: Old vs Old+OSKM
summary_MML_OO <- prom_wide_MML %>%
  mutate(cat_OO_MML = factor(cat_OO_MML, levels = dir_levels)) %>%
  group_by(H3_class, cat_OO_MML) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(H3_class) %>%
  mutate(
    n_total = sum(n),
    prop    = 100 * n / n_total
  ) %>%
  ungroup()

p_MML_OO <- ggplot(summary_MML_OO,
                   aes(x = H3_class, y = prop,
                       fill = cat_OO_MML)) +
  geom_col() +
  scale_fill_manual(values = dir_colors,
                    name   = "Δ direction\n(Old − Old+OSKM)") +
  geom_text(
    aes(label = sprintf("%.1f%%\n%d/%d", prop, n, n_total)),
    position = position_stack(vjust = 0.5),
    size     = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 11) +
  labs(
    title = "MML, Old vs Old+OSKM",
    x     = "H3K27me3 status",
    y     = "Promoters (%)"
  )

#NME: Old vs Young
summary_NME_OY <- prom_wide_NME %>%
  mutate(cat_OY_NME = factor(cat_OY_NME, levels = dir_levels)) %>%
  group_by(H3_class, cat_OY_NME) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(H3_class) %>%
  mutate(
    n_total = sum(n),
    prop    = 100 * n / n_total
  ) %>%
  ungroup()

p_NME_OY <- ggplot(summary_NME_OY,
                   aes(x = H3_class, y = prop,
                       fill = cat_OY_NME)) +
  geom_col() +
  scale_fill_manual(values = dir_colors,
                    name   = "Δ direction\n(Old − Young)") +
  geom_text(
    aes(label = sprintf("%.1f%%\n%d/%d", prop, n, n_total)),
    position = position_stack(vjust = 0.5),
    size     = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 11) +
  labs(
    title = "NME, Old vs Young",
    x     = "H3K27me3 status",
    y     = "Promoters (%)"
  )

#NME: Old vs Old+OSKM
summary_NME_OO <- prom_wide_NME %>%
  mutate(cat_OO_NME = factor(cat_OO_NME, levels = dir_levels)) %>%
  group_by(H3_class, cat_OO_NME) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(H3_class) %>%
  mutate(
    n_total = sum(n),
    prop    = 100 * n / n_total
  ) %>%
  ungroup()

p_NME_OO <- ggplot(summary_NME_OO,
                   aes(x = H3_class, y = prop,
                       fill = cat_OO_NME)) +
  geom_col() +
  scale_fill_manual(values = dir_colors,
                    name   = "Δ direction\n(Old − Old+OSKM)") +
  geom_text(
    aes(label = sprintf("%.1f%%\n%d/%d", prop, n, n_total)),
    position = position_stack(vjust = 0.5),
    size     = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_size = 11) +
  labs(
    title = "NME, Old vs Old+OSKM",
    x     = "H3K27me3 status",
    y     = "Promoters (%)"
  )

#Combine into 4-panel figure
p_four_panel_stacked <- (p_MML_OY | p_MML_OO) /
  (p_NME_OY | p_NME_OO) +
  plot_annotation(tag_levels = "A")

p_four_panel_stacked
