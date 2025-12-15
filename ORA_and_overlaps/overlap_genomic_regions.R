# Overlap of genomic regions

# Required R packages: 
- GenomicRanges
- regioneR

library(GenomicRanges)
library(regioneR)


# Read H3K9me2 consensus peaks:
H3K9me2_peaks <- read.csv("/Users/oscarcamacho/Documents/Feinberg_lab/ChIP_aging_epidermis/H3K9me2/diffbind/EDD_spikein_norm/summit_FALSE_filter0_spikein_norm_no_blacklist.csv", header=TRUE)
H3K9me2_peaks <- H3K9me2_peaks[ !H3K9me2_peaks$seqnames %in% c("chrX", "chrY"), ]

#Convert peak coordinates to granges object
H3K9me2_gr <- GRanges(
  seqnames = H3K9me2_peaks$seqnames,
  ranges   = IRanges(start = H3K9me2_peaks$start, end = H3K9me2_peaks$end)
)
H3K9me2_gr




# From dmrseq (Block settings):
dmrseq_blocks <- read.csv("/Users/oscarcamacho/Documents/Feinberg_lab/ChIP_aging_epidermis/WGBS_aging_epidermis/dmrseq_Kenna/OldvYoung_block_regions_dmrseq_default_maxPerms100.csv")
#dmrseq_blocks <- dmrseq_blocks[ !dmrseq_blocks$seqnames %in% c("chrX", "chrY"), ]

# Subset blocks with qval less than 0.05 and more than 100 Cpgs:
dmrseq_blocks <- dmrseq_blocks[dmrseq_blocks$qval < 0.05, ]
dmrseq_blocks <- dmrseq_blocks[dmrseq_blocks$L > 100, ]
dmrseq_blocks_gr <- GRanges(
  seqnames = dmrseq_blocks$seqnames,
  ranges   = IRanges(start = dmrseq_blocks$start, end = dmrseq_blocks$end)
)

# Check proportion of shared elements (dmrseq blocks):
hits <- findOverlaps(H3K9me2_gr,dmrseq_blocks_gr)
# Proportion of H3K9me2 peaks that overlap any DMR
n_H3K9me2        <- length(H3K9me2_gr)
n_H3K9me2_olap   <- length(unique(queryHits(hits)))
prop_H3K9me2_olap <- n_H3K9me2_olap / n_H3K9me2
prop_H3K9me2_olap # 69.31% with all chr. 68.49% without sex chr
# Proportion of DMRs that overlap any H3K9me2 peak
n_dmr        <- length(dmrseq_blocks_gr)
n_dmr_olap   <- length(unique(subjectHits(hits)))
prop_dmr_olap <- n_dmr_olap / n_dmr
n_dmr_olap
n_dmr
prop_dmr_olap # 87.23 % with all chr. 91.14% without sex chr



mm10_genome <- filterChromosomes(getGenome("mm10"), keep.chr=c("chr1", "chr2", "chr3","chr4", "chr5", "chr6","chr7",
                                 "chr8", "chr9","chr10", "chr11", "chr12","chr13", "chr14", "chr15",
                                 "chr16", "chr17","chr18", "chr19", "chrX","chrY")) 



set.seed(123)

pt <- permTest(
  A = dmrseq_blocks_gr,
  B = H3K9me2_gr,
  randomize.function = randomizeRegions,
  genome = mm10_genome,      # leave this as the full genome object
  ntimes = 1000,
  evaluate.function = numOverlaps,
  alternative = "greater", count.once=TRUE
)
pt

