# Overlap of genomic regions

# Required R packages: 
- GenomicRanges
- regioneR

library(GenomicRanges)
library(regioneR)

# Upload coordinates for DNA methylation blocks:
dmrseq_blocks <- read.csv("Dataset_EV8.xlsx", skip=2)
# Get block DMRs:
dmrseq_blocks <- dmrseq_blocks[dmrseq_blocks$qval < 0.05, ]
dmrseq_blocks <- dmrseq_blocks[dmrseq_blocks$L > 100, ]
dmrseq_blocks_gr <- GRanges(
  seqnames = dmrseq_blocks$seqnames,
  ranges   = IRanges(start = dmrseq_blocks$start, end = dmrseq_blocks$end)
)

# Upload coordinates for H3K9me2 LOCKs:
H3K9me2_peaks <- read.csv("Dataset_EV9.xlsx", skip=2)
H3K9me2_gr <- GRanges(
  seqnames = H3K9me2_peaks$seqnames,
  ranges   = IRanges(start = H3K9me2_peaks$start, end = H3K9me2_peaks$end)
)

# Assess significance of the overlap of genomic regions using the R package regioneR, based on permutationt testing:
# Define the mouse genome as canonical chromosomes:
mm10_genome <- filterChromosomes(getGenome("mm10"), keep.chr=c("chr1", "chr2", "chr3","chr4", "chr5", "chr6","chr7",
                                 "chr8", "chr9","chr10", "chr11", "chr12","chr13", "chr14", "chr15",
                                 "chr16", "chr17","chr18", "chr19", "chrX","chrY")) 

set.seed(123)
pt <- permTest(
  A = dmrseq_blocks_gr,
  B = H3K9me2_gr,
  randomize.function = randomizeRegions,
  genome = mm10_genome,      
  ntimes = 1000,
  evaluate.function = numOverlaps,
  alternative = "greater", count.once=TRUE
)
