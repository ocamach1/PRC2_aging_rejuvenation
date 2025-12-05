library("DiffBind")
library(csaw)

# Reading sample sheet with metadata of our samples. See DiffBind documentation to prepare metadata sheet. 
# Briefly, the metadata sheet contains the paths for:
# - IP bam file (containing reads uniquely aligned to the mouse genome, derived from the IP sample - DNA bound to the target of interest).
# - Input bam file (containing reads uniquely aligned to the mouse genome, derived from the Input sample - Total chromatin DNA with no antibody pull-down).
# - Spike-in bam file (containing reads uniquely aligning to Drosophila genome, derived from the IP sample). 
# - BED file containing the coordinates of the peaks defined by SICER2.

samples <- read.csv("peaksets_H3K27me3.csv")

# Creating a dba object
aging_H3K27me3 <- dba(sampleSheet=samples)

# Creating binding matrix with scores based on read counts for every sample (affinity scores)
aging_H3K27me3 <- dba.count(aging_H3K27me3, summit=FALSE)

# Normalization using spike-in normalization
aging_H3K27me3 <- dba.normalize(aging_H3K27me3, spikein=TRUE) 

# Defining the contrast we are interested (Old vs Young)
aging_H3K27me3 <- dba.contrast(aging_H3K27me3, minMembers=2,
                          reorderMeta=list(Condition="Young")) # Here we define 'Young' as the control condition

# Performing the differential analysis 
aging_H3K27me3 <- dba.analyze(aging_H3K27me3, method=DBA_DESEQ2)   
aging_H3K27me3.DB <- dba.report(aging_H3K27me3, th=1)

# Count peaks < 0.05 FDR
count_below_0.05 <- sum(aging_H3K27me3.DB$`FDR` < 0.05)

# Save the dataframe with differentially enriched peaks
write.csv(aging_H3K27me3.DB, file = "H3K27me3_aging_differential_peaks.csv", row.names = FALSE)
