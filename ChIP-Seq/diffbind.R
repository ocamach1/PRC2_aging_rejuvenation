library("DiffBind")
library(csaw)

# Read sample sheet with metadata of our samples 
# (See diffbind documentation to prepare metadata sheet. We include paths for ChIP, input and spikein (reads aligned to Drosophila) bam files, and bed files with the coordinates of the peaks defined by sicer2)
samples <- read.csv("peaksets_H3K27me3.csv")

# Creating a dba object
aging_H3K27me3 <- dba(sampleSheet=samples)

# Creating binding matrix with scores based on read counts for every sample (affinity scores)
aging_H3K27me3 <- dba.count(aging_H3K27me3, summit=2000)

# Normalization with spikein (using bam files with reads that uniquely aligned to the Drosophila genome, specified in the metadata sheet)
aging_H3K27me3 <- dba.normalize(aging_H3K27me3, spikein=TRUE) 

# Defining the comparison we are interested
aging_H3K27me3 <- dba.contrast(aging_H3K27me3, minMembers=2,
                          reorderMeta=list(Condition="Young")) # Here we define the control condition

# Performing the differential analysis 
aging_H3K27me3 <- dba.analyze(aging_H3K27me3, method=DBA_DESEQ2)   
aging_H3K27me3.DB <- dba.report(aging_H3K27me3, th=1)

# Count peaks < 0.05 FDR
count_below_0.05 <- sum(aging_H3K27me3.DB$`FDR` < 0.05)

# Save the dataframe with differentially enriched peaks
write.csv(aging_H3K27me3.DB, file = "H3K27me3_aging_differential_peaks.csv", row.names = FALSE)
